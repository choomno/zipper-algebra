"""Deterministic parameter samples backed by validated solver-run records."""
from __future__ import annotations

from collections import Counter
from collections.abc import Collection, Iterator
from concurrent.futures import ThreadPoolExecutor
from typing import TypedDict, cast

from ..cofactor.atoms import AtomSet
from ..cofactor.proof import Basis, Proof, Word
from ..cofactor.skeleton import Anchor, Skeleton
from ..presentation import Presentation
from .config import RunConfig
from .runner import SolverRun, run_point


EXCLUDED_Q = frozenset({2, 4})
Point = tuple[int, int]
ProofShape = tuple[tuple[Word, str, int | str, Word], ...]


class ShapeCensus(TypedDict):
    keys: list[str]
    anchor_shape: dict[str, ProofShape]
    census: dict[str, Counter[ProofShape | None]]
    per_point: dict[Point, dict[str, ProofShape | None]]


def spiral() -> Iterator[Point]:
    """Positive prime pairs in the historical square-spiral order, off the diagonal."""
    primes: list[int] = []
    value = 2
    while True:
        if all(value % prime for prime in primes if prime * prime <= value):
            primes.append(value)
            for other in primes[:-1]:
                yield value, other
            for other in reversed(primes[:-1]):
                yield other, value
        value += 1


def admissible(q: int, t: int, atoms: AtomSet,
               excluded_q: Collection[int] = EXCLUDED_Q) -> bool:
    """Sampling exclusions and zeros of the specified denominator factors."""
    return q not in excluded_q and all(atom.eval(q, t) != 0 for atom in atoms.atoms)


class SampleSet:
    """One presentation/target with parsed proofs from validated run bundles.

    Points are retained in deterministic submission order. Calling the same
    discovery command again reuses the runner's validated content-addressed
    records; arbitrary proof filenames are never interpreted as completed runs.
    """

    def __init__(self, pres: Presentation, target: str, anchor: Anchor, config: RunConfig,
                 atoms: AtomSet | None = None, *, producer: dict | None = None) -> None:
        self.pres, self.target, self.config = pres, target, config
        # tuple preserves the supplied three-coordinate anchor; its inferred
        # variadic result otherwise loses the fixed arity of the input type.
        self.anchor = cast(Anchor, tuple(anchor))
        self.atoms = AtomSet([]) if atoms is None else atoms
        self.atoms.validate()
        self._runs: dict[Point, SolverRun] = {}
        self._bases: dict[Point, Basis] = {}
        self._proofs: dict[Point, Proof] = {}
        self.producer = producer

    @property
    def points(self) -> list[Point]:
        return list(self._runs)

    @property
    def char(self) -> int:
        return self.anchor[2]

    def _remember(self, point: Point, result: SolverRun) -> None:
        self._runs[point] = result
        self._bases[point] = Basis.from_text(result.files["basis.txt"].decode(), self.char,
                                            self.pres.var_names())
        self._proofs[point] = Proof.from_text(result.files["proof.txt"].decode(), self.pres.var_names())

    def run_point(self, q: int, t: int, *, refresh: bool = False) -> SolverRun:
        result = run_point(self.pres, self.target, (q, t, self.char), self.config,
                           refresh=refresh, producer=self.producer)
        self._remember((q, t), result)
        return result

    def run_anchor(self, *, refresh: bool = False) -> SolverRun:
        if not admissible(*self.anchor[:2], self.atoms):
            raise ValueError("anchor is excluded by the sampling policy or lies on a denominator zero")
        return self.run_point(*self.anchor[:2], refresh=refresh)

    def basis(self, q: int, t: int) -> Basis:
        return self._bases[q, t]

    def proof(self, q: int, t: int) -> Proof:
        return self._proofs[q, t]

    def anchor_basis(self) -> Basis:
        return self.basis(*self.anchor[:2])

    def anchor_proof(self) -> Proof:
        return self.proof(*self.anchor[:2])

    def reanchor(self, point: Point) -> None:
        self.basis(*point)  # A replacement anchor must already have validated proof data.
        self.anchor = (*point, self.char)

    def next_points(self, amount: int) -> list[Point]:
        if type(amount) is not int or amount < 0:
            raise ValueError("sample growth must be a nonnegative integer")
        selected: list[Point] = []
        points = spiral()
        while len(selected) < amount:
            point = next(points)
            if point not in self._runs and admissible(*point, self.atoms):
                selected.append(point)
        return selected

    def grow(self, amount: int, *, jobs: int = 1, refresh: bool = False) -> list[Point]:
        if type(jobs) is not int or jobs < 1:
            raise ValueError("jobs must be a positive integer")
        points = self.next_points(amount)
        if points:
            def one(point: Point) -> SolverRun:
                return run_point(self.pres, self.target, (*point, self.char), self.config,
                                 refresh=refresh, producer=self.producer)
            with ThreadPoolExecutor(max_workers=jobs) as executor:
                results = list(executor.map(one, points))
            for point, result in zip(points, results):
                self._remember(point, result)
        return points

    def shape_census(self, skeleton: Skeleton) -> ShapeCensus:
        keys = skeleton.basis_keys
        basis, proof = self.anchor_basis(), self.anchor_proof()
        anchor_shapes = {key: proof.shape(basis.index_by_key[key], basis) for key in keys}
        census: dict[str, Counter[ProofShape | None]] = {key: Counter() for key in keys}
        per_point: dict[Point, dict[str, ProofShape | None]] = {}
        for point in self.points:
            basis, proof = self.basis(*point), self.proof(*point)
            shapes: dict[str, ProofShape | None] = {}
            for key in keys:
                index = basis.index_by_key.get(key)
                shape = None if index is None else proof.shape(index, basis)
                census[key][shape] += 1
                shapes[key] = shape
            per_point[point] = shapes
        return {"keys": keys, "anchor_shape": anchor_shapes,
                "census": census, "per_point": per_point}

    def modal_shapes(self, skeleton: Skeleton) -> tuple[dict[str, ProofShape], list[str]]:
        census = self.shape_census(skeleton)
        modal: dict[str, ProofShape] = {}
        off: list[str] = []
        for key in census["keys"]:
            counts = census["census"][key]
            best = max(counts, key=lambda shape: (counts[shape], shape is not None))
            modal[key] = census["anchor_shape"][key] if best is None else best
            if modal[key] != census["anchor_shape"][key]:
                off.append(key)
        return modal, off

    def stable(self, skeleton: Skeleton) -> list[Point]:
        census = self.shape_census(skeleton)
        modal = self.modal_shapes(skeleton)[0]
        return [point for point in self.points
                if all(census["per_point"][point][key] == modal[key]
                       for key in census["keys"])]
