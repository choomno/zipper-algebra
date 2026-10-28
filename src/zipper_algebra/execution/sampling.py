"""Deterministic parameter samples backed by validated solver-run records."""
from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor

from ..cofactor.atoms import AtomSet
from ..cofactor.proof import Basis, Proof
from .runner import run_point


EXCLUDED_Q = frozenset({2, 4})


def spiral():
    """Positive prime pairs in the historical square-spiral order, off the diagonal."""
    primes = []
    value = 2
    while True:
        if all(value % prime for prime in primes if prime * prime <= value):
            primes.append(value)
            for other in primes[:-1]:
                yield value, other
            for other in reversed(primes[:-1]):
                yield other, value
        value += 1


def admissible(q, t, atoms, excluded_q=EXCLUDED_Q):
    """Sampling exclusions and zeros of the specified denominator factors."""
    return q not in excluded_q and all(atom.eval(q, t) != 0 for atom in atoms.atoms)


class SampleSet:
    """One presentation/target with parsed proofs from validated run bundles.

    Points are retained in deterministic submission order. Calling the same
    discovery command again reuses the runner's validated content-addressed
    records; arbitrary proof filenames are never interpreted as completed runs.
    """

    def __init__(self, pres, target, anchor, config, atoms=None, *, producer=None):
        self.pres, self.target, self.anchor, self.config = pres, target, tuple(anchor), config
        self.atoms = AtomSet.from_json(AtomSet.DEFAULT.to_json()) if atoms is None else atoms
        self.atoms.validate()
        self._runs, self._bases, self._proofs = {}, {}, {}
        self.producer = producer

    @property
    def points(self):
        return list(self._runs)

    @property
    def char(self):
        return self.anchor[2]

    def _remember(self, point, result):
        self._runs[point] = result
        self._bases[point] = Basis.from_text(result.files["basis.txt"].decode(), self.char,
                                            self.pres.var_names())
        self._proofs[point] = Proof.from_text(result.files["proof.txt"].decode(), self.pres.var_names())

    def run_point(self, q, t, *, refresh=False):
        result = run_point(self.pres, self.target, (q, t, self.char), self.config,
                           refresh=refresh, producer=self.producer)
        self._remember((q, t), result)
        return result

    def run_anchor(self, *, refresh=False):
        if not admissible(*self.anchor[:2], self.atoms):
            raise ValueError("anchor is excluded by the sampling policy or lies on a denominator zero")
        return self.run_point(*self.anchor[:2], refresh=refresh)

    def basis(self, q, t):
        return self._bases[q, t]

    def proof(self, q, t):
        return self._proofs[q, t]

    def anchor_basis(self):
        return self.basis(*self.anchor[:2])

    def anchor_proof(self):
        return self.proof(*self.anchor[:2])

    def reanchor(self, point):
        self.basis(*point)  # A replacement anchor must already have validated proof data.
        self.anchor = (*point, self.char)

    def next_points(self, amount):
        if type(amount) is not int or amount < 0:
            raise ValueError("sample growth must be a nonnegative integer")
        selected = []
        points = spiral()
        while len(selected) < amount:
            point = next(points)
            if point not in self._runs and admissible(*point, self.atoms):
                selected.append(point)
        return selected

    def grow(self, amount, *, jobs=1, refresh=False):
        if type(jobs) is not int or jobs < 1:
            raise ValueError("jobs must be a positive integer")
        points = self.next_points(amount)
        if points:
            def one(point):
                return run_point(self.pres, self.target, (*point, self.char), self.config,
                                 refresh=refresh, producer=self.producer)
            with ThreadPoolExecutor(max_workers=jobs) as executor:
                results = list(executor.map(one, points))
            for point, result in zip(points, results):
                self._remember(point, result)
        return points

    def shape_census(self, skeleton):
        keys = skeleton.basis_keys
        basis, proof = self.anchor_basis(), self.anchor_proof()
        anchor_shapes = {key: proof.shape(basis.index_by_key[key], basis) for key in keys}
        census = {key: Counter() for key in keys}
        per_point = {}
        for point in self.points:
            basis, proof = self.basis(*point), self.proof(*point)
            shapes = {}
            for key in keys:
                index = basis.index_by_key.get(key)
                shape = None if index is None else proof.shape(index, basis)
                census[key][shape] += 1
                shapes[key] = shape
            per_point[point] = shapes
        return {"keys": keys, "anchor_shape": anchor_shapes,
                "census": census, "per_point": per_point}

    def modal_shapes(self, skeleton):
        census = self.shape_census(skeleton)
        modal, off = {}, []
        for key in census["keys"]:
            counts = census["census"][key]
            best = max(counts, key=lambda shape: (counts[shape], shape is not None))
            modal[key] = census["anchor_shape"][key] if best is None else best
            if modal[key] != census["anchor_shape"][key]:
                off.append(key)
        return modal, off

    def stable(self, skeleton):
        census = self.shape_census(skeleton)
        modal = self.modal_shapes(skeleton)[0]
        return [point for point in self.points
                if all(census["per_point"][point][key] == modal[key]
                       for key in census["keys"])]
