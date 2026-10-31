"""Derive contextual word equations from a solver proof and certify a cofactor.

The proof graph supplies candidate contexts only. Its numeric coefficients are
ignored, including cancellations between distinct derivation paths. A full
kernel certificate and a separate free-algebra identity decide acceptance over
Q(q,t); a finite graph is not a completeness theorem for this search.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from graphlib import CycleError, TopologicalSorter
from math import lcm

from ..fields import RationalFunctionField
from ..matrices import ExactMatrix
from ..presentation import Presentation
from ..reconstruction import reconstruct_kernel
from ..verification.linear import KernelCertificate, verify_kernel
from .atoms import Atom, AtomSet
from .cofactor_io import Cofactor, symbolic_target, validate_word
from .proof import Basis, Proof, Word
from .skeleton import Skeleton
from .verify import check

Context = tuple[Word, int, Word]
_FIELD = RationalFunctionField(("q", "t"))


@dataclass(frozen=True)
class ContextSystem:
    """Column j is left*relation[index]*right; the final column is -target."""

    presentation_hash: str
    target: str
    contexts: tuple[Context, ...]
    words: tuple[Word, ...]
    matrix: ExactMatrix
    statistics: dict

    def to_json(self):
        return {"schema": "zipper-context-system-v1",
                "presentation_hash": self.presentation_hash, "target": self.target,
                "contexts": [[list(left), index, list(right)]
                             for left, index, right in self.contexts],
                "word_rows": [list(word) for word in self.words],
                "target_column": len(self.contexts),
                "matrix_hash": self.matrix.content_hash(), "statistics": self.statistics}


@dataclass(frozen=True)
class Discovery:
    system: ContextSystem
    kernel: KernelCertificate
    cofactor: Cofactor
    exact: dict


def build_system(pres: Presentation, target: str, skeleton: Skeleton,
                 basis: Basis, proof: Proof, *, max_contexts=100000,
                 max_products=10000000) -> ContextSystem:
    """Expand a finite acyclic proof graph structurally and build [A|-h].

    ``max_contexts`` bounds the sum of stored intermediate context sets and the
    output set. ``max_products`` bounds contextual substitutions plus nonzero
    word-coefficient placements in the augmented matrix. Exhaustion is an
    inconclusive search failure, never evidence of nonmembership.
    """
    for value in (max_contexts, max_products):
        if type(value) is not int or value < 1:
            raise ValueError("context and product budgets must be positive integers")
    if (skeleton.presentation_hash != pres.content_hash() or skeleton.target != target
            or not skeleton.in_ideal):
        raise ValueError("a zero-remainder skeleton for this presentation and target is required")
    if (basis.char != 0 or len(proof) != len(basis)
            or basis.var_names != pres.var_names()):
        raise ValueError("proof and characteristic-zero basis must match the presentation")
    relations, _ = pres.symbolic_ncpolys()

    def reference(left, kind, index, right):
        validate_word(left, pres.n)
        validate_word(right, pres.n)
        if (kind not in ("basis", "input") or type(index) is not int
                or not 0 <= index < (len(basis) if kind == "basis" else len(relations))):
            raise ValueError("proof reference is outside this presentation or basis")
        return tuple(left), kind, index, tuple(right)

    rules = [tuple(reference(left, kind, index, right)
                   for _coefficient, left, kind, index, right in line)
             for line in proof.terms]
    try:
        roots = [reference(step["left"], "basis", basis.index(step["basis_key"]), step["right"])
                 for step in skeleton.steps]
    except (KeyError, TypeError) as error:
        raise ValueError("skeleton steps must name existing basis keys and word contexts") from error
    dependencies = {i: {index for _, kind, index, _ in line if kind == "basis"}
                    for i, line in enumerate(rules)}
    needed: set[int] = set()
    pending = [index for _, _, index, _ in roots]
    while pending:
        index = pending.pop()
        if index not in needed:
            needed.add(index)
            pending.extend(dependencies[index])
    try:
        order = tuple(TopologicalSorter({index: dependencies[index] for index in needed}).static_order())
    except CycleError as error:
        raise ValueError("cyclic proof graph") from error
    expanded: dict[int, set[Context]] = {}
    work = retained = 0

    def consume(amount):
        nonlocal work
        work += amount
        if work > max_products:
            raise RuntimeError("context-system product budget exhausted")

    def substitute(line):
        out: set[Context] = set()
        for left, kind, index, right in line:
            children = {((), index, ())} if kind == "input" else expanded[index]
            for u, relation, v in children:
                consume(1)
                out.add((left + u, relation, v + right))
                if retained + len(out) > max_contexts:
                    raise RuntimeError("context-system storage budget exhausted")
        return out

    for index in order:
        expanded[index] = substitute(rules[index])
        retained += len(expanded[index])
    contexts = tuple(sorted(substitute(roots), key=lambda term: (len(term[0]) + len(term[2]), term)))
    target_poly = symbolic_target(pres, target)
    vectors = []
    for left, index, right in contexts:
        vector = {left + word + right: coefficient
                  for word, coefficient in relations[index].items()}
        consume(sum(bool(coefficient) for coefficient in vector.values()))
        vectors.append(vector)
    consume(sum(bool(coefficient) for coefficient in target_poly.values()))
    words = tuple(sorted(set(target_poly).union(*(set(vector) for vector in vectors)),
                         key=lambda word: (len(word), word)))
    row_indices = {word: index for index, word in enumerate(words)}
    entries = [(row_indices[word], column, _FIELD.polynomial(coefficient))
               for column, vector in enumerate(vectors) for word, coefficient in vector.items()]
    entries.extend((row_indices[word], len(contexts), -_FIELD.polynomial(coefficient))
                   for word, coefficient in target_poly.items())
    matrix = ExactMatrix(_FIELD, len(words), len(contexts) + 1, tuple(entries))
    statistics = {"needed_basis_nodes": len(needed), "distinct_contexts": len(contexts),
                  "stored_contexts": retained + len(contexts), "products": work,
                  "grammar_edges": sum(len(rules[index]) for index in needed) + len(roots)}
    return ContextSystem(pres.content_hash(), target, contexts, words, matrix, statistics)


def _cofactor(pres, target, contexts, coefficients, anchor):
    """Clear a true polynomial LCM, without an anticipated denominator family."""
    from sympy import QQ
    from sympy.polys.rings import ring

    polynomial_ring, _, _ = ring("q,t", QQ)
    fractions = [polynomial_ring.from_dict(dict(value.numerator)).cancel(
                     polynomial_ring.from_dict(dict(value.denominator))) for value in coefficients]
    numerators = [numerator for numerator, _ in fractions]
    denominators = [bottom for _, bottom in fractions]
    denominator = polynomial_ring.one
    for numerator, bottom in zip(numerators, denominators):
        if numerator:
            denominator = denominator.lcm(bottom)
    denominator *= lcm(*(int(value.denominator) for value in denominator.values()))
    integral = {monomial: int(value) for monomial, value in denominator.items()}
    atoms = AtomSet([], []) if denominator == polynomial_ring.one else AtomSet([Atom("LCM", integral)], [1])
    terms = []
    labels = pres.labels()
    for (left, index, right), numerator, bottom in zip(contexts, numerators, denominators):
        if numerator:
            scaled = numerator * denominator.exquo(bottom)
            scaled_numerator = {monomial: Fraction(int(value.numerator), int(value.denominator))
                            for monomial, value in scaled.items()}
            terms.append((left, labels[index], right, scaled_numerator, atoms.zero() if not len(atoms) else (1,)))
    return Cofactor.from_presentation(pres, target, atoms, terms,
                                     "proof-graph FiniteFlow reconstruction", anchor=anchor)


def reconstruct_cofactor(pres: Presentation, target: str, system: ContextSystem, *,
                         anchor=None, max_primes=5, max_degree=24, threads=1,
                         start_prime=0, shift=None, max_products=10000000) -> Discovery:
    """Reconstruct a kernel, normalize a target-coordinate vector, and check it.

    A zero target row throughout the complete kernel means the finite supplied
    context family cannot express the target. It makes no claim about the ideal
    generated by all contextual multiples of the defining relations.
    """
    if (system.presentation_hash != pres.content_hash() or system.target != target
            or system.matrix.field != _FIELD
            or system.matrix.ncols != len(system.contexts) + 1):
        raise ValueError("context system does not match this presentation and target")
    certificate = reconstruct_kernel(system.matrix, max_primes=max_primes, max_degree=max_degree,
                                     threads=threads, start_prime=start_prime, shift=shift,
                                     max_products=max_products)
    try:
        verify_kernel(system.matrix, certificate, max_products=max_products)
    except (ValueError, ZeroDivisionError, TypeError) as error:
        raise RuntimeError("reconstructed kernel failed independent exact verification") from error
    target_row = len(system.contexts)
    values = {(row, column): value for row, column, value in certificate.kernel.entries}
    candidates = [(column, value) for (row, column), value in values.items() if row == target_row]
    if not candidates:
        raise RuntimeError("the finite proof-graph context family cannot express this target")
    column, scale = min(candidates, key=lambda pair: pair[0])
    coefficients = [values.get((row, column), _FIELD.zero) / scale for row in range(target_row)]
    cofactor = _cofactor(pres, target, system.contexts, coefficients, anchor)
    exact = check(cofactor, pres, target, verbose=False)
    if not exact["pass"]:
        raise RuntimeError("independent exact free-algebra identity check failed")
    return Discovery(system, certificate, cofactor, exact)
