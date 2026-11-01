"""Reconstruct shared generic Gröbner derivations from numerical support proposals.

Numerical words propose supports; scalar values select compatible nodes at the
anchor. Local generic kernel solves recover a polynomial
and its derivation together over the original field, retaining references to
previously reconstructed nodes. Only independent exact completion checking can
turn this finite support search into a certified generic reduction.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from graphlib import TopologicalSorter

from .algebra import AlgebraPresentation, AlgebraElement
from .matrices import ExactMatrix
from .fields import RationalFunctionField
from .reconstruction import reconstruct_kernel
from .verification.derivations import (DerivationGraph, DerivationTerm,
                                       _Budget, _positive_limit, verify_derivations)
from .verification.groebner import (GroebnerCertificate, ProofTerm, VerifiedReduction,
                                   verify_groebner)
from .verification.linear import verify_kernel
from .verification.working_field import WorkingField

__all__ = ["GroebnerProduction", "SupportFailure", "reconstruct_groebner"]


class SupportFailure(RuntimeError):
    """A numerical proof proposal cannot be lifted with the supplied supports."""


@dataclass(frozen=True, slots=True)
class GroebnerProduction:
    """A reconstructed candidate that passed exact generic completion checks."""

    certificate: GroebnerCertificate
    reduction: VerifiedReduction
    statistics: dict


def _graph(certificate):
    return DerivationGraph(certificate.presentation_hash, certificate.basis,
                           tuple(tuple(DerivationTerm(term.coefficient, term.left,
                                                     "input" if term.kind == "input" else "node",
                                                     term.index, term.right) for term in row)
                                 for row in certificate.derivations))


def _select_node(matrix, kernel, support, contexts, numeric, point, options, max_products):
    """Select a generic identity compatible with the observed monic node.

    A local kernel may be underdetermined. Choosing an arbitrary vector can
    change a parent node at the anchor and make later shared references
    unusable. Numerical coefficients select a constant combination of exact
    generic kernel vectors; they are never interpolated as generic functions.
    """
    field = matrix.field
    arithmetic = WorkingField(field)
    entries = {(row, column): arithmetic.from_element(value)
               for row, column, value in kernel.kernel.entries}
    leading = max(support, key=lambda word: (len(word), word))
    anchor = len(contexts) + support.index(leading)
    observed = dict(numeric.terms)
    lead_value = observed[leading]

    def rational(value):
        return Fraction(value.numerator[0][1]) if value.numerator else Fraction(0)

    target = [rational(observed[word] / lead_value) if word in observed else Fraction(0)
              for word in support]
    evaluated = {}
    available = []
    for column in range(kernel.kernel.ncols):
        try:
            values = {(row, column): rational(value.specialize(point))
                      for row, col, value in kernel.kernel.entries
                      if col == column and row >= len(contexts)}
        except ZeroDivisionError:
            continue
        available.append(column)
        evaluated.update(values)
    if not available:
        raise SupportFailure("no generic kernel column has regular node coefficients at the anchor")
    for column in available:
        scale = evaluated.get((anchor, column), 0)
        if scale and all(evaluated.get((len(contexts) + row, column), 0) == scale * value
                         for row, value in enumerate(target)):
            inverse = entries[anchor, column]
            return [arithmetic.to_element(entries.get((row, column), arithmetic.zero) / inverse)
                    for row in range(matrix.ncols)]
    ground = RationalFunctionField()
    indices = {column: index for index, column in enumerate(available)}
    combination = ExactMatrix(ground, len(support), len(available) + 1,
                       tuple((row - len(contexts), indices[column], ground.scalar(value))
                             for (row, column), value in evaluated.items())
                       + tuple((row, len(available), ground.scalar(-value))
                               for row, value in enumerate(target)))
    solution = reconstruct_kernel(combination, max_products=max_products, **options)
    verify_kernel(combination, solution, max_products=max_products)
    weights = {(row, column): rational(value) for row, column, value in solution.kernel.entries}
    choices = [(column, value) for (row, column), value in weights.items()
               if row == len(available)]
    if not choices:
        raise SupportFailure("numerical node has no anchor-compatible generic proof in this support")
    column, scale = min(choices)
    selected = [(available[j], arithmetic.from_element(field.scalar(value / scale)))
                for (j, col), value in weights.items() if col == column and j < len(available)]
    coefficients = [sum((entries.get((row, j), arithmetic.zero) * weight
                         for j, weight in selected), arithmetic.zero)
                    for row in range(matrix.ncols)]
    return [arithmetic.to_element(value / coefficients[anchor]) for value in coefficients]


def reconstruct_groebner(presentation: AlgebraPresentation, sample: GroebnerCertificate, *,
                         point: Mapping[str, int | Fraction], max_nodes=100000,
                         max_operations=10000000, max_products=10000000,
                         max_reductions=2000000, max_words=100000,
                         max_primes=5, max_degree=24, threads=1, start_prime=0,
                         shift=None) -> GroebnerProduction:
    """Lift a numerical proof shape using exact generic equations.

    ``sample`` must be bound to the full literal specialization at ``point``;
    compacted solver input references must already have been remapped. Every
    numerical derivation is independently verified. Its words and contextual
    references propose finite local systems, never a completeness assumption.

    Each system contains contextual input/parent columns followed by negative
    monomials drawn from those columns up to the numerical leading word. This
    restores tails that cancelled at the sample, while requiring the proposed
    leading word to survive generically. A kernel vector with nonzero leading-word
    coordinate is normalized to produce a monic generic node. Numerical node
    coefficients select an anchor-compatible combination when a kernel is
    underdetermined; generic coefficients come from exact generic equations.
    Inadequate support,
    exhausted bounds and failed completion checks give no dimension verdict.

    ``max_operations`` separately caps numerical graph verification, aggregate
    generic matrix placements, and final generic graph verification. Kernel
    product/reconstruction bounds apply per local solve; ``max_reductions`` and
    ``max_words`` bound the final completion checker. These are not process RAM
    or wall-clock limits. A completed reduction can have dimension=None when
    finite standard-word enumeration exceeds its bound.
    """
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(sample, GroebnerCertificate):
        raise TypeError("generic production requires an AlgebraPresentation and GroebnerCertificate")
    field = presentation.algebra.field
    if field.characteristic:
        raise ValueError("generic f4ncgb production currently requires characteristic zero")
    if (not isinstance(point, Mapping) or set(point) != set(field.parameters)
            or any(type(value) is not int and not isinstance(value, Fraction) for value in point.values())):
        raise ValueError("point must assign every declared parameter an exact rational value")
    _positive_limit(max_nodes, "max_nodes")
    budget = _Budget(max_operations)
    for name, value in (("max_products", max_products), ("max_reductions", max_reductions),
                        ("max_words", max_words)):
        _positive_limit(value, name)
    for name, value, minimum in (("max_primes", max_primes, 1), ("max_degree", max_degree, 1),
                                 ("threads", threads, 1), ("start_prime", start_prime, 0)):
        if type(value) is not int or value < minimum:
            raise ValueError(f"{name} must be an integer at least {minimum}")
    if max_primes > 2047 or start_prime + max_primes >= 2048 or max_degree >= 2**16 or threads >= 2**32:
        raise ValueError("reconstruction options exceed FiniteFlow native limits")
    if shift is not None:
        shift = tuple(shift)
        if (len(shift) != len(field.parameters)
                or any(type(x) is not int or not 0 <= x < 2**64 for x in shift)):
            raise ValueError("the shift needs one unsigned exact integer per parameter")
    specialized = presentation.specialize(point)
    verify_derivations(specialized, _graph(sample), max_nodes=max_nodes,
                       max_operations=max_operations)
    dependencies = {i: {term.index for term in row if term.kind == "basis"}
                    for i, row in enumerate(sample.derivations)}
    order = tuple(TopologicalSorter(dependencies).static_order())
    nodes: dict[int, AlgebraElement] = {}
    derivations = {}
    counts = {"nodes": len(sample.basis), "systems": 0, "matrix_rows": 0,
              "matrix_columns": 0, "matrix_entries": 0}
    algebra = presentation.algebra
    for index in order:
        leading = max((word for word, _ in sample.basis[index].terms), key=lambda word: (len(word), word))
        # Repeated raw references propose the same column. Preserve first-use
        # order without retaining duplicate unknowns or numerical coefficients.
        contexts = tuple(dict.fromkeys((term.left, term.kind, term.index, term.right)
                                        for term in sample.derivations[index]))
        vectors = []
        for left, kind, reference, right in contexts:
            budget.consume()
            source = presentation.relations[reference] if kind == "input" else nodes[reference]
            vector = {}
            for word, coefficient in source.terms:
                budget.consume()
                vector[left + word + right] = coefficient
            vectors.append(vector)
        sampled_support = tuple(word for word, _ in sample.basis[index].terms)
        expanded_support = tuple(sorted({word for vector in vectors for word in vector
                                if (len(word), word) <= (len(leading), leading)}
                               | set(sampled_support)))
        def solve_node(support):
            budget.consume(len(support) + sum(map(len, vectors)))
            columns = vectors + [{word: -field.one} for word in support]
            words = tuple(sorted(set().union(*(set(vector) for vector in columns)),
                                 key=lambda word: (len(word), word)))
            rows = {word: row for row, word in enumerate(words)}
            matrix = ExactMatrix(field, len(words), len(columns),
                                 tuple((rows[word], column, coefficient)
                                       for column, vector in enumerate(columns)
                                       for word, coefficient in vector.items()))
            kernel = reconstruct_kernel(matrix, max_primes=max_primes, max_degree=max_degree,
                                         threads=threads, start_prime=start_prime, shift=shift,
                                         max_products=max_products)
            verify_kernel(matrix, kernel, max_products=max_products)
            selection_options = {"max_primes": max_primes, "max_degree": max_degree,
                                 "threads": threads, "start_prime": start_prime}
            counts["systems"] += 1
            counts["matrix_rows"] += matrix.nrows
            counts["matrix_columns"] += matrix.ncols
            counts["matrix_entries"] += len(matrix.entries)
            return _select_node(matrix, kernel, support, contexts, sample.basis[index],
                                point, selection_options, max_products)

        support = sampled_support
        try:
            coefficients = solve_node(support)
        except SupportFailure:
            if expanded_support == sampled_support:
                raise
            support = expanded_support
            coefficients = solve_node(support)
        nodes[index] = algebra.element(dict(zip(support, coefficients[len(contexts):])))
        derivations[index] = tuple(ProofTerm(value, left, kind, reference, right)
                                   for (left, kind, reference, right), value
                                   in zip(contexts, coefficients[:len(contexts)]) if not value.is_zero())
    certificate = GroebnerCertificate(presentation.content_hash(),
                                      tuple(nodes[i] for i in range(len(nodes))),
                                      tuple(derivations[i] for i in range(len(nodes))))
    verify_derivations(presentation, _graph(certificate), max_nodes=max_nodes,
                       max_operations=max_operations)
    reduction = verify_groebner(presentation, certificate, max_reductions=max_reductions,
                                max_words=max_words, prune=True)
    counts["matrix_placements"] = budget.used
    return GroebnerProduction(certificate, reduction, counts)
