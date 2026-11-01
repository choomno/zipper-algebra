"""Reconstruct one explicit original-input ideal witness with FiniteFlow.

A derivation graph proposes word contexts without propagating its coefficients.
Native sparse affine solving reconstructs one particular vector. Acceptance is
an independent exact identity against every original defining polynomial; no
kernel, rank, basis-completeness, or shortest-proof assertion is made.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import importlib
from typing import Any, Iterable, TypeAlias

from ..algebra import AlgebraPresentation, CofactorIdentity
from ..fields import FieldElement, RationalFunctionField
from ..reconstruction import _NATIVE_LOCK, _native_errors
from ..verification.derivations import DerivationGraph, VerificationLimit, _Budget, _positive_limit
from ..verification.working_field import WorkingField
from ._contexts import Context, _collect_contexts, _validate_graph

__all__ = ["CofactorReconstruction", "reconstruct_derivation"]

Rows: TypeAlias = list[list[tuple[int, FieldElement]]]


@dataclass(frozen=True, slots=True)
class CofactorReconstruction:
    """A checked original-input witness and bounded integer construction counts."""

    identity: CofactorIdentity
    _statistics: tuple[tuple[str, int], ...]

    @property
    def statistics(self) -> dict[str, int]:
        """Counts: candidate_contexts before solving, contexts in the checked witness."""
        return dict(self._statistics)


def _options(field: RationalFunctionField, max_primes: int, max_degree: int,
             threads: int, start_prime: int, shift: Iterable[int] | None) -> tuple[int, ...] | None:
    for name, value in (("max_primes", max_primes), ("max_degree", max_degree), ("threads", threads)):
        _positive_limit(value, name)
    if type(start_prime) is not int or start_prime < 0:
        raise ValueError("start_prime must be a nonnegative integer")
    if max_primes > 2047 or start_prime + max_primes >= 2048 or max_degree >= 2**16 or threads >= 2**32:
        raise ValueError("reconstruction options exceed FiniteFlow native limits")
    if field.characteristic:
        raise ValueError("FiniteFlow cofactor reconstruction requires characteristic zero")
    if shift is None:
        return None
    result = tuple(shift)
    if (len(result) != len(field.parameters)
            or any(type(value) is not int or not 0 <= value < 2**64 for value in result)):
        raise ValueError("the shift needs one unsigned exact integer per parameter")
    return result


def _native_dimensions(nrows: int, nvars: int) -> None:
    if (type(nrows) is not int or type(nvars) is not int or nrows < 0 or nvars < 0
            or nrows >= 2**32 or nvars >= 2**32):
        raise ValueError("cofactor matrix dimensions exceed FiniteFlow native limits")


def _matrix(presentation: AlgebraPresentation, graph: DerivationGraph, node: int,
            contexts: tuple[Context, ...], budget: _Budget, maximum: int) -> Rows:
    target = graph.nodes[node]
    count = len(target.terms)
    if count > maximum:
        raise VerificationLimit("cofactor reconstruction exceeds its sparse matrix-entry bound")
    for _, index, _ in contexts:
        budget.consume()
        count += len(presentation.relations[index].terms)
        if count > maximum:
            raise VerificationLimit("cofactor reconstruction exceeds its sparse matrix-entry bound")
    words = set(dict(target.terms))
    for left, index, right in contexts:
        for word, _ in presentation.relations[index].terms:
            budget.consume(len(left) + len(word) + len(right) + 1)
            words.add(left + word + right)
    _native_dimensions(len(words), len(contexts))
    indices = {word: i for i, word in enumerate(sorted(words, key=lambda w: (len(w), w), reverse=True))}
    rows: Rows = [[] for _ in words]
    for column, (left, index, right) in enumerate(contexts):
        for word, coefficient in presentation.relations[index].terms:
            budget.consume()
            rows[indices[left + word + right]].append((column, coefficient))
    for word, coefficient in target.terms:
        budget.consume()
        rows[indices[word]].append((len(contexts), coefficient))  # C*x=h, positive affine RHS.
    return rows


def _rational(value: Any) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (str, int, Fraction)):
        raise ValueError("native affine coefficients must be exact rational values")
    return Fraction(value)


def _decode_analytic(ff: Any, output: Any, field: RationalFunctionField,
                     count: int, arithmetic: WorkingField) -> list[FieldElement]:
    if (not isinstance(output, ff.RatFunList) or len(output) != count
            or type(output.nvars()) is not int or output.nvars() != len(field.parameters)):
        raise ValueError("native affine output count/type/parameter arity mismatch")
    values = []
    zero = field.zero
    for index in range(count):
        top, bottom = output.monomials(index)
        for terms in (top, bottom):
            if any(len(powers) != len(field.parameters)
                   or any(type(power) is not int or not 0 <= power < 2**16 for power in powers)
                   for _, powers in terms):
                raise ValueError("native affine exponent arity or value mismatch")
            if len({tuple(powers) for _, powers in terms}) != len(terms):
                raise ValueError("native affine output repeats a monomial")
        value = FieldElement(field, {tuple(e): _rational(c) for c, e in top},
                            {tuple(e): _rational(c) for c, e in bottom})
        values.append(zero if value.is_zero() else arithmetic.to_element(arithmetic.from_element(value)))
    return values


def _native_candidate(field: RationalFunctionField, rows: Rows, nvars: int, budget: _Budget, *,
                      max_primes: int, max_degree: int, threads: int, start_prime: int,
                      shift: tuple[int, ...] | None) -> dict[int, FieldElement]:
    # Encode unique literal coefficients before entering native error handling:
    # unsupported input exponents remain invalid-data errors, not native failures.
    sparse_rows = []
    coefficients: list[FieldElement] = []
    for row in rows:
        sparse_rows.append([j for j, _ in row])
        coefficients.extend(c for _, c in row)
    unique: list[tuple[tuple[tuple[str, tuple[int, ...]], ...], ...]] = []
    indexes = []
    lookup = {}
    if field.parameters:
        for coefficient in coefficients:
            budget.consume()
            key = coefficient.numerator, coefficient.denominator
            if key not in lookup:
                for terms in key:
                    budget.consume(len(terms))
                    if any(power >= 2**16 for powers, _ in terms for power in powers):
                        raise ValueError("FiniteFlow monomial exponents must be below 2**16")
                lookup[key] = len(unique)
                unique.append(tuple(tuple((str(c), powers) for powers, c in terms) for terms in key))
            indexes.append(lookup[key])
        constants: list[str] = []
    else:
        budget.consume(len(coefficients))
        constants = [str(Fraction(c.numerator[0][1]) / Fraction(c.denominator[0][1])) for c in coefficients]
    del coefficients, lookup
    arithmetic = WorkingField(field, cache_size=0)
    try:
        ff = importlib.import_module("fflow")
    except ImportError as error:
        raise ImportError("cofactor reconstruction requires the FiniteFlow bindings; see docs/RECONSTRUCTION.md") from error
    with _NATIVE_LOCK:
        if not hasattr(ff, "LSolveValidateFullSystem"):
            raise ImportError("cofactor reconstruction requires the maintained FiniteFlow fork")
        if shift is not None and not hasattr(ff, "SetReconstructionShift"):
            raise ImportError("the FiniteFlow fork needs SetReconstructionShift for a custom shift")
        available = ff.NAvailablePrimes()
        if type(available) is not int or available < 1:
            raise RuntimeError("FiniteFlow returned an invalid prime-table size")
        if start_prime + max_primes >= available:
            raise ValueError("the prime budget must leave a final unused checking prime")
        if shift is not None and any(value >= ff.PrimeNo(available - 1) for value in shift):
            raise ValueError("shift values must be below the smallest FiniteFlow prime")
        with _native_errors(ff):
            graph, input_node = ff.NewGraphWithInput(len(field.parameters))
            functions = output = None
            try:
                if shift is not None:
                    ff.SetReconstructionShift(graph, list(shift))
                if field.parameters:
                    functions = ff.NewIdxRatFunList(len(field.parameters), unique, indexes)
                    node = ff.AlgAnalyticSparseLSolve(graph, input_node, nvars, sparse_rows, functions)
                else:
                    node = ff.AlgNumericSparseLSolve(graph, nvars, sparse_rows, constants)
                ff.LSolveValidateFullSystem(graph, node, False)
                ff.LSolveOnlyNonHomogeneous(graph, node)
                ff.LSolveSparseOutput(graph, node, False)
                ff.LSolveEqWeight(graph, node, list(range(len(rows))))
                ff.SetOutputNode(graph, node)
                ff.LearnEx(graph, prime_no=start_prime)
                if ff.LSolveIsImpossible(graph, node):
                    raise ValueError("native learned an inconsistent affine proposal")
                if ff.LSolveNEqsNVars(graph, node) != (len(rows), nvars):
                    raise ValueError("native affine system shape mismatch")
                count = ff.LSolveNDepVars(graph, node)
                if type(count) is not int or count != nvars:
                    raise ValueError("native affine variable layout/coverage mismatch")
                dependent = tuple(ff.LSolveDepVars(graph, node)) if count else ()
                if (len(dependent) != nvars
                        or any(type(j) is not int for j in dependent)
                        or set(dependent) != set(range(nvars))
                        or ff.LSolveNIndepVars(graph, node) != 0
                        or ff.GraphNParsOut(graph) != nvars or ff.LSolveOutputIsSparse(graph, node)):
                    raise ValueError("native affine variable layout/coverage mismatch")
                if field.parameters:
                    output = ff.ReconstructFunction(graph, n_threads=threads, max_primes=max_primes,
                                                    max_deg=max_degree, start_mod=start_prime)
                    values = _decode_analytic(ff, output, field, nvars, arithmetic)
                else:
                    output = ff.ReconstructNumeric(graph, n_threads=threads, max_primes=max_primes,
                                                   start_mod=start_prime)
                    if not isinstance(output, list) or len(output) != nvars:
                        raise ValueError("native numeric affine output count/type mismatch")
                    values = [field.scalar(_rational(value)) for value in output]
                return {j: value for j, value in zip(dependent, values) if not value.is_zero()}
            except (TypeError, IndexError, OverflowError) as error:
                raise RuntimeError("FiniteFlow returned malformed affine reconstruction data") from error
            finally:
                try:
                    ff.DeleteGraph(graph)
                finally:
                    # Retain coefficient owners until graph deletion; destruct
                    # native lists while the global session lock is still held.
                    functions = output = None


def reconstruct_derivation(presentation: AlgebraPresentation, graph: DerivationGraph, node: int, *,
                           max_operations: int = 10_000_000, max_contexts: int = 100_000,
                           max_nodes: int = 100_000, max_matrix_entries: int = 100_000,
                           max_primes: int = 5, max_degree: int = 24, threads: int = 1,
                           start_prime: int = 0, shift: Iterable[int] | None = None) -> CofactorReconstruction:
    """Reconstruct and check one explicit cofactor against the original inputs.

    Full raw graph validation precedes context collection, including zero and
    unrelated references. The graph proposes support only; its identities need
    not be trusted. Only nonzero paths contribute structural contexts. A failed
    native proposal establishes neither membership nor nonmembership.

    max_contexts bounds the sum of live context-set entries, including a new
    output before its children are released. max_matrix_entries bounds sparse
    coefficient entries including the affine RHS, checked before matrix/word
    allocation. Construction and the final independent exact check separately
    use max_operations. Native prime/degree controls do not bound time, RAM, or
    coefficient growth. This optional backend requires characteristic zero.
    """
    _positive_limit(max_contexts, "max_contexts")
    _positive_limit(max_matrix_entries, "max_matrix_entries")
    budget = _Budget(max_operations)
    order = _validate_graph(presentation, graph, budget, max_nodes)
    if type(node) is not int or not 0 <= node < len(graph.nodes):
        raise ValueError("the selected derivation node is out of range")
    field = presentation.algebra.field
    shift_values = _options(field, max_primes, max_degree, threads, start_prime, shift)
    target = graph.nodes[node]
    stats = {"graph_nodes": len(graph.nodes), "needed_nodes": 0, "expanded_nodes": 0,
             "contexts": 0, "candidate_contexts": 0, "peak_contexts": 0, "structural_operations": 0,
             "matrix_rows": 0, "matrix_entries": 0, "output_slots": 0}
    if target.is_zero():
        identity = CofactorIdentity(presentation.content_hash(), target, ())
    else:
        contexts, context_stats = _collect_contexts(presentation, graph, node, order, budget, max_contexts)
        context_stats["candidate_contexts"] = context_stats.pop("contexts")
        stats.update(context_stats)
        # Prefer short contexts and sparse defining polynomials, independently
        # of any expected certificate or specialized coefficient values.
        budget.consume(len(contexts))
        contexts = tuple(sorted(contexts, key=lambda c: (len(c[0]) + len(c[2]),
                            len(presentation.relations[c[1]].terms), c[1], c[0], c[2])))
        if not contexts:
            raise RuntimeError("the structural context proposal is empty for this nonzero target")
        rows = _matrix(presentation, graph, node, contexts, budget, max_matrix_entries)
        stats.update(matrix_rows=len(rows), matrix_entries=sum(map(len, rows)), output_slots=len(contexts))
        coefficients = _native_candidate(field, rows, len(contexts), budget, max_primes=max_primes,
                                         max_degree=max_degree, threads=threads, start_prime=start_prime,
                                         shift=shift_values)
        del rows
        algebra = presentation.algebra
        identity = CofactorIdentity(presentation.content_hash(), target,
            tuple((algebra.element({contexts[j][0]: coefficient}), contexts[j][1],
                   algebra.element({contexts[j][2]: 1}))
                  for j, coefficient in sorted(coefficients.items(),
                      key=lambda item: (len(contexts[item[0]][0]) + len(contexts[item[0]][2]), contexts[item[0]]))))
    # Import locally: general's optional backend dispatch also imports this module.
    from .general import verify_cofactor
    if not verify_cofactor(presentation, identity, max_operations=max_operations):
        raise RuntimeError("FiniteFlow candidate fails its complete exact original-input identity")
    stats.update(contexts=len(identity.terms), cofactor_terms=len(identity.terms), operations=budget.used)
    return CofactorReconstruction(identity, tuple(sorted(stats.items())))
