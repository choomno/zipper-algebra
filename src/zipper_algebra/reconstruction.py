"""FiniteFlow reconstruction of exact right kernels over rational function fields.

FiniteFlow learns and reconstructs a candidate. ``reconstruct_kernel`` checks
its full generic identity and a lower rank witness before returning a certificate.
``reconstruct_kernel_candidate`` returns unverified vectors for producers whose
final mathematical certificate is checked independently.
The external FiniteFlow Python bindings are optional and loaded only on use.
"""
from __future__ import annotations

from contextlib import contextmanager
from fractions import Fraction
import importlib
import random
import threading

from .fields import FieldElement
from .matrices import ExactMatrix
from .verification.linear import KernelCertificate, _evaluate, verify_kernel

_PRIMES = (2147483647, 2147483629, 2147483587)
# FiniteFlow graphs share a native global Session. Internal native probe threads
# are supported, but separate Python callers must not mutate that Session at once.
_NATIVE_LOCK = threading.Lock()


@contextmanager
def _native_errors(ff):
    """Normalize native status and decoded-coefficient failures."""
    failures = tuple(getattr(ff, name) for name in
                     ("FFlowError", "Failed", "MissingPoints", "MissingPrimes")
                     if hasattr(ff, name))
    try:
        yield
    except failures + (ValueError, ZeroDivisionError) as error:
        raise RuntimeError(f"FiniteFlow reconstruction failed: {error}") from error


def _lower_inverse(matrix, rows, columns):
    """Produce a small modular witness using SymPy's exact matrix routines."""
    from sympy import GF
    from sympy.polys.matrices import DomainMatrix
    from sympy.polys.matrices.exceptions import DMNonInvertibleMatrixError
    random_points = random.Random(0)
    for prime in _PRIMES:
        for _ in range(8 if matrix.field.parameters else 1):
            point = tuple(random_points.randrange(prime) for _ in matrix.field.parameters)
            try:
                values = {(i, j): _evaluate(c, point, prime) for i, j, c in matrix.entries}
            except ValueError:
                continue
            if not columns:
                return prime, point, ()
            minor = [[values.get((i, j), 0) for j in columns] for i in rows]
            try:
                inverse = DomainMatrix.from_list(minor, GF(prime)).inv().to_list()
            except DMNonInvertibleMatrixError:
                continue
            return prime, point, tuple(tuple(int(x) % prime for x in row) for row in inverse)
    raise RuntimeError("no denominator-safe nonzero minor found within the bounded witness search")



def _reconstruct_kernel_data(matrix: ExactMatrix, *, max_primes=5, max_degree=24,
                             threads=1, max_products=10000000, start_prime=0, shift=None):
    """Run the shared native producer; returned vectors and indices are proposals."""
    if not isinstance(matrix, ExactMatrix):
        raise TypeError("kernel reconstruction requires an exact matrix")
    if matrix.field.characteristic:
        raise ValueError("FiniteFlow reconstruction currently requires characteristic zero")
    for name, value, minimum in (("max_primes", max_primes, 1),
                                 ("max_degree", max_degree, 1), ("threads", threads, 1),
                                 ("max_products", max_products, 0), ("start_prime", start_prime, 0)):
        if type(value) is not int or value < minimum:
            raise ValueError(f"{name} must be an integer at least {minimum}")
    if max_primes > 2047 or start_prime + max_primes >= 2048 or max_degree >= 2**16 or threads >= 2**32:
        raise ValueError("reconstruction options exceed FiniteFlow native limits")
    if matrix.nrows >= 2**32 or matrix.ncols >= 2**32:
        raise ValueError("matrix dimensions exceed FiniteFlow native limits")
    field = matrix.field
    if shift is not None:
        shift = tuple(shift)
        if (len(shift) != len(field.parameters)
                or any(type(x) is not int or not 0 <= x < 2**64 for x in shift)):
            raise ValueError("the shift needs one unsigned exact integer per parameter")
    # Empty and zero systems have a complete immediate proof and need no native
    # graph. There are no literal denominator factors hidden in canonical zero.
    if not matrix.entries:
        kernel = ExactMatrix(field, matrix.ncols, matrix.ncols,
                             tuple((i, i, field.one) for i in range(matrix.ncols)))
        return kernel, tuple(range(matrix.ncols)), (), ()
    if any(power >= 2**16 for _, _, c in matrix.entries
           for terms in (c.numerator, c.denominator) for exponents, _ in terms for power in exponents):
        raise ValueError("FiniteFlow monomial exponents must be below 2**16")
    try:
        ff = importlib.import_module("fflow")
    except ImportError as error:
        raise ImportError("kernel reconstruction requires the FiniteFlow Python bindings; "
                          "see docs/RECONSTRUCTION.md") from error
    if not hasattr(ff, "LSolveValidateFullSystem"):
        raise ImportError("kernel reconstruction requires the zipper-algebra FiniteFlow fork "
                          "with LSolveValidateFullSystem; see docs/RECONSTRUCTION.md")
    if start_prime + max_primes >= ff.NAvailablePrimes():
        raise ValueError("the prime budget must leave a final unused checking prime")
    if shift is not None:
        if not hasattr(ff, "SetReconstructionShift"):
            raise ImportError("the FiniteFlow fork needs SetReconstructionShift for a custom shift")
        if any(x >= ff.PrimeNo(ff.NAvailablePrimes() - 1) for x in shift):
            raise ValueError("shift values must be below the smallest FiniteFlow prime")
    with _NATIVE_LOCK, _native_errors(ff):
        graph, input_node = ff.NewGraphWithInput(len(field.parameters))
        try:
            if shift is not None:
                ff.SetReconstructionShift(graph, list(shift))
            coefficients = []
            sparse_rows: list[list[int]] = [[] for _ in range(matrix.nrows)]
            for i, j, coefficient in matrix.entries:
                sparse_rows[i].append(j)
                coefficients.append(tuple(tuple((str(c), powers) for powers, c in terms)
                                          for terms in (coefficient.numerator, coefficient.denominator)))
            if field.parameters:
                functions = ff.NewIdxRatFunList(len(field.parameters), coefficients)
                node = ff.AlgAnalyticSparseLSolve(graph, input_node, matrix.ncols, sparse_rows, functions)
            else:
                numeric_coefficients = [str(Fraction(c.numerator[0][1]) / Fraction(c.denominator[0][1]))
                                        for _, _, c in matrix.entries]
                node = ff.AlgNumericSparseLSolve(graph, matrix.ncols, sparse_rows, numeric_coefficients)
            ff.LSolveValidateFullSystem(graph, node, True)
            ff.LSolveOnlyHomogeneous(graph, node)
            ff.LSolveSparseOutput(graph, node, False)
            ff.SetOutputNode(graph, node)
            ff.LearnEx(graph, prime_no=start_prime)
            dependent = tuple(ff.LSolveDepVars(graph, node)) if ff.LSolveNDepVars(graph, node) else ()
            free = tuple(ff.LSolveIndepVars(graph, node)) if ff.LSolveNIndepVars(graph, node) else ()
            rows = tuple(ff.LSolveIndepEqs(graph, node)) if ff.LSolveNIndepEqs(graph, node) else ()
            if (len(rows) != len(dependent) or len(dependent) + len(free) != matrix.ncols
                    or any(type(i) is not int for i in dependent + free + rows)
                    or set(dependent + free) != set(range(matrix.ncols))
                    or len(set(rows)) != len(rows) or any(not 0 <= i < matrix.nrows for i in rows)):
                raise RuntimeError("FiniteFlow returned inconsistent learned matrix indices")
            entries = [(i, j, field.one) for j, i in enumerate(free)]
            if dependent and free:
                values = []
                if field.parameters:
                    output = ff.ReconstructFunction(graph, n_threads=threads,
                                                    max_primes=max_primes, max_deg=max_degree, start_mod=start_prime)
                    if (not isinstance(output, ff.RatFunList) or len(output) != len(dependent) * len(free)
                            or output.nvars() != len(field.parameters)):
                        raise RuntimeError("FiniteFlow did not reconstruct the complete kernel coefficients")
                    for index in range(len(output)):
                        top, bottom = output.monomials(index)
                        # Validate exponent arity/exact scalars and never silently
                        # merge repeated native monomials at the public boundary.
                        if any(len({tuple(e) for _, e in terms}) != len(terms) for terms in (top, bottom)):
                            raise RuntimeError("FiniteFlow returned repeated monomials")
                        values.append(FieldElement(field, {tuple(e): Fraction(c) for c, e in top},
                                                   {tuple(e): Fraction(c) for c, e in bottom}))
                else:
                    numeric = ff.ReconstructNumeric(graph, n_threads=threads, max_primes=max_primes,
                                                    start_mod=start_prime)
                    if not isinstance(numeric, list) or len(numeric) != len(dependent) * len(free):
                        raise RuntimeError("FiniteFlow did not reconstruct the complete numeric kernel")
                    values = [field.scalar(Fraction(value)) for value in numeric]
                entries.extend((dependent[index // len(free)], index % len(free), value)
                               for index, value in enumerate(values))
            kernel = ExactMatrix(field, matrix.ncols, len(free), tuple(entries))
        finally:
            ff.DeleteGraph(graph)
    return kernel, free, rows, dependent


def reconstruct_kernel_candidate(matrix: ExactMatrix, *, max_primes=5, max_degree=24,
                                 threads=1, start_prime=0, shift=None) -> ExactMatrix:
    """Return unverified right-kernel vectors proposed by FiniteFlow.

    No rank, completeness or exact-kernel assertion is established here. Native
    full-system probe checks and wire validation remain enabled, but candidates
    can still be wrong or incomplete. A downstream checker must establish the
    actual identity or algebraic statement before accepting a result. In
    particular, a returned zero-column matrix does not prove full rank.

    This route skips modular-minor construction and independent kernel checking.
    It is intended for ideal-identity producers whose final exact certificate
    checks all local identities. Use ``reconstruct_kernel`` for a complete-kernel
    certificate, as needed by rank and minimal-polynomial computations.

    Native prime, degree, thread and shift controls match the certified route.
    There is no independent verification or product-check budget on this path.
    Neither route imposes process RAM or wall-clock limits.
    """
    kernel, _, _, _ = _reconstruct_kernel_data(matrix, max_primes=max_primes,
                    max_degree=max_degree, threads=threads, start_prime=start_prime, shift=shift)
    return kernel


def reconstruct_kernel(matrix: ExactMatrix, *, max_primes=5, max_degree=24,
                       threads=1, max_products=10000000, start_prime=0, shift=None) -> KernelCertificate:
    """Reconstruct and certify the entire right kernel of ``matrix``.

    The producer currently supports Q and Q(parameters). Parameters and entries
    are passed as exact monomial arrays, without parsing strings as expressions.
    Parameter-field arithmetic requires the symbolic extra, including zero
    systems; nonzero ground systems also need it to produce the minor witness.
    Ground-field certificate checks with small primes are dependency-free.
    ``max_primes`` (1..2047) and ``max_degree`` (1..65535) bound searches; they are
    not wall-clock or memory limits. No partial or merely sampled answer is
    returned when reconstruction or independent verification fails. ``start_prime``
    selects the first FiniteFlow prime-table index for learning and reconstruction;
    ``shift`` optionally supplies an exact coordinate vector for degree learning.
    """
    kernel, free, rows, dependent = _reconstruct_kernel_data(matrix, max_primes=max_primes,
                    max_degree=max_degree, threads=threads, max_products=max_products,
                    start_prime=start_prime, shift=shift)
    if not matrix.entries:
        certificate = KernelCertificate(matrix.content_hash(), kernel, free,
                                        (), (), _PRIMES[0], (2,) * len(matrix.field.parameters), ())
        return verify_kernel(matrix, certificate, max_products=max_products)
    prime, point, inverse = _lower_inverse(matrix, rows, dependent)
    certificate = KernelCertificate(matrix.content_hash(), kernel, free, rows, dependent,
                                    prime, point, inverse)
    try:
        return verify_kernel(matrix, certificate, max_products=max_products)
    except (ValueError, ZeroDivisionError) as error:
        raise RuntimeError(f"FiniteFlow candidate failed independent exact verification: {error}") from error
