"""Produce exactly checked invariants of certified finite-dimensional algebras."""
from __future__ import annotations

from .matrices import ExactMatrix
from .reconstruction import reconstruct_kernel
from .verification.derivations import VerificationLimit, _positive_limit
from .verification.finite import FiniteElement, VerifiedFiniteAlgebra
from .verification.krylov import MinimalPolynomialCertificate, verify_minimal_polynomial
from .verification.linear import KernelCertificate
from .verification.working_field import WorkingField

__all__ = ["minimal_polynomial"]


def _orbit_matrix(algebra, orbit, degree):
    return ExactMatrix(algebra.field, algebra.dimension, degree + 1,
                       tuple((row, column, value)
                             for column in range(degree + 1)
                             for row, value in enumerate(orbit[column].coordinates)
                             if not value.is_zero()))


def minimal_polynomial(algebra: VerifiedFiniteAlgebra, element: FiniteElement, *,
                       seed: FiniteElement | None = None, max_degree: int | None = None,
                       max_operations=10000000, max_products=10000000, max_primes=5,
                       reconstruction_degree=24, threads=1, start_prime=0,
                       shift=None) -> MinimalPolynomialCertificate:
    """Certify the monic polynomial of least degree with ``p(element)*seed=0``.

    The default seed is the unit, giving the element's minimal polynomial.
    An explicit seed instead gives the minimal annihilator of that orbit;
    it need not annihilate the whole regular representation. The zero seed,
    including the zero algebra, has polynomial 1.

    Exact prefix kernels find the first dependency without a numerical degree
    guess. Prefix lengths grow geometrically, then the dependent prefix is
    truncated to the certified cyclic rank. Production uses at most
    ``min(max_degree, dim A)`` multiplication steps; independent verification
    recomputes the required orbit. An exhausted degree budget is inconclusive.
    ``max_operations`` applies separately to each multiplication in production
    and independent verification, and ``max_products`` to each kernel check.
    They are not wall-time or resident-memory limits.

    Nonzero orbits currently use FiniteFlow over Q or Q(parameters). Trivial
    zero orbits also work in positive characteristic without a native backend.
    Native search options have the meanings of ``reconstruct_kernel``.
    """
    if not isinstance(algebra, VerifiedFiniteAlgebra):
        raise TypeError("minimal polynomials require a verified finite algebra")
    if not isinstance(element, FiniteElement) or element.algebra is not algebra:
        raise TypeError("the element must belong to this verified algebra")
    if seed is None:
        seed = algebra.one
    if not isinstance(seed, FiniteElement) or seed.algebra is not algebra:
        raise TypeError("the seed must belong to this verified algebra")
    if max_degree is not None and (type(max_degree) is not int or max_degree < 0):
        raise ValueError("max_degree must be a nonnegative integer or None")
    _positive_limit(max_operations, "max_operations")
    for name, value, minimum in (("max_products", max_products, 0),
                                 ("max_primes", max_primes, 1),
                                 ("reconstruction_degree", reconstruction_degree, 1),
                                 ("threads", threads, 1), ("start_prime", start_prime, 0)):
        if type(value) is not int or value < minimum:
            raise ValueError(f"{name} must be an integer at least {minimum}")
    if (max_primes > 2047 or start_prime + max_primes >= 2048
            or reconstruction_degree >= 2**16 or threads >= 2**32):
        raise ValueError("reconstruction options exceed FiniteFlow native limits")
    if shift is not None:
        shift = tuple(shift)
        if (len(shift) != len(algebra.field.parameters)
                or any(type(value) is not int or not 0 <= value < 2**64 for value in shift)):
            raise ValueError("the shift needs one unsigned exact integer per parameter")

    field = algebra.field
    orbit = [seed]
    if all(value.is_zero() for value in seed.coordinates):
        matrix = _orbit_matrix(algebra, orbit, 0)
        identity = ExactMatrix(field, 1, 1, ((0, 0, field.one),))
        witness = KernelCertificate(matrix.content_hash(), identity, (0,), (), (),
                                    field.characteristic or 2147483647,
                                    (2,) * len(field.parameters), ())
        coefficients = (field.one,)
    else:
        if field.characteristic:
            raise ValueError("nonzero-orbit production currently requires characteristic zero")
        limit = algebra.dimension if max_degree is None else min(max_degree, algebra.dimension)
        if not limit:
            raise VerificationLimit("minimal-polynomial degree budget exhausted")
        degree = 1
        options = {"max_primes": max_primes, "max_degree": reconstruction_degree,
                   "threads": threads, "max_products": max_products,
                   "start_prime": start_prime, "shift": shift}
        while True:
            while len(orbit) <= degree:
                orbit.append(algebra.multiply(element, orbit[-1], max_operations=max_operations))
            matrix = _orbit_matrix(algebra, orbit, degree)
            witness = reconstruct_kernel(matrix, **options)
            if witness.nullity:
                # If an orbit has its first relation in degree d, invariance
                # under multiplication forces every later vector into the
                # first d-dimensional span. A dependent prefix therefore has
                # rank exactly d; no specialization or pivot guess is needed.
                if witness.rank != degree:
                    degree = witness.rank
                    matrix = _orbit_matrix(algebra, orbit, degree)
                    witness = reconstruct_kernel(matrix, **options)
                if witness.rank != degree or witness.nullity != 1:
                    raise RuntimeError("reconstruction did not return a unique orbit relation")
                values = {row: value for row, _, value in witness.kernel.entries}
                leading = values.get(degree, field.zero)
                if leading.is_zero():
                    raise RuntimeError("reconstruction did not return a monic orbit relation")
                arithmetic = WorkingField(field)
                scale = arithmetic.from_element(leading)
                coefficients = tuple(arithmetic.to_element(
                    arithmetic.from_element(values.get(row, field.zero)) / scale)
                    for row in range(degree)) + (field.one,)
                break
            if degree == limit:
                raise VerificationLimit("minimal-polynomial degree budget exhausted")
            degree = min(2 * degree, limit)

    certificate = MinimalPolynomialCertificate(
        presentation_hash=algebra.presentation.content_hash(), basis_words=algebra.basis_words,
        element=element.coordinates, seed=seed.coordinates, coefficients=coefficients,
        kernel=witness)
    return verify_minimal_polynomial(algebra, certificate, max_operations=max_operations,
                                     max_products=max_products)
