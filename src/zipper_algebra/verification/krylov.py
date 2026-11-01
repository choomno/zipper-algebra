"""Exact certificates for minimal polynomials on cyclic algebra orbits.

For an element a and a seed v in a verified finite algebra, a certificate
proves the monic polynomial of least degree m with m(a)v = 0. The seed 1 gives
the minimal polynomial of a; another seed concerns only that cyclic orbit.
No reconstruction library is loaded by this verifier.
"""
from __future__ import annotations

from dataclasses import dataclass
import re

from ..fields import FieldElement
from ..matrices import ExactMatrix
from .derivations import _positive_limit
from .finite import FiniteElement, VerifiedFiniteAlgebra
from .linear import KernelCertificate, verify_kernel
from .working_field import WorkingField


@dataclass(frozen=True, slots=True)
class MinimalPolynomialCertificate:
    """A monic orbit relation, with an exact complete Krylov-kernel witness.

    Coordinates refer to the ordered named basis of the supplied verified
    algebra. JSON stores the element and seed sparsely. Coefficients are in
    ascending order, so coefficients[j] multiplies the j-th power.
    """

    # Coefficient equality is algebraic; literal encodings are not equality hashes.
    __hash__ = None  # type: ignore[assignment]

    presentation_hash: str
    basis_words: tuple[tuple[int, ...], ...]
    element: tuple[FieldElement, ...]
    seed: tuple[FieldElement, ...]
    coefficients: tuple[FieldElement, ...]
    kernel: KernelCertificate

    def __post_init__(self):
        if (not isinstance(self.presentation_hash, str)
                or re.fullmatch(r"[0-9a-f]{64}", self.presentation_hash) is None):
            raise ValueError("a minimal-polynomial certificate needs the full presentation hash")
        if (not isinstance(self.basis_words, tuple)
                or any(not isinstance(word, tuple)
                       or any(type(letter) is not int or letter < 1 for letter in word)
                       for word in self.basis_words)
                or len(set(self.basis_words)) != len(self.basis_words)):
            raise ValueError("the named basis must contain distinct words of positive indices")
        if not isinstance(self.kernel, KernelCertificate):
            raise TypeError("a minimal-polynomial certificate needs a kernel certificate")
        dimension = len(self.basis_words)
        for values in (self.element, self.seed):
            if not isinstance(values, tuple) or len(values) != dimension:
                raise ValueError("element and seed coordinates must match the named basis")
        if (not isinstance(self.coefficients, tuple)
                or not 1 <= len(self.coefficients) <= dimension + 1):
            raise ValueError("the polynomial degree must be between zero and the algebra dimension")
        field = self.kernel.kernel.field
        if any(not isinstance(value, FieldElement) or value.field != field
               for values in (self.element, self.seed, self.coefficients) for value in values):
            raise TypeError("all certificate coefficients must belong to the declared field")

    @property
    def degree(self):
        return len(self.coefficients) - 1

    def to_json(self):
        def sparse(values):
            return [[row, value.to_json()] for row, value in enumerate(values)
                    if not value.is_zero()]
        return {"schema": "zipper-minimal-polynomial-certificate-v1",
                "presentation_hash": self.presentation_hash,
                "basis_words": [list(word) for word in self.basis_words],
                "element": sparse(self.element), "seed": sparse(self.seed),
                "coefficients": [value.to_json() for value in self.coefficients],
                "kernel": self.kernel.to_json()}

    @classmethod
    def from_json(cls, value):
        keys = {"schema", "presentation_hash", "basis_words", "element", "seed",
                "coefficients", "kernel"}
        if (not isinstance(value, dict) or set(value) != keys
                or value["schema"] != "zipper-minimal-polynomial-certificate-v1"
                or any(not isinstance(value[key], list)
                       for key in ("basis_words", "element", "seed", "coefficients"))
                or any(not isinstance(word, list) for word in value["basis_words"])):
            raise ValueError("invalid minimal-polynomial certificate schema")
        kernel = KernelCertificate.from_json(value["kernel"])
        words = tuple(tuple(word) for word in value["basis_words"])
        field = kernel.kernel.field

        def coordinates(entries):
            values = [field.zero] * len(words)
            seen = set()
            for entry in entries:
                if (not isinstance(entry, list) or len(entry) != 2
                        or type(entry[0]) is not int or not 0 <= entry[0] < len(words)
                        or entry[0] in seen):
                    raise ValueError("invalid sparse minimal-polynomial coordinates")
                row, coefficient = entry
                seen.add(row)
                values[row] = FieldElement.from_json(coefficient)
            return tuple(values)

        return cls(value["presentation_hash"], words, coordinates(value["element"]),
                   coordinates(value["seed"]),
                   tuple(FieldElement.from_json(c) for c in value["coefficients"]), kernel)


def krylov_matrix(algebra: VerifiedFiniteAlgebra, element: FiniteElement,
                  seed: FiniteElement, degree: int, *, max_operations=10000000) -> ExactMatrix:
    """Return the columns v, av, ..., a**degree v, in the certified basis.

    The operation limit applies to each multiplication, as in algebra.multiply.
    Degree is bounded by the verified dimension. No full regular operator or
    multiplication table is materialized.
    """
    if not isinstance(algebra, VerifiedFiniteAlgebra):
        raise TypeError("Krylov computation requires a verified finite algebra")
    if any(not isinstance(value, FiniteElement) or value.algebra is not algebra
           for value in (element, seed)):
        raise TypeError("the element and seed must belong to the supplied verified algebra")
    if type(degree) is not int or not 0 <= degree <= algebra.dimension:
        raise ValueError("Krylov degree must be between zero and the verified dimension")
    _positive_limit(max_operations, "max_operations")
    entries: list[tuple[int, int, FieldElement]] = []
    vector = seed
    for column in range(degree + 1):
        entries.extend((row, column, value) for row, value in enumerate(vector.coordinates)
                       if not value.is_zero())
        if column < degree:
            vector = algebra.multiply(element, vector, max_operations=max_operations)
    return ExactMatrix(algebra.field, algebra.dimension, degree + 1, tuple(entries))


def verify_minimal_polynomial(algebra: VerifiedFiniteAlgebra,
                              certificate: MinimalPolynomialCertificate, *,
                              max_operations=10000000, max_products=10000000):
    """Check generic annihilation and minimality; return the supplied certificate.

    Recompute the cyclic orbit from the verified algebra, then check the complete
    kernel and its rank. A rank-d prefix with d+1 columns and a relation whose
    last coefficient is one proves independence of the first d powers. This
    rules out every lower-degree annihilator. The zero seed has polynomial 1.
    max_operations is per multiplication; max_products bounds the kernel check.
    These are counted-work limits, not wall-time or memory limits.
    """
    if (not isinstance(algebra, VerifiedFiniteAlgebra)
            or not isinstance(certificate, MinimalPolynomialCertificate)):
        raise TypeError("minimal-polynomial verification requires a verified algebra and certificate")
    _positive_limit(max_operations, "max_operations")
    if type(max_products) is not int or max_products < 0:
        raise ValueError("the product budget must be a nonnegative integer")
    if (certificate.presentation_hash != algebra.presentation.content_hash()
            or certificate.basis_words != algebra.basis_words):
        raise ValueError("the minimal-polynomial certificate is bound to a different algebra or basis")
    if certificate.kernel.kernel.field != algebra.field:
        raise TypeError("the certificate and verified algebra have different coefficient fields")
    degree = certificate.degree
    if (certificate.kernel.rank != degree or certificate.kernel.nullity != 1
            or certificate.coefficients[-1] != algebra.field.one):
        raise ValueError("a minimal-polynomial certificate needs a monic relation and rank-degree witness")
    matrix = krylov_matrix(algebra, algebra.element(certificate.element),
                           algebra.element(certificate.seed), degree,
                           max_operations=max_operations)
    verify_kernel(matrix, certificate.kernel, max_products=max_products)
    working = WorkingField(algebra.field)
    vector = {row: working.from_element(value)
              for row, _column, value in certificate.kernel.kernel.entries}
    scale = vector.get(degree, working.zero)
    for row, coefficient in enumerate(certificate.coefficients):
        if working.from_element(coefficient) * scale - vector.get(row, working.zero) != 0:
            raise ValueError("the polynomial coefficients differ from the exact Krylov relation")
    return certificate
