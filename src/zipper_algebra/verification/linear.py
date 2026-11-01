"""Independent exact verification of reconstructed right kernels.

A free identity block proves independence of the proposed kernel vectors. A
nonzero specialized minor gives the opposite rank bound; together with A*N=0
this proves that the proposed vectors form the entire kernel over the declared
function field. Specialization is used only for the lower bound.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import re
from typing import Any

from ..matrices import ExactMatrix
from ..fields import FieldElement
from ..linear import validate_field_characteristic
from .working_field import WorkingField
from .derivations import VerificationLimit


def _prime(prime):
    # Bound existing trial division to small moduli. Large FF-sized primes use
    # SymPy's definitive (<2**64) test without constructing a coefficient field.
    # https://docs.sympy.org/latest/modules/ntheory.html#sympy.ntheory.primetest.isprime
    if type(prime) is not int or not 2 <= prime < 2**64:
        raise ValueError("the minor witness needs a prime below 2**64")
    if prime < 2**32:
        validate_field_characteristic(prime)
        return
    try:
        from sympy import isprime
    except ImportError as error:
        raise ImportError("large witness primes require the symbolic extra; "
                          "install 'zipper-algebra[symbolic]'") from error
    if not isprime(prime):
        raise ValueError("the minor witness needs a prime below 2**64")


def _indices(values, size, count):
    if (not isinstance(values, tuple) or len(values) != count
            or any(type(i) is not int or not 0 <= i < size for i in values)
            or len(set(values)) != count):
        raise ValueError("witness indices must be distinct and inside the matrix")


def _rational_mod(value, prime):
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("specialization values must be exact integers or Fractions")
    value = Fraction(value)
    try:
        return value.numerator * pow(value.denominator, -1, prime) % prime
    except ValueError:
        raise ValueError("a scalar denominator vanishes at the witness prime") from None


def _evaluate(value: FieldElement, point, prime):
    def polynomial(terms):
        total = 0
        for powers, coefficient in terms:
            term = _rational_mod(coefficient, prime)
            for x, power in zip(point, powers):
                term = term * pow(x, power, prime) % prime
            total = (total + term) % prime
        return total
    numerator, denominator = polynomial(value.numerator), polynomial(value.denominator)
    if not denominator:
        raise ValueError("a matrix denominator vanishes at the witness point")
    return numerator * pow(denominator, -1, prime) % prime


@dataclass(frozen=True, slots=True)
class KernelCertificate:
    """A complete kernel basis and a denominator-safe modular minor witness."""

    # Coefficient equality is algebraic; literal encodings are not equality hashes.
    __hash__ = None  # type: ignore[assignment]

    matrix_hash: str
    kernel: ExactMatrix
    identity_rows: tuple[int, ...]
    pivot_rows: tuple[int, ...]
    pivot_columns: tuple[int, ...]
    prime: int
    point: tuple[int | Fraction, ...]
    minor_inverse: tuple[tuple[int, ...], ...]

    def __post_init__(self):
        if (not isinstance(self.matrix_hash, str) or len(self.matrix_hash) != 64
                or any(c not in "0123456789abcdef" for c in self.matrix_hash)):
            raise ValueError("a kernel certificate needs the complete matrix SHA-256")
        if not isinstance(self.kernel, ExactMatrix):
            raise TypeError("a kernel certificate requires an exact matrix")
        for name in ("identity_rows", "pivot_rows", "pivot_columns", "point"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        object.__setattr__(self, "minor_inverse", tuple(tuple(row) for row in self.minor_inverse))

    @property
    def rank(self):
        return len(self.pivot_columns)

    @property
    def nullity(self):
        return self.kernel.ncols

    def to_json(self):
        return {"schema": "zipper-kernel-certificate-v1", "matrix_hash": self.matrix_hash,
                "kernel": self.kernel.to_json(), "identity_rows": list(self.identity_rows),
                "pivot_rows": list(self.pivot_rows), "pivot_columns": list(self.pivot_columns),
                "prime": self.prime, "point": [str(x) for x in self.point],
                "minor_inverse": [list(row) for row in self.minor_inverse]}

    @classmethod
    def from_json(cls, value):
        keys = {"schema", "matrix_hash", "kernel", "identity_rows", "pivot_rows",
                "pivot_columns", "prime", "point", "minor_inverse"}
        if (not isinstance(value, dict) or set(value) != keys
                or value["schema"] != "zipper-kernel-certificate-v1"
                or any(not isinstance(value[key], list) for key in
                       ("identity_rows", "pivot_rows", "pivot_columns", "point", "minor_inverse"))
                or any(not isinstance(x, str) or not re.fullmatch(r"[+-]?[0-9]+(?:/[0-9]+)?", x)
                       for x in value["point"])
                or any(not isinstance(row, list) for row in value["minor_inverse"])):
            raise ValueError("invalid kernel certificate schema")
        try:
            point = tuple(Fraction(x) for x in value["point"])
        except ZeroDivisionError as error:
            raise ValueError("the minor witness point contains a zero scalar denominator") from error
        return cls(value["matrix_hash"], ExactMatrix.from_json(value["kernel"]),
                   tuple(value["identity_rows"]), tuple(value["pivot_rows"]),
                   tuple(value["pivot_columns"]), value["prime"],
                   point,
                   tuple(tuple(row) for row in value["minor_inverse"]))


def verify_kernel(matrix: ExactMatrix, certificate: KernelCertificate, *, max_products=10000000):
    """Return the certificate after proving both rank bounds and A*N=0.

    ``max_products`` bounds scalar products, not elapsed time or memory. The
    verifier does not invoke FiniteFlow or trust its learning/reconstruction.
    Parameter-field arithmetic and witness primes >=2**32 require the symbolic
    extra; missing optional arithmetic raises ImportError.
    """
    if not isinstance(matrix, ExactMatrix) or not isinstance(certificate, KernelCertificate):
        raise TypeError("kernel verification requires an exact matrix and certificate")
    if type(max_products) is not int or max_products < 0:
        raise ValueError("the product budget must be a nonnegative integer")
    kernel, rank = certificate.kernel, certificate.rank
    if certificate.matrix_hash != matrix.content_hash():
        raise ValueError("the certificate is bound to a different matrix")
    if (kernel.field != matrix.field or kernel.nrows != matrix.ncols
            or rank > min(matrix.nrows, matrix.ncols) or kernel.ncols != matrix.ncols - rank):
        raise ValueError("kernel field or dimensions disagree with the rank witness")
    _indices(certificate.identity_rows, matrix.ncols, kernel.ncols)
    _indices(certificate.pivot_rows, matrix.nrows, rank)
    _indices(certificate.pivot_columns, matrix.ncols, rank)
    prime = certificate.prime
    if matrix.field.characteristic not in (0, prime):
        raise ValueError("the minor witness cannot change the source characteristic")
    _prime(prime)
    if len(certificate.point) != len(matrix.field.parameters):
        raise ValueError("the witness needs one assignment per parameter")
    point = tuple(_rational_mod(x, prime) for x in certificate.point)
    inverse = certificate.minor_inverse
    if (len(inverse) != rank or any(len(row) != rank for row in inverse)
            or any(type(x) is not int or not 0 <= x < prime for row in inverse for x in row)):
        raise ValueError("the modular minor inverse has invalid dimensions or entries")
    used = 0

    def tick():
        nonlocal used
        used += 1
        if used > max_products:
            raise VerificationLimit("kernel verification exceeded its scalar-product budget")

    # Evaluate every stored coefficient, including those outside the selected
    # minor. This certifies a well-defined specialization of the input matrix.
    specialized = {(i, j): _evaluate(c, point, prime) for i, j, c in matrix.entries}
    minor = [[specialized.get((i, j), 0) for j in certificate.pivot_columns]
             for i in certificate.pivot_rows]
    for i in range(rank):
        for j in range(rank):
            total = 0
            for k in range(rank):
                tick()
                total = (total + inverse[i][k] * minor[k][j]) % prime
            if total != int(i == j):
                raise ValueError("the selected modular minor inverse identity fails")
    working = WorkingField(matrix.field)
    native_kernel = {(i, j): working.from_element(c) for i, j, c in kernel.entries}
    for i, row in enumerate(certificate.identity_rows):
        for j in range(kernel.ncols):
            if native_kernel.get((row, j), working.zero) - int(i == j) != 0:
                raise ValueError("the kernel free-coordinate block is not the identity")
    rows: dict[int, list[tuple[int, Any]]] = {}
    for (i, j), c in native_kernel.items():
        rows.setdefault(i, []).append((j, c))
    residual: dict[int, Any] = {}
    current_row = None
    for i, k, coefficient in matrix.entries:
        if current_row is not None and i != current_row:
            if any(value != 0 for value in residual.values()):
                raise ValueError("the proposed kernel fails A*N=0")
            residual.clear()
        current_row = i
        native = working.from_element(coefficient)
        for j, other in rows.get(k, ()):
            tick()
            residual[j] = residual.get(j, working.zero) + native * other
    if any(value != 0 for value in residual.values()):
        raise ValueError("the proposed kernel fails A*N=0")
    return certificate
