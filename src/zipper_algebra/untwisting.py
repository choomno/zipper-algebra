"""Exact polynomial input for independently specified untwisting relations.

Coefficients are rational and listed in ascending powers. The data retain the
literal polynomial: accepting an input does not assert monicity, an invertible
constant term, or finite dimensionality of the resulting presentation.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from fractions import Fraction
import re

from .coefficients import NumericRing, SymbolicRing
from .polynomial import Poly

__all__ = ["UntwistingPolynomial", "UntwistingFamily"]

_RATIONAL = re.compile(r"[+-]?[0-9]+(?:/[0-9]+)?")


@dataclass(frozen=True, slots=True, init=False)
class UntwistingPolynomial:
    """A nonzero polynomial over Q, with coefficients in ascending order.

    Trailing zero coefficients are removed. Nonmonic polynomials and zero
    constant terms remain available for raw quotient presentations; the
    hypotheses of a finiteness theorem must be checked separately.
    """

    coefficients: tuple[Fraction, ...]

    def __init__(self, coefficients: Iterable[int | Fraction]) -> None:
        values = tuple(coefficients)
        if any(isinstance(value, bool) or not isinstance(value, (int, Fraction))
               for value in values):
            raise TypeError("polynomial coefficients must be exact integers or Fractions")
        exact = tuple(Fraction(value) for value in values)
        last = len(exact)
        while last and not exact[last - 1]:
            last -= 1
        if not last:
            raise ValueError("an untwisting polynomial must be nonzero")
        object.__setattr__(self, "coefficients", exact[:last])

    @property
    def degree(self) -> int:
        """Degree of the rational input polynomial, before specialization."""
        return len(self.coefficients) - 1

    @property
    def constant_term(self) -> Fraction:
        return self.coefficients[0]

    @property
    def is_monic(self) -> bool:
        return self.coefficients[-1] == 1

    def evaluate(self, operator: Poly, ring: NumericRing | SymbolicRing) -> Poly:
        """Evaluate by Horner's rule without rewriting free words.

        Rational coefficients are coerced by ``ring.scalar``. In positive
        characteristic this requires their denominators to be units. Polynomial
        degree and monicity can change under that specialization.
        """
        if not isinstance(operator, Poly):
            raise TypeError("the polynomial argument must be a free-algebra polynomial")
        if operator.char != ring.char:
            raise TypeError("the polynomial argument and coefficient ring have different characteristics")
        result = Poly.scalar(ring.scalar(self.coefficients[-1]))
        for coefficient in reversed(self.coefficients[:-1]):
            result = result * operator + Poly.scalar(ring.scalar(coefficient))
        return result

    def to_json(self) -> list[int | str]:
        """Encode exact coefficients using integers and rational strings."""
        return [value.numerator if value.denominator == 1 else str(value)
                for value in self.coefficients]

    @classmethod
    def from_json(cls, value: object) -> UntwistingPolynomial:
        """Read a JSON coefficient list; decimal and exponent strings are rejected."""
        if not isinstance(value, list):
            raise ValueError("polynomial data must be a JSON coefficient list")
        coefficients: list[Fraction] = []
        for item in value:
            if type(item) is int:
                coefficients.append(Fraction(item))
            elif isinstance(item, str) and _RATIONAL.fullmatch(item):
                try:
                    coefficients.append(Fraction(item))
                except (ValueError, ZeroDivisionError) as error:
                    raise ValueError("coefficient strings must represent exact rational numbers") from error
            else:
                raise ValueError("JSON coefficients must be integers or rational strings")
        return cls(coefficients)


@dataclass(frozen=True, slots=True, init=False)
class UntwistingFamily:
    """Immutable per-level polynomial input; absent levels impose no relation.

    The enclosing presentation checks that the selected levels fit its strand
    count. Empty families are valid and do not add untwisting relations.
    """

    items: tuple[tuple[int, UntwistingPolynomial], ...]

    def __init__(self, polynomials: Mapping[int, UntwistingPolynomial]) -> None:
        if not isinstance(polynomials, Mapping):
            raise TypeError("an untwisting family must map levels to polynomials")
        for level, polynomial in polynomials.items():
            if type(level) is not int or level < 2:
                raise ValueError("untwisting levels must be integers at least two")
            if not isinstance(polynomial, UntwistingPolynomial):
                raise TypeError("each untwisting level must specify an UntwistingPolynomial")
        object.__setattr__(self, "items", tuple(sorted(polynomials.items())))

    def to_json(self) -> dict[str, list[int | str]]:
        return {str(level): polynomial.to_json() for level, polynomial in self.items}

    @classmethod
    def from_json(cls, value: object) -> UntwistingFamily:
        """Read canonical decimal level names and exact coefficient lists."""
        if not isinstance(value, dict):
            raise ValueError("untwisting-family data must be a JSON object")
        polynomials: dict[int, UntwistingPolynomial] = {}
        for key, coefficients in value.items():
            if not isinstance(key, str) or not key.isdecimal():
                raise ValueError("untwisting level names must be canonical decimal integers")
            level = int(key)
            if level < 2 or str(level) != key:
                raise ValueError("untwisting level names must be canonical integers at least two")
            polynomials[level] = UntwistingPolynomial.from_json(coefficients)
        return cls(polynomials)
