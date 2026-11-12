"""Normalized coefficients for repeated exact verification arithmetic.

Public ``FieldElement`` values retain literal numerator/denominator data. This
adapter leaves those data intact and uses a separate working representation:
``Fraction`` over Q, existing scalars over a prime field, and SymPy's cancelling
sparse fraction fields when parameters are present. See SymPy's domain guide,
https://docs.sympy.org/latest/modules/polys/domainsintro.html#rational-function-fields.
"""
from __future__ import annotations

from collections import OrderedDict
from fractions import Fraction
from threading import Lock
from typing import Any

from ..fields import FieldElement, RationalFunctionField, Terms

__all__ = ["WorkingField"]

_CACHE_LOCK = Lock()


class WorkingField:
    """Internal exact arithmetic adapter for one declared coefficient field.

    Native values stay inside a single computation and are not certificate
    encodings. Conversion back can remove denominator factors. A converted
    coefficient expression can therefore be defined where the literal input
    or its derivation is undefined. This does not justify specialization of a
    quotient-algebra statement; the proof may invert additional factors.
    At most ``cache_size`` parameter conversion entries are retained, with least
    recently used entries evicted first. A size of zero disables this cache.
    """

    def __init__(self, field: RationalFunctionField, *, cache_size: int = 16384) -> None:
        if not isinstance(field, RationalFunctionField):
            raise TypeError("working arithmetic requires a declared coefficient field")
        if type(cache_size) is not int or cache_size < 0:
            raise ValueError("working coefficient cache size must be a nonnegative integer")
        self.field = field
        self._native: Any = None
        self._cache_size = cache_size
        self._cache: OrderedDict[tuple[Terms, Terms], Any] = OrderedDict()
        if field.parameters:
            try:
                from sympy import GF, QQ, Symbol, field as fraction_field
            except ImportError as error:
                raise ImportError(
                    "parameter-field verification requires the symbolic extra; "
                    "install 'zipper-algebra[symbolic]'") from error
            ground = GF(field.characteristic) if field.characteristic else QQ
            self._native = fraction_field([Symbol(name) for name in field.parameters], ground)[0]

    @property
    def zero(self):
        if self._native is not None:
            return self._native.zero
        return self.field.zero if self.field.characteristic else Fraction(0)

    @property
    def one(self):
        if self._native is not None:
            return self._native.one
        return self.field.one if self.field.characteristic else Fraction(1)

    def from_element(self, value: FieldElement):
        """Convert an exact value from this field without modifying its data."""
        if not isinstance(value, FieldElement) or value.field != self.field:
            raise TypeError("coefficient belongs to a different working field")
        if self._native is None:
            if self.field.characteristic:
                return value
            # FieldElement normalizes constant denominators to one at construction.
            return Fraction(value.numerator[0][1]) if value.numerator else Fraction(0)
        key = (value.numerator, value.denominator)
        with _CACHE_LOCK:
            if key in self._cache:
                self._cache.move_to_end(key)
                return self._cache[key]
        ring = self._native.ring
        ground = self._native.domain
        top = ring.from_dict({monomial: ground.convert(coefficient)
                              for monomial, coefficient in value.numerator})
        bottom = ring.from_dict({monomial: ground.convert(coefficient)
                                 for monomial, coefficient in value.denominator})
        native = self._native.new(top, bottom)
        if self._cache_size:
            with _CACHE_LOCK:
                self._cache[key] = native
                if len(self._cache) > self._cache_size:
                    self._cache.popitem(last=False)
        return native

    def to_element(self, value) -> FieldElement:
        """Convert a native working value to exact public coefficient data."""
        if isinstance(value, (int, Fraction)):
            return self.field.scalar(value)
        if self._native is None:
            if self.field.characteristic and isinstance(value, FieldElement):
                return self.field.scalar(value)
            raise TypeError("expected an exact scalar from this working field")
        if type(value) is not type(self._native.zero) or value.field != self._native:
            raise TypeError("expected an exact scalar from this working field")

        def coefficients(polynomial):
            if self.field.characteristic:
                return {monomial: int(coefficient) for monomial, coefficient in polynomial.items()}
            ground = self._native.domain
            return {monomial: Fraction(int(ground.numer(coefficient)), int(ground.denom(coefficient)))
                    for monomial, coefficient in polynomial.items()}

        return FieldElement(self.field, coefficients(value.numer), coefficients(value.denom))
