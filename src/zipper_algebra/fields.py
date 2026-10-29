"""Exact rational functions over Q or a checked prime field.

Parameter order is part of the coefficient domain. Fractions retain polynomial
denominators, without a multivariate gcd calculation. Equality is checked by
cross multiplication. Specialization evaluates the stored expressions wherever
their denominators survive; it is not a map from an entire function field.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
import re

from .linear import validate_field_characteristic

__all__ = ["RationalFunctionField", "FieldElement"]

Exponent = tuple[int, ...]
Scalar = int | Fraction
Terms = tuple[tuple[Exponent, Scalar], ...]
_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*")
_RATIONAL = re.compile(r"[+-]?[0-9]+(?:/[0-9]+)?")


def _base(value: Scalar, characteristic: int) -> Scalar:
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise TypeError("field scalars must be exact integers or Fractions")
    rational = Fraction(value)
    if not characteristic:
        return rational
    try:
        inverse = pow(rational.denominator, -1, characteristic)
    except ValueError:
        raise ZeroDivisionError("scalar denominator vanishes in the prime field") from None
    return rational.numerator * inverse % characteristic


def _clean(field: RationalFunctionField, coefficients: Mapping[Exponent, Scalar]) -> Terms:
    if not isinstance(coefficients, Mapping):
        raise TypeError("a sparse polynomial must map exponent tuples to exact scalars")
    terms = []
    for exponents, coefficient in coefficients.items():
        if (not isinstance(exponents, tuple) or len(exponents) != len(field.parameters)
                or any(type(power) is not int or power < 0 for power in exponents)):
            raise ValueError("monomials must have one nonnegative integer exponent per parameter")
        value = _base(coefficient, field.characteristic)
        if value:
            terms.append((exponents, value))
    return tuple(sorted(terms))


def _add(field: RationalFunctionField, left: Terms, right: Terms) -> Terms:
    result = dict(left)
    for monomial, coefficient in right:
        result[monomial] = result.get(monomial, 0) + coefficient
    return _clean(field, result)


def _mul(field: RationalFunctionField, left: Terms, right: Terms) -> Terms:
    result: dict[Exponent, Scalar] = {}
    for monomial, coefficient in left:
        for other, scalar in right:
            product = tuple(a + b for a, b in zip(monomial, other))
            result[product] = result.get(product, 0) + coefficient * scalar
    return _clean(field, result)


def _encode_scalar(value: Scalar) -> int | str:
    rational = Fraction(value)
    return rational.numerator if rational.denominator == 1 else str(rational)


def _decode_scalar(value: object) -> Scalar:
    if type(value) is int:
        return value
    if isinstance(value, str) and _RATIONAL.fullmatch(value):
        try:
            return Fraction(value)
        except (ValueError, ZeroDivisionError):
            pass
    raise ValueError("JSON scalars must be integers or exact rational strings")


@dataclass(frozen=True, slots=True)
class RationalFunctionField:
    """Q(parameters) or GF(p)(parameters), with an ordered parameter tuple.

    With no parameters this is Q or GF(p). Positive characteristic is verified
    by exact trial division; very large proposed characteristics can therefore
    require substantial validation time.
    """

    parameters: tuple[str, ...] = ()
    characteristic: int = 0

    def __post_init__(self) -> None:
        validate_field_characteristic(self.characteristic)
        if (not isinstance(self.parameters, tuple)
                or any(not isinstance(name, str) or not _NAME.fullmatch(name)
                       for name in self.parameters)):
            raise ValueError("parameters must be a tuple of ASCII identifiers")
        if len(set(self.parameters)) != len(self.parameters):
            raise ValueError("parameter names must be distinct")

    def scalar(self, value: Scalar | FieldElement) -> FieldElement:
        if isinstance(value, FieldElement):
            if value.field != self:
                raise TypeError("cannot mix coefficient fields")
            return value
        return FieldElement(self, {(0,) * len(self.parameters): value})

    def polynomial(self, coefficients: Mapping[Exponent, Scalar]) -> FieldElement:
        return FieldElement(self, coefficients)

    def parameter(self, name: str) -> FieldElement:
        if name not in self.parameters:
            raise ValueError(f"unknown parameter {name!r}")
        index = self.parameters.index(name)
        return self.polynomial({tuple(int(i == index) for i in range(len(self.parameters))): 1})

    @property
    def zero(self) -> FieldElement:
        return self.polynomial({})

    @property
    def one(self) -> FieldElement:
        return self.scalar(1)

    def to_json(self) -> dict:
        return {"schema": "zipper-field-v1", "characteristic": self.characteristic,
                "parameters": list(self.parameters)}

    @classmethod
    def from_json(cls, value: object) -> RationalFunctionField:
        if (not isinstance(value, dict)
                or set(value) != {"schema", "characteristic", "parameters"}
                or value["schema"] != "zipper-field-v1"
                or not isinstance(value["parameters"], list)):
            raise ValueError("invalid coefficient-field declaration")
        return cls(tuple(value["parameters"]), value["characteristic"])


@dataclass(frozen=True, slots=True, init=False, eq=False)
class FieldElement:
    """An immutable exact rational function with a declared coefficient field.

    No polynomial gcd is used. Equal fractions can retain different serialized
    representatives and exceptional loci. Values are unhashable because equality
    is algebraic, rather than equality of numerator/denominator tuples.
    """

    field: RationalFunctionField
    numerator: Terms
    denominator: Terms
    # The runtime spelling for an unhashable value conflicts with object's
    # method annotation in mypy; retaining identity hashing would be incorrect.
    __hash__ = None  # type: ignore[assignment]

    def __init__(self, field: RationalFunctionField,
                 numerator: Mapping[Exponent, Scalar],
                 denominator: Mapping[Exponent, Scalar] | None = None) -> None:
        if not isinstance(field, RationalFunctionField):
            raise TypeError("a rational function requires a coefficient field")
        unit = (((0,) * len(field.parameters), _base(1, field.characteristic)),)
        top = _clean(field, numerator)
        bottom = unit if denominator is None else _clean(field, denominator)
        if not bottom:
            raise ZeroDivisionError("a rational-function denominator must be nonzero")
        if not top:
            bottom = unit
        else:
            leading = bottom[-1][1]
            inverse = (Fraction(1, leading) if not field.characteristic
                       else pow(int(leading), -1, field.characteristic))
            top = _clean(field, {monomial: coefficient * inverse for monomial, coefficient in top})
            bottom = _clean(field, {monomial: coefficient * inverse for monomial, coefficient in bottom})
        object.__setattr__(self, "field", field)
        object.__setattr__(self, "numerator", top)
        object.__setattr__(self, "denominator", bottom)

    def _lift(self, other: object) -> FieldElement | None:
        if isinstance(other, FieldElement):
            if other.field != self.field:
                raise TypeError("cannot mix coefficient fields")
            return other
        if isinstance(other, (int, Fraction)):
            return self.field.scalar(other)
        return None

    def __add__(self, other):
        value = self._lift(other)
        if value is None:
            return NotImplemented
        if self.is_zero():
            return value
        if value.is_zero():
            return self
        if self.denominator == value.denominator:
            top = _add(self.field, self.numerator, value.numerator)
            bottom = self.denominator
        else:
            top = _add(self.field, _mul(self.field, self.numerator, value.denominator),
                       _mul(self.field, value.numerator, self.denominator))
            bottom = _mul(self.field, self.denominator, value.denominator)
        return FieldElement(self.field, dict(top), dict(bottom))

    __radd__ = __add__

    def __neg__(self) -> FieldElement:
        return FieldElement(self.field, {m: -c for m, c in self.numerator}, dict(self.denominator))

    def __sub__(self, other):
        value = self._lift(other)
        return NotImplemented if value is None else self + (-value)

    def __rsub__(self, other):
        value = self._lift(other)
        return NotImplemented if value is None else value + (-self)

    def __mul__(self, other):
        value = self._lift(other)
        if value is None:
            return NotImplemented
        if self.is_zero() or value.is_zero():
            return self.field.zero
        unit = self.field.one.denominator
        if value.numerator == value.denominator == unit:
            return self
        if self.numerator == self.denominator == unit:
            return value
        return FieldElement(self.field,
                            dict(_mul(self.field, self.numerator, value.numerator)),
                            dict(_mul(self.field, self.denominator, value.denominator)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        value = self._lift(other)
        if value is None:
            return NotImplemented
        if value.is_zero():
            raise ZeroDivisionError("division by zero in the coefficient field")
        if self.is_zero():
            return self.field.zero
        return FieldElement(self.field,
                            dict(_mul(self.field, self.numerator, value.denominator)),
                            dict(_mul(self.field, self.denominator, value.numerator)))

    def __rtruediv__(self, other):
        value = self._lift(other)
        return NotImplemented if value is None else value / self

    def __pow__(self, exponent: int) -> FieldElement:
        if type(exponent) is not int:
            raise TypeError("rational-function powers must be integers")
        base = self
        if exponent < 0:
            base = self.field.one / self
            exponent = -exponent
        result = self.field.one
        while exponent:
            if exponent % 2:
                result = result * base
            exponent //= 2
            if exponent:
                base = base * base
        return result

    def __eq__(self, other: object) -> bool:
        try:
            value = self._lift(other)
        except TypeError:
            return False
        return (value is not None and
                ((self.numerator == value.numerator and self.denominator == value.denominator)
                 or _mul(self.field, self.numerator, value.denominator)
                 == _mul(self.field, value.numerator, self.denominator)))

    def is_zero(self) -> bool:
        return not self.numerator

    def specialize(self, values: Mapping[str, Scalar], *,
                   characteristic: int | None = None) -> FieldElement:
        """Evaluate this stored expression at a complete exact parameter point.

        Characteristic-zero expressions can also be reduced to a prime field.
        The rational constants and the stored denominator must be defined there.
        Parameterized output and changes between positive characteristics are
        not part of this operation.
        """
        if not isinstance(values, Mapping) or set(values) != set(self.field.parameters):
            raise ValueError("specialization requires exactly the declared parameter names")
        target = RationalFunctionField((), self.field.characteristic
                                       if characteristic is None else characteristic)
        if self.field.characteristic and target.characteristic != self.field.characteristic:
            raise ValueError("specialization cannot change positive characteristic")
        point = [_base(values[name], target.characteristic) for name in self.field.parameters]

        def evaluate(terms: Terms) -> Scalar:
            total: Scalar = 0
            for monomial, coefficient in terms:
                term = _base(coefficient, target.characteristic)
                for value, exponent in zip(point, monomial):
                    if target.characteristic:
                        term = (int(term) * pow(int(value), exponent, target.characteristic)) % target.characteristic
                    else:
                        term *= value ** exponent
                total = _base(total + term, target.characteristic)
            return total

        top = evaluate(self.numerator)
        bottom = evaluate(self.denominator)
        return target.scalar(top) / target.scalar(bottom)

    def to_json(self) -> dict:
        def encode(terms):
            return [[list(monomial), _encode_scalar(coefficient)]
                    for monomial, coefficient in terms]
        return {"schema": "zipper-field-element-v1", "field": self.field.to_json(),
                "numerator": encode(self.numerator), "denominator": encode(self.denominator)}

    @classmethod
    def from_json(cls, value: object) -> FieldElement:
        if (not isinstance(value, dict)
                or set(value) != {"schema", "field", "numerator", "denominator"}
                or value["schema"] != "zipper-field-element-v1"):
            raise ValueError("invalid rational-function encoding")
        field = RationalFunctionField.from_json(value["field"])

        def decode(terms):
            if not isinstance(terms, list):
                raise ValueError("sparse polynomial data must be a list")
            result = {}
            for item in terms:
                if (not isinstance(item, list) or len(item) != 2
                        or not isinstance(item[0], list)
                        or any(type(power) is not int or power < 0 for power in item[0])):
                    raise ValueError("invalid sparse polynomial term")
                monomial = tuple(item[0])
                if monomial in result:
                    raise ValueError("duplicate sparse polynomial monomial")
                result[monomial] = _decode_scalar(item[1])
            return result

        return cls(field, decode(value["numerator"]), decode(value["denominator"]))

    def __str__(self) -> str:
        def render(terms):
            parts = []
            for monomial, coefficient in terms:
                powers = [name if exponent == 1 else f"{name}^{exponent}"
                          for name, exponent in zip(self.field.parameters, monomial) if exponent]
                parts.append("*".join([str(coefficient)] + powers))
            return "+".join(parts) if parts else "0"
        if self.denominator == self.field.one.denominator:
            return render(self.numerator)
        return f"({render(self.numerator)})/({render(self.denominator)})"
