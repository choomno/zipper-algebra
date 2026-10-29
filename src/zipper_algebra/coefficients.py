"""Exact rational, residue-ring and bivariate symbolic coefficient arithmetic."""
from __future__ import annotations

from fractions import Fraction

__all__ = ["Coef", "QTCoef", "NumericRing", "SymbolicRing"]


def _validate_characteristic(char):
    if type(char) is not int or char < 0 or char == 1:
        raise ValueError("characteristic must be zero or an integer at least two")


class Coef:
    """Exact scalars over Q or Z/mZ; the latter is GF(m) when m is prime."""
    __slots__ = ("v", "char")
    v: int | Fraction

    def __init__(self, v: int | Fraction, char: int) -> None:
        _validate_characteristic(char)
        if isinstance(v, bool) or not isinstance(v, (int, Fraction)):
            raise TypeError("coefficients must be exact integers or Fractions")
        self.char = char
        value = Fraction(v)
        if char == 0:
            self.v = value
        else:
            try:
                inverse = pow(value.denominator, -1, char)
            except ValueError:
                raise ZeroDivisionError(
                    f"coefficient denominator is not invertible modulo {char}") from None
            self.v = value.numerator * inverse % char

    def _check_parent(self, other):
        if self.char != other.char:
            raise TypeError(f"cannot mix coefficients of characteristics {self.char} and {other.char}")

    @staticmethod
    def zero(char: int) -> "Coef":
        return Coef(0, char)

    @staticmethod
    def one(char: int) -> "Coef":
        return Coef(1, char)

    def __add__(self, other: "Coef") -> "Coef":
        if not isinstance(other, Coef):
            return NotImplemented
        self._check_parent(other)
        return Coef(self.v + other.v, self.char)

    def __sub__(self, other: "Coef") -> "Coef":
        if not isinstance(other, Coef):
            return NotImplemented
        self._check_parent(other)
        return Coef(self.v - other.v, self.char)

    def __mul__(self, other: "Coef") -> "Coef":
        if not isinstance(other, Coef):
            return NotImplemented
        self._check_parent(other)
        return Coef(self.v * other.v, self.char)

    def __neg__(self) -> "Coef":
        return Coef(-self.v, self.char)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Coef):
            return NotImplemented
        return self.v == other.v and self.char == other.char

    def is_zero(self) -> bool:
        return self.v == 0

    def __str__(self) -> str:
        if self.char == 0:
            assert isinstance(self.v, Fraction)
            return (str(self.v.numerator) if self.v.denominator == 1
                    else f"{self.v.numerator}/{self.v.denominator}")
        return str(self.v)


class QTCoef:
    """Scalar backend with q, t SYMBOLIC: an element of Q[q,t] as {(deg_q, deg_t): Fraction}.

    Same operator protocol as Coef, so Poly and build_relations_ring work unchanged over it.
    Interoperates with a characteristic-0 Coef by promoting it to a constant. Exists so the
    trusted emitter can be handed symbols directly; before 2026-09-17 the presentation over
    Q(q,t) was hand-copied in symbolic_identity_check.py and reconstructed by sampling in
    papers/generic/fit_symbolic.py, because this class did not."""
    __slots__ = ("v", "char")
    v: dict[tuple[int, int], Fraction]

    def __init__(self, v):
        self.char = 0
        if isinstance(v, dict):
            self.v = {}
            for monomial, coefficient in v.items():
                if (not isinstance(monomial, tuple) or len(monomial) != 2
                        or any(type(e) is not int or e < 0 for e in monomial)):
                    raise ValueError("symbolic monomials need two nonnegative integer exponents")
                if isinstance(coefficient, bool) or not isinstance(coefficient, (int, Fraction)):
                    raise TypeError("symbolic coefficients must be exact integers or Fractions")
                if coefficient:
                    self.v[monomial] = Fraction(coefficient)
        else:
            if isinstance(v, bool) or not isinstance(v, (int, Fraction)):
                raise TypeError("symbolic coefficients must be exact integers or Fractions")
            f = Fraction(v)
            self.v = {(0, 0): f} if f != 0 else {}

    @staticmethod
    def _lift(other):
        if isinstance(other, QTCoef):
            return other
        if isinstance(other, Coef):
            if other.char != 0:
                raise TypeError("cannot mix a GF(p) Coef with symbolic q,t")
            return QTCoef(other.v)
        if isinstance(other, (int, Fraction)):
            return QTCoef(other)
        return None

    def __add__(self, other):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        out = dict(self.v)
        for k, c in o.v.items():
            out[k] = out.get(k, Fraction(0)) + c
        return QTCoef(out)

    __radd__ = __add__

    def __neg__(self):
        return QTCoef({k: -c for k, c in self.v.items()})

    def __sub__(self, other):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        return self + (-o)

    def __rsub__(self, other):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        return o + (-self)

    def __mul__(self, other):
        o = self._lift(other)
        if o is None:
            return NotImplemented
        out: dict[tuple[int, int], Fraction] = {}
        for (a, b), c in self.v.items():
            for (i, j), d in o.v.items():
                out[(a + i, b + j)] = out.get((a + i, b + j), Fraction(0)) + c * d
        return QTCoef(out)

    __rmul__ = __mul__

    def __eq__(self, other):
        try:
            o = self._lift(other)
        except TypeError:
            return False
        return o is not None and self.v == o.v

    def is_zero(self) -> bool:
        return not self.v

    def bidegree(self):
        return (max((i for i, _ in self.v), default=0), max((j for _, j in self.v), default=0))

    def __str__(self) -> str:
        """Singular-readable: (c)*q^i*t^j terms joined by '+' -- generic_dim.coeff_str's format."""
        if not self.v:
            return "0"
        parts = []
        for (i, j), c in sorted(self.v.items()):
            if c.denominator != 1:
                raise ValueError(f"non-integer coefficient {c}: emission assumes Z[q,t]")
            m = []
            if i:
                m.append("q" if i == 1 else f"q^{i}")
            if j:
                m.append("t" if j == 1 else f"t^{j}")
            body = "*".join(m)
            parts.append(f"({c.numerator})" + ("*" + body if body else ""))
        return "+".join(parts)


class NumericRing:
    """Exact parameter values over Q (char 0) or the residue ring Z/mZ.

    A positive modulus describes a field only when it is prime. Callers using
    field algorithms or emitting input for a prime-field solver must supply a prime.
    """
    def __init__(self, q_int: int, t_int: int, char: int):
        self.char = char
        self.q = Coef(q_int, char)
        self.t = Coef(t_int, char)
        self.one = Coef.one(char)
        self.q_int, self.t_int = q_int, t_int

    def scalar(self, v) -> Coef:
        return Coef(v, self.char)

    def describe(self) -> str:
        return f"(q,t)=({self.q_int},{self.t_int}) over " + ("Q" if self.char == 0 else f"Z/{self.char}Z")


class SymbolicRing:
    """q, t transcendental: coefficients in Z[q,t] (QTCoef). Results over Q(q,t) come from here."""
    char = 0

    def __init__(self):
        self.q = QTCoef({(1, 0): 1})
        self.t = QTCoef({(0, 1): 1})
        self.one = QTCoef(1)

    def scalar(self, v) -> QTCoef:
        return QTCoef(v)

    def describe(self) -> str:
        return "q,t symbolic over Q"
