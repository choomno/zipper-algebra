"""Sparse free-associative polynomials over the exact coefficient backends."""
from __future__ import annotations

from .coefficients import Coef, QTCoef, _validate_characteristic

__all__ = ["Poly"]


class Poly:
    """A non-commutative polynomial. monomials are tuples of variable indices
    (1-based). The empty tuple () is the identity. Coefficients drop on collision."""
    __slots__ = ("terms", "char")

    def __init__(self, char: int):
        _validate_characteristic(char)
        self.terms: dict[tuple[int, ...], Coef] = {}
        self.char = char

    def add_term(self, mon: tuple[int, ...], coef: Coef) -> None:
        if not isinstance(coef, (Coef, QTCoef)) or coef.char != self.char:
            raise TypeError("a polynomial term must have the polynomial's coefficient characteristic")
        if not isinstance(mon, tuple) or any(type(i) is not int or i < 1 for i in mon):
            raise ValueError("monomials must be tuples of positive integer variable indices")
        if coef.is_zero():
            return
        existing = self.terms.get(mon)
        if existing is None:
            self.terms[mon] = coef
        else:
            new = existing + coef
            if new.is_zero():
                del self.terms[mon]
            else:
                self.terms[mon] = new

    def __add__(self, other: "Poly") -> "Poly":
        if not isinstance(other, Poly):
            return NotImplemented
        if self.char != other.char:
            raise TypeError("cannot mix polynomials of different characteristics")
        result = Poly(self.char)
        for mon, c in self.terms.items():
            result.add_term(mon, c)
        for mon, c in other.terms.items():
            result.add_term(mon, c)
        return result

    def __sub__(self, other: "Poly") -> "Poly":
        if not isinstance(other, Poly):
            return NotImplemented
        if self.char != other.char:
            raise TypeError("cannot mix polynomials of different characteristics")
        result = Poly(self.char)
        for mon, c in self.terms.items():
            result.add_term(mon, c)
        for mon, c in other.terms.items():
            result.add_term(mon, -c)
        return result

    def __neg__(self) -> "Poly":
        result = Poly(self.char)
        for mon, c in self.terms.items():
            result.add_term(mon, -c)
        return result

    def __mul__(self, other: "Poly") -> "Poly":
        if not isinstance(other, Poly):
            return NotImplemented
        if self.char != other.char:
            raise TypeError("cannot mix polynomials of different characteristics")
        result = Poly(self.char)
        for m1, c1 in self.terms.items():
            for m2, c2 in other.terms.items():
                result.add_term(m1 + m2, c1 * c2)
        return result

    def scalar_mul(self, c: Coef) -> "Poly":
        if not isinstance(c, (Coef, QTCoef)) or c.char != self.char:
            raise TypeError("a scalar must have the polynomial's coefficient characteristic")
        result = Poly(self.char)
        for mon, coef in self.terms.items():
            result.add_term(mon, coef * c)
        return result

    @staticmethod
    def scalar(c: Coef) -> "Poly":
        p = Poly(c.char)
        p.add_term((), c)
        return p

    @staticmethod
    def generator(idx: int, char: int) -> "Poly":
        p = Poly(char)
        p.add_term((idx,), Coef.one(char))
        return p

    @staticmethod
    def generator_over(idx: int, one) -> "Poly":
        """The generator with the given multiplicative identity (Coef or QTCoef)."""
        p = Poly(one.char)
        p.add_term((idx,), one)
        return p

    def to_msolve_string(self, var_names: list[str]) -> str:
        """Render in msolve syntax. Variable names are the strings the user
        chose (e.g. 'x1', 'y2'); 'mon' tuples are 1-based indices into them."""
        if not self.terms:
            return "0"
        # Sort terms: longer monomials first, then lex.
        items = sorted(self.terms.items(), key=lambda kv: (-len(kv[0]), kv[0]))
        parts: list[str] = []
        for mon, coef in items:
            cs = str(coef)
            # mon is a tuple of 1-based var indices
            mon_str = "*".join(var_names[i - 1] for i in mon) if mon else ""
            if mon_str == "":
                term = cs
            elif cs == "1":
                term = mon_str
            elif cs == "-1":
                term = "-" + mon_str
            else:
                term = f"{cs}*{mon_str}"
            # Connect with sign
            if parts and not term.startswith("-"):
                parts.append("+" + term)
            else:
                parts.append(term)
        return "".join(parts)
