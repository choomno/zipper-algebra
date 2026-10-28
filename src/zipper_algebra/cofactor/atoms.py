#!/usr/bin/env python3
"""Structured denominator factors for exact rational reconstruction.

An Atom is a nonzero polynomial in Z[q,t]. An AtomSet is an ordered family of these
polynomials together with exponent search bounds. Its factors need not be irreducible.
A represented coefficient has the form

    N(q,t) / product_i atom_i(q,t)**e_i, with N in Q[q,t].

The atom set and bounds specify a finite search family for denominators. Inference from
sampled coefficients proposes a family; it does not prove that the family contains the
coefficient functions of a generic basis. Reconstructed candidates must pass held-out
checks, and the final mathematical conclusion comes from exact certificate verification.

Exact univariate reconstruction is implemented in univariate.py; optional sampled
factor inference is isolated in inference.py. Both preserve exact sample validation
and supply candidates for a separate final identity check.

AtomSet.DEFAULT retains the historical six factors, positional order and search bounds,
so preserved certificate exponent vectors keep their interpretation. AtomSet.UNITS is the
historical q,t seed for Laurent coefficient presentations, rather than a universal required
factor set. Arbitrary nonzero polynomial factors are accepted by certificate verification.
"""
from __future__ import annotations

import json
import re
from fractions import Fraction
from itertools import product
from typing import Iterable

# ----------------------------------------------------------------------------
# Atom: a polynomial in Z[q,t] as {(i, j): c} meaning sum c * q^i * t^j
# ----------------------------------------------------------------------------

class Atom:
    __slots__ = ("name", "coefs")

    def __init__(self, name: str, coefs: dict):
        self.name = name
        self.coefs = dict(coefs)
        self.validate()
        self.coefs = {k: int(v) for k, v in self.coefs.items() if v != 0}

    def validate(self):
        """A denominator factor must be a nonzero polynomial in Z[q,t].

        Irreducibility is useful for reconstruction, but is not needed for verification.
        Check again at verification boundaries: these objects are mutable.
        """
        if not isinstance(self.name, str) or not self.name.strip() or self.name == "1":
            raise ValueError("an atom needs a nonempty name other than the reserved name '1'")
        depth = 0
        for character in self.name:
            depth += (character == "(") - (character == ")")
            if depth < 0:
                raise ValueError("atom names must have balanced parentheses")
        if depth:
            raise ValueError("atom names must have balanced parentheses")
        if not isinstance(self.coefs, dict):
            raise ValueError("atom coefficients must be a polynomial dictionary")
        for mon, c in self.coefs.items():
            if (not isinstance(mon, tuple) or len(mon) != 2
                    or any(type(e) is not int or e < 0 for e in mon)):
                raise ValueError("atom monomials must be pairs of nonnegative integers")
            if (isinstance(c, bool) or not isinstance(c, (int, Fraction))
                    or Fraction(c).denominator != 1):
                raise ValueError("atom coefficients must be exact integers")
        if not any(self.coefs.values()):
            raise ValueError("a denominator atom must be a nonzero polynomial")

    def eval(self, q, t):
        return sum(c * q ** i * t ** j for (i, j), c in self.coefs.items())

    def bidegree(self):
        return (max((i for i, _ in self.coefs), default=0),
                max((j for _, j in self.coefs), default=0))

    def to_json(self):
        return {"name": self.name, "coefs": [[i, j, c] for (i, j), c in sorted(self.coefs.items())]}

    @classmethod
    def from_json(cls, d):
        return cls(d["name"], {(i, j): c for i, j, c in d["coefs"]})

    def _signature(self):
        """Coefficients up to overall sign: q-t and t-q are the same atom (a unit apart)."""
        if not self.coefs:
            return ()
        lead = max(self.coefs)
        sgn = 1 if self.coefs[lead] > 0 else -1
        return tuple(sorted((m, sgn * c) for m, c in self.coefs.items()))

    def __eq__(self, other):
        return isinstance(other, Atom) and self._signature() == other._signature()

    def __hash__(self):
        return hash(self._signature())

    def __repr__(self):
        return f"Atom({self.name})"


def poly_str(coefs: dict) -> str:
    """Canonical text for {(i,j): c}: descending in q then t, e.g. 't+q^2' -> 'q^2+t'."""
    terms = []
    for (i, j), c in sorted(coefs.items(), reverse=True):
        mon = "*".join(s for s in (("q" if i == 1 else f"q^{i}" if i else ""),
                                   ("t" if j == 1 else f"t^{j}" if j else "")) if s)
        if mon:
            body = mon if abs(c) == 1 else f"{abs(c)}*{mon}"
        else:
            body = str(abs(c))
        terms.append(("-" if c < 0 else "+") + body)
    s = "".join(terms).lstrip("+")
    return s or "0"


# ----------------------------------------------------------------------------
# AtomSet: an ordered list of atoms + exponent box; denominators are exponent tuples
# ----------------------------------------------------------------------------

class AtomSet:
    def __init__(self, atoms: Iterable[Atom], max_exp: Iterable[int] | None = None,
                 provenance: str = ""):
        self.atoms = list(atoms)
        self.max_exp = list(max_exp) if max_exp is not None else [3] * len(self.atoms)
        self.provenance = provenance
        self.validate()

    def validate(self):
        for atom in self.atoms:
            if not isinstance(atom, Atom):
                raise ValueError("an atom set must contain Atom objects")
            atom.validate()
        if len(set(self.names())) != len(self.atoms):
            raise ValueError("denominator atom names must be distinct")
        self.validate_exponents(self.max_exp)

    def validate_exponents(self, exps):
        if not isinstance(exps, (tuple, list)) or len(exps) != len(self.atoms):
            raise ValueError(f"denominator exponent vector must have length {len(self.atoms)}")
        if any(type(e) is not int or e < 0 for e in exps):
            raise ValueError("denominator exponents must be nonnegative integers")

    # --- the denominator algebra used by every phase -------------------------
    def __len__(self):
        return len(self.atoms)

    def names(self):
        return [a.name for a in self.atoms]

    def eval_denom(self, exps, q, t):
        self.validate_exponents(exps)
        out = 1
        for a, e in zip(self.atoms, exps):
            if e:
                out *= a.eval(q, t) ** e
        return out

    def denom_str(self, exps) -> str:
        self.validate_exponents(exps)
        parts = []
        for a, e in zip(self.atoms, exps):
            if not e:
                continue
            base = a.name if re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", a.name) \
                else f"({a.name})"
            parts.append(base if e == 1 else f"{base}^{e}")
        return "*".join(parts) if parts else "1"

    @staticmethod
    def lcm(e1, e2):
        AtomSet._validate_pair(e1, e2)
        return tuple(max(a, b) for a, b in zip(e1, e2))

    @staticmethod
    def mul(e1, e2):
        AtomSet._validate_pair(e1, e2)
        return tuple(a + b for a, b in zip(e1, e2))

    @staticmethod
    def _validate_pair(e1, e2):
        if (not isinstance(e1, (tuple, list)) or not isinstance(e2, (tuple, list))
                or len(e1) != len(e2)
                or any(type(e) is not int or e < 0 for e in (*e1, *e2))):
            raise ValueError("denominator vectors must have equal lengths and nonnegative integer powers")

    def zero(self):
        return tuple(0 for _ in self.atoms)

    def parse_denom_str(self, s: str):
        """Inverse of denom_str.  Accepts 'q^2*t*(q-1)' etc.  Unknown factor -> ValueError."""
        exps = [0] * len(self.atoms)
        s = s.strip()
        if s in ("", "1"):
            return tuple(exps)
        by_name = {a.name: k for k, a in enumerate(self.atoms)}
        # Atom names themselves can contain products, e.g. (2*t+1).
        factors, start, depth = [], 0, 0
        for i, ch in enumerate(s):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth < 0:
                    raise ValueError(f"unbalanced denominator {s!r}")
            elif ch == "*" and depth == 0:
                factors.append(s[start:i])
                start = i + 1
        if depth:
            raise ValueError(f"unbalanced denominator {s!r}")
        factors.append(s[start:])
        for piece in factors:
            piece = piece.strip()
            e = 1
            # split a trailing ^e that is not inside parentheses
            if "^" in piece and not piece.endswith(")"):
                base, _, ex = piece.rpartition("^")
                piece, e = base, int(ex)
                if e < 0:
                    raise ValueError("denominator exponents must be nonnegative integers")
            if piece.startswith("(") and piece.endswith(")"):
                piece = piece[1:-1]
            if piece not in by_name:
                raise ValueError(f"unknown denominator factor {piece!r}; atoms are {self.names()}")
            exps[by_name[piece]] += e
        self.validate_exponents(exps)
        return tuple(exps)

    def candidates(self):
        """All exponent tuples inside the box, ascending total degree -- the search order of
        fit_rational.  For the DEFAULT set this reproduces the historical enumeration exactly."""
        out = [(exps, sum(exps)) for exps in product(*(range(m + 1) for m in self.max_exp))]
        out.sort(key=lambda x: x[1])
        return [e for e, _ in out]

    # --- persistence ---------------------------------------------------------
    def to_json(self):
        return {"atoms": [a.to_json() for a in self.atoms],
                "max_exp": self.max_exp,
                "provenance": self.provenance}

    @classmethod
    def from_json(cls, d):
        return cls([Atom.from_json(a) for a in d["atoms"]], d.get("max_exp"), d.get("provenance", ""))

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.to_json(), f, indent=1)

    @classmethod
    def load(cls, path):
        with open(path) as f:
            return cls.from_json(json.load(f))

    def extended(self, new_atoms: Iterable[Atom], max_exp_new: int = 3, provenance=""):
        """A copy with new atoms appended -- positions of existing atoms are unchanged, so
        exponent tuples written under the old set stay valid under the new one (zero-padded)."""
        atoms, box = list(self.atoms), list(self.max_exp)
        for a in new_atoms:
            if a not in atoms:
                atoms.append(a)
                box.append(max_exp_new)
        return AtomSet(atoms, box, provenance or self.provenance)

    def __repr__(self):
        return f"AtomSet({self.names()}, box={self.max_exp})"


# Historical seed for presentations over a Laurent base where q and t are units.
# The general AtomSet API also accepts other factor families, including the empty family.
AtomSet.UNITS = AtomSet([Atom("q", {(1, 0): 1}), Atom("t", {(0, 1): 1})], max_exp=[5, 4],
                        provenance="units of Z[q^,t^]")

# The six atoms the X_3 pipeline has always used, in the positional order every existing
# lift/cofactor file encodes, with the historical exponent box of gen_denoms_sorted().
AtomSet.DEFAULT = AtomSet(
    [Atom("q",     {(1, 0): 1}),
     Atom("t",     {(0, 1): 1}),
     Atom("q-1",   {(1, 0): 1, (0, 0): -1}),
     Atom("t-q",   {(0, 1): 1, (1, 0): -1}),
     Atom("t+q^2", {(0, 1): 1, (2, 0): 1}),
     Atom("q+1",   {(1, 0): 1, (0, 0): 1})],
    max_exp=[5, 4, 2, 1, 2, 1],
    provenance="historical six (prop:atoms, X_3 at n=4)")


def infer(values_by_slot: dict, seed: AtomSet | None = None, max_bideg=(2, 1), coef_bound=2,
          min_lines=3, verbose=True):
    """Propose denominator factors; compatibility entry point for inference.infer."""
    from .inference import infer as infer_denominators
    return infer_denominators(values_by_slot, seed=seed, max_bideg=max_bideg,
                              coef_bound=coef_bound, min_lines=min_lines, verbose=verbose)
