#!/usr/bin/env python3
"""AtomSet: the denominator atoms of a presentation, as DATA rather than a constant.

`reconstruction_theory.tex`, prop:atoms.  Run Buchberger--Mora on the presentation over
K = Q(q,t); let d in Z[q,t] be a common denominator of every basis coefficient and of every
leading coefficient inverted (used as a pivot) during the run.  The ATOMS are the distinct
irreducible factors of d.  Then every coefficient the lift has to recover is

        c(q,t) = N(q,t) / prod_i a_i^{e_i},      N in Z[q,t],

so the atoms ARE the hypothesis class of `fit_rational`, and the exponent box is its size.

They are a property of the presentation AND the monomial order.  A hard-coded set is therefore
wrong for any other quotient, and a wrong set cannot be repaired by sampling more points --
`lift_until_verified.py` failure mode (b).  This module makes the set a computed input.

How they are obtained.  No vetted engine runs the elimination over Q(q,t) and exposes its
pivots: f4ncgb has them (linear_algebra.cpp) but only at a specialization, where a pivot is an
integer; Singular letterplace clears content and discards them, and factoring the leading
coefficients of its reduced basis provably misses atoms (INDEX.md, measured 2026-09-16).  So
the atoms are inferred from what the pipeline already produces -- exactly the fallback the
planned-feature note prescribes.  At a grid point p every proof coefficient is a rational whose
denominator is a product of specialized pivots, hence divisible by a(p) for each atom a that
occurs.  Across many points that pins a(q,t) down: a polynomial of small bidegree and small
integer coefficients that divides the residual denominator at EVERY point is an atom, and a
wrong polynomial fails divisibility somewhere with overwhelming probability.

This inference is SEARCH.  The PROOF is what the pipeline already does with the result: the
lift built from the atom set must pass Phase 5's symbolic identity check over Q(q,t) with q,t
transcendental.  Discover by specialization, prove by transcendence.

Backward compatibility.  `AtomSet.DEFAULT` is the six atoms the pipeline has always used, in
their historical positional order with the historical exponent box, so that every existing
lift file, cofactor file and denominator string keeps its meaning byte for byte.
"""
from __future__ import annotations

import json
import math
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
        self.coefs = {k: int(v) for k, v in coefs.items() if int(v) != 0}

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

    def __eq__(self, other):
        return isinstance(other, Atom) and self.coefs == other.coefs

    def __hash__(self):
        return hash(tuple(sorted(self.coefs.items())))

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
        assert len(self.max_exp) == len(self.atoms)
        self.provenance = provenance

    # --- the denominator algebra used by every phase -------------------------
    def __len__(self):
        return len(self.atoms)

    def names(self):
        return [a.name for a in self.atoms]

    def eval_denom(self, exps, q, t):
        out = 1
        for a, e in zip(self.atoms, exps):
            if e:
                out *= a.eval(q, t) ** e
        return out

    def denom_str(self, exps) -> str:
        parts = []
        for a, e in zip(self.atoms, exps):
            if not e:
                continue
            base = a.name if len(a.coefs) == 1 and abs(next(iter(a.coefs.values()))) == 1 \
                else f"({a.name})"
            parts.append(base if e == 1 else f"{base}^{e}")
        return "*".join(parts) if parts else "1"

    @staticmethod
    def lcm(e1, e2):
        return tuple(max(a, b) for a, b in zip(e1, e2))

    @staticmethod
    def mul(e1, e2):
        return tuple(a + b for a, b in zip(e1, e2))

    def zero(self):
        return tuple(0 for _ in self.atoms)

    def parse_denom_str(self, s: str):
        """Inverse of denom_str.  Accepts 'q^2*t*(q-1)' etc.  Unknown factor -> ValueError."""
        exps = [0] * len(self.atoms)
        s = s.strip()
        if s in ("", "1"):
            return tuple(exps)
        by_name = {a.name: k for k, a in enumerate(self.atoms)}
        for piece in s.split("*"):
            piece = piece.strip()
            e = 1
            # split a trailing ^e that is not inside parentheses
            if "^" in piece and not piece.endswith(")"):
                base, _, ex = piece.rpartition("^")
                piece, e = base, int(ex)
            if piece.startswith("(") and piece.endswith(")"):
                piece = piece[1:-1]
            if piece not in by_name:
                raise ValueError(f"unknown denominator factor {piece!r}; atoms are {self.names()}")
            exps[by_name[piece]] += e
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


# ----------------------------------------------------------------------------
# Inference from sampled denominators
# ----------------------------------------------------------------------------

def _strip(r: int, vals: Iterable[int]) -> int:
    """Divide r by each value as many times as it divides (values of +-1 or 0 are skipped)."""
    for v in vals:
        v = abs(v)
        if v <= 1:
            continue
        while r % v == 0:
            r //= v
    return r


def _candidate_polys(max_bideg=(2, 1), coef_bound=2):
    """Primitive polynomials of bidegree <= max_bideg with |coefficients| <= coef_bound,
    non-constant, sign-normalised (leading monomial in (q-degree, t-degree) lex has c > 0),
    in ascending (total degree, coefficient size)."""
    mons = [(i, j) for i in range(max_bideg[0] + 1) for j in range(max_bideg[1] + 1)]
    mons.sort(reverse=True)                      # leading monomial first
    rng = range(-coef_bound, coef_bound + 1)
    seen = set()
    out = []
    for cs in product(rng, repeat=len(mons)):
        coefs = {m: c for m, c in zip(mons, cs) if c}
        if not coefs or list(coefs) == [(0, 0)]:
            continue
        lead = max(coefs)                        # lex on (i, j)
        if coefs[lead] < 0:
            continue                             # its negative is enumerated instead
        if math.gcd(*[abs(c) for c in coefs.values()]) != 1:
            continue                             # not primitive
        key = tuple(sorted(coefs.items()))
        if key in seen:
            continue
        seen.add(key)
        deg = max(i + j for i, j in coefs)
        size = sum(abs(c) for c in coefs.values())
        out.append((deg, size, len(coefs), coefs))
    out.sort(key=lambda x: (x[0], x[1], x[2]))
    return [c for _, _, _, c in out]


def infer(values, seed: AtomSet | None = None, max_bideg=(2, 1), coef_bound=2,
          min_support=3, verbose=True):
    """Infer the atom set from sampled coefficient values.

    values : iterable of ((q, t), Fraction) -- every coefficient the lift will have to fit,
             at every grid point.  Only the DENOMINATORS are used.
    seed   : atoms assumed present (their values are stripped first).  None = start from
             the empty set, which is the honest test of the method.

    Returns (atomset, report).  atomset = seed + discovered, in discovery order, so
    positions of seed atoms are unchanged.  report['residual'] lists the points where a
    denominator remained unexplained; non-empty means the class must be widened
    (larger max_bideg / coef_bound) and is exactly the situation more sampling cannot fix.
    """
    seed = seed if seed is not None else AtomSet([], [])
    # Pool residual denominators by point: R(p) = lcm over slots of the stripped denominator.
    pooled: dict[tuple, int] = {}
    for (q, t), val in values:
        den = Fraction(val).denominator
        r = _strip(den, (a.eval(q, t) for a in seed.atoms))
        if r > 1:
            pooled[(q, t)] = (pooled.get((q, t), 1) * r) // math.gcd(pooled.get((q, t), 1), r)
    pts = sorted(pooled)
    found: list[Atom] = []
    report = {"seed": seed.names(), "points_with_residual": len(pts), "found": [],
              "residual": {}, "search": {"max_bideg": max_bideg, "coef_bound": coef_bound}}
    if not pts:
        if verbose:
            print(f"[atoms] seed {seed.names()} explains every denominator at every point")
        return seed, report

    cands = _candidate_polys(max_bideg, coef_bound)
    if verbose:
        print(f"[atoms] {len(pts)} points carry residual denominators after stripping "
              f"{seed.names() or 'nothing'}; searching {len(cands)} candidate polynomials "
              f"of bidegree <= {max_bideg}, |coef| <= {coef_bound}")
    known = set(seed.atoms)
    while True:
        best, best_score = None, 0.0
        for coefs in cands:
            a = Atom(poly_str(coefs), coefs)
            if a in known:
                continue
            support, score, ok = 0, 0.0, True
            for p in pts:
                v = abs(a.eval(*p))
                if v <= 1:
                    continue                     # uninformative at this point
                if pooled[p] % v != 0:
                    ok = False
                    break
                support += 1
                score += math.log(v)
            if ok and support >= min_support and score > best_score:
                best, best_score = a, score
        if best is None:
            break
        found.append(best)
        known.add(best)
        if verbose:
            print(f"[atoms]   found {best.name}")
        # strip it and drop points now fully explained
        for p in list(pts):
            pooled[p] = _strip(pooled[p], [best.eval(*p)])
            if pooled[p] == 1:
                pts.remove(p)
        if not pts:
            break

    report["found"] = [a.name for a in found]
    report["residual"] = {f"{p[0]},{p[1]}": pooled[p] for p in pts}
    if verbose:
        if pts:
            print(f"[atoms] WIDEN THE CLASS: {len(pts)} points keep unexplained denominators, "
                  f"e.g. {pts[0]} -> {pooled[pts[0]]}. Raise max_bideg / coef_bound.")
        else:
            print(f"[atoms] every denominator explained by {seed.names() + report['found']}")
    return seed.extended(found, provenance=f"inferred from {len(values) if hasattr(values, '__len__') else '?'} "
                                             f"sampled coefficients; seed {seed.names()}"), report


if __name__ == "__main__":
    # Smoke test on synthetic data: denominators built from the historical six at random points.
    import random
    random.seed(1)
    six = AtomSet.DEFAULT
    pts = [(q, t) for q in (3, 5, 7, 11, 13, 17, 19, 23) for t in (3, 5, 7, 11, 13, 17) if q != t]
    vals = []
    for p in pts:
        for exps in [(1, 2, 1, 0, 0, 0), (0, 0, 0, 1, 1, 0), (2, 0, 0, 0, 1, 1), (0, 1, 1, 1, 0, 1)]:
            vals.append((p, Fraction(random.randint(1, 97), six.eval_denom(exps, *p))))
    got, rep = infer(vals, seed=None)
    print(got, rep["residual"] == {})
    assert set(got.names()) == set(six.names()), got.names()
    print("OK: rediscovered the six from denominators alone")
