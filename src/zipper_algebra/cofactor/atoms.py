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

Univariate reconstruction uses exact rational arithmetic, rejects underdetermined fits,
cancels polynomial gcds, and checks the reduced candidate at every supplied location.
Repeated locations contribute one constraint; inconsistent repeated values are invalid.

AtomSet.DEFAULT retains the historical six factors, positional order and search bounds,
so preserved certificate exponent vectors keep their interpretation. AtomSet.UNITS is the
historical q,t seed for Laurent coefficient presentations, rather than a universal required
factor set. Arbitrary nonzero polynomial factors are accepted by certificate verification.
"""
from __future__ import annotations

import json
import math
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


# ----------------------------------------------------------------------------
# Inference: reconstruct each coefficient's denominator as a POLYNOMIAL along grid lines,
# then read the atoms off by polynomial divisibility.
#
# Why not integer divisibility of the sampled denominators?  Because the values of
# different atoms at a point share small prime factors (at (q,t)=(3,5): q-1 = 2, t-q = 2,
# q+1 = 4), so "strip a(p) from den(p) as often as it divides" over-strips and a wrong or
# reducible polynomial can look like an atom.  A univariate rational function through enough
# points is unique, and polynomial divisibility has no coincidences, so we work there.
# ----------------------------------------------------------------------------

def _poly_eval(coefs, x):
    return sum(c * x ** i for i, c in enumerate(coefs))


def _poly_divmod(num, den):
    """Exact univariate division over Fraction; polynomials as coefficient lists, low->high."""
    num = [Fraction(c) for c in num]
    den = [Fraction(c) for c in den]
    while den and den[-1] == 0:
        den.pop()
    if not den:
        raise ZeroDivisionError
    q = [Fraction(0)] * max(0, len(num) - len(den) + 1)
    r = num[:]
    while len(r) >= len(den) and any(r):
        if r[-1] == 0:
            r.pop(); continue
        k = len(r) - len(den)
        c = r[-1] / den[-1]
        q[k] = c
        for i, d in enumerate(den):
            r[k + i] -= c * d
        r.pop()
    while r and r[-1] == 0:
        r.pop()
    return q, r


def _solve(aug, ncols):
    """Exact Gauss-Jordan on an augmented matrix; returns a solution vector or None."""
    n = len(aug)
    rank, piv = 0, []
    for col in range(ncols):
        pv = next((r for r in range(rank, n) if aug[r][col] != 0), None)
        if pv is None:
            continue
        aug[rank], aug[pv] = aug[pv], aug[rank]
        inv = 1 / aug[rank][col]
        aug[rank] = [x * inv for x in aug[rank]]
        for r in range(n):
            if r != rank and aug[r][col] != 0:
                f = aug[r][col]
                aug[r] = [a - f * b for a, b in zip(aug[r], aug[rank])]
        piv.append(col); rank += 1
    for r in range(n):
        if all(aug[r][c] == 0 for c in range(ncols)) and aug[r][ncols] != 0:
            return None
    if rank != ncols:
        return None
    sol = [Fraction(0)] * ncols
    for r, col in enumerate(piv):
        sol[col] = aug[r][ncols]
    return sol


def exact_fraction(value):
    """An exact rational input, excluding bool and implicit binary-float coercion."""
    if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
        raise ValueError("sample coordinates and values must be exact rational numbers")
    return Fraction(value)


def normalise_points(points, dimensions):
    """Sorted, distinct exact samples; consistent repeats do not add information."""
    values = {}
    for point in points:
        if not isinstance(point, (tuple, list)) or len(point) != 2:
            raise ValueError("each sample must contain a coordinate and a value")
        location, value = point
        if dimensions == 1:
            location = exact_fraction(location)
        else:
            if not isinstance(location, (tuple, list)) or len(location) != dimensions:
                raise ValueError(f"sample locations need {dimensions} coordinates")
            location = tuple(exact_fraction(v) for v in location)
        value = exact_fraction(value)
        if location in values and values[location] != value:
            raise ValueError(f"conflicting values at sample location {location!r}")
        values[location] = value
    return sorted(values.items())


def _trim(poly):
    poly = [Fraction(c) for c in poly]
    while poly and poly[-1] == 0:
        poly.pop()
    return poly


def _reduce_univariate_fraction(numerator, denominator):
    """Cancel a polynomial gcd over Q and normalize the denominator to monic."""
    numerator, denominator = _trim(numerator), _trim(denominator)
    if not denominator:
        raise ValueError("a rational function needs a nonzero denominator")
    if not numerator:
        return [Fraction(0)], [Fraction(1)]
    a, b = numerator, denominator
    while b:
        _, remainder = _poly_divmod(a, b)
        a, b = b, _trim(remainder)
    common = [c / a[-1] for c in a]
    numerator, nr = _poly_divmod(numerator, common)
    denominator, dr = _poly_divmod(denominator, common)
    assert not _trim(nr) and not _trim(dr)
    numerator, denominator = _trim(numerator), _trim(denominator)
    lead = denominator[-1]
    return [c / lead for c in numerator], [c / lead for c in denominator]


def reconstruct_univariate(points, holdout=2):
    """Reconstruct a reduced rational candidate from distinct exact rational samples.

    Fit on all but `holdout` locations and predict the remaining ones. Every supplied
    location must be defined and correct after gcd cancellation. Return (N, D), low to
    high coefficients with D monic, or None when no determined candidate succeeds.
    Agreement on samples is a reconstruction candidate, not a symbolic identity proof.
    """
    if type(holdout) is not int or holdout < 1:
        raise ValueError("holdout must be a positive integer")
    points = normalise_points(points, 1)
    if len(points) < holdout + 3:
        return None
    fit = points[:-holdout]
    for total in range(len(fit) - 1):
        for dd in range(total + 1):
            dn = total - dd
            ncols = dn + 1 + dd
            aug = []
            for x, value in fit:
                row = [x ** i for i in range(dn + 1)]
                row += [-value * x ** j for j in range(dd)]
                row.append(value * x ** dd)
                aug.append(row)
            sol = _solve(aug, ncols)
            if sol is None:
                continue
            N, D = _reduce_univariate_fraction(sol[:dn + 1], sol[dn + 1:] + [Fraction(1)])
            if all((d := _poly_eval(D, x)) != 0 and _poly_eval(N, x) / d == value
                   for x, value in points):
                return N, D
    return None


def _specialise(coefs: dict, fixed: str, x0) -> list:
    """Bivariate {(i,j): c} restricted to q = x0 ('q') or t = x0 ('t'); returns a univariate
    coefficient list in the other variable, low->high."""
    out = {}
    for (i, j), c in coefs.items():
        if fixed == "t":
            out[i] = out.get(i, 0) + c * Fraction(x0) ** j
        else:
            out[j] = out.get(j, 0) + c * Fraction(x0) ** i
    deg = max(out) if out else 0
    poly = [out.get(k, Fraction(0)) for k in range(deg + 1)]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def _candidate_polys(max_bideg=(2, 1), coef_bound=2):
    """Primitive polynomials of bidegree <= max_bideg with |coefficients| <= coef_bound,
    non-constant, sign-normalised, in ascending (total degree, coefficient size)."""
    mons = [(i, j) for i in range(max_bideg[0] + 1) for j in range(max_bideg[1] + 1)]
    mons.sort(reverse=True)
    rng = range(-coef_bound, coef_bound + 1)
    out = []
    for cs in product(rng, repeat=len(mons)):
        coefs = {m: c for m, c in zip(mons, cs) if c}
        if not coefs or list(coefs) == [(0, 0)]:
            continue
        if coefs[max(coefs)] < 0:
            continue
        if math.gcd(*[abs(c) for c in coefs.values()]) != 1:
            continue
        out.append((max(i + j for i, j in coefs), sum(abs(c) for c in coefs.values()), coefs))
    out.sort(key=lambda x: (x[0], x[1]))
    return [c for _, _, c in out]



# ----------------------------------------------------------------------------
# Residual factors: turn "widen the class" into "here is the atom you are missing"
# ----------------------------------------------------------------------------

def _primitive(coefs: list) -> tuple:
    """A rational coefficient list -> the primitive integer one, sign-normalised on the lead."""
    from fractions import Fraction as F
    cs = [F(c) for c in coefs]
    den = 1
    for c in cs:
        den = den * c.denominator // math.gcd(den, c.denominator)
    ints = [int(c * den) for c in cs]
    g = 0
    for v in ints:
        g = math.gcd(g, abs(v))
    if g:
        ints = [v // g for v in ints]
    if ints and ints[-1] < 0:
        ints = [-v for v in ints]
    return tuple(ints)


def residual_factors(lines, verbose=True):
    """Factor every leftover line denominator over Q and report the irreducible factors.

    `lines` is the [(fixed_var, x0, D)] left by infer() after stripping the atoms it accepted.
    The brute-force candidate search cannot be widened far -- it enumerates all coefficient
    vectors in a bidegree box, which is exponential -- so when it fails the useful thing is not a
    bigger box but the ACTUAL factors, which sympy gets by factoring each univariate residual.
    Returns {(fixed_var, primitive coefficient tuple): [x0, ...]}.
    """
    import sympy as sp
    x = sp.Symbol("x")
    out = {}
    for fixed, x0, D in lines:
        if len(D) < 2:
            continue
        poly = sp.Poly(list(reversed([sp.Rational(c) for c in D])), x)
        for f, _e in sp.factor_list(poly)[1]:
            cs = list(reversed(f.all_coeffs()))
            out.setdefault((fixed, _primitive(cs)), []).append(x0)
    if verbose and out:
        print(f"[atoms] residual irreducible factors ({len(out)} distinct):")
        for (fixed, cs), xs in sorted(out.items(), key=lambda kv: -len(kv[1]))[:8]:
            var = "q" if fixed == "t" else "t"
            print(f"[atoms]   {_fmt_univariate(cs, var)}   on {len(xs)} lines "
                  f"({'t' if fixed == 't' else 'q'} = {sorted(set(xs))[:5]}...)")
    return out


def _fmt_univariate(cs, var) -> str:
    parts = []
    for i, c in enumerate(cs):
        if not c:
            continue
        mon = "" if i == 0 else (var if i == 1 else f"{var}^{i}")
        body = f"{abs(c)}" if not mon else (mon if abs(c) == 1 else f"{abs(c)}*{mon}")
        parts.append(("-" if c < 0 else "+") + body)
    return "".join(parts).lstrip("+") or "0"


def atoms_from_residual_factors(factors, max_other_deg=4, min_lines=3, verbose=True):
    """Lift line factors to bivariate atoms.

    A factor of the residual on the line t = t0 is the restriction a(q, t0) of an atom a(q,t).
    Group the factors by their degree in q, and for each coefficient position interpolate its
    values across the lines as a polynomial in t (exactly, verified on held-out lines).  This
    is the same finite-difference argument the rest of the pipeline uses, one variable at a time,
    and it needs no bidegree box: the degree comes out of the data.
    """
    from fractions import Fraction as F
    by_deg = {}
    for (fixed, cs), xs in factors.items():
        by_deg.setdefault((fixed, len(cs)), []).append((cs, xs))
    found = []
    for (fixed, ln), group in sorted(by_deg.items()):
        pts = {}
        for cs, xs in group:
            for x0 in xs:
                pts.setdefault(x0, set()).add(cs)
        usable = {x0: next(iter(v)) for x0, v in pts.items() if len(v) == 1}
        if len(usable) < min_lines + 2:
            continue
        xs = sorted(usable)
        coefs = {}
        ok = True
        for pos in range(ln):
            data = [(x0, F(usable[x0][pos])) for x0 in xs]
            fit = _fit_univariate_exact(data, max_other_deg)
            if fit is None:
                ok = False
                break
            coefs[pos] = fit
        if not ok:
            continue
        # assemble: coefficient of q^pos is a polynomial in t (or vice versa)
        bi = {}
        for pos, poly in coefs.items():
            for d, c in poly.items():
                key = (pos, d) if fixed == "t" else (d, pos)
                bi[key] = bi.get(key, F(0)) + c
        bi = {k: v for k, v in bi.items() if v}
        if not bi or all(k == (0, 0) for k in bi):
            continue
        den = 1
        for c in bi.values():
            den = den * F(c).denominator // math.gcd(den, F(c).denominator)
        bi = {k: int(F(v) * den) for k, v in bi.items()}
        a = Atom(poly_str(bi), bi)
        if a not in found:
            found.append(a)
            if verbose:
                print(f"[atoms]   lifted to a bivariate atom: {a.name}")
    return found


def _fit_univariate_exact(data, max_deg):
    """Fit {x: value} by a polynomial in x of degree <= max_deg, exactly, with two held out."""
    from fractions import Fraction as F
    data = sorted(data)
    if len(data) < 3:
        return None
    fit, hold = data[:-2], data[-2:]
    for deg in range(max_deg + 1):
        if len(fit) < deg + 1:
            break
        rows = [[F(x) ** i for i in range(deg + 1)] + [F(v)] for x, v in fit[: deg + 1]]
        sol = _solve(rows, deg + 1)
        if sol is None:
            continue
        good = all(sum(sol[i] * F(x) ** i for i in range(deg + 1)) == F(v)
                   for x, v in data)
        if good:
            return {i: sol[i] for i in range(deg + 1) if sol[i]}
    return None


def infer(values_by_slot: dict, seed: AtomSet | None = None, max_bideg=(2, 1), coef_bound=2,
          min_lines=3, verbose=True):
    """Infer the atom set from sampled coefficient values.

    COST, measured 2026-09-17.  This is the expensive step of the pipeline at scale and the
    numbers should be known before it is called on something large.  At n=4 bigelow: 231 slots,
    924 line denominators, seconds.  At n=3 square-zero: 834 slots, 20066 line denominators, and
    the greedy candidate search is O(#candidates x #lines) -- 16k x 20k divisibility tests -- so
    it runs for tens of minutes.  Two consequences: (1) the result is persisted per
    (presentation, target) and re-read, never recomputed; (2) the candidate box cannot be widened
    to fix a miss, because the enumeration is (2*coef_bound+1)^((dq+1)*(dt+1)) -- that is what
    residual_factors/atoms_from_residual_factors are for.  If this needs to scale further, the
    fix is to run the greedy search on a SAMPLE of the lines and verify the accepted atoms on the
    rest, not to enlarge the box.

    values_by_slot : {slot: [((q, t), Fraction), ...]} -- for every coefficient the lift must
                     fit, its value at every grid point.  Grid points sharing a q or a t form
                     the lines along which denominators are reconstructed.
    seed           : atoms assumed present.  None = start empty (the honest test).

    Returns (atomset, report).  Seed positions are unchanged; discovered atoms are appended.
    report['unexplained'] > 0 means some reconstructed denominator has a factor outside the
    search box. This is a diagnostic for the sampled reconstruction, not a proof that
    additional samples cannot help; the search class may also need revision.
    """
    seed = seed if seed is not None else AtomSet([], [])
    # 1. reconstruct denominators along lines: list of (fixed_var, x0, D_univariate)
    lines = []
    n_rec = n_try = 0
    for slot, pts in values_by_slot.items():
        by_t, by_q = {}, {}
        for (q, t), v in pts:
            by_t.setdefault(t, []).append((q, v))
            by_q.setdefault(q, []).append((t, v))
        for fixed, groups in (("t", by_t), ("q", by_q)):
            for x0, line in groups.items():
                if len(line) < 6:
                    continue
                n_try += 1
                rec = reconstruct_univariate(line)
                if rec is None:
                    continue
                n_rec += 1
                D = rec[1]
                if len(D) > 1:
                    lines.append((fixed, x0, D))
    if verbose:
        print(f"[atoms] reconstructed {n_rec}/{n_try} line denominators over "
              f"{len(values_by_slot)} slots; {len(lines)} non-trivial")
    report = {"seed": seed.names(), "lines_reconstructed": n_rec, "lines_tried": n_try,
              "found": [], "multiplicity": {}, "unexplained": 0,
              "search": {"max_bideg": list(max_bideg), "coef_bound": coef_bound}}

    # 2. strip the seed atoms (polynomial division, recording multiplicities)
    mult = {a.name: 0 for a in seed.atoms}

    def strip(atom):
        m_max = 0
        for idx, (fixed, x0, D) in enumerate(lines):
            a1 = _specialise(atom.coefs, fixed, x0)
            if len(a1) < 2:
                continue
            m = 0
            while len(D) >= len(a1):
                qq, r = _poly_divmod(D, a1)
                if r:
                    break
                D, m = qq, m + 1
            lines[idx] = (fixed, x0, D)
            m_max = max(m_max, m)
        return m_max

    for a in seed.atoms:
        mult[a.name] = strip(a)

    # 3. search
    found = []
    known = set(seed.atoms)
    cands = _candidate_polys(max_bideg, coef_bound)
    while True:
        best, best_key = None, None
        for coefs in cands:
            a = Atom(poly_str(coefs), coefs)
            if a in known:
                continue
            hits, tvals, qvals = 0, set(), set()
            for fixed, x0, D in lines:
                if len(D) < 2:
                    continue
                a1 = _specialise(coefs, fixed, x0)
                if len(a1) < 2 or len(a1) > len(D):
                    continue
                if not _poly_divmod(D, a1)[1]:
                    hits += 1
                    (tvals if fixed == "t" else qvals).add(x0)
            dep_t = any(j for _, j in coefs)
            dep_q = any(i for i, _ in coefs)
            if hits < min_lines:
                continue
            if dep_t and len(tvals) < 2 and len(qvals) < 2:
                continue
            if dep_q and len(qvals) < 2 and len(tvals) < 2:
                continue
            key = (hits,)
            if best_key is None or key > best_key:
                best, best_key = a, key
        if best is None:
            break
        found.append(best); known.add(best)
        mult[best.name] = strip(best)
        if verbose:
            print(f"[atoms]   found {best.name}  (multiplicity <= {mult[best.name]}, "
                  f"{best_key[0]} lines)")
        if all(len(D) < 2 for _, _, D in lines):
            break

    left = sum(1 for _, _, D in lines if len(D) >= 2)
    report["found"] = [a.name for a in found]
    report["multiplicity"] = mult
    report["unexplained"] = left
    if left:
        # Do not just say "widen the box" -- the box cannot be widened far (the candidate search
        # is exponential in it). Factor what is left and name the missing atoms.
        try:
            rf = residual_factors(lines, verbose=verbose)
            lifted = atoms_from_residual_factors(rf, verbose=verbose)
        except Exception as e:                      # sympy missing, or a factorisation failure
            rf, lifted = {}, []
            if verbose:
                print(f"[atoms]   (residual factorisation unavailable: {e})")
        report["residual_factors"] = [_fmt_univariate(cs, "q" if f == "t" else "t")
                                      for (f, cs) in rf]
        if lifted:
            # Lifting returns distinct residual factors after all found atoms were stripped.
            for a in lifted:
                found.append(a)
                mult[a.name] = strip(a)
            left = sum(1 for _, _, D in lines if len(D) >= 2)
            report["found"] = [a.name for a in found]
            report["multiplicity"] = mult
        report["unexplained"] = left
    if verbose:
        if left:
            ex = next(D for _, _, D in lines if len(D) >= 2)
            print(f"[atoms] WIDEN THE CLASS: {left} reconstructed denominators keep a factor "
                  f"outside bidegree {max_bideg} / |coef| <= {coef_bound}, e.g. {[str(c) for c in ex]}")
        else:
            print(f"[atoms] every reconstructed denominator is a product of "
                  f"{seed.names() + report['found']}")
    out = seed.extended(found, provenance=f"inferred over {len(values_by_slot)} slots; seed {seed.names()}")
    # exponent box from observed multiplicities (+1 slack), never below the seed's box
    out.max_exp = [max(m0, mult.get(a.name, 0) + 1) for a, m0 in zip(out.atoms, out.max_exp)]
    return out, report


if __name__ == "__main__":
    # Smoke test on synthetic data: denominators built from the historical six at random points.
    import random
    random.seed(1)
    six = AtomSet.DEFAULT
    pts = [(q, t) for q in (3, 5, 7, 8, 9, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
           for t in (3, 5, 7, 11, 13, 17) if q != t]
    slots = {}
    for si, exps in enumerate([(1, 2, 1, 0, 0, 0), (0, 0, 0, 1, 1, 0), (2, 0, 0, 0, 1, 1), (0, 1, 1, 1, 0, 1)]):
        num = {(i, j): random.randint(-5, 5) for i in range(3) for j in range(2)}
        slots[si] = [(p, Fraction(sum(c * p[0] ** i * p[1] ** j for (i, j), c in num.items()),
                                  six.eval_denom(exps, *p))) for p in pts]
    got, rep = infer(slots, seed=AtomSet.UNITS)
    print(got, "unexplained:", rep["unexplained"])
    assert set(got.atoms) == set(six.atoms), got.names()
    assert rep["unexplained"] == 0
    print("OK: from the units alone, rediscovered the six from denominators")
