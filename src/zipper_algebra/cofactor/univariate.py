"""Exact sampled reconstruction over rational univariate polynomials.

Coefficient lists are ordered from constant to highest degree. These helpers are
shared by denominator inference and bivariate rational reconstruction; sampled
agreement produces candidates that still need exact certificate verification.
"""
from __future__ import annotations

from fractions import Fraction


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
