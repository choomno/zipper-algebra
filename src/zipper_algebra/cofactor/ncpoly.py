#!/usr/bin/env python3
"""Exact arithmetic for free-algebra polynomials with coefficients in Q[q,t].

    QTPoly = {(deg_q, deg_t): Fraction}                 an element of Q[q,t]
    NCPoly = {word: QTPoly}, word = tuple of 1-based variable indices

This independent sparse-polynomial implementation supplies the exact verifier's
accumulation arithmetic. It does not construct presentations or choose targets.
"""
from __future__ import annotations

from fractions import Fraction
from typing import Dict, Tuple

QTPoly = Dict[Tuple[int, int], Fraction]


def qt_zero() -> QTPoly:
    return {}


def qt_const(c) -> QTPoly:
    f = Fraction(c)
    if f == 0:
        return {}
    return {(0, 0): f}


def qt_q() -> QTPoly:
    return {(1, 0): Fraction(1)}


def qt_t() -> QTPoly:
    return {(0, 1): Fraction(1)}


def qt_monomial(a: int, b: int, c: Fraction) -> QTPoly:
    if c == 0:
        return {}
    return {(a, b): Fraction(c)}


def qt_add(p: QTPoly, q: QTPoly) -> QTPoly:
    out = dict(p)
    for k, v in q.items():
        s = out.get(k, Fraction(0)) + v
        if s == 0:
            out.pop(k, None)
        else:
            out[k] = s
    return out


def qt_sub(p: QTPoly, q: QTPoly) -> QTPoly:
    out = dict(p)
    for k, v in q.items():
        s = out.get(k, Fraction(0)) - v
        if s == 0:
            out.pop(k, None)
        else:
            out[k] = s
    return out


def qt_neg(p: QTPoly) -> QTPoly:
    return {k: -v for k, v in p.items()}


def qt_mul(p: QTPoly, q: QTPoly) -> QTPoly:
    out: QTPoly = {}
    for (a1, b1), c1 in p.items():
        for (a2, b2), c2 in q.items():
            k = (a1 + a2, b1 + b2)
            s = out.get(k, Fraction(0)) + c1 * c2
            if s == 0:
                out.pop(k, None)
            else:
                out[k] = s
    return out


def qt_pow(p: QTPoly, n: int) -> QTPoly:
    r = qt_const(1)
    for _ in range(n):
        r = qt_mul(r, p)
    return r


def qt_eq(p: QTPoly, q: QTPoly) -> bool:
    return p == q


def qt_is_zero(p: QTPoly) -> bool:
    return len(p) == 0


def qt_str(p: QTPoly) -> str:
    if not p:
        return "0"
    items = sorted(p.items())
    parts = []
    for (a, b), c in items:
        cs = (str(c.numerator) if c.denominator == 1
              else f"{c.numerator}/{c.denominator}")
        ms = []
        if a:
            ms.append(f"q^{a}")
        if b:
            ms.append(f"t^{b}")
        mstr = "*".join(ms) if ms else ""
        if mstr == "":
            parts.append(cs)
        elif cs == "1":
            parts.append(mstr)
        elif cs == "-1":
            parts.append("-" + mstr)
        else:
            parts.append(f"{cs}*{mstr}")
    s = " + ".join(parts).replace("+ -", "- ")
    return s


# ===========================================================================
# Non-commutative polynomial in (Q[q,t]) < x_1, x_2, x_3, y_1, y_2, y_3 >
#   represented as dict (mon_tuple of 1-based var indices) -> QTPoly
# ===========================================================================

NCMon = Tuple[int, ...]
NCPoly = Dict[NCMon, QTPoly]


def nc_zero() -> NCPoly:
    return {}


def nc_const(c: QTPoly) -> NCPoly:
    if not c:
        return {}
    return {(): dict(c)}


def nc_var(idx: int) -> NCPoly:
    """The variable with this 1-based index (the index is into the presentation's var_names)."""
    return {(idx,): qt_const(1)}


def nc_var_named(name: str, var_index: dict) -> NCPoly:
    """The variable by name; var_index maps a name to its 1-based index (no module globals:
    which variables exist is a property of the presentation, not of the arithmetic)."""
    return nc_var(var_index[name])


def _set_or_drop(d: NCPoly, m: NCMon, p: QTPoly) -> None:
    if qt_is_zero(p):
        d.pop(m, None)
    else:
        d[m] = p


def nc_add(a: NCPoly, b: NCPoly) -> NCPoly:
    out = {m: dict(p) for m, p in a.items()}
    for m, p in b.items():
        cur = out.get(m)
        if cur is None:
            out[m] = dict(p)
        else:
            s = qt_add(cur, p)
            _set_or_drop(out, m, s)
    return out


def nc_sub(a: NCPoly, b: NCPoly) -> NCPoly:
    out = {m: dict(p) for m, p in a.items()}
    for m, p in b.items():
        cur = out.get(m)
        if cur is None:
            out[m] = qt_neg(p)
        else:
            s = qt_sub(cur, p)
            _set_or_drop(out, m, s)
    return out


def nc_neg(a: NCPoly) -> NCPoly:
    return {m: qt_neg(p) for m, p in a.items()}


def nc_mul(a: NCPoly, b: NCPoly) -> NCPoly:
    out: NCPoly = {}
    for m1, c1 in a.items():
        for m2, c2 in b.items():
            mk = m1 + m2
            prod = qt_mul(c1, c2)
            if not prod:
                continue
            cur = out.get(mk)
            if cur is None:
                out[mk] = prod
            else:
                s = qt_add(cur, prod)
                _set_or_drop(out, mk, s)
    return out


def nc_scalar(a: NCPoly, c: QTPoly) -> NCPoly:
    if not c:
        return {}
    out: NCPoly = {}
    for m, p in a.items():
        prod = qt_mul(p, c)
        if prod:
            out[m] = prod
    return out


def nc_is_zero(a: NCPoly) -> bool:
    return len(a) == 0


def nc_size(a: NCPoly) -> int:
    return len(a)
