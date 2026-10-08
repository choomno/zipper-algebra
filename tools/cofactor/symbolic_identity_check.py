
#!/usr/bin/env python3
"""
Pure-polynomial-identity verifier for the X_3 cofactor at n=4 over the
function field Q(q,t).

Goal
----
Check, IN THE FREE ALGEBRA Q(q,t)<x_1,x_2,x_3,y_1,y_2,y_3>, the identity

    sum_{j=1}^{695}  c_j(q,t) * L_j * R_{i_j} * R'_j  ==  X_3                 (*)

where:
  * The R_i are the 31 input relations produced by
    tools/zipper_in_ms.py :: build_relations(n=4, ...) with q,t
    promoted to formal variables (no specialisation).
  * L_j, R'_j are non-commutative monomials in x_1,x_2,x_3,y_1,y_2,y_3.
  * c_j(q,t) = N_j(q,t)/D_j(q,t) is read literally from
    logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt.
  * X_3 = S[3][1] = (q^2*y_1*y_2 - x_1*x_2) * (q*y_1 + (1-q) - x_1).

The script never specialises q,t to numbers, and uses only exact rational
(Fraction) arithmetic throughout.

Strategy for clearing denominators
----------------------------------
Inspection of the cofactor file shows every D_j is a monomial in the set
    { q, t, (q-1), (t+q^2) }
with degrees q^a*t^b*(q-1)^1*(t+q^2)^1 where 0<=a<=5, 1<=b<=4.  Define
    M = q^5 * t^4 * (q-1) * (t+q^2).
Then (M/D_j) is a polynomial in Z[q,t] for every j.  Multiplying both
sides of (*) by M removes all denominators and reduces (*) to a
polynomial identity in Z[q,t]<x_1,...,y_3>:

    sum_j (M/D_j) * N_j * L_j * R_{i_j} * R'_j  ==  M * X_3.                  (**)

We verify (**) exactly by NC polynomial subtraction.  Since the original
(*) and (**) differ only by the nonzero scalar M(q,t) in the (commutative)
coefficient ring, they are equivalent identities in the free algebra.

NOTE on sympy
-------------
The task description suggested sympy for the (q,t) coefficient ring.
sympy is not installed in this sandbox and pip is unavailable, so we use
Python's stdlib `fractions.Fraction` together with a small bivariate
polynomial class.  All arithmetic is exact; no floating point is used.
This is mathematically equivalent: Z[q,t] <-> ZZ[q,t] in sympy, and the
chosen denominator-clearing makes Q(q,t) unnecessary.
"""

from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py lives
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)


import re
import time
_TOOLS_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _TOOLS_DIR not in sys.path:
    sys.path.insert(0, _TOOLS_DIR)
from fractions import Fraction
from typing import Dict, List, Tuple

ROOT = "/home/choom/braid-research/zipper-algebra"
COF_FILE = os.path.join(
    ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt"
)

# Free-algebra variables.  Order MUST match zipper_in_ms.build_relations:
#   x1, x2, x3, y1, y2, y3   (1-based indices 1..6)
VAR_NAMES = ["x1", "x2", "x3", "y1", "y2", "y3"]
VAR_INDEX = {v: i + 1 for i, v in enumerate(VAR_NAMES)}  # 1-based

N_STRANDS = 4  # so x_1..x_3 and y_1..y_3


# ===========================================================================
# Bivariate polynomial in Q[q,t]:   dict  (deg_q, deg_t) -> Fraction
# ===========================================================================

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
    """Variable by 1-based index in VAR_NAMES."""
    return {(idx,): qt_const(1)}


def nc_var_named(name: str) -> NCPoly:
    return nc_var(VAR_INDEX[name])


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


# ===========================================================================
# Build the 31 input relations symbolically over Q[q,t]
#   (mirrors zipper_in_ms.build_relations with n=4, untw='bigelow')
# ===========================================================================

def build_inputs_and_X3(presentation=None, target=(3, 1)) -> Tuple[List[NCPoly], NCPoly]:
    """The input relations and the target, EXACT over Q[q,t].

    Until 2026-09-17 this function was a 130-line hand copy of zipper_in_ms.build_relations
    for n=4/bigelow (GENERALITY_CENSUS.md §3.1: "a second, divergable copy of the
    presentation"). It now delegates to tools/presentation.py, the single builder, over the
    symbolic coefficient ring. Default: Presentation.DEFAULT and X_3 = S[3][1].
    """
    from presentation import Presentation
    pres = presentation or Presentation(N_STRANDS)
    relations, S = pres.symbolic_ncpolys()
    if target not in S:
        raise ValueError(f"target S{list(target)} not available at n={pres.n}")
    return relations, S[target]


# ===========================================================================
# Parse a numerator polynomial of the cofactor file:
#   "1*q^4*t^0 + -2*q^5*t^0 + 1*q^6*t^0"
# Tokens with an explicit q^a and t^b factor.  Coefficient is a signed
# integer.
# ===========================================================================

_NUM_TERM_RE = re.compile(
    r"\s*(?P<sign>[+\-])?\s*"
    r"(?P<coef>-?\d+)\s*\*\s*q\^(?P<aq>\d+)\s*\*\s*t\^(?P<bt>\d+)\s*"
)


def parse_numerator(s: str) -> QTPoly:
    s = s.strip()
    # strip outer parens if any
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()
    out: QTPoly = {}
    # tokenize by splitting at top-level ' + ' (file has explicit signs)
    # but signs may be '+' or '-'.  We scan with the regex above sequentially.
    pos = 0
    first = True
    while pos < len(s):
        m = _NUM_TERM_RE.match(s, pos)
        if not m:
            # skip a stray whitespace or sign char if present
            if s[pos].isspace() or s[pos] == "+":
                pos += 1
                continue
            raise ValueError(
                f"Cannot parse numerator at position {pos}: "
                f"{s[pos:pos+40]!r} in full string {s!r}"
            )
        sign_g = m.group("sign")
        coef = int(m.group("coef"))
        if sign_g == "-":
            coef = -coef
        aq = int(m.group("aq"))
        bt = int(m.group("bt"))
        s_new = out.get((aq, bt), Fraction(0)) + coef
        if s_new == 0:
            out.pop((aq, bt), None)
        else:
            out[(aq, bt)] = s_new
        pos = m.end()
        first = False
    return out


def parse_denominator(s: str) -> Tuple[int, int, int, int]:
    """Return (a_q, b_t, k_qm1, k_tpq2) for D = q^a * t^b * (q-1)^k_qm1 * (t+q^2)^k_tpq2."""
    s = s.strip()
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()
    a_q = 0
    b_t = 0
    k_qm1 = 0
    k_tpq2 = 0
    parts = [p.strip() for p in s.split("*") if p.strip()]
    for p in parts:
        if p == "q":
            a_q += 1
        elif p == "t":
            b_t += 1
        elif p.startswith("q^"):
            a_q += int(p[2:])
        elif p.startswith("t^"):
            b_t += int(p[2:])
        elif p == "(q-1)":
            k_qm1 += 1
        elif p == "(t+q^2)":
            k_tpq2 += 1
        else:
            raise ValueError(f"Unknown denominator factor: {p!r} in {s!r}")
    return a_q, b_t, k_qm1, k_tpq2


def parse_monomial(s: str) -> NCMon:
    s = s.strip()
    if s == "" or s == "()":
        return ()
    pieces = [p.strip() for p in s.split("*") if p.strip()]
    out: List[int] = []
    for p in pieces:
        if p not in VAR_INDEX:
            raise ValueError(f"Unknown var in monomial: {p!r}")
        out.append(VAR_INDEX[p])
    return tuple(out)


# Line format:
#   <left> * <label> * <right>  ::  <numerator> / <denominator>
# where <label> is either 'R_<idx>' (legacy) or a Z/T label from
# relations_naming.RELATION_LABEL (e.g. 'Z_2', 'T_3', 'INV_1_p').
# left/right may be '()' (empty).
_LINE_RE = re.compile(
    r"^\s*(?P<left>.+?)\s+\*\s+(?P<lbl>[A-Z][A-Za-z0-9_]*)\s+\*\s+(?P<right>.+?)\s*::\s*"
    r"(?P<num>\([^()]*\))\s*/\s*(?P<den>\(.+\))\s*$"
)


def parse_cofactor_line(ln: str):
    from relations_naming import parse_label
    m = _LINE_RE.match(ln)
    if not m:
        raise ValueError(f"Cannot parse cofactor line: {ln!r}")
    left = parse_monomial(m.group("left"))
    idx = parse_label(m.group("lbl"))
    right = parse_monomial(m.group("right"))
    num = parse_numerator(m.group("num"))
    den = parse_denominator(m.group("den"))
    return left, idx, right, num, den


# ===========================================================================
# Main verification.
#
# We clear the common denominator M = q^A * t^B * (q-1)^K1 * (t+q^2)^K2,
# where (A, B, K1, K2) is the coordinatewise max of all (a, b, k1, k2).
# Then we check
#     sum_j (M / D_j) * N_j * L_j * R_{i_j} * R'_j  ==  M * X_3
# as an exact polynomial identity in Z[q,t]<x1..y3>.
# ===========================================================================

def main() -> int:
    t_start = time.time()

    # --- 1. Build inputs and X_3 symbolically over Q[q,t] ---------------
    relations, X3 = build_inputs_and_X3()
    t_built = time.time()
    print(f"[build] {len(relations)} input relations and X_3 built in "
          f"{t_built - t_start:.2f}s.  X_3 has {nc_size(X3)} NC monomials.")
    assert len(relations) == 31, f"Expected 31 inputs, got {len(relations)}"

    # --- 2. Parse the cofactor file --------------------------------------
    with open(COF_FILE) as f:
        raw = f.readlines()
    lines = [ln.rstrip("\n") for ln in raw
             if ln.strip() and not ln.lstrip().startswith("#")]
    print(f"[cof]   {len(lines)} cofactor tuples to parse.")

    parsed = [parse_cofactor_line(ln) for ln in lines]

    # --- 3. Compute LCM denominator (A, B, K1, K2) -----------------------
    A = max(d[4][0] for d in parsed)  # max a_q
    B = max(d[4][1] for d in parsed)  # max b_t
    K1 = max(d[4][2] for d in parsed)  # max k_{q-1}
    K2 = max(d[4][3] for d in parsed)  # max k_{t+q^2}
    print(f"[lcm]   M = q^{A} * t^{B} * (q-1)^{K1} * (t+q^2)^{K2}")

    # Precompute (q-1) and (t+q^2) and their powers, plus (M / D_j) for j.
    qm1 = qt_sub(qt_q(), qt_const(1))
    tpq2 = qt_add(qt_t(), qt_pow(qt_q(), 2))

    # Cache q^k, t^k, (q-1)^k, (t+q^2)^k
    qpow = {0: qt_const(1)}
    for k in range(1, A + 1):
        qpow[k] = qt_mul(qpow[k - 1], qt_q())
    tpow = {0: qt_const(1)}
    for k in range(1, B + 1):
        tpow[k] = qt_mul(tpow[k - 1], qt_t())
    qm1pow = {0: qt_const(1)}
    for k in range(1, K1 + 1):
        qm1pow[k] = qt_mul(qm1pow[k - 1], qm1)
    tpq2pow = {0: qt_const(1)}
    for k in range(1, K2 + 1):
        tpq2pow[k] = qt_mul(tpq2pow[k - 1], tpq2)

    # M = q^A * t^B * (q-1)^K1 * (t+q^2)^K2 as a Q[q,t] polynomial.
    M_poly = qt_mul(qt_mul(qt_mul(qpow[A], tpow[B]), qm1pow[K1]), tpq2pow[K2])
    print(f"[lcm]   M has {len(M_poly)} terms when expanded.")

    # --- 4. Accumulate LHS = sum (M/D_j) * N_j * L_j * R_{i_j} * R'_j ----
    total: NCPoly = nc_zero()
    t_pre = time.time()
    for n_idx, (left, idx, right, num, den) in enumerate(parsed):
        # (M / D_j) factor = q^(A-a) * t^(B-b) * (q-1)^(K1-k1) * (t+q^2)^(K2-k2)
        a_q, b_t, k1, k2 = den
        # Coefficient (as Q[q,t] polynomial): num * (M/D_j)
        cofactor_coef = qt_mul(
            qt_mul(qpow[A - a_q], tpow[B - b_t]),
            qt_mul(qm1pow[K1 - k1], tpq2pow[K2 - k2])
        )
        coef_poly = qt_mul(num, cofactor_coef)
        if not coef_poly:
            continue
        # Build L_j * R_{i_j} * R'_j (3-fold NC product).
        Rij = relations[idx]
        # left, right are bare monomials -> turn into NC polys
        L_nc: NCPoly = {left: qt_const(1)}
        R_nc: NCPoly = {right: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, Rij), R_nc)
        term = nc_scalar(term, coef_poly)
        total = nc_add(total, term)
        if (n_idx + 1) % 100 == 0:
            print(f"[acc]   {n_idx + 1}/{len(parsed)} tuples accumulated "
                  f"({time.time() - t_pre:.1f}s).  total has {nc_size(total)} NC monomials.")

    t_lhs = time.time()
    print(f"[sum]   LHS built in {t_lhs - t_built:.2f}s with "
          f"{nc_size(total)} NC monomials.")

    # --- 5. RHS = M * X_3 ------------------------------------------------
    rhs = nc_scalar(X3, M_poly)
    print(f"[rhs]   M * X_3 has {nc_size(rhs)} NC monomials.")

    # --- 6. Compare ------------------------------------------------------
    diff = nc_sub(total, rhs)
    t_end = time.time()

    if nc_is_zero(diff):
        print("\n" + "=" * 72)
        print("PASS: sum_j c_j(q,t) * L_j * R_{i_j} * R'_j  ==  X_3")
        print(f"      identity verified in the FREE algebra")
        print(f"      Q(q,t)<x_1,x_2,x_3,y_1,y_2,y_3>.")
        print(f"      Tuples used:           {len(parsed)}")
        print(f"      Distinct R_i used:     "
              f"{sorted({d[1] for d in parsed})}")
        print(f"      LCM cleared:           M = q^{A}*t^{B}*(q-1)^{K1}*(t+q^2)^{K2}")
        print(f"      Total wall time:       {t_end - t_start:.2f}s")
        print("=" * 72)
        return 0

    # Non-zero residual: emit informative diagnostics.
    print("\n" + "=" * 72)
    print(f"FAIL: residual has {nc_size(diff)} non-zero NC monomial(s).")
    print("=" * 72)
    items = sorted(diff.items(), key=lambda kv: (len(kv[0]), kv[0]))
    for mon, coef in items[:20]:
        if mon:
            mstr = "*".join(VAR_NAMES[i - 1] for i in mon)
        else:
            mstr = "1"
        print(f"  [{mstr}]  coef = {qt_str(coef)}")
    if len(items) > 20:
        print(f"  ... {len(items) - 20} more non-zero monomials suppressed.")
    print(f"Total wall time: {t_end - t_start:.2f}s")
    return 1


if __name__ == "__main__":
    sys.exit(main())
