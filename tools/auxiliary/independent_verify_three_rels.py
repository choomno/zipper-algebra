#!/usr/bin/env python3
"""
Independent verification of cofactor-X3-n4-Qqt-three-rels.txt.

Parses the 161 data lines from scratch (no use of
symbolic_identity_check.parse_cofactor_line) and verifies that

    sum_j  N_j(q,t) * L_j * relation(label_j) * R_j   ==   D * X_3

as an element of Q(q,t)[B_3], where D = t^4 (q-1)(t+q^2) is the common
denominator.

Uses ONLY the NC polynomial machinery (nc_*, qt_*) plus
build_inputs_and_X3() and the Garside normalization helpers from
simplify_cofactor.py.  Parsing and summation are written inline.
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py lives
COFACTOR = os.path.join(PARENT, "cofactor")  # these are auxiliary to the pipeline, not in it
for _p in (HERE, PARENT, COFACTOR):
    if _p not in sys.path:
        sys.path.insert(0, _p)


import re
import time
from fractions import Fraction
from typing import Dict, List, Tuple

HERE = "/home/choom/braid-research/zipper-algebra/tools"
sys.path.insert(0, HERE)

from symbolic_identity_check import (
    VAR_NAMES, VAR_INDEX,
    qt_const, qt_q, qt_t, qt_add, qt_sub, qt_mul, qt_pow, qt_is_zero, qt_str,
    nc_zero, nc_add, nc_sub, nc_mul, nc_scalar, nc_is_zero, nc_size,
    build_inputs_and_X3,
)
from simplify_cofactor import normalize_word, canonical_word
COF_FILE = "/home/choom/braid-research/zipper-algebra/logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-three-rels.txt"


# ---------------------------------------------------------------------------
# Independent inline parsers (do NOT import parse_cofactor_line).
# ---------------------------------------------------------------------------

LABEL_TO_RELIDX = {"Z_2": 22, "T_2": 25, "T_3": 28}

# Each numerator term: optional leading sign + signed integer + * q^a * t^b.
# The file actually emits tokens like " + -2*q^2*t^2" or " + 1*q^1*t^2"; the
# coefficient already carries its own sign, so the explicit ' + ' is only a
# separator.  We tokenize by splitting on ' + ' boundaries.
def parse_numer_inline(s: str) -> Dict[Tuple[int, int], Fraction]:
    """Parse '(1*q^1*t^2 + -2*q^2*t^2 + ...)' into Q[q,t] dict."""
    s = s.strip()
    assert s.startswith("(") and s.endswith(")"), f"bad numer: {s!r}"
    inner = s[1:-1].strip()
    if inner == "":
        return {}
    # Split on top-level ' + ' (there are no parens inside).
    tokens = inner.split(" + ")
    out: Dict[Tuple[int, int], Fraction] = {}
    term_re = re.compile(r"^\s*(-?\d+)\s*\*\s*q\^(\d+)\s*\*\s*t\^(\d+)\s*$")
    for tok in tokens:
        m = term_re.match(tok)
        if not m:
            raise ValueError(f"cannot parse numer token {tok!r}")
        c = int(m.group(1))
        a = int(m.group(2))
        b = int(m.group(3))
        key = (a, b)
        out[key] = out.get(key, Fraction(0)) + c
        if out[key] == 0:
            out.pop(key, None)
    return out


def parse_denom_inline(s: str) -> Tuple[int, int, int, int]:
    """Parse '(t^4*(q-1)*(t+q^2))' into (a_q, b_t, k_qm1, k_tpq2)."""
    s = s.strip()
    assert s.startswith("(") and s.endswith(")"), f"bad denom: {s!r}"
    inner = s[1:-1].strip()
    # Split on '*' carefully (factors include "(q-1)" and "(t+q^2)").
    parts = []
    depth = 0
    cur = ""
    for ch in inner:
        if ch == "(":
            depth += 1
            cur += ch
        elif ch == ")":
            depth -= 1
            cur += ch
        elif ch == "*" and depth == 0:
            parts.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        parts.append(cur.strip())
    a_q = b_t = k_qm1 = k_tpq2 = 0
    for p in parts:
        if p == "(q-1)":
            k_qm1 += 1
        elif p == "(t+q^2)":
            k_tpq2 += 1
        elif p == "q":
            a_q += 1
        elif p == "t":
            b_t += 1
        elif p.startswith("q^"):
            a_q += int(p[2:])
        elif p.startswith("t^"):
            b_t += int(p[2:])
        else:
            raise ValueError(f"unknown denom factor {p!r}")
    return a_q, b_t, k_qm1, k_tpq2


def parse_monomial_inline(s: str) -> Tuple[int, ...]:
    """Parse '()' or 'x1*x2*y1' into a tuple of 1-based var indices."""
    s = s.strip()
    if s == "()" or s == "":
        return ()
    out = []
    for tok in s.split("*"):
        tok = tok.strip()
        if tok not in VAR_INDEX:
            raise ValueError(f"unknown variable {tok!r}")
        out.append(VAR_INDEX[tok])
    return tuple(out)


def parse_line_inline(ln: str):
    """
    Format:   <L_str> * <label> * <R_str>  ::  <num> / <denom>
    Each side around ' * ' may be '()' or a star-separated var list.
    Returns (L_mon, label_str, R_mon, num_qt, denom_4t).
    """
    # Split off the coefficient with the '::' delimiter.
    left_side, right_side = ln.split("::", 1)
    # Coefficient half: '<num> / <denom>'
    # Both num and denom are parenthesized; split at the first ' / ' between
    # the closing ) of num and opening ( of denom.
    coef = right_side.strip()
    # Find ') / (' as the natural split.
    m = re.match(r"^(\([^()]*\))\s*/\s*(\(.+\))\s*$", coef)
    if not m:
        raise ValueError(f"cannot split coef in {ln!r}")
    num_str = m.group(1)
    den_str = m.group(2)
    num = parse_numer_inline(num_str)
    den = parse_denom_inline(den_str)
    # Word half: '<L_str> * <label> * <R_str>' where label is one of
    # Z_2, T_2, T_3.  L_str or R_str may contain '*' as well.
    # We locate ' * <label> * ' as the unique splitter.
    word_str = left_side.strip()
    found = None
    for lab in LABEL_TO_RELIDX:
        sep = f" * {lab} * "
        idx = word_str.find(sep)
        if idx >= 0:
            L_str = word_str[:idx]
            R_str = word_str[idx + len(sep):]
            found = (L_str, lab, R_str)
            break
    if found is None:
        raise ValueError(f"no label found in {ln!r}")
    L_str, lab, R_str = found
    L_mon = parse_monomial_inline(L_str)
    R_mon = parse_monomial_inline(R_str)
    return L_mon, lab, R_mon, num, den


# ---------------------------------------------------------------------------
# Garside-normalization on NC polynomials (replays
# reduce_to_three_rels.nc_normalize_mod_B3 verbatim; uses simplify_cofactor
# helpers).
# ---------------------------------------------------------------------------

def nc_normalize_mod_B3(p):
    out = {}
    for mon, coef in p.items():
        k, atoms = normalize_word(mon)
        cmon = canonical_word(k, atoms)
        cur = out.get(cmon)
        if cur is None:
            out[cmon] = dict(coef)
        else:
            s = qt_add(cur, coef)
            if qt_is_zero(s):
                out.pop(cmon, None)
            else:
                out[cmon] = s
    return out


# ---------------------------------------------------------------------------
# Main.
# ---------------------------------------------------------------------------

def main() -> int:
    t0 = time.time()

    # 1. Parse all 161 data lines.
    with open(COF_FILE) as f:
        raw = f.readlines()
    data_lines = [ln.rstrip("\n") for ln in raw
                  if ln.strip() and not ln.lstrip().startswith("#")]
    print(f"[parse] {len(data_lines)} non-comment lines")
    if len(data_lines) != 161:
        print(f"[WARN]  expected 161 lines, got {len(data_lines)}")

    parsed = [parse_line_inline(ln) for ln in data_lines]

    # Sanity: only Z_2, T_2, T_3 labels appear.
    label_hist: Dict[str, int] = {}
    for _, lab, _, _, _ in parsed:
        label_hist[lab] = label_hist.get(lab, 0) + 1
    print(f"[hist]  labels: {label_hist}")
    assert set(label_hist).issubset(set(LABEL_TO_RELIDX)), \
        f"unexpected label(s): {set(label_hist) - set(LABEL_TO_RELIDX)}"

    # Sanity: all denominators are exactly t^4 (q-1)(t+q^2) -> (0, 4, 1, 1).
    denom_set = {den for _, _, _, _, den in parsed}
    print(f"[denom] distinct denominators: {denom_set}")
    assert denom_set == {(0, 4, 1, 1)}, \
        f"expected single denom (0,4,1,1) = t^4 (q-1)(t+q^2); got {denom_set}"

    # 2. Build the 31 input relations and X_3.  Pick relations[22], [25], [28].
    relations, X3 = build_inputs_and_X3()
    assert len(relations) == 31, len(relations)
    Z_2 = relations[22]
    T_2 = relations[25]
    T_3 = relations[28]
    relation_by_label = {"Z_2": Z_2, "T_2": T_2, "T_3": T_3}
    print(f"[build] Z_2: {nc_size(Z_2)} mon, T_2: {nc_size(T_2)} mon, T_3: {nc_size(T_3)} mon")
    print(f"[build] X_3: {nc_size(X3)} mon")

    # 3. Accumulate sum_j N_j * L_j * relation(lab_j) * R_j (no /D factor yet
    #    -- since all denominators are equal to D, we just sum numerators).
    total = nc_zero()
    for n_idx, (L_mon, lab, R_mon, num, den) in enumerate(parsed):
        Rij = relation_by_label[lab]
        L_nc = {L_mon: qt_const(1)}
        R_nc = {R_mon: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, Rij), R_nc)
        term = nc_scalar(term, num)
        total = nc_add(total, term)
        if (n_idx + 1) % 50 == 0:
            print(f"[acc]   {n_idx + 1}/{len(parsed)}  total NC monomials: {nc_size(total)}")
    print(f"[sum]   total has {nc_size(total)} NC monomials")

    # 4. Build D * X_3 where D = t^4 (q-1)(t+q^2).
    qm1 = qt_sub(qt_q(), qt_const(1))
    tpq2 = qt_add(qt_t(), qt_pow(qt_q(), 2))
    D = qt_mul(qt_mul(qt_pow(qt_t(), 4), qm1), tpq2)
    print(f"[D]     D = t^4*(q-1)*(t+q^2) has {len(D)} terms when expanded: {qt_str(D)}")
    rhs = nc_scalar(X3, D)
    print(f"[rhs]   D*X_3 has {nc_size(rhs)} NC monomials")

    # 5. Compare in Q(q,t)[B_3] via Garside normalization.
    print("[norm]  normalizing LHS ...")
    total_B3 = nc_normalize_mod_B3(total)
    print(f"[norm]  |total_B3| = {nc_size(total_B3)}")
    print("[norm]  normalizing RHS ...")
    rhs_B3 = nc_normalize_mod_B3(rhs)
    print(f"[norm]  |rhs_B3|   = {nc_size(rhs_B3)}")

    diff = nc_sub(total_B3, rhs_B3)
    t_end = time.time()
    if nc_is_zero(diff):
        print("=" * 72)
        print(f"PASS: residual is zero.  |residual| = 0 NC monomials.")
        print(f"      161 tuples * (Z_2, T_2, T_3) verify D*X_3 in Q(q,t)[B_3].")
        print(f"      wall time: {t_end - t0:.2f}s")
        print("=" * 72)
        return 0
    else:
        print("=" * 72)
        print(f"FAIL: residual has {nc_size(diff)} non-zero NC monomials.")
        items = sorted(diff.items(), key=lambda kv: (len(kv[0]), kv[0]))
        for mon, coef in items[:10]:
            mstr = "*".join(VAR_NAMES[i - 1] for i in mon) if mon else "1"
            print(f"  [{mstr}]  {qt_str(coef)}")
        print("=" * 72)
        return 1


if __name__ == "__main__":
    sys.exit(main())
