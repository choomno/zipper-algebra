#!/usr/bin/env python3
"""Cancel the spurious (q-1) factor from the 161-tuple three-relation cofactor.

Empirical fact: every one of the 161 numerators N_j(q,t) is divisible by (q-1).
The pipeline accumulated (q-1) into the LCM denominator at an earlier stage
but never re-tested whether it survives the final simplifications. It does not:
the cofactor identity holds over the strictly smaller common denominator
    D' := t^4 (t + q^2)
i.e.,
    D' * X_3  =  sum_{j=1}^{161}  M_j'(q,t) * L_j * R_{i_j} * R'_j
with M_j' = M_j / (q-1) in Z[q,t].

This script verifies and produces the minimal-denominator file
    cofactor-X3-n4-Qqt-minimal.txt.

It then RE-VERIFIES the minimal-denominator identity in Q(q,t)[B_3] using the
same Garside-normalization route as simplify_cofactor.py.
"""
from __future__ import annotations
import os
import sys
from fractions import Fraction
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from symbolic_identity_check import (
    parse_cofactor_line, build_inputs_and_X3,
    nc_zero, nc_add, nc_sub, nc_mul, nc_scalar, nc_is_zero, nc_size,
    qt_const, qt_q, qt_t, qt_add, qt_sub, qt_mul, qt_pow, qt_is_zero,
    VAR_NAMES,
)
from simplify_cofactor import normalize_word, canonical_word
from relations_naming import label

ROOT = "/home/choom/gap_scripts/zipper-algebra"
SRC = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-three-rels.txt")
OUT = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-minimal.txt")


def divide_by_qm1(num):
    """Exact polynomial division of num (dict (deg_q, deg_t) -> Fraction) by (q-1).

    Raises ValueError if not exact (i.e. if N(q=1, t) is not identically 0).
    Method: for each t-power b, extract the univariate polynomial p_b(q) =
    sum_a num[(a, b)] * q^a, then divide p_b(q) by (q-1) using synthetic
    division (Horner's method).  p_b(q) / (q-1) = sum_{a>=0} s_a q^a where
    s_a = sum_{i > a} num[(i, b)].
    """
    out = {}
    # Group by t-degree
    by_t = defaultdict(dict)
    max_q = -1
    for (a, b), c in num.items():
        by_t[b][a] = c
        if a > max_q:
            max_q = a
    for b, pb in by_t.items():
        # Check pb(1) = 0
        if sum(pb.values()) != 0:
            raise ValueError(f"Not divisible by (q-1): t-degree {b}, sum = {sum(pb.values())}")
        # synthetic division: p_b(q) / (q-1)
        # If p_b(q) = sum_a c_a q^a, then p_b / (q-1) = sum_{a >= 0} s_a q^a
        # with s_a = sum_{i > a} c_i (Horner from the top down).
        s = Fraction(0)
        coeffs = []
        for a in range(max_q, -1, -1):
            # s_{a-1} = s_a * 1 + c_a   (synthetic division by (q-1), root q=1)
            # Wait: if p(q) = (q-1) Q(q), with Q(q) = sum q_k q^k, then
            # c_a = q_{a-1} - q_a  (for a >= 1)
            # c_0 = -q_0
            # Equivalently, Horner: starting from top, q_{max-1} = c_max,
            # then q_{a-1} = q_a + c_a for a from max-1 down to 1, and
            # the remainder is q_0 + c_0 which must be zero.
            pass
        # Use a different formulation: divide as q^k coefficient extraction
        # p_b(q) = c_max q^max + ... + c_0
        # p_b/(q-1) = q_{max-1} q^{max-1} + ... + q_0
        # with q_{max-1} = c_max, q_{a} = q_{a+1} + c_{a+1} (a from max-2 down)
        # actually the standard formula: if p(q)/(q-1) = sum q_k q^k, then
        # q_k = sum_{i > k} c_i for k = 0..max-1. (Plus remainder c_0 + q_0 = 0.)
        # Let's just do it.
        max_b = max(pb.keys()) if pb else -1
        for k in range(max_b - 1, -1, -1):
            qk = Fraction(0)
            for i in range(k + 1, max_b + 1):
                qk += pb.get(i, Fraction(0))
            if qk != 0:
                out[(k, b)] = out.get((k, b), Fraction(0)) + qk
    return out


def main():
    # 1) Parse the 161-tuple file
    tuples = []
    for ln in open(SRC):
        s = ln.strip()
        if not s or s.startswith('#'): continue
        tuples.append(parse_cofactor_line(ln))
    print(f"[parse] {len(tuples)} tuples from {os.path.basename(SRC)}.")

    # 2) Divide every numerator by (q-1)
    new_tuples = []
    for (left, idx, right, num, den) in tuples:
        a_q, b_t, k1, k2 = den
        # Verify denom has (q-1)^1
        if k1 != 1:
            raise SystemExit(f"Unexpected denom: {den}")
        new_num = divide_by_qm1(num)
        new_den = (a_q, b_t, 0, k2)  # drop (q-1)
        new_tuples.append((left, idx, right, new_num, new_den))
    print(f"[cancel] divided every numerator by (q-1) exactly.")

    # 3) Verify the new cofactor identity in Q(q,t)[B_3]
    relations, X3 = build_inputs_and_X3()
    A, B, K1, K2 = 0, 4, 0, 1   # new (a_q, b_t, k1, k2) for D' = t^4 (t+q^2)
    qm1 = qt_sub(qt_q(), qt_const(1))
    tpq2 = qt_add(qt_t(), qt_pow(qt_q(), 2))
    qpow = {0: qt_const(1)}
    tpow = {0: qt_const(1)}
    for k in range(1, 6): tpow[k] = qt_mul(tpow[k-1], qt_t())
    tpq2pow = {0: qt_const(1), 1: tpq2}

    # rhs = D' * X_3
    D_prime = qt_mul(tpow[B], tpq2pow[K2])
    rhs = nc_scalar(X3, D_prime)

    total = nc_zero()
    for (left, idx, right, num, den) in new_tuples:
        a_q, b_t, k1, k2 = den
        # factor = D' / den_j = t^{B-b_t} * (t+q^2)^{K2-k2}
        factor = qt_mul(tpow[B - b_t], tpq2pow[K2 - k2])
        coef = qt_mul(num, factor)
        Rij = relations[idx]
        L_nc = {left: qt_const(1)}
        R_nc = {right: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, Rij), R_nc)
        total = nc_add(total, nc_scalar(term, coef))

    def normalize_modB3(p):
        out = {}
        for mon, coef in p.items():
            k, atoms = normalize_word(mon)
            cmon = canonical_word(k, atoms)
            cur = out.get(cmon)
            if cur is None:
                out[cmon] = dict(coef)
            else:
                s = qt_add(cur, coef)
                if qt_is_zero(s): out.pop(cmon, None)
                else: out[cmon] = s
        return out

    total_B3 = normalize_modB3(total)
    rhs_B3 = normalize_modB3(rhs)
    diff = nc_sub(total_B3, rhs_B3)
    print(f"[verif] |sum_B3| = {nc_size(total_B3)}  |D'*X_3_B3| = {nc_size(rhs_B3)}  |diff| = {nc_size(diff)}")
    if not nc_is_zero(diff):
        print("[FAIL] residual is not zero")
        return 1
    print("[verif] PASS: D' * X_3 = sum_j M_j' * L_j * R_{i_j} * R'_j in Q(q,t)[B_3].")
    print(f"        with D' = t^4 (t+q^2), the strictly smaller common denominator.")

    # 4) Write output file
    histogram = Counter(t[1] for t in new_tuples)
    items = sorted(new_tuples, key=lambda x: (x[1], len(x[0]), x[0], len(x[2]), x[2]))

    def mon_to_str(m):
        if not m: return "()"
        return "*".join(VAR_NAMES[i - 1] for i in m)

    def num_to_str(p):
        if not p: return "(0)"
        terms = []
        for (a, b), c in sorted(p.items()):
            cs = f"{c.numerator}" if c.denominator == 1 else f"{c.numerator}/{c.denominator}"
            terms.append(f"{cs}*q^{a}*t^{b}")
        return "(" + " + ".join(terms) + ")"

    den_str = "t^4*(t+q^2)"

    with open(OUT, "w") as f:
        f.write("# Minimal-denominator three-relation Q(q,t) X_3 cofactor at n=4.\n")
        f.write("# Obtained from cofactor-X3-n4-Qqt-three-rels.txt by canceling\n")
        f.write("# the spurious (q-1) factor from every numerator and from the\n")
        f.write("# common denominator. Every numerator in the source file is\n")
        f.write("# divisible by (q-1) exactly; the cancellation is exact in Z[q,t].\n")
        f.write("#\n")
        f.write("# Relation labels:\n")
        f.write("#   Z_2 = level-2 zipper      (q*y_2 + (1-q) - x_2 - q*y_1*y_2 + x_1*x_2) * S[2][1]\n")
        f.write("#   T_2 = level-2 untwisting  (x_1 - t) * S[2][1]\n")
        f.write("#   T_3 = level-3 untwisting  (x_2*x_1 - q^2*y_2*y_1) * S[3][1]\n")
        f.write("#\n")
        f.write(f"# Tuple count: {len(items)}  (same as three-rels.txt; only denominator changed)\n")
        labeled_hist = {label(i): c for i, c in sorted(histogram.items())}
        f.write(f"# Relation histogram: {labeled_hist}\n")
        f.write(f"# Global denominator: D' = {den_str}\n")
        f.write("#\n")
        f.write("# IDENTITY (in Q(q,t)[B_3]):\n")
        f.write("#   sum_j  c_j(q,t) * L_j * <label_j> * R'_j  ==  X_3 = S[3][1]\n")
        f.write("# with c_j(q,t) = M_j'(q,t) / D', M_j' in Z[q,t].\n")
        f.write("# Verified by Garside normalization of the free-algebra residual.\n")
        f.write("#\n")
        for (left, idx, right, num, den) in items:
            f.write(f"{mon_to_str(left)} * {label(idx)} * {mon_to_str(right)}  ::  "
                    f"{num_to_str(num)} / ({den_str})\n")

    print(f"[write] Wrote {OUT}")
    print(f"[stats] Total tuples: {len(items)}")
    print(f"[stats] Histogram: {labeled_hist}")
    print(f"[stats] Global denominator: D' = {den_str}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
