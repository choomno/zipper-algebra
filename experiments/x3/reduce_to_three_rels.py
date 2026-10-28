#!/usr/bin/env python3
"""
Reduce the 166-tuple simplified X_3 cofactor at n=4 from 4 input relations
{R_22, R_25, R_26, R_28} to 3 input relations {R_22, R_25, R_28} by using
the conjugacy
    R_26 = w * R_25 * w^{-1}      with  w = sigma_1 * sigma_2 = x_1 * x_2.
This is verified by direct computation:
    w (x_1 - t) w^{-1} = (x_2 - t),   w S[2][1] w^{-1} = S[2][2],
so w * R_25 * w^{-1} = (x_2 - t) * S[2][2] = R_26.

Substitution rule on tuples (L, R_26, R):
    L * R_26 * R   =   (L * w) * R_25 * (w^{-1} * R)
                  =   (L * x_1*x_2) * R_25 * (y_2*y_1 * R).
We replace each R_26 tuple by the corresponding R_25 tuple, then re-canonicalize
left and right multipliers in B_3 (Delta-extraction + left-greedy Garside form),
and merge buckets with the existing R_25 tuples.

Verification: the resulting sum equals X_3 in Q(q,t)[B_3] (group algebra),
checked by Garside-normalization mod B_3 of the residual against X_3.
"""

from __future__ import annotations
import argparse
from experiments.output import prepare_output_dir
import os
from pathlib import Path
import sys

import time
from typing import Dict, Tuple
from collections import Counter
from .relations_naming import label
from .legacy_format import (
    VAR_NAMES,
    qt_zero, qt_const, qt_q, qt_t, qt_add, qt_sub, qt_neg, qt_mul, qt_pow,
    qt_is_zero, qt_str,
    nc_zero, nc_add, nc_sub, nc_mul, nc_scalar, nc_is_zero, nc_size,
    build_inputs_and_X3,
    parse_cofactor_line,
)
from .simplify_cofactor import (
    SIGMA1, SIGMA2, Y1, Y2,
    normalize_word, canonical_word,
)


# w = sigma_1 sigma_2 = x_1 x_2   (free-algebra word)
W_WORD: Tuple[int, ...] = (SIGMA1, SIGMA2)
# w^{-1} = sigma_2^{-1} sigma_1^{-1} = y_2 y_1
W_INV_WORD: Tuple[int, ...] = (Y2, Y1)


def parse_lines(path: str):
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            yield parse_cofactor_line(ln)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="generated certificate from the preceding stage")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    prepare_output_dir(args.output.parent)
    t_start = time.time()
    print("=" * 72)
    print("Reduction: eliminate R_26 via R_26 = (sigma_1 sigma_2) R_25 (sigma_1 sigma_2)^{-1}")
    print("=" * 72)

    relations, X3 = build_inputs_and_X3()
    R25 = relations[25]
    R26 = relations[26]

    # --- 0. Sanity-check the conjugacy in the free algebra mod B_3 ---------
    from .simplify_cofactor import normalize_word
    # Compute w * R_25 * w^{-1} as an NC poly, and compare to R_26 mod B_3.
    L_nc = {W_WORD: qt_const(1)}
    R_nc = {W_INV_WORD: qt_const(1)}
    conj_R25 = nc_mul(nc_mul(L_nc, R25), R_nc)
    # Normalize both mod B_3
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
    conj_R25_B3 = nc_normalize_mod_B3(conj_R25)
    R26_B3 = nc_normalize_mod_B3(R26)
    diff_conj = nc_sub(conj_R25_B3, R26_B3)
    if not nc_is_zero(diff_conj):
        print(f"[FAIL]  conjugacy check failed: residual has {nc_size(diff_conj)} monomials")
        for mon, coef in list(diff_conj.items())[:5]:
            mstr = "*".join(VAR_NAMES[x - 1] for x in mon) if mon else "1"
            print(f"  {mstr}  : {qt_str(coef)}")
        return 1
    print("[check] PASS: w * R_25 * w^{-1} == R_26 in Q(q,t)[B_3].")

    # --- 1. Parse simplified cofactor --------------------------------------
    tuples = list(parse_lines(args.input))
    print(f"[parse] {len(tuples)} input tuples from {os.path.basename(args.input)}.")
    src_hist = Counter(idx for (_, idx, _, _, _) in tuples)
    print(f"[hist]  src R_i histogram: {dict(sorted(src_hist.items()))}")

    # --- 2. Determine global denominator -----------------------------------
    A = max(d[4][0] for d in tuples)
    B = max(d[4][1] for d in tuples)
    K1 = max(d[4][2] for d in tuples)
    K2 = max(d[4][3] for d in tuples)
    print(f"[lcm]   D = q^{A} * t^{B} * (q-1)^{K1} * (t+q^2)^{K2}")
    assert (A, B, K1, K2) == (0, 4, 1, 1), \
        f"expected D = t^4 (q-1)(t+q^2), got exponents {(A,B,K1,K2)}"

    qm1 = qt_sub(qt_q(), qt_const(1))
    tpq2 = qt_add(qt_t(), qt_pow(qt_q(), 2))
    qpow = {0: qt_const(1)}
    for k in range(1, max(A, 6) + 1):
        qpow[k] = qt_mul(qpow[k - 1], qt_q())
    tpow = {0: qt_const(1)}
    for k in range(1, max(B, 6) + 1):
        tpow[k] = qt_mul(tpow[k - 1], qt_t())
    qm1pow = {0: qt_const(1)}
    for k in range(1, K1 + 1):
        qm1pow[k] = qt_mul(qm1pow[k - 1], qm1)
    tpq2pow = {0: qt_const(1)}
    for k in range(1, K2 + 1):
        tpq2pow[k] = qt_mul(tpq2pow[k - 1], tpq2)

    # --- 3. Substitute and rebucket ----------------------------------------
    bucket: Dict[Tuple[Tuple[int, ...], int, Tuple[int, ...]], "QTPoly"] = {}
    for (left, idx, right, num, den) in tuples:
        a_q, b_t, k1, k2 = den
        factor = qt_mul(qt_mul(qpow[A - a_q], tpow[B - b_t]),
                        qt_mul(qm1pow[K1 - k1], tpq2pow[K2 - k2]))
        coef_poly = qt_mul(num, factor)
        if not coef_poly:
            continue
        if idx == 26:
            # (L, 26, R)  ->  (L * w, 25, w^{-1} * R)
            new_left = left + W_WORD
            new_right = W_INV_WORD + right
            new_idx = 25
        else:
            new_left = left
            new_right = right
            new_idx = idx
        # Canonicalize new_left and new_right
        kL, atomsL = normalize_word(new_left)
        kR, atomsR = normalize_word(new_right)
        L_word = canonical_word(kL, atomsL)
        R_word = canonical_word(kR, atomsR)
        key = (L_word, new_idx, R_word)
        if key in bucket:
            bucket[key] = qt_add(bucket[key], coef_poly)
            if qt_is_zero(bucket[key]):
                del bucket[key]
        else:
            bucket[key] = coef_poly

    print(f"[bucket] {len(bucket)} buckets after R_26 -> R_25 substitution.")
    new_hist = Counter(k[1] for k in bucket.keys())
    print(f"[hist]  new R_i histogram: {dict(sorted(new_hist.items()))}")

    # --- 4. Verify the new cofactor equals X_3 mod B_3 ---------------------
    rhs = nc_scalar(X3, qt_mul(qt_mul(qt_mul(qpow[A], tpow[B]), qm1pow[K1]), tpq2pow[K2]))
    total = nc_zero()
    for (L_word, idx, R_word), coef_poly in bucket.items():
        Rij = relations[idx]
        L_nc = {L_word: qt_const(1)}
        R_nc = {R_word: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, Rij), R_nc)
        total = nc_add(total, nc_scalar(term, coef_poly))
    print(f"[verif] |sum| = {nc_size(total)} NC monomials (pre-normalize).")
    total_B3 = nc_normalize_mod_B3(total)
    rhs_B3 = nc_normalize_mod_B3(rhs)
    diff_B3 = nc_sub(total_B3, rhs_B3)
    print(f"[verif] |sum_B3|={nc_size(total_B3)} |rhs_B3|={nc_size(rhs_B3)} |diff|={nc_size(diff_B3)}")
    if not nc_is_zero(diff_B3):
        print("[FAIL]  reduced cofactor does NOT equal X_3 mod B_3. First residuals:")
        for i, (mon, coef) in enumerate(sorted(diff_B3.items(), key=lambda kv: (len(kv[0]), kv[0]))):
            if i >= 8: break
            mstr = "*".join(VAR_NAMES[x - 1] for x in mon) if mon else "1"
            print(f"  {mstr}  : {qt_str(coef)}")
        return 2
    print("[verif] PASS: reduced cofactor sum == X_3 in Q(q,t)[B_3].")

    # --- 5. Write output ---------------------------------------------------
    items = sorted(bucket.items(), key=lambda kv: (kv[0][1], len(kv[0][0]), kv[0][0], len(kv[0][2]), kv[0][2]))

    def mon_to_str(m: Tuple[int, ...]) -> str:
        if not m:
            return "()"
        return "*".join(VAR_NAMES[i - 1] for i in m)

    def num_to_str(p) -> str:
        if not p:
            return "(0)"
        terms = []
        for (a, b), c in sorted(p.items()):
            if c.denominator == 1:
                cs = f"{c.numerator}"
            else:
                cs = f"{c.numerator}/{c.denominator}"
            terms.append(f"{cs}*q^{a}*t^{b}")
        return "(" + " + ".join(terms) + ")"

    parts = []
    if A > 0:
        parts.append(f"q^{A}" if A != 1 else "q")
    if B > 0:
        parts.append(f"t^{B}" if B != 1 else "t")
    if K1 > 0:
        parts.append("(q-1)" if K1 == 1 else "(q-1)" + f"^{K1}")
    if K2 > 0:
        parts.append("(t+q^2)" if K2 == 1 else "(t+q^2)" + f"^{K2}")
    den_str = "*".join(parts) if parts else "1"

    final_hist = Counter(k[1] for k in bucket.keys())
    with open(args.output, "w") as f:
        f.write("# Three-relation Q(q,t) X_3 cofactor at n=4 (Bigelow zipper algebra).\n")
        f.write("# 166-tuple simplified cofactor further reduced to 3 input relations\n")
        f.write("# {Z_2, T_2, T_3} by eliminating T_2_j2 via the conjugacy\n")
        f.write("#   T_2_j2 = (sigma_1 sigma_2) * T_2 * (sigma_1 sigma_2)^{-1}\n")
        f.write("# i.e. each tuple (L, T_2_j2, R) is rewritten as\n")
        f.write("#   (L * x_1*x_2) * T_2 * (y_2*y_1 * R),\n")
        f.write("# then re-canonicalized in B_3 and bucketed.\n")
        f.write("#\n")
        f.write("# Identity holds in the group algebra Q(q,t)[B_3].\n")
        f.write("#\n")
        f.write(f"# Tuple count: {len(items)}  (vs 166 four-relation, vs 695 original)\n")
        labeled_hist = {label(i): c for i, c in sorted(final_hist.items())}
        f.write(f"# Relation histogram: {labeled_hist}\n")
        f.write(f"# Global denominator: D = {den_str}\n")
        f.write("#\n")
        f.write("# IDENTITY (in Q(q,t)[B_3]):\n")
        f.write("#   sum_j  c_j(q,t) * L_j * <label_j> * R'_j  ==  X_3 = S[3][1]\n")
        f.write("# where L_j, R'_j are canonical B_3 words in {x_1,x_2,y_1,y_2}.\n")
        f.write("# Verified by Garside normalization of the free-algebra residual.\n")
        f.write("#\n")
        for (L_word, idx, R_word), coef_poly in items:
            L_s = mon_to_str(L_word)
            R_s = mon_to_str(R_word)
            num_s = num_to_str(coef_poly)
            f.write(f"{L_s} * {label(idx)} * {R_s}  ::  {num_s} / ({den_str})\n")

    print(f"[write] Wrote {args.output}")
    print(f"[stats] Total tuples:        {len(items)}")
    print(f"[stats] vs 166 (4-rel):      {len(items)/166:.2f}x")
    print(f"[stats] vs 695 (original):   {len(items)/695:.2f}x")
    print(f"[stats] By R_i:              {dict(sorted(final_hist.items()))}")
    print(f"[time]  total wall: {time.time()-t_start:.1f}s")

    # Re-parse and re-verify
    print("[verify] Re-parsing output and verifying in Q(q,t)[B_3]...")
    out_tuples = list(parse_lines(args.output))
    print(f"[verify] Re-parsed {len(out_tuples)} tuples.")
    Aout = max(d[4][0] for d in out_tuples)
    Bout = max(d[4][1] for d in out_tuples)
    K1out = max(d[4][2] for d in out_tuples)
    K2out = max(d[4][3] for d in out_tuples)
    assert (Aout, Bout, K1out, K2out) == (A, B, K1, K2)
    total_v = nc_zero()
    for (left, idx, right, num, den) in out_tuples:
        a_q, b_t, k1, k2 = den
        factor = qt_mul(qt_mul(qpow[A - a_q], tpow[B - b_t]),
                        qt_mul(qm1pow[K1 - k1], tpq2pow[K2 - k2]))
        coef_poly = qt_mul(num, factor)
        if not coef_poly:
            continue
        Rij = relations[idx]
        L_nc = {left: qt_const(1)}
        R_nc = {right: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, Rij), R_nc)
        total_v = nc_add(total_v, nc_scalar(term, coef_poly))
    total_v_B3 = nc_normalize_mod_B3(total_v)
    diff_v = nc_sub(total_v_B3, nc_normalize_mod_B3(rhs))
    if nc_is_zero(diff_v):
        print("[verify] PASS: sum == X_3 in Q(q,t)[B_3].")
        return 0
    else:
        print(f"[verify] FAIL: residual {nc_size(diff_v)} monomials.")
        return 3


if __name__ == "__main__":
    sys.exit(main())
