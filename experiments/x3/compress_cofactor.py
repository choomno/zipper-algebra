#!/usr/bin/env python3
"""
Compress the 695-tuple X_3 cofactor onto the chosen relation set
{R_0, R_1, R_2, R_3, R_6, R_22, R_25, R_28} for the n=4 Bigelow zipper algebra.

Strategy:
  1. For each dependent input R_d in {R_7, R_14, R_15, R_16, R_17, R_26},
     express R_d as a sum sum_k a_k * L_k * R_{i_k} * R'_k over the chosen
     generators. Verify each derivation as a polynomial identity in the
     free algebra Q(q,t)<x_1,x_2,x_3,y_1,y_2,y_3>.
  2. Substitute these into the 695-tuple cofactor.
  3. Re-verify the compressed cofactor (sum equals X_3 exactly).
  4. Write out cofactor-X3-n4-Qqt-compressed.txt.
"""
from __future__ import annotations
import argparse
from experiments.output import prepare_output_dir, capture_experiment_provenance
from zipper_algebra.provenance import write_output_provenance
import os
from pathlib import Path
import sys

import re
import time
from fractions import Fraction
from typing import Dict, List, Tuple
from .relations_naming import label
from .legacy_format import (
    VAR_NAMES, VAR_INDEX, N_STRANDS,
    qt_zero, qt_const, qt_q, qt_t, qt_add, qt_sub, qt_neg, qt_mul, qt_pow,
    qt_eq, qt_is_zero, qt_str,
    nc_zero, nc_const, nc_var, nc_var_named, nc_add, nc_sub, nc_neg,
    nc_mul, nc_scalar, nc_is_zero, nc_size,
    build_inputs_and_X3,
    parse_cofactor_line, read_cofactor_input, parse_numerator, parse_denominator,
    COF_FILE, ROOT,
)


# --------- Convenience constructors -----------------------------------
def mon_nc(mon: Tuple[int, ...]) -> "NCPoly":
    return {mon: qt_const(1)}


def mon_from_str(s: str) -> Tuple[int, ...]:
    if not s or s == "()":
        return ()
    return tuple(VAR_INDEX[v] for v in s.split("*"))


def mon_to_str(m: Tuple[int, ...]) -> str:
    if not m:
        return "()"
    return "*".join(VAR_NAMES[i - 1] for i in m)


# --------- Derivation table ------------------------------------------
# Each dependent relation is represented as a list of (coef_qt, L_mon,
# R_idx, R'_mon) tuples, where the sum equals the dependent relation
# exactly in the free algebra. Verified at runtime.

# Helpers
def Q(c) -> "QTPoly":
    return qt_const(c)

# (q-1)
def qm1() -> "QTPoly":
    return qt_sub(qt_q(), Q(1))


def derivations(relations):
    """Explicit free-algebra identities for the inverse braid and mixed relations.

    Each tuple is (coefficient, left word, defining-relation index, right word).
    The caller checks every identity against the ordered defining presentation
    before substitution. R_26 remains a defining relation at this stage.
    """
    KEEP = {0, 1, 2, 3, 6, 22, 25, 28}

    derivs: Dict[int, List[Tuple["QTPoly", Tuple[int, ...], int, Tuple[int, ...]]]] = {}

    # Inverse braid relation.
    derivs[7] = [
        (Q(-1), mon_from_str("y2*y1*y2"), 6, mon_from_str("y1*y2*y1")),
        (Q(1),  mon_from_str("y2*y1*y2"), 0, ()),
        (Q(1),  mon_from_str("y2*y1*y2*x1"), 2, mon_from_str("y1")),
        (Q(1),  mon_from_str("y2*y1*y2*x1*x2"), 0, mon_from_str("y2*y1")),
        (Q(-1), (), 3, mon_from_str("y1*y2*y1")),
        (Q(-1), mon_from_str("y2"), 1, mon_from_str("x2*y1*y2*y1")),
        (Q(-1), mon_from_str("y2*y1"), 3, mon_from_str("x1*x2*y1*y2*y1")),
    ]

    # Mixed relation MIX_1_a.
    derivs[14] = [
        (Q(-1), mon_from_str("y2"), 6, mon_from_str("y1")),
        (Q(1),  mon_from_str("y2*x1*x2"), 0, ()),
        (Q(-1), (), 3, mon_from_str("x1*x2*y1")),
    ]

    # Mixed relation MIX_1_b.
    derivs[15] = [
        (Q(1),  mon_from_str("y1*y2"), 6, mon_from_str("y1*y2")),
        (Q(-1), mon_from_str("y1*y2*x1"), 2, ()),
        (Q(-1), mon_from_str("y1*y2*x1*x2"), 0, mon_from_str("y2")),
        (Q(1),  (), 1, mon_from_str("x2*y1*y2")),
        (Q(1),  mon_from_str("y1"), 3, mon_from_str("x1*x2*y1*y2")),
    ]

    # Mixed relation MIX_1_c.
    derivs[16] = [
        (Q(1),  mon_from_str("y2*y1"), 6, mon_from_str("y2*y1")),
        (Q(-1), (), 3, mon_from_str("x1*y2*y1")),
        (Q(-1), mon_from_str("y2"), 1, mon_from_str("x2*x1*y2*y1")),
        (Q(1),  mon_from_str("y2*y1*x2"), 0, ()),
        (Q(1),  mon_from_str("y2*y1*x2*x1"), 2, mon_from_str("y1")),
    ]

    # Mixed relation MIX_1_d.
    derivs[17] = [
        (Q(-1), mon_from_str("y1"), 6, mon_from_str("y2")),
        (Q(1),  (), 1, mon_from_str("x2*x1*y2")),
        (Q(-1), mon_from_str("y1*x2*x1"), 2, ()),
    ]

    # Retained adjacent untwisting relation.
    derivs[26] = None

    return derivs, KEEP


# --------- Verify a derivation ----------------------------------------
def eval_derivation(deriv, relations) -> "NCPoly":
    """Evaluate sum_k a_k * L_k * R_{i_k} * R'_k as NCPoly."""
    out = nc_zero()
    for (coef, L_mon, idx, R_mon) in deriv:
        L_nc = {L_mon: qt_const(1)}
        R_nc = {R_mon: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, relations[idx]), R_nc)
        term = nc_scalar(term, coef)
        out = nc_add(out, term)
    return out




# --------- Main compression -------------------------------------------
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=COF_FILE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    prepare_output_dir(args.output.parent)
    producer = capture_experiment_provenance(__file__)
    tuples, input_record = read_cofactor_input(args.input)
    inputs = [input_record]
    t_start = time.time()
    relations, X3 = build_inputs_and_X3()
    print(f"[build] {len(relations)} relations, X_3 with {nc_size(X3)} mons")

    derivs, KEEP = derivations(relations)

    # Verify each derivation
    print("\n[verify derivations]")
    for d_idx in [7, 14, 15, 16, 17]:
        deriv = derivs[d_idx]
        if deriv is None:
            print(f"  R_{d_idx}: SKIPPED (no derivation provided)")
            continue
        lhs = eval_derivation(deriv, relations)
        rhs = relations[d_idx]
        diff = nc_sub(lhs, rhs)
        ok = nc_is_zero(diff)
        # Check all source indices are in KEEP
        src_idxs = sorted({t[2] for t in deriv})
        in_keep = all(i in KEEP for i in src_idxs)
        print(f"  R_{d_idx}: {'OK' if ok else 'FAIL'}  uses {src_idxs} "
              f"({'all in KEEP' if in_keep else 'MISSING'})  {len(deriv)} terms")
        if not ok:
            # print first few residual monomials for debug
            items = sorted(diff.items(), key=lambda kv: (len(kv[0]), kv[0]))
            for mon, coef in items[:5]:
                mstr = "*".join(VAR_NAMES[i - 1] for i in mon) if mon else "1"
                print(f"    residual [{mstr}] = {qt_str(coef)}")
            return 1

    print("  R_26: retained (no derivation)")
    KEEP.add(26)

    parsed = tuples
    print(f"\n[cof] {len(parsed)} input tuples")

    # --- Substitute dependents
    # Each tuple is (L, idx, R, num, den)  with num/den in QTPoly form.
    # If idx is a dependent, replace by the derivation, multiplying out:
    #   c * L * R_idx * R'  ->  sum_k (c * a_k) * (L * L_k) * R_{i_k} * (R'_k * R')
    #
    # We need to convert (num, den) into a single coefficient over Q[q,t] by
    # multiplying through by den (which is a monomial in q, t, (q-1), (t+q^2)).
    #
    # Actually it's easier to keep the (num, den) representation. Denominators
    # are MONOMIALS in {q, t, (q-1), (t+q^2)} so we can track them as 4-tuples.
    # The compressed cofactor will have the same denominator structure (since
    # we're adding coefficients of the original tuples + new ones that
    # introduce no new denominators).

    # Compressed tuple format: {(L, idx, R) -> (num_dict, (a_q, b_t, k1, k2))}
    # where num_dict: QTPoly (in Q[q,t]).

    # For derivations: a_k is in Q[q,t] (constants), so multiplying by a_k just
    # multiplies the numerator. Denominator stays the same.

    # Step 1: convert (num, den) to a uniform denominator representation.
    # Actually, let me NOT unify the LCM. Instead, for each output (L, idx, R),
    # accumulate a list of (num, den) pairs and at the end combine to common
    # denominator.

    # Simpler approach: represent each tuple's coefficient as (num: QTPoly, den_tup),
    # and for each output (L, idx, R), keep a list. At end, find LCM of dens
    # and combine.

    from collections import defaultdict
    # accumulator: (L, idx, R) -> list of (num: QTPoly, den_tup)
    acc: Dict[Tuple[Tuple[int, ...], int, Tuple[int, ...]],
              List[Tuple["QTPoly", Tuple[int, int, int, int]]]] = defaultdict(list)

    for (L, idx, R, num, den) in parsed:
        deriv = derivs.get(idx)
        if deriv is None:
            # idx is in KEEP (or treated as kept). Pass through.
            acc[(L, idx, R)].append((num, den))
        else:
            # idx is dependent. Substitute.
            for (a_k, L_k, src_idx, R_k) in deriv:
                new_L = L + L_k
                new_R = R_k + R
                # New num = num * a_k
                new_num = qt_mul(num, a_k)
                if not new_num:
                    continue
                acc[(new_L, src_idx, new_R)].append((new_num, den))

    # Step 2: combine list of (num, den) pairs into a single (num, den) per key,
    # by finding the LCM of all dens. Drop zero results.
    def lcm_den(dens):
        return tuple(max(d[i] for d in dens) for i in range(4))

    def den_str(d):
        a, b, k1, k2 = d
        parts = []
        if a == 1: parts.append("q")
        elif a > 1: parts.append(f"q^{a}")
        if b == 1: parts.append("t")
        elif b > 1: parts.append(f"t^{b}")
        for _ in range(k1): parts.append("(q-1)")
        for _ in range(k2): parts.append("(t+q^2)")
        if not parts:
            return "(1)"
        return "(" + "*".join(parts) + ")"

    def num_str(n):
        # Render QTPoly in the file format: "(1*q^a*t^b + -2*q^c*t^d + ...)"
        if not n:
            return "(0)"
        items = sorted(n.items(), key=lambda kv: kv[0])
        parts = []
        first = True
        for (a, b), c in items:
            assert c.denominator == 1, f"Non-integer numerator coefficient: {c}"
            cs = str(c.numerator)
            term = f"{cs}*q^{a}*t^{b}"
            if first:
                parts.append(term)
                first = False
            else:
                if cs.startswith("-"):
                    parts.append(" + " + term)
                else:
                    parts.append(" + " + term)
        return "(" + "".join(parts) + ")"

    compressed: List[Tuple[Tuple[int, ...], int, Tuple[int, ...], "QTPoly",
                          Tuple[int, int, int, int]]] = []

    for key, entries in acc.items():
        # Each entry is (num, den). Combine with LCM denominator.
        if len(entries) == 1:
            num, den = entries[0]
        else:
            tgt_den = lcm_den([e[1] for e in entries])
            total_num: "QTPoly" = {}
            for num, den in entries:
                # multiplier = q^(tgt-d) * t^(...) * (q-1)^(...) * (t+q^2)^(...)
                a_diff = tgt_den[0] - den[0]
                b_diff = tgt_den[1] - den[1]
                k1_diff = tgt_den[2] - den[2]
                k2_diff = tgt_den[3] - den[3]
                mult = qt_const(1)
                if a_diff: mult = qt_mul(mult, qt_pow(qt_q(), a_diff))
                if b_diff: mult = qt_mul(mult, qt_pow(qt_t(), b_diff))
                if k1_diff:
                    qm1_p = qt_sub(qt_q(), qt_const(1))
                    mult = qt_mul(mult, qt_pow(qm1_p, k1_diff))
                if k2_diff:
                    tpq2_p = qt_add(qt_t(), qt_pow(qt_q(), 2))
                    mult = qt_mul(mult, qt_pow(tpq2_p, k2_diff))
                contrib = qt_mul(num, mult)
                total_num = qt_add(total_num, contrib)
            num = total_num
            den = tgt_den

        # Simplify: cancel common q, t, (q-1), (t+q^2) factors.
        # For q and t factors: find min q-degree and min t-degree across num.
        if num:
            min_q = min(k[0] for k in num)
            min_t = min(k[1] for k in num)
            # Don't strip below 0
            strip_q = min(min_q, den[0])
            strip_t = min(min_t, den[1])
            if strip_q > 0 or strip_t > 0:
                new_num: "QTPoly" = {}
                for (a, b), c in num.items():
                    new_num[(a - strip_q, b - strip_t)] = c
                num = new_num
                den = (den[0] - strip_q, den[1] - strip_t, den[2], den[3])

        if not num:
            continue
        compressed.append((key[0], key[1], key[2], num, den))

    # --- Verify ---
    print(f"\n[compressed] {len(compressed)} tuples")
    distinct_idxs = sorted({c[1] for c in compressed})
    print(f"[compressed] distinct R_i: {distinct_idxs}")

    # Build LHS sum
    print("\n[verify compressed] computing LHS sum...")
    total: "NCPoly" = nc_zero()
    # Compute global LCM denominator
    A = max(c[4][0] for c in compressed)
    B = max(c[4][1] for c in compressed)
    K1 = max(c[4][2] for c in compressed)
    K2 = max(c[4][3] for c in compressed)
    print(f"[lcm] M = q^{A}*t^{B}*(q-1)^{K1}*(t+q^2)^{K2}")

    qm1_poly = qt_sub(qt_q(), qt_const(1))
    tpq2_poly = qt_add(qt_t(), qt_pow(qt_q(), 2))
    M_poly = qt_const(1)
    M_poly = qt_mul(M_poly, qt_pow(qt_q(), A))
    M_poly = qt_mul(M_poly, qt_pow(qt_t(), B))
    M_poly = qt_mul(M_poly, qt_pow(qm1_poly, K1))
    M_poly = qt_mul(M_poly, qt_pow(tpq2_poly, K2))

    for k, (L, idx, R, num, den) in enumerate(compressed):
        a_q, b_t, k1, k2 = den
        cofactor_mult = qt_const(1)
        cofactor_mult = qt_mul(cofactor_mult, qt_pow(qt_q(), A - a_q))
        cofactor_mult = qt_mul(cofactor_mult, qt_pow(qt_t(), B - b_t))
        cofactor_mult = qt_mul(cofactor_mult, qt_pow(qm1_poly, K1 - k1))
        cofactor_mult = qt_mul(cofactor_mult, qt_pow(tpq2_poly, K2 - k2))
        cp = qt_mul(num, cofactor_mult)
        if not cp:
            continue
        Rij = relations[idx]
        L_nc = {L: qt_const(1)}
        R_nc = {R: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, Rij), R_nc)
        term = nc_scalar(term, cp)
        total = nc_add(total, term)
        if (k + 1) % 200 == 0:
            print(f"  accumulated {k+1}/{len(compressed)} tuples (total NC mons={nc_size(total)})")

    rhs = nc_scalar(X3, M_poly)
    diff = nc_sub(total, rhs)

    if not nc_is_zero(diff):
        print("\nFAIL: compressed cofactor does not equal M*X_3.")
        items = sorted(diff.items(), key=lambda kv: (len(kv[0]), kv[0]))
        for mon, coef in items[:10]:
            mstr = "*".join(VAR_NAMES[i - 1] for i in mon) if mon else "1"
            print(f"  [{mstr}] = {qt_str(coef)}")
        return 1

    print("\nPASS: compressed cofactor verified.")

    # --- Write output ---
    OUT = args.output
    # Sort: by idx, then by L mon, then by R mon
    compressed.sort(key=lambda c: (c[1], c[0], c[2]))
    with open(OUT, "w") as f:
        f.write("# Compressed symbolic Q(q,t) cofactor of X_3 over n=4 Bigelow zipper.\n")
        f.write(f"# Original: 695 tuples using indices {{0,1,2,3,6,7,14,15,16,17,22,25,26,28}}.\n")
        f.write(f"# Compressed: {len(compressed)} tuples using indices {distinct_idxs}.\n")
        f.write("# Dependent inputs eliminated via derivations:\n")
        f.write("#   BR_1_y = y_1 y_2 y_1 - y_2 y_1 y_2 derived from BR_1 + INV_{1,2}_{p,m}\n")
        f.write("#   MIX_1_a = x_1 x_2 y_1 - y_2 x_1 x_2 derived from BR_1, INV_1_p, INV_2_m\n")
        f.write("#   MIX_1_b = -x_2 y_1 y_2 + y_1 y_2 x_1 derived from BR_1 + INV_{1,2}_{p,m}\n")
        f.write("#   MIX_1_c = x_1 y_2 y_1 - y_2 y_1 x_2 derived from BR_1 + INV_{1,2}_{p,m}\n")
        f.write("#   MIX_1_d = -x_2 x_1 y_2 + y_1 x_2 x_1 derived from BR_1, INV_1_m, INV_2_p\n")
        if 26 in KEEP:
            f.write("#   T_2_j2 KEPT (not eliminated; derivation by braid conjugation not implemented)\n")
        else:
            f.write("#   T_2_j2 derived from T_2, BR_1 + INV_{1,2}_{p,m}\n")
        f.write("#\n")
        for L, idx, R, num, den in compressed:
            f.write(f"{mon_to_str(L)} * {label(idx)} * {mon_to_str(R)}  ::  "
                    f"{num_str(num)} / {den_str(den)}\n")

    write_output_provenance(args.output, producer=producer, inputs=inputs)
    print(f"\n[output] wrote {OUT}")
    print(f"  tuples:               {len(compressed)} (was 695)")
    print(f"  distinct R_i used:    {distinct_idxs} ({len(distinct_idxs)} total)")
    print(f"  denominator atoms:    q^{A} * t^{B} * (q-1)^{K1} * (t+q^2)^{K2}")
    print(f"  total wall time:      {time.time() - t_start:.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
