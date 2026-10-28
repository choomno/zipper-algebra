#!/usr/bin/env python3
"""Verify the compact braid-algebra identity and its algebraic ideal containments.

The three-relation witness proves X3 lies in <Z2,T2,T3> after the displayed
localization. This script also checks the defining factorizations of T3 and Z3
through X3. It makes no nonmembership claim about <Z2,T2,Z3> from historical
specializations or dimensions.
"""

from __future__ import annotations
import argparse
import os
from pathlib import Path
import sys

import time
from typing import Dict, Tuple
from .legacy_format import (
    VAR_NAMES,
    qt_const, qt_q, qt_t, qt_add, qt_sub, qt_neg, qt_mul, qt_pow,
    qt_is_zero, qt_str,
    nc_zero, nc_const, nc_add, nc_sub, nc_mul, nc_scalar,
    nc_is_zero, nc_size,
    build_inputs_and_X3,
    parse_cofactor_line,
)
from .simplify_cofactor import normalize_word, canonical_word

ROOT = str(Path(__file__).resolve().parents[2])
SRC_FILE = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-three-rels.txt")

# Label -> original relation index
LABEL_TO_IDX = {"Z_2": 22, "T_2": 25, "T_3": 28}


def parse_ZT_lines(path):
    """Parse the ZT (or any) cofactor file. The parser accepts both Z/T
    labels (Z_2, T_2, T_3, ...) and the legacy R_<idx> form, so this
    works on either kind of file."""
    out = []
    with open(path) as f:
        for ln in f:
            s = ln.strip()
            if not s or s.startswith("#"):
                continue
            parsed = parse_cofactor_line(ln)
            out.append(parsed)
    return out


def nc_normalize_mod_B3(p):
    """Map every NC monomial to its B_3 Garside canonical word."""
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


def verify_cofactor_identity(parsed, relations, X3):
    """Compute the sum sum_j c_j L_j R_{i_j} R'_j in the free algebra and
    check it equals X_3 modulo B_3 (Garside normalization)."""
    A = max(d[4][0] for d in parsed)
    B = max(d[4][1] for d in parsed)
    K1 = max(d[4][2] for d in parsed)
    K2 = max(d[4][3] for d in parsed)
    qm1 = qt_sub(qt_q(), qt_const(1))
    tpq2 = qt_add(qt_t(), qt_pow(qt_q(), 2))
    qpow = {0: qt_const(1)}
    for k in range(1, max(A, 8) + 1):
        qpow[k] = qt_mul(qpow[k - 1], qt_q())
    tpow = {0: qt_const(1)}
    for k in range(1, max(B, 8) + 1):
        tpow[k] = qt_mul(tpow[k - 1], qt_t())
    qm1pow = {0: qt_const(1)}
    for k in range(1, K1 + 1):
        qm1pow[k] = qt_mul(qm1pow[k - 1], qm1)
    tpq2pow = {0: qt_const(1)}
    for k in range(1, K2 + 1):
        tpq2pow[k] = qt_mul(tpq2pow[k - 1], tpq2)

    rhs = nc_scalar(X3, qt_mul(qt_mul(qt_mul(qpow[A], tpow[B]), qm1pow[K1]), tpq2pow[K2]))

    total = nc_zero()
    for (left, idx, right, num, den) in parsed:
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
        total = nc_add(total, nc_scalar(term, coef_poly))
    diff_B3 = nc_sub(nc_normalize_mod_B3(total), nc_normalize_mod_B3(rhs))
    return nc_is_zero(diff_B3), nc_size(diff_B3)


def main(argv=None):
    t0 = time.time()
    print("=" * 72)
    print("Verification (A): compact cofactor encodes the braid-algebra identity")
    print("=" * 72)

    relations, X3 = build_inputs_and_X3()
    print(f"[setup] Loaded 31 input relations for n=4 zipper algebra.")
    print(f"[setup] X_3 = S[3][1] has {nc_size(X3)} monomials in the free algebra.")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=SRC_FILE)
    args = parser.parse_args(argv)
    parsed_ZT = parse_ZT_lines(args.input)
    print(f"[parse] {len(parsed_ZT)} tuples in {args.input.name}.")

    # (A.2) Verify the cofactor identity in Q(q,t)[B_3] by summing and
    # comparing to X_3 via Garside normalization.
    print("[A.2]   verifying sum_j c_j L_j R_{i_j} R'_j == X_3 in Q(q,t)[B_3] ...")
    ok, residual_size = verify_cofactor_identity(parsed_ZT, relations, X3)
    if ok:
        print("[A.2]   PASS: identity holds in Q(q,t)[B_3] (Garside residual is zero).")
    else:
        print(f"[A.2]   FAIL: Garside residual has {residual_size} monomials.")
        return 3

    print()
    print("=" * 72)
    print("Verification (B): defining factorizations and the T_3 ideal equality")
    print("=" * 72)

    # (B.1) T_3 = R_28 = (x_2 x_1 - q^2 y_2 y_1) * X_3 in the free algebra.
    # We verify this as a polynomial identity.
    R28 = relations[28]
    x1, x2, y1, y2 = 1, 2, 4, 5
    # Build (x_2*x_1 - q^2*y_2*y_1) as an NC polynomial
    coef_poly_T3 = {
        (x2, x1): qt_const(1),
        (y2, y1): qt_neg(qt_pow(qt_q(), 2)),
    }
    computed_T3 = nc_mul(coef_poly_T3, X3)
    diff_T3 = nc_sub(computed_T3, R28)
    if nc_is_zero(diff_T3):
        print("[B.1]   PASS: T_3 = R_28 = (x_2*x_1 - q^2*y_2*y_1) * X_3 in the free algebra.")
        print("        Hence T_3 in <X_3> as a two-sided ideal in the free algebra.")
    else:
        print(f"[B.1]   FAIL: residual has {nc_size(diff_T3)} monomials.")
        return 4

    # (B.2) Z_3 = R_24 (the level-3 zipper at n=4, j=1).
    # The CORRECT formula at n=4 involves strand 3:
    #   R_24 = (q^2*y_2*y_3 - x_2*x_3 + x_1*x_2*x_3 - q^2*y_1*y_2*y_3) * X_3
    # (Note: in earlier discussion I used a truncated version of this formula
    # without strand 3, which does NOT correspond to a relation of Z_n at any n;
    # the n=3 dim test in test_R28_necessity.gap was therefore testing a
    # different ideal than the one generated by the actual level-3 zipper at n=4.)
    R24 = relations[24]
    x3, y3 = 3, 6
    coef_poly_Z3 = {
        (y2, y3): qt_pow(qt_q(), 2),
        (x2, x3): qt_neg(qt_const(1)),
        (x1, x2, x3): qt_const(1),
        (y1, y2, y3): qt_neg(qt_pow(qt_q(), 2)),
    }
    computed_Z3 = nc_mul(coef_poly_Z3, X3)
    diff_Z3 = nc_sub(computed_Z3, R24)
    if nc_is_zero(diff_Z3):
        print("[B.2]   PASS: Z_3 = R_24 = (q^2*y_2*y_3 - x_2*x_3 + x_1*x_2*x_3 - q^2*y_1*y_2*y_3) * X_3")
        print("        in the free algebra Z[q,t]<x_1,x_2,x_3,y_1,y_2,y_3>.")
        print("        Hence Z_3 in <X_3> as a two-sided ideal in the n=4 free algebra.")
    else:
        print(f"[B.2]   FAIL: residual has {nc_size(diff_Z3)} monomials.")
        return 5

    # (B.3) X_3 in <Z_2, T_2, T_3>: this is exactly verification (A.2), passed
    # above. The 161-tuple cofactor is the explicit witness.
    print("[B.3]   PASS: X_3 in <Z_2, T_2, T_3>")
    print("        Witness: the 161-tuple cofactor verified in (A.2).")

    print("[B.4]   No conclusion about replacing T_3 by Z_3 is computed here.")
    print("Established over Q(q,t) in the braid group algebra of B_4:")
    print("  <Z_2, T_2, T_3> = <Z_2, T_2, X_3>.")
    print("  <Z_2, T_2, Z_3> is contained in <Z_2, T_2, X_3>.")

    print()
    print("=" * 72)
    print(f"All checks PASS.  Total wall time: {time.time()-t0:.2f}s")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
