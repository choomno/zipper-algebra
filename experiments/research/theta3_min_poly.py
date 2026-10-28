#!/usr/bin/env python3
r"""Probe three-strand lasso-twist annihilators at specified parameter points.

The historical candidate is the reflected two-strand polynomial p*. Agreement
of Krylov calculations at finitely many points does not establish its generic
validity or minimality. The supplied solver output has no basis-completeness
certificate in this recipe. The independent exact route is retained as
 gap_scripts/theta3_exact.gap in the repository root.

Run as python -m experiments.research.theta3_min_poly after installing the package.
The presentation is emitted by the maintained Presentation constructor.
"""
from __future__ import annotations
import argparse, os, subprocess, sys

from zipper_algebra import core, Presentation
from zipper_algebra.linear import validate_field_characteristic
from experiments.variants import sigma_min_poly as smp  # min_poly_on_vector, poly utils

F4_DEFAULT = os.environ.get("F4NCGB", "f4ncgb")
THETA3_WORD = ("x1", "x2", "x2", "x1")          # Theta_3 = sigma_1 sigma_2^2 sigma_1

# p_2 = (sigma-t)^2 : --untwist square-zero ;  p_2 = sigma^2-t : --untwist bigelow --framing-power 2
CASES = {"linear": ("bigelow", 1), "square-zero": ("square-zero", 1), "framing": ("bigelow", 2)}


def ncmul(A, B, p):
    out = {}
    for ca, wa in A:
        for cb, wb in B:
            w = wa + wb
            out[w] = (out.get(w, 0) + ca * cb) % p
    return [(c % p, w) for w, c in out.items() if c % p]


def build_X3(q, p):
    """X_3 = (q^2 y1 y2 - x1 x2)(q y1 + (1-q) - x1) as list[(coef, word)] mod p."""
    q %= p
    f1 = [(q * q % p, ("y1", "y2")), ((-1) % p, ("x1", "x2"))]
    f2 = [(q, ("y1",)), ((1 - q) % p, ()), ((-1) % p, ("x1",))]
    return ncmul(f1, f2, p)


def pstar_expected(case, t, p):
    """p*(lambda) coefficients low->high at the specialization t (mod p)."""
    t %= p
    return {"linear": [t, 1],
            "square-zero": [t * t % p, (2 * t) % p, 1],
            "framing": [(-t) % p, 0, 1]}[case]


def compute(case, q, t, p, f4, scratch, maxiter=5000):
    """Return (dim, min_poly_coeffs_or_None). None min poly means X_3 = 0."""
    validate_field_characteristic(p)
    if not p:
        raise ValueError("this Krylov recipe requires a positive prime characteristic")
    untw, fp = CASES[case]
    ms = os.path.join(scratch, f"theta3-{case}-{q}-{t}.ms")
    bs = os.path.join(scratch, f"theta3-{case}-{q}-{t}.basis")
    pres = Presentation(3, untwist=untw, framing_power=fp, max_untwist_level=2)
    pres.emit_ms(ms, q, t, p)
    subprocess.run([f4, "-m", str(maxiter), "-i", ms, "-o", bs], check=True, capture_output=True)
    polys = core.parse_basis_Fp(bs)
    lms = [lm for lm, _ in polys]
    variables = pres.var_names()
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(delta, is_match, variables, max_total=200000)
    if standard is None:
        return "enumeration-limit", None
    N = len(standard)
    mon_idx = {m: i for i, m in enumerate(standard)}
    lm_to_tail = core.lm_to_tail_Fp(polys, p)
    v0 = core.reduce_to_basis_Fp(build_X3(q, p), mon_idx, lm_to_tail, delta, is_match, match_pat, p)
    if not v0:
        return N, None
    M = {}
    for j, m_j in enumerate(standard):
        red = core.reduce_to_basis_Fp([(1, THETA3_WORD + m_j)], mon_idx, lm_to_tail,
                                      delta, is_match, match_pat, p)
        if red:
            M[j] = red
    return N, smp.min_poly_on_vector(M, v0, N, p)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", choices=list(CASES) + ["all"], default="all")
    ap.add_argument("--grid", default="5,7 11,13 3,10 13,17 2,9 7,4 10,3",
                    help="space-separated q,t points")
    ap.add_argument("--char", type=int, default=2**31 - 1)
    ap.add_argument("--f4ncgb", default=F4_DEFAULT)
    ap.add_argument("--scratch", default="/tmp")
    args = ap.parse_args()
    p = args.char
    grid = [tuple(int(x) for x in pt.split(",")) for pt in args.grid.split()]
    cases = list(CASES) if args.case == "all" else [args.case]
    for case in cases:
        print(f"\n=== {case}: min poly of Theta_3 on X_3 (p={p}) ===")
        allok = True
        for (q, t) in grid:
            N, mp = compute(case, q, t, p, args.f4ncgb, args.scratch)
            if N == "enumeration-limit":
                print("  Standard-word enumeration budget exhausted; no dimension or vanishing result.")
                allok = False
                continue
            if mp is None:
                print(f"  (q,t)=({q},{t}): dim={N}  X_3 = 0"); allok = False; continue
            exp = pstar_expected(case, t, p)
            ok = (mp == exp); allok = allok and ok
            tag = "== p*" if ok else "MISMATCH exp " + smp.fmt_poly(exp, p, "L")
            print(f"  (q,t)=({q:>2},{t:>2}): dim={N}  minpoly = {smp.fmt_poly(mp, p, 'L'):22s} {tag}")
        print(f"  --> min poly = p*(t) at all grid points: {allok}")


if __name__ == "__main__":
    main()
