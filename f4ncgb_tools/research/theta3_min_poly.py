#!/usr/bin/env python3
"""Derived three-strand untwisting: minimal polynomial of the lasso twist
Theta_3 on the slider X_3 in  Z_3 / (zipper relations + two-strand untwisting p_2),
with NO three-strand relation imposed.

This operationalizes the finitization program: impose only the least needed to
finitize (here the two-strand untwisting alone suffices — the algebra is finite
even though the slider *module* is infinite), then read the higher untwisting off
as the min poly of the lasso twist.

FINDING (this script, mod-p at a (q,t) grid; confirmed exactly over Q(q,t) for the
linear case in ../../gap_scripts/theta3_exact.gap): for
    p_2 = sigma - t      -> min poly = lambda + t         = p*
    p_2 = (sigma - t)^2   -> min poly = (lambda + t)^2     = p*
    p_2 = sigma^2 - t     -> min poly = lambda^2 - t       = p*
i.e. the algebra's derived relation is EXACTLY p* (the reflection of p), dropping
the rigid (lambda - q^2) factor the module retains (thm:modelp gives only the
divisibility ceiling (lambda - q^2) p*). This confirms conj:recursion's sharpening
and explains Conjecture 1.

TRUST: mod-p at a grid is a SCOUT (speculative). The load-bearing statement is the
exact Q(q,t) confirmation in gap_scripts/theta3_exact.gap. Use --char 0 grid +
reconstruction (TODO: corelib.reconstruct) for a self-contained exact result.

Pipeline: zipper_in_ms (emit two-strand-only presentation) -> f4ncgb (GB) ->
core (parse/AC/reduce) -> min poly of left-mult-by-Theta_3 on the reduced X_3.
The min-poly / F_p-linear-algebra helpers are reused from variants.sigma_min_poly
(TODO: hoist into corelib.fieldlinalg per the refactor plan).
"""
from __future__ import annotations
import argparse, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
if REPO not in sys.path:
    sys.path.insert(0, REPO)
from f4ncgb_tools import core
from f4ncgb_tools.variants import sigma_min_poly as smp  # min_poly_on_vector, poly utils

F4_DEFAULT = "/home/choom/src/f4ncgb/build-release/f4ncgb"
EMIT = os.path.join(REPO, "f4ncgb_tools", "zipper_in_ms.py")
THETA3_WORD = ("x1", "x2", "x2", "x1")          # Theta_3 = sigma_1 sigma_2^2 sigma_1

# p_2 = (sigma-t)^2 : --untwist jordan ;  p_2 = sigma^2-t : --untwist bigelow --framing-power 2
CASES = {"linear": ("bigelow", 1), "jordan": ("jordan", 1), "framing": ("bigelow", 2)}


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
            "jordan": [t * t % p, (2 * t) % p, 1],
            "framing": [(-t) % p, 0, 1]}[case]


def compute(case, q, t, p, f4, scratch, maxiter=5000):
    """Return (dim, min_poly_coeffs_or_None). None min poly means X_3 = 0."""
    untw, fp = CASES[case]
    ms = os.path.join(scratch, f"theta3-{case}-{q}-{t}.ms")
    bs = os.path.join(scratch, f"theta3-{case}-{q}-{t}.basis")
    subprocess.run(["python3", EMIT, "--n", "3", "--untwist", untw, "--framing-power", str(fp),
                    "--q", str(q), "--t", str(t), "--char", str(p),
                    "--max-untwist-level", "2", "--out", ms], check=True, capture_output=True)
    subprocess.run([f4, "-m", str(maxiter), "-i", ms, "-o", bs], check=True, capture_output=True)
    polys = core.parse_basis_Fp(bs)
    lms = [lm for lm, _ in polys]
    variables = core.collect_variables(polys)
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(delta, is_match, variables, max_total=200000)
    if standard is None:
        return "infinite", None
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
            if mp is None:
                print(f"  (q,t)=({q},{t}): dim={N}  X_3 = 0"); allok = False; continue
            exp = pstar_expected(case, t, p)
            ok = (mp == exp); allok = allok and ok
            tag = "== p*" if ok else "MISMATCH exp " + smp.fmt_poly(exp, p, "L")
            print(f"  (q,t)=({q:>2},{t:>2}): dim={N}  minpoly = {smp.fmt_poly(mp, p, 'L'):22s} {tag}")
        print(f"  --> min poly = p*(t) at all grid points: {allok}")


if __name__ == "__main__":
    main()
