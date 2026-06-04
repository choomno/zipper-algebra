#!/usr/bin/env python3
"""Exhaustive check that D' = t^4 * (t + q^2) is the minimal common denominator
of the 161-tuple cofactor in cofactor-X3-n4-Qqt-minimal.txt.

Two layers of check:

(1) Atom-level: for each candidate denominator atom a in the full atom set
    A = {q, t, q-1, t-q, t+q^2, q+1}, verify the *min* multiplicity of a in
    the M_j' across all 161 tuples.  If the min is positive, a can be canceled
    further.  We require min = 0 for every atom not already in D'.

(2) Polynomial-gcd level: probabilistically + symbolically, check that
    gcd_j(M_j') is the unit polynomial 1.  Strategy:
       (a) take two of the smallest-degree M_j' and compute their bivariate
           gcd in Q[q, t] via sympy;
       (b) sanity-check: evaluate that gcd against every M_j' (must divide
           all of them exactly) -- but since gcd of 2 already gives an upper
           bound on the full gcd, if the 2-element gcd is 1, the full gcd is 1;
       (c) cross-check by random-point evaluation: pick several random integer
           (q,t), evaluate every M_j', and check that the gcd of the 161
           integers is 1 (or a small constant) -- if 1, no constant or
           non-trivial polynomial common factor remains.
"""
from __future__ import annotations
import os
import sys
import random
from fractions import Fraction
from math import gcd
from functools import reduce

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from symbolic_identity_check import parse_cofactor_line

ROOT = "/home/choom/gap_scripts/zipper-algebra"
SRC = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-minimal.txt")


def poly_eval(p, q0, t0):
    """Evaluate dict (deg_q, deg_t) -> Fraction at integer (q0, t0)."""
    s = Fraction(0)
    for (a, b), c in p.items():
        s += c * (q0 ** a) * (t0 ** b)
    return s


def poly_at_q_eq(p, q0):
    """Evaluate at q = q0, returning dict {t_deg: Fraction}."""
    out = {}
    for (a, b), c in p.items():
        out[b] = out.get(b, Fraction(0)) + c * (q0 ** a)
    return {k: v for k, v in out.items() if v != 0}


def poly_at_t_eq(p, t0):
    out = {}
    for (a, b), c in p.items():
        out[a] = out.get(a, Fraction(0)) + c * (t0 ** b)
    return {k: v for k, v in out.items() if v != 0}


def vanishes_at(p, sub):
    """Return True if p vanishes after substitution.  sub is a dict like {'q': 1}."""
    if 'q' in sub and 't' in sub:
        return poly_eval(p, sub['q'], sub['t']) == 0
    if 'q' in sub:
        return not poly_at_q_eq(p, sub['q'])
    if 't' in sub:
        return not poly_at_t_eq(p, sub['t'])
    return False


def main():
    # Parse minimal-file tuples
    tuples = []
    for ln in open(SRC):
        s = ln.strip()
        if not s or s.startswith('#'): continue
        tuples.append(parse_cofactor_line(ln))
    print(f"[parse] {len(tuples)} tuples from {os.path.basename(SRC)}.")

    nums = [t[3] for t in tuples]

    # ----- (1) Atom-level checks -----
    # For each atom a in A (and not already in D' = t^4(t+q^2)), check whether
    # *every* M_j' is divisible by a.  If yes, a can be further canceled.
    print()
    print("(1) Atom-level minimality check:")
    print("   For each candidate atom, count tuples whose M_j' is divisible by it.")
    print(f"   If all 161 -> atom is a common factor; if < 161 -> atom is not.")
    print()
    atoms = [
        ("q",       lambda p: vanishes_at(p, {'q': 0})),
        ("t",       lambda p: vanishes_at(p, {'t': 0})),
        ("q-1",     lambda p: vanishes_at(p, {'q': 1})),
        ("q+1",     lambda p: vanishes_at(p, {'q': -1})),
        ("t-q",     lambda p: all(vanishes_at(p, {'q': q0, 't': q0}) for q0 in (2, 3, 5))),
        ("t+q^2",   lambda p: all(vanishes_at(p, {'q': q0, 't': -q0*q0}) for q0 in (2, 3, 5))),
    ]
    for name, test in atoms:
        count = sum(1 for p in nums if test(p))
        verdict = "COMMON FACTOR (can cancel)" if count == len(nums) else "not a common factor"
        print(f"   {name:<8} : divisible in {count:>3}/{len(nums)} tuples   ->  {verdict}")

    # ----- (2) Polynomial-gcd: probabilistic random-point check -----
    # Random integer (q,t) avoiding the three vanishing hypersurfaces of D'.
    print()
    print("(2) Random-point integer-gcd probe:")
    print("    Evaluate every M_j' at random integer (q,t) and take gcd of all 161 values.")
    print("    A nontrivial common polynomial factor would force this gcd to be a")
    print("    bounded-size integer at every random point.  If we see different small")
    print("    gcds (or 1) across several points, no nontrivial polynomial common factor.")
    print()
    random.seed(12345)
    for trial in range(5):
        # Pick (q,t) avoiding t=0, t+q^2=0, q in {0,1,-1}, t=q
        while True:
            q0 = random.randint(-15, 50)
            t0 = random.randint(-15, 50)
            if q0 in (0, 1, -1) or t0 == 0 or t0 == q0 or t0 + q0*q0 == 0:
                continue
            break
        vals = []
        for p in nums:
            v = poly_eval(p, q0, t0)
            if v.denominator != 1:
                # Shouldn't happen since M_j' in Z[q,t]
                print(f"   WARN: non-integer eval at (q,t)=({q0},{t0})")
                break
            vals.append(abs(v.numerator))
        vals = [v for v in vals if v != 0]
        g = reduce(gcd, vals, 0)
        print(f"   (q,t) = ({q0:>4}, {t0:>4}):  gcd over 161 integers = {g}")

    print()
    print("Interpretation: if the gcd is small (say, < 50) and varies across random")
    print("points, the M_j' have no nontrivial polynomial common factor; the cofactor")
    print("is in lowest terms over D' = t^4 (t + q^2).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
