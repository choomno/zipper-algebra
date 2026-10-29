#!/usr/bin/env python3
"""Selected-factor and sampled-gcd diagnostics for a fixed cofactor expression.

The program records vanishing at selected substitutions and integer gcds at
five reproducibly sampled points. Positive probe results do not certify a
polynomial common factor or minimality among all certificates.
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path
import sys

import random
from fractions import Fraction
from math import gcd
from functools import reduce
from .legacy_format import parse_cofactor_line



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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="candidate expression to examine")
    args = parser.parse_args(argv)
    # Parse minimal-file tuples
    tuples = []
    for ln in open(args.input):
        s = ln.strip()
        if not s or s.startswith('#'): continue
        tuples.append(parse_cofactor_line(ln))
    print(f"[parse] {len(tuples)} tuples from {os.path.basename(args.input)}.")

    nums = [t[3] for t in tuples]

    # ----- (1) Atom-level checks -----
    # Record selected substitutions. Nonzero values refute the candidate divisor;
    # vanishing at every selected point is not a general divisibility certificate.
    print()
    print("(1) Selected-substitution factor probes:")
    print("   For each candidate atom, count numerators vanishing at its selected substitutions.")
    print("   All-zero probes are inconclusive; a nonzero probe refutes that candidate factor.")
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
        verdict = "all probes vanish (not a divisibility proof)" if count == len(nums) else "candidate refuted by a nonzero probe"
        print(f"   {name:<8} : vanishing probes in {count:>3}/{len(nums)} tuples   ->  {verdict}")

    # ----- (2) Polynomial-gcd: probabilistic random-point check -----
    # Random integer (q,t) avoiding the three vanishing hypersurfaces of D'.
    print()
    print("(2) Random-point integer-gcd probe:")
    print("    Evaluate every M_j' at random integer (q,t) and take gcd of all 161 values.")
    print("    These values record sampled arithmetic only; no symbolic gcd is computed.")
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
    print("The sampled integer gcds are diagnostics, not a polynomial gcd certificate.")
    print("The atom checks concern this fixed expression, not all possible certificates.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
