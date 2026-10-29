"""Infer the historical denominator family from deterministic synthetic samples.

Run after installing the package with its optional inference extra:
    python examples/infer_denominators.py
"""
from fractions import Fraction
import random

from zipper_algebra.cofactor.atoms import AtomSet
from zipper_algebra.cofactor.inference import infer


def main():
    random.seed(1)
    six = AtomSet.DEFAULT
    pts = [(q, t) for q in (3, 5, 7, 8, 9, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
           for t in (3, 5, 7, 11, 13, 17) if q != t]
    slots = {}
    for si, exps in enumerate([(1, 2, 1, 0, 0, 0), (0, 0, 0, 1, 1, 0), (2, 0, 0, 0, 1, 1), (0, 1, 1, 1, 0, 1)]):
        num = {(i, j): random.randint(-5, 5) for i in range(3) for j in range(2)}
        slots[si] = [(p, Fraction(sum(c * p[0] ** i * p[1] ** j for (i, j), c in num.items()),
                                  six.eval_denom(exps, *p))) for p in pts]
    got, rep = infer(slots, seed=AtomSet.UNITS)
    print(got, "unexplained:", rep["unexplained"])
    assert set(got.atoms) == set(six.atoms), got.names()
    assert rep["unexplained"] == 0
    print("OK: from the units alone, rediscovered the six from denominators")


if __name__ == "__main__":
    main()
