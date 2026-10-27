"""Exercise the optional exact polynomial-factorization dependency."""
from fractions import Fraction

from zipper_algebra.cofactor.atoms import residual_factors


if __name__ == "__main__":
    # q^2-1 factors into q-1 and q+1 over the rationals.
    result = residual_factors([("t", 3, [Fraction(-1), Fraction(0), Fraction(1)])],
                              verbose=False)
    assert set(result) == {("t", (-1, 1)), ("t", (1, 1))}, result
    print("Optional exact factorization passed")
