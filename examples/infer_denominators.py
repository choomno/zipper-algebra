"""Discover unfamiliar denominator factors from exact synthetic observations.

Run after installing the package with its optional inference extra:
    python examples/infer_denominators.py

The formula constructs test data. Only its sampled rational values are passed
into inference: no denominator family or seed is supplied. Both factors lie
outside the small initial candidate box, so their discovery uses reconstruction
and factorization of the sampled line denominators. This is a heuristic proposal;
production certificates still require their separate exact identity check.
"""
from fractions import Fraction
import importlib.util

from zipper_algebra.cofactor.atoms import Atom
from zipper_algebra.cofactor.inference import infer


def main():
    if importlib.util.find_spec("sympy") is None:
        raise SystemExit("Install the inference extra from the project root: "
                         "python -m pip install '.[inference]'")
    observations = {
        0: [((q, t), Fraction(q + 2 * t, (q + t + 5) * (q * q + t + 7)))
            for q in range(2, 15) for t in range(17, 30)]
    }
    family, report = infer(observations, max_bideg=(1, 1), coef_bound=1, verbose=False)

    # These expected factors check the synthetic experiment; they are not inputs
    # to inference. The initial box excludes their constants and the q^2 term.
    expected = {
        Atom("q+t+5", {(1, 0): 1, (0, 1): 1, (0, 0): 5}),
        Atom("q^2+t+7", {(2, 0): 1, (0, 1): 1, (0, 0): 7}),
    }
    assert set(family.atoms) == expected, family.names()
    assert report["seed"] == []
    assert report["unexplained"] == 0
    print("Discovered denominator factors:", ", ".join(family.names()))
    print("No denominator family or seed was supplied.")
    return family, report


if __name__ == "__main__":
    main()
