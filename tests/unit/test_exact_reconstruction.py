"""Exact synthetic regressions for reconstruction, without solvers or research files."""
from __future__ import annotations

import contextlib
from fractions import Fraction
import io
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import pytest
from unittest.mock import patch

from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.univariate import (_poly_eval, _reduce_univariate_fraction,
                                               reconstruct_univariate)
from zipper_algebra.cofactor.lift import (NeedMoreSamples, fit_poly_qt, fit_rational,
                   fit_rational_known, is_reduced, lift, line_denominator, poly_eval_qt)
from zipper_algebra.coefficients import NumericRing
from zipper_algebra.polynomial import Poly
from zipper_algebra.presentation import Presentation
from zipper_algebra.cofactor.proof import Basis, Proof


class TestExactReconstruction:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    def test_underdetermined_polynomial_box_is_not_filled_with_arbitrary_zeros(self):
        points = [((q, q), Fraction(q)) for q in range(1, 9)]
        assert fit_poly_qt(points, 1, 1) is None
        identified = points + [((1, 2), Fraction(1)), ((2, 1), Fraction(2))]
        assert fit_poly_qt(identified, 1, 1) == {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 0}

    def test_constant_rational_numerator_is_allowed(self):
        points = [((q, t), Fraction(1, 2)) for q, t in ((2, 3), (3, 4), (4, 5))]
        assert fit_rational(points, AtomSet([]), 0, 0) == ((), {(0, 0): Fraction(1, 2)})

    def test_rational_coefficients_with_nontrivial_denominator(self):
        atoms = AtomSet([Atom("q", {(1, 0): 1})], [1])
        points = [((q, 2), Fraction(1, 2 * q)) for q in range(2, 9)]
        for function in (lambda: fit_rational(points, atoms, 0, 0),
                         lambda: fit_rational_known(points, atoms, (1,), 0, 0)):
            assert function() == ((1,), {(0, 0): Fraction(1, 2)})

    def test_composite_atom_shortcut_only_checks_whole_atom_cancellation(self):
        atoms = AtomSet([Atom("q*t", {(1, 1): 1})], [1])
        points = [((q, t), Fraction(1, t)) for q, t in ((2, 3), (3, 5), (5, 7), (7, 11))]
        denominator, numerator = fit_rational_known(points, atoms, (1,), 1, 0)
        assert numerator == {(0, 0): 0, (1, 0): 1}
        assert is_reduced(atoms, denominator, numerator)
        # The whole atom q*t is not removable, although its proper factor q is.
        assert poly_eval_qt(numerator, 13, 17) / atoms.eval_denom(denominator, 13, 17) == Fraction(1, 17)

    def test_specialization_collision_is_a_heuristic_not_an_exponent_bound(self):
        atoms = AtomSet([Atom("q", {(1, 0): 1}),
                         Atom("q+t", {(1, 0): 1, (0, 1): 1})], [2, 2])
        points = [((q, t), Fraction(1, q * (q + t)))
                  for q in range(2, 12) for t in (0, 1)]
        proposed = line_denominator(points, atoms)
        assert proposed == (2, 1)
        # At t=0 the two distinct factors coincide. The q exponent is overestimated.
        assert fit_rational_known(points, atoms, proposed, 1, 0) is None
        recovered = fit_rational(points, atoms, 1, 0)
        assert recovered == ((1, 1), {(0, 0): 1})

    def test_every_sample_can_supply_rank_beyond_a_degenerate_prefix(self):
        points = [((1, t), Fraction(1)) for t in range(1, 7)]
        points += [((2, 1), Fraction(2)), ((3, 1), Fraction(3))]
        candidate = fit_rational(points, AtomSet([]), 1, 0)
        assert candidate == ((), {(0, 0): 0, (1, 0): 1})

    def test_duplicate_locations_count_once_and_conflicts_are_rejected(self):
        assert fit_poly_qt([((1, 2), 3)] * 5, 1, 0) is None
        with pytest.raises(ValueError, match='conflicting'):
            fit_poly_qt([((1, 2), 3), ((1, 2), 4)], 0, 0)
        with pytest.raises(ValueError, match='conflicting'):
            reconstruct_univariate([(1, 2), (1, 3)])
        points = [(x, Fraction(1, x + 1)) for x in range(1, 9)]
        assert reconstruct_univariate(points) == reconstruct_univariate(points * 2)

    def test_inexact_or_malformed_samples_and_degree_bounds_are_rejected(self):
        for points in ([((1.0, 2), 3)], [((1, 2), 0.5)], [((True, 2), 3)],
                       [((1,), 3)], [((1, 2), "3")]):
            with self.subtests.test(points=points):
                with pytest.raises(ValueError):
                    fit_rational(points, AtomSet([]))
        for bounds in ((-1, 0), (True, 0), (1.5, 0)):
            with self.subtests.test(bounds=bounds):
                with pytest.raises(ValueError):
                    fit_poly_qt([], *bounds)
        with pytest.raises(ValueError):
            reconstruct_univariate([(1.0, 2)])

    def test_univariate_rational_reconstruction_is_reduced_and_defined(self):
        points = [(x, Fraction(2 * x + 1, (x + 1) * (x + 2))) for x in range(1, 13)]
        numerator, denominator = reconstruct_univariate(points)
        assert numerator == [1, 2]
        assert denominator == [2, 3, 1]
        for x, expected in points:
            assert _poly_eval(denominator, x) != 0
            assert _poly_eval(numerator, x) / _poly_eval(denominator, x) == expected
        assert _reduce_univariate_fraction([-1, 1], [-1, 1]) == ([1], [1])
        assert _reduce_univariate_fraction([0, 0], [-1, 1]) == ([0], [1])

    def test_a_removable_pole_cannot_hide_a_wrong_fitting_value(self):
        # The old cross-multiplied fit N=D=x-1 accepted the point (1,0) via 0=0,
        # while every other point and both holdouts had value 1.
        assert reconstruct_univariate([(1, 0)] + [(x, 1) for x in range(2, 9)]) is None

    @pytest.mark.parametrize("holdout", [0, -1, True, 1.5])
    @pytest.mark.parametrize("operation", [
        pytest.param(lambda holdout: reconstruct_univariate([], holdout=holdout), id="univariate"),
        pytest.param(lambda holdout: lift(SimpleNamespace(anchor=(1, 2, 0)), None,
                                          holdout=holdout), id="lift"),
    ])
    def test_invalid_holdouts_are_rejected(self, holdout, operation):
        with pytest.raises(ValueError):
            operation(holdout)

    def test_pole_in_supplied_holdout_rejects_the_lift(self):
        points = [(q, 0) for q in range(1, 7)] + [(0, 0)]
        samples = synthetic_samples(points, lambda q, t: Fraction(1, q) if q else Fraction(99))
        with pytest.raises(NeedMoreSamples) as context:
            lift(synthetic_skeleton(), samples, AtomSet([Atom("q", {(1, 0): 1})], [1]),
                 holdout=1, max_pq=0, max_pt=0, verbose=False)
        assert context.value.slots == [('x1', 0, 'holdout')]

    def test_failed_holdout_is_inconclusive_and_more_data_can_succeed(self):
        points = [(q, 1) for q in range(1, 7)] + [(1, 2), (2, 2)]
        function = lambda q, t: Fraction(q + (t - 1) * q * q)
        samples = synthetic_samples(points, function)
        with pytest.raises(NeedMoreSamples) as context:
            lift(synthetic_skeleton(), samples, AtomSet([]), holdout=2,
                 max_pq=2, max_pt=1, verbose=False)
        assert context.value.slots == [('x1', 0, 'holdout')]
        expanded = points + [(q, t) for q in range(3, 9) for t in (2, 3)]
        result = lift(synthetic_skeleton(), synthetic_samples(expanded, function),
                      AtomSet([]), holdout=2, max_pq=2, max_pt=1, verbose=False)
        numerator = result.entries["x1"][0][1]
        assert poly_eval_qt(numerator, 10, 4) == function(10, 4)

    def test_inexact_holdout_value_is_rejected(self):
        samples = synthetic_samples([(q, 0) for q in range(1, 9)],
                                    lambda q, t: Fraction(1) if q < 8 else 1.0)
        with pytest.raises(ValueError, match='exact rational'):
            lift(synthetic_skeleton(), samples, AtomSet([]), holdout=1, verbose=False)

    def test_repeated_samples_do_not_supply_independent_holdout_data(self):
        samples = synthetic_samples([(q, 0) for q in range(1, 7)] * 2,
                                    lambda q, t: Fraction(1))
        with pytest.raises(NeedMoreSamples):
            lift(synthetic_skeleton(), samples, AtomSet([]), verbose=False)

    def test_atom_inference_uses_only_fitting_coefficients(self):
        points = [(q, 0) for q in range(1, 9)]
        samples = synthetic_samples(points, lambda q, t: Fraction(1))
        with patch.object(sys.modules[lift.__module__], "infer_atoms",
                          return_value=(AtomSet([]), {"unexplained": 0})) as inference:
            lift(synthetic_skeleton(), samples, holdout=2, verbose=False)
        values = inference.call_args.args[0]
        for slot_points in values.values():
            assert [tuple(point) for point, value in slot_points] == points[:-2]

    def test_finite_characteristic_lift_is_not_rational_reconstruction(self):
        with pytest.raises(ValueError, match='characteristic-zero'):
            lift(SimpleNamespace(anchor=(1, 2, 5)), None)
    def test_known_sample_characteristic_must_match_rational_domain(self):
        samples = synthetic_samples([(q, 0) for q in range(1, 9)], lambda q, t: Fraction(1))
        samples.char = 5
        with pytest.raises(ValueError, match='characteristic-zero'):
            lift(synthetic_skeleton(), samples, AtomSet([]), verbose=False)


def synthetic_skeleton():
    return SimpleNamespace(anchor=(1, 0, 0), steps=[{"basis_key": "x1"}])


def synthetic_samples(points, function):
    basis = SimpleNamespace(key=["x1"], index_by_key={"x1": 0})
    return SimpleNamespace(points=points, basis=lambda q, t: basis,
        proof=lambda q, t: SimpleNamespace(terms=[[(function(q, t), (), "input", 0, ())]]))


