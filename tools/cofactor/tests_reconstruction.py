"""Exact synthetic regressions for reconstruction, without solvers or research files."""
from __future__ import annotations

import contextlib
from fractions import Fraction
import io
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch

from .atoms import (Atom, AtomSet, _poly_eval, _reduce_univariate_fraction,
                    reconstruct_univariate)
from .lift import (NeedMoreSamples, WrongClass, fit_poly_qt, fit_rational,
                   fit_rational_known, is_reduced, lift, line_denominator, poly_eval_qt)
from . import pipeline
from ..zipper_in_ms import Poly, NumericRing
from ..presentation import Presentation
from .proof import Basis, Proof


class TestExactReconstruction(unittest.TestCase):
    def test_underdetermined_polynomial_box_is_not_filled_with_arbitrary_zeros(self):
        points = [((q, q), Fraction(q)) for q in range(1, 9)]
        self.assertIsNone(fit_poly_qt(points, 1, 1))
        identified = points + [((1, 2), Fraction(1)), ((2, 1), Fraction(2))]
        self.assertEqual(fit_poly_qt(identified, 1, 1),
                         {(0, 0): 0, (0, 1): 0, (1, 0): 1, (1, 1): 0})

    def test_constant_rational_numerator_is_allowed(self):
        points = [((q, t), Fraction(1, 2)) for q, t in ((2, 3), (3, 4), (4, 5))]
        self.assertEqual(fit_rational(points, AtomSet([]), 0, 0),
                         ((), {(0, 0): Fraction(1, 2)}))

    def test_rational_coefficients_with_nontrivial_denominator(self):
        atoms = AtomSet([Atom("q", {(1, 0): 1})], [1])
        points = [((q, 2), Fraction(1, 2 * q)) for q in range(2, 9)]
        for function in (lambda: fit_rational(points, atoms, 0, 0),
                         lambda: fit_rational_known(points, atoms, (1,), 0, 0)):
            self.assertEqual(function(), ((1,), {(0, 0): Fraction(1, 2)}))

    def test_composite_atom_shortcut_only_checks_whole_atom_cancellation(self):
        atoms = AtomSet([Atom("q*t", {(1, 1): 1})], [1])
        points = [((q, t), Fraction(1, t)) for q, t in ((2, 3), (3, 5), (5, 7), (7, 11))]
        denominator, numerator = fit_rational_known(points, atoms, (1,), 1, 0)
        self.assertEqual(numerator, {(0, 0): 0, (1, 0): 1})
        self.assertTrue(is_reduced(atoms, denominator, numerator))
        # The whole atom q*t is not removable, although its proper factor q is.
        self.assertEqual(poly_eval_qt(numerator, 13, 17) /
                         atoms.eval_denom(denominator, 13, 17), Fraction(1, 17))

    def test_specialization_collision_is_a_heuristic_not_an_exponent_bound(self):
        atoms = AtomSet([Atom("q", {(1, 0): 1}),
                         Atom("q+t", {(1, 0): 1, (0, 1): 1})], [2, 2])
        points = [((q, t), Fraction(1, q * (q + t)))
                  for q in range(2, 12) for t in (0, 1)]
        proposed = line_denominator(points, atoms)
        self.assertEqual(proposed, (2, 1))
        # At t=0 the two distinct factors coincide. The q exponent is overestimated.
        self.assertIsNone(fit_rational_known(points, atoms, proposed, 1, 0))
        recovered = fit_rational(points, atoms, 1, 0)
        self.assertEqual(recovered, ((1, 1), {(0, 0): 1}))

    def test_every_sample_can_supply_rank_beyond_a_degenerate_prefix(self):
        points = [((1, t), Fraction(1)) for t in range(1, 7)]
        points += [((2, 1), Fraction(2)), ((3, 1), Fraction(3))]
        candidate = fit_rational(points, AtomSet([]), 1, 0)
        self.assertEqual(candidate, ((), {(0, 0): 0, (1, 0): 1}))

    def test_duplicate_locations_count_once_and_conflicts_are_rejected(self):
        self.assertIsNone(fit_poly_qt([((1, 2), 3)] * 5, 1, 0))
        with self.assertRaisesRegex(ValueError, "conflicting"):
            fit_poly_qt([((1, 2), 3), ((1, 2), 4)], 0, 0)
        with self.assertRaisesRegex(ValueError, "conflicting"):
            reconstruct_univariate([(1, 2), (1, 3)])
        points = [(x, Fraction(1, x + 1)) for x in range(1, 9)]
        self.assertEqual(reconstruct_univariate(points), reconstruct_univariate(points * 2))

    def test_inexact_or_malformed_samples_and_degree_bounds_are_rejected(self):
        for points in ([((1.0, 2), 3)], [((1, 2), 0.5)], [((True, 2), 3)],
                       [((1,), 3)], [((1, 2), "3")]):
            with self.subTest(points=points), self.assertRaises(ValueError):
                fit_rational(points, AtomSet([]))
        for bounds in ((-1, 0), (True, 0), (1.5, 0)):
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                fit_poly_qt([], *bounds)
        with self.assertRaises(ValueError):
            reconstruct_univariate([(1.0, 2)])

    def test_univariate_rational_reconstruction_is_reduced_and_defined(self):
        points = [(x, Fraction(2 * x + 1, (x + 1) * (x + 2))) for x in range(1, 13)]
        numerator, denominator = reconstruct_univariate(points)
        self.assertEqual(numerator, [1, 2])
        self.assertEqual(denominator, [2, 3, 1])
        for x, expected in points:
            self.assertNotEqual(_poly_eval(denominator, x), 0)
            self.assertEqual(_poly_eval(numerator, x) / _poly_eval(denominator, x), expected)
        self.assertEqual(_reduce_univariate_fraction([-1, 1], [-1, 1]), ([1], [1]))
        self.assertEqual(_reduce_univariate_fraction([0, 0], [-1, 1]), ([0], [1]))

    def test_a_removable_pole_cannot_hide_a_wrong_fitting_value(self):
        # The old cross-multiplied fit N=D=x-1 accepted the point (1,0) via 0=0,
        # while every other point and both holdouts had value 1.
        self.assertIsNone(reconstruct_univariate([(1, 0)] + [(x, 1) for x in range(2, 9)]))

    def test_invalid_holdouts_are_rejected(self):
        for holdout in (0, -1, True, 1.5):
            with self.subTest(holdout=holdout), self.assertRaises(ValueError):
                reconstruct_univariate([], holdout=holdout)
            with self.subTest(lift_holdout=holdout), self.assertRaises(ValueError):
                lift(SimpleNamespace(anchor=(1, 2, 0)), None, holdout=holdout)

    def test_pole_in_supplied_holdout_rejects_the_lift(self):
        points = [(q, 0) for q in range(1, 7)] + [(0, 0)]
        samples = synthetic_samples(points, lambda q, t: Fraction(1, q) if q else Fraction(99))
        with self.assertRaises(NeedMoreSamples) as context:
            lift(synthetic_skeleton(), samples, AtomSet([Atom("q", {(1, 0): 1})], [1]),
                 holdout=1, max_pq=0, max_pt=0, verbose=False)
        self.assertEqual(context.exception.slots, [("x1", 0, "holdout")])

    def test_failed_holdout_is_inconclusive_and_more_data_can_succeed(self):
        points = [(q, 1) for q in range(1, 7)] + [(1, 2), (2, 2)]
        function = lambda q, t: Fraction(q + (t - 1) * q * q)
        samples = synthetic_samples(points, function)
        with self.assertRaises(NeedMoreSamples) as context:
            lift(synthetic_skeleton(), samples, AtomSet([]), holdout=2,
                 max_pq=2, max_pt=1, verbose=False)
        self.assertIsInstance(context.exception, WrongClass)
        expanded = points + [(q, t) for q in range(3, 9) for t in (2, 3)]
        result = lift(synthetic_skeleton(), synthetic_samples(expanded, function),
                      AtomSet([]), holdout=2, max_pq=2, max_pt=1, verbose=False)
        numerator = result.entries["x1"][0][1]
        self.assertEqual(poly_eval_qt(numerator, 10, 4), function(10, 4))

    def test_inexact_holdout_value_is_rejected(self):
        samples = synthetic_samples([(q, 0) for q in range(1, 9)],
                                    lambda q, t: Fraction(1) if q < 8 else 1.0)
        with self.assertRaisesRegex(ValueError, "exact rational"):
            lift(synthetic_skeleton(), samples, AtomSet([]), holdout=1, verbose=False)

    def test_repeated_samples_do_not_supply_independent_holdout_data(self):
        samples = synthetic_samples([(q, 0) for q in range(1, 7)] * 2,
                                    lambda q, t: Fraction(1))
        with self.assertRaises(NeedMoreSamples):
            lift(synthetic_skeleton(), samples, AtomSet([]), verbose=False)

    def test_atom_inference_uses_only_fitting_coefficients(self):
        points = [(q, 0) for q in range(1, 9)]
        samples = synthetic_samples(points, lambda q, t: Fraction(1))
        with patch.object(sys.modules[lift.__module__], "infer_atoms",
                          return_value=(AtomSet([]), {"unexplained": 0})) as inference:
            lift(synthetic_skeleton(), samples, holdout=2, verbose=False)
        values = inference.call_args.args[0]
        for slot_points in values.values():
            self.assertEqual([tuple(point) for point, value in slot_points], points[:-2])

    def test_finite_characteristic_lift_is_not_rational_reconstruction(self):
        with self.assertRaisesRegex(ValueError, "characteristic-zero"):
            lift(SimpleNamespace(anchor=(1, 2, 5)), None)
    def test_known_sample_characteristic_must_match_rational_domain(self):
        samples = synthetic_samples([(q, 0) for q in range(1, 9)], lambda q, t: Fraction(1))
        samples.char = 5
        with self.assertRaisesRegex(ValueError, "characteristic-zero"):
            lift(synthetic_skeleton(), samples, AtomSet([]), verbose=False)


def synthetic_skeleton():
    return SimpleNamespace(anchor=(1, 0, 0), steps=[{"basis_key": "x1"}])


def synthetic_samples(points, function):
    basis = SimpleNamespace(key=["x1"], index_by_key={"x1": 0})
    return SimpleNamespace(points=points, basis=lambda q, t: basis,
        proof=lambda q, t: SimpleNamespace(terms=[[(function(q, t), (), "input", 0, ())]]))


class TestPipelineReconstructionOutcome(unittest.TestCase):
    def fake_run(self, lift_effect, verification=None, multiplier_failure=False, unstable_shapes=False,
                 unresolved_target=None):
        paths = SimpleNamespace(target="X_3", anchor_ms=Path("input.ms"),
                                skeleton=Path("skeleton.json"), lift=Path("lift.json"),
                                root=Path("out"), ensure_dirs=lambda: None,
                                cofactor=lambda phase: Path(f"cofactor-{phase}.txt"))
        disk = [0, 1, 2]
        paths.points_on_disk = lambda: disk
        growth = []
        def grow(count, jobs=None):
            growth.append(count)
            disk.extend(range(count))
        empty_basis = Basis([], 0, Presentation(3).var_names())
        samples = SimpleNamespace(run_anchor=lambda **kwargs: SimpleNamespace(
            hit=True, key="test", seconds=0), anchor_basis=lambda: empty_basis,
            anchor_proof=lambda: Proof([]),
            modal_shapes=lambda skeleton: ({}, ["x1"] if unstable_shapes else []),
            stable=lambda skeleton: [],
            grow=grow, _atoms=AtomSet([]))
        remainder_terms = len(unresolved_target.terms) if unresolved_target is not None else 0
        skel = SimpleNamespace(target="X_3", steps=[], basis_keys=[], remainder_terms=remainder_terms,
                               in_ideal=not remainder_terms, save=lambda path: None,
                               step_keys=lambda: [])
        pres = SimpleNamespace(slug=lambda: "test", content_hash=lambda: "hash",
                               n_relations=lambda: 0, emit_ms=lambda *args: 0)
        sampler_module = ModuleType(f"{pipeline.__package__}.samples")
        sampler_module.SampleSet = lambda *args, **kwargs: samples
        candidate = SimpleNamespace(terms=[], save=lambda path: None)
        valid_lift = SimpleNamespace(entries={}, save=lambda path: None,
                                    samples_used=4, atom_report=None, atoms=AtomSet([]))
        effects = [lift_effect] * 2 if isinstance(lift_effect, Exception) else [valid_lift] * 2
        actual_reduce = pipeline.reduce_target
        with patch.dict(sys.modules, {sampler_module.__name__: sampler_module}), \
             patch.object(pipeline, "Paths", return_value=paths), \
             patch.object(pipeline.compose_mod, "resolve_target",
                          return_value=unresolved_target if unresolved_target is not None else Poly(0)), \
             patch.object(pipeline, "reduce_target", return_value=([], Poly(0)),
                          side_effect=(lambda target, basis: actual_reduce(target, basis))
                          if unresolved_target is not None else None), \
             patch.object(pipeline.Skeleton, "from_reduction", return_value=skel), \
             patch.object(pipeline, "do_lift", side_effect=effects), \
             patch.object(pipeline, "_multiplier_samples", return_value={}), \
             patch.object(pipeline.compose_mod, "compose", return_value=candidate,
                          side_effect=[NeedMoreSamples("multiplier fit"), candidate]
                          if multiplier_failure else None), \
             patch.object(pipeline.verify_mod, "check", side_effect=verification), \
             contextlib.redirect_stdout(io.StringIO()):
            result = pipeline.run(pres, "X_3", (1, 2, 0), start=3, step=2, maxpts=4)
        return result, growth

    def test_failed_reconstruction_grows_to_cap_and_returns_inconclusive(self):
        result, growth = self.fake_run(WrongClass("failed held-out fit"))
        self.assertEqual(growth, [1])
        self.assertEqual(result["status"], "inconclusive")
        self.assertFalse(result["verified"])
        self.assertIsNone(result["in_ideal"])

    def test_failed_exact_candidate_can_retry_with_more_data(self):
        result, growth = self.fake_run(None, [{"pass": False}, {"pass": True}])
        self.assertEqual(growth, [1])
        self.assertEqual(result["status"], "verified")
        self.assertTrue(result["verified"])

    def test_failed_multiplier_fit_can_retry_with_more_data(self):
        result, growth = self.fake_run(None, [{"pass": True}], multiplier_failure=True)
        self.assertEqual(growth, [1])
        self.assertEqual(result["status"], "verified")

    def test_no_verified_candidate_at_cap_is_inconclusive(self):
        result, growth = self.fake_run(None, [{"pass": False}, {"pass": False}])
        self.assertEqual(result["status"], "inconclusive")
        self.assertFalse(result["verified"])
        self.assertIn("failed exact verification", result["reason"])

    def test_cli_returns_nonzero_for_inconclusive_reconstruction(self):
        with patch.object(pipeline, "run", return_value={"status": "inconclusive",
                                                        "reason": "sample cap"}), \
             patch.object(sys, "argv", ["pipeline"]), \
             contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(pipeline.main(), 3)

    def test_unavailable_modal_shape_grows_to_cap_and_is_inconclusive(self):
        result, growth = self.fake_run(None, unstable_shapes=True)
        self.assertEqual(growth, [1])
        self.assertEqual(result["status"], "inconclusive")
        self.assertIn("modal proof shape", result["reason"])

    def test_defining_relation_nonzero_against_truncated_basis_is_inconclusive(self):
        # It belongs to the defining ideal by construction. An empty, truncated basis
        # leaves the relation unchanged, so its remainder cannot prove non-membership.
        relation = Presentation(3).relations(NumericRing(1, 2, 0))[0]
        self.assertTrue(relation.terms)
        result, growth = self.fake_run(None, unresolved_target=relation)
        self.assertEqual(growth, [])
        self.assertEqual(result["status"], "inconclusive")
        self.assertIsNone(result["in_ideal"])
        self.assertEqual(result["remainder_terms"], len(relation.terms))
        self.assertEqual(result["basis_completeness"], "uncertified")

    def test_invalid_sampling_limits_fail_before_execution(self):
        for parameters in ({"start": 0}, {"step": 0}, {"maxpts": 0},
                           {"start": 6, "maxpts": 5}):
            with self.subTest(parameters=parameters), self.assertRaises(ValueError):
                pipeline.run(None, None, (1, 2, 0), **parameters)


if __name__ == "__main__":
    unittest.main(verbosity=2)
