"""Exact reconstruction examples, failure diagnostics and saved proof data."""
from __future__ import annotations

import contextlib
from fractions import Fraction as F
import io
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from zipper_algebra.cofactor import lift as L
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.proof import Basis, Proof
from zipper_algebra.polynomial import Poly


class Samples:
    """A small exact two-element proof graph with coefficients specified independently."""

    def __init__(self, points, function=lambda q, t: F(1), graph=False, stable=False):
        self.points = list(points)
        self.function = function
        self.graph = graph
        self.char = 0
        if stable:
            self.stable = lambda skeleton: self.points + self.points[:1]

    def basis(self, q, t):
        names = ["a", "b"] if self.graph else ["a"]
        return SimpleNamespace(key=names, index_by_key={k: i for i, k in enumerate(names)})

    def proof(self, q, t):
        terms = [[(self.function(q, t), (), "input", 0, ())]]
        if self.graph:
            terms.append([(F(q + t), (1,), "basis", 0, (2,)),
                          (F(q * t), (), "input", 1, ())])
        return SimpleNamespace(terms=terms)


def skeleton(graph=False, **overrides):
    result = dict(anchor=(1, 1, 0), steps=[{"basis_key": "b" if graph else "a"}])
    result.update(overrides)
    return SimpleNamespace(**result)


class TestExactFittingDetails(unittest.TestCase):
    def test_bivariate_rational_coefficients_and_fractional_locations(self):
        def polynomial(q, t):
            return F(1, 3) + 2 * q - F(1, 2) * t + q * t
        points = [((q, t), polynomial(q, t))
                  for q in (F(1, 2), F(3, 2), F(5, 2)) for t in (F(1, 3), F(2, 3))]
        result = L.fit_poly_qt(points, 1, 1)
        self.assertEqual(result, {(0, 0): F(1, 3), (0, 1): F(-1, 2),
                                  (1, 0): F(2), (1, 1): F(1)})
        self.assertEqual(L.poly_eval_qt(result, F(7, 2), F(4, 3)),
                         polynomial(F(7, 2), F(4, 3)))
        self.assertEqual(L.poly_eval_qt({}, 2, 3), 0)

    def test_rank_consistency_and_small_sample_constraints(self):
        for bound in (-1, True, F(1, 2)):
            with self.subTest(bound=bound), self.assertRaises(ValueError):
                L.fit_poly_qt([], bound, 0)
        self.assertIsNone(L.fit_poly_qt([((q, q), F(q)) for q in range(1, 7)], 1, 1))
        self.assertIsNone(L.fit_poly_qt([((1, 1), 2)], 1, 1))
        self.assertIsNone(L.fit_poly_qt([((1, 1), 2), ((2, 1), 3)], 0, 0))
        self.assertIsNone(L._fit_numerator([((1, 1), F(2))], 1, 1))
        data = [((q, 0), F(q)) for q in range(1, 7)]
        self.assertEqual(L._fit_numerator(data, 1, 0), {(0, 0): 0, (1, 0): 1})
        self.assertIsNone(L._fit_numerator(data, 0, 0))

    def test_denominator_poles_bounds_and_unsuccessful_candidates(self):
        q = AtomSet([Atom("q", {(1, 0): 1})], [1])
        self.assertIsNone(L._scale_by_denominator([((0, 2), F(3))], q, (1,)))
        points = [((i, 1), F(1, i)) for i in range(1, 7)]
        self.assertEqual(L.fit_rational(points, q, 0, 0), ((1,), {(0, 0): 1}))
        self.assertIsNone(L.fit_rational(points, q, 0, 0, candidates=[(0,)]))
        self.assertIsNone(L.fit_rational_known(points, q, (2,), 0, 0))
        self.assertIsNone(L.fit_rational_known([((0, 1), F(1))], q, (1,), 0, 0))
        self.assertIsNone(L.fit_rational_known(points, q, (0,), 0, 0))
        for bounds in ((-1, 0), (0, True), (1.5, 0)):
            for fitter in (lambda: L.fit_rational([], q, *bounds),
                           lambda: L.fit_rational_known([], q, (0,), *bounds)):
                with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                    fitter()

    def test_exact_polynomial_division_and_atom_cancellation(self):
        q_minus_t = {(1, 0): 1, (0, 1): -1}
        self.assertTrue(L._divisible({(2, 0): 1, (0, 2): -1}, q_minus_t))
        self.assertFalse(L._divisible({(2, 0): 1, (0, 2): 1}, q_minus_t))
        self.assertTrue(L._divisible({}, q_minus_t))
        self.assertFalse(L._divisible({(0, 0): 1}, q_minus_t))
        atoms = AtomSet([Atom("q-t", q_minus_t)], [1])
        self.assertFalse(L.is_reduced(atoms, (1,), {(1, 0): 2, (0, 1): -2}))
        self.assertTrue(L.is_reduced(atoms, (0,), {(1, 0): 2, (0, 1): -2}))
        self.assertEqual(L._poly_divmod([-1, 0, 1, 0], [-1, 1, 0]),
                         ([1, 1, 0], []))
        self.assertEqual(L._poly_divmod([1], [0, 1]), ([], [1]))
        self.assertEqual(L._poly_divmod([0, 0], [1]), ([0, 0], []))
        with self.assertRaises(ZeroDivisionError):
            L._poly_divmod([1], [0, 0])

    def test_specialization_both_directions_and_zero_leading_coefficients(self):
        polynomial = {(2, 0): 1, (1, 1): -2, (0, 2): 1}
        self.assertEqual(L._specialise(polynomial, "t", 3), [9, -6, 1])
        self.assertEqual(L._specialise(polynomial, "q", 2), [4, -4, 1])
        self.assertEqual(L._specialise({(2, 1): 1}, "t", 0), [0])
        self.assertEqual(L._specialise({}, "q", 2), [0])

    def test_line_denominator_limits_rejected_lines_and_increasing_powers(self):
        atoms = AtomSet([Atom("q", {(1, 0): 1}), Atom("t", {(0, 1): 1})], [2, 2])
        points = [((q, t), F(1, q * t * t)) for q in range(1, 9) for t in range(1, 9)]
        self.assertEqual(L.line_denominator(points, atoms, max_lines=1), (1, 2))
        self.assertIsNone(L.line_denominator(points, atoms, max_lines=0))
        self.assertIsNone(L.line_denominator([((q, 1), int(q == 6))
                                              for q in range(1, 7)], atoms))
        self.assertIsNone(L.line_denominator([((1, 1), 2)], atoms))
        # The reconstructed denominator is q+1, so q itself is not a divisor.
        self.assertEqual(L.line_denominator([((q, 2), F(1, q + 1))
                                              for q in range(1, 9)], atoms), (0, 0))


class TestProofStabilityAndLift(unittest.TestCase):
    GRID = [(q, t) for q in range(1, 5) for t in range(1, 4)]

    def test_closure_follows_basis_references_once_and_preserves_anchor_order(self):
        samples = Samples(self.GRID, graph=True)
        b, p = samples.basis(1, 1), samples.proof(1, 1)
        self.assertEqual(L.closure(["b", "a", "b"], p, b), ["a", "b"])
        self.assertEqual(L.term_structure(p.terms[1], b),
                         (((1,), "basis", "a", (2,)), ((), "input", 1, ())))
        self.assertEqual(L.holdout_size(8), 4)
        self.assertEqual(L.holdout_size(100), 10)

    def test_unstable_missing_short_or_different_proofs_are_excluded(self):
        anchor = SimpleNamespace(key=["a"], index_by_key={"a": 0})
        proof = SimpleNamespace(terms=[[(F(1), (), "input", 0, ())]])
        bases = {(1, 1): anchor, (2, 1): SimpleNamespace(key=[], index_by_key={}),
                 (3, 1): anchor, (4, 1): anchor}
        proofs = {(1, 1): proof, (2, 1): SimpleNamespace(terms=[]),
                  (3, 1): SimpleNamespace(terms=[]),
                  (4, 1): SimpleNamespace(terms=[[(F(1), (1,), "input", 0, ())]])}
        samples = SimpleNamespace(points=list(bases), basis=lambda q, t: bases[q, t],
                                  proof=lambda q, t: proofs[q, t])
        self.assertEqual(L.stable_points(["a"], proof, anchor, samples), [(1, 1)])

    def test_complete_lift_known_denominators_roundtrip_and_verbose_reports(self):
        samples = Samples(self.GRID, graph=True, stable=True)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            result = L.lift(skeleton(True, basis_keys=["a", "b"]), samples,
                            AtomSet([]), holdout=2,
                            known={("a", 0): (), ("b", 0): (), ("b", 1): ()},
                            max_pq=1, max_pt=1)
        self.assertEqual(result.keys(), ["a", "b"])
        self.assertEqual(result.samples_used, len(self.GRID))
        self.assertEqual(result.holdout, 2)
        self.assertEqual(result.timing["routes"], {"known": 3, "enum": 0})
        self.assertIn("fitting on 10, holding out 2", out.getvalue())
        self.assertEqual(result.entries["b"][0][4], "a")
        self.assertEqual(result.entries["b"][1][4], 1)
        self.assertEqual(L.poly_eval_qt(result.entries["b"][0][1], 11, 12), 23)
        self.assertEqual(L.poly_eval_qt(result.entries["b"][1][1], 11, 10), 110)
        encoded = result.to_json()
        self.assertEqual(encoded["entries"]["b"]["terms"][0]["ref_anchor_index"], 0)
        self.assertEqual(L.Lift.from_json(encoded), result)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lift.json"
            result.save(path)
            self.assertEqual(L.Lift.load(path), result)
        self.assertNotEqual(result, object())
        self.assertNotEqual(result, L.Lift(AtomSet([]), {}, 1, 1, (1, 1, 0), {}))

    def test_lift_equality_uses_the_ordered_atom_definitions(self):
        entries = {"a": [((1,), {(0, 0): F(1)}, (), "input", 0, ())]}
        def result(atoms):
            return L.Lift(atoms, entries, 8, 2, (1, 1, 0), {"a": 0})
        original = result(AtomSet([Atom("a", {(1, 0): 1})], [1], provenance="first fit"))
        changed = result(AtomSet([Atom("a", {(1, 0): 1, (0, 0): 1})], [1]))
        self.assertNotEqual(original, changed)
        same = result(AtomSet([Atom("a", {(1, 0): 1})], [5], provenance="independent fit"))
        self.assertEqual(original, same)
        self.assertNotEqual(original, result(AtomSet([Atom("renamed", {(1, 0): 1})], [1])))
        # Exponent positions have meaning: reordering the same named factors
        # changes the represented coefficients when the entries stay fixed.
        first = Atom("q", {(1, 0): 1})
        second = Atom("t", {(0, 1): 1})
        ordered_entries = {"a": [((1, 0), {(0, 0): F(1)}, (), "input", 0, ())]}
        self.assertNotEqual(
            L.Lift(AtomSet([first, second]), ordered_entries, 8, 2, (1, 1, 0), {"a": 0}),
            L.Lift(AtomSet([second, first]), ordered_entries, 8, 2, (1, 1, 0), {"a": 0}))

    def test_known_line_route_and_fallback_to_enumeration(self):
        points = [(q, 1) for q in range(1, 10)]
        samples = Samples(points, lambda q, t: F(1, q))
        atoms = AtomSet([Atom("q", {(1, 0): 1})], [1])
        result = L.lift(skeleton(), samples, atoms, holdout=1, known="lines",
                        max_pq=0, max_pt=0, verbose=False)
        self.assertEqual(result.timing["routes"], {"known": 1, "enum": 0})
        self.assertEqual(result.entries["a"][0][:2], ((1,), {(0, 0): 1}))
        # An unsuitable supplied denominator is rejected; enumeration still succeeds.
        fallback = L.lift(skeleton(), samples, atoms, holdout=1, known={("a", 0): (0,)},
                          max_pq=0, max_pt=0, verbose=False)
        self.assertEqual(fallback.timing["routes"], {"known": 0, "enum": 1})
        self.assertEqual(fallback, result)
        # No long line supplies a proposal, which also falls back to enumeration.
        scattered = Samples([(q, q) for q in range(1, 8)], lambda q, t: F(2))
        result = L.lift(skeleton(), scattered, AtomSet([]), holdout=1, known="lines",
                        max_pq=0, max_pt=0, verbose=False)
        self.assertEqual(result.timing["routes"], {"known": 0, "enum": 1})

    def test_validation_and_inconclusive_result_types(self):
        samples = Samples(self.GRID)
        for bad in (0, -1, True, F(1, 2)):
            with self.subTest(holdout=bad), self.assertRaises(ValueError):
                L.lift(skeleton(), samples, AtomSet([]), holdout=bad)
        for anchor in ((1, 1), (1, 1, 5), (1, 1, False)):
            with self.subTest(anchor=anchor), self.assertRaises(ValueError):
                L.lift(skeleton(anchor=anchor), samples, AtomSet([]))
        for char in (5, False):
            samples.char = char
            with self.assertRaises(ValueError):
                L.lift(skeleton(), samples, AtomSet([]))
        samples.char = 0
        for bounds in ((-1, 0), (0, True), (F(1, 2), 1)):
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                L.lift(skeleton(), samples, AtomSet([]), max_pq=bounds[0], max_pt=bounds[1])
        with self.assertRaisesRegex(ValueError, "closure"):
            L.lift(skeleton(basis_keys=["b"]), samples, AtomSet([]), verbose=False)
        with self.assertRaises(L.NeedMoreSamples):
            L.lift(skeleton(), Samples(self.GRID[:4]), AtomSet([]), verbose=False)
        with self.assertRaises(L.NeedMoreSamples) as failed:
            L.lift(skeleton(), Samples(self.GRID, lambda q, t: F(q)), AtomSet([]),
                   holdout=1, max_pq=0, max_pt=0, verbose=False)
        self.assertEqual(failed.exception.slots, [("a", 0)])
        error = L.WrongClass("inconclusive", [("a", 0)])
        self.assertIsInstance(error, L.NeedMoreSamples)
        self.assertEqual(error.slots, [("a", 0)])

    def test_holdout_rejects_wrong_prediction_and_a_pole(self):
        locations = [(q, 1) for q in range(1, 9)]
        with self.assertRaises(L.WrongClass) as failed:
            L.lift(skeleton(), Samples(locations, lambda q, t: F(2 if q == 8 else 1)),
                   AtomSet([]), holdout=1, max_pq=0, max_pt=0, verbose=False)
        self.assertEqual(failed.exception.slots, [("a", 0, "holdout")])
        # Sample order reserves q=0 for holdout, after all nonzero fitting locations.
        pole_samples = Samples(locations + [(0, 1)], lambda q, t: F(1, q) if q else F(0))
        with self.assertRaises(L.WrongClass):
            L.lift(skeleton(), pole_samples, AtomSet([Atom("q", {(1, 0): 1})], [1]),
                   holdout=1, max_pq=0, max_pt=0, verbose=False)

    def test_actual_atom_inference_and_its_verbose_diagnostics(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            result = L.lift(skeleton(), Samples([(q, 0) for q in range(1, 9)]),
                            holdout=2, max_pq=0, max_pt=0)
        self.assertEqual(result.entries["a"][0][1], {(0, 0): 1})
        self.assertEqual(result.atom_report["unexplained"], 0)
        silent = L.lift(skeleton(), Samples([(q, 0) for q in range(1, 9)]),
                        holdout=2, max_pq=0, max_pt=0, verbose=False)
        self.assertEqual(silent, result)
        self.assertIn("Inferring the atom set", output.getvalue())
        # A single line can expose a residual factor while offering too little
        # evidence to infer its bivariate origin. That is an inconclusive fit.
        samples = Samples([(q, 0) for q in range(1, 15)], lambda q, t: F(1, q**4 + 1))
        output = io.StringIO()
        with contextlib.redirect_stdout(output), self.assertRaises(L.NeedMoreSamples):
            L.lift(skeleton(), samples, holdout=2, max_pq=0, max_pt=0)
        self.assertIn("WARNING", output.getvalue())


class TestHistoricalSyntheticRegressions(unittest.TestCase):
    def test_stability_reidentifies_consistently_renumbered_proof_references(self):
        x, y = Poly.generator(1, 0), Poly.generator(2, 0)
        anchor = Basis([x, y], 0, ["x", "y"])
        reordered = Basis([y, x], 0, ["x", "y"])
        proof = Proof([[(F(5), (), "input", 0, ())],
                       [(F(7), (1,), "basis", 0, (2,))]])
        renumbered = Proof([[(F(13), (1,), "basis", 1, (2,))],
                            [(F(11), (), "input", 0, ())]])
        wrong_context = Proof([[(F(17), (2,), "basis", 1, (1,))],
                               [(F(19), (), "input", 0, ())]])
        bases = {(5, 7): anchor, (11, 13): reordered, (19, 17): reordered}
        proofs = {(5, 7): proof, (11, 13): renumbered, (19, 17): wrong_context}
        samples = SimpleNamespace(points=list(bases), basis=lambda q, t: bases[q, t],
                                  proof=lambda q, t: proofs[q, t])
        self.assertEqual(L.closure(["y"], proof, anchor), ["x", "y"])
        self.assertEqual(L.stable_points(["x", "y"], proof, anchor, samples),
                         [(5, 7), (11, 13)])

    def test_multifactor_rational_lift_preserves_multiplicity_and_known_route(self):
        points = [(q, t) for q in (5, 3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
                  for t in (7, 3, 5, 11, 13)]
        def coefficient(q, t):
            return F(q**3-t, q*(q-1)**2)
        samples = Samples(points, coefficient)
        skel = skeleton(anchor=(5, 7, 0))
        enumerated = L.lift(skel, samples, AtomSet.DEFAULT, verbose=False)
        known = L.lift(skel, samples, AtomSet.DEFAULT, known="lines", verbose=False)
        self.assertEqual(enumerated, known)
        self.assertEqual(enumerated.entries["a"],
                         [((1, 0, 2, 0, 0, 0), {(3, 0): F(1), (0, 1): F(-1)},
                           (), "input", 0, ())])
        self.assertEqual(known.timing["routes"], {"known": 1, "enum": 0})
        for q, t in ((43, 17), (47, 19)):
            den, num, *_ = known.entries["a"][0]
            self.assertEqual(L.poly_eval_qt(num, q, t)/known.atoms.eval_denom(den, q, t),
                             coefficient(q, t))
