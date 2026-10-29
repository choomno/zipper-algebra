"""Exact denominator arithmetic and small functional reconstruction examples."""
from contextlib import redirect_stdout
from fractions import Fraction as F
from io import StringIO
import importlib.util
import random
import runpy
import sys
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from zipper_algebra.cofactor import atoms as a
from zipper_algebra.cofactor.inference import infer
from zipper_algebra.cofactor import inference as i
from zipper_algebra.cofactor import univariate as u

HAS_SYMPY = importlib.util.find_spec("sympy") is not None


def grid_values(denominator, qs=range(2, 10), ts=range(11, 19)):
    return [( (q, t), F(1, denominator(q, t))) for q in qs for t in ts]


class TestAtomArithmetic(unittest.TestCase):
    def test_integer_polynomial_validation(self):
        for name in (None, "", " ", "1"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                a.Atom(name, {(0, 0): 1})
        for coefs in ({(1,): 1}, {(1, -1): 1}, {(True, 0): 1},
                      {("q", 0): 1}, {(1, 0): True}, {(1, 0): 1.0},
                      {(1, 0): F(1, 2)}, {}):
            with self.subTest(coefs=coefs), self.assertRaises(ValueError):
                a.Atom("bad", coefs)
        atom = a.Atom("q", {(1, 0): F(1), (0, 0): 0})
        self.assertEqual(atom.coefs, {(1, 0): 1})
        atom.coefs = []
        with self.assertRaises(ValueError):
            atom.validate()

    def test_atom_labels_reject_unbalanced_parentheses(self):
        for name in ("q)", "q(", ")("):
            with self.subTest(name=name), self.assertRaises(ValueError):
                a.Atom(name, {(1, 0): 1})

    def test_balanced_atom_labels_round_trip_as_denominator_factors(self):
        for name in ("g(t)", "(q*t+1)", "q^2", "q*t"):
            family = a.AtomSet([a.Atom(name, {(1, 0): 1, (0, 1): 1})])
            for exponent in (1, 2):
                with self.subTest(name=name, exponent=exponent):
                    written = family.denom_str((exponent,))
                    self.assertEqual(family.parse_denom_str(written), (exponent,))
                    self.assertEqual(family.eval_denom((exponent,), 2, 3), 5**exponent)

    def test_polynomial_values_persistence_and_sign_units(self):
        atom = a.Atom("difference", {(1, 0): 1, (0, 1): -1})
        negative = a.Atom("opposite", {(1, 0): -1, (0, 1): 1})
        self.assertEqual(atom, negative)
        self.assertEqual(hash(atom), hash(negative))
        self.assertNotEqual(atom, "q-t")
        self.assertEqual(atom.eval(F(3, 2), F(1, 2)), 1)
        self.assertEqual(atom.bidegree(), (1, 1))
        self.assertEqual(a.Atom.from_json(atom.to_json()), atom)
        self.assertEqual(repr(atom), "Atom(difference)")
        self.assertEqual(a.poly_str({(2, 3): 2, (1, 1): -1, (0, 1): 3,
                                     (0, 0): -4}), "2*q^2*t^3-q*t+3*t-4")
        self.assertEqual(a.poly_str({}), "0")
        self.assertEqual(a.poly_str({(1, 0): 1, (0, 2): -1}), "q-t^2")

    def test_atom_sets_preserve_denominator_positions(self):
        q = a.Atom("q", {(1, 0): 1})
        t = a.Atom("t", {(0, 1): 1})
        compound = a.Atom("q*t+1", {(1, 1): 1, (0, 0): 1})
        family = a.AtomSet([q, t, compound], [2, 1, 2], "initial")
        self.assertEqual(len(family), 3)
        self.assertEqual(family.zero(), (0, 0, 0))
        self.assertEqual(family.eval_denom((2, 0, 1), 2, 3), 28)
        self.assertEqual(family.denom_str((2, 0, 1)), "q^2*(q*t+1)")
        self.assertEqual(family.denom_str((0, 1, 2)), "t*(q*t+1)^2")
        self.assertEqual(family.denom_str((0, 0, 0)), "1")
        self.assertEqual(family.parse_denom_str("q*q*(q*t+1)"), (2, 0, 1))
        self.assertEqual(family.parse_denom_str("q^2*(q*t+1)^2"), (2, 0, 2))
        self.assertEqual(family.parse_denom_str("q^0"), (0, 0, 0))
        self.assertEqual(family.parse_denom_str(""), (0, 0, 0))
        self.assertEqual(family.parse_denom_str("1"), (0, 0, 0))
        self.assertEqual(a.AtomSet.lcm((2, 0, 1), (1, 3, 0)), (2, 3, 1))
        self.assertEqual(a.AtomSet.mul((2, 0, 1), (1, 3, 0)), (3, 3, 1))
        self.assertEqual(a.AtomSet([q], [2]).candidates(), [(0,), (1,), (2,)])
        extended = family.extended([a.Atom("minus_q", {(1, 0): -1}),
                                    a.Atom("q+1", {(1, 0): 1, (0, 0): 1})],
                                   4, "extended")
        self.assertEqual(extended.names(), ["q", "t", "q*t+1", "q+1"])
        self.assertEqual(extended.max_exp, [2, 1, 2, 4])
        self.assertEqual(extended.provenance, "extended")
        self.assertEqual(family.extended([]).provenance, "initial")
        self.assertIn("box=[2, 1, 2]", repr(family))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "atoms.json"
            family.save(path)
            restored = a.AtomSet.load(path)
        self.assertEqual(restored.to_json(), family.to_json())
        self.assertEqual(a.AtomSet.from_json({"atoms": [q.to_json()]}).max_exp, [3])

    def test_polynomial_display_uses_coefficients_and_preserves_historical_spelling(self):
        defaults = a.AtomSet.DEFAULT
        exponents = (2, 1, 2, 1, 2, 1)
        expected = "q^2*t*(q-1)^2*(t-q)*(t+q^2)^2*(q+1)"
        self.assertEqual(defaults.denom_polynomial_str(exponents), expected)
        self.assertEqual(defaults.denom_polynomial_str(exponents), defaults.denom_str(exponents))
        self.assertEqual(defaults.denom_polynomial_str(defaults.zero()), "1")
        self.assertEqual(a.AtomSet([]).denom_polynomial_str(()), "1")
        family = a.AtomSet([a.Atom("t+5", {(1, 0): 1, (0, 0): -1}),
                            a.Atom("q", {(1, 1): -2, (0, 0): 3}),
                            a.Atom("constant", {(0, 0): -2})])
        self.assertEqual(family.denom_polynomial_str((1, 2, 1)),
                         "(q-1)*(-2*q*t+3)^2*(-2)")
        self.assertEqual(family.denom_str((1, 2, 1)), "(t+5)*q^2*constant")
        self.assertEqual(family.eval_denom((1, 2, 1), 2, 3), -162)
        with self.assertRaises(ValueError):
            family.denom_polynomial_str((1,))

    def test_malformed_families_and_denominators_are_rejected(self):
        q = a.Atom("q", {(1, 0): 1})
        for values in (["q"], [q, q]):
            with self.assertRaises(ValueError):
                a.AtomSet(values)
        family = a.AtomSet([q])
        for exponents in ((-1,), (True,), (1.0,), (), "1"):
            with self.subTest(exponents=exponents), self.assertRaises(ValueError):
                family.validate_exponents(exponents)
        for pair in [((1,), (1, 2)), ((-1,), (0,)), ("1", (1,))]:
            with self.subTest(pair=pair), self.assertRaises(ValueError):
                a.AtomSet.lcm(*pair)
        for text in ("q)", "(q", "(q)^-1", "q^-1", "unknown", "q**q"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                family.parse_denom_str(text)


class TestUnivariateDenominatorArithmetic(unittest.TestCase):
    def test_division_exact_and_with_remainder(self):
        quotient, remainder = u._poly_divmod([-1, 0, 1, 0], [-1, 1, 0])
        self.assertEqual(quotient, [1, 1, 0])
        self.assertEqual(remainder, [])
        self.assertEqual(u._poly_divmod([1], [1, 1]), ([], [1]))
        self.assertEqual(u._poly_divmod([1, 0, 1], [1, 1]), ([-1, 1], [2]))
        with self.assertRaises(ZeroDivisionError):
            u._poly_divmod([1], [0, 0])
        self.assertEqual(u._poly_eval([F(1, 2), 2, 3], 2), F(33, 2))

    def test_linear_systems_detect_missing_information(self):
        self.assertEqual(u._solve([[F(0), F(1), F(3)],
                                  [F(1), F(0), F(2)]], 2), [2, 3])
        self.assertIsNone(u._solve([[F(0), F(1)]], 1))
        self.assertIsNone(u._solve([[F(0), F(0)]], 1))
        self.assertEqual(u._solve([[F(1), F(2)], [F(2), F(4)]], 1), [2])

    def test_samples_are_exact_and_distinct(self):
        for value in (True, 0.5, "1"):
            with self.assertRaises(ValueError):
                u.exact_fraction(value)
        self.assertEqual(u.normalise_points([(2, 4), (1, 1), (2, 4)], 1),
                         [(F(1), F(1)), (F(2), F(4))])
        self.assertEqual(u.normalise_points([((1, 2), 3)], 2),
                         [((F(1), F(2)), F(3))])
        for points, dimension in [([1], 1), ([(1, 2, 3)], 1),
                                  ([(1, 2)], 2), ([((1,), 2)], 2),
                                  ([(1, 2), (1, 3)], 1)]:
            with self.subTest(points=points), self.assertRaises(ValueError):
                u.normalise_points(points, dimension)

    def test_fraction_reduction_and_rational_reconstruction(self):
        self.assertEqual(u._reduce_univariate_fraction([0, 0], [2, 2]), ([0], [1]))
        self.assertEqual(u._reduce_univariate_fraction([-2, 0, 2], [-2, 2]),
                         ([1, 1], [1]))
        with self.assertRaises(ValueError):
            u._reduce_univariate_fraction([1], [0])
        for holdout in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                u.reconstruct_univariate([], holdout)
        self.assertIsNone(u.reconstruct_univariate([(1, 1)]))
        points = [(x, F(x + 2, x - 1)) for x in range(2, 11)]
        numerator, denominator = u.reconstruct_univariate(points)
        self.assertEqual(numerator, [2, 1])
        self.assertEqual(denominator, [-1, 1])
        self.assertEqual(u.reconstruct_univariate([(x, 0) for x in range(7)]), ([0], [1]))
        self.assertIsNone(u.reconstruct_univariate([(x, int(x == 6)) for x in range(7)]))

    def test_specialization_and_candidate_order(self):
        polynomial = {(2, 1): 1, (1, 0): 2, (0, 2): -1}
        self.assertEqual(u._specialise(polynomial, "t", 3), [-9, 2, 3])
        self.assertEqual(u._specialise(polynomial, "q", 2), [4, 4, -1])
        self.assertEqual(u._specialise({(1, 1): 1, (1, 0): -1}, "t", 1), [0])
        self.assertEqual(u._specialise({}, "q", 2), [0])
        candidates = i._candidate_polys((1, 0), 2)
        expected = {"q", "q-1", "q+1", "q-2", "q+2", "2*q-1", "2*q+1"}
        self.assertEqual({a.poly_str(p) for p in candidates}, expected)
        self.assertEqual([sum(abs(c) for c in p.values()) for p in candidates],
                         sorted(sum(abs(c) for c in p.values()) for p in candidates))
        self.assertEqual(i._candidate_polys((0, 0), 2), [])

    def test_primitive_and_readable_univariate_polynomials(self):
        self.assertEqual(i._primitive([]), ())
        self.assertEqual(i._primitive([0, 0]), (0, 0))
        self.assertEqual(i._primitive([F(1, 2), -1]), (-1, 2))
        self.assertEqual(i._fmt_univariate([0, -1, 2, -3], "q"), "-q+2*q^2-3*q^3")
        self.assertEqual(i._fmt_univariate([0], "q"), "0")
        self.assertEqual(i._fmt_univariate([2, 1], "t"), "2+t")

    def test_polynomial_line_fits_have_independent_holdouts(self):
        self.assertEqual(i._fit_univariate_exact([(x, x*x-1) for x in range(7)], 3),
                         {0: F(-1), 2: F(1)})
        self.assertIsNone(i._fit_univariate_exact([(0, 1), (1, 2)], 3))
        self.assertIsNone(i._fit_univariate_exact([(1, 0), (1, 0), (1, 0), (2, 1), (3, 2)], 4))
        self.assertIsNone(i._fit_univariate_exact([(x, x*x) for x in range(5)], 1))


class TestAtomInference(unittest.TestCase):
    def test_empty_and_unreconstructible_lines(self):
        with redirect_stdout(StringIO()) as output:
            family, report = infer({}, max_bideg=(1, 0), coef_bound=1)
        self.assertEqual(family.names(), [])
        self.assertEqual(report["lines_tried"], 0)
        self.assertIn("every reconstructed denominator", output.getvalue())
        points = [((x, 1), int(x == 6)) for x in range(1, 7)]
        family, report = infer({"bad": points}, max_bideg=(1, 0), coef_bound=1, verbose=False)
        self.assertEqual(report["lines_tried"], 1)
        self.assertEqual(report["lines_reconstructed"], 0)

    def test_polynomial_coefficients_need_no_denominator_atoms(self):
        points = [((q, t), q + t) for q in range(2, 9) for t in range(11, 18)]
        family, report = infer({0: points}, max_bideg=(1, 0), coef_bound=1,
                                  verbose=False)
        self.assertEqual(family.names(), [])
        self.assertEqual(report["lines_reconstructed"], 14)
        self.assertEqual(report["unexplained"], 0)

    def test_quiet_search_has_the_same_exact_factor_result(self):
        family, report = infer({0: grid_values(lambda q, t: q+t)},
                                  max_bideg=(1, 1), coef_bound=1, verbose=False)
        self.assertEqual(family.atoms, [a.Atom("q+t", {(1, 0): 1, (0, 1): 1})])
        self.assertEqual(report["unexplained"], 0)

    def test_seed_divisibility_and_multiplicity(self):
        seed = a.AtomSet([a.Atom("q+t", {(1, 0): 1, (0, 1): 1})], [0])
        family, report = infer({0: grid_values(lambda q, t: (q+t)**2)}, seed,
                                  max_bideg=(1, 1), coef_bound=1, verbose=False)
        self.assertEqual(family.names(), ["q+t"])
        self.assertEqual(report["multiplicity"], {"q+t": 2})
        self.assertEqual(family.max_exp, [3])
        self.assertEqual(report["unexplained"], 0)

    def test_greedy_search_finds_independent_factors(self):
        with redirect_stdout(StringIO()) as output:
            family, report = infer({0: grid_values(lambda q, t: (q+t)*(q-t))},
                                      max_bideg=(1, 1), coef_bound=1)
        self.assertEqual(set(family.atoms), {a.Atom("q+t", {(1, 0): 1, (0, 1): 1}),
                                            a.Atom("q-t", {(1, 0): 1, (0, 1): -1})})
        self.assertEqual(report["unexplained"], 0)
        self.assertIn("found", output.getvalue())

    def test_one_line_cannot_support_a_bivariate_factor(self):
        with mock.patch.dict(sys.modules, {"sympy": None}):
            family, report = infer({0: grid_values(lambda q, t: q, range(1, 8), [0])},
                                      max_bideg=(1, 1), coef_bound=1, min_lines=1, verbose=False)
        self.assertEqual(family.names(), [])
        self.assertEqual(report["unexplained"], 1)

    def test_unavailable_factorization_remains_inconclusive(self):
        points = grid_values(lambda q, t: q+t+5, range(1, 8), range(1, 4))
        with mock.patch.dict(sys.modules, {"sympy": None}):
            with redirect_stdout(StringIO()) as output:
                family, report = infer({0: points}, max_bideg=(0, 0), coef_bound=1)
        self.assertEqual(family.names(), [])
        self.assertEqual(report["unexplained"], 3)
        self.assertEqual(report["residual_factors"], [])
        self.assertIn("factorisation unavailable", output.getvalue())
        self.assertIn("WIDEN THE CLASS", output.getvalue())

    def test_line_factor_lifting_and_ambiguity(self):
        factors = {("t", (x, 1)): [x] for x in range(1, 8)}
        factors.update({("q", (x, 1)): [x] for x in range(1, 8)})
        with redirect_stdout(StringIO()) as output:
            result = i.atoms_from_residual_factors(factors)
        self.assertEqual(result, [a.Atom("q+t", {(1, 0): 1, (0, 1): 1})])
        self.assertIn("lifted to a bivariate atom", output.getvalue())
        odd = range(1, 12, 2)
        rational = {("t", ((x+1)//2, 1)): [x] for x in odd}
        self.assertEqual(i.atoms_from_residual_factors(rational, verbose=False),
                         [a.Atom("2*q+t+1", {(1, 0): 2, (0, 1): 1, (0, 0): 1})])
        self.assertEqual(i.atoms_from_residual_factors({("t", (1, 1)): [1]}, verbose=False), [])
        conflicting = {("t", (1, 1)): list(range(7)), ("t", (2, 1)): list(range(7))}
        self.assertEqual(i.atoms_from_residual_factors(conflicting, verbose=False), [])
        nonpolynomial = {("t", (2**x, 1)): [x] for x in range(7)}
        self.assertEqual(i.atoms_from_residual_factors(nonpolynomial, max_other_deg=1, verbose=False), [])
        constant = {("t", (1,)): list(range(7))}
        self.assertEqual(i.atoms_from_residual_factors(constant, verbose=False), [])

    def test_standalone_example_rediscovers_its_synthetic_denominator_family(self):
        state = random.getstate()
        try:
            with redirect_stdout(StringIO()) as output:
                runpy.run_path(Path(__file__).resolve().parents[2] / "examples/infer_denominators.py",
                               run_name="__main__")
        finally:
            random.setstate(state)
        self.assertIn("rediscovered the six from denominators", output.getvalue())

    @unittest.skipUnless(HAS_SYMPY, "optional Sympy inference dependency")
    def test_actual_residual_factorization(self):
        # (x-1)^2*(x^2+1) has exactly these two irreducible rational factors.
        lines = [("t", 2, [1, -2, 2, -2, 1]), ("q", 3, [-1, 1]), ("t", 4, [1])]
        with redirect_stdout(StringIO()) as output:
            factors = i.residual_factors(lines)
        self.assertEqual(factors, {("t", (-1, 1)): [2], ("t", (1, 0, 1)): [2],
                                   ("q", (-1, 1)): [3]})
        self.assertIn("residual irreducible factors", output.getvalue())
        self.assertEqual(i.residual_factors([], verbose=False), {})

    @unittest.skipUnless(HAS_SYMPY, "optional Sympy inference dependency")
    def test_missing_large_constant_is_identified_by_residual_lifting(self):
        points = grid_values(lambda q, t: q+t+5)
        with redirect_stdout(StringIO()):
            family, report = infer({0: points}, max_bideg=(1, 1), coef_bound=1)
        self.assertEqual(family.atoms, [a.Atom("q+t+5", {(1, 0): 1, (0, 1): 1, (0, 0): 5})])
        self.assertEqual(report["unexplained"], 0)
        self.assertTrue(report["residual_factors"])


if __name__ == "__main__":
    unittest.main()
