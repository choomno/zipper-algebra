"""Exact arithmetic, signed basis grammar and the zero-dimensional unit quotient."""
from fractions import Fraction
import operator
from pathlib import Path
import tempfile
import unittest

from tools import core
from tools.zipper_in_ms import Coef, NumericRing, Poly, QTCoef


class TestCoefficientArithmetic(unittest.TestCase):
    def test_rationals_reduce_by_inverting_the_denominator(self):
        self.assertEqual(Coef(Fraction(1, 2), 5).v, 3)
        self.assertEqual(Coef(Fraction(-2, 3), 5).v, 1)
        self.assertEqual((Coef(Fraction(1, 2), 5) * Coef(2, 5)).v, 1)
        self.assertEqual((Coef(Fraction(1, 2), 0) + Coef(Fraction(1, 3), 0)).v,
                         Fraction(5, 6))

    def test_noninvertible_denominator_is_rejected(self):
        with self.assertRaises(ZeroDivisionError):
            Coef(Fraction(1, 5), 5)

    def test_mixed_characteristics_are_rejected_in_both_orders(self):
        for left, right in ((Coef(1, 0), Coef(2, 5)), (Coef(1, 5), Coef(2, 7))):
            for operation in (operator.add, operator.sub, operator.mul):
                for first, second in ((left, right), (right, left)):
                    with self.subTest(operation=operation, first=first.char, second=second.char):
                        with self.assertRaises(TypeError):
                            operation(first, second)
            self.assertNotEqual(left, right)

    def test_coefficients_do_not_coerce_inexact_inputs_or_invalid_characteristics(self):
        for value in (0.5, True, "1/2"):
            with self.subTest(value=value), self.assertRaises(TypeError):
                Coef(value, 0)
        for characteristic in (-1, 1, 3.0, True):
            with self.subTest(characteristic=characteristic), self.assertRaises(ValueError):
                Coef(1, characteristic)
            with self.subTest(characteristic=characteristic), self.assertRaises(ValueError):
                Poly(characteristic)

    def test_rational_constants_interoperate_with_symbolic_polynomials(self):
        constant = Coef(Fraction(1, 2), 0)
        symbol = QTCoef({(1, 0): 1})
        self.assertEqual((constant * symbol).v, {(1, 0): Fraction(1, 2)})
        self.assertEqual((symbol * constant).v, {(1, 0): Fraction(1, 2)})
        self.assertEqual(QTCoef(Fraction(1, 2)), constant)
        self.assertEqual(constant, QTCoef(Fraction(1, 2)))

    def test_symbolic_coefficients_obey_the_exact_polynomial_domain(self):
        for value in (0.5, True, "1/2", {(0, 0): 0.5}, {(0, 0): True}):
            with self.subTest(value=value), self.assertRaises(TypeError):
                QTCoef(value)
        for value in ({(-1, 0): 1}, {(0.0, 1): 1}, {(True, 0): 1}, {(0,): 1}):
            with self.subTest(value=value), self.assertRaises(ValueError):
                QTCoef(value)
        for operation in (operator.add, operator.sub, operator.mul):
            for left, right in ((QTCoef(1), 0.5), (0.5, QTCoef(1)),
                                (QTCoef(1), Coef(1, 5)), (Coef(1, 5), QTCoef(1))):
                with self.subTest(operation=operation), self.assertRaises(TypeError):
                    operation(left, right)
        self.assertNotEqual(QTCoef(1), Coef(1, 5))
        self.assertNotEqual(Coef(1, 5), QTCoef(1))

    def test_composite_modulus_is_described_as_a_residue_ring(self):
        ring = NumericRing(2, 3, 4)
        self.assertIn("Z/4Z", ring.describe())
        self.assertNotIn("GF", ring.describe())
        self.assertEqual((ring.q * ring.q).v, 0)

    def test_zero_polynomials_cannot_hide_a_parent_mismatch(self):
        left, right = Poly(0), Poly(5)
        for operation in (operator.add, operator.sub, operator.mul):
            for first, second in ((left, right), (right, left)):
                with self.subTest(operation=operation), self.assertRaises(TypeError):
                    operation(first, second)
        with self.assertRaises(TypeError):
            left.add_term((), Coef.zero(5))
        with self.assertRaises(TypeError):
            left.scalar_mul(Coef.zero(5))

    def test_polynomial_words_have_positive_integer_indices(self):
        for word in ((0,), (-1,), (True,), (1.0,), [1]):
            with self.subTest(word=word), self.assertRaises(ValueError):
                Poly(0).add_term(word, Coef.zero(0))


class TestBasisGrammar(unittest.TestCase):
    def parse(self, text, parser):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "basis.txt"
            path.write_text(text)
            return parser(str(path))

    def test_negative_bare_words_are_signs_not_variable_names(self):
        expected = [(('x1', 'x2'), [(-1, ('x1', 'x2')), (1, ('y1',)), (-3, ())])]
        text = "-x1*x2 + y1 - 3,\n"
        self.assertEqual(self.parse(text, core.parse_basis_Fp), expected)
        self.assertEqual(self.parse(text, core.parse_basis_Q), expected)
        self.assertEqual(self.parse(text, core.parse_leading_monomials), [('x1', 'x2')])

    def test_compact_and_spaced_rational_leading_terms_agree(self):
        expected = [(('x1',), [(Fraction(-2, 3), ('x1',)),
                               (Fraction(5, 7), ('y1',)), (Fraction(-1, 2), ())])]
        for text in ("-2/3*x1+5/7*y1-1/2\n", "-2/3*x1 + 5/7*y1 - 1/2,\n"):
            with self.subTest(text=text):
                self.assertEqual(self.parse(text, core.parse_basis_Q), expected)
                self.assertEqual(self.parse(text, core.parse_leading_monomials), [('x1',)])
        with self.assertRaises(ValueError):
            self.parse("1/2*x1\n", core.parse_basis_Fp)

    def test_nonzero_signed_rational_scalars_generate_the_unit_ideal(self):
        for text in ("-1\n", "-2/3,\n", "+5/7\n"):
            with self.subTest(text=text):
                self.assertEqual(self.parse(text, core.parse_leading_monomials), [()])
                self.assertEqual(self.parse(text, core.parse_basis_Q)[0][0], ())

    def test_zero_terms_are_not_leading_words(self):
        text = "0,\n0*x1 + x2,\n"
        self.assertEqual(self.parse(text, core.parse_basis_Q), [(('x2',), [(Fraction(1), ('x2',))])])
        self.assertEqual(self.parse(text, core.parse_leading_monomials), [('x2',)])
        cancelled = "x1 - x1 + x2,\n"
        self.assertEqual(self.parse(cancelled, core.parse_basis_Q),
                         [(('x2',), [(Fraction(1), ('x2',))])])
        self.assertEqual(self.parse(cancelled, core.parse_leading_monomials), [('x2',)])

    def test_malformed_terms_are_rejected(self):
        for text in ("3x1\n", "x1**x2\n", "x1+\n", "1/0*x1\n", "x1^2\n",
                     "1 2*x1\n", "x1 x2\n", "x1,,,\n", ",\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.parse(text, core.parse_basis_Q)

    def test_preserved_basis_and_leading_reader_agree(self):
        path = Path(__file__).resolve().parents[1] / "logs_and_traces/f4ncgb/basis-n4-Q.txt"
        parsed = core.parse_basis_Q(str(path))
        self.assertEqual(len(parsed), 137)
        self.assertEqual([leading for leading, _ in parsed], core.parse_leading_monomials(str(path)))
        self.assertTrue(all(isinstance(coefficient, Fraction)
                            for _, terms in parsed for coefficient, _ in terms))


class TestUnitIdeal(unittest.TestCase):
    def test_zero_residue_does_not_produce_a_sparse_vector_entry(self):
        delta, matches, patterns, _ = core.build_aho_corasick([], [])
        self.assertEqual(core.reduce_to_basis_Fp([(5, ())], {(): 0}, {},
                                                delta, matches, patterns, 5), {})
    def test_a_leading_coefficient_that_vanishes_modulo_the_characteristic_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "not invertible"):
            core.lm_to_tail_Fp([(("x1",), [(5, ("x1",)), (1, ())])], 5)

    def test_empty_leading_word_leaves_no_standard_monomials(self):
        delta, matches, patterns, _ = core.build_aho_corasick([()], ["x1"])
        self.assertEqual(core.enumerate_standard_monomials(delta, matches, ["x1"]), [])
        self.assertEqual(core.find_leading_match((), delta, matches, patterns), (0, 0, ()))
        self.assertEqual(core.find_leading_match(("x1",), delta, matches, patterns), (0, 0, ()))

    def test_every_polynomial_reduces_to_zero_in_the_unit_quotient(self):
        delta, matches, patterns, _ = core.build_aho_corasick([()], ["x1"])
        for text in ("-1\n", "2\n"):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "unit.basis"
                path.write_text(text)
                rational = core.parse_basis_Q(str(path))
                finite = core.parse_basis_Fp(str(path))
            self.assertEqual(core.reduce_to_basis_Q(
                [(Fraction(3, 2), ()), (Fraction(-7), ("x1", "x1"))], {},
                core.lm_to_tail_Q(rational), delta, matches, patterns), {})
            self.assertEqual(core.reduce_to_basis_Fp(
                [(3, ()), (2, ("x1",))], {}, core.lm_to_tail_Fp(finite, 5),
                delta, matches, patterns, 5), {})
