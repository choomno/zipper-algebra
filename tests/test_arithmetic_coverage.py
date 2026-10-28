"""Algebraic laws, presentation variants and exact solver-input emission."""
from contextlib import redirect_stdout
from fractions import Fraction
import io
import operator
from pathlib import Path
import runpy
import sys
import tempfile
import unittest
from unittest.mock import patch

from zipper_algebra import zipper_in_ms as z
from zipper_algebra.cofactor import ncpoly as nc


class TestScalarAndFreePolynomialArithmetic(unittest.TestCase):
    def test_exact_scalar_signs_equality_and_rendering(self):
        half = z.Coef(Fraction(1, 2), 0)
        third = z.Coef(Fraction(1, 3), 0)
        self.assertEqual((half - third).v, Fraction(1, 6))
        self.assertEqual((-half).v, Fraction(-1, 2))
        self.assertEqual(str(half), '1/2')
        self.assertEqual(str(z.Coef(2, 0)), '2')
        self.assertEqual(str(z.Coef(-1, 5)), '4')
        self.assertTrue(z.Coef.zero(5).is_zero())
        self.assertFalse(half.is_zero())
        self.assertNotEqual(half, object())
        for operation in (operator.add, operator.sub, operator.mul):
            with self.subTest(operation=operation), self.assertRaises(TypeError):
                operation(half, object())

    def test_symbolic_arithmetic_with_rational_constants_and_reflection(self):
        q = z.QTCoef({(1, 0): 1})
        t = z.QTCoef({(0, 1): 1})
        self.assertEqual(((q + t) * (q - t)).v,
                         {(2, 0): Fraction(1), (0, 2): Fraction(-1)})
        self.assertEqual((2 - q).v, {(0, 0): Fraction(2), (1, 0): Fraction(-1)})
        self.assertEqual((z.Coef(2, 0) - q).v,
                         {(0, 0): Fraction(2), (1, 0): Fraction(-1)})
        self.assertEqual((z.Coef(2, 0) + q).v,
                         {(0, 0): Fraction(2), (1, 0): Fraction(1)})
        self.assertEqual((Fraction(1, 2) * q).v, {(1, 0): Fraction(1, 2)})
        self.assertEqual((q + Fraction(1, 3)).v,
                         {(1, 0): Fraction(1), (0, 0): Fraction(1, 3)})
        self.assertTrue((q - q).is_zero())
        self.assertEqual((q + t * t).bidegree(), (1, 2))
        self.assertEqual(z.QTCoef(0).bidegree(), (0, 0))
        self.assertNotEqual(q, object())
        self.assertNotEqual(q, z.Coef(1, 5))
        for operation in (operator.add, operator.sub, operator.mul):
            for left, right in ((q, object()), (object(), q)):
                with self.subTest(operation=operation), self.assertRaises(TypeError):
                    operation(left, right)
        with self.assertRaises(TypeError):
            q + z.Coef(1, 7)

    def test_symbolic_rendering_distinguishes_powers_constants_and_rationals(self):
        self.assertEqual(str(z.QTCoef(0)), '0')
        polynomial = z.QTCoef({(0, 0): -2, (1, 0): 1, (0, 1): 3,
                               (2, 3): -4})
        self.assertEqual(str(polynomial), '(-2)+(3)*t+(1)*q+(-4)*q^2*t^3')
        with self.assertRaisesRegex(ValueError, 'non-integer coefficient'):
            str(z.QTCoef(Fraction(1, 2)))
        self.assertEqual(z.NumericRing(2, 3, 0).describe(), '(q,t)=(2,3) over Q')
        self.assertEqual(z.NumericRing(2, 3, 5).scalar(Fraction(1, 2)).v, 3)
        symbolic = z.SymbolicRing()
        self.assertEqual(symbolic.scalar(2).v, {(0, 0): Fraction(2)})
        self.assertEqual(symbolic.describe(), 'q,t symbolic over Q')

    def test_free_products_preserve_order_and_cancel_coefficient_collisions(self):
        x, y = z.Poly.generator(1, 0), z.Poly.generator(2, 0)
        square = (x + y) * (x + y)
        mixed = square - x*x - y*y
        self.assertEqual({w: c.v for w, c in mixed.terms.items()},
                         {(1, 2): Fraction(1), (2, 1): Fraction(1)})
        self.assertEqual({w: c.v for w, c in (-mixed).terms.items()},
                         {(1, 2): Fraction(-1), (2, 1): Fraction(-1)})
        self.assertFalse((mixed - mixed).terms)
        self.assertEqual({w: c.v for w, c in mixed.scalar_mul(z.Coef(3, 0)).terms.items()},
                         {(1, 2): Fraction(3), (2, 1): Fraction(3)})
        with self.assertRaises(TypeError):
            mixed.scalar_mul(object())
        with self.assertRaises(TypeError):
            mixed.add_term((1,), z.Coef(1, 5))
        for operation in (operator.add, operator.sub, operator.mul):
            with self.subTest(operation=operation), self.assertRaises(TypeError):
                operation(x, object())
        p = z.Poly(0)
        p.add_term((1,), z.Coef(2, 0))
        p.add_term((1,), z.Coef(3, 0))
        p.add_term((), z.Coef(0, 0))
        self.assertEqual(p.terms[(1,)].v, 5)
        self.assertEqual(str(p.terms[(1,)]), '5')
        self.assertEqual((z.Poly.scalar(z.Coef(0, 0))).terms, {})

    def test_msolve_rendering_is_unambiguous_about_signs_and_empty_words(self):
        p = z.Poly(0)
        for word, coefficient in [((1, 2), 1), ((2, 1), -1), ((1,), 2), ((), -3)]:
            p.add_term(word, z.Coef(coefficient, 0))
        self.assertEqual(p.to_msolve_string(['x', 'y']), 'x*y-y*x+2*x-3')
        self.assertEqual(z.Poly(0).to_msolve_string(['x']), '0')
        self.assertEqual(z.Poly.scalar(z.Coef(2, 0)).to_msolve_string(['x']), '2')


class TestSparsePolynomialLaws(unittest.TestCase):
    def test_commuting_parameter_arithmetic_and_text(self):
        q, t = nc.qt_q(), nc.qt_t()
        self.assertEqual(nc.qt_zero(), {})
        self.assertEqual(nc.qt_monomial(2, 3, 0), {})
        self.assertEqual(nc.qt_monomial(2, 3, Fraction(1, 2)),
                         {(2, 3): Fraction(1, 2)})
        self.assertEqual(nc.qt_mul(nc.qt_add(q, t), nc.qt_sub(q, t)),
                         {(2, 0): Fraction(1), (0, 2): Fraction(-1)})
        self.assertEqual(nc.qt_sub(q, q), {})
        self.assertEqual(nc.qt_add(q, nc.qt_neg(q)), {})
        self.assertEqual(nc.qt_sub(q, t), {(1, 0): Fraction(1), (0, 1): Fraction(-1)})
        self.assertEqual(nc.qt_pow(q, 0), {(0, 0): Fraction(1)})
        self.assertEqual(nc.qt_pow(q, 2), {(2, 0): Fraction(1)})
        self.assertTrue(nc.qt_eq(q, {(1, 0): Fraction(1)}))
        self.assertFalse(nc.qt_eq(q, t))
        self.assertTrue(nc.qt_is_zero(nc.qt_const(0)))
        self.assertFalse(nc.qt_is_zero(q))
        self.assertEqual(nc.qt_str({}), '0')
        self.assertEqual(nc.qt_str({(0, 0): Fraction(-1, 2), (0, 1): Fraction(-1),
                                    (1, 0): Fraction(1), (2, 0): Fraction(3, 2)}),
                         '-1/2 - t^1 + q^1 + 3/2*q^2')

    def test_noncommuting_polynomials_and_zero_coefficient_normalization(self):
        one = {(0, 0): Fraction(1)}
        minus = {(0, 0): Fraction(-1)}
        x, y = nc.nc_var(1), nc.nc_var_named('y', {'y': 2})
        self.assertEqual(nc.nc_zero(), {})
        self.assertEqual(nc.nc_const({}), {})
        self.assertEqual(nc.nc_const(one), {(): one})
        with self.assertRaises(KeyError):
            nc.nc_var_named('z', {'x': 1})
        self.assertEqual(nc.nc_sub(nc.nc_add(x, y), x), y)
        self.assertEqual(nc.nc_sub(x, y), {(1,): one, (2,): minus})
        self.assertEqual(nc.nc_add(x, nc.nc_neg(x)), {})
        self.assertEqual(nc.nc_sub(x, x), {})
        self.assertEqual(nc.nc_add(x, x), {(1,): {(0, 0): Fraction(2)}})
        self.assertEqual(nc.nc_sub(nc.nc_scalar(x, nc.qt_const(2)), x), x)
        self.assertFalse(nc.nc_is_zero(x))
        self.assertTrue(nc.nc_is_zero(nc.nc_zero()))
        self.assertEqual(nc.nc_size(nc.nc_add(x, y)), 2)
        self.assertEqual(nc.nc_mul(nc.nc_add(x, y), nc.nc_sub(x, y)),
                         {(1, 1): one, (1, 2): minus, (2, 1): one, (2, 2): minus})
        # Cross terms of (x+x^2)(x-x^2) cancel, despite arising from different factors.
        left = {(1,): one, (1, 1): one}
        right = {(1,): one, (1, 1): minus}
        self.assertEqual(nc.nc_mul(left, right), {(1, 1): one, (1, 1, 1, 1): minus})
        self.assertEqual(nc.nc_scalar(x, {}), {})
        self.assertEqual(nc.nc_scalar(x, nc.qt_const(2)), {(1,): {(0, 0): Fraction(2)}})
        # Explicit zero coefficients represent the zero polynomial as well.
        self.assertEqual(nc.nc_mul({(1,): {}}, y), {})
        self.assertEqual(nc.nc_scalar({(1,): {}}, one), {})


class TestPresentationBuildersAndEmitter(unittest.TestCase):
    def test_all_relation_labels_encode_anchor_and_variant(self):
        cases = [
            (z.RelTag('INV', 2, None, 'm'), 5, 'INV_2_m'),
            (z.RelTag('BR', 1, 2, 'y'), 5, 'BR_1_y'),
            (z.RelTag('BR', 1, 2, 'x'), 5, 'BR_1'),
            (z.RelTag('FC', 1, 3, 'xy'), 4, 'FC_xy'),
            (z.RelTag('FC', 1, 4, 'xy'), 5, 'FC_14_xy'),
            (z.RelTag('MIX', 1, 2, 'a'), 5, 'MIX_1_a'),
            (z.RelTag('Z', 3, 1), 5, 'Z_3'),
            (z.RelTag('T', 3, 2), 5, 'T_3_j2'),
            (z.RelTag('EXTRA', 3, 2), 5, 'EXTRA_3_2'),
            (z.RelTag('EXTRA', 3, 2, 'power'), 5, 'EXTRA_3_2_power'),
        ]
        for tag, n, expected in cases:
            with self.subTest(expected=expected):
                self.assertEqual(tag.label(n), expected)
        self.assertEqual(repr(z.RelTag('T', 3, 1)), "RelTag(T,3,1,'')")

    def test_slider_and_two_strand_relations_against_explicit_coefficients(self):
        _, _, sliders, _ = z.build_relations_ring(3, z.NumericRing(5, 7, 0), 'bigelow')
        self.assertEqual({w: c.v for w, c in sliders[(3, 1)].terms.items()}, {
            (3, 4, 3): Fraction(125), (3, 4): Fraction(-100),
            (3, 4, 1): Fraction(-25), (1, 2, 3): Fraction(-5),
            (1, 2): Fraction(4), (1, 2, 1): Fraction(1),
        })
        modes = {
            'bigelow': {(1, 3): Fraction(5), (1, 1): Fraction(-1),
                        (1,): Fraction(3), (3,): Fraction(-35), (): Fraction(28)},
            'silly_question': {(1, 1, 1): Fraction(-1), (1, 1, 3): Fraction(5),
                               (1, 1): Fraction(-4), (3,): Fraction(-245),
                               (1,): Fraction(49), (): Fraction(196)},
            'square-zero': {(1, 1, 1): Fraction(-1), (1, 1, 3): Fraction(5),
                            (1, 1): Fraction(10), (1, 3): Fraction(-70),
                            (1,): Fraction(7), (3,): Fraction(245), (): Fraction(-196)},
        }
        # Keep inverse letters literal: no x*y cancellation is permitted in the emitter.
        for mode, expected in modes.items():
            _, relations, _, mode_tags = z.build_relations_ring(3, z.NumericRing(5, 7, 0), mode)
            rel = relations[next(i for i, tag in enumerate(mode_tags) if tag.kind == 'T' and tag.level == 2)]
            self.assertEqual({w: c.v for w, c in rel.terms.items()}, expected)

    def test_twist_and_framing_powers_have_the_stated_literal_words(self):
        _, relations, sliders, tags = z.build_relations_ring(
            4, z.NumericRing(2, 3, 0), 'bigelow', framing_power=3, theta_power=2)
        rels = {tag.label(4): p for tag, p in zip(tags, relations)}
        # Level two is (x1^3-3)(2*y1-1-x1), with no braid cancellation.
        self.assertEqual({w: c.v for w, c in rels['T_2'].terms.items()}, {
            (1, 1, 1, 4): Fraction(2), (1, 1, 1): Fraction(-1),
            (1, 1, 1, 1): Fraction(-1), (4,): Fraction(-6),
            (): Fraction(3), (1,): Fraction(3),
        })
        theta_square = (1, 2, 2, 1, 1, 2, 2, 1)
        expected = {theta_square + w: c.v for w, c in sliders[(3, 1)].terms.items()}
        expected.update({w: -16*c.v for w, c in sliders[(3, 1)].terms.items()})
        self.assertEqual({w: c.v for w, c in rels['T_3'].terms.items()}, expected)
        # At level four the parameter still selects the canonical reversed-word relation.
        self.assertTrue(all(len(w) <= 9 for w in rels['T_4'].terms))
        _, limited, _, limited_tags = z.build_relations_ring(
            5, z.SymbolicRing(), 'none', max_zipper_level=3, max_untwist_level=2)
        self.assertFalse(any(tag.kind == 'T' for tag in limited_tags))
        self.assertFalse(any(tag.kind == 'Z' and tag.level > 3 for tag in limited_tags))
        self.assertTrue(any(tag.kind == 'FC' for tag in limited_tags))
        self.assertTrue(limited)

    def test_unimplemented_and_invalid_powers_are_explicit_errors(self):
        for kwargs in ({'framing_power': 0}, {'theta_power': 0}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                z.build_relations_ring(3, z.NumericRing(2, 3, 0), 'bigelow', **kwargs)
        with self.assertRaises(NotImplementedError):
            z.build_relations_ring(3, z.SymbolicRing(), 'indeterminate')

    def test_target_and_zero_slider_append_order_and_missing_targets(self):
        names, relations, sliders = z.build_relations(3, 5, 7, 0, 'none')
        self.assertEqual(names, ['x1', 'x2', 'y1', 'y2'])
        before = len(relations)
        self.assertEqual(z.append_zero_sliders(relations, sliders, 3, 3), 1)
        self.assertIs(relations[-1], sliders[(3, 1)])
        self.assertEqual(z.append_target(relations, sliders, 2, 2), before+1)
        self.assertIs(relations[-1], sliders[(2, 2)])
        self.assertEqual(z.append_zero_sliders([], {}, 2, 3), 0)
        with self.assertRaisesRegex(ValueError, 'not in slider dict'):
            z.append_target([], sliders, 4)

    def test_emitter_without_optional_appends_has_a_complete_msolve_body(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)/'input.ms'
            argv = ['zipper_in_ms', '--n', '2', '--q', '2', '--t', '3',
                    '--untwist', 'none', '--out', str(out)]
            stdout = io.StringIO()
            with patch.object(sys, 'argv', argv), redirect_stdout(stdout):
                self.assertEqual(z.main(), 0)
            self.assertEqual(out.read_text(), 'x1,y1\n0\nx1*y1-1,\ny1*x1-1\n')
            self.assertIn('Wrote 2 relations in 2 variables', stdout.getvalue())

    def test_emitter_module_entry_point_supports_target_and_zero_sliders(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)/'input.ms'
            argv = ['zipper_in_ms', '--n', '2', '--out', str(out),
                    '--zero-from', '2', '--target', '2']
            stdout = io.StringIO()
            with patch.object(sys, 'argv', argv), redirect_stdout(stdout):
                with self.assertRaises(SystemExit) as exit_:
                    runpy.run_path(z.__file__, run_name='__main__')
            self.assertEqual(exit_.exception.code, 0)
            lines = out.read_text().splitlines()
            self.assertEqual(lines[-2:], ['-x1+5*y1-4,', '-x1+5*y1-4'])
            self.assertIn('Appended 1 zero-slider relations', stdout.getvalue())
            self.assertIn('Appended S[2][1]', stdout.getvalue())
