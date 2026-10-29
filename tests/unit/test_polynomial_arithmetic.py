"""Algebraic laws, presentation variants and exact solver-input emission."""
from contextlib import redirect_stdout
from fractions import Fraction
import io
import operator
from pathlib import Path
import subprocess
import sys
import tempfile
import pytest
from unittest.mock import patch

from zipper_algebra import msolve
from zipper_algebra.coefficients import Coef, QTCoef, NumericRing, SymbolicRing
from zipper_algebra.polynomial import Poly
from zipper_algebra.relations import (RelTag, build_relations, build_relations_ring,
                                     append_zero_sliders, append_target)
from zipper_algebra.cofactor import ncpoly as nc


class TestScalarAndFreePolynomialArithmetic:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    def test_exact_scalar_signs_equality_and_rendering(self):
        half = Coef(Fraction(1, 2), 0)
        third = Coef(Fraction(1, 3), 0)
        assert (half - third).v == Fraction(1, 6)
        assert (-half).v == Fraction(-1, 2)
        assert str(half) == '1/2'
        assert str(Coef(2, 0)) == '2'
        assert str(Coef(-1, 5)) == '4'
        assert Coef.zero(5).is_zero()
        assert not half.is_zero()
        assert half != object()
        for operation in (operator.add, operator.sub, operator.mul):
            with self.subtests.test(operation=operation):
                with pytest.raises(TypeError):
                    operation(half, object())

    def test_symbolic_arithmetic_with_rational_constants_and_reflection(self):
        q = QTCoef({(1, 0): 1})
        t = QTCoef({(0, 1): 1})
        assert ((q + t) * (q - t)).v == {(2, 0): Fraction(1), (0, 2): Fraction(-1)}
        assert (2 - q).v == {(0, 0): Fraction(2), (1, 0): Fraction(-1)}
        assert (Coef(2, 0) - q).v == {(0, 0): Fraction(2), (1, 0): Fraction(-1)}
        assert (Coef(2, 0) + q).v == {(0, 0): Fraction(2), (1, 0): Fraction(1)}
        assert (Fraction(1, 2) * q).v == {(1, 0): Fraction(1, 2)}
        assert (q + Fraction(1, 3)).v == {(1, 0): Fraction(1), (0, 0): Fraction(1, 3)}
        assert (q - q).is_zero()
        assert (q + t * t).bidegree() == (1, 2)
        assert QTCoef(0).bidegree() == (0, 0)
        assert q != object()
        assert q != Coef(1, 5)
        for operation in (operator.add, operator.sub, operator.mul):
            for left, right in ((q, object()), (object(), q)):
                with self.subtests.test(operation=operation):
                    with pytest.raises(TypeError):
                        operation(left, right)
        with pytest.raises(TypeError):
            q + Coef(1, 7)

    def test_symbolic_rendering_distinguishes_powers_constants_and_rationals(self):
        assert str(QTCoef(0)) == '0'
        polynomial = QTCoef({(0, 0): -2, (1, 0): 1, (0, 1): 3,
                               (2, 3): -4})
        assert str(polynomial) == '(-2)+(3)*t+(1)*q+(-4)*q^2*t^3'
        with pytest.raises(ValueError, match='non-integer coefficient'):
            str(QTCoef(Fraction(1, 2)))
        assert NumericRing(2, 3, 0).describe() == '(q,t)=(2,3) over Q'
        assert NumericRing(2, 3, 5).scalar(Fraction(1, 2)).v == 3
        symbolic = SymbolicRing()
        assert symbolic.scalar(2).v == {(0, 0): Fraction(2)}
        assert symbolic.describe() == 'q,t symbolic over Q'

    def test_free_products_preserve_order_and_cancel_coefficient_collisions(self):
        x, y = Poly.generator(1, 0), Poly.generator(2, 0)
        square = (x + y) * (x + y)
        mixed = square - x*x - y*y
        assert {w: c.v for w, c in mixed.terms.items()} == {(1, 2): Fraction(1), (2, 1): Fraction(1)}
        assert {w: c.v for w, c in (-mixed).terms.items()} == {(1, 2): Fraction(-1), (2, 1): Fraction(-1)}
        assert not (mixed - mixed).terms
        assert {w: c.v for w, c in mixed.scalar_mul(Coef(3, 0)).terms.items()} == {(1, 2): Fraction(3), (2, 1): Fraction(3)}
        with pytest.raises(TypeError):
            mixed.scalar_mul(object())
        with pytest.raises(TypeError):
            mixed.add_term((1,), Coef(1, 5))
        for operation in (operator.add, operator.sub, operator.mul):
            with self.subtests.test(operation=operation):
                with pytest.raises(TypeError):
                    operation(x, object())
        p = Poly(0)
        p.add_term((1,), Coef(2, 0))
        p.add_term((1,), Coef(3, 0))
        p.add_term((), Coef(0, 0))
        assert p.terms[1,].v == 5
        assert str(p.terms[1,]) == '5'
        assert Poly.scalar(Coef(0, 0)).terms == {}

    def test_msolve_rendering_is_unambiguous_about_signs_and_empty_words(self):
        p = Poly(0)
        for word, coefficient in [((1, 2), 1), ((2, 1), -1), ((1,), 2), ((), -3)]:
            p.add_term(word, Coef(coefficient, 0))
        assert p.to_msolve_string(['x', 'y']) == 'x*y-y*x+2*x-3'
        assert Poly(0).to_msolve_string(['x']) == '0'
        assert Poly.scalar(Coef(2, 0)).to_msolve_string(['x']) == '2'


class TestSparsePolynomialLaws:
    def test_commuting_parameter_arithmetic_and_text(self):
        q, t = nc.qt_q(), nc.qt_t()
        assert nc.qt_zero() == {}
        assert nc.qt_monomial(2, 3, 0) == {}
        assert nc.qt_monomial(2, 3, Fraction(1, 2)) == {(2, 3): Fraction(1, 2)}
        assert nc.qt_mul(nc.qt_add(q, t), nc.qt_sub(q, t)) == {(2, 0): Fraction(1), (0, 2): Fraction(-1)}
        assert nc.qt_sub(q, q) == {}
        assert nc.qt_add(q, nc.qt_neg(q)) == {}
        assert nc.qt_sub(q, t) == {(1, 0): Fraction(1), (0, 1): Fraction(-1)}
        assert nc.qt_pow(q, 0) == {(0, 0): Fraction(1)}
        assert nc.qt_pow(q, 2) == {(2, 0): Fraction(1)}
        assert nc.qt_eq(q, {(1, 0): Fraction(1)})
        assert not nc.qt_eq(q, t)
        assert nc.qt_is_zero(nc.qt_const(0))
        assert not nc.qt_is_zero(q)
        assert nc.qt_str({}) == '0'
        assert nc.qt_str({(0, 0): Fraction(-1, 2), (0, 1): Fraction(-1), (1, 0): Fraction(1), (2, 0): Fraction(3, 2)}) == '-1/2 - t^1 + q^1 + 3/2*q^2'

    def test_noncommuting_polynomials_and_zero_coefficient_normalization(self):
        one = {(0, 0): Fraction(1)}
        minus = {(0, 0): Fraction(-1)}
        x, y = nc.nc_var(1), nc.nc_var_named('y', {'y': 2})
        assert nc.nc_zero() == {}
        assert nc.nc_const({}) == {}
        assert nc.nc_const(one) == {(): one}
        with pytest.raises(KeyError):
            nc.nc_var_named('z', {'x': 1})
        assert nc.nc_sub(nc.nc_add(x, y), x) == y
        assert nc.nc_sub(x, y) == {(1,): one, (2,): minus}
        assert nc.nc_add(x, nc.nc_neg(x)) == {}
        assert nc.nc_sub(x, x) == {}
        assert nc.nc_add(x, x) == {(1,): {(0, 0): Fraction(2)}}
        assert nc.nc_sub(nc.nc_scalar(x, nc.qt_const(2)), x) == x
        assert not nc.nc_is_zero(x)
        assert nc.nc_is_zero(nc.nc_zero())
        assert nc.nc_size(nc.nc_add(x, y)) == 2
        assert nc.nc_mul(nc.nc_add(x, y), nc.nc_sub(x, y)) == {(1, 1): one, (1, 2): minus, (2, 1): one, (2, 2): minus}
        # Cross terms of (x+x^2)(x-x^2) cancel, despite arising from different factors.
        left = {(1,): one, (1, 1): one}
        right = {(1,): one, (1, 1): minus}
        assert nc.nc_mul(left, right) == {(1, 1): one, (1, 1, 1, 1): minus}
        assert nc.nc_scalar(x, {}) == {}
        assert nc.nc_scalar(x, nc.qt_const(2)) == {(1,): {(0, 0): Fraction(2)}}
        # Explicit zero coefficients represent the zero polynomial as well.
        assert nc.nc_mul({(1,): {}}, y) == {}
        assert nc.nc_scalar({(1,): {}}, one) == {}


class TestPresentationBuildersAndEmitter:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    def test_all_relation_labels_encode_anchor_and_variant(self):
        cases = [
            (RelTag('INV', 2, None, 'm'), 5, 'INV_2_m'),
            (RelTag('BR', 1, 2, 'y'), 5, 'BR_1_y'),
            (RelTag('BR', 1, 2, 'x'), 5, 'BR_1'),
            (RelTag('FC', 1, 3, 'xy'), 4, 'FC_xy'),
            (RelTag('FC', 1, 4, 'xy'), 5, 'FC_14_xy'),
            (RelTag('MIX', 1, 2, 'a'), 5, 'MIX_1_a'),
            (RelTag('Z', 3, 1), 5, 'Z_3'),
            (RelTag('T', 3, 2), 5, 'T_3_j2'),
            (RelTag('EXTRA', 3, 2), 5, 'EXTRA_3_2'),
            (RelTag('EXTRA', 3, 2, 'power'), 5, 'EXTRA_3_2_power'),
        ]
        for tag, n, expected in cases:
            with self.subtests.test(expected=expected):
                assert tag.label(n) == expected
        assert repr(RelTag('T', 3, 1)) == "RelTag(T,3,1,'')"

    def test_slider_and_two_strand_relations_against_explicit_coefficients(self):
        _, _, sliders, _ = build_relations_ring(3, NumericRing(5, 7, 0), 'bigelow')
        assert {w: c.v for w, c in sliders[3, 1].terms.items()} == {(3, 4, 3): Fraction(125), (3, 4): Fraction(-100), (3, 4, 1): Fraction(-25), (1, 2, 3): Fraction(-5), (1, 2): Fraction(4), (1, 2, 1): Fraction(1)}
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
            _, relations, _, mode_tags = build_relations_ring(3, NumericRing(5, 7, 0), mode)
            rel = relations[next(i for i, tag in enumerate(mode_tags) if tag.kind == 'T' and tag.level == 2)]
            assert {w: c.v for w, c in rel.terms.items()} == expected

    def test_twist_and_framing_powers_have_the_stated_literal_words(self):
        _, relations, sliders, tags = build_relations_ring(
            4, NumericRing(2, 3, 0), 'bigelow', framing_power=3, theta_power=2)
        rels = {tag.label(4): p for tag, p in zip(tags, relations)}
        # Level two is (x1^3-3)(2*y1-1-x1), with no braid cancellation.
        assert {w: c.v for w, c in rels['T_2'].terms.items()} == {(1, 1, 1, 4): Fraction(2), (1, 1, 1): Fraction(-1), (1, 1, 1, 1): Fraction(-1), (4,): Fraction(-6), (): Fraction(3), (1,): Fraction(3)}
        theta_square = (1, 2, 2, 1, 1, 2, 2, 1)
        expected = {theta_square + w: c.v for w, c in sliders[(3, 1)].terms.items()}
        expected.update({w: -16*c.v for w, c in sliders[(3, 1)].terms.items()})
        assert {w: c.v for w, c in rels['T_3'].terms.items()} == expected
        # At level four the parameter still selects the canonical reversed-word relation.
        assert all((len(w) <= 9 for w in rels['T_4'].terms))
        _, limited, _, limited_tags = build_relations_ring(
            5, SymbolicRing(), 'none', max_zipper_level=3, max_untwist_level=2)
        assert not any((tag.kind == 'T' for tag in limited_tags))
        assert not any((tag.kind == 'Z' and tag.level > 3 for tag in limited_tags))
        assert any((tag.kind == 'FC' for tag in limited_tags))
        assert limited

    def test_unimplemented_and_invalid_powers_are_explicit_errors(self):
        for kwargs in ({'framing_power': 0}, {'theta_power': 0}):
            with self.subtests.test(kwargs=kwargs):
                with pytest.raises(ValueError):
                    build_relations_ring(3, NumericRing(2, 3, 0), 'bigelow', **kwargs)
        with pytest.raises(NotImplementedError):
            build_relations_ring(3, SymbolicRing(), 'indeterminate')

    def test_target_and_zero_slider_append_order_and_missing_targets(self):
        names, relations, sliders = build_relations(3, 5, 7, 0, 'none')
        assert names == ['x1', 'x2', 'y1', 'y2']
        before = len(relations)
        assert append_zero_sliders(relations, sliders, 3, 3) == 1
        assert relations[-1] is sliders[3, 1]
        assert append_target(relations, sliders, 2, 2) == before + 1
        assert relations[-1] is sliders[2, 2]
        assert append_zero_sliders([], {}, 2, 3) == 0
        with pytest.raises(ValueError, match='not in slider dict'):
            append_target([], sliders, 4)

    def test_emitter_without_optional_appends_has_a_complete_msolve_body(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)/'input.ms'
            argv = ['msolve', '--n', '2', '--q', '2', '--t', '3',
                    '--untwist', 'none', '--out', str(out)]
            stdout = io.StringIO()
            with patch.object(sys, 'argv', argv), redirect_stdout(stdout):
                assert msolve.main() == 0
            assert out.read_text() == 'x1,y1\n0\nx1*y1-1,\ny1*x1-1\n'
            assert 'Wrote 2 relations in 2 variables' in stdout.getvalue()

    def test_emitter_module_entry_point_supports_explicit_augmentation_and_zero_sliders(self, tmp_path):
        out = tmp_path / "input.ms"
        result = subprocess.run(
            [sys.executable, "-m", "zipper_algebra.msolve", "--n", "2", "--out", str(out),
             "--zero-from", "2", "--append-target", "2"],
            text=True, capture_output=True, check=True)
        assert result.stderr == ""
        lines = out.read_text().splitlines()
        assert lines[-2:] == ["-x1+5*y1-4,", "-x1+5*y1-4"]
        assert "Appended 1 zero-slider relations" in result.stdout
        assert "Augmented presentation: added S[2][1] as defining relation" in result.stdout
