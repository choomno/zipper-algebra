"""Certificate serialization, specialization and malformed-input regressions."""
from __future__ import annotations

from contextlib import redirect_stdout
from fractions import Fraction
from io import StringIO
import json
from pathlib import Path
import runpy
import sys
import tempfile
import pytest
from unittest.mock import patch
import warnings

from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.cofactor_io import (
    Cofactor, _strip_outer, num_str, parse_num, parse_word, resolve_target,
    symbolic_target, target_name, validate_qt_polynomial, validate_word,
    var_names, word_from_json, word_str,
)
from zipper_algebra.presentation import Presentation
from zipper_algebra.coefficients import Coef, NumericRing, QTCoef, SymbolicRing
from zipper_algebra.polynomial import Poly


class TestCertificateTextAndWords:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    def test_word_representations_preserve_order_and_unit(self):
        names = var_names(3)
        index = {name: i + 1 for i, name in enumerate(names)}
        assert names == ['x1', 'x2', 'y1', 'y2']
        assert word_str((4, 1, 3), names) == 'y2*x1*y1'
        assert word_str((), names) == '()'
        assert parse_word(' y2 * x1 * y1 ', index) == (4, 1, 3)
        for text in ('', '()'):
            assert parse_word(text, index) == ()
        assert word_from_json([], index) == ()
        assert word_from_json([4, 1], index) == (4, 1)
        assert word_from_json(['y2', 'x1'], index) == (4, 1)
        with pytest.raises(ValueError, match='unknown variable'):
            parse_word('x9', index)
        for word in ('x1', (0,), (True,), (5,)):
            with self.subtests.test(word=word):
                with pytest.raises(ValueError):
                    validate_word(word, 3)
        validate_word([1, 4], 3)

    def test_polynomial_text_round_trip_and_cancellation(self):
        polynomial = {(2, 1): Fraction(-2, 3), (0, 0): Fraction(5, 7)}
        assert parse_num(num_str(polynomial)) == polynomial
        assert num_str({}) == '0'
        assert parse_num('') == {}
        assert parse_num('q^2 - t + 3/2') == {(2, 0): Fraction(1), (0, 1): Fraction(-1), (0, 0): Fraction(3, 2)}
        assert parse_num('q + -q + 0*t') == {}
        assert parse_num(' +q') == {(1, 0): Fraction(1)}
        assert _strip_outer(' (q*t) ') == 'q*t'
        assert _strip_outer('(q)*(t)') == '(q)*(t)'
        assert _strip_outer('q*t') == 'q*t'
        for text in ('x', 'q2', 't2', 'q12', 't12'):
            with self.subtests.test(text=text):
                with pytest.raises(ValueError):
                    parse_num(text)
        for polynomial in ([], {(0,): 1}, {(0, 0): True}, {(0, 0): .5}):
            with self.subtests.test(polynomial=polynomial):
                with pytest.raises(ValueError):
                    validate_qt_polynomial(polynomial)
        validate_qt_polynomial({(0, 1): Fraction(1, 2)})


    def test_incomplete_numerator_expressions_are_rejected(self):
        # Missing factors must not silently change the expression. In particular,
        # a sign with no operand is not a constant coefficient of magnitude one.
        malformed = ('q**t', '*q', 'q*', '-', '+', 'q +  + t', 'q + -')
        for text in malformed:
            with self.subtests.test(text=text):
                with pytest.raises(ValueError):
                    parse_num(text)
        assert parse_num('') == {}


class TestCertificateTargets:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    def test_slider_spellings_and_direct_polynomial(self):
        pres, ring = Presentation(3), NumericRing(2, 3, 0)
        expected = {(3,): Fraction(2), (): Fraction(-1), (1,): Fraction(-1)}
        for target in ((2, 1), [2, 1], ' S[2][1] ', 'x_2'):
            with self.subtests.test(target=target):
                p = resolve_target(pres, target, ring)
                assert {w: c.v for w, c in p.terms.items()} == expected
        p = Poly.scalar(Coef(Fraction(3, 2), 0))
        assert resolve_target(pres, p, ring) is p
        assert symbolic_target(pres, p) == {(): {(0, 0): Fraction(3, 2)}}
        assert target_name(' S[2][1] ') == 'X_2'
        assert target_name('S[2][2]') == 'S[2][2]'
        assert target_name(17) == '17'
        for target in ((2.0, 1), 3, 'S[3][1]suffix', 'missing'):
            with self.subtests.test(target=target):
                with pytest.raises((ValueError, TypeError)):
                    resolve_target(pres, target, ring)

    def test_extra_targets_are_evaluated_over_the_selected_ring(self):
        pres, ring = Presentation(4), NumericRing(2, 3, 0)
        # The named lasso relation is an explicit polynomial operation on X3.
        slider = pres.target(3, 1, ring)
        theta = Poly(0)
        theta.add_term((1, 2, 2, 1), Coef.one(0))
        expected = (theta - Poly.scalar(Coef(4, 0))) * slider
        for name in ('(Theta_3-q^2)X_3', 'EXTRA[(Theta_3-q^2)X_3]'):
            got = resolve_target(pres, name, ring)
            assert {w: c.v for w, c in got.terms.items()} == {w: c.v for w, c in expected.terms.items()}

    def test_symbolic_targets_revalidate_mutable_coefficients(self):
        pres = Presentation(2)
        finite = Poly.scalar(Coef.one(5))
        with pytest.raises(ValueError, match='characteristic zero'):
            symbolic_target(pres, finite)
        p = Poly(0)
        p.add_term((1,), Coef.one(0))
        p.terms[(1,)].char = 5
        with pytest.raises(ValueError, match='characteristic zero'):
            symbolic_target(pres, p)
        p.terms[(1,)] = Coef.one(0)
        p.terms[(1,)].v = .25
        with pytest.raises(ValueError, match='exact polynomials'):
            symbolic_target(pres, p)
        p.terms[(1,)] = QTCoef({(2, 0): 2})
        assert symbolic_target(pres, p) == {(1,): {(2, 0): Fraction(2)}}
        assert symbolic_target(pres, {(): {(0, 0): 0}, (1,): {(1, 0): 2}}) == {(1,): {(1, 0): Fraction(2)}}


class TestCertificateData:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    @pytest.fixture(autouse=True)
    def setup_certificate(self):
        self.pres = Presentation(2)
        self.atoms = AtomSet([Atom('t', {(0, 1): 1})])

    def certificate(self, terms=None):
        if terms is None:
            terms = [((1,), 'INV_1_p', (2,), {(1, 0): Fraction(1)}, (1,))]
        return Cofactor.from_presentation(self.pres, 'X_2', self.atoms, terms,
                                         'serialization-regression', anchor=(2, 3, 0),
                                         notes=['an exact contextual multiple'])

    def test_header_views_and_contextual_specialization(self):
        cof = self.certificate()
        assert cof.header['anchor'] == (2, 3, 0)
        assert cof.header['notes'] == ['an exact contextual multiple']
        assert cof.n_terms == 1
        assert cof.labels() == ['INV_1_p']
        assert cof.relation_indices(self.pres) == [0]
        assert cof.lcm_denominator() == (1,)
        assert '1 terms' in repr(cof)
        actual = cof.evaluate(2, 3, self.pres)
        expected = {(1, 1, 2, 2): Fraction(2, 3), (1, 2): Fraction(-2, 3)}
        assert {w: c.v for w, c in actual.terms.items()} == expected
        target = Poly(0)
        for word, value in expected.items():
            target.add_term(word, Coef(value, 0))
        assert cof.numeric_check(2, 3, self.pres, target) == (True, 0)
        assert cof.numeric_check(2, 3, self.pres, Poly(0)) == (False, 2)
        assert cof.numeric_check(2, 3, self.pres)[0] == False
        with pytest.raises(ZeroDivisionError, match='denominator atom'):
            cof.evaluate(2, 0, self.pres)
        zero = self.certificate([((), 'INV_1_p', (), {}, (0,))])
        assert not zero.evaluate(2, 3, self.pres).terms

    def test_data_boundaries_and_legacy_relation_indices(self):
        for n in (True, 1, 2.0):
            with self.subtests.test(n=n):
                with pytest.raises(ValueError):
                    Cofactor([], n, self.atoms)
        with pytest.raises(ValueError, match='AtomSet'):
            Cofactor([], 2, [])
        with pytest.raises(ValueError, match='hash must be a string'):
            Cofactor([], 2, self.atoms, presentation_hash=3)
        for term in ((), ((), '', (), {}, (0,)), ((), 3, (), {}, (0,))):
            with self.subtests.test(term=term):
                with pytest.raises(ValueError):
                    self.certificate([term])
        pytest.raises(ValueError, self.certificate().validate, Presentation(3))
        wrong = self.certificate()
        wrong.presentation_hash = 'different-content'
        with pytest.raises(ValueError, match='does not match'):
            wrong.validate(self.pres)
        bare = Cofactor([], 2, self.atoms)
        assert 'anchor' not in bare.header
        assert 'notes' not in bare.header
        assert not any((line.startswith('# anchor:') for line in bare.header_lines()))
        cof = self.certificate([((), 'R_0', (), {}, (0,))])
        assert cof.relation_indices(self.pres) == [0]
        cof.terms[0] = ((), 'R_999', (), {}, (0,))
        with pytest.raises(ValueError, match='relation index'):
            cof.relation_indices(self.pres)
        cof.terms[0] = ((), 'not_a_relation', (), {}, (0,))
        with pytest.raises(ValueError, match='relation label'):
            cof.relation_indices(self.pres)
        with pytest.raises(ValueError, match='unknown relation label'):
            cof.validate(self.pres)

    def test_current_file_round_trip_and_explicit_target_override(self):
        cof = self.certificate()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'certificate.txt'
            cof.save(path)
            actual = Cofactor.load(path)
            assert actual.terms == cof.terms
            assert actual.header == cof.header
            override = Cofactor.load(path, pres=self.pres, target='S[2][1]')
            assert override.target == 'S[2][1]'
            assert Cofactor.strip_header(path.read_text()) == '\n'.join(cof.body_lines()) + '\n'
            assert Cofactor.strip_header('# comment\nbody\n') == 'body\n'

    def test_legacy_numeric_file_tracks_anchor_and_notes(self):
        text = '# (q, t) = (-2, 3)\n#\n# note: retained note\n\n+1/2 * x1 * [0] * y1\n-1 * [1]\n'
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'legacy.txt'; path.write_text(text)
            with pytest.raises(ValueError, match='pass pres='):
                Cofactor.load(path)
            cof = Cofactor.load(path, pres=self.pres, atoms=self.atoms, target='X_2')
            assert cof.anchor == (-2, 3, 0)
            assert cof.target == 'X_2'
            assert cof.terms[0] == ((1,), 'INV_1_p', (2,), {(0, 0): Fraction(1, 2)}, (0,))
            assert cof.terms[1][1] == 'INV_1_m'
            assert 'retained note' in cof.notes
            no_anchor = Cofactor.load(path, pres=self.pres)
            assert no_anchor.target == ''
            path.write_text('# provenance without an anchor\n+1 * [0]\n')
            assert Cofactor.load(path, pres=self.pres).anchor is None

    def test_legacy_symbolic_and_header_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'certificate.txt'
            path.write_text('() * INV_1_p * () :: (1/2*q^1*t^0) / (q)\n')
            with pytest.raises(ValueError, match='unknown denominator factor'):
                Cofactor.load(path, pres=self.pres)
            cof = Cofactor.load(path, pres=self.pres, atoms=AtomSet([Atom('q', {(1, 0): 1})]))
            assert cof.terms[0][3] == {(1, 0): Fraction(1, 2)}
            path.write_text('# presentation: wrong-slug\n')
            with pytest.raises(ValueError, match='cannot read n'):
                Cofactor.load(path)
            path.write_text('# presentation: n2-bigelow\n# unknown: provenance\n# anchor: 1,2,0\nnot a cofactor\n')
            with pytest.raises(ValueError, match='cannot parse cofactor line'):
                Cofactor.load(path)
            path.write_text('# presentation: n2-bigelow\n# unknown: provenance\n+0 * [0]\n')
            cof = Cofactor.load(path)
            assert cof.terms == [((), 'R_0', (), {}, ())]
            assert 'unknown: provenance' in cof.notes
            path.write_text('# presentation: n2-bigelow\n# atoms: '+json.dumps(self.atoms.to_json())+'\n')
            cof = Cofactor.load(path, atoms=AtomSet([]))
            assert cof.atoms.names() == ['t']

    def test_malformed_cofactor_line_boundaries(self):
        for text in ('+1 * [0] * [1]', '+1 * x1', '+1 * [-1]', '+1 * [999]'):
            with self.subtests.test(text=text):
                with pytest.raises(ValueError):
                    Cofactor._parse_extract_line(text, {'x1': 1}, self.atoms, self.pres.labels())
        term = Cofactor._parse_extract_line('+1 * x1 * [0] * y1', {'x1': 1, 'y1': 2}, self.atoms, None)
        assert term[:3] == ((1,), 'R_0', (2,))
        with pytest.raises(ValueError, match='expected'):
            Cofactor._parse_symbolic_line('INV_1_p :: 1 / 1', {}, self.atoms)

    def test_rewriter_command_describes_and_writes_without_solver(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'source.txt'; output = Path(directory) / 'copy.txt'
            self.certificate().save(source)
            for flags in ([], ['--n', '2', '--target', 'S[2][1]', '--out', str(output)]):
                stream = StringIO()
                with warnings.catch_warnings(), patch.object(sys, 'argv', ['cofactor_io', str(source), *flags]), redirect_stdout(stream):
                    warnings.simplefilter('ignore', RuntimeWarning)
                    runpy.run_module('zipper_algebra.cofactor.cofactor_io', run_name='__main__')
                assert 'LCM denominator: t' in stream.getvalue()
            assert Cofactor.load(output).target == 'S[2][1]'
