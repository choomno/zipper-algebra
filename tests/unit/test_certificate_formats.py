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
import unittest
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


class TestCertificateTextAndWords(unittest.TestCase):
    def test_word_representations_preserve_order_and_unit(self):
        names = var_names(3)
        index = {name: i + 1 for i, name in enumerate(names)}
        self.assertEqual(names, ['x1', 'x2', 'y1', 'y2'])
        self.assertEqual(word_str((4, 1, 3), names), 'y2*x1*y1')
        self.assertEqual(word_str((), names), '()')
        self.assertEqual(parse_word(' y2 * x1 * y1 ', index), (4, 1, 3))
        for text in ('', '()'):
            self.assertEqual(parse_word(text, index), ())
        self.assertEqual(word_from_json([], index), ())
        self.assertEqual(word_from_json([4, 1], index), (4, 1))
        self.assertEqual(word_from_json(['y2', 'x1'], index), (4, 1))
        with self.assertRaisesRegex(ValueError, 'unknown variable'):
            parse_word('x9', index)
        for word in ('x1', (0,), (True,), (5,)):
            with self.subTest(word=word), self.assertRaises(ValueError):
                validate_word(word, 3)
        validate_word([1, 4], 3)

    def test_polynomial_text_round_trip_and_cancellation(self):
        polynomial = {(2, 1): Fraction(-2, 3), (0, 0): Fraction(5, 7)}
        self.assertEqual(parse_num(num_str(polynomial)), polynomial)
        self.assertEqual(num_str({}), '0')
        self.assertEqual(parse_num(''), {})
        self.assertEqual(parse_num('q^2 - t + 3/2'),
                         {(2, 0): Fraction(1), (0, 1): Fraction(-1), (0, 0): Fraction(3, 2)})
        self.assertEqual(parse_num('q + -q + 0*t'), {})
        self.assertEqual(parse_num(' +q'), {(1, 0): Fraction(1)})
        self.assertEqual(_strip_outer(' (q*t) '), 'q*t')
        self.assertEqual(_strip_outer('(q)*(t)'), '(q)*(t)')
        self.assertEqual(_strip_outer('q*t'), 'q*t')
        for text in ('x', 'q2', 't2', 'q12', 't12'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_num(text)
        for polynomial in ([], {(0,): 1}, {(0, 0): True}, {(0, 0): .5}):
            with self.subTest(polynomial=polynomial), self.assertRaises(ValueError):
                validate_qt_polynomial(polynomial)
        validate_qt_polynomial({(0, 1): Fraction(1, 2)})


    def test_incomplete_numerator_expressions_are_rejected(self):
        # Missing factors must not silently change the expression. In particular,
        # a sign with no operand is not a constant coefficient of magnitude one.
        malformed = ('q**t', '*q', 'q*', '-', '+', 'q +  + t', 'q + -')
        for text in malformed:
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_num(text)
        self.assertEqual(parse_num(''), {})


class TestCertificateTargets(unittest.TestCase):
    def test_slider_spellings_and_direct_polynomial(self):
        pres, ring = Presentation(3), NumericRing(2, 3, 0)
        expected = {(3,): Fraction(2), (): Fraction(-1), (1,): Fraction(-1)}
        for target in ((2, 1), [2, 1], ' S[2][1] ', 'x_2'):
            with self.subTest(target=target):
                p = resolve_target(pres, target, ring)
                self.assertEqual({w: c.v for w, c in p.terms.items()}, expected)
        p = Poly.scalar(Coef(Fraction(3, 2), 0))
        self.assertIs(resolve_target(pres, p, ring), p)
        self.assertEqual(symbolic_target(pres, p), {(): {(0, 0): Fraction(3, 2)}})
        self.assertEqual(target_name(' S[2][1] '), 'X_2')
        self.assertEqual(target_name('S[2][2]'), 'S[2][2]')
        self.assertEqual(target_name(17), '17')
        for target in ((2.0, 1), 3, 'S[3][1]suffix', 'missing'):
            with self.subTest(target=target), self.assertRaises((ValueError, TypeError)):
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
            self.assertEqual({w: c.v for w, c in got.terms.items()},
                             {w: c.v for w, c in expected.terms.items()})

    def test_symbolic_targets_revalidate_mutable_coefficients(self):
        pres = Presentation(2)
        finite = Poly.scalar(Coef.one(5))
        with self.assertRaisesRegex(ValueError, 'characteristic zero'):
            symbolic_target(pres, finite)
        p = Poly(0)
        p.add_term((1,), Coef.one(0))
        p.terms[(1,)].char = 5
        with self.assertRaisesRegex(ValueError, 'characteristic zero'):
            symbolic_target(pres, p)
        p.terms[(1,)] = Coef.one(0)
        p.terms[(1,)].v = .25
        with self.assertRaisesRegex(ValueError, 'exact polynomials'):
            symbolic_target(pres, p)
        p.terms[(1,)] = QTCoef({(2, 0): 2})
        self.assertEqual(symbolic_target(pres, p), {(1,): {(2, 0): Fraction(2)}})
        self.assertEqual(symbolic_target(pres, {(): {(0, 0): 0}, (1,): {(1, 0): 2}}),
                         {(1,): {(1, 0): Fraction(2)}})


class TestCertificateData(unittest.TestCase):
    def setUp(self):
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
        self.assertEqual(cof.header['anchor'], (2, 3, 0))
        self.assertEqual(cof.header['notes'], ['an exact contextual multiple'])
        self.assertEqual(cof.n_terms, 1)
        self.assertEqual(cof.labels(), ['INV_1_p'])
        self.assertEqual(cof.relation_indices(self.pres), [0])
        self.assertEqual(cof.lcm_denominator(), (1,))
        self.assertIn('1 terms', repr(cof))
        actual = cof.evaluate(2, 3, self.pres)
        expected = {(1, 1, 2, 2): Fraction(2, 3), (1, 2): Fraction(-2, 3)}
        self.assertEqual({w: c.v for w, c in actual.terms.items()}, expected)
        target = Poly(0)
        for word, value in expected.items():
            target.add_term(word, Coef(value, 0))
        self.assertEqual(cof.numeric_check(2, 3, self.pres, target), (True, 0))
        self.assertEqual(cof.numeric_check(2, 3, self.pres, Poly(0)), (False, 2))
        self.assertEqual(cof.numeric_check(2, 3, self.pres)[0], False)
        with self.assertRaisesRegex(ZeroDivisionError, 'denominator atom'):
            cof.evaluate(2, 0, self.pres)
        zero = self.certificate([((), 'INV_1_p', (), {}, (0,))])
        self.assertFalse(zero.evaluate(2, 3, self.pres).terms)

    def test_data_boundaries_and_legacy_relation_indices(self):
        for n in (True, 1, 2.0):
            with self.subTest(n=n), self.assertRaises(ValueError):
                Cofactor([], n, self.atoms)
        with self.assertRaisesRegex(ValueError, 'AtomSet'):
            Cofactor([], 2, [])
        with self.assertRaisesRegex(ValueError, 'hash must be a string'):
            Cofactor([], 2, self.atoms, presentation_hash=3)
        for term in ((), ((), '', (), {}, (0,)), ((), 3, (), {}, (0,))):
            with self.subTest(term=term), self.assertRaises(ValueError):
                self.certificate([term])
        self.assertRaises(ValueError, self.certificate().validate, Presentation(3))
        wrong = self.certificate()
        wrong.presentation_hash = 'different-content'
        with self.assertRaisesRegex(ValueError, 'does not match'):
            wrong.validate(self.pres)
        bare = Cofactor([], 2, self.atoms)
        self.assertNotIn('anchor', bare.header)
        self.assertNotIn('notes', bare.header)
        self.assertFalse(any(line.startswith('# anchor:') for line in bare.header_lines()))
        cof = self.certificate([((), 'R_0', (), {}, (0,))])
        self.assertEqual(cof.relation_indices(self.pres), [0])
        cof.terms[0] = ((), 'R_999', (), {}, (0,))
        with self.assertRaisesRegex(ValueError, 'relation index'):
            cof.relation_indices(self.pres)
        cof.terms[0] = ((), 'not_a_relation', (), {}, (0,))
        with self.assertRaisesRegex(ValueError, 'relation label'):
            cof.relation_indices(self.pres)
        with self.assertRaisesRegex(ValueError, 'unknown relation label'):
            cof.validate(self.pres)

    def test_current_file_round_trip_and_explicit_target_override(self):
        cof = self.certificate()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'certificate.txt'
            cof.save(path)
            actual = Cofactor.load(path)
            self.assertEqual(actual.terms, cof.terms)
            self.assertEqual(actual.header, cof.header)
            override = Cofactor.load(path, pres=self.pres, target='S[2][1]')
            self.assertEqual(override.target, 'S[2][1]')
            self.assertEqual(Cofactor.strip_header(path.read_text()), '\n'.join(cof.body_lines())+'\n')
            self.assertEqual(Cofactor.strip_header('# comment\nbody\n'), 'body\n')

    def test_legacy_numeric_file_tracks_anchor_and_notes(self):
        text = '# (q, t) = (-2, 3)\n#\n# note: retained note\n\n+1/2 * x1 * [0] * y1\n-1 * [1]\n'
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'legacy.txt'; path.write_text(text)
            with self.assertRaisesRegex(ValueError, 'pass pres='):
                Cofactor.load(path)
            cof = Cofactor.load(path, pres=self.pres, atoms=self.atoms, target='X_2')
            self.assertEqual(cof.anchor, (-2, 3, 0))
            self.assertEqual(cof.target, 'X_2')
            self.assertEqual(cof.terms[0], ((1,), 'INV_1_p', (2,), {(0, 0): Fraction(1, 2)}, (0,)))
            self.assertEqual(cof.terms[1][1], 'INV_1_m')
            self.assertIn('retained note', cof.notes)
            no_anchor = Cofactor.load(path, pres=self.pres)
            self.assertEqual(no_anchor.target, '')
            path.write_text('# provenance without an anchor\n+1 * [0]\n')
            self.assertIsNone(Cofactor.load(path, pres=self.pres).anchor)

    def test_legacy_symbolic_and_header_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'certificate.txt'
            path.write_text('() * INV_1_p * () :: (1/2*q^1*t^0) / (q)\n')
            cof = Cofactor.load(path, pres=self.pres)
            self.assertEqual(cof.terms[0][3], {(1, 0): Fraction(1, 2)})
            path.write_text('# presentation: wrong-slug\n')
            with self.assertRaisesRegex(ValueError, 'cannot read n'):
                Cofactor.load(path)
            path.write_text('# presentation: n2-bigelow\n# unknown: provenance\n# anchor: 1,2,0\nnot a cofactor\n')
            with self.assertRaisesRegex(ValueError, 'cannot parse cofactor line'):
                Cofactor.load(path)
            path.write_text('# presentation: n2-bigelow\n# unknown: provenance\n+0 * [0]\n')
            cof = Cofactor.load(path)
            self.assertEqual(cof.terms, [((), 'R_0', (), {}, AtomSet.DEFAULT.zero())])
            self.assertIn('unknown: provenance', cof.notes)
            path.write_text('# presentation: n2-bigelow\n# atoms: '+json.dumps(self.atoms.to_json())+'\n')
            cof = Cofactor.load(path, atoms=AtomSet([]))
            self.assertEqual(cof.atoms.names(), ['t'])

    def test_malformed_cofactor_line_boundaries(self):
        for text in ('+1 * [0] * [1]', '+1 * x1', '+1 * [-1]', '+1 * [999]'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                Cofactor._parse_extract_line(text, {'x1': 1}, self.atoms, self.pres.labels())
        term = Cofactor._parse_extract_line('+1 * x1 * [0] * y1', {'x1': 1, 'y1': 2}, self.atoms, None)
        self.assertEqual(term[:3], ((1,), 'R_0', (2,)))
        with self.assertRaisesRegex(ValueError, 'expected'):
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
                self.assertIn('LCM denominator: t', stream.getvalue())
            self.assertEqual(Cofactor.load(output).target, 'S[2][1]')
