"""Standalone GAP identities preserve exact coefficients and ordered words."""
from contextlib import redirect_stdout, redirect_stderr
from fractions import Fraction
from copy import deepcopy
import io
import hashlib
import json
import shutil
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from zipper_algebra import Cofactor, Presentation
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.verification.gap import gap_script, main
from zipper_algebra.provenance import validate_provenance

ROOT = Path(__file__).resolve().parents[1]


class TestGAPCertificateExport(unittest.TestCase):
    def setUp(self):
        self.presentation = Presentation(2)
        self.atoms = AtomSet([Atom('display_name_is_not_executable', {(1, 0): 1})])
        self.cofactor = Cofactor.from_presentation(
            self.presentation, 'X_2', self.atoms,
            [((1,), 'INV_1_p', (2,), {(0, 0): Fraction(1, 2)}, (3,))], 'test')

    def test_exact_fraction_denominator_and_noncommuting_contexts(self):
        result = gap_script(self.cofactor, self.presentation)
        self.assertIn('FunctionField(Rationals, ["q", "t"])', result)
        self.assertIn('FreeAssociativeAlgebraWithOne(K, 2)', result)
        self.assertIn('Length(g) <> 3 or g[1] <> One(A)', result)
        self.assertIn('((1/2)*q^0*t^0)', result)
        self.assertIn('((1)*q^1*t^0)^3', result)
        self.assertIn('*g[2]*r[1]*g[3]', result)
        self.assertNotIn('display_name_is_not_executable', result)
        self.assertIn('residual = Zero(A)', result)
        self.assertIn('QUIT_GAP(1)', result)

    def test_zero_polynomial_empty_certificate_and_explicit_zero_target(self):
        zero = Cofactor.from_presentation(self.presentation, '', AtomSet([]), [], 'test')
        result = gap_script(zero, self.presentation, target={})
        self.assertIn('residual := Zero(A)', result)
        self.assertNotIn('residual := residual -', result)
        with self.assertRaisesRegex(ValueError, 'explicit target'):
            gap_script(zero, self.presentation)

    def test_constant_coefficients_and_zero_numerator_are_rendered(self):
        self.cofactor.terms = [((), 'INV_1_p', (), {(0, 0): Fraction(-2, 3)}, (0,)),
                               ((), 'INV_1_p', (), {}, (0,))]
        result = gap_script(self.cofactor, self.presentation)
        self.assertIn('((-2/3)*q^0*t^0)/(1)', result)
        self.assertIn('((0)/(1))*One(A)*r[1]*One(A)', result)

    def test_mutated_invalid_relation_and_denominator_are_rejected(self):
        self.cofactor.terms[0] = ((1,), 'not-a-relation', (), {(0, 0): 1}, (0,))
        with self.assertRaisesRegex(ValueError, 'relation label'):
            gap_script(self.cofactor, self.presentation)
        self.cofactor.terms[0] = ((1,), 'INV_1_p', (), {(0, 0): 1}, (-1,))
        with self.assertRaisesRegex(ValueError, 'nonnegative'):
            gap_script(self.cofactor, self.presentation)

    def test_invalid_word_indices_denominator_lengths_and_scalar_types_are_rejected(self):
        left, label, right, numerator, denominator = self.cofactor.terms[0]
        cases = [((index,), label, right, numerator, denominator)
                 for index in (0, 3, '1', True)]
        cases += [(left, label, right, numerator, exponents) for exponents in ((), (0, 0), (True,))]
        cases += [(left, label, right, {(0, 0): value}, denominator)
                  for value in ('1); QUIT_GAP(0)', 0.5, True)]
        for term in cases:
            with self.subTest(term=term), self.assertRaises(ValueError):
                self.cofactor.terms = [term]
                gap_script(self.cofactor, self.presentation)

    def test_claim_labels_and_reversed_word_are_auditable(self):
        left, label, _, numerator, denominator = self.cofactor.terms[0]
        self.cofactor.terms = [(left, label, (2, 1), numerator, denominator)]
        script = gap_script(self.cofactor, self.presentation, target='X_2')
        line = next(line for line in script.splitlines() if line.startswith('# claim: '))
        claim = json.loads(line[len('# claim: '):])
        self.assertEqual(claim, {'n': 2, 'untwist': 'bigelow', 'target': 'X_2',
                               'presentation_hash': self.presentation.content_hash(), 'terms': 1,
                               'generators': ['x1', 'y1']})
        self.assertIn('*g[3]*g[2];;', script)
        self.assertIn('# "INV_1_p"', script)

    def test_producer_comment_is_json_encoded(self):
        result = gap_script(self.cofactor, self.presentation,
                            producer={'value': 'untrusted\nQUIT_GAP(0);'})
        line = next(line for line in result.splitlines() if line.startswith('# producer: '))
        self.assertEqual(json.loads(line[len('# producer: '):]),
                         {'value': 'untrusted\nQUIT_GAP(0);'})
        self.assertNotIn('\nQUIT_GAP(0);', result)

    def test_current_certificate_cli_emits_source_and_input_provenance(self):
        with redirect_stdout(io.StringIO()) as output:
            self.assertEqual(main([str(ROOT / 'examples/x3/certificate.txt'), '--n', '4',
                                   '--target', 'X_3']), 0)
        script = output.getvalue()
        self.assertEqual(script.count('residual := residual -'), 129)
        producer_line = next(line for line in script.splitlines() if line.startswith('# producer: '))
        producer = json.loads(producer_line[len('# producer: '):])
        self.assertEqual(validate_provenance(producer), producer)
        input_line = next(line for line in script.splitlines() if line.startswith('# input-certificate: '))
        identity = json.loads(input_line[len('# input-certificate: '):])
        raw = (ROOT / 'examples/x3/certificate.txt').read_bytes()
        self.assertEqual(identity, {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()})
        self.assertIn('src/zipper_algebra/verification/gap.py', producer['source_files'])

    def test_module_command_emits_an_executable_check_without_running_gap(self):
        result = subprocess.run(
            [sys.executable, '-m', 'zipper_algebra.verification.gap',
             str(ROOT / 'examples/x3/certificate.txt'), '--n', '4'],
            text=True, capture_output=True, check=True)
        self.assertEqual(result.stderr, '')
        self.assertIn('CheckCertificate := function()', result.stdout)
        self.assertEqual(result.stdout.count('residual := residual -'), 129)

    def test_cli_reports_missing_invalid_and_undecodable_input(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.txt'
            for payload in (None, b'not a cofactor', b'\xff',
                            b'# atoms: {}\n# target: X_2\n() * INV_1_p * () :: (1) / (1)\n'):
                if payload is not None:
                    path.write_bytes(payload)
                with self.subTest(payload=payload), redirect_stderr(io.StringIO()), \
                     self.assertRaises(SystemExit) as error:
                    main([str(path), '--n', '2'])
                self.assertEqual(error.exception.code, 2)


@unittest.skipUnless(shutil.which('gap'), 'GAP is an optional independent checker')
class TestNativeGAPCertificateCheck(unittest.TestCase):
    def run_script(self, script, flags=('-T', '--quitonbreak')):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'check.g'
            path.write_text(script)
            return subprocess.run([shutil.which('gap'), '-A', '-q', '-r', *flags, str(path)],
                                  stdin=subprocess.DEVNULL, text=True, capture_output=True, timeout=30)

    def test_generic_certificate_rejects_changes_to_coefficients_words_denominators_and_target(self):
        presentation = Presentation(4)
        original = Cofactor.load(str(ROOT / 'examples/x3/certificate.txt'), pres=presentation)
        for change in ('none', 'numerator', 'denominator', 'word_order', 'target'):
            cofactor = deepcopy(original)
            target = 'X_2' if change == 'target' else 'X_3'
            if change in ('numerator', 'denominator'):
                left, label, right, numerator, denominator = cofactor.terms[0]
                numerator, denominator = dict(numerator), list(denominator)
                if change == 'numerator':
                    numerator[next(iter(numerator))] += 1
                else:
                    denominator[0] += 1
                cofactor.terms[0] = (left, label, right, numerator, tuple(denominator))
            if change == 'word_order':
                index = next(i for i, term in enumerate(cofactor.terms)
                             if len(term[2]) > 1 and term[2][::-1] != term[2])
                left, label, right, numerator, denominator = cofactor.terms[index]
                cofactor.terms[index] = (left, label, right[::-1], numerator, denominator)
            result = self.run_script(gap_script(cofactor, presentation, target))
            with self.subTest(change=change):
                self.assertEqual(result.returncode, int(change != 'none'), result.stdout + result.stderr)
                self.assertEqual(result.stderr, '')
                self.assertEqual(result.stdout, 'PASS: exact cofactor identity over Q(q,t)\n'
                                 if change == 'none' else 'FAIL: nonzero free-algebra residual\n')

    def test_nondefault_presentation_and_explicit_rational_target(self):
        presentation = Presentation(2, untwist='none')
        atoms = AtomSet([Atom('q', {(1, 0): 1})])
        cofactor = Cofactor.from_presentation(
            presentation, '', atoms,
            [((1,), 'INV_1_p', (2,), {(3, 0): Fraction(1, 2)}, (3,))], 'test')
        # (q^3/(2*q^3))*x1*(x1*y1-1)*y1 is this literal free polynomial.
        target = {(1, 1, 2, 2): {(0, 0): Fraction(1, 2)},
                  (1, 2): {(0, 0): Fraction(-1, 2)}}
        script = gap_script(cofactor, presentation, target)
        claim = json.loads(next(line[len('# claim: '):] for line in script.splitlines()
                                if line.startswith('# claim: ')))
        self.assertEqual(claim['untwist'], 'none')
        self.assertEqual(claim['target'], str(target))
        result = self.run_script(script)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(result.stdout, 'PASS: exact cofactor identity over Q(q,t)\n')
        self.assertEqual(result.stderr, '')

    def test_runtime_error_cannot_be_mistaken_for_success_without_command_line_error_flags(self):
        presentation = Presentation(4)
        cofactor = Cofactor.load(str(ROOT / 'examples/x3/certificate.txt'), pres=presentation)
        script = gap_script(cofactor, presentation).replace('Length(g) <> 7', 'Length(g) <> 99')
        for flags in ((), ('-T',)):
            result = self.run_script(script, flags)
            with self.subTest(flags=flags):
                self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
                self.assertNotIn('PASS:', result.stdout)
                self.assertIn('Unexpected free-algebra generator ordering', result.stderr)


if __name__ == '__main__':
    unittest.main()
