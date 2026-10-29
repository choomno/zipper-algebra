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
import pytest

from zipper_algebra import Cofactor, Presentation
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.verification.gap import gap_script, main
from zipper_algebra.provenance import validate_provenance

ROOT = Path(__file__).resolve().parents[1]


class TestGAPCertificateExport:
    @pytest.fixture(autouse=True)
    def resources(self, request):
        self.presentation = Presentation(2)
        self.atoms = AtomSet([Atom('display_name_is_not_executable', {(1, 0): 1})])
        self.cofactor = Cofactor.from_presentation(
            self.presentation, 'X_2', self.atoms,
            [((1,), 'INV_1_p', (2,), {(0, 0): Fraction(1, 2)}, (3,))], 'test')

    def test_exact_fraction_denominator_and_noncommuting_contexts(self):
        result = gap_script(self.cofactor, self.presentation)
        assert 'FunctionField(Rationals, ["q", "t"])' in result
        assert 'FreeAssociativeAlgebraWithOne(K, 2)' in result
        assert 'Length(g) <> 3 or g[1] <> One(A)' in result
        assert '((1/2)*q^0*t^0)' in result
        assert '((1)*q^1*t^0)^3' in result
        assert '*g[2]*r[1]*g[3]' in result
        assert 'display_name_is_not_executable' not in result
        assert 'residual = Zero(A)' in result
        assert 'QUIT_GAP(1)' in result

    def test_zero_polynomial_empty_certificate_and_explicit_zero_target(self):
        zero = Cofactor.from_presentation(self.presentation, '', AtomSet([]), [], 'test')
        result = gap_script(zero, self.presentation, target={})
        assert 'residual := Zero(A)' in result
        assert 'residual := residual -' not in result
        with pytest.raises(ValueError, match='explicit target'):
            gap_script(zero, self.presentation)

    def test_constant_coefficients_and_zero_numerator_are_rendered(self):
        self.cofactor.terms = [((), 'INV_1_p', (), {(0, 0): Fraction(-2, 3)}, (0,)),
                               ((), 'INV_1_p', (), {}, (0,))]
        result = gap_script(self.cofactor, self.presentation)
        assert '((-2/3)*q^0*t^0)/(1)' in result
        assert '((0)/(1))*One(A)*r[1]*One(A)' in result

    def test_mutated_invalid_relation_and_denominator_are_rejected(self):
        self.cofactor.terms[0] = ((1,), 'not-a-relation', (), {(0, 0): 1}, (0,))
        with pytest.raises(ValueError, match='relation label'):
            gap_script(self.cofactor, self.presentation)
        self.cofactor.terms[0] = ((1,), 'INV_1_p', (), {(0, 0): 1}, (-1,))
        with pytest.raises(ValueError, match='nonnegative'):
            gap_script(self.cofactor, self.presentation)

    def test_invalid_word_indices_denominator_lengths_and_scalar_types_are_rejected(self, subtests):
        left, label, right, numerator, denominator = self.cofactor.terms[0]
        cases = [((index,), label, right, numerator, denominator)
                 for index in (0, 3, '1', True)]
        cases += [(left, label, right, numerator, exponents) for exponents in ((), (0, 0), (True,))]
        cases += [(left, label, right, {(0, 0): value}, denominator)
                  for value in ('1); QUIT_GAP(0)', 0.5, True)]
        for term in cases:
            with subtests.test(term=term), pytest.raises(ValueError):
                self.cofactor.terms = [term]
                gap_script(self.cofactor, self.presentation)

    def test_claim_labels_and_reversed_word_are_auditable(self):
        left, label, _, numerator, denominator = self.cofactor.terms[0]
        self.cofactor.terms = [(left, label, (2, 1), numerator, denominator)]
        script = gap_script(self.cofactor, self.presentation, target='X_2')
        line = next(line for line in script.splitlines() if line.startswith('# claim: '))
        claim = json.loads(line[len('# claim: '):])
        assert claim == {'n': 2, 'untwist': 'bigelow', 'target': 'X_2', 'presentation_hash': self.presentation.content_hash(), 'terms': 1, 'generators': ['x1', 'y1']}
        assert '*g[3]*g[2];;' in script
        assert '# "INV_1_p"' in script

    def test_producer_comment_is_json_encoded(self):
        result = gap_script(self.cofactor, self.presentation,
                            producer={'value': 'untrusted\nQUIT_GAP(0);'})
        line = next(line for line in result.splitlines() if line.startswith('# producer: '))
        assert json.loads(line[len('# producer: '):]) == {'value': 'untrusted\nQUIT_GAP(0);'}
        assert '\nQUIT_GAP(0);' not in result

    def test_current_certificate_cli_emits_source_and_input_provenance(self):
        with redirect_stdout(io.StringIO()) as output:
            assert main([str(ROOT / 'examples/x3/certificate.txt'), '--n', '4', '--target', 'X_3']) == 0
        script = output.getvalue()
        assert script.count('residual := residual -') == 129
        producer_line = next(line for line in script.splitlines() if line.startswith('# producer: '))
        producer = json.loads(producer_line[len('# producer: '):])
        assert validate_provenance(producer) == producer
        input_line = next(line for line in script.splitlines() if line.startswith('# input-certificate: '))
        identity = json.loads(input_line[len('# input-certificate: '):])
        raw = (ROOT / 'examples/x3/certificate.txt').read_bytes()
        assert identity == {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}
        assert 'src/zipper_algebra/verification/gap.py' in producer['source_files']

    def test_module_command_emits_an_executable_check_without_running_gap(self):
        result = subprocess.run(
            [sys.executable, '-m', 'zipper_algebra.verification.gap',
             str(ROOT / 'examples/x3/certificate.txt'), '--n', '4'],
            text=True, capture_output=True, check=True)
        assert result.stderr == ''
        assert 'CheckCertificate := function()' in result.stdout
        assert result.stdout.count('residual := residual -') == 129

    @pytest.mark.parametrize("payload", [
        None,
        b'not a cofactor',
        b'\xff',
        b'# atoms: {}\n# target: X_2\n() * INV_1_p * () :: (1) / (1)\n',
    ], ids=["missing", "malformed", "invalid-utf8", "malformed-atoms"])
    def test_cli_reports_missing_invalid_and_undecodable_input(self, tmp_path, payload):
        path = tmp_path / 'input.txt'
        if payload is not None:
            path.write_bytes(payload)
        with redirect_stderr(io.StringIO()), pytest.raises(SystemExit) as error:
            main([str(path), '--n', '2'])
        assert error.value.code == 2


@pytest.mark.skipif(not shutil.which('gap'), reason='GAP is an optional independent checker')
class TestNativeGAPCertificateCheck:
    def run_script(self, script, flags=('-T', '--quitonbreak')):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'check.g'
            path.write_text(script)
            return subprocess.run([shutil.which('gap'), '-A', '-q', '-r', *flags, str(path)],
                                  stdin=subprocess.DEVNULL, text=True, capture_output=True, timeout=30)

    @pytest.mark.parametrize("change", ["none", "numerator", "denominator", "word_order", "target"])
    def test_generic_certificate_rejects_changes_to_coefficients_words_denominators_and_target(self, change):
        presentation = Presentation(4)
        original = Cofactor.load(str(ROOT / 'examples/x3/certificate.txt'), pres=presentation)
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
        assert result.returncode == int(change != 'none'), result.stdout + result.stderr
        assert result.stderr == ''
        assert result.stdout == ('PASS: exact cofactor identity over Q(q,t)\n' if change == 'none' else 'FAIL: nonzero free-algebra residual\n')

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
        assert claim['untwist'] == 'none'
        assert claim['target'] == str(target)
        result = self.run_script(script)
        assert result.returncode == 0, result.stdout + result.stderr
        assert result.stdout == 'PASS: exact cofactor identity over Q(q,t)\n'
        assert result.stderr == ''

    @pytest.mark.parametrize("flags", [(), ("-T",)], ids=["default-error-handling", "no-break-loop"])
    def test_runtime_error_cannot_be_mistaken_for_success_without_command_line_error_flags(self, flags):
        presentation = Presentation(4)
        cofactor = Cofactor.load(str(ROOT / 'examples/x3/certificate.txt'), pres=presentation)
        script = gap_script(cofactor, presentation).replace('Length(g) <> 7', 'Length(g) <> 99')
        result = self.run_script(script, flags)
        assert result.returncode == 1, result.stdout + result.stderr
        assert 'PASS:' not in result.stdout
        assert 'Unexpected free-algebra generator ordering' in result.stderr
