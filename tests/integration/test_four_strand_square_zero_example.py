"""Bounded producer/parser checks; the generic research computation stays local."""
from fractions import Fraction
import importlib.util
import json
import os
import shutil
import subprocess
from pathlib import Path
import sys

import pytest

from zipper_algebra.verification.groebner import GroebnerCertificate, verify_groebner


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT/'examples'/'four_strand_square_zero.py'


def load_example():
    spec = importlib.util.spec_from_file_location('four_strand_square_zero', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mock_environment(module, tmp_path, monkeypatch):
    library = tmp_path/'freegb.lib'
    library.write_text('version="version freegb.lib test";\n')
    metadata = {'singular_version_id': '44105', 'freegb_version': 'version freegb.lib test',
                'freegb_source': {'name': library.name, **module.file_identity(library)}}
    monkeypatch.setattr(module, 'solver_environment', lambda *args, **kwargs: (library, metadata, 0.01))
    return metadata


def test_full_generic_input_retains_all_untwistings_and_parameter_order():
    module = load_example()
    full = module.presentation()
    source, subset = module.header(full, degree=24)
    assert full.algebra.field.parameters == ('q', 't')
    assert full.algebra.generators == module.NAMES
    assert len(full.relations) == 31
    assert full.labels[-1] == 'T_4'
    assert [full.labels[i] for i in subset] == list(full.labels[:-1])
    assert len(subset) == 30
    assert 'freeAlgebra(r,24)' in source
    q, t = (full.algebra.field.parameter(name) for name in ('q', 't'))
    assert module.parse_coefficient('(q^2-t^2)/(q-t)') == q+t
    assert module.parse_coefficient('-2*q + 3/2*t') == -2*q + full.algebra.field.scalar(3)/2*t


@pytest.mark.parametrize('coefficient', [
    '__import__("os")', 'q.real', 'q(t)', 'q[0]', 'q**t', 'q**2', 'q^129', 'True', '1.5',
    'q//2', 'q%2', 'q^-1', 'q**(1+1)', '', 'q\n+t', 'q if t else 1',
])
def test_coefficient_parser_rejects_code_and_unbounded_exponents(coefficient):
    with pytest.raises(ValueError):
        load_example().parse_coefficient(coefficient)


@pytest.mark.parametrize('text', [
    'G|1|0|1,|0|',                         # zero terms are not a canonical sparse record
    'G|1|1|7,|0|',                         # outside the six-generator alphabet
    'G|1|1|1|0|',                          # missing delimiter in word
    'G|1|1|1,|1|',                         # module component is prohibited
    'G|1|1|1,|0|\nG|1|2|1,|0|',          # duplicate polynomial word
    'G|2|1|1,|0|',                         # a truncated stream with missing node 1
    'G|1|1|1,|0|\nD|1|1||1|1|',          # cyclic reference
    'G|1|1|1,|0|\nD|1|1||0|32|',         # undeclared original relation
    'G|1|1|1,|0|\nD|2|1||0|1|',          # derivation without its polynomial
    'G|0|1|1,|0|', 'G|1|1|1,|0|extra',
])
def test_malformed_solver_records_cannot_become_certificates(text):
    with pytest.raises(ValueError):
        load_example().parse_output(text)


def test_runner_exports_replayable_dag_and_provenance_without_claiming_completeness(tmp_path, monkeypatch):
    subprocess.run(['git', 'init', '-q', str(tmp_path)], check=True)
    module = load_example()
    # A tiny true identity is enough to exercise orchestration. It deliberately
    # does not complete the full presentation: the independent verifier must
    # reject it rather than accepting producer success as a dimension proof.
    polynomial = 'G|1|1|1,4,|0|\nG|1|-1||0|\n'
    outputs = {
        'candidate': polynomial+'CANDIDATE_DONE|1\n',
        'derivation': polynomial+'D|1|1||0|1|\nFINISH|1|0|0\n',
        'bridge': 'G|2|1|1,4,|0|\nG|2|-1||0|\nD|2|1||1|1|\nBRIDGE_DONE|1\n',
    }
    stages = []
    def fake_run(binary, source, directory, stage, **limits):
        stages.append(stage)
        assert len(module.presentation().relations) == 31
        assert limits == {'timeout': 180, 'memory_mib': 6144}
        (directory/(stage+'.sing')).write_text(source)
        (directory/(stage+'.log')).write_text(outputs[stage])
        return outputs[stage], 0.01
    monkeypatch.setattr(module, 'run_singular', fake_run)
    runtime = mock_environment(module, tmp_path, monkeypatch)
    output = tmp_path/'run'
    report = module.produce(Path(sys.executable), output)
    assert stages == ['candidate', 'derivation', 'bridge']
    assert report['input_relations'] == 31
    assert report['runtime'] == runtime
    assert report['candidate_polynomials'] == 1
    assert report['proof_nodes'] == 2
    assert report['verification'] == 'candidate produced; independent verification required'
    certificate = GroebnerCertificate.from_json(json.loads((output/'certificate.json').read_text()))
    full = module.presentation()
    assert certificate.presentation_hash == full.content_hash()
    assert certificate.derivations[0][0].kind == 'basis'
    assert certificate.derivations[0][0].index == 1
    assert certificate.derivations[1][0].kind == 'input'
    assert certificate.derivations[1][0].index == 0
    with pytest.raises(ValueError, match='input relation'):
        verify_groebner(full, certificate, max_reductions=1000)
    proof = (output/'proof.txt').read_text()
    assert proof.startswith('META|31|1|2\n')
    assert proof.endswith('END|31|1|2\n')
    assert any(line.startswith('I|31|') for line in proof.splitlines())
    provenance = json.loads((output/'provenance.json').read_text())
    assert {'presentation.json', 'certificate.json', 'proof.txt', 'report.json'} <= set(provenance['outputs'])
    assert [entry['name'] for entry in provenance['inputs']][:2] == [SCRIPT.name, SCRIPT.with_suffix('.sing').name]
    assert (output/'.gitignore').read_text() == '*\n!.gitignore\n'
    ignored = subprocess.run(['git', '-C', str(tmp_path), 'check-ignore',
                              str(output/'certificate.json')],
                             capture_output=True, text=True, check=True)
    assert ignored.stdout.strip() == str(output/'certificate.json')
    with pytest.raises(ValueError, match='empty'):
        module.produce(Path(sys.executable), output)


def test_missing_stage_marker_is_rejected(tmp_path, monkeypatch):
    module = load_example()
    mock_environment(module, tmp_path, monkeypatch)
    monkeypatch.setattr(module, 'run_singular', lambda *args, **kwargs: ('G|1|1|1,|0|\n', 0))
    with pytest.raises(RuntimeError, match='completion marker'):
        module.produce(Path(sys.executable), tmp_path/'run')


@pytest.mark.parametrize('limits', [dict(timeout=0), dict(memory_mib=1), dict(max_nodes=1), dict(degree=11)])
def test_resource_limits_are_validated_before_creating_outputs(tmp_path, limits):
    with pytest.raises(ValueError, match='resource limits'):
        load_example().produce(Path(sys.executable), tmp_path/'run', **limits)
    assert not (tmp_path/'run').exists()


def test_actual_subprocess_wrapper_captures_failures_and_bounds_runtime(tmp_path):
    module = load_example()
    text, seconds = module.run_singular(Path(sys.executable), 'print("bounded check")', tmp_path,
                                        'success', timeout=3, memory_mib=256)
    assert text.strip() == 'bounded check'
    assert 0 <= seconds < 3
    with pytest.raises(RuntimeError, match='status 7'):
        module.run_singular(Path(sys.executable), 'raise SystemExit(7)', tmp_path,
                            'failure', timeout=3, memory_mib=256)
    with pytest.raises(RuntimeError, match='Singular error'):
        module.run_singular(Path(sys.executable), 'print("   ? syntax error")', tmp_path,
                            'text-error', timeout=3, memory_mib=256)
    with pytest.raises(RuntimeError, match='exceeded'):
        module.run_singular(Path(sys.executable), 'import time; time.sleep(10)', tmp_path,
                            'timeout', timeout=0.05, memory_mib=256)


def test_three_term_square_zero_identity_is_a_literal_free_polynomial_equality():
    from zipper_algebra.algebra import CofactorIdentity
    module = load_example()
    full = module.presentation()
    identity = module.square_zero_identity(full)
    assert len(identity.terms) == 3
    assert [full.labels[index] for _, index, _ in identity.terms] == ['T_2', 'INV_1_m', 'INV_1_p']
    assert identity.verify(full)
    restored = CofactorIdentity.from_json(identity.to_json())
    assert restored == identity
    # Perturb the target, preserving all true cofactor terms.
    false = CofactorIdentity(identity.presentation_hash, identity.target+full.algebra.one, identity.terms)
    assert not false.verify(full)


@pytest.mark.parametrize('record', ['G|1|1|2,|0|', 'D|1|1||0|1|'])
def test_bridge_cannot_modify_an_existing_node(record):
    module = load_example()
    proof = 'G|1|1|1,|0|\nD|1|1||0|1|\n'
    with pytest.raises(ValueError, match='modify an existing'):
        module.parse_bridge(proof, record, proved_count=1, inputs=30)
    with pytest.raises(ValueError, match='contiguous'):
        module.parse_output('G|2|1|2,|0|\n')


def test_auxiliary_input_parser_rejects_the_omitted_relation():
    with pytest.raises(ValueError, match='out-of-range'):
        load_example().parse_output('G|1|1|1,|0|\nD|1|1||0|31|\n', inputs=30)


def test_cofactor_failure_precedes_any_solver_execution(tmp_path, monkeypatch):
    module = load_example()
    def fail_identity(full):
        raise ArithmeticError('cofactor test failure')
    def forbidden_solver(*args, **kwargs):
        raise AssertionError('solver must not start')
    monkeypatch.setattr(module, 'square_zero_identity', fail_identity)
    monkeypatch.setattr(module, 'solver_environment', forbidden_solver)
    with pytest.raises(ArithmeticError, match='cofactor test failure'):
        module.produce(Path(sys.executable), tmp_path/'run')


def test_timeout_race_is_reported_and_stdin_is_closed(tmp_path, monkeypatch):
    module = load_example()
    class ExitingProcess:
        pid = 12345
        def __init__(self):
            self.waits = 0
        def wait(self, timeout=None):
            self.waits += 1
            if self.waits == 1:
                raise subprocess.TimeoutExpired('fake Singular', timeout)
            return 0
    def launch(*args, **kwargs):
        assert kwargs['stdin'] is subprocess.DEVNULL
        return ExitingProcess()
    def kill(*args):
        raise ProcessLookupError('already exited')
    monkeypatch.setattr(module.subprocess, 'Popen', launch)
    monkeypatch.setattr(module.os, 'killpg', kill)
    with pytest.raises(RuntimeError, match='exceeded'):
        module.run_singular(Path(sys.executable), 'quit;', tmp_path, 'race', timeout=1, memory_mib=256)


@pytest.mark.parametrize('error', [ArithmeticError('arithmetic'), ZeroDivisionError('division'),
                                 UnicodeError('encoding'), OSError('filesystem')])
def test_cli_reports_expected_io_encoding_and_arithmetic_failures_without_tracebacks(monkeypatch, capsys, error):
    module = load_example()
    monkeypatch.setattr(module.shutil, 'which', lambda command: sys.executable)
    def fail(*args, **kwargs):
        raise error
    monkeypatch.setattr(module, 'produce', fail)
    with pytest.raises(SystemExit) as stopped:
        module.main([])
    assert stopped.value.code == 1
    assert str(error) in capsys.readouterr().err


def test_solver_library_is_resolved_hashed_and_pinned(tmp_path, monkeypatch):
    module = load_example()
    library = tmp_path/'freegb.lib'
    library.write_text('version="version freegb.lib 4.1.2.0 Feb_2019 ";\n')
    def fake_run(binary, source, directory, stage, **limits):
        assert stage == 'environment'
        assert limits['timeout'] <= 15
        text = 'SOLVER_VERSION|44105\nLIBRARY_PATH|'+str(tmp_path)+'\n'
        (directory/'environment.log').write_text(text)
        return text, 0.01
    monkeypatch.setattr(module, 'run_singular', fake_run)
    path, metadata, seconds = module.solver_environment(Path(sys.executable), tmp_path,
                                                       timeout=180, memory_mib=256)
    assert path == library.resolve()
    assert metadata['singular_version_id'] == '44105'
    assert metadata['freegb_version'] == 'version freegb.lib 4.1.2.0 Feb_2019'
    assert metadata['freegb_source']['sha256'] == module.file_identity(library)['sha256']
    source, _ = module.header(module.presentation(), degree=24, library=path,
                              library_version=metadata['freegb_version'])
    assert 'LIB "'+str(path)+'";' in source
    assert 'FREEGB_VERSION|version freegb.lib 4.1.2.0 Feb_2019' in source
    assert seconds == 0.01


def test_genuine_fraction_serialization_round_trips_through_singular(tmp_path):
    binary = os.environ.get('SINGULAR_TEST_BIN') or shutil.which('Singular')
    if binary is None:
        pytest.skip('optional bounded check needs a local Singular executable')
    module = load_example()
    field = module.FIELD
    q, t = (field.parameter(name) for name in field.parameters)
    terms = {(1,): field.scalar(Fraction(1, 2)), (2,): field.scalar(Fraction(-3, 2)),
             (): (q**2-t)/(q-1)}
    source = ('ring r=(0,q,t),(x1,x2,x3,y1,y2,y3),Dp;\nshort=0;\n'
              'poly p='+module.polynomial_source(terms)+';\n'
              'poly expected=number(1)/number(2)*x1-number(3)/number(2)*x2+(q^2-t)/(q-1);\n'
              'print("ROUNDTRIP|"+string(p-expected));\n'
              'print("HALF|"+string(leadcoef(p)));\nquit;\n')
    text, _ = module.run_singular(Path(binary), source, tmp_path, 'fraction', timeout=5, memory_mib=6144)
    assert 'ROUNDTRIP|0' in text.splitlines()
    assert 'HALF|1/2' in text.splitlines()


def test_candidate_records_need_not_arrive_in_index_order(tmp_path, monkeypatch):
    module = load_example()
    mock_environment(module, tmp_path, monkeypatch)
    forward = 'G|1|1|1,4,|0|\nG|1|-1||0|\n'
    reverse = 'G|2|1|4,1,|0|\nG|2|-1||0|\n'
    proof = forward+reverse+'D|1|1||0|1|\nD|2|1||0|2|\nFINISH|2|0|0\n'
    bridge = (forward.replace('G|1|', 'G|3|')+reverse.replace('G|2|', 'G|4|')+
              'D|3|1||1|1|\nD|4|1||1|2|\nBRIDGE_DONE|2\n')
    def run(binary, source, directory, stage, **limits):
        if stage == 'candidate':
            return reverse+forward+'CANDIDATE_DONE|2\n', 0.01
        if stage == 'derivation':
            assert 'Missing=list("63","36");' in source
            return proof, 0.01
        parsed, _ = module.parse_output(forward+reverse)
        expected = ',\n'.join(module.polynomial_source(parsed[index]) for index in (1, 2))
        assert 'ideal Targets='+expected+';' in source
        assert 'B[1]=' in source and 'B[2]=' in source
        return bridge, 0.01
    monkeypatch.setattr(module, 'run_singular', run)
    report = module.produce(Path(sys.executable), tmp_path/'run')
    assert report['candidate_polynomials'] == 2
    assert report['proof_nodes'] == 4
    certificate = GroebnerCertificate.from_json(
        json.loads((tmp_path/'run'/'certificate.json').read_text()))
    assert certificate.basis[:2] == module.presentation().relations[:2]
    assert [(row[0].kind, row[0].index) for row in certificate.derivations] == [
        ('basis', 2), ('basis', 3), ('input', 0), ('input', 1)]
