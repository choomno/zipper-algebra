"""Bounded, independent controls for the finite-basis reproduction adapter."""
from dataclasses import replace
from fractions import Fraction
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.finite import WordBasisCertificate, verify_word_basis
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT/'examples'/'four_strand_finite_basis.py'


def load_example():
    spec = importlib.util.spec_from_file_location('four_strand_finite_basis', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def toy_source(module):
    """The exact one-dimensional quotient Q<x1,...,y3>/(x1,...,y3)."""
    helper = module.load_helper()
    field = RationalFunctionField()
    algebra = FreeAlgebra(field, helper.NAMES)
    relations = tuple(algebra.generator(name) for name in algebra.generators)
    full = AlgebraPresentation(algebra, relations, tuple(f'zero_{name}' for name in algebra.generators))
    certificate = GroebnerCertificate(full.content_hash(), relations,
                                     tuple((ProofTerm(field.one, (), 'input', index, ()),)
                                           for index in range(6)))
    # The actual parser is exercised over exact QQ values in this small control;
    # the public four-strand producer always supplies its QQ(q,t) helper.
    helper.parse_coefficient = lambda text: field.scalar(Fraction(text))
    helper.presentation = lambda: full
    return helper, full, certificate


def frames():
    lines = []
    for generator in range(1, 7):
        node = 6+generator
        lines.extend((f'ACTION|{generator}|1|{node}', f'G|{node}|1|{generator},|0|',
                      f'D|{node}|1||1|{generator}|'))
    return '\n'.join(lines)+'\nACTIONS_DONE|6\n'


def write_source(directory, full, certificate):
    directory.mkdir()
    (directory/'presentation.json').write_text(json.dumps(full.to_json()))
    (directory/'certificate.json').write_text(json.dumps(certificate.to_json()))
    records = ['META|6|6|6']
    for kind in ('I', 'G'):
        records.extend(f'{kind}|{index}|1|{index},|0|' for index in range(1, 7))
    records.extend(f'D|{index}|1||0|{index}|' for index in range(1, 7))
    records.append('END|6|6|6')
    (directory/'proof.txt').write_text('\n'.join(records)+'\n')


def test_parser_creates_a_genuinely_verified_finite_basis_from_local_identities():
    module = load_example()
    helper, full, source = toy_source(module)
    certificate = module.parse_actions(frames(), full, source, 6, ((),), helper)
    assert len(certificate.proof.nodes) == 12
    assert certificate.action_nodes == tuple((index,) for index in range(6, 12))
    assert all(action.columns == ((),) for action in certificate.actions)
    restored = WordBasisCertificate.from_json(certificate.to_json())
    assert restored == certificate
    verified = verify_word_basis(full, restored, max_operations=10_000)
    assert verified.dimension == 1
    assert verified.basis_words == ((),)
    for relation in full.relations:
        assert verified.coordinates(relation) == (full.algebra.field.zero,)


@pytest.mark.parametrize('mutate,match', [
    (lambda text: text.replace('ACTION|1|1|7\n', ''), 'missing a column'),
    (lambda text: text.replace('G|7|1|1,|0|', 'G|7|2|1,|0|'), 'target disagrees'),
    (lambda text: text.replace('G|7|1|1,|0|', 'Z|7'), 'target disagrees'),
    (lambda text: text+'A|1|1|1|1,\n', 'outside the basis'),
    (lambda text: text+'A|1|1|0|\n', 'coefficient is zero'),
    (lambda text: text+'A|1|1|1|\nA|1|1|2|\n', 'duplicate sparse action'),
    (lambda text: text+'ACTION|1|1|7\n', 'duplicate or mismatched'),
    (lambda text: text.replace('ACTION|1|1|7', 'ACTION|1|1|8'), 'mismatched'),
    (lambda text: text+'G|1|1|1,|0|\n', 'modify a source'),
    (lambda text: text.replace('D|7|1||1|1|', 'D|7|1||1|7|'), 'outside its declared'),
    (lambda text: text.replace('D|7|1||1|1|', 'D|7|1||0|1|'), 'reference proved candidate'),
    (lambda text: text.replace('ACTIONS_DONE|6', 'ACTIONS_DONE|5'), 'completion marker'),
])
def test_action_parser_rejects_missing_columns_and_false_or_malformed_targets(mutate, match):
    module = load_example()
    helper, full, source = toy_source(module)
    with pytest.raises(ValueError, match=match):
        module.parse_actions(mutate(frames()), full, source, 6, ((),), helper)


def test_false_local_witness_is_rejected_by_the_finite_verifier():
    module = load_example()
    helper, full, source = toy_source(module)
    text = frames().replace('D|7|1||1|1|', 'D|7|2||1|1|')
    candidate = module.parse_actions(text, full, source, 6, ((),), helper)
    with pytest.raises(ValueError, match='identity'):
        verify_word_basis(full, candidate, max_operations=10_000)


def test_zero_target_is_retained_as_an_explicit_graph_node():
    module = load_example()
    helper, full, source = toy_source(module)
    # The proposed singleton basis {x1} supplies an x1*x1 action target, which
    # disagrees with an exported zero polynomial; reject it before verification.
    text = frames().replace('G|7|1|1,|0|', 'Z|7').replace('D|7|1||1|1|\n', '')
    word = (1,)
    with pytest.raises(ValueError, match='target disagrees'):
        module.parse_actions(text, full, source, 6, (word,), helper)
    # A proper independent control uses basis {1,x1} for Q[x1]/(x1²), embedded
    # in the six-generator alphabet with all other generators zero.
    algebra = full.algebra
    x = algebra.generator('x1')
    relations = (x*x, *full.relations[1:])
    full = replace(full, relations=relations)
    source = GroebnerCertificate(full.content_hash(), relations,
                                tuple((ProofTerm(algebra.field.one, (), 'input', index, ()),)
                                      for index in range(6)))
    data = []
    for g in range(1, 7):
        for j in range(1, 3):
            k = 6+(g-1)*2+j
            data.append(f'ACTION|{g}|{j}|{k}')
            if (g, j) == (1, 1):
                data.extend((f'Z|{k}', 'A|1|1|1|1,'))
            else:
                word = (g,) if j == 1 else (g, 1)
                data.append(f'G|{k}|1|'+''.join(f'{letter},' for letter in word)+'|0|')
                data.append(f'D|{k}|1||1|{g}|'+('1,' if j == 2 and g != 1 else ''))
    data.append('ACTIONS_DONE|12')
    certificate = module.parse_actions('\n'.join(data), full, source, 6, ((), (1,)), helper)
    assert certificate.proof.nodes[6].is_zero()
    assert certificate.proof.derivations[6] == ()
    assert verify_word_basis(full, certificate, max_operations=20_000).dimension == 2


def test_source_encodings_must_agree_exactly(tmp_path):
    module = load_example()
    helper, full, source = toy_source(module)
    directory = tmp_path/'source'
    write_source(directory, full, source)
    assert module.load_source(directory, helper) == (full, source, 6)
    proof = directory/'proof.txt'
    original = proof.read_text()
    proof.write_text(original.replace('G|1|1|1,|0|', 'G|1|2|1,|0|'))
    with pytest.raises(ValueError, match='polynomials disagree'):
        module.load_source(directory, helper)
    proof.write_text(original.replace('D|1|1||0|1|', 'D|1|2||0|1|'))
    with pytest.raises(ValueError, match='derivations disagree'):
        module.load_source(directory, helper)
    proof.write_text(original.replace('END|6|6|6', 'END|6|5|6'))
    with pytest.raises(ValueError, match='counts or terminal'):
        module.load_source(directory, helper)


def test_duplicate_json_keys_and_a_different_presentation_are_rejected(tmp_path):
    module = load_example()
    path = tmp_path/'duplicate.json'
    path.write_text('{"schema":1,"schema":2}')
    with pytest.raises(ValueError, match='duplicate JSON key'):
        module.read_json(path)
    helper, full, source = toy_source(module)
    directory = tmp_path/'source'
    write_source(directory, full, source)
    helper.presentation = lambda: replace(full, labels=tuple(label+'changed' for label in full.labels))
    with pytest.raises(ValueError, match='not the full generic'):
        module.load_source(directory, helper)


def test_orchestration_saves_one_replayable_certificate_and_binds_all_sources(tmp_path, monkeypatch):
    module = load_example()
    subprocess.run(['git', 'init', str(tmp_path)], capture_output=True, text=True, check=True)
    helper, full, source = toy_source(module)
    source_run = tmp_path/'source'
    write_source(source_run, full, source)
    monkeypatch.setattr(module, 'load_helper', lambda: helper)
    library = tmp_path/'freegb.lib'
    library.write_text('version="test freegb";')
    runtime = {'singular_version_id': '44105', 'freegb_version': 'test freegb',
               'freegb_source': {'name': 'freegb.lib', **module.file_identity(library)}}
    monkeypatch.setattr(helper, 'solver_environment', lambda *args, **kwargs: (library, runtime, 0.01))
    def fake_run(binary, text, directory, stage, **limits):
        assert stage == 'actions'
        assert 'list Heads;' in text
        assert 'Heads[i]' in text
        assert limits == {'timeout': 180, 'memory_mib': 6144}
        (directory/'actions.sing').write_text(text)
        (directory/'actions.log').write_text(frames())
        return frames(), 0.01
    monkeypatch.setattr(helper, 'run_singular', fake_run)
    output = tmp_path/'custom-output'
    report = module.produce(Path(sys.executable), source_run, output)
    assert report['action_columns'] == 6
    assert report['proof_nodes'] == 12
    assert report['matrix_nonzero_entries'] == 0
    assert report['verification'] == 'candidate produced; finite-basis verification required'
    certificate = WordBasisCertificate.from_json(module.read_json(output/'word-basis-certificate.json'))
    assert verify_word_basis(full, certificate, max_operations=10_000).dimension == 1
    provenance = json.loads((output/'provenance.json').read_text())
    assert set(provenance['outputs']) == {'actions.sing', 'actions.log', 'actions.txt', 'presentation.json',
                                         'word-basis-certificate.json', 'report.json', '.gitignore'}
    inputs = {entry['name'] for entry in provenance['inputs']}
    assert {'presentation.json', 'certificate.json', 'proof.txt', SCRIPT.name, SCRIPT.with_suffix('.sing').name,
            'four_strand_square_zero.py', 'freegb.lib'} <= inputs
    assert (output/'.gitignore').read_text() == '*\n!.gitignore\n'
    ignored = subprocess.run(['git', '-C', str(tmp_path), 'check-ignore',
                              str(output/'word-basis-certificate.json')],
                             capture_output=True, text=True, check=True)
    assert ignored.stdout.strip() == str(output/'word-basis-certificate.json')
    with pytest.raises(ValueError, match='empty'):
        module.produce(Path(sys.executable), source_run, output)


def test_infinite_candidate_word_language_is_bounded():
    module = load_example()
    _, _, source = toy_source(module)
    with pytest.raises(ValueError, match='enumeration budget'):
        module.standard_words(source.basis[:1], maximum=10)


@pytest.mark.parametrize('token', ['NaN', 'Infinity', '-Infinity'])
def test_nonfinite_json_is_not_valid_certificate_input(tmp_path, token):
    module = load_example()
    path = tmp_path/'nonfinite.json'
    path.write_text('{"value":'+token+'}')
    with pytest.raises(ValueError, match='nonfinite'):
        module.read_json(path)


@pytest.mark.parametrize('limits', [dict(timeout=0), dict(memory_mib=1), dict(degree=11)])
def test_invalid_resource_limits_do_not_create_a_run(tmp_path, limits):
    with pytest.raises(ValueError, match='resource limits'):
        load_example().produce(Path(sys.executable), tmp_path/'source', tmp_path/'output', **limits)
    assert not (tmp_path/'output').exists()


def production_case(tmp_path, monkeypatch, *, provenance=False, during_solver=None):
    """One exact finite quotient, with mutable temporary computation inputs."""
    module = load_example()
    helper, full, source = toy_source(module)
    source_run = tmp_path/'source'
    write_source(source_run, full, source)
    if provenance:
        module.write_bundle_provenance(source_run, module.SOURCE_NAMES)
    engine = tmp_path/'engine.sing'
    engine.write_text(module.ENGINE.read_text())
    helper_path = tmp_path/'helper.py'
    helper_path.write_text('# bounded test helper identity\n')
    binary = tmp_path/'solver'
    binary.write_text('bounded test solver identity\n')
    library = tmp_path/'freegb.lib'
    library.write_text('version="test freegb";')
    runtime = {'singular_version_id': '44105', 'freegb_version': 'test freegb',
               'freegb_source': {'name': 'freegb.lib', **module.file_identity(library)}}
    monkeypatch.setattr(module, 'ENGINE', engine)
    monkeypatch.setattr(module, 'HELPER', helper_path)
    monkeypatch.setattr(module, 'load_helper', lambda: helper)
    monkeypatch.setattr(helper, 'solver_environment', lambda *args, **kwargs: (library, runtime, 0.01))
    paths = {'engine': engine, 'helper': helper_path, 'binary': binary, 'library': library,
             'presentation': source_run/'presentation.json', 'certificate': source_run/'certificate.json',
             'proof': source_run/'proof.txt', 'provenance': source_run/'provenance.json'}
    def fake_run(binary, text, directory, stage, **limits):
        (directory/'actions.sing').write_text(text)
        (directory/'actions.log').write_text(frames())
        if during_solver is not None:
            during_solver(module, paths)
        return frames(), 0.01
    monkeypatch.setattr(helper, 'run_singular', fake_run)
    return module, binary, source_run, tmp_path/'output', paths


@pytest.mark.parametrize('changed', ['engine', 'helper', 'binary', 'library', 'presentation',
                                     'certificate', 'proof', 'provenance'])
def test_changed_computation_inputs_cannot_receive_output_provenance(tmp_path, monkeypatch, changed):
    def mutate(module, paths):
        with paths[changed].open('a') as stream:
            stream.write('\n')
    module, binary, source_run, output, paths = production_case(
        tmp_path, monkeypatch, provenance=True, during_solver=mutate)
    with pytest.raises(ValueError, match='input changed during'):
        module.produce(binary, source_run, output)
    assert not (output/'provenance.json').exists()
    assert not (output/'word-basis-certificate.json').exists()


def test_package_source_snapshot_is_revalidated_before_outputs_are_published(tmp_path, monkeypatch):
    from zipper_algebra.execution import artifacts
    def mutate_package_snapshot(module, paths):
        monkeypatch.setattr(artifacts, 'current_source_files', lambda names: {})
    module, binary, source_run, output, paths = production_case(
        tmp_path, monkeypatch, during_solver=mutate_package_snapshot)
    with pytest.raises(ValueError, match='source files changed'):
        module.produce(binary, source_run, output)
    assert not (output/'provenance.json').exists()
    assert not (output/'word-basis-certificate.json').exists()


def test_removed_input_and_new_source_provenance_are_detected(tmp_path, monkeypatch):
    first = tmp_path/'removed'
    first.mkdir()
    module, binary, source_run, output, paths = production_case(
        first, monkeypatch, during_solver=lambda module, paths: paths['engine'].unlink())
    with pytest.raises(ValueError, match='input became unreadable'):
        module.produce(binary, source_run, output)
    second = tmp_path/'appeared'
    second.mkdir()
    module, binary, source_run, output, paths = production_case(
        second, monkeypatch, during_solver=lambda module, paths: paths['provenance'].write_text('{}'))
    with pytest.raises(ValueError, match='provenance appeared'):
        module.produce(binary, source_run, output)


@pytest.mark.parametrize('supplied,status', [(True, 'present_verified'), (False, 'absent')])
def test_source_provenance_status_and_actual_json_wire_round_trip(tmp_path, monkeypatch, supplied, status):
    module, binary, source_run, output, paths = production_case(tmp_path, monkeypatch, provenance=supplied)
    parsed_texts = []
    strict = module.strict_loads
    def observe_wire(text):
        parsed_texts.append(text)
        return strict(text)
    monkeypatch.setattr(module, 'strict_loads', observe_wire)
    report = module.produce(binary, source_run, output)
    assert report['source_provenance_status'] == status
    saved = (output/'word-basis-certificate.json').read_text()
    assert saved in parsed_texts
    candidate = WordBasisCertificate.from_json(module.read_json(output/'word-basis-certificate.json'))
    assert candidate.basis_words == ((),)
    source_inputs = {entry['name'] for entry in module.read_json(output/'provenance.json')['inputs']}
    assert ('provenance.json' in source_inputs) is supplied


def test_source_provenance_must_match_the_actual_consumed_files(tmp_path, monkeypatch):
    module, binary, source_run, output, paths = production_case(tmp_path, monkeypatch, provenance=True)
    with paths['certificate'].open('a') as stream:
        stream.write('\n')
    with pytest.raises(ValueError, match='source provenance does not match consumed file: certificate.json'):
        module.produce(binary, source_run, output)
    assert not output.exists()


def test_read_limit_applies_to_actual_bytes_instead_of_file_stat(tmp_path, monkeypatch):
    module = load_example()
    monkeypatch.setattr(module, 'INPUT_BUDGET', 16)
    path = tmp_path/'grew.json'
    path.write_bytes(b'"'+b'x'*20+b'"')
    actual_stat = Path.stat
    def earlier_smaller_stat(current, *args, **kwargs):
        result = actual_stat(current, *args, **kwargs)
        if current == path:
            fields = list(result)
            fields[6] = 0
            return os.stat_result(fields)
        return result
    monkeypatch.setattr(Path, 'stat', earlier_smaller_stat)
    assert path.stat().st_size == 0
    with pytest.raises(ValueError, match='input budget'):
        module.read_json(path)
    path.write_bytes(b'{"value":1}')
    assert module.read_json(path) == {'value': 1}


def test_trace_normal_form_disagreement_has_its_own_diagnostic():
    module = load_example()
    helper, full, source = toy_source(module)
    with pytest.raises(ValueError, match=r'traced normal form.*TRACE_MISMATCH\|1\|1'):
        module.parse_actions('TRACE_MISMATCH|1|1\n', full, source, 6, ((),), helper)


def test_provenance_link_can_verify_named_finite_certificate_outputs(tmp_path):
    module = load_example()
    names = ('presentation.json', 'word-basis-certificate.json')
    for name in names:
        (tmp_path/name).write_text('{"bounded":"test output"}')
    module.write_bundle_provenance(tmp_path, names)
    snapshots = [(tmp_path/name, module.file_identity(tmp_path/name), True) for name in names]
    assert module.source_provenance_status(tmp_path, snapshots, names=names) == 'present_verified'
    record = module.read_json(tmp_path/'provenance.json')
    record['outputs']['word-basis-certificate.json']['sha256'] = '0'*64
    (tmp_path/'provenance.json').write_text(json.dumps(record))
    with pytest.raises(ValueError, match='does not match consumed file: word-basis-certificate.json'):
        module.source_provenance_status(tmp_path, snapshots, names=names)
