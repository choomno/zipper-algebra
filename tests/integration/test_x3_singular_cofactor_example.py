"""Bounded protocol and exact-acceptance checks, without a research computation."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

from zipper_algebra.algebra import AlgebraPresentation, CofactorIdentity
from zipper_algebra.cofactor.general import verify_cofactor


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT/'examples'/'x3_singular_cofactor.py'


def load_example():
    spec = importlib.util.spec_from_file_location('x3_singular_cofactor', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def framed(polynomial, *, coefficient='1', marker='TARGET_DONE|1\n'):
    words = lambda word: ''.join(f'{letter},' for letter in word)
    return (''.join(f'G|1|{value}|{words(word)}|0|\n' for word, value in polynomial.terms)
            + f'D|1|{coefficient}||0|1|\n'+marker)


def fake_native(module, tmp_path, monkeypatch, *, coefficient='1', marker='TARGET_DONE|1\n'):
    helper = module.load_helper()
    full, _ = module.presentation()
    # Prove one actual defining relation, with all 31 original inputs retained.
    # This makes frontend/acceptance tests small without priming real discovery.
    target = full.relations[0]
    monkeypatch.setattr(module, 'presentation', lambda: (full, target))
    library = tmp_path/'freegb.lib'
    library.write_text('version="version freegb.lib test";\n')
    metadata = {'singular_version_id': '44105', 'freegb_version': 'version freegb.lib test',
                'freegb_source': {'name': library.name, **module.file_identity(library)}}

    def environment(binary, output, **limits):
        assert limits == {'timeout': 120, 'memory_mib': 4096}
        (output/'environment.sing').write_text('quit;\n')
        (output/'environment.log').write_text('SOLVER_VERSION|44105\n')
        return library, metadata, 0.001

    def run(binary, source, output, stage, **limits):
        assert stage == 'producer'
        assert limits == {'timeout': 120, 'memory_mib': 4096}
        ideal = source.split('ideal I=', 1)[1].split(';\npoly Target=', 1)[0]
        assert ideal == ',\n'.join(helper.polynomial_source(dict(r.terms)) for r in full.relations)
        assert len(full.relations) == 31
        result = framed(target, coefficient=coefficient, marker=marker)
        (output/'producer.sing').write_text(source)
        (output/'producer.log').write_text(result)
        return result, 0.001

    monkeypatch.setattr(helper, 'solver_environment', environment)
    monkeypatch.setattr(helper, 'run_singular', run)
    monkeypatch.setattr(module, 'load_helper', lambda: helper)
    return full, target


def test_literal_original_input_and_fresh_target_search():
    module = load_example()
    full, target = module.presentation()
    assert full.content_hash() == 'c155b61e62039db7d7aa7b710c5bbc38e35495067b472dd08efcc6acabf8b7df'
    assert len(full.relations) == 31 and full.labels[-1] == 'T_4'
    assert target not in full.relations
    assert full.algebra.generators == ('x1', 'x2', 'x3', 'y1', 'y2', 'y3')
    assert full.algebra.field.parameters == ('q', 't')
    source = module.source_for(full, target, Path('/usr/share/singular/LIB/freegb.lib'),
                               degree=24, max_nodes=1000)
    helper = module.load_helper()
    ideal = source.split('ideal I=', 1)[1].split(';\npoly Target=', 1)[0]
    assert ideal == ',\n'.join(helper.polynomial_source(dict(r.terms)) for r in full.relations)
    assert 'poly Target='+helper.polynomial_source(dict(target.terms))+';' in source
    assert 'targetDone=checkTarget();' in source
    assert 'freeAlgebra(r,24)' in source
    assert 'if(res[1]!=0){return(0);}' in source
    assert module.load_helper() is module.load_helper()


def test_publish_reload_and_exactly_verify_original_relation_identity(tmp_path, monkeypatch):
    module = load_example()
    full, target = fake_native(module, tmp_path, monkeypatch)
    output = tmp_path/'run'
    report = module.produce(Path(sys.executable), output)
    reloaded = CofactorIdentity.from_json(json.loads((output/'cofactor.json').read_text()))
    original = AlgebraPresentation.from_json(json.loads((output/'presentation.json').read_text()))
    assert original == full and reloaded.target == target
    assert verify_cofactor(original, reloaded)
    assert report['proof_nodes'] == report['cofactor_terms'] == 1
    assert report['input_relations'] == 31
    provenance = json.loads((output/'provenance.json').read_text())
    assert {'presentation.json', 'target.json', 'cofactor.json', 'graph.json',
            'producer.sing', 'producer.log', 'report.json'} <= set(provenance['outputs'])
    assert {SCRIPT.name, 'four_strand_square_zero.py', 'four_strand_square_zero.sing',
            Path(sys.executable).resolve().name, 'freegb.lib'} <= {
                identity['name'] for identity in provenance['inputs']}
    assert all(identity['sha256'] for identity in provenance['inputs'])
    assert (output/'.gitignore').read_text() == '*\n!.gitignore\n'
    with pytest.raises(ValueError, match='empty'):
        module.produce(Path(sys.executable), output)


def test_native_success_and_target_marker_cannot_substitute_for_exact_proof(tmp_path, monkeypatch):
    module = load_example()
    fake_native(module, tmp_path, monkeypatch, coefficient='2')
    with pytest.raises(ValueError, match='expanded cofactor'):
        module.produce(Path(sys.executable), tmp_path/'run')
    assert not (tmp_path/'run'/'report.json').exists()
    assert not (tmp_path/'run'/'provenance.json').exists()


def test_changed_package_source_cannot_receive_a_verified_report(tmp_path, monkeypatch):
    module = load_example()
    fake_native(module, tmp_path, monkeypatch)

    def changed(producer):
        raise ValueError('package source changed during computation')

    monkeypatch.setattr(module, 'validate_source_snapshot', changed)
    with pytest.raises(ValueError, match='package source changed'):
        module.produce(Path(sys.executable), tmp_path/'run')
    assert not (tmp_path/'run'/'report.json').exists()
    assert not (tmp_path/'run'/'provenance.json').exists()
    assert not (tmp_path/'run'/'cofactor.json').exists()


def test_shared_node_contexts_are_expanded_in_order():
    module = load_example()
    full, _ = module.presentation()
    target = full.algebra.generator('x2')*full.relations[0]*full.algebra.generator('y3')
    first = framed(full.relations[0], marker='').split('D|')[0]
    second = ''.join(f'G|2|{c}|'+''.join(f'{i},' for i in w)+'|0|\n' for w, c in target.terms)
    text = first+'D|1|1||0|1|\n'+second+'D|2|1|2,|1|1|6,\nTARGET_DONE|2\n'
    graph, node = module.parse_target(full, target, text, max_nodes=1000)
    identity = module.expand_derivation(full, graph, node)
    assert verify_cofactor(full, identity)
    left, index, right = identity.terms[0]
    assert index == 0 and left == full.algebra.generator('x2') and right == full.algebra.generator('y3')


@pytest.mark.parametrize('document', ['target.json', 'presentation.json'])
def test_modified_serialized_input_cannot_be_published(tmp_path, monkeypatch, document):
    module = load_example()
    fake_native(module, tmp_path, monkeypatch)
    helper = module.load_helper()
    native = helper.run_singular

    def changed(binary, source, output, stage, **limits):
        text, elapsed = native(binary, source, output, stage, **limits)
        full, _ = module.presentation()
        value = (full.algebra.zero.to_json() if document == 'target.json' else
                 AlgebraPresentation(full.algebra, full.relations[:-1], full.labels[:-1]).to_json())
        (output/document).write_text(json.dumps(value))
        return text, elapsed

    monkeypatch.setattr(helper, 'run_singular', changed)
    with pytest.raises(ValueError, match='published original-input cofactor'):
        module.produce(Path(sys.executable), tmp_path/'run')
    assert not (tmp_path/'run'/'cofactor.json').exists()


@pytest.mark.parametrize('marker', ['', 'TARGET_DONE|0\n', 'TARGET_DONE|2\n',
                                    'TARGET_DONE|1\nTARGET_DONE|1\n', 'TARGET_DONE|bad\n'])
def test_missing_or_wrong_target_marker_is_inconclusive(tmp_path, monkeypatch, marker):
    module = load_example()
    fake_native(module, tmp_path, monkeypatch, marker=marker)
    with pytest.raises((RuntimeError, ValueError), match='target'):
        module.produce(Path(sys.executable), tmp_path/'run')


@pytest.mark.parametrize('text', [
    'G|1|1|1,|0|\nTARGET_DONE|1\n',
    'G|1|1|1,|0|\nD|1|1||0|32|\nTARGET_DONE|1\n',
    'G|1|1|1,|0|\nD|1|1||1|1|\nTARGET_DONE|1\n',
])
def test_malformed_or_unproved_native_nodes_are_rejected(text):
    module = load_example()
    full, target = module.presentation()
    with pytest.raises(ValueError):
        module.parse_target(full, target, text, max_nodes=1000)


def test_wrong_target_and_node_bound_are_rejected():
    module = load_example()
    full, target = module.presentation()
    with pytest.raises(ValueError, match='final target'):
        module.parse_target(full, target, framed(full.relations[0]), max_nodes=1000)
    with pytest.raises(ValueError, match='node bound'):
        module.parse_target(full, target, framed(full.relations[0]), max_nodes=-1)


@pytest.mark.parametrize('limits', [
    {'timeout': 0}, {'timeout': float('nan')}, {'timeout': float('inf')}, {'timeout': True},
    {'memory_mib': 255}, {'degree': 11}, {'degree': 25}, {'max_nodes': 30},
    {'max_operations': 0}, {'max_contexts': 0}, {'max_nodes': True},
])
def test_invalid_limits_do_not_create_outputs(tmp_path, limits):
    with pytest.raises(ValueError, match='resource limits'):
        load_example().produce(Path(sys.executable), tmp_path/'run', **limits)
    assert not (tmp_path/'run').exists()


def test_engine_adaptation_fails_closed_on_source_changes(tmp_path, monkeypatch):
    module = load_example()
    engine = tmp_path/'engine.sing'
    engine.write_text(module.ENGINE.read_text().replace('proc completeToHeads()', 'proc renamed()'))
    monkeypatch.setattr(module, 'ENGINE', engine)
    with pytest.raises(ValueError, match='source shape'):
        module.target_engine()


@pytest.mark.parametrize('path', ['unsafe".lib', 'unsafe\\.lib', 'unsafe\n.lib'])
def test_singular_library_path_cannot_inject_source(path):
    module = load_example()
    full, target = module.presentation()
    with pytest.raises(ValueError, match='library path'):
        module.source_for(full, target, path, degree=24, max_nodes=1000)


def test_cli_help_does_not_launch_research_computation():
    result = subprocess.run([sys.executable, str(SCRIPT), '--help'], capture_output=True,
                            text=True, check=True, timeout=10)
    assert '--singular' in result.stdout and '--max-contexts' in result.stdout


def test_cli_reports_missing_executable(monkeypatch, capsys):
    module = load_example()
    monkeypatch.setattr(module.shutil, 'which', lambda name: None)
    with pytest.raises(SystemExit) as caught:
        module.main(['--singular', 'absent'])
    assert caught.value.code == 2
    assert 'executable not found' in capsys.readouterr().err


@pytest.mark.parametrize('error', [ZeroDivisionError('invalid exact scalar'),
                                 subprocess.SubprocessError('native resource setup failed')])
def test_cli_reports_exact_arithmetic_and_native_setup_failures(monkeypatch, capsys, error):
    module = load_example()
    monkeypatch.setattr(module.shutil, 'which', lambda name: sys.executable)

    def failed(*args, **kwargs):
        raise error

    monkeypatch.setattr(module, 'produce', failed)
    with pytest.raises(SystemExit) as caught:
        module.main([])
    assert caught.value.code == 1
    assert capsys.readouterr().err == f'{error}\n'
