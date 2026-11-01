"""Publication boundaries around independently specified finite-algebra proofs."""
import json
from pathlib import Path
import shutil
import sys

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.execution import production as pipeline
from zipper_algebra.execution.artifacts import ArtifactError, canonical_json, digest_bytes, file_identity
from zipper_algebra.execution.config import RunConfig
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.production import GroebnerProduction
from zipper_algebra.verification.finite import WordBasisCertificate, verify_word_basis
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm, verify_groebner


# This process emits only the literal input relations and their direct proofs.
# The chosen three relations already are a complete basis for M2(Q).
PRODUCER = '''
import argparse
from pathlib import Path
p=argparse.ArgumentParser()
for flag in ("-i", "-o", "-p"): p.add_argument(flag, required=True)
a,_=p.parse_known_args()
relations=Path(a.i).read_text().splitlines()[2:]
Path(a.o).write_text("\\n".join(r.rstrip(",").replace("-", " -").replace("+", " +") for r in relations)+"\\n")
Path(a.p).write_text("\\n".join("[i%d]"%i for i in range(len(relations)))+"\\n")
'''


@pytest.fixture
def setup(tmp_path, monkeypatch):
    field = RationalFunctionField(())
    algebra = FreeAlgebra(field, ('u', 'v'))
    u, v = (algebra.generator(name) for name in algebra.generators)
    presentation = AlgebraPresentation(algebra, (u*u, v*v, v*u+u*v-1), ('uu', 'vv', 'vu'))
    source = tmp_path/'literal_producer.py'
    source.write_text(PRODUCER)
    config = RunConfig(Path(sys.executable), tmp_path/'out', command_prefix=(source,), flags=())
    certificate = GroebnerCertificate(presentation.content_hash(), presentation.relations,
        tuple((ProofTerm(field.one, (), 'input', i, ()),) for i in range(3)))
    # Reconstruction itself is tested separately. This boundary supplies an
    # honest complete generic candidate whose source identities are literal.
    result = GroebnerProduction(certificate, verify_groebner(presentation, certificate),
                               {'systems': 3})
    monkeypatch.setattr(pipeline, 'reconstruction_provenance', lambda: {'engine': 'test-boundary'})
    monkeypatch.setattr(pipeline, 'reconstruct_groebner', lambda *a, **k: result)
    return presentation, config, result


def read(path):
    return json.loads(Path(path).read_text())


def assert_inconclusive(report, stage):
    assert report['status'] == 'inconclusive'
    assert not report['verified'] and report['dimension'] is None
    assert not report['basis_complete'] and report['stage'] == stage
    assert 'generation' not in report and 'artifacts' not in report


def test_production_publishes_only_byte_bound_certified_algebra(setup):
    presentation, config, _ = setup
    report = pipeline.produce_algebra(presentation, {}, config, refresh=True, shift=())
    assert report['status'] == 'verified' and report['verified']
    assert report['basis_complete'] and report['dimension'] == 4
    assert read(config.output_root/'result.json') == report
    assert (config.output_root/'samples.sqlite3').is_file()
    assert (config.output_root/'.gitignore').read_text() == '*\n!.gitignore\n'
    artifacts = report['artifacts']
    assert set(artifacts) == {'presentation.json', 'groebner-certificate.json',
                             'word-basis-certificate.json', 'production.json', 'provenance.json'}
    restored_presentation = AlgebraPresentation.from_json(read(artifacts['presentation.json']))
    assert restored_presentation.content_hash() == presentation.content_hash()
    restored_groebner = GroebnerCertificate.from_json(read(artifacts['groebner-certificate.json']))
    assert verify_groebner(presentation, restored_groebner).dimension == 4
    restored_basis = WordBasisCertificate.from_json(read(artifacts['word-basis-certificate.json']))
    algebra = verify_word_basis(presentation, restored_basis)
    assert algebra.dimension == 4
    assert restored_basis.basis_words == ((), (1,), (2,), (1, 2))
    metadata = read(artifacts['production.json'])
    assert metadata['completion']['basis_complete']
    assert metadata['dimension'] == 4
    assert metadata['reconstruction']['options']['shift'] == []
    assert metadata['solver']['store'] == 'samples.sqlite3'
    assert not metadata['solver']['hit']
    provenance = read(artifacts['provenance.json'])
    assert provenance['producer'] == report['producer']
    assert set(provenance['outputs']) == set(artifacts)-{'provenance.json'}
    assert all(file_identity(Path(path)) == report['files'][name] for name, path in artifacts.items())
    assert digest_bytes(canonical_json(report['files'])) == report['generation']


def test_cache_reuse_and_reconstruction_control_forwarding(setup, monkeypatch):
    presentation, config, result = setup
    calls = []
    def reconstruct(pres, sample, **options):
        assert pres is presentation
        assert sample.presentation_hash == presentation.content_hash()
        calls.append(options)
        return result
    monkeypatch.setattr(pipeline, 'reconstruct_groebner', reconstruct)
    first = pipeline.produce_algebra(presentation, {}, config)
    second = pipeline.produce_algebra(presentation, {}, config, max_primes=3, max_degree=9,
                                     threads=2, start_prime=1, max_products=12345,
                                     max_operations=12345, max_reductions=12345,
                                     max_words=30, max_nodes=100)
    assert first['verified'] and second['verified']
    # Different validation bounds deliberately change the numerical request key.
    third = pipeline.produce_algebra(presentation, {}, config)
    assert third['verified'] and third['solver']['hit']
    assert calls[1] == {'point': {}, 'max_nodes': 100, 'max_operations': 12345,
                       'max_products': 12345, 'max_reductions': 12345, 'max_words': 30,
                       'max_primes': 3, 'max_degree': 9, 'threads': 2, 'start_prime': 1,
                       'shift': None}


@pytest.mark.parametrize(('name', 'failure', 'stage'), [
    ('run_algebra_point', pipeline.ExecutionError('solver timed out'), 'solver'),
    ('solver_certificate', ArtifactError('retained proof changed'), 'reconstruction'),
    ('reconstruct_groebner', RuntimeError('support cannot lift'), 'reconstruction'),
    ('reconstruction_provenance', ImportError('native dependency unavailable'), 'reconstruction'),
    ('word_basis_from_groebner', RuntimeError('word bound exhausted'), 'actions'),
    ('_publish', OSError('output is full'), 'publication'),
])
def test_failed_stages_retain_an_inconclusive_result(setup, monkeypatch, name, failure, stage):
    presentation, config, _ = setup
    def fail(*args, **kwargs):
        raise failure
    monkeypatch.setattr(pipeline, name, fail)
    report = pipeline.produce_algebra(presentation, {}, config)
    assert_inconclusive(report, stage)
    assert report['reason'] == str(failure)
    assert read(config.output_root/'result.json') == report


def test_dependency_or_source_changes_prevent_publication(setup, monkeypatch):
    presentation, config, _ = setup
    dependencies = iter(({'version': 1}, {'version': 2}))
    monkeypatch.setattr(pipeline, 'reconstruction_provenance', lambda: next(dependencies))
    report = pipeline.produce_algebra(presentation, {}, config)
    assert_inconclusive(report, 'publication')
    assert 'dependencies changed' in report['reason']
    monkeypatch.setattr(pipeline, 'reconstruction_provenance', lambda: {})
    def changed(producer):
        raise ArtifactError('source changed')
    monkeypatch.setattr(pipeline, 'validate_source_snapshot', changed)
    report = pipeline.produce_algebra(presentation, {}, config)
    assert_inconclusive(report, 'publication')
    assert 'source changed' in report['reason']


def test_wrong_public_types_are_rejected_before_output(setup):
    presentation, config, _ = setup
    for pres, cfg in ((None, config), (presentation, None)):
        with pytest.raises(TypeError, match='AlgebraPresentation and RunConfig'):
            pipeline.produce_algebra(pres, {}, cfg)
    assert not config.output_root.exists()


def documents_from_report(report):
    return {name: read(path) for name, path in report['artifacts'].items() if name != 'provenance.json'}


def test_identical_immutable_generation_is_reused_but_corruption_is_rejected(setup):
    presentation, config, _ = setup
    report = pipeline.produce_algebra(presentation, {}, config)
    documents = documents_from_report(report)
    repeated = pipeline._publish(config.output_root, documents, report['producer'])
    assert repeated['generation'] == report['generation']
    certificate = Path(report['artifacts']['word-basis-certificate.json'])
    certificate.write_text('{}')
    with pytest.raises(ArtifactError, match='existing algebra generation'):
        pipeline._publish(config.output_root, documents, report['producer'])


def test_changed_serialized_document_is_rejected_before_selection(setup, monkeypatch):
    presentation, config, _ = setup
    actual_write = pipeline.write_json_atomic
    def corrupt(path, document):
        actual_write(path, {} if Path(path).name == 'presentation.json' else document)
    monkeypatch.setattr(pipeline, 'write_json_atomic', corrupt)
    report = pipeline.produce_algebra(presentation, {}, config)
    assert_inconclusive(report, 'publication')
    assert 'serialized algebra production differs' in report['reason']
    assert not any((config.output_root/'algebras').iterdir())


@pytest.mark.parametrize('winner', [False, True])
def test_atomic_generation_rename_requires_a_complete_winner(setup, monkeypatch, winner):
    presentation, config, _ = setup
    report = pipeline.produce_algebra(presentation, {}, config)
    documents = documents_from_report(report)
    shutil.rmtree(Path(report['artifacts']['production.json']).parent)
    original = pipeline.os.replace
    def competing(source, destination):
        if Path(source).name == 'generation':
            if winner:
                shutil.copytree(source, destination)
            raise OSError('competing rename')
        return original(source, destination)
    monkeypatch.setattr(pipeline.os, 'replace', competing)
    if winner:
        assert pipeline._publish(config.output_root, documents, report['producer'])['generation'] == report['generation']
    else:
        with pytest.raises(OSError, match='competing rename'):
            pipeline._publish(config.output_root, documents, report['producer'])


def test_interrupted_attempt_cannot_leave_an_old_verified_verdict(setup, monkeypatch):
    presentation, config, _ = setup
    completed = pipeline.produce_algebra(presentation, {}, config)
    assert completed['verified']
    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt
    monkeypatch.setattr(pipeline, 'reconstruct_groebner', interrupt)
    with pytest.raises(KeyboardInterrupt):
        pipeline.produce_algebra(presentation, {}, config)
    unfinished = read(config.output_root/'result.json')
    assert unfinished['status'] == 'running'
    assert not unfinished['verified'] and unfinished['dimension'] is None
    assert 'artifacts' not in unfinished and 'generation' not in unfinished
    # Previously selected immutable evidence survives, but is no longer the
    # current attempt's result.
    assert Path(completed['artifacts']['word-basis-certificate.json']).is_file()


def test_serialized_presentation_must_match_both_certificates(setup):
    presentation, config, _ = setup
    completed = pipeline.produce_algebra(presentation, {}, config)
    documents = documents_from_report(completed)
    changed = AlgebraPresentation(presentation.algebra, presentation.relations,
                                   ('new-label', 'vv', 'vu'))
    documents['presentation.json'] = changed.to_json()
    with pytest.raises(ArtifactError, match='saved presentation'):
        pipeline._publish(config.output_root, documents, completed['producer'])


def test_literal_parameter_pole_is_reported_as_inconclusive(setup):
    _, config, _ = setup
    field = RationalFunctionField(('p',))
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    p = field.parameter('p')
    presentation = AlgebraPresentation(parent, (x*x-parent.scalar(1/(p-1)),), ('pole',))
    report = pipeline.produce_algebra(presentation, {'p': 1}, config)
    assert_inconclusive(report, 'solver')
    assert 'division by zero' in report['reason']
    assert read(config.output_root/'result.json') == report
