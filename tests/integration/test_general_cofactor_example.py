"""A saved general cofactor proves the literal target-minus-representative identity."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from zipper_algebra.action_production import word_basis_from_rewrites
from zipper_algebra.algebra import AlgebraElement, AlgebraPresentation, CofactorIdentity, FreeAlgebra
from zipper_algebra.cofactor.general import verify_cofactor
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm

ROOT = Path(__file__).resolve().parents[2]


def inputs(tmp_path, parameters=()):
    field = RationalFunctionField(parameters)
    parent = FreeAlgebra(field, ('u', 'v'))
    u, v = parent.generator('u'), parent.generator('v')
    a = field.parameter('alpha') if parameters else field.one
    presentation = AlgebraPresentation(parent, (u*u-a, v*v, v*u+u*v-1), ('u2', 'v2', 'mixed'))
    proposal = GroebnerCertificate(presentation.content_hash(), presentation.relations,
        tuple((ProofTerm(field.one, (), 'input', index, ()),) for index in range(3)))
    certificate = word_basis_from_rewrites(presentation, proposal)
    paths = {name: tmp_path / name for name in ('presentation.json', 'word-basis-certificate.json', 'target.json')}
    paths['presentation.json'].write_text(json.dumps(presentation.to_json()))
    paths['word-basis-certificate.json'].write_text(json.dumps(certificate.to_json()))
    return presentation, u, v, paths


def run(paths, output, *options):
    command = [sys.executable, str(ROOT / 'examples/general_cofactor.py'),
               '--presentation', str(paths['presentation.json']),
               '--word-basis', str(paths['word-basis-certificate.json']),
               '--target', str(paths['target.json']), '--output', str(output), *options]
    result = subprocess.run(command, text=True, capture_output=True, timeout=30)
    return result, json.loads(result.stdout)


@pytest.mark.parametrize('parameters', [(), ('alpha', 'unused')])
@pytest.mark.parametrize('membership', [False, True])
def test_saved_cofactor_proves_noncommuting_reduction(tmp_path, parameters, membership):
    if parameters:
        pytest.importorskip('sympy')
    presentation, u, v, paths = inputs(tmp_path, parameters)
    target = v*u*v-v if membership else v*u*v+u
    expected = presentation.algebra.zero if membership else v+u
    paths['target.json'].write_text(json.dumps(target.to_json()))
    result, report = run(paths, tmp_path/'out')
    assert result.returncode == 0, result.stderr + result.stdout
    assert report['status'] == 'verified' and report['verified']
    assert report['in_ideal'] is membership
    assert report['presentation_hash'] == presentation.content_hash()
    representative = AlgebraElement.from_json(json.loads(Path(report['artifacts']['representative.json']).read_text()))
    identity = CofactorIdentity.from_json(json.loads(Path(report['artifacts']['cofactor.json']).read_text()))
    assert representative == expected
    assert identity.target == target-expected
    assert identity.terms and verify_cofactor(presentation, identity)
    provenance = json.loads((tmp_path/'out/provenance.json').read_text())
    assert provenance['producer'] == report['producer']
    for name, path in paths.items():
        source = path.read_bytes()
        assert {'name': name, 'bytes': len(source), 'sha256': hashlib.sha256(source).hexdigest()} in report['inputs']
    for name, identity in provenance['outputs'].items():
        content = (tmp_path/'out'/name).read_bytes()
        assert identity == {'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()}
    assert (tmp_path/'out/.gitignore').read_text() == '*\n!.gitignore\n'


def test_exhausted_expansion_budget_does_not_claim_membership(tmp_path):
    presentation, u, v, paths = inputs(tmp_path)
    paths['target.json'].write_text(json.dumps((v*u*v-v).to_json()))
    result, report = run(paths, tmp_path/'out', '--max-operations', '1')
    assert result.returncode == 3
    assert report['status'] == 'inconclusive' and not report['verified']
    assert 'in_ideal' not in report and 'artifacts' not in report
    assert not (tmp_path/'out/cofactor.json').exists()
    assert 'budget' in report['reason']
    assert json.loads((tmp_path/'out/report.json').read_text()) == report
