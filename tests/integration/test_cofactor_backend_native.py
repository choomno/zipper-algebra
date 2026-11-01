"""Small real-native queries and the installed original-input export recipe."""
import hashlib
import json

import pytest

from zipper_algebra.algebra import AlgebraElement, CofactorIdentity
from zipper_algebra.cofactor.general import cofactor_representation, verify_cofactor
from zipper_algebra.verification.finite import WordBasisCertificate
from tests.integration.test_general_cofactor_example import inputs, run


@pytest.fixture(autouse=True)
def native_dependencies():
    pytest.importorskip('sympy')
    pytest.importorskip('fflow')


@pytest.mark.parametrize('membership', [False, True])
def test_native_query_and_saved_problem_have_an_exact_original_input_cofactor(tmp_path, membership):
    p, u, v, paths = inputs(tmp_path, ('alpha', 'unused'))
    a = p.algebra.field.parameter('alpha')
    target = u*u*u+v*u*v
    expected = a*u+v
    if membership:
        target = target-expected
        expected = p.algebra.zero
    direct = cofactor_representation(p,
        WordBasisCertificate.from_json(json.loads(paths['word-basis-certificate.json'].read_text())),
        target, backend='finiteflow', max_primes=3, max_degree=6, start_prime=1, shift=(3, 5))
    assert direct.remainder == expected
    assert direct.identity.target == target-expected and verify_cofactor(p, direct.identity)
    paths['target.json'].write_text(json.dumps(target.to_json()))
    output = tmp_path/'native-output'
    result, report = run(paths, output, '--backend', 'finiteflow', '--max-primes', '3',
                         '--max-degree', '6', '--start-prime', '1', '--shift', '3', '5')
    assert result.returncode == 0, result.stderr+result.stdout
    assert report['verified'] and report['backend'] == 'finiteflow'
    assert report['in_ideal'] is membership
    representative = AlgebraElement.from_json(json.loads((output/'representative.json').read_text()))
    identity = CofactorIdentity.load(output/'cofactor.json')
    assert representative == expected and identity.target == target-expected
    assert verify_cofactor(p, identity)
    reconstruction = report['reconstruction']
    assert reconstruction['options']['max_primes'] == 3
    assert reconstruction['options']['shift'] == [3, 5]
    assert reconstruction['dependencies']['implementation'] == 'FiniteFlow'
    assert set(reconstruction['dependencies']['files']) == {'binding', 'extension', 'library'}
    provenance = json.loads((output/'provenance.json').read_text())
    assert set(provenance['outputs']) == {'presentation.json', 'target.json', 'representative.json',
                                          'cofactor.json', 'report.json'}
    for name, recorded in provenance['outputs'].items():
        assert hashlib.sha256((output/name).read_bytes()).hexdigest() == recorded['sha256']
