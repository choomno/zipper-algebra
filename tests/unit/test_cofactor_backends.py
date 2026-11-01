"""Backend selection preserves the checked representative and original-input identity."""
from dataclasses import FrozenInstanceError
from fractions import Fraction
import json

import pytest

from zipper_algebra.cofactor.general import cofactor_representation, expand_derivation, verify_cofactor
from zipper_algebra.cofactor.reconstruction import CofactorReconstruction
from zipper_algebra.algebra import AlgebraPresentation, CofactorIdentity, FreeAlgebra
from zipper_algebra.fields import FieldElement, RationalFunctionField
from zipper_algebra.verification.working_field import WorkingField
from tests.unit.test_finite_algebra import matrix_algebra_certificate
from tests.unit.test_general_cofactor import simple_graph


def test_optional_backend_receives_limits_and_returns_a_checked_input_identity(monkeypatch):
    import zipper_algebra.cofactor.reconstruction as reconstruction
    p, graph, _, _ = simple_graph()
    expected = expand_derivation(p, graph, 0)
    calls = []

    def reconstruct(presentation, proposal, node, **options):
        calls.append((presentation, proposal, node, options))
        return CofactorReconstruction(expected, (("contexts", len(expected.terms)),))

    monkeypatch.setattr(reconstruction, 'reconstruct_derivation', reconstruct)
    actual = expand_derivation(p, graph, 0, backend='finiteflow', max_operations=500,
        max_contexts=200, max_nodes=100, max_matrix_entries=300, max_primes=3,
        max_degree=8, threads=2, start_prime=1, shift=())
    assert actual == expected and verify_cofactor(p, actual)
    assert calls == [(p, graph, 0, {'max_operations': 500, 'max_contexts': 200, 'max_nodes': 100,
        'max_matrix_entries': 300, 'max_primes': 3, 'max_degree': 8,
        'threads': 2, 'start_prime': 1, 'shift': ()})]


def test_query_builds_the_same_target_minus_normal_form_for_either_backend(monkeypatch):
    import zipper_algebra.cofactor.reconstruction as reconstruction
    p, certificate, x, y = matrix_algebra_certificate()
    target = y*x*y+x
    expected = cofactor_representation(p, certificate, target)
    calls = []

    def reconstruct(presentation, graph, node, **options):
        assert graph.nodes[node] == target-(y+x)
        calls.append(options)
        identity = expand_derivation(presentation, graph, node)
        return CofactorReconstruction(identity, (("contexts", len(identity.terms)),))

    monkeypatch.setattr(reconstruction, 'reconstruct_derivation', reconstruct)
    result = cofactor_representation(p, certificate, target, backend='finiteflow',
                                    max_matrix_entries=300, max_primes=2)
    assert result.remainder == expected.remainder == y+x
    assert result.identity.target == target-result.remainder
    assert verify_cofactor(p, result.identity)
    assert calls[0]['max_matrix_entries'] == 300 and calls[0]['max_primes'] == 2
    assert result.statistics['basis_dimension'] == 4


@pytest.mark.parametrize('backend', [None, False, 'unknown', 'FINITEFLOW'])
def test_unknown_backend_is_rejected_before_computation(backend):
    p, graph, _, _ = simple_graph()
    with pytest.raises(ValueError, match='backend'):
        expand_derivation(p, graph, 0, backend=backend)
    p, certificate, x, _ = matrix_algebra_certificate()
    with pytest.raises(ValueError, match='backend'):
        cofactor_representation(p, certificate, x, backend=backend)


def test_optional_backend_checks_an_immediate_zero_without_importing_native_code(monkeypatch):
    import zipper_algebra.cofactor.reconstruction as reconstruction
    p, certificate, x, y = matrix_algebra_certificate()

    def unexpected(*args, **kwargs):
        pytest.fail('a zero ideal difference must not need native solving')

    monkeypatch.setattr(reconstruction, '_native_candidate', unexpected)
    result = cofactor_representation(p, certificate, x+y, backend='finiteflow')
    assert result.remainder == x+y
    assert result.identity.target.is_zero() and not result.identity.terms
    assert verify_cofactor(p, result.identity)


def test_reconstruction_result_is_immutable_and_statistics_are_snapshots():
    p, graph, _, _ = simple_graph()
    identity = expand_derivation(p, graph, 0)
    result = CofactorReconstruction(identity, (("contexts", 2),))
    with pytest.raises(FrozenInstanceError):
        result.identity = identity
    result.statistics['contexts'] = 999
    assert result.statistics == {'contexts': 2}


def test_constant_denominator_json_keeps_exact_rational_input_cofactors():
    field = RationalFunctionField()
    half = FieldElement(field, {(): 1}, {(): 2})
    assert half.denominator == field.one.denominator
    assert WorkingField(field).from_element(half) == Fraction(1, 2)
    encoded = half.to_json()
    # This valid input encoding has not yet undergone constructor normalization.
    encoded['numerator'], encoded['denominator'] = [[[], 1]], [[[], 2]]
    loaded = FieldElement.from_json(json.loads(json.dumps(encoded)))
    assert loaded.denominator == field.one.denominator
    assert WorkingField(field).from_element(loaded) == Fraction(1, 2)
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    p = AlgebraPresentation(parent, (x,), ('x',))
    false = CofactorIdentity(p.content_hash(), x, ((parent.scalar(loaded), 0, parent.one),))
    assert not verify_cofactor(p, CofactorIdentity.from_json(json.loads(json.dumps(false.to_json()))))
    true = CofactorIdentity(p.content_hash(), parent.scalar(loaded)*x, false.terms)
    assert verify_cofactor(p, true)
