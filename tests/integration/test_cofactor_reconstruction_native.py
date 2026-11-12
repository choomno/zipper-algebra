"""Tiny native particular-witness reconstructions with exact original-input checks."""
from fractions import Fraction
import json

import pytest

from zipper_algebra.algebra import AlgebraElement, AlgebraPresentation, CofactorIdentity, FreeAlgebra
from zipper_algebra.cofactor.general import verify_cofactor
from zipper_algebra.cofactor.reconstruction import reconstruct_derivation
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.matrices import ExactMatrix
from zipper_algebra.reconstruction import reconstruct_kernel
from zipper_algebra.verification.derivations import DerivationGraph, DerivationTerm
from zipper_algebra.verification.linear import verify_kernel


@pytest.fixture
def native():
    pytest.importorskip('sympy', reason='parameter arithmetic is optional')
    return pytest.importorskip('fflow', reason='FiniteFlow native Python bindings are optional')


def inputs(field, relations, target, contexts):
    parent = target.algebra
    presentation = AlgebraPresentation(parent, tuple(relations),
        tuple(f'r{i}' for i in range(len(relations))))
    # These coefficients deliberately do not encode the desired answer.
    graph = DerivationGraph(presentation.content_hash(), (target,),
        (tuple(DerivationTerm(field.scalar(7), left, 'input', index, right)
               for left, index, right in contexts),))
    return presentation, graph


def collected(identity):
    result = {}
    for left, index, right in identity.terms:
        assert len(left.terms) == len(right.terms) == 1
        (u, a), (v, b) = left.terms[0], right.terms[0]
        result[u, index, v] = a*b
    return result


def test_native_rational_coefficient_uses_positive_affine_rhs(native):
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    presentation, graph = inputs(field, (Fraction(2, 3)*x,), x, (((), 0, ()),))
    result = reconstruct_derivation(presentation, graph, 0)
    assert collected(result.identity) == {((), 0, ()): field.scalar(Fraction(3, 2))}
    assert verify_cofactor(presentation, result.identity)


def test_native_parameter_noncommuting_identity_survives_json_save_reload(native, tmp_path):
    field = RationalFunctionField(('s', 't'))
    s, t = map(field.parameter, field.parameters)
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = map(parent.generator, parent.generators)
    relations = (x*y-s*y*x, y*y)
    a, b = s/(t+1), Fraction(2, 3)*(t+1)
    target = a*x*relations[0]*y+b*y*relations[1]*x
    presentation, graph = inputs(field, relations, target,
        (((1,), 0, (2,)), ((2,), 1, (1,))))
    result = reconstruct_derivation(presentation, graph, 0, max_primes=3,
        max_degree=8, start_prime=1, shift=(2, 3), threads=1)
    assert collected(result.identity) == {((1,), 0, (2,)): a, ((2,), 1, (1,)): b}
    paths = {name: tmp_path/f'{name}.json' for name in ('presentation', 'target', 'identity')}
    for name, value in [('presentation', presentation), ('target', target), ('identity', result.identity)]:
        paths[name].write_text(json.dumps(value.to_json()))
    saved_presentation = AlgebraPresentation.from_json(json.loads(paths['presentation'].read_text()))
    saved_target = AlgebraElement.from_json(json.loads(paths['target'].read_text()))
    saved_identity = CofactorIdentity.from_json(json.loads(paths['identity'].read_text()))
    assert saved_identity.target == saved_target == target
    assert verify_cofactor(saved_presentation, saved_identity)


def test_native_underdetermined_vector_maps_unsorted_ids_and_zero_free_coordinates(native, monkeypatch):
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = map(parent.generator, parent.generators)
    presentation, graph = inputs(field, (x, y, x+y), x+y,
        (((), 0, ()), ((), 1, ()), ((), 2, ())))
    observed = []
    original = native.LSolveDepVars
    def record(*args):
        ids = original(*args)
        observed.append(list(ids))
        return ids
    monkeypatch.setattr(native, 'LSolveDepVars', record)
    result = reconstruct_derivation(presentation, graph, 0)
    assert verify_cofactor(presentation, result.identity)
    assert observed and set(observed[0]) == {0, 1, 2}
    assert observed[0] != sorted(observed[0])
    assert collected(result.identity) == {((), 0, ()): field.one, ((), 1, ()): field.one}
    assert result.statistics['candidate_contexts'] == 3
    assert result.statistics['contexts'] == result.statistics['cofactor_terms'] == 2


def test_candidate_mode_does_not_change_certified_complete_kernel_in_same_process(native):
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = map(parent.generator, parent.generators)
    presentation, graph = inputs(field, (x, y), x+y, (((), 0, ()), ((), 1, ())))
    first = reconstruct_derivation(presentation, graph, 0)
    matrix = ExactMatrix.from_rows(field, [[1, 1]])
    kernel = reconstruct_kernel(matrix)
    assert verify_kernel(matrix, kernel) is kernel
    assert kernel.kernel.nrows == 2 and kernel.kernel.ncols == 1
    second = reconstruct_derivation(presentation, graph, 0)
    assert verify_cofactor(presentation, first.identity)
    assert verify_cofactor(presentation, second.identity)
    assert collected(first.identity) == collected(second.identity)


def test_native_learning_prime_pole_fails_and_next_prime_recovers(native):
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    prime = native.PrimeNo(0)
    presentation, graph = inputs(field, (Fraction(1, prime)*x,), x, (((), 0, ()),))
    with pytest.raises(RuntimeError):
        reconstruct_derivation(presentation, graph, 0)
    recovered = reconstruct_derivation(presentation, graph, 0, start_prime=1, max_primes=3)
    assert verify_cofactor(presentation, recovered.identity)
    assert collected(recovered.identity) == {((), 0, ()): field.scalar(prime)}
