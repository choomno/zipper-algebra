"""Independent matrix-unit controls for certified finite word bases.

The matrices and every column proof below are specified directly. No Groebner
basis or reduction producer supplies the expected coordinates.
"""
from copy import deepcopy
from dataclasses import FrozenInstanceError
from fractions import Fraction
import json

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.derivations import (
    DerivationGraph, DerivationTerm, verify_derivations,
)
from zipper_algebra.verification.finite import (
    GeneratorAction, WordBasisCertificate, verify_word_basis,
)


def presentation(parent, relations):
    return AlgebraPresentation(parent, tuple(relations),
                               tuple(f'r{i}' for i in range(len(relations))))


def sparse_action(field, columns):
    return GeneratorAction(field, tuple(tuple((row, field.scalar(value)) for row, value in col)
                                        for col in columns))


def term(field, index, *, left=(), right=(), value=1, kind='input'):
    return DerivationTerm(field.scalar(value), left, kind, index, right)


def matrix_algebra_certificate(characteristic=0):
    field = RationalFunctionField((), characteristic)
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = (parent.generator(name) for name in parent.generators)
    relations = (x*x, y*y, y*x+x*y-1)
    pres = presentation(parent, relations)
    # x=E12 and y=E21. In B=(1,x,y,xy), left multiplication has columns:
    # xB=(x,0,xy,0), yB=(y,1-xy,0,y).
    mx = sparse_action(field, (((1, 1),), (), ((3, 1),), ()))
    my = sparse_action(field, (((2, 1),), ((0, 1), (3, -1)), (), ((2, 1),)))
    nodes = (parent.zero, x*x, x*x*y, relations[2], y*y, y*x*y-y)
    graph = DerivationGraph(pres.content_hash(), nodes, (
        (), (term(field, 0),), (term(field, 0, right=(2,)),),
        (term(field, 2),), (term(field, 1),),
        (term(field, 2, right=(2,)), term(field, 1, left=(1,), value=-1)),
    ))
    certificate = WordBasisCertificate(pres.content_hash(), ((), (1,), (2,), (1, 2)),
                                       (mx, my), graph, ((0, 1, 0, 2), (0, 3, 4, 5)))
    return pres, certificate, x, y


def quadratic_certificate(characteristic=0):
    field = RationalFunctionField(('s',), characteristic)
    parent = FreeAlgebra(field, ('x',))
    x, s = parent.generator('x'), field.parameter('s')
    relation = x*x-s
    pres = presentation(parent, (relation,))
    matrix = sparse_action(field, (((1, 1),), ((0, s),)))
    graph = DerivationGraph(pres.content_hash(), (parent.zero, relation),
                            ((), (term(field, 0),)))
    certificate = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                       (matrix,), graph, ((0, 1),))
    return pres, certificate, x, s


@pytest.mark.parametrize('characteristic', [0, 2, 7])
def test_matrix_unit_algebra_has_an_exact_noncommutative_four_word_basis(characteristic):
    pres, certificate, x, y = matrix_algebra_certificate(characteristic)
    result = verify_word_basis(pres, certificate)
    field = pres.algebra.field
    assert result.dimension == 4
    assert result.basis_words == ((), (1,), (2,), (1, 2))
    assert result.coordinates(x*y) == tuple(field.scalar(value) for value in (0, 0, 0, 1))
    assert result.coordinates(y*x) == tuple(field.scalar(value) for value in (1, 0, 0, -1))
    assert result.coordinates(x*y*x) == tuple(field.scalar(value) for value in (0, 1, 0, 0))
    assert result.coordinates(x*x) == (field.zero,)*4
    assert result.coordinates(pres.algebra.one) == (field.one, field.zero, field.zero, field.zero)
    assert result.coordinates(pres.algebra.zero) == (field.zero,)*4
    assert result.generator_action('x') == certificate.actions[0]
    assert result.generator_action('y') == certificate.actions[1]


@pytest.mark.parametrize('characteristic', [0, 3])
def test_symbolic_quadratic_word_basis_retains_the_declared_parameter_field(characteristic):
    pytest.importorskip('sympy')
    pres, certificate, x, s = quadratic_certificate(characteristic)
    result = verify_word_basis(pres, certificate)
    field = pres.algebra.field
    assert result.dimension == 2
    assert result.coordinates(x**5) == (field.zero, s*s)
    assert result.coordinates(x**4) == (s*s, field.zero)
    assert result.coordinates(x*x-s) == (field.zero, field.zero)
    assert result.coordinates(Fraction(1, 2)*x) == (field.zero, field.scalar(Fraction(1, 2)))


def test_column_proofs_are_essential_even_for_a_cyclic_representation():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    free = presentation(parent, ())
    # The one-dimensional representation x->0 is cyclic and satisfies every
    # defining relation of the free algebra. Its upper-bound claim x=0 is false.
    graph = DerivationGraph(free.content_hash(), (x,), ((),))
    false = WordBasisCertificate(free.content_hash(), ((),),
                                 (sparse_action(field, ((),)),), graph, ((0,),))
    with pytest.raises(ValueError):
        verify_word_basis(free, false)


def test_valid_spanning_proofs_and_a_representation_do_not_establish_independence():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    pres = presentation(parent, (x,))
    # Q<x>/(x) is one-dimensional. The named words 1,x are distinct free words
    # but dependent in the quotient. Both action-column identities are honest.
    graph = DerivationGraph(pres.content_hash(), (x, x*x),
                            ((term(field, 0),), (term(field, 0, left=(1,)),)))
    false = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                 (sparse_action(field, ((), ())),), graph, ((0, 1),))
    with pytest.raises(ValueError):
        verify_word_basis(pres, false)


def test_corrupting_one_action_column_invalidates_the_bound_proof_target():
    pres, cert, _, _ = matrix_algebra_certificate()
    field = pres.algebra.field
    wrong = sparse_action(field, (((0, 1), (1, 1)), (), ((3, 1),), ()))
    false = WordBasisCertificate(cert.presentation_hash, cert.basis_words,
                                 (wrong, cert.actions[1]), cert.proof, cert.action_nodes)
    with pytest.raises(ValueError):
        verify_word_basis(pres, false)


def test_transposed_or_reversed_actions_cannot_pass_the_noncommutative_control():
    pres, cert, _, _ = matrix_algebra_certificate()
    reversed_actions = WordBasisCertificate(cert.presentation_hash, cert.basis_words,
                                            cert.actions[::-1], cert.proof, cert.action_nodes)
    with pytest.raises(ValueError):
        verify_word_basis(pres, reversed_actions)


def test_every_declared_generator_including_an_inverse_requires_an_action():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = (parent.generator(name) for name in parent.generators)
    pres = presentation(parent, (x*x-1, x*y-1, y*x-1))
    action = sparse_action(field, (((1, 1),), ((0, 1),)))
    graph = DerivationGraph(pres.content_hash(), (parent.zero, x*x-1, y-x, y*x-1), (
        (), (term(field, 0),),
        (term(field, 0, left=(2,), value=-1), term(field, 2, right=(1,))),
        (term(field, 2),),
    ))
    certificate = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                       (action, action), graph, ((0, 1), (2, 3)))
    assert verify_word_basis(pres, certificate).dimension == 2
    with pytest.raises((TypeError, ValueError)):
        missing = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                       (action,), graph, ((0, 1),))
        verify_word_basis(pres, missing)


def test_certificate_is_bound_to_literal_relations_and_parent_even_for_zero():
    pres, cert, _, _ = matrix_algebra_certificate()
    other = AlgebraPresentation(pres.algebra, (2*pres.relations[0],)+pres.relations[1:], pres.labels)
    with pytest.raises(ValueError):
        verify_word_basis(other, cert)
    result = verify_word_basis(pres, cert)
    foreign = FreeAlgebra(RationalFunctionField(('s',)), pres.algebra.generators)
    with pytest.raises(TypeError):
        result.coordinates(foreign.zero)
    with pytest.raises((TypeError, ValueError)):
        result.coordinates(0)
    with pytest.raises(ValueError):
        result.generator_action('unknown')


def test_zero_algebra_has_a_separate_verified_unit_witness():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    pres = presentation(parent, (parent.scalar(2),))
    graph = DerivationGraph(pres.content_hash(), (parent.one,),
                            ((term(field, 0, value=Fraction(1, 2)),),))
    cert = WordBasisCertificate(pres.content_hash(), (),
                                (sparse_action(field, ()),), graph, ((),), unit_node=0)
    result = verify_word_basis(pres, cert)
    assert result.dimension == 0
    assert result.coordinates(parent.generator('x')) == ()
    assert result.coordinates(parent.one) == ()
    assert result.basis_words == ()
    with pytest.raises((TypeError, ValueError)):
        false = WordBasisCertificate(pres.content_hash(), (),
                                     (sparse_action(field, ()),), graph, ((),))
        verify_word_basis(pres, false)


def test_scalar_algebra_requires_no_generator_actions_but_keeps_the_unit_word():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ())
    pres = presentation(parent, ())
    graph = DerivationGraph(pres.content_hash(), (), ())
    cert = WordBasisCertificate(pres.content_hash(), ((),), (), graph, ())
    result = verify_word_basis(pres, cert)
    assert result.dimension == 1
    assert result.coordinates(parent.scalar(7)) == (field.scalar(7),)


def test_derivation_nodes_share_forward_references_without_expanding_the_proof():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    pres = presentation(parent, (x*x,))
    graph = DerivationGraph(pres.content_hash(), (x*x, 2*x*x, parent.zero), (
        (term(field, 1, value=Fraction(1, 2), kind='node'),),
        (term(field, 0, value=2),), (),
    ))
    verify_derivations(pres, graph)
    cert = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                (sparse_action(field, (((1, 1),), ())),), graph, ((2, 0),))
    assert verify_word_basis(pres, cert).dimension == 2


@pytest.mark.parametrize('zero_reference', [False, True])
def test_graph_cycles_and_dangling_references_are_checked_before_zero_cancellation(zero_reference):
    pres, cert, _, _ = matrix_algebra_certificate()
    field = pres.algebra.field
    rows = list(cert.proof.derivations)
    if zero_reference:
        rows[0] = (term(field, 999, value=0, kind='node'),)
    else:
        rows[0] = (term(field, 0, kind='node'),)
    graph = DerivationGraph(pres.content_hash(), cert.proof.nodes, tuple(rows))
    false = WordBasisCertificate(cert.presentation_hash, cert.basis_words,
                                 cert.actions, graph, cert.action_nodes)
    with pytest.raises(ValueError):
        verify_word_basis(pres, false)


def test_multiplication_uses_left_generator_actions_in_the_written_word_order():
    pres, certificate, _, _ = matrix_algebra_certificate()
    result = verify_word_basis(pres, certificate)
    unit = result.element([1, 0, 0, 0])
    x = result.element([0, 1, 0, 0])
    y = result.element([0, 0, 1, 0])
    xy = result.element([0, 0, 0, 1])
    yx = result.element([1, 0, 0, -1])
    assert result.multiply(x, y) == xy
    assert result.multiply(y, x) == yx
    assert result.multiply(xy, xy) == xy
    assert result.multiply(unit, x) == x
    assert result.multiply(y, unit) == y
    assert result.multiply(result.element([0]*4), unit) == result.element([0]*4)
    other, other_certificate, _, _ = matrix_algebra_certificate(characteristic=7)
    foreign = verify_word_basis(other, other_certificate).element([0]*4)
    with pytest.raises(TypeError):
        result.multiply(x, foreign)
    with pytest.raises(TypeError):
        result.multiply(foreign, x)


@pytest.mark.parametrize('words', [
    ((1,),),
    ((), (), (2,), (1, 2)),
    ((), (1,), (2,), (3,)),
    [(), (1,), (2,), (1, 2)],
])
def test_invalid_basis_word_lists_do_not_define_a_certificate(words):
    pres, certificate, _, _ = matrix_algebra_certificate()
    with pytest.raises((TypeError, ValueError)):
        false = WordBasisCertificate(certificate.presentation_hash, words,
                                     certificate.actions, certificate.proof, certificate.action_nodes)
        verify_word_basis(pres, false)


def test_reordering_basis_words_requires_reordering_both_actions_and_proofs():
    pres, certificate, _, _ = matrix_algebra_certificate()
    false = WordBasisCertificate(certificate.presentation_hash, ((), (2,), (1,), (1, 2)),
                                 certificate.actions, certificate.proof, certificate.action_nodes)
    with pytest.raises(ValueError):
        verify_word_basis(pres, false)


@pytest.mark.parametrize('option,value', [
    ('max_operations', 0), ('max_operations', True), ('max_operations', 1.5),
    ('max_dimension', 0), ('max_dimension', True), ('max_dimension', 1.5),
    ('max_nodes', 0), ('max_nodes', True), ('max_nodes', 1.5),
])
def test_verification_limits_are_exact_positive_integers(option, value):
    pres, certificate, _, _ = matrix_algebra_certificate()
    with pytest.raises((TypeError, ValueError)):
        verify_word_basis(pres, certificate, **{option: value})


def test_word_basis_and_derivation_graph_json_round_trip_exact_parameter_data():
    pytest.importorskip('sympy')
    pres, certificate, x, s = quadratic_certificate()
    restored = WordBasisCertificate.from_json(json.loads(json.dumps(certificate.to_json())))
    assert restored.to_json() == certificate.to_json()
    assert verify_word_basis(pres, restored).coordinates(x**4) == (s*s, pres.algebra.field.zero)
    graph = DerivationGraph.from_json(json.loads(json.dumps(certificate.proof.to_json())))
    assert graph.to_json() == certificate.proof.to_json()
    with pytest.raises(FrozenInstanceError):
        graph.nodes = ()
    with pytest.raises(FrozenInstanceError):
        restored.basis_words = ()


@pytest.mark.parametrize('value', [None, [], 'certificate', {}])
def test_json_certificate_boundaries_require_objects_with_a_schema(value):
    with pytest.raises((TypeError, ValueError)):
        WordBasisCertificate.from_json(value)
    with pytest.raises((TypeError, ValueError)):
        DerivationGraph.from_json(value)
    with pytest.raises((TypeError, ValueError)):
        GeneratorAction.from_json(value)


@pytest.mark.parametrize('mutation', [
    lambda data: data.update(schema='other'),
    lambda data: data.update(extra='ignored?'),
    lambda data: data.pop('presentation_hash'),
    lambda data: data.update(basis_words={}),
    lambda data: data.update(actions={}),
    lambda data: data.update(action_nodes={}),
    lambda data: data.update(proof={}),
])
def test_corrupt_word_basis_json_is_not_silently_coerced(mutation):
    _, certificate, _, _ = matrix_algebra_certificate()
    data = deepcopy(certificate.to_json())
    mutation(data)
    with pytest.raises((TypeError, ValueError)):
        WordBasisCertificate.from_json(data)


def test_zero_dimensional_matrices_without_a_unit_ideal_proof_do_not_certify_zero_algebra():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    pres = presentation(parent, ())
    graph = DerivationGraph(pres.content_hash(), (parent.zero,), ((),))
    false = WordBasisCertificate(pres.content_hash(), (),
                                 (sparse_action(field, ()),), graph, ((),), unit_node=0)
    with pytest.raises(ValueError):
        verify_word_basis(pres, false)


def test_finite_elements_of_different_quotients_cannot_mix_even_at_zero():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    quotients = []
    for constant in (0, 1):
        relation = x*x-constant
        pres = presentation(parent, (relation,))
        graph = DerivationGraph(pres.content_hash(), (parent.zero, relation),
                                ((), (term(field, 0),)))
        action = sparse_action(field, (((1, 1),), ((0, constant),)))
        cert = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                    (action,), graph, ((0, 1),))
        quotients.append(verify_word_basis(pres, cert))
    left, right = quotients
    with pytest.raises(TypeError):
        left.multiply(left.element([0, 0]), right.element([1, 0]))


def test_cyclic_actions_with_honest_spanning_proofs_must_still_satisfy_input_relations():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    pres = presentation(parent, (parent.one,))
    # In the zero algebra every proposed column has an ideal proof. The swap
    # matrix is cyclic on 1,x, but it cannot represent the relation 1=0.
    graph = DerivationGraph(pres.content_hash(), (parent.zero, x*x-1), (
        (), (term(field, 0, left=(1,), right=(1,)), term(field, 0, value=-1)),
    ))
    action = sparse_action(field, (((1, 1),), ((0, 1),)))
    false = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                 (action,), graph, ((0, 1),))
    with pytest.raises(ValueError):
        verify_word_basis(pres, false)


def test_finite_elements_preserve_their_verified_parent_and_exact_coordinates():
    pres, certificate, _, _ = matrix_algebra_certificate()
    result = verify_word_basis(pres, certificate)
    x = result.element([0, 1, 0, 0])
    y = result.element([0, 0, 1, 0])
    assert x.algebra is result
    assert x.coordinates == tuple(pres.algebra.field.scalar(value) for value in (0, 1, 0, 0))
    assert x+y == result.element([0, 1, 1, 0])
    assert (x+y)-x == y
    assert -x == result.element([0, -1, 0, 0])
    assert x*y == result.element([0, 0, 0, 1])
    assert 2*x == x*2 == result.element([0, 2, 0, 0])
    assert x+1 == 1+x == result.element([1, 1, 0, 0])
    assert 1-x == result.element([1, -1, 0, 0])
    assert x-x == result.zero
    assert x*result.one == x
    assert result.element(iter([0, 1, 0, 0])) == x
    with pytest.raises((TypeError, ValueError)):
        result.element([0, 1, 0])
    with pytest.raises(TypeError):
        result.element([False, 0, 0, 0])
    with pytest.raises(TypeError):
        result.element([0.0, 0, 0, 0])
    with pytest.raises(TypeError):
        result.element([RationalFunctionField(('s',)).zero, 0, 0, 0])
    with pytest.raises(FrozenInstanceError):
        x.coordinates = ()
    with pytest.raises(FrozenInstanceError):
        result.basis_words = ()
    report = result.report
    report['dimension'] = 1000
    assert result.report['dimension'] == 4
    with pytest.raises(TypeError):
        result.coordinates(FreeAlgebra(pres.algebra.field, ('y', 'x')).zero)


@pytest.mark.parametrize('option,value', [('max_dimension', 3), ('max_nodes', 1), ('max_operations', 1)])
def test_resource_exhaustion_cannot_return_a_verified_dimension(option, value):
    pres, certificate, _, _ = matrix_algebra_certificate()
    with pytest.raises(RuntimeError):
        verify_word_basis(pres, certificate, **{option: value})


def test_later_coordinate_and_product_operations_have_explicit_fresh_budgets():
    pres, certificate, x, y = matrix_algebra_certificate()
    result = verify_word_basis(pres, certificate)
    with pytest.raises(RuntimeError):
        result.coordinates(x+y, max_operations=1)
    with pytest.raises(RuntimeError):
        result.multiply(result.element([0, 1, 1, 0]), result.element([0, 1, 1, 0]),
                        max_operations=1)
    assert result.coordinates(x) == result.element([0, 1, 0, 0]).coordinates


@pytest.mark.parametrize('columns', [
    [], ([],), (((0,),),), (((-1, 1),),), (((True, 1),),), (((1, 1),),),
    (((0, 1), (0, 2)),),
])
def test_generator_action_shape_and_row_indices_are_checked(columns):
    field = RationalFunctionField()
    converted = columns
    if isinstance(columns, tuple):
        converted = tuple(tuple((pair[0], field.scalar(pair[1])) if len(pair) == 2 else pair
                                for pair in col) if isinstance(col, tuple) else col for col in columns)
    with pytest.raises((TypeError, ValueError)):
        GeneratorAction(field, converted)


def test_generator_action_coefficients_are_exact_and_parented_before_dropping_zeros():
    field = RationalFunctionField()
    with pytest.raises(TypeError):
        GeneratorAction('Q', ())
    with pytest.raises(TypeError):
        GeneratorAction(field, (((0, RationalFunctionField(('s',)).zero),),))
    with pytest.raises((TypeError, ValueError)):
        GeneratorAction(field, (((0, 0.5),),))
    with pytest.raises((TypeError, ValueError)):
        GeneratorAction(field, (((False, field.zero),),))


@pytest.mark.parametrize('coefficient,left,kind,index,right', [
    (1, (), 'input', 0, ()),
    (None, (), 'input', 0, ()),
])
def test_derivation_coefficients_are_explicit_field_elements(coefficient, left, kind, index, right):
    with pytest.raises(TypeError):
        DerivationTerm(coefficient, left, kind, index, right)


@pytest.mark.parametrize('left,kind,index,right', [
    ([], 'input', 0, ()), ((0,), 'input', 0, ()), ((True,), 'input', 0, ()),
    ((), 'basis', 0, ()), ((), 'input', -1, ()), ((), 'input', True, ()),
    ((), 'input', 0.5, ()), ((), 'input', 0, []), ((), 'input', 0, (-1,)),
])
def test_raw_derivation_contexts_and_references_are_exact(left, kind, index, right):
    with pytest.raises((TypeError, ValueError)):
        DerivationTerm(RationalFunctionField().one, left, kind, index, right)


@pytest.mark.parametrize('hash_value', [None, '0'*16, 'G'*64])
def test_derivation_graphs_need_full_lowercase_presentation_hashes(hash_value):
    with pytest.raises((TypeError, ValueError)):
        DerivationGraph(hash_value, (), ())


def test_derivation_graph_shape_and_coefficient_parents_are_validated():
    pres, certificate, _, _ = matrix_algebra_certificate()
    parent = pres.algebra
    foreign = FreeAlgebra(RationalFunctionField(('s',)), parent.generators)
    cases = [
        ([], ()), ((None,), ((),)), ((parent.zero,), []),
        ((parent.zero,), ()), ((parent.zero,), ([],)), ((parent.zero,), ((None,),)),
        ((parent.zero, foreign.zero), ((), ())),
        ((parent.zero,), ((term(foreign.field, 0),),)),
        ((parent.zero,), ((term(parent.field, 0, left=(3,)),),)),
    ]
    for nodes, rows in cases:
        with pytest.raises((TypeError, ValueError)):
            DerivationGraph(certificate.presentation_hash, nodes, rows)


def test_false_or_foreign_local_derivations_never_become_ideal_witnesses():
    pres, certificate, x, _ = matrix_algebra_certificate()
    field = pres.algebra.field
    false = DerivationGraph(pres.content_hash(), (x,), ((term(field, 0),),))
    with pytest.raises(ValueError):
        verify_derivations(pres, false)
    other = AlgebraPresentation(pres.algebra, (2*pres.relations[0],)+pres.relations[1:], pres.labels)
    with pytest.raises(ValueError):
        verify_derivations(other, certificate.proof)
    foreign = FreeAlgebra(RationalFunctionField(('s',)), pres.algebra.generators)
    foreign_graph = DerivationGraph(pres.content_hash(), (foreign.zero,), ((),))
    with pytest.raises(TypeError):
        verify_derivations(pres, foreign_graph)
    with pytest.raises(TypeError):
        verify_derivations(None, certificate.proof)
    with pytest.raises(TypeError):
        verify_derivations(pres, None)
    for option in ('max_operations', 'max_nodes'):
        with pytest.raises(ValueError):
            verify_derivations(pres, certificate.proof, **{option: 0})
    with pytest.raises(RuntimeError):
        verify_derivations(pres, certificate.proof, max_nodes=1)
    with pytest.raises(RuntimeError):
        verify_derivations(pres, certificate.proof, max_operations=1)


@pytest.mark.parametrize('kind', ['input', 'node'])
def test_zero_derivation_terms_cannot_hide_out_of_range_references(kind):
    pres, certificate, _, _ = matrix_algebra_certificate()
    graph = DerivationGraph(pres.content_hash(), (pres.algebra.zero,),
                            ((term(pres.algebra.field, 999, value=0, kind=kind),),))
    with pytest.raises(ValueError):
        verify_derivations(pres, graph)


@pytest.mark.parametrize('mutation', [
    lambda data: data.update(schema='other'),
    lambda data: data.update(extra=1),
    lambda data: data.pop('presentation_hash'),
    lambda data: data.update(nodes={}),
    lambda data: data.update(derivations={}),
    lambda data: data['derivations'].__setitem__(0, {}),
    lambda data: data['derivations'][1].__setitem__(0, None),
    lambda data: data['derivations'][1][0].update(left=(1,)),
    lambda data: data['derivations'][1][0].update(right=(1,)),
    lambda data: data['derivations'][1][0].update(index=True),
])
def test_derivation_graph_json_rejects_malformed_data(mutation):
    _, certificate, _, _ = matrix_algebra_certificate()
    data = deepcopy(certificate.proof.to_json())
    mutation(data)
    with pytest.raises((TypeError, ValueError)):
        DerivationGraph.from_json(data)


def test_defining_relations_must_vanish_on_every_column_not_only_the_cyclic_vector():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = (parent.generator(name) for name in parent.generators)
    pres = presentation(parent, (x, y*y))
    # The actual quotient has basis 1,y. Add the dependent word xy as a false
    # third basis word. Both input relations annihilate e0, and all proposed
    # words are cyclic, but x acts nontrivially on y's coordinate e1.
    mx = sparse_action(field, ((), ((2, 1),), ()))
    my = sparse_action(field, (((1, 1),), (), ()))
    graph = DerivationGraph(pres.content_hash(), (parent.zero, x, x*x*y, y*y, y*x*y), (
        (), (term(field, 0),), (term(field, 0, right=(1, 2)),),
        (term(field, 1),), (term(field, 0, left=(2,), right=(2,)),),
    ))
    false = WordBasisCertificate(pres.content_hash(), ((), (2,), (1, 2)),
                                 (mx, my), graph, ((1, 0, 2), (0, 3, 4)))
    with pytest.raises(ValueError):
        verify_word_basis(pres, false)


@pytest.mark.parametrize('mutation', [
    lambda data: data.update(schema='other'),
    lambda data: data.update(extra=1),
    lambda data: data.update(columns={}),
    lambda data: data['columns'].__setitem__(0, {}),
    lambda data: data['columns'][0].__setitem__(0, {}),
    lambda data: data['columns'][0].__setitem__(0, [1]),
    lambda data: data['columns'][0][0].__setitem__(0, True),
])
def test_generator_action_json_checks_entry_shapes_and_exact_row_indices(mutation):
    _, certificate, _, _ = matrix_algebra_certificate()
    data = deepcopy(certificate.actions[0].to_json())
    mutation(data)
    with pytest.raises((TypeError, ValueError)):
        GeneratorAction.from_json(data)


@pytest.mark.parametrize('hash_value', [None, '0'*16, 'G'*64])
def test_word_basis_certificate_hashes_are_exact_full_sha256(hash_value):
    _, certificate, _, _ = matrix_algebra_certificate()
    with pytest.raises((TypeError, ValueError)):
        WordBasisCertificate(hash_value, certificate.basis_words, certificate.actions,
                              certificate.proof, certificate.action_nodes)


def test_word_basis_certificate_rejects_unbound_proofs_and_malformed_action_references():
    pres, certificate, _, _ = matrix_algebra_certificate()
    wrong_graph = DerivationGraph('0'*64, certificate.proof.nodes, certificate.proof.derivations)
    cases = [
        (list(certificate.actions), certificate.proof, certificate.action_nodes),
        ((None, None), certificate.proof, certificate.action_nodes),
        (certificate.actions, None, certificate.action_nodes),
        (certificate.actions, wrong_graph, certificate.action_nodes),
        (certificate.actions, certificate.proof, list(certificate.action_nodes)),
        (certificate.actions, certificate.proof, ([0, 1, 0, 2], (0, 3, 4, 5))),
        (certificate.actions, certificate.proof, ((True, 1, 0, 2), (0, 3, 4, 5))),
        (certificate.actions, certificate.proof, ((-1, 1, 0, 2), (0, 3, 4, 5))),
        (certificate.actions, certificate.proof, ((0,), (0, 3, 4, 5))),
        (certificate.actions, certificate.proof, (certificate.action_nodes[0],)),
    ]
    for actions, graph, references in cases:
        with pytest.raises((TypeError, ValueError)):
            WordBasisCertificate(pres.content_hash(), certificate.basis_words, actions, graph, references)
    with pytest.raises(ValueError):
        WordBasisCertificate(pres.content_hash(), certificate.basis_words, certificate.actions,
                              certificate.proof, certificate.action_nodes, unit_node=0)


def test_verified_values_require_a_real_verifier_and_reject_unsupported_operators():
    from zipper_algebra.verification.finite import FiniteElement, VerifiedFiniteAlgebra
    with pytest.raises(TypeError):
        VerifiedFiniteAlgebra()
    with pytest.raises(TypeError):
        FiniteElement('unverified', ())
    pres, certificate, _, _ = matrix_algebra_certificate()
    result = verify_word_basis(pres, certificate)
    x = result.element([0, 1, 0, 0])
    other = verify_word_basis(pres, certificate)
    foreign_zero = other.zero
    with pytest.raises(TypeError):
        x+foreign_zero
    assert result.zero != foreign_zero
    assert result.one != 1
    for operation in (lambda: x+object(), lambda: object()-x, lambda: x-object(),
                      lambda: x*object(), lambda: object()*x):
        with pytest.raises(TypeError):
            operation()
    with pytest.raises(TypeError):
        verify_word_basis(None, certificate)
    with pytest.raises(TypeError):
        verify_word_basis(pres, None)


def test_action_fields_and_proof_node_bounds_are_checked_without_coercing_zero_data():
    pres, certificate, _, _ = matrix_algebra_certificate()
    foreign = RationalFunctionField((), 7)
    actions = tuple(GeneratorAction(foreign, tuple(tuple((row, foreign.scalar(value.numerator[0][1]))
                                                       for row, value in column)
                                                 for column in action.columns))
                    for action in certificate.actions)
    bad_field = WordBasisCertificate(certificate.presentation_hash, certificate.basis_words,
                                     actions, certificate.proof, certificate.action_nodes)
    with pytest.raises(TypeError):
        verify_word_basis(pres, bad_field)
    bad_reference = WordBasisCertificate(certificate.presentation_hash, certificate.basis_words,
                                         certificate.actions, certificate.proof,
                                         ((999, 1, 0, 2), certificate.action_nodes[1]))
    with pytest.raises(ValueError):
        verify_word_basis(pres, bad_reference)


def rational_quadratic_data(characteristic, parameters):
    """Nonmonomial denominators, with a separately scaled defining relation."""
    field = RationalFunctionField(parameters, characteristic)
    s = field.parameter('s')
    if len(parameters) == 1:
        c = 1/(s+1)
        scale = Fraction(5, 7)*(s+2)/(s*s+s+1)
    else:
        t = field.parameter('t')
        c = (s+t+1)/(s*t+s+1)
        scale = Fraction(5, 7)*(s+t+2)/(s*s+t*t+1)
    return field, c, scale


@pytest.mark.parametrize('characteristic', [0, 3])
@pytest.mark.parametrize('parameters', [('s',), ('s', 't')])
def test_nonmonomial_denominators_and_nonmonic_defining_coefficients(characteristic, parameters):
    pytest.importorskip('sympy')
    field, c, scale = rational_quadratic_data(characteristic, parameters)
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    monic = x*x-c
    pres = presentation(parent, (scale*monic, parent.zero))
    action = sparse_action(field, (((1, 1),), ((0, c),)))
    graph = DerivationGraph(pres.content_hash(), (parent.zero, monic),
                            ((), (term(field, 0, value=1/scale),)))
    cert = WordBasisCertificate(pres.content_hash(), ((), (1,)), (action,), graph, ((0, 1),))
    result = verify_word_basis(pres, cert)
    assert result.dimension == 2
    assert result.coordinates(x*x) == (c, field.zero)
    assert result.coordinates(x**5) == (field.zero, c*c)
    assert result.coordinates(scale*monic) == (field.zero, field.zero)
    finite_x = result.element((0, 1))
    assert finite_x*finite_x == result.element((c, 0))
    assert result.multiply(result.element((1, c)), finite_x) == result.element((c*c, 1))
    assert result.multiply(finite_x, result.element((1, c))) == result.element((c*c, 1))
    assert result.coordinates(parent.one) == (field.one, field.zero)
    assert result.coordinates(parent.zero) == (field.zero, field.zero)
    restored = WordBasisCertificate.from_json(json.loads(json.dumps(cert.to_json())))
    assert verify_word_basis(pres, restored).coordinates(x**4) == (c*c, field.zero)


@pytest.mark.parametrize('characteristic', [0, 3])
def test_large_exact_rational_function_coefficients_are_not_rounded(characteristic):
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s',), characteristic)
    s = field.parameter('s')
    c = ((2**80+7)*s+11)/(s*s+s+1)
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    relation = x*x-c
    pres = presentation(parent, (relation,))
    graph = DerivationGraph(pres.content_hash(), (parent.zero, relation),
                            ((), (term(field, 0),)))
    cert = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                (sparse_action(field, (((1, 1),), ((0, c),))),), graph, ((0, 1),))
    result = verify_word_basis(pres, cert)
    assert result.coordinates(x*x) == (c, field.zero)
    assert result.coordinates(x**4) == (c*c, field.zero)


@pytest.mark.parametrize('characteristic', [0, 3])
def test_rational_constant_relation_cannot_vanish_when_clearing_matrix_denominators(characteristic):
    pytest.importorskip('sympy')
    field, c, scale = rational_quadratic_data(characteristic, ('s', 't'))
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    pres = presentation(parent, (parent.scalar(scale),))
    # Honest spanning proofs exist because the true algebra is zero. The
    # proposed companion representation is cyclic but cannot satisfy scale*1=0.
    graph = DerivationGraph(pres.content_hash(), (parent.zero, x*x-c), (
        (), (term(field, 0, left=(1,), right=(1,), value=1/scale),
             term(field, 0, value=-c/scale)),
    ))
    cert = WordBasisCertificate(pres.content_hash(), ((), (1,)),
                                (sparse_action(field, (((1, 1),), ((0, c),))),), graph, ((0, 1),))
    with pytest.raises(ValueError, match='defining relation'):
        verify_word_basis(pres, cert)


@pytest.mark.parametrize('characteristic', [0, 3])
def test_denominator_clearing_does_not_replace_cyclic_independence_by_a_scaled_vector(characteristic):
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s',), characteristic)
    s = field.parameter('s')
    u = 1/(s+1)
    c = u*u
    scale = (s+2)/(s*s+s+1)
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    f = x*x-c
    pres = presentation(parent, (scale*f,))
    # B=(1,x,x²) is dependent in the actual quadratic quotient. The extra
    # coordinate is a one-dimensional module of eigenvalue u, so all relation
    # matrices and spanning identities are honest. x²*e0=c*e0 is not e2.
    action = sparse_action(field, (((1, 1),), ((0, c),),
                                   ((0, -u*c), (1, c), (2, u))))
    last_identity = (x-u)*f
    graph = DerivationGraph(pres.content_hash(), (parent.zero, f, last_identity), (
        (), (term(field, 0, value=1/scale),),
        (term(field, 0, left=(1,), value=1/scale), term(field, 0, value=-u/scale)),
    ))
    cert = WordBasisCertificate(pres.content_hash(), ((), (1,), (1, 1)),
                                (action,), graph, ((0, 1, 2),))
    with pytest.raises(ValueError, match='cyclic coordinate basis'):
        verify_word_basis(pres, cert)


@pytest.mark.parametrize('characteristic', [0, 3])
def test_parameter_fields_with_no_actions_or_zero_dimensional_seeds(characteristic):
    pytest.importorskip('sympy')
    field, c, scale = rational_quadratic_data(characteristic, ('s', 't'))
    scalar_parent = FreeAlgebra(field, ())
    scalar_pres = presentation(scalar_parent, (scalar_parent.zero,))
    scalar_graph = DerivationGraph(scalar_pres.content_hash(), (), ())
    scalar_cert = WordBasisCertificate(scalar_pres.content_hash(), ((),), (), scalar_graph, ())
    scalar_result = verify_word_basis(scalar_pres, scalar_cert)
    assert scalar_result.coordinates(scalar_parent.scalar(c)) == (c,)
    assert scalar_result.coordinates(scalar_parent.zero) == (field.zero,)
    assert scalar_result.multiply(scalar_result.element((c,)), scalar_result.one) == scalar_result.element((c,))
    zero_parent = FreeAlgebra(field, ('x',))
    zero_pres = presentation(zero_parent, (zero_parent.scalar(scale),))
    zero_graph = DerivationGraph(zero_pres.content_hash(), (zero_parent.one,),
                                 ((term(field, 0, value=1/scale),),))
    zero_cert = WordBasisCertificate(zero_pres.content_hash(), (),
                                     (sparse_action(field, ()),), zero_graph, ((),), unit_node=0)
    zero_result = verify_word_basis(zero_pres, zero_cert)
    assert zero_result.dimension == 0
    assert zero_result.coordinates(zero_parent.scalar(c)) == ()
    assert zero_result.one == zero_result.zero


def test_coordinate_iterables_are_rejected_after_at_most_dimension_plus_one_entries():
    pres, certificate, _, _ = matrix_algebra_certificate()
    result = verify_word_basis(pres, certificate)
    consumed = []

    def longer_than_dimension():
        for index in range(result.dimension+1):
            consumed.append(index)
            yield 0
        raise AssertionError('coordinate reader exceeded its dimension-plus-one bound')

    with pytest.raises(ValueError, match='coordinate count'):
        result.element(longer_than_dimension())
    assert consumed == list(range(result.dimension+1))


def test_zero_algebra_coordinate_iterables_are_bounded_too():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    pres = presentation(parent, (parent.one,))
    graph = DerivationGraph(pres.content_hash(), (parent.one,), ((term(field, 0),),))
    certificate = WordBasisCertificate(pres.content_hash(), (),
                                       (sparse_action(field, ()),), graph, ((),), unit_node=0)
    result = verify_word_basis(pres, certificate)
    consumed = []

    def nonempty():
        consumed.append(0)
        yield 0
        raise AssertionError('zero-algebra coordinates consumed more than one item')

    with pytest.raises(ValueError, match='coordinate count'):
        result.element(nonempty())
    assert consumed == [0]
    assert result.element(iter(())).coordinates == ()


@pytest.mark.parametrize('characteristic', [0, 3])
def test_two_noncommuting_actions_can_have_distinct_polynomial_denominators(characteristic):
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s', 't'), characteristic)
    s, t = (field.parameter(name) for name in field.parameters)
    a = 1/(s+1)
    b = (s+t+1)/(t*t+t+1)
    scales = ((s+2)/(s*s+1), (t+2)/(t*t+1), (s+t+2)/(s*t+1))
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = (parent.generator(name) for name in parent.generators)
    monic = (x*x-a, y*y-b, y*x+x*y)
    pres = presentation(parent, tuple(scale*relation for scale, relation in zip(scales, monic)))
    mx = sparse_action(field, (((1, 1),), ((0, a),), ((3, 1),), ((2, a),)))
    my = sparse_action(field, (((2, 1),), ((3, -1),), ((0, b),), ((1, -b),)))
    graph = DerivationGraph(pres.content_hash(),
                            (parent.zero, monic[0], monic[0]*y, monic[2], monic[1], y*x*y+b*x), (
        (), (term(field, 0, value=1/scales[0]),),
        (term(field, 0, right=(2,), value=1/scales[0]),),
        (term(field, 2, value=1/scales[2]),),
        (term(field, 1, value=1/scales[1]),),
        (term(field, 2, right=(2,), value=1/scales[2]),
         term(field, 1, left=(1,), value=-1/scales[1])),
    ))
    certificate = WordBasisCertificate(pres.content_hash(), ((), (1,), (2,), (1, 2)),
                                       (mx, my), graph, ((0, 1, 0, 2), (0, 3, 4, 5)))
    result = verify_word_basis(pres, certificate)
    zero = field.zero
    assert result.dimension == 4
    assert result.coordinates((x*y)**2) == (-a*b, zero, zero, zero)
    assert result.coordinates(y*x) == (zero, zero, zero, field.scalar(-1))
    assert result.coordinates(x*x*y) == (zero, zero, a, zero)
    assert result.multiply(result.element((0, 1, 0, 0)), result.element((0, 1, 0, 0))) == result.element((a, 0, 0, 0))
    assert result.multiply(result.element((0, 0, 1, 0)), result.element((0, 0, 1, 0))) == result.element((b, 0, 0, 0))
    u, v, w = 1/(s+1), 1/(t+1), 1/(s*t+1)
    left = result.element((0, 1, u, 0))
    right = result.element((0, v, w, 0))
    # (x+u*y)(v*x+w*y) = v*a+u*w*b+(w-u*v)*xy. The right
    # vector has two different nonmonomial denominators to clear together.
    assert result.multiply(left, right) == result.element((v*a+u*w*b, 0, 0, w-u*v))
    assert result.multiply(left, result.zero) == result.zero
