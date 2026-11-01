"""Independent finite-algebra controls for shared reduction-proof production."""
from dataclasses import replace
from fractions import Fraction

import pytest

from zipper_algebra.action_production import word_basis_from_groebner
from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.finite import WordBasisCertificate, verify_word_basis
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm, VerificationLimit


def inputs(algebra, relations):
    pres = AlgebraPresentation(algebra, tuple(relations), tuple(f'r{i}' for i in range(len(relations))))
    cert = GroebnerCertificate(pres.content_hash(), pres.relations,
                               tuple((ProofTerm(algebra.field.one, (), 'input', i, ()),)
                                     for i in range(len(relations))))
    return pres, cert


def matrix(action):
    field = action.field
    columns = [dict(column) for column in action.columns]
    return tuple(tuple(column.get(row, field.zero) for column in columns)
                 for row in range(action.dimension))


@pytest.mark.parametrize('characteristic', [0, 2, 7])
def test_exterior_actions_have_expected_noncommuting_orientation(characteristic):
    A = FreeAlgebra(RationalFunctionField((), characteristic), ('x', 'y'))
    x, y = map(A.generator, A.generators)
    p, g = inputs(A, [x*x, y*y, y*x+x*y])
    c = word_basis_from_groebner(p, g)
    assert c.basis_words == ((), (1,), (2,), (1, 2))
    assert matrix(c.actions[0]) == tuple(tuple(A.field.scalar(v) for v in row) for row in (
        (0,0,0,0), (1,0,0,0), (0,0,0,0), (0,0,1,0)))
    assert matrix(c.actions[1]) == tuple(tuple(A.field.scalar(v) for v in row) for row in (
        (0,0,0,0), (0,0,0,0), (1,0,0,0), (0,-1,0,0)))
    assert verify_word_basis(p, WordBasisCertificate.from_json(c.to_json())).dimension == 4
    assert c.proof.nodes[:3] == g.basis
    assert len(c.proof.nodes) < 3 + 1 + 2 * len(c.basis_words)


def test_generic_clifford_family_recovers_independently_specified_matrices():
    pytest.importorskip('sympy')
    f = RationalFunctionField(('alpha', 'beta', 'gamma'))
    A = FreeAlgebra(f, ('x', 'y'))
    x, y = map(A.generator, A.generators)
    a, b, c = map(f.parameter, f.parameters)
    p, g = inputs(A, [x*x-a, y*y-b, y*x+x*y-c])
    result = word_basis_from_groebner(p, g)
    zero, one = f.zero, f.one
    assert matrix(result.actions[0]) == ((zero,a,zero,zero), (one,zero,zero,zero),
                                          (zero,zero,zero,a), (zero,zero,one,zero))
    assert matrix(result.actions[1]) == ((zero,c,b,zero), (zero,zero,zero,-b),
                                          (one,zero,zero,c), (zero,-one,zero,zero))
    finite = verify_word_basis(p, result)
    xv, yv = finite.element((0,1,0,0)), finite.element((0,0,1,0))
    assert xv*yv != yv*xv
    assert xv*xv == finite.scalar(a)
    assert yv*xv + xv*yv == finite.scalar(c)
    # yxy needs both a contextual input and a nontrivial child reduction.
    yxy = next(node for node in result.proof.nodes if dict(node.terms).get((2,1,2)) == one)
    assert yxy == y*x*y + b*x-c*y


def test_characteristic_two_tail_cancellation_is_exact():
    A = FreeAlgebra(RationalFunctionField((), 2), ('x', 'y'))
    x, y = map(A.generator, A.generators)
    p, g = inputs(A, [x*x, y*y-y, y*x-x*y-x])
    result = word_basis_from_groebner(p, g)
    assert result.basis_words == ((), (1,), (2,), (1,2))
    assert result.actions[1].columns[-1] == ()  # yxy=xyy+xy=0
    assert verify_word_basis(p, result).dimension == 4


def test_nonmonic_fractional_relation_uses_inverse_leading_coefficient():
    A = FreeAlgebra(RationalFunctionField(), ('x',))
    x = A.generator('x')
    p, g = inputs(A, [Fraction(2,3)*x*x-2*x])
    result = word_basis_from_groebner(p, g)
    assert result.actions[0].columns == (((1,A.field.one),), ((1,A.field.scalar(3)),))
    assert result.proof.derivations[-1][0].coefficient == A.field.scalar(Fraction(3,2))


def test_forward_and_shared_original_proof_references_are_retained():
    A = FreeAlgebra(RationalFunctionField(), ('x',))
    x = A.generator('x')
    p, _ = inputs(A, [x*x])
    f = A.field
    g = GroebnerCertificate(p.content_hash(), (x*x, 2*x*x), (
        (ProofTerm(f.scalar(Fraction(1,2)), (), 'basis', 1, ()),),
        (ProofTerm(f.scalar(2), (), 'input', 0, ()),)))
    result = word_basis_from_groebner(p, g)
    assert result.proof.derivations[0][0].kind == 'node'
    assert result.proof.derivations[0][0].index == 1
    assert verify_word_basis(p, result).dimension == 2


@pytest.mark.parametrize('names', [(), ('x', 'y')])
def test_unit_ideal_has_explicit_unit_proof_and_empty_actions(names):
    A = FreeAlgebra(RationalFunctionField(), names)
    p, g = inputs(A, [A.scalar(3)])
    result = word_basis_from_groebner(p, g)
    assert result.basis_words == ()
    assert len(result.actions) == len(names)
    assert all(action.columns == () for action in result.actions)
    assert result.proof.nodes[result.unit_node] == A.one
    assert result.proof.derivations[result.unit_node][0].coefficient == A.field.scalar(Fraction(1,3))
    assert verify_word_basis(p, result).dimension == 0


def test_ground_field_with_no_generators_or_relations_is_one_dimensional():
    A = FreeAlgebra(RationalFunctionField(), ())
    p, g = inputs(A, [])
    result = word_basis_from_groebner(p, g)
    assert result.basis_words == ((),)
    assert result.actions == () and result.action_nodes == ()
    assert verify_word_basis(p, result).dimension == 1


def test_unbounded_free_algebra_does_not_get_a_finite_action_certificate():
    A = FreeAlgebra(RationalFunctionField(), ('x',))
    p, g = inputs(A, [])
    with pytest.raises(VerificationLimit, match='enumeration'):
        word_basis_from_groebner(p, g, max_words=5)


def test_finite_enumeration_budget_is_inconclusive():
    A = FreeAlgebra(RationalFunctionField(), ('x',))
    x = A.generator('x')
    p, g = inputs(A, [x**4])
    with pytest.raises(VerificationLimit, match='enumeration'):
        word_basis_from_groebner(p, g, max_words=3)


@pytest.mark.parametrize('option,value', [('max_nodes',0), ('max_nodes',True),
                                         ('max_operations',0), ('max_words',0), ('max_reductions',0)])
def test_invalid_bounds_are_rejected(option, value):
    A = FreeAlgebra(RationalFunctionField(), ('x',))
    p,g = inputs(A, [A.generator('x')**2])
    with pytest.raises(ValueError):
        word_basis_from_groebner(p,g, **{option:value})


@pytest.mark.parametrize('limit', [1,2])
def test_initial_and_generated_proof_node_bounds(limit):
    A = FreeAlgebra(RationalFunctionField(), ('x',))
    x = A.generator('x')
    p,g = inputs(A, [x*x, 2*x*x])
    with pytest.raises(VerificationLimit, match='node bound'):
        word_basis_from_groebner(p,g,max_nodes=limit)


def test_construction_operation_exhaustion_is_inconclusive():
    A = FreeAlgebra(RationalFunctionField(), ('x',))
    p,g = inputs(A, [A.generator('x')**2])
    with pytest.raises(VerificationLimit):
        word_basis_from_groebner(p,g,max_operations=1)


def test_incomplete_or_forged_groebner_input_fails_before_action_production():
    A = FreeAlgebra(RationalFunctionField(), ('x','y'))
    x,y = map(A.generator,A.generators)
    p,g = inputs(A, [x*y-x,y*x-y])
    with pytest.raises(ValueError,match='ambiguity'):
        word_basis_from_groebner(p,g)
    p,g = inputs(A, [x*x,y*y,y*x+x*y])
    forged = replace(g,basis=(x,g.basis[1],g.basis[2]))
    with pytest.raises(ValueError,match='identity'):
        word_basis_from_groebner(p,forged)


def test_generic_three_strand_hecke_algebra_has_six_regular_basis_words():
    pytest.importorskip('sympy')
    f = RationalFunctionField(('q',))
    A = FreeAlgebra(f, ('x', 'y'))
    x,y = map(A.generator,A.generators)
    q = f.parameter('q')
    p,g = inputs(A, [x*x-(1-q)*x-q, y*y-(1-q)*y-q, y*x*y-x*y*x])
    result = word_basis_from_groebner(p,g)
    assert result.basis_words == ((),(1,),(2,),(1,2),(2,1),(1,2,1))
    finite = verify_word_basis(p,result)
    xv,yv = finite.element((0,1,0,0,0,0)), finite.element((0,0,1,0,0,0))
    assert xv*yv*xv == yv*xv*yv
    assert xv*yv != yv*xv
    assert (xv-finite.one)*(xv+finite.scalar(q)) == finite.zero
    assert finite.dimension == 6
