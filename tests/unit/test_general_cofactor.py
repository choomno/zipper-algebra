"""Independent identities for explicit cofactors over general coefficient fields."""
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
import json

import pytest

from zipper_algebra.algebra import AlgebraPresentation, CofactorIdentity, FreeAlgebra
from zipper_algebra.cofactor.general import (
    cofactor_representation, expand_derivation, verify_cofactor,
)
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.derivations import DerivationGraph, DerivationTerm, VerificationLimit
from zipper_algebra.verification.finite import GeneratorAction, WordBasisCertificate
from tests.unit.test_finite_algebra import matrix_algebra_certificate, quadratic_certificate


def presentation(parent, relations):
    return AlgebraPresentation(parent, tuple(relations), tuple(f'r{i}' for i in range(len(relations))))


def term(field, index, *, left=(), right=(), value=1, kind='input'):
    return DerivationTerm(field.scalar(value), left, kind, index, right)


def collected(identity):
    result = {}
    field = identity.target.algebra.field
    for left, index, right in identity.terms:
        assert len(left.terms) == len(right.terms) == 1
        (u,a), (v,b) = left.terms[0], right.terms[0]
        key = u,index,v
        result[key] = result.get(key, field.zero)+a*b
    return result


def simple_graph(field=None):
    parent = FreeAlgebra(field or RationalFunctionField(), ('x','y'))
    x,y = map(parent.generator,parent.generators)
    p = presentation(parent,(x*y-y*x,))
    # Forward shared reference. Its contexts cannot be commuted or reversed.
    target = x*x*y*y-x*y*x*y+2*y*x*y*x-2*y*y*x*x
    graph = DerivationGraph(p.content_hash(), (target,p.relations[0]), (
        (term(parent.field,1,left=(1,),right=(2,),kind='node'),
         term(parent.field,1,left=(2,),right=(1,),value=2,kind='node')),
        (term(parent.field,0),),
    ))
    return p,graph,x,y


def test_forward_shared_graph_expands_into_original_noncommuting_contexts():
    p,graph,x,y = simple_graph()
    cofactor = expand_derivation(p,graph,0)
    assert cofactor.target == x*x*y*y-x*y*x*y+2*y*x*y*x-2*y*y*x*x
    assert collected(cofactor) == {((1,),0,(2,)):p.algebra.field.one,
                                   ((2,),0,(1,)):p.algebra.field.scalar(2)}
    assert verify_cofactor(p,cofactor)
    assert cofactor.verify(p)
    encoded = CofactorIdentity.from_json(json.loads(json.dumps(cofactor.to_json())))
    assert verify_cofactor(p,encoded)


def test_repeated_contexts_collect_and_cancel_before_output():
    p,_,x,y = simple_graph()
    f=p.algebra.field
    graph=DerivationGraph(p.content_hash(), (3*p.relations[0],p.relations[0]), (
        (term(f,1,left=(1,),right=(2,),kind='node'),
         term(f,1,left=(1,),right=(2,),value=-1,kind='node'),
         term(f,1,value=3,kind='node')),
        (term(f,0),),
    ))
    result=expand_derivation(p,graph,0)
    assert collected(result)=={((),0,()):f.scalar(3)}
    assert result.target==3*(x*y-y*x)
    assert verify_cofactor(p,result)


@pytest.mark.parametrize('characteristic',[2,7])
def test_prime_field_coefficient_cancellation_and_zero_nodes(characteristic):
    f=RationalFunctionField((),characteristic)
    p,_,x,y=simple_graph(f)
    graph=DerivationGraph(p.content_hash(),(p.algebra.zero,),(
        (term(f,0,value=3),term(f,0,value=-3)),))
    result=expand_derivation(p,graph,0)
    assert result.target.is_zero() and result.terms==()
    assert verify_cofactor(p,result)


def test_generic_fraction_cancellation_preserves_original_field():
    pytest.importorskip('sympy')
    f=RationalFunctionField(('s',))
    p,_,_,_=simple_graph(f)
    s=f.parameter('s')
    value=(s*s-1)/(s-1)
    graph=DerivationGraph(p.content_hash(),(value*p.relations[0],),
                          ((DerivationTerm(value,(),'input',0,()),),))
    result=expand_derivation(p,graph,0)
    assert collected(result)=={((),0,()):s+1}
    assert result.target.algebra.field==f
    assert verify_cofactor(p,result)


def test_unselected_arithmetic_is_not_a_claim_about_the_returned_identity():
    p,_,x,_=simple_graph()
    graph=DerivationGraph(p.content_hash(),(p.relations[0],x),((term(p.algebra.field,0),),()))
    result=expand_derivation(p,graph,0)
    assert result.target==p.relations[0]
    assert verify_cofactor(p,result)
    with pytest.raises(ValueError):
        expand_derivation(p,graph,1)


@pytest.mark.parametrize('bad_kind',['out-of-range','zero-out-of-range','cycle'])
def test_every_raw_reference_is_validated_including_unselected_zero_terms(bad_kind):
    p,_,_,_=simple_graph();f=p.algebra.field
    last=term(f,99,value=0 if bad_kind=='zero-out-of-range' else 1)
    if bad_kind=='cycle': last=term(f,1,value=0,kind='node')
    graph=DerivationGraph(p.content_hash(),(p.relations[0],p.algebra.zero),
                          ((term(f,0),),(last,)))
    with pytest.raises(ValueError):
        expand_derivation(p,graph,0)


def test_selected_forged_polynomial_fails_final_exact_identity():
    p,graph,x,_=simple_graph()
    with pytest.raises(ValueError):
        expand_derivation(p,replace(graph,nodes=(graph.nodes[0]+x,graph.nodes[1])),0)


@pytest.mark.parametrize('bad_node',[-1,2,True,1.5,None])
def test_selected_node_must_be_an_existing_integer(bad_node):
    p,graph,_,_=simple_graph()
    with pytest.raises((ValueError,TypeError)):
        expand_derivation(p,graph,bad_node)


def test_graph_literal_hash_and_parent_are_bound_even_for_zero():
    p,graph,_,_=simple_graph()
    with pytest.raises(ValueError):
        expand_derivation(p,replace(graph,presentation_hash='0'*64),0)
    other=FreeAlgebra(RationalFunctionField(),('a','b'))
    foreign=DerivationGraph(p.content_hash(),(other.zero,), ((),))
    with pytest.raises(TypeError):
        expand_derivation(p,foreign,0)


@pytest.mark.parametrize('which',['presentation','graph'])
def test_expansion_rejects_wrong_public_types(which):
    p,graph,_,_=simple_graph()
    with pytest.raises(TypeError):
        expand_derivation(None if which=='presentation' else p,
                          None if which=='graph' else graph,0)


@pytest.mark.parametrize('option',['max_operations','max_contexts','max_nodes'])
@pytest.mark.parametrize('value',[0,True,-1])
def test_expansion_bounds_are_positive_exact_integers(option,value):
    p,graph,_,_=simple_graph()
    with pytest.raises(ValueError):
        expand_derivation(p,graph,0,**{option:value})


@pytest.mark.parametrize('options',[{'max_operations':1},{'max_contexts':1},{'max_nodes':1}])
def test_expansion_exhaustion_is_inconclusive(options):
    p,graph,_,_=simple_graph()
    with pytest.raises(VerificationLimit):
        expand_derivation(p,graph,0,**options)


def test_verifier_accepts_polynomial_factors_and_rejects_a_false_identity():
    A=FreeAlgebra(RationalFunctionField(),('x','y'))
    x,y=map(A.generator,A.generators)
    p=presentation(A,(x*y,))
    target=x*x*y*x-x*x*y*y+y*x*y*x-y*x*y*y
    identity=CofactorIdentity(p.content_hash(),target,((x+y,0,x-y),))
    assert verify_cofactor(p,identity)
    assert not verify_cofactor(p,replace(identity,target=target+1))


@pytest.mark.parametrize('zero_factor',[False,True])
def test_verifier_rejects_raw_unknown_relation_even_when_factor_is_zero(zero_factor):
    p,_,_,_=simple_graph();A=p.algebra
    identity=CofactorIdentity(p.content_hash(),A.zero,((A.zero if zero_factor else A.one,99,A.one),))
    with pytest.raises(ValueError): verify_cofactor(p,identity)


def test_verifier_rejects_wrong_hash_parent_public_types_and_bad_budget():
    p,graph,_,_=simple_graph();identity=expand_derivation(p,graph,1)
    with pytest.raises(ValueError): verify_cofactor(p,replace(identity,presentation_hash='0'*64))
    foreign=FreeAlgebra(RationalFunctionField(),('z',))
    bad=CofactorIdentity(p.content_hash(),foreign.zero,())
    with pytest.raises(TypeError): verify_cofactor(p,bad)
    for args in [(None,identity),(p,None)]:
        with pytest.raises(TypeError): verify_cofactor(*args)
    for value in [0,True,-1]:
        with pytest.raises(ValueError): verify_cofactor(p,identity,max_operations=value)
    with pytest.raises(VerificationLimit): verify_cofactor(p,identity,max_operations=1)


@pytest.mark.parametrize('characteristic',[0,2,7])
def test_matrix_algebra_elements_get_representatives_and_explicit_original_cofactors(characteristic):
    p,certificate,x,y=matrix_algebra_certificate(characteristic)
    result=cofactor_representation(p,certificate,y*x)
    assert result.remainder==1-x*y
    assert result.identity.target==y*x+x*y-1
    assert collected(result.identity)=={((),2,()):p.algebra.field.one}
    assert verify_cofactor(p,result.identity)
    assert result.identity.verify(p)
    reverse=cofactor_representation(p,certificate,x*y)
    assert reverse.remainder==x*y
    assert reverse.identity.target.is_zero()
    assert result.remainder!=reverse.remainder


def test_higher_noncommuting_words_use_both_closure_proofs_and_left_contexts():
    p,certificate,x,y=matrix_algebra_certificate()
    target=y*x*y*x+2*x*x*y+3*y*x+4
    result=cofactor_representation(p,certificate,target)
    assert result.remainder==8-4*x*y
    assert result.identity.target==target-(8-4*x*y)
    assert verify_cofactor(p,result.identity) and result.identity.verify(p)
    zero=cofactor_representation(p,certificate,x*x*y+y*y*x+x*y+y*x-1)
    assert zero.remainder.is_zero()
    assert zero.identity.target==x*x*y+y*y*x+x*y+y*x-1
    assert zero.identity.verify(p)
    with pytest.raises(FrozenInstanceError): result.remainder=p.algebra.zero


@pytest.mark.parametrize('characteristic',[0,3])
def test_parameter_field_reduction_produces_generic_cofactor(characteristic):
    pytest.importorskip('sympy')
    p,certificate,x,s=quadratic_certificate(characteristic)
    result=cofactor_representation(p,certificate,x**5)
    assert result.remainder==s*s*x
    assert result.identity.target==x**5-s*s*x
    assert verify_cofactor(p,result.identity)


def test_two_parameter_quadratic_cofactor_uses_the_declared_field():
    pytest.importorskip('sympy')
    f=RationalFunctionField(('a','b'));A=FreeAlgebra(f,('x',))
    a,b=map(f.parameter,f.parameters);x=A.generator('x');r=x*x-a*x-b
    p=presentation(A,(r,))
    action=GeneratorAction(f,(((1,f.one),),((0,b),(1,a))))
    graph=DerivationGraph(p.content_hash(),(A.zero,r),((),(term(f,0),)))
    certificate=WordBasisCertificate(p.content_hash(),((),(1,)),(action,),graph,((0,1),))
    result=cofactor_representation(p,certificate,x**3)
    assert result.remainder==a*b+(a*a+b)*x
    assert result.identity.target.algebra.field==f
    assert verify_cofactor(p,result.identity)


@pytest.mark.parametrize('target_kind',['zero','scalar'])
def test_zero_or_scalar_target_needs_no_nonzero_cofactor_terms(target_kind):
    p,certificate,_,_=matrix_algebra_certificate()
    target=p.algebra.zero if target_kind=='zero' else p.algebra.scalar(Fraction(2,3))
    result=cofactor_representation(p,certificate,target)
    assert result.remainder==target
    assert result.identity.target.is_zero() and result.identity.terms==()


@pytest.mark.parametrize('names',[(),('x','y')])
def test_zero_algebra_supplies_explicit_unit_ideal_cofactors(names):
    f=RationalFunctionField();A=FreeAlgebra(f,names)
    p=presentation(A,(A.scalar(3),))
    graph=DerivationGraph(p.content_hash(),(A.one,),((term(f,0,value=Fraction(1,3)),),))
    certificate=WordBasisCertificate(p.content_hash(),(),tuple(GeneratorAction(f,()) for _ in names),
                                     graph,tuple(() for _ in names),0)
    target=A.scalar(2) if not names else A.generator('x')*A.generator('y')+A.generator('y')+2
    result=cofactor_representation(p,certificate,target)
    assert result.remainder.is_zero() and result.identity.target==target
    assert verify_cofactor(p,result.identity) and result.identity.verify(p)


def test_field_without_generators_or_relations_has_zero_difference_identity():
    f=RationalFunctionField();A=FreeAlgebra(f,());p=presentation(A,())
    graph=DerivationGraph(p.content_hash(),(),())
    certificate=WordBasisCertificate(p.content_hash(),((),),(),graph,())
    result=cofactor_representation(p,certificate,A.scalar(Fraction(7,3)))
    assert result.remainder==A.scalar(Fraction(7,3))
    assert result.identity.terms==()


@pytest.mark.parametrize('which',['presentation','certificate','target'])
def test_general_reduction_rejects_wrong_public_types(which):
    p,c,x,_=matrix_algebra_certificate()
    args=[p,c,x]
    args[('presentation','certificate','target').index(which)]=None
    with pytest.raises(TypeError): cofactor_representation(*args)


def test_reduction_checks_literal_certificate_and_target_parent_before_zero_shortcut():
    p,c,_,_=matrix_algebra_certificate()
    with pytest.raises(ValueError):
        cofactor_representation(p,replace(c,presentation_hash='0'*64,
                                         proof=replace(c.proof,presentation_hash='0'*64)),p.algebra.zero)
    foreign=FreeAlgebra(RationalFunctionField(('z',)),p.algebra.generators)
    with pytest.raises(TypeError): cofactor_representation(p,c,foreign.zero)
    wrong=replace(c,actions=c.actions[::-1])
    with pytest.raises(ValueError): cofactor_representation(p,wrong,p.algebra.zero)


@pytest.mark.parametrize('option',['max_operations','max_contexts','max_nodes','max_dimension'])
@pytest.mark.parametrize('value',[0,True,-1])
def test_general_reduction_requires_positive_limits(option,value):
    p,c,x,_=matrix_algebra_certificate()
    with pytest.raises(ValueError): cofactor_representation(p,c,x,**{option:value})


@pytest.mark.parametrize('options',[{'max_operations':1},{'max_contexts':1},{'max_nodes':1},{'max_dimension':1}])
def test_general_reduction_bounds_are_inconclusive(options):
    p,c,x,y=matrix_algebra_certificate()
    with pytest.raises(VerificationLimit): cofactor_representation(p,c,y*x*y*x+x*x*y,**options)


def test_nonempty_cofactor_terms_can_cancel_to_the_zero_identity():
    p,_,x,y=simple_graph()
    identity=CofactorIdentity(p.content_hash(),p.algebra.zero,((x,0,y),(-x,0,y)))
    assert verify_cofactor(p,identity)


def test_scalar_nonzero_target_with_no_terms_is_not_an_identity():
    p,_,_,_=simple_graph()
    identity=CofactorIdentity(p.content_hash(),p.algebra.one,())
    assert not verify_cofactor(p,identity)


def test_valid_zero_coefficient_reference_does_not_create_a_context():
    p,_,_,_=simple_graph();f=p.algebra.field
    graph=DerivationGraph(p.content_hash(),(p.relations[0],),
                         ((term(f,0),term(f,0,value=0,left=(1,),right=(2,))),))
    result=expand_derivation(p,graph,0)
    assert collected(result)=={((),0,()):f.one}


def test_reduction_statistics_are_a_copy_and_count_the_explicit_output():
    p,c,x,y=matrix_algebra_certificate()
    result=cofactor_representation(p,c,y*x)
    assert result.statistics['basis_dimension']==4
    assert result.statistics['contexts']==len(result.identity.terms)
    statistics=result.statistics
    statistics['basis_dimension']=999
    assert result.statistics['basis_dimension']==4


def test_generated_proof_nodes_respect_the_bound_after_initial_verification():
    p,c,x,y=matrix_algebra_certificate()
    with pytest.raises(VerificationLimit,match='proof-node bound'):
        cofactor_representation(p,c,y*x,max_nodes=len(c.proof.nodes))


def test_certified_basis_need_not_be_prefix_closed_or_a_standard_word_basis():
    f=RationalFunctionField();A=FreeAlgebra(f,('x',));x=A.generator('x')
    relation=x*x-x-1
    p=presentation(A,(relation,))
    # The coordinates are B=(1,x^2), not the usual B=(1,x). Its generator
    # action and both closure proofs are literal independent input data.
    action=GeneratorAction(f,(((0,f.scalar(-1)),(1,f.one)),
                             ((0,f.scalar(-1)),(1,f.scalar(2)))))
    graph=DerivationGraph(p.content_hash(),(-relation,x**3-2*x*x+1),(
        (term(f,0,value=-1),),
        (term(f,0,left=(1,)),term(f,0,value=-1)),
    ))
    certificate=WordBasisCertificate(p.content_hash(),((),(1,1)),(action,),graph,((0,1),))
    result=cofactor_representation(p,certificate,x**5)
    assert result.remainder==5*x*x-2
    assert result.identity.target==x**5-5*x*x+2
    assert result.identity.target==(x**3+x*x+2*x-2)*relation
    assert verify_cofactor(p,result.identity) and result.identity.verify(p)
