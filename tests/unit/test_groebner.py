"""Independent finite-algebra controls and deliberately incomplete rewrite systems."""
from copy import deepcopy
from dataclasses import FrozenInstanceError
from fractions import Fraction
import json

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.groebner import (
    GroebnerCertificate, ProofTerm, VerificationLimit, verify_groebner,
)


def algebra(names=('x', 'y'), parameters=(), characteristic=0):
    return FreeAlgebra(RationalFunctionField(parameters, characteristic), names)


def presentation(parent, relations):
    return AlgebraPresentation(parent, tuple(relations),
                               tuple(f'r{i}' for i in range(len(relations))))


def input_certificate(pres, *, indices=None):
    indices = tuple(range(len(pres.relations))) if indices is None else tuple(indices)
    return GroebnerCertificate(pres.content_hash(), tuple(pres.relations[i] for i in indices),
                              tuple((ProofTerm(pres.algebra.field.one, (), 'input', i, ()),)
                                    for i in indices))


@pytest.mark.parametrize('characteristic', [0, 2, 7])
def test_exterior_algebra_has_four_independently_known_standard_words(characteristic):
    parent = algebra(characteristic=characteristic)
    x, y = (parent.generator(name) for name in parent.generators)
    pres = presentation(parent, [x*x, y*y, y*x+x*y])
    result = verify_groebner(pres, input_certificate(pres))
    assert result.dimension == 4
    assert result.standard_words == ((), (1,), (2,), (1, 2))
    assert result.normal_form(y*x) == -x*y
    assert result.normal_form((x+y)**2).is_zero()
    assert result.normal_form(y*x*y).is_zero()
    assert result.normal_form(x*y) == x*y
    assert result.report['standard_words_complete'] is True


def test_commutative_truncated_polynomial_algebra_and_nonnormal_input():
    parent = algebra()
    x, y = (parent.generator(name) for name in parent.generators)
    pres = presentation(parent, [x*x, y*y, y*x-x*y])
    result = verify_groebner(pres, input_certificate(pres))
    assert result.dimension == 4
    assert result.normal_form((x+y)**2) == 2*x*y
    assert result.normal_form(y*x+x*y) == 2*x*y


@pytest.mark.parametrize('characteristic', [0, 3])
def test_symbolic_quadratic_relation_over_an_explicit_parameter_field(characteristic):
    parent = algebra(('x',), ('s',), characteristic)
    x, s = parent.generator('x'), parent.field.parameter('s')
    pres = presentation(parent, [x*x-s])
    result = verify_groebner(pres, input_certificate(pres))
    assert result.dimension == 2
    assert result.standard_words == ((), (1,))
    assert result.normal_form(x**5) == s*s*x
    assert result.normal_form((x*x-s)*(x+1)).is_zero()


def test_nonmonic_rules_are_normalized_without_changing_the_original_ideal():
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [Fraction(2, 3)*x*x-2*x])
    result = verify_groebner(pres, input_certificate(pres))
    assert result.dimension == 2
    assert result.normal_form(x*x) == 3*x
    assert result.normal_form(x**3) == 9*x


@pytest.mark.parametrize('case', ['self', 'inclusion', 'equal-leading', 'distinct-overlap'])
def test_every_kind_of_unresolved_ambiguity_prevents_completion(case):
    parent = algebra(('a', 'b', 'c'))
    a, b, c = (parent.generator(name) for name in parent.generators)
    relations = {
        'self': [a*b*a-b],
        'inclusion': [a*b*a-a, b-1],
        'equal-leading': [a*b-a, a*b-b],
        'distinct-overlap': [a*b-a, b*c-b],
    }[case]
    pres = presentation(parent, relations)
    # Each candidate has an honest trivial input proof and reduces its own
    # defining relation to zero. These failures concern completeness itself.
    with pytest.raises(ValueError):
        verify_groebner(pres, input_certificate(pres), prune=False)


def test_pruning_rechecks_reverse_inclusion_and_preserves_a_redundant_system():
    parent = algebra()
    x, y = (parent.generator(name) for name in parent.generators)
    bad = presentation(parent, [x*y-x, x*y-y])
    with pytest.raises(ValueError):
        verify_groebner(bad, input_certificate(bad), prune=True)
    good = presentation(parent, [x*y-x, y-1])
    for prune in (False, True):
        result = verify_groebner(good, input_certificate(good), max_words=5, prune=prune)
        assert result.normal_form(y) == parent.one
        assert result.normal_form(x*y) == x
        assert result.dimension is None


def test_basis_proofs_alone_do_not_establish_reverse_ideal_inclusion():
    parent = algebra()
    x, y = (parent.generator(name) for name in parent.generators)
    pres = presentation(parent, [x*x, y*y])
    with pytest.raises(ValueError):
        verify_groebner(pres, input_certificate(pres, indices=(0,)))


def test_a_forged_local_identity_cannot_enlarge_the_input_ideal():
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [x*x])
    forged = GroebnerCertificate(pres.content_hash(), (x,),
                                ((ProofTerm(parent.field.one, (), 'input', 0, ()),),))
    with pytest.raises(ValueError):
        verify_groebner(pres, forged)


def test_forward_derivation_references_and_nonmonic_duplicate_rules_are_valid():
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [x*x])
    certificate = GroebnerCertificate(pres.content_hash(), (x*x, 2*x*x), (
        (ProofTerm(parent.field.scalar(Fraction(1, 2)), (), 'basis', 1, ()),),
        (ProofTerm(parent.field.scalar(2), (), 'input', 0, ()),),
    ))
    for prune in (False, True):
        result = verify_groebner(pres, certificate, prune=prune)
        assert result.dimension == 2
        assert result.normal_form(x*x).is_zero()


@pytest.mark.parametrize('zero_cycle', [False, True])
def test_cyclic_derivations_are_rejected_even_if_local_equalities_hold(zero_cycle):
    parent = algebra(('x',))
    x = parent.generator('x')
    if zero_cycle:
        pres = presentation(parent, [x])
        certificate = GroebnerCertificate(pres.content_hash(), (x,), ((
            ProofTerm(parent.field.one, (), 'input', 0, ()),
            ProofTerm(parent.field.zero, (), 'basis', 0, ()),
        ),))
    else:
        pres = presentation(parent, [])
        certificate = GroebnerCertificate(pres.content_hash(), (x, x), (
            (ProofTerm(parent.field.one, (), 'basis', 1, ()),),
            (ProofTerm(parent.field.one, (), 'basis', 0, ()),),
        ))
    with pytest.raises(ValueError):
        verify_groebner(pres, certificate)


def test_derivation_contexts_preserve_the_order_of_noncommuting_words():
    parent = algebra()
    x, y = (parent.generator(name) for name in parent.generators)
    relation = x*y-x
    pres = presentation(parent, [relation])
    proof = ((ProofTerm(parent.field.one, (), 'input', 0, ()),),
             (ProofTerm(parent.field.one, (2,), 'input', 0, (1,)),))
    good = GroebnerCertificate(pres.content_hash(), (relation, y*relation*x), proof)
    assert verify_groebner(pres, good, prune=True, max_words=6).normal_form(y*relation*x).is_zero()
    bad = GroebnerCertificate(pres.content_hash(), (relation, x*relation*y), proof)
    with pytest.raises(ValueError):
        verify_groebner(pres, bad, prune=True, max_words=6)


def test_zero_algebra_and_scalar_only_free_algebra_are_distinguished():
    parent = algebra(('x',))
    x = parent.generator('x')
    unit_ideal = presentation(parent, [parent.scalar(3)])
    zero = verify_groebner(unit_ideal, input_certificate(unit_ideal))
    assert zero.dimension == 0
    assert zero.standard_words == ()
    assert zero.normal_form(x+1).is_zero()
    scalar_parent = algebra(())
    scalar_pres = presentation(scalar_parent, [])
    scalars = verify_groebner(scalar_pres, input_certificate(scalar_pres))
    assert scalars.dimension == 1
    assert scalars.standard_words == ((),)
    assert scalars.normal_form(scalar_parent.scalar(7)) == scalar_parent.scalar(7)


def test_empty_basis_and_finite_enumeration_budget_do_not_assert_infinite_dimension():
    parent = algebra(('x',))
    x = parent.generator('x')
    free = presentation(parent, [])
    result = verify_groebner(free, input_certificate(free), max_words=4)
    assert result.dimension is None
    assert result.standard_words is None
    assert result.report['standard_words_complete'] is False
    assert result.normal_form(x**20) == x**20
    finite = presentation(parent, [x*x])
    limited = verify_groebner(finite, input_certificate(finite), max_words=1)
    assert limited.dimension is None
    assert limited.normal_form(x*x).is_zero()


def test_zero_defining_relations_and_wrong_parent_zero_elements():
    parent = algebra(('x',))
    pres = presentation(parent, [parent.zero])
    cert = GroebnerCertificate(pres.content_hash(), (), ())
    result = verify_groebner(pres, cert, max_words=3)
    assert result.normal_form(parent.zero) == parent.zero
    with pytest.raises(TypeError):
        result.normal_form(algebra(('x',), ('s',)).zero)
    with pytest.raises((TypeError, ValueError)):
        result.normal_form(0)


def test_budget_exhaustion_produces_no_verified_reduction():
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [x*x, x**8])
    with pytest.raises(VerificationLimit):
        verify_groebner(pres, input_certificate(pres, indices=(0,)), max_reductions=1)


@pytest.mark.parametrize('option', ['max_words', 'max_reductions'])
@pytest.mark.parametrize('bad', [0, -1, True, 1.5, None])
def test_verification_budgets_are_positive_exact_integers(option, bad):
    parent = algebra(())
    pres = presentation(parent, [])
    with pytest.raises((TypeError, ValueError)):
        verify_groebner(pres, input_certificate(pres), **{option: bad})


def test_certificate_is_bound_to_the_literal_presentation_and_full_parent():
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [x*x])
    cert = input_certificate(pres)
    other = presentation(parent, [2*x*x])
    with pytest.raises(ValueError):
        verify_groebner(other, cert)
    foreign = algebra(('x',), ('s',))
    foreign_x = foreign.generator('x')
    false_parent = GroebnerCertificate(pres.content_hash(), (foreign_x*foreign_x,),
                                      ((ProofTerm(foreign.field.one, (), 'input', 0, ()),),))
    with pytest.raises((TypeError, ValueError)):
        verify_groebner(pres, false_parent)


def test_certificate_json_round_trip_keeps_parameter_coefficients_and_derivations():
    parent = algebra(('x',), ('s',))
    x, s = parent.generator('x'), parent.field.parameter('s')
    pres = presentation(parent, [x*x-s])
    cert = input_certificate(pres)
    restored = GroebnerCertificate.from_json(json.loads(json.dumps(cert.to_json())))
    assert restored.to_json() == cert.to_json()
    assert verify_groebner(pres, restored).normal_form(x**4) == parent.scalar(s*s)
    with pytest.raises(FrozenInstanceError):
        restored.basis = ()
    with pytest.raises(FrozenInstanceError):
        restored.derivations[0][0].index = 99


@pytest.mark.parametrize('mutation', [
    lambda data: data.update(schema='other'),
    lambda data: data.update(extra='ignored?'),
    lambda data: data.pop('presentation_hash'),
    lambda data: data.update(basis={}),
    lambda data: data.update(derivations={}),
    lambda data: data['derivations'].__setitem__(0, {}),
    lambda data: data['derivations'][0].__setitem__(0, []),
    lambda data: data['derivations'][0][0].update(extra=1),
    lambda data: data['derivations'][0][0].update(kind='other'),
    lambda data: data['derivations'][0][0].update(index=True),
    lambda data: data['derivations'][0][0].update(left=[0]),
])
def test_malformed_certificate_json_is_rejected(mutation):
    parent = algebra(('x',))
    pres = presentation(parent, [parent.generator('x')**2])
    data = deepcopy(input_certificate(pres).to_json())
    mutation(data)
    with pytest.raises((TypeError, ValueError)):
        GroebnerCertificate.from_json(data)


@pytest.mark.parametrize('coefficient,left,kind,index,right', [
    (1, (), 'input', 0, ()),
    (None, (), 'input', 0, ()),
    ('one', (), 'input', 0, ()),
])
def test_proof_terms_require_declared_field_elements(coefficient, left, kind, index, right):
    with pytest.raises(TypeError):
        ProofTerm(coefficient, left, kind, index, right)


@pytest.mark.parametrize('left,kind,index,right', [
    ([], 'input', 0, ()), ((0,), 'input', 0, ()), ((True,), 'input', 0, ()),
    ((), 'other', 0, ()), ((), 'input', -1, ()), ((), 'input', True, ()),
    ((), 'input', 0.5, ()), ((), 'input', 0, []), ((), 'input', 0, (-1,)),
])
def test_proof_terms_validate_raw_contexts_kinds_and_reference_indices(left, kind, index, right):
    with pytest.raises((TypeError, ValueError)):
        ProofTerm(RationalFunctionField().one, left, kind, index, right)


@pytest.mark.parametrize('reference_kind', ['input', 'basis'])
def test_zero_coefficients_do_not_hide_invalid_references(reference_kind):
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [x*x])
    bad = GroebnerCertificate(pres.content_hash(), (x*x,), ((
        ProofTerm(parent.field.one, (), 'input', 0, ()),
        ProofTerm(parent.field.zero, (), reference_kind, 999, ()),
    ),))
    with pytest.raises(ValueError):
        verify_groebner(pres, bad)


@pytest.mark.parametrize('hash_value', [None, 1, '0'*16, 'A'*64, 'g'*64])
def test_certificate_hashes_are_full_lowercase_sha256(hash_value):
    with pytest.raises((TypeError, ValueError)):
        GroebnerCertificate(hash_value, (), ())


def test_certificate_container_parent_and_derivation_validation():
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [x*x])
    one = ProofTerm(parent.field.one, (), 'input', 0, ())
    foreign = algebra(('x',), ('s',))
    fx = foreign.generator('x')
    cases = [
        ([x*x], ((one,),)),
        ((parent.zero,), ((one,),)),
        ((x*x,), []),
        ((x*x,), ()),
        ((x*x,), ([one],)),
        ((x*x,), ((None,),)),
        ((x*x, fx), ((one,), (one,))),
        ((x*x,), ((ProofTerm(foreign.field.one, (), 'input', 0, ()),),)),
        ((x*x,), ((ProofTerm(parent.field.one, (2,), 'input', 0, ()),),)),
        ((x*x,), ((ProofTerm(parent.field.one, (), 'input', 0, (2,)),),)),
    ]
    for basis, proofs in cases:
        with pytest.raises((TypeError, ValueError)):
            GroebnerCertificate(pres.content_hash(), basis, proofs)


@pytest.mark.parametrize('value', [None, [], 'certificate'])
def test_certificate_json_requires_an_object(value):
    with pytest.raises((TypeError, ValueError)):
        GroebnerCertificate.from_json(value)


def test_budget_limit_on_a_later_reduction_does_not_become_a_false_normal_form():
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [x*x-x])
    result = verify_groebner(pres, input_certificate(pres), max_reductions=10)
    with pytest.raises(VerificationLimit):
        result.normal_form(x**50)
    # Previously checked irreducible words remain exact after the interruption.
    assert result.normal_form(x) == x


def test_matrix_algebra_control_has_noncommuting_idempotent_products():
    parent = algebra()
    x, y = (parent.generator(name) for name in parent.generators)
    # x=E12,y=E21 in M2(Q): x²=y²=0 and xy+yx=1.
    pres = presentation(parent, [x*x, y*y, y*x+x*y-1])
    result = verify_groebner(pres, input_certificate(pres))
    assert result.dimension == 4
    assert result.normal_form(x*y) == x*y
    assert result.normal_form(y*x) == parent.one-x*y
    assert result.normal_form(x*y) != result.normal_form(y*x)
    assert result.normal_form(x*y*x) == x
    assert result.normal_form((x*y)**2) == x*y
    assert result.normal_form(x*y*y*x).is_zero()


def test_commutative_polynomial_ring_has_exact_normal_forms_without_dimension_claim():
    parent = algebra()
    x, y = (parent.generator(name) for name in parent.generators)
    pres = presentation(parent, [y*x-x*y])
    result = verify_groebner(pres, input_certificate(pres), max_words=12)
    assert result.dimension is None
    assert result.normal_form(y*y*x*x) == x*x*y*y
    assert result.normal_form((y*x)*(y*x)) == x*x*y*y
    assert result.normal_form(x+y) == x+y


def test_derived_redundant_rule_can_have_a_nonzero_inclusion_ambiguity():
    parent = algebra(('a', 'b', 'c'))
    a, b, c = (parent.generator(name) for name in parent.generators)
    r = b*c*b-1
    # The self-overlap r*cb-bc*r equals bc-cb. Add a derived relation
    # whose leading word abb contains the input leading word a. The defining
    # inputs still reduce to zero, but the inclusion ambiguity is bc-cb.
    derived = a*b*b+b*c-c*b
    pres = presentation(parent, [a, r])
    cert = GroebnerCertificate(pres.content_hash(), (a, derived, r), (
        (ProofTerm(parent.field.one, (), 'input', 0, ()),),
        (ProofTerm(parent.field.one, (), 'input', 0, (2, 2)),
         ProofTerm(parent.field.one, (), 'input', 1, (3, 2)),
         ProofTerm(parent.field.scalar(-1), (2, 3), 'input', 1, ())),
        (ProofTerm(parent.field.one, (), 'input', 1, ()),),
    ))
    with pytest.raises(ValueError, match='inclusion ambiguity'):
        verify_groebner(pres, cert)


def test_verification_result_cannot_be_instantiated_or_relabelled_by_callers():
    from zipper_algebra.verification.groebner import VerifiedReduction
    with pytest.raises(TypeError):
        VerifiedReduction()
    parent = algebra(())
    pres = presentation(parent, [])
    result = verify_groebner(pres, input_certificate(pres))
    with pytest.raises(FrozenInstanceError):
        result.standard_words = ()
    report = result.report
    report['dimension'] = 9000
    assert result.report['dimension'] == 1
    with pytest.raises(TypeError):
        verify_groebner(None, input_certificate(pres))
    with pytest.raises(TypeError):
        verify_groebner(pres, None)


def test_branching_rewrites_share_normal_forms_of_repeated_descendants():
    parent = algebra(('x',))
    x = parent.generator('x')
    pres = presentation(parent, [x*x-x-1])
    result = verify_groebner(pres, input_certificate(pres))
    # Fibonacci recurrence x^n=F_n*x+F_{n-1} in Q[x]/(x²-x-1).
    assert result.normal_form(x**8) == 21*x+13
    assert result.normal_form(x**10) == 55*x+34
