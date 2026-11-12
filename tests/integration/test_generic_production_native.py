"""Small installed-native generic proof reconstruction with literal sample proofs."""
from fractions import Fraction

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.production import SupportFailure, reconstruct_groebner
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm, verify_groebner


@pytest.fixture
def native_finiteflow():
    pytest.importorskip('sympy', reason='generic reconstruction requires symbolic arithmetic')
    return pytest.importorskip('fflow', reason='native reconstruction is optional')


def trivial_sample(pres, point):
    numeric = pres.specialize(point)
    return GroebnerCertificate(numeric.content_hash(), numeric.relations,
        tuple((ProofTerm(numeric.algebra.field.one, (), 'input', index, ()),)
              for index in range(len(numeric.relations))))


def test_native_generic_noncommutative_four_dimensional_control(native_finiteflow):
    field = RationalFunctionField(('alpha', 'beta', 'gamma'))
    a, b, c = (field.parameter(name) for name in field.parameters)
    parent = FreeAlgebra(field, ('u', 'v'))
    u, v = parent.generator('u'), parent.generator('v')
    pres = AlgebraPresentation(parent, (u*u-a, v*v-b, v*u+u*v-c), ('uu', 'vv', 'vu'))
    point = {'alpha': Fraction(2, 3), 'beta': 5, 'gamma': 7}
    result = reconstruct_groebner(pres, trivial_sample(pres, point), point=point,
                                 max_primes=3, max_degree=4, threads=2,
                                 start_prime=1, shift=(3, 5, 7))
    assert result.reduction.standard_words == ((), (1,), (2,), (1, 2))
    assert result.reduction.normal_form(v*u) == -u*v+c
    assert result.reduction.normal_form(u*u*v) == a*v
    assert verify_groebner(pres, GroebnerCertificate.from_json(result.certificate.to_json())).dimension == 4


def test_native_forward_reference_reconstructs_a_shared_generic_proof(native_finiteflow):
    field = RationalFunctionField(('coefficient', 'constant'))
    coefficient, constant = field.parameter('coefficient'), field.parameter('constant')
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    pres = AlgebraPresentation(parent, (coefficient*x*x-constant,), ('quadratic',))
    point = {'coefficient': 3, 'constant': 2}
    numeric = pres.specialize(point)
    relation = numeric.relations[0]
    one = numeric.algebra.field.one
    sample = GroebnerCertificate(numeric.content_hash(), (2*relation, relation), (
        (ProofTerm(numeric.algebra.field.scalar(2), (), 'basis', 1, ()),),
        (ProofTerm(one, (), 'input', 0, ()),),
    ))
    result = reconstruct_groebner(pres, sample, point=point, max_primes=3, max_degree=4)
    assert result.reduction.normal_form(x*x) == parent.scalar(constant/coefficient)
    assert any(term.kind == 'basis' and term.index == 1
               for term in result.certificate.derivations[0])


def test_native_structural_tail_recovery_uses_generic_contexts(native_finiteflow):
    field = RationalFunctionField(('a',))
    parent = FreeAlgebra(field, ('x',))
    x, a = parent.generator('x'), field.parameter('a')
    pres = AlgebraPresentation(parent, (x*x+a*x,), ('quadratic',))
    result = reconstruct_groebner(pres, trivial_sample(pres, {'a': 0}), point={'a': 0},
                                 max_primes=3, max_degree=4)
    assert result.certificate.basis == (x*x+a*x,)
    assert result.reduction.normal_form(x*x) == -a*x
    assert result.reduction.dimension == 2


def test_native_loss_of_sampled_leading_word_remains_inconclusive(native_finiteflow):
    field = RationalFunctionField(('a',))
    parent = FreeAlgebra(field, ('x',))
    x, a = parent.generator('x'), field.parameter('a')
    pres = AlgebraPresentation(parent, (a*x*x+x,), ('quadratic',))
    with pytest.raises(SupportFailure):
        reconstruct_groebner(pres, trivial_sample(pres, {'a': 0}), point={'a': 0},
                             max_primes=3, max_degree=4)


def test_native_underdetermined_node_preserves_the_numeric_monic_tail(native_finiteflow):
    field = RationalFunctionField(('a',))
    parent = FreeAlgebra(field, ('y', 'x'))
    y, x, a = parent.generator('y'), parent.generator('x'), field.parameter('a')
    pres = AlgebraPresentation(parent, (x-a*y, y), ('mixed', 'y'))
    numeric = pres.specialize({'a': 3})
    one = numeric.algebra.field.one
    sample = GroebnerCertificate(numeric.content_hash(),
        (numeric.relations[0]+2*numeric.relations[1], numeric.relations[1]), (
        (ProofTerm(one, (), 'input', 0, ()),
         ProofTerm(numeric.algebra.field.scalar(2), (), 'input', 1, ())),
        (ProofTerm(one, (), 'input', 1, ()),)))
    result = reconstruct_groebner(pres, sample, point={'a': 3}, max_primes=3,
                                 max_degree=4, threads=2, shift=(7,))
    assert result.certificate.basis[0] == x-y
    assert result.reduction.dimension == 1
