"""Exact generic reconstruction from independently specified numerical proof DAGs."""
from dataclasses import replace
from fractions import Fraction
import json

import pytest

from zipper_algebra import production
from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.matrices import ExactMatrix
from zipper_algebra.production import SupportFailure, reconstruct_groebner
from zipper_algebra.verification.groebner import GroebnerCertificate, ProofTerm, VerificationLimit, verify_groebner
from tests.unit.test_minimal_polynomial_producer import exact_kernel


def presentation(parent, relations):
    return AlgebraPresentation(parent, tuple(relations), tuple(f'r{i}' for i in range(len(relations))))


def input_sample(pres, point, indices=None):
    """Trivial local proofs, specified directly from the ordered literal relations."""
    specialized = pres.specialize(point)
    indices = tuple(range(len(specialized.relations))) if indices is None else tuple(indices)
    return GroebnerCertificate(specialized.content_hash(),
        tuple(specialized.relations[index] for index in indices),
        tuple((ProofTerm(specialized.algebra.field.one, (), 'input', index, ()),)
              for index in indices))


def clifford():
    field = RationalFunctionField(('alpha', 'beta', 'gamma'))
    a, b, c = (field.parameter(name) for name in field.parameters)
    parent = FreeAlgebra(field, ('u', 'v'))
    u, v = (parent.generator(name) for name in parent.generators)
    pres = presentation(parent, (u*u-a, v*v-b, v*u+u*v-c))
    return pres, u, v, (a, b, c)


@pytest.fixture
def oracle(monkeypatch):
    pytest.importorskip('sympy')
    calls = []
    def solve(matrix, **options):
        calls.append((matrix, options))
        return exact_kernel(matrix, **options)
    monkeypatch.setattr(production, 'reconstruct_kernel', solve)
    return calls


def test_three_parameter_noncommutative_family_has_four_standard_words(oracle):
    pres, u, v, (a, b, c) = clifford()
    point = {'gamma': Fraction(7, 2), 'alpha': 2, 'beta': 3}
    result = reconstruct_groebner(pres, input_sample(pres, point), point=point)
    assert result.reduction.dimension == 4
    assert result.reduction.standard_words == ((), (1,), (2,), (1, 2))
    assert result.reduction.normal_form(u*u) == pres.algebra.scalar(a)
    assert result.reduction.normal_form(v*v) == pres.algebra.scalar(b)
    assert result.reduction.normal_form(v*u) == -u*v+c
    assert result.reduction.normal_form(u*v) != result.reduction.normal_form(v*u)
    assert all(matrix.field.parameters == ('alpha', 'beta', 'gamma') for matrix, _ in oracle)
    assert len(oracle) == 3
    restored = GroebnerCertificate.from_json(json.loads(json.dumps(result.certificate.to_json())))
    assert verify_groebner(pres, restored).dimension == 4


def test_nonmonic_numerical_node_reconstructs_rational_parameter_coefficients(oracle):
    field = RationalFunctionField(('scale', 'tail'))
    scale, tail = field.parameter('scale'), field.parameter('tail')
    parent = FreeAlgebra(field, ('letter',))
    x = parent.generator('letter')
    pres = presentation(parent, (scale*x*x-tail,))
    point = {'scale': 2, 'tail': Fraction(5, 3)}
    sample = input_sample(pres, point)
    result = reconstruct_groebner(pres, sample, point=point,
        max_primes=3, max_degree=8, threads=2, start_prime=1, shift=[7, 9], max_products=4567)
    assert result.reduction.normal_form(x*x) == parent.scalar(tail/scale)
    assert result.certificate.basis[0] == x*x-tail/scale
    assert result.certificate.derivations[0][0].coefficient == field.one/scale
    assert all(options == {'max_primes': 3, 'max_degree': 8, 'threads': 2,
                           'start_prime': 1, 'shift': (7, 9), 'max_products': 4567}
               for _, options in oracle)


def test_forward_reference_retains_noncommuting_contexts_without_cofactor_expansion(oracle):
    pres, u, v, _ = clifford()
    point = {'alpha': 2, 'beta': 3, 'gamma': 5}
    numeric = pres.specialize(point)
    r0, r1, r2 = numeric.relations
    nu, nv = (numeric.algebra.generator(name) for name in numeric.algebra.generators)
    one = numeric.algebra.field.one
    sample = GroebnerCertificate(numeric.content_hash(), (nv*r0*nu, r0, r1, r2), (
        (ProofTerm(one, (2,), 'basis', 1, (1,)),),
        (ProofTerm(one, (), 'input', 0, ()),),
        (ProofTerm(one, (), 'input', 1, ()),),
        (ProofTerm(one, (), 'input', 2, ()),),
    ))
    result = reconstruct_groebner(pres, sample, point=point)
    assert result.reduction.dimension == 4
    assert result.certificate.basis[0] == v*pres.relations[0]*u
    refs = result.certificate.derivations[0]
    assert any(term.kind == 'basis' and term.index == 1 and term.left == (2,)
               and term.right == (1,) and not term.coefficient.is_zero() for term in refs)
    assert not any(term.kind == 'input' for term in refs)
    assert len(result.certificate.basis) == 4


def test_shared_graph_stays_linear_in_nodes_instead_of_expanding_all_paths(oracle):
    parent = FreeAlgebra(RationalFunctionField(('a',)), ('x',))
    x, a = parent.generator('x'), parent.field.parameter('a')
    pres = presentation(parent, (x*x-a,))
    numeric = pres.specialize({'a': 2})
    one = numeric.algebra.field.one
    count = 12
    nodes = tuple((2**index)*numeric.relations[0] for index in range(count))
    rows = ((ProofTerm(one, (), 'input', 0, ()),),) + tuple(
        (ProofTerm(one, (), 'basis', index-1, ()), ProofTerm(one, (), 'basis', index-1, ()))
        for index in range(1, count))
    sample = GroebnerCertificate(numeric.content_hash(), nodes, rows)
    result = reconstruct_groebner(pres, sample, point={'a': 2})
    assert result.reduction.dimension == 2
    assert len(result.certificate.basis) == count
    assert sum(map(len, result.certificate.derivations)) <= 2*count
    assert all(any(term.kind == 'basis' and term.index == index-1 for term in row)
               for index, row in enumerate(result.certificate.derivations) if index)


def test_structural_context_support_recovers_a_tail_canceled_at_the_sample(oracle):
    parent = FreeAlgebra(RationalFunctionField(('a',)), ('x',))
    x, a = parent.generator('x'), parent.field.parameter('a')
    pres = presentation(parent, (x*x+a*x,))
    sample = input_sample(pres, {'a': 0})
    assert sample.basis[0].terms == (((1, 1), sample.basis[0].algebra.field.one),)
    result = reconstruct_groebner(pres, sample, point={'a': 0})
    assert result.certificate.basis[0] == x*x+a*x
    assert result.reduction.normal_form(x*x) == -a*x
    assert result.reduction.dimension == 2


@pytest.mark.parametrize('leading_loss', ['degree', 'lexicographic'])
def test_structural_tail_expansion_never_inserts_words_above_the_sampled_leading_word(oracle, leading_loss):
    parent = FreeAlgebra(RationalFunctionField(('a',)), ('x', 'y'))
    x, y = parent.generator('x'), parent.generator('y')
    a = parent.field.parameter('a')
    relation = a*x*x+x if leading_loss == 'degree' else a*y+x
    pres = presentation(parent, (relation,))
    sample = input_sample(pres, {'a': 0})
    assert sample.basis[0] == sample.basis[0].algebra.generator('x')
    with pytest.raises(SupportFailure):
        reconstruct_groebner(pres, sample, point={'a': 0})
    # The unwanted higher word is an equation row, not a proposed target
    # coefficient column that could silently replace the sampled leading word.
    generic_systems = [matrix for matrix, _ in oracle if matrix.field.parameters]
    assert len(generic_systems) == 1
    assert generic_systems[0].ncols == 2
    assert generic_systems[0].nrows == 2


def test_canceled_noncommuting_context_tails_are_recovered_through_a_shared_parent(oracle):
    pres, u, v, (a, _b, c) = clifford()
    point = {'alpha': 0, 'beta': 3, 'gamma': 0}
    numeric = pres.specialize(point)
    r0, r1, r2 = numeric.relations
    nu, nv = numeric.algebra.generator('u'), numeric.algebra.generator('v')
    one = numeric.algebra.field.one
    sample = GroebnerCertificate(numeric.content_hash(), (nv*r0*nu, r0, r1, r2), (
        (ProofTerm(one, (2,), 'basis', 1, (1,)),),
        (ProofTerm(one, (), 'input', 0, ()),),
        (ProofTerm(one, (), 'input', 1, ()),),
        (ProofTerm(one, (), 'input', 2, ()),),
    ))
    result = reconstruct_groebner(pres, sample, point=point)
    assert result.certificate.basis[0] == v*u*u*u-a*v*u
    assert result.reduction.normal_form(v*u) == -u*v+c
    assert result.reduction.dimension == 4
    assert all(term.kind == 'basis' and term.index == 1
               for term in result.certificate.derivations[0])


def test_empty_specialized_ideal_does_not_certify_the_generic_free_algebra(oracle):
    parent = FreeAlgebra(RationalFunctionField(('a',)), ('x',))
    x, a = parent.generator('x'), parent.field.parameter('a')
    pres = presentation(parent, (a*x,))
    with pytest.raises((RuntimeError, ValueError)):
        reconstruct_groebner(pres, input_sample(pres, {'a': 0}, ()), point={'a': 0}, max_words=3)


def test_local_input_membership_does_not_replace_generic_completion(oracle):
    parent = FreeAlgebra(RationalFunctionField(), ('a', 'b'))
    a, b = (parent.generator(name) for name in parent.generators)
    pres = presentation(parent, (a*b*a-b,))
    with pytest.raises(ValueError):
        reconstruct_groebner(pres, input_sample(pres, {}), point={}, max_words=8)


def test_omitted_relation_does_not_pass_reverse_ideal_inclusion(oracle):
    pres, *_ = clifford()
    point = {'alpha': 2, 'beta': 3, 'gamma': 5}
    with pytest.raises(ValueError):
        reconstruct_groebner(pres, input_sample(pres, point, (0, 1)), point=point, max_words=8)


def test_unknown_word_enumeration_is_not_a_false_finiteness_verdict(oracle):
    parent = FreeAlgebra(RationalFunctionField(), ('x',))
    pres = presentation(parent, ())
    result = reconstruct_groebner(pres, input_sample(pres, {}), point={}, max_words=3)
    assert result.reduction.dimension is None
    assert result.reduction.standard_words is None
    finite = presentation(parent, (parent.generator('x')**2,))
    bounded = reconstruct_groebner(finite, input_sample(finite, {}), point={}, max_words=1)
    assert bounded.reduction.dimension is None


def test_unit_ideal_and_scalar_free_algebra_have_distinct_exact_dimensions(oracle):
    parent = FreeAlgebra(RationalFunctionField(), ('x',))
    unit = presentation(parent, (parent.scalar(3),))
    zero = reconstruct_groebner(unit, input_sample(unit, {}), point={})
    assert zero.reduction.dimension == 0
    scalar_parent = FreeAlgebra(RationalFunctionField(), ())
    scalar = presentation(scalar_parent, ())
    assert reconstruct_groebner(scalar, input_sample(scalar, {}), point={}).reduction.dimension == 1


@pytest.mark.parametrize('kind', ['bad-identity', 'cycle', 'invalid-zero-reference'])
def test_false_numeric_proofs_fail_before_reconstruction(monkeypatch, kind):
    parent = FreeAlgebra(RationalFunctionField(), ('x',))
    x = parent.generator('x')
    pres = presentation(parent, (x*x,))
    cert = input_sample(pres, {})
    if kind == 'bad-identity':
        cert = replace(cert, basis=(x,))
    elif kind == 'cycle':
        cert = replace(cert, derivations=((ProofTerm(parent.field.one, (), 'basis', 0, ()),),))
    else:
        cert = replace(cert, derivations=((ProofTerm(parent.field.one, (), 'input', 0, ()),
                                          ProofTerm(parent.field.zero, (), 'input', 1, ())),))
    monkeypatch.setattr(production, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises(ValueError):
        reconstruct_groebner(pres, cert, point={})


def test_sample_is_bound_to_exact_named_specialization_before_native_work(monkeypatch):
    pres, *_ = clifford()
    point = {'alpha': 2, 'beta': 3, 'gamma': 5}
    sample = input_sample(pres, point)
    monkeypatch.setattr(production, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises(ValueError):
        reconstruct_groebner(pres, sample, point={**point, 'alpha': 7})
    foreign_parent = FreeAlgebra(RationalFunctionField(), ('u', 'other'))
    foreign = presentation(foreign_parent, (foreign_parent.generator('u')**2,))
    with pytest.raises((TypeError, ValueError)):
        reconstruct_groebner(pres, input_sample(foreign, {}), point=point)


def test_positive_characteristic_is_explicitly_unsupported(monkeypatch):
    parent = FreeAlgebra(RationalFunctionField((), 7), ('x',))
    pres = presentation(parent, (parent.generator('x')**2,))
    monkeypatch.setattr(production, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises(ValueError, match='characteristic zero'):
        reconstruct_groebner(pres, input_sample(pres, {}), point={})


@pytest.mark.parametrize('limit,value', [('max_nodes', 0), ('max_nodes', True),
    ('max_operations', 0), ('max_operations', True), ('max_reductions', 0),
    ('max_words', 0), ('max_primes', 0), ('max_degree', 0), ('threads', 0),
    ('start_prime', -1), ('shift', (-1,))])
def test_invalid_limits_fail_before_native_work(monkeypatch, limit, value):
    parent = FreeAlgebra(RationalFunctionField(), ('x',))
    pres = presentation(parent, (parent.generator('x')**2,))
    monkeypatch.setattr(production, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises((ValueError, TypeError)):
        reconstruct_groebner(pres, input_sample(pres, {}), point={}, **{limit: value})


@pytest.mark.parametrize('option,value', [('max_nodes', 1), ('max_operations', 1), ('max_reductions', 1)])
def test_actual_work_limits_remain_inconclusive(oracle, option, value):
    pres, *_ = clifford()
    point = {'alpha': 2, 'beta': 3, 'gamma': 5}
    with pytest.raises(VerificationLimit):
        reconstruct_groebner(pres, input_sample(pres, point), point=point, **{option: value})


def test_backend_failure_is_not_reported_as_nonmembership(monkeypatch):
    parent = FreeAlgebra(RationalFunctionField(), ('x',))
    pres = presentation(parent, (parent.generator('x')**2,))
    def fail(*args, **kwargs):
        raise RuntimeError('bounded native search exhausted')
    monkeypatch.setattr(production, 'reconstruct_kernel', fail)
    with pytest.raises(RuntimeError, match='bounded native'):
        reconstruct_groebner(pres, input_sample(pres, {}), point={})


def test_stale_valid_kernel_cannot_certify_a_different_local_system(oracle, monkeypatch):
    parent = FreeAlgebra(RationalFunctionField(), ('x',))
    pres = presentation(parent, (parent.generator('x')**2,))
    def stale(matrix, **options):
        # An honestly certified zero-system kernel, returned for a different
        # nonzero local identity system, models incorrect cache reuse.
        return exact_kernel(ExactMatrix(matrix.field, matrix.nrows, matrix.ncols))
    monkeypatch.setattr(production, 'reconstruct_kernel', stale)
    with pytest.raises((ValueError, RuntimeError)):
        reconstruct_groebner(pres, input_sample(pres, {}), point={})


@pytest.mark.parametrize('point', [None, {}, {'alpha': 1, 'beta': 2},
    {'alpha': 1, 'beta': 2, 'gamma': 3, 'extra': 4},
    {'alpha': True, 'beta': 2, 'gamma': 3},
    {'alpha': 1.5, 'beta': 2, 'gamma': 3}])
def test_parameter_assignment_requires_every_named_exact_rational(monkeypatch, point):
    pres, *_ = clifford()
    sample = input_sample(pres, {'alpha': 2, 'beta': 3, 'gamma': 5})
    monkeypatch.setattr(production, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises((TypeError, ValueError)):
        reconstruct_groebner(pres, sample, point=point)


@pytest.mark.parametrize('bad_pres,bad_sample', [(None, None), (None, 'bad'), ('bad', None)])
def test_untyped_input_is_rejected_before_any_native_work(monkeypatch, bad_pres, bad_sample):
    monkeypatch.setattr(production, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises(TypeError):
        reconstruct_groebner(bad_pres, bad_sample, point={})


def test_literal_coefficient_poles_are_not_silently_canceled(monkeypatch):
    parent = FreeAlgebra(RationalFunctionField(('a',)), ('x',))
    a, x = parent.field.parameter('a'), parent.generator('x')
    pres = presentation(parent, (x*x-parent.scalar(parent.field.one/a),))
    safe = input_sample(pres, {'a': 2})
    monkeypatch.setattr(production, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises(ZeroDivisionError):
        reconstruct_groebner(pres, safe, point={'a': 0})


@pytest.mark.parametrize('options', [
    {'max_primes': 2048}, {'max_primes': 1, 'start_prime': 2047},
    {'max_degree': 2**16}, {'threads': 2**32}, {'threads': True},
    {'shift': (True,)}, {'shift': (2**64,)}, {'shift': (-1,)},
    {'shift': ()}, {'shift': (1, 2)}, {'shift': 3},
])
def test_native_options_are_validated_even_for_an_empty_proof_graph(monkeypatch, options):
    parent = FreeAlgebra(RationalFunctionField(('a',)), ())
    pres = presentation(parent, ())
    sample = input_sample(pres, {'a': 2})
    monkeypatch.setattr(production, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises((ValueError, TypeError)):
        reconstruct_groebner(pres, sample, point={'a': 2}, **options)


def test_aggregate_generic_placement_budget_can_exhaust_after_numeric_proof_passes(oracle):
    pres, *_ = clifford()
    point = {'alpha': 2, 'beta': 3, 'gamma': 5}
    # The numerical graph costs ten units. Generic systems also place the
    # candidate support monomials, exhausting the separately bounded stage.
    with pytest.raises(VerificationLimit, match='operation budget'):
        reconstruct_groebner(pres, input_sample(pres, point), point=point, max_operations=10)
    assert 0 < len(oracle) < 3


def test_underdetermined_node_uses_an_exact_constant_combination_to_preserve_the_anchor(oracle):
    field = RationalFunctionField(('a',))
    parent = FreeAlgebra(field, ('y', 'x'))
    y, x, a = parent.generator('y'), parent.generator('x'), field.parameter('a')
    pres = presentation(parent, (x-a*y, y))
    numeric = pres.specialize({'a': 3})
    one = numeric.algebra.field.one
    sample = GroebnerCertificate(numeric.content_hash(),
        (numeric.relations[0]+2*numeric.relations[1], numeric.relations[1]), (
        (ProofTerm(one, (), 'input', 0, ()),
         ProofTerm(numeric.algebra.field.scalar(2), (), 'input', 1, ())),
        (ProofTerm(one, (), 'input', 1, ()),)))
    result = reconstruct_groebner(pres, sample, point={'a': 3}, threads=2,
                                 max_primes=3, max_degree=8, start_prime=1, shift=[7])
    # The generic kernel has y and x vectors, neither of which separately
    # specializes to x-y. Their exact constant combination preserves that node.
    assert result.certificate.basis[0] == x-y
    assert result.reduction.dimension == 1
    assert {term.index: term.coefficient for term in result.certificate.derivations[0]} == {
        0: field.one, 1: a-1}
    ground_solves = [(matrix, options) for matrix, options in oracle if not matrix.field.parameters]
    assert len(ground_solves) == 1
    matrix, options = ground_solves[0]
    assert (matrix.nrows, matrix.ncols) == (2, 3)
    assert options == {'max_products': 10000000, 'max_primes': 3, 'max_degree': 8,
                       'threads': 2, 'start_prime': 1}
    assert all(options.get('shift') == (7,) for matrix, options in oracle if matrix.field.parameters)


def test_a_valid_kernel_basis_with_an_anchor_pole_fails_without_a_generic_verdict(oracle, monkeypatch):
    from zipper_algebra.verification.linear import KernelCertificate, verify_kernel
    field = RationalFunctionField(('a',))
    parent = FreeAlgebra(field, ('y', 'x'))
    y, x, a = parent.generator('y'), parent.generator('x'), field.parameter('a')
    pres = presentation(parent, (x+a*y,))
    alternate = []
    def different_kernel_basis(matrix, **options):
        if matrix.field.parameters and matrix.ncols == 3:
            # An honest complete kernel with identity coordinate on the y
            # target, rather than on x, has x coefficient 1/a. Its rank witness
            # is taken at a=2, where the chosen minor is invertible.
            kernel = ExactMatrix(field, 3, 1, ((0, 0, field.one/a),
                (1, 0, field.one), (2, 0, field.one/a)))
            certificate = KernelCertificate(matrix.content_hash(), kernel, (1,), (0, 1),
                (0, 2), 101, (2,), ((51, 0), (51, 100)))
            alternate.append(verify_kernel(matrix, certificate))
            return certificate
        return exact_kernel(matrix, **options)
    monkeypatch.setattr(production, 'reconstruct_kernel', different_kernel_basis)
    with pytest.raises(SupportFailure, match='regular'):
        reconstruct_groebner(pres, input_sample(pres, {'a': 0}), point={'a': 0})
    assert len(alternate) == 1


def test_one_pole_column_does_not_hide_another_regular_generic_proof(oracle, monkeypatch):
    from zipper_algebra.verification.linear import KernelCertificate, verify_kernel
    field = RationalFunctionField(('a',))
    parent = FreeAlgebra(field, ('y', 'x'))
    y, x, a = parent.generator('y'), parent.generator('x'), field.parameter('a')
    pres = presentation(parent, (x+a*y, y, x))
    numeric = pres.specialize({'a': 0})
    one, zero = numeric.algebra.field.one, numeric.algebra.field.zero
    sample = GroebnerCertificate(numeric.content_hash(), numeric.relations, (
        (ProofTerm(one, (), 'input', 0, ()), ProofTerm(zero, (), 'input', 1, ()),
         ProofTerm(zero, (), 'input', 2, ())),
        (ProofTerm(one, (), 'input', 1, ()),),
        (ProofTerm(one, (), 'input', 2, ()),)))
    alternate = []
    def mixed_kernel(matrix, **options):
        if matrix.field.parameters and matrix.ncols == 4:
            # One complete-kernel column has target coefficient -1/a;
            # another is the regular direct proof of x from input relation 2.
            kernel = ExactMatrix(field, 4, 2, ((0, 0, -field.one/a),
                (1, 0, field.one), (3, 0, -field.one/a),
                (2, 1, field.one), (3, 1, field.one)))
            certificate = KernelCertificate(matrix.content_hash(), kernel, (1, 2),
                (0, 1), (0, 3), 101, (2,), ((51, 0), (51, 100)))
            alternate.append(verify_kernel(matrix, certificate))
            return certificate
        return exact_kernel(matrix, **options)
    monkeypatch.setattr(production, 'reconstruct_kernel', mixed_kernel)
    result = reconstruct_groebner(pres, sample, point={'a': 0})
    assert result.certificate.basis[0] == x
    assert result.reduction.dimension == 1
    assert len(alternate) == 1
    assert result.certificate.derivations[0] == (ProofTerm(field.one, (), 'input', 2, ()),)
