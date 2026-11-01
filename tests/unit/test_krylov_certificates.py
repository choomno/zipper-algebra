"""Literal cyclic-orbit witnesses, independent of the reconstruction producer."""
from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
import json

import pytest

from tests.unit.test_finite_algebra import matrix_algebra_certificate
from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.matrices import ExactMatrix
from zipper_algebra.verification.derivations import (
    DerivationGraph, DerivationTerm, VerificationLimit,
)
from zipper_algebra.verification.finite import (
    GeneratorAction, WordBasisCertificate, verify_word_basis,
)
from zipper_algebra.verification.krylov import (
    MinimalPolynomialCertificate, krylov_matrix, verify_minimal_polynomial,
)
from zipper_algebra.verification.linear import KernelCertificate, verify_kernel


def _monogenic(lower, field=None):
    """K[x]/(x**d + sum(lower[i]*x**i)), with literal column proofs."""
    field = RationalFunctionField() if field is None else field
    if field.parameters:
        pytest.importorskip('sympy')
    degree = len(lower)
    parent = FreeAlgebra(field, ('x',))
    relation = parent.element({(1,) * degree: field.one,
                               **{(1,) * i: field.scalar(c) for i, c in enumerate(lower)}})
    presentation = AlgebraPresentation(parent, (relation,), ('polynomial',))
    columns = tuple(((i + 1, field.one),) for i in range(degree - 1)) + (
        tuple((i, -field.scalar(c)) for i, c in enumerate(lower)),)
    action = GeneratorAction(field, columns)
    graph = DerivationGraph(presentation.content_hash(), (parent.zero, relation),
                            ((), (DerivationTerm(field.one, (), 'input', 0, ()),)))
    certificate = WordBasisCertificate(
        presentation.content_hash(), tuple((1,) * i for i in range(degree)),
        (action,), graph, ((0,) * (degree - 1) + (1,),))
    return verify_word_basis(presentation, certificate)


def _matrix_algebra(characteristic=0):
    presentation, certificate, _, _ = matrix_algebra_certificate(characteristic)
    return verify_word_basis(presentation, certificate)


def _triangular_algebra():
    """UT_2(K) in the basis (1,p,n), p=E11 and n=E12."""
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('p', 'n'))
    p, n = (parent.generator(name) for name in parent.generators)
    relations = (p*p-p, n*n, p*n-n, n*p)
    presentation = AlgebraPresentation(parent, relations, ('idempotent', 'nilpotent', 'pn', 'np'))
    def action(columns):
        return GeneratorAction(field, tuple(tuple((row, field.scalar(value)) for row, value in column)
                                            for column in columns))
    p_action = action((((1, 1),), ((1, 1),), ((2, 1),)))
    n_action = action((((2, 1),), (), ()))
    graph = DerivationGraph(presentation.content_hash(), (parent.zero,) + relations,
                            ((),) + tuple((DerivationTerm(field.one, (), 'input', i, ()),)
                                          for i in range(4)))
    certificate = WordBasisCertificate(presentation.content_hash(), ((), (1,), (2,)),
                                       (p_action, n_action), graph, ((0, 1, 3), (0, 4, 2)))
    return verify_word_basis(presentation, certificate)


def _zero_algebra():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    presentation = AlgebraPresentation(parent, (parent.one,), ('unit',))
    graph = DerivationGraph(presentation.content_hash(), (parent.one,),
                            ((DerivationTerm(field.one, (), 'input', 0, ()),),))
    certificate = WordBasisCertificate(presentation.content_hash(), (),
                                       (GeneratorAction(field, ()),), graph, ((),), unit_node=0)
    return verify_word_basis(presentation, certificate)


def _matrix(field, columns, dimension):
    return ExactMatrix(field, dimension, len(columns), tuple(
        (row, column, field.scalar(value)) for column, vector in enumerate(columns)
        for row, value in enumerate(vector)))


def _certificate(algebra, element, seed, coefficients, columns, rows, *,
                 inverse=None, pivot_columns=None, identity_row=None, vector=None):
    """Use a hand-selected minor, with no elimination or reconstruction routine."""
    field = algebra.field
    degree = len(coefficients) - 1
    original = _matrix(field, columns, algebra.dimension)
    polynomial = tuple(field.scalar(c) for c in coefficients)
    vector = polynomial if vector is None else tuple(field.scalar(c) for c in vector)
    kernel = ExactMatrix(field, degree + 1, 1,
                         tuple((row, 0, value) for row, value in enumerate(vector)))
    witness = KernelCertificate(original.content_hash(), kernel,
        (degree if identity_row is None else identity_row,), tuple(rows),
        tuple(range(degree)) if pivot_columns is None else tuple(pivot_columns),
        field.characteristic or 101, (1,) * len(field.parameters),
        tuple(tuple(int(i == j) for j in range(degree)) for i in range(degree))
        if inverse is None else inverse)
    certificate = MinimalPolynomialCertificate(algebra.presentation.content_hash(), algebra.basis_words,
        tuple(field.scalar(c) for c in element), tuple(field.scalar(c) for c in seed),
        polynomial, witness)
    return original, certificate


def _dual_certificate():
    algebra = _monogenic((0, 0))
    matrix, certificate = _certificate(algebra, (0, 1), (1, 0), (0, 0, 1),
                                        ((1, 0), (0, 1), (0, 0)), (0, 1))
    return algebra, matrix, certificate


def test_noncommuting_krylov_columns_use_left_multiplication():
    algebra = _matrix_algebra()
    x, y = algebra.element((0, 1, 0, 0)), algebra.element((0, 0, 1, 0))
    expected = _matrix(algebra.field, ((0, 0, 1, 0), (0, 0, 0, 1), (0, 0, 0, 0)), 4)
    assert krylov_matrix(algebra, x, y, 2).to_json() == expected.to_json()
    # y*x=1-xy would give a different second column.
    assert algebra.multiply(y, x).coordinates == tuple(algebra.field.scalar(c) for c in (1, 0, 0, -1))


@pytest.mark.parametrize('characteristic', [0, 2, 7])
def test_matrix_unit_relative_nilpotent_polynomial(characteristic):
    algebra = _matrix_algebra(characteristic)
    matrix, certificate = _certificate(algebra, (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 1),
        ((0, 0, 1, 0), (0, 0, 0, 1), (0, 0, 0, 0)), (2, 3))
    assert verify_kernel(matrix, certificate.kernel) is certificate.kernel
    assert verify_minimal_polynomial(algebra, certificate) is certificate
    assert certificate.degree == 2


def test_relative_polynomial_and_element_polynomial_are_different():
    algebra = _matrix_algebra()
    # e=xy is a proper idempotent. On x it acts as the identity; on 1 its
    # minimal polynomial is lambda*(lambda-1).
    _, relative = _certificate(algebra, (0, 0, 0, 1), (0, 1, 0, 0), (-1, 1),
                                ((0, 1, 0, 0), (0, 1, 0, 0)), (1,))
    _, full = _certificate(algebra, (0, 0, 0, 1), (1, 0, 0, 0), (0, -1, 1),
        ((1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 1)), (0, 3))
    assert verify_minimal_polynomial(algebra, relative) is relative
    assert verify_minimal_polynomial(algebra, full) is full
    assert relative.degree == 1 and full.degree == 2
    # Replacing the seed, while keeping its old valid relation, cannot turn a
    # relative assertion into an assertion about the element itself.
    with pytest.raises(ValueError):
        verify_minimal_polynomial(algebra, replace(relative, seed=full.seed))


def test_triangular_left_action_is_distinguished_from_right_action():
    algebra = _triangular_algebra()
    # a=p+n. Then a*p=p while p*a=p+n, so swapping the multiplication order
    # would invalidate this literal relative polynomial.
    matrix, certificate = _certificate(algebra, (0, 1, 1), (0, 1, 0), (-1, 1),
                                        ((0, 1, 0), (0, 1, 0)), (1,))
    assert verify_kernel(matrix, certificate.kernel) is certificate.kernel
    assert verify_minimal_polynomial(algebra, certificate) is certificate
    assert krylov_matrix(algebra, algebra.element(certificate.element),
                         algebra.element(certificate.seed), 1).to_json() == matrix.to_json()


def test_nonminimal_annihilator_is_rejected_despite_a_complete_valid_kernel():
    algebra = _matrix_algebra()
    field = algebra.field
    # For a=1, (lambda-1)^2 annihilates every seed but is not minimal.
    original = _matrix(field, ((1, 0, 0, 0),) * 3, 4)
    kernel = ExactMatrix.from_rows(field, [[-1, -1], [1, 0], [0, 1]])
    witness = KernelCertificate(original.content_hash(), kernel, (1, 2), (0,), (0,),
                                 101, (), ((1,),))
    assert verify_kernel(original, witness) is witness
    false = MinimalPolynomialCertificate(algebra.presentation.content_hash(), algebra.basis_words,
        algebra.one.coordinates, algebra.one.coordinates,
        tuple(field.scalar(c) for c in (1, -2, 1)), witness)
    with pytest.raises(ValueError):
        verify_minimal_polynomial(algebra, false)


def test_dual_numbers_relation_has_a_genuine_nonzero_lower_power():
    algebra, matrix, certificate = _dual_certificate()
    assert certificate.element != algebra.zero.coordinates
    assert matrix.nrows == 2 and matrix.ncols == 3
    assert verify_minimal_polynomial(algebra, certificate) is certificate


def test_degree_three_and_shorter_seed_orbits():
    algebra = _monogenic((0, 0, 0))
    _, full = _certificate(algebra, (0, 1, 0), (1, 0, 0), (0, 0, 0, 1),
        ((1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)), (0, 1, 2))
    _, relative = _certificate(algebra, (0, 1, 0), (0, 1, 0), (0, 0, 1),
        ((0, 1, 0), (0, 0, 1), (0, 0, 0)), (1, 2))
    assert verify_minimal_polynomial(algebra, full) is full
    assert verify_minimal_polynomial(algebra, relative) is relative


@pytest.mark.parametrize('characteristic,parameters', [(0, ('a', 'b', 'c')), (2, ('s',)), (3, ('s',))])
def test_generic_parameter_polynomials_remain_over_the_declared_field(characteristic, parameters):
    pytest.importorskip('sympy')
    field = RationalFunctionField(parameters, characteristic)
    if len(parameters) == 3:
        a, b, c = (field.parameter(name) for name in parameters)
        lower = (c, a / b)
    else:
        s = field.parameter('s')
        lower = (field.one, s)
    algebra = _monogenic(lower, field)
    matrix, certificate = _certificate(algebra, (0, 1), (1, 0), (*lower, 1),
                                        ((1, 0), (0, 1), (-lower[0], -lower[1])), (0, 1))
    assert verify_kernel(matrix, certificate.kernel) is certificate.kernel
    assert verify_minimal_polynomial(algebra, certificate) is certificate
    assert certificate.coefficients == tuple(field.scalar(c) for c in (*lower, 1))


def test_irreducible_quadratic_over_f2_is_a_two_dimensional_algebra():
    field = RationalFunctionField((), 2)
    algebra = _monogenic((1, 1), field)
    _, certificate = _certificate(algebra, (0, 1), (1, 0), (1, 1, 1),
                                   ((1, 0), (0, 1), (1, 1)), (0, 1))
    assert verify_minimal_polynomial(algebra, certificate) is certificate
    assert certificate.degree == 2


def test_nonunit_last_kernel_coordinate_is_normalized_exactly():
    algebra = _monogenic((0, 0))
    # a=2+epsilon has polynomial (lambda-2)^2. A valid kernel can be
    # normalized at its first coordinate, rather than its last coordinate.
    matrix, certificate = _certificate(algebra, (2, 1), (1, 0), (4, -4, 1),
        ((1, 0), (2, 1), (4, 4)), (0, 1), inverse=((1, 100), (25, 51)),
        pivot_columns=(1, 2), identity_row=0, vector=(1, -1, Fraction(1, 4)))
    assert verify_kernel(matrix, certificate.kernel) is certificate.kernel
    assert verify_minimal_polynomial(algebra, certificate) is certificate


def test_nonunit_seed_coordinates_and_a_nonidentity_rank_minor():
    algebra = _monogenic((0, 0))
    # v=2+epsilon, a=1+epsilon: K has columns (2,1),(2,3),(2,5).
    matrix, certificate = _certificate(algebra, (1, 1), (2, 1), (1, -2, 1),
        ((2, 1), (2, 3), (2, 5)), (0, 1), inverse=((26, 50), (25, 51)))
    assert verify_kernel(matrix, certificate.kernel) is certificate.kernel
    assert verify_minimal_polynomial(algebra, certificate) is certificate


@pytest.mark.parametrize('zero_parent', [False, True])
def test_zero_seed_has_polynomial_one_including_the_zero_algebra(zero_parent):
    algebra = _zero_algebra() if zero_parent else _monogenic((0, 0))
    field = algebra.field
    matrix, certificate = _certificate(algebra, algebra.one.coordinates,
        algebra.zero.coordinates, (1,), (algebra.zero.coordinates,), ())
    assert matrix.ncols == 1 and matrix.nrows == algebra.dimension
    assert certificate.degree == 0 and certificate.coefficients == (field.one,)
    assert verify_minimal_polynomial(algebra, certificate, max_products=0) is certificate


def test_certificate_can_move_between_independent_verifications_of_the_same_basis():
    algebra, _, certificate = _dual_certificate()
    second = _monogenic((0, 0))
    assert second is not algebra
    assert verify_minimal_polynomial(second, certificate) is certificate
    with pytest.raises(TypeError):
        krylov_matrix(second, algebra.one, second.one, 1)
    with pytest.raises(TypeError):
        krylov_matrix(second, second.one, algebra.one, 1)


def test_certificate_binding_rejects_a_different_presentation_or_basis():
    algebra, _, certificate = _dual_certificate()
    with pytest.raises(ValueError):
        verify_minimal_polynomial(algebra, replace(certificate, presentation_hash='f' * 64))
    matrix_algebra = _matrix_algebra()
    _, certificate = _certificate(matrix_algebra, (0, 0, 0, 1), (1, 0, 0, 0), (0, -1, 1),
        ((1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 1)), (0, 3))
    words = ((), (2,), (1,), (1, 2))
    with pytest.raises(ValueError):
        verify_minimal_polynomial(matrix_algebra, replace(certificate, basis_words=words))


def test_wrong_polynomial_or_false_native_kernel_is_rejected():
    algebra, matrix, certificate = _dual_certificate()
    field = algebra.field
    with pytest.raises(ValueError):
        verify_minimal_polynomial(algebra, replace(certificate,
            coefficients=(field.one, field.zero, field.one)))
    false_vector = ExactMatrix.from_rows(field, [[1], [0], [1]])
    with pytest.raises(ValueError):
        verify_minimal_polynomial(algebra, replace(certificate,
            kernel=replace(certificate.kernel, kernel=false_vector)))
    assert matrix.content_hash() == certificate.kernel.matrix_hash


def test_unrelated_matrix_hash_and_corrupt_rank_witness_are_rejected():
    algebra, _, certificate = _dual_certificate()
    for wrong in [replace(certificate.kernel, matrix_hash='a' * 64),
                  replace(certificate.kernel, minor_inverse=((1, 0), (0, 2)))]:
        with pytest.raises(ValueError):
            verify_minimal_polynomial(algebra, replace(certificate, kernel=wrong))


def test_monicity_is_checked_as_a_mathematical_requirement():
    algebra, _, certificate = _dual_certificate()
    field = algebra.field
    for leading in (field.zero, field.scalar(2)):
        with pytest.raises(ValueError):
            verify_minimal_polynomial(algebra, replace(certificate,
                coefficients=(field.zero, field.zero, leading)))


@pytest.mark.parametrize('degree', [-1, True, Fraction(1, 2), 3])
def test_invalid_krylov_degrees_are_rejected(degree):
    algebra, _, _ = _dual_certificate()
    with pytest.raises(ValueError):
        krylov_matrix(algebra, algebra.one, algebra.one, degree)


@pytest.mark.parametrize('limit', [0, True, -1, 1.5])
def test_operation_limits_are_exact_positive_integers(limit):
    algebra, _, certificate = _dual_certificate()
    with pytest.raises(ValueError):
        krylov_matrix(algebra, algebra.one, algebra.one, 1, max_operations=limit)
    with pytest.raises(ValueError):
        verify_minimal_polynomial(algebra, certificate, max_operations=limit)


@pytest.mark.parametrize('limit', [-1, True, 1.5])
def test_kernel_product_limits_are_exact_nonnegative_integers(limit):
    algebra, _, certificate = _dual_certificate()
    with pytest.raises(ValueError):
        verify_minimal_polynomial(algebra, certificate, max_products=limit)


def test_exhausted_counted_work_never_returns_a_verified_polynomial():
    algebra, _, certificate = _dual_certificate()
    with pytest.raises(VerificationLimit):
        krylov_matrix(algebra, algebra.one, algebra.one, 1, max_operations=1)
    with pytest.raises(VerificationLimit):
        verify_minimal_polynomial(algebra, certificate, max_operations=1)
    with pytest.raises(RuntimeError, match='budget'):
        verify_minimal_polynomial(algebra, certificate, max_products=0)


def test_wrong_api_types_and_foreign_parents_fail_before_arithmetic():
    algebra, _, certificate = _dual_certificate()
    with pytest.raises(TypeError):
        krylov_matrix(object(), algebra.one, algebra.one, 1)
    with pytest.raises(TypeError):
        krylov_matrix(algebra, algebra.one.coordinates, algebra.one, 1)
    with pytest.raises(TypeError):
        verify_minimal_polynomial(object(), certificate)
    with pytest.raises(TypeError):
        verify_minimal_polynomial(algebra, certificate.kernel)


def test_model_is_frozen_and_preserves_literal_coordinate_order():
    _, _, certificate = _dual_certificate()
    with pytest.raises(FrozenInstanceError):
        certificate.element = ()
    wire = certificate.to_json()
    assert wire['schema'] == 'zipper-minimal-polynomial-certificate-v1'
    assert wire['element'][0][0] == 1 and wire['seed'][0][0] == 0
    assert MinimalPolynomialCertificate.from_json(json.loads(json.dumps(wire))) == certificate


@pytest.mark.parametrize('change', [
    {'presentation_hash': 'short'}, {'basis_words': [(), (1,)]},
    {'basis_words': ((), ())}, {'basis_words': ((), (True,))},
    {'element': ()}, {'seed': []}, {'coefficients': []}, {'coefficients': ()},
    {'element': (0, 1)}, {'coefficients': (0, 0, 1)},
    {'kernel': object()},
])
def test_malformed_certificate_model_data_is_rejected(change):
    _, _, certificate = _dual_certificate()
    with pytest.raises((TypeError, ValueError)):
        replace(certificate, **change)


def test_certificate_field_mismatch_is_rejected_including_zero():
    _, _, certificate = _dual_certificate()
    foreign = RationalFunctionField(('s',))
    for change in ({'element': (foreign.zero, foreign.one)},
                   {'coefficients': (foreign.zero, foreign.zero, foreign.one)}):
        with pytest.raises(TypeError):
            replace(certificate, **change)


def test_well_formed_certificate_in_another_field_cannot_change_the_algebra_field():
    algebra, _, certificate = _dual_certificate()
    foreign = RationalFunctionField((), 7)
    def convert(values):
        return tuple(foreign.scalar(Fraction(dict(value.numerator).get((), 0)) /
                                    Fraction(dict(value.denominator)[()])) for value in values)
    matrix = certificate.kernel.kernel
    foreign_kernel = ExactMatrix(foreign, matrix.nrows, matrix.ncols,
        tuple((i, j, convert((value,))[0]) for i, j, value in matrix.entries))
    # The wire shape is valid. Its asserted presentation hash and coordinate
    # field disagree, and that discrepancy must be checked before arithmetic.
    forged = replace(certificate, element=convert(certificate.element),
                     seed=convert(certificate.seed), coefficients=convert(certificate.coefficients),
                     kernel=replace(certificate.kernel, kernel=foreign_kernel))
    with pytest.raises(TypeError):
        verify_minimal_polynomial(algebra, forged)


@pytest.mark.parametrize('mutation', [
    lambda wire: wire.update(extra=True),
    lambda wire: wire.pop('seed'),
    lambda wire: wire.update(schema='unknown'),
    lambda wire: wire.update(basis_words='not an array'),
    lambda wire: wire.update(basis_words=[None]),
    lambda wire: wire.update(element={}),
    lambda wire: wire.update(seed={}),
    lambda wire: wire.update(coefficients={}),
    lambda wire: wire.update(element=[[True, wire['element'][0][1]]]),
    lambda wire: wire.update(element=[[-1, wire['element'][0][1]]]),
    lambda wire: wire.update(element=[[2, wire['element'][0][1]]]),
    lambda wire: wire.update(element=[wire['element'][0], wire['element'][0]]),
    lambda wire: wire.update(element=[None]),
    lambda wire: wire.update(element=[[1]]),
    lambda wire: wire.update(element=[[1, 1.25]]),
    lambda wire: wire.update(coefficients=[]),
    lambda wire: wire['kernel'].update(extra=True),
])
def test_strict_json_shape_scalar_and_coordinate_validation(mutation):
    _, _, certificate = _dual_certificate()
    wire = deepcopy(certificate.to_json())
    mutation(wire)
    with pytest.raises((TypeError, ValueError)):
        MinimalPolynomialCertificate.from_json(wire)


def test_json_reader_rejects_nonobject_and_foreign_field_coefficients():
    with pytest.raises(ValueError):
        MinimalPolynomialCertificate.from_json([])
    _, _, certificate = _dual_certificate()
    wire = certificate.to_json()
    wire['element'][0][1] = RationalFunctionField(('s',)).zero.to_json()
    with pytest.raises(TypeError):
        MinimalPolynomialCertificate.from_json(wire)


def test_zero_orbit_roundtrip_keeps_the_declared_field_in_the_kernel():
    algebra = _zero_algebra()
    _, certificate = _certificate(algebra, (), (), (1,), ((),), ())
    wire = certificate.to_json()
    assert wire['element'] == wire['seed'] == []
    restored = MinimalPolynomialCertificate.from_json(wire)
    assert restored == certificate
    assert verify_minimal_polynomial(algebra, restored) is restored


def test_minimal_polynomial_certificates_have_no_generated_equality_hash():
    assert MinimalPolynomialCertificate.__hash__ is None
