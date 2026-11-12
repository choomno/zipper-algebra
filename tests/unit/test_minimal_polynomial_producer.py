"""Exact hand-specified orbit answers and independent producer-boundary checks."""
from dataclasses import replace
from fractions import Fraction

import pytest

from zipper_algebra import invariants
from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import FieldElement, RationalFunctionField
from zipper_algebra.invariants import minimal_polynomial
from zipper_algebra.matrices import ExactMatrix
from zipper_algebra.verification.derivations import DerivationGraph, DerivationTerm, VerificationLimit
from zipper_algebra.verification.finite import GeneratorAction, WordBasisCertificate, verify_word_basis
from zipper_algebra.verification.linear import KernelCertificate, _evaluate, verify_kernel
from tests.unit.test_finite_algebra import matrix_algebra_certificate, quadratic_certificate


def truncated_algebra(degree=4, characteristic=0):
    """The literal basis 1,x,...,x^(degree-1) in K[x]/(x^degree)."""
    field = RationalFunctionField((), characteristic)
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    presentation = AlgebraPresentation(parent, (x**degree,), ('nilpotence',))
    proof = DerivationGraph(presentation.content_hash(), (parent.zero, x**degree),
                            ((), (DerivationTerm(field.one, (), 'input', 0, ()),)))
    action = GeneratorAction(field, tuple(((index + 1, field.one),) if index + 1 < degree else ()
                                         for index in range(degree)))
    certificate = WordBasisCertificate(presentation.content_hash(),
        tuple((1,) * index for index in range(degree)), (action,), proof,
        (tuple(0 if index + 1 < degree else 1 for index in range(degree)),))
    algebra = verify_word_basis(presentation, certificate)
    return algebra, algebra.element((0, 1) + (0,) * (degree - 2))


@pytest.fixture
def nilpotent():
    return truncated_algebra()


def exact_kernel(matrix, **options):
    """Independent SymPy elimination for tiny controls, never the production backend."""
    sp = pytest.importorskip('sympy')
    parameters = sp.symbols(matrix.field.parameters)

    def expression(value):
        def polynomial(terms):
            return sum(sp.Rational(coefficient.numerator, coefficient.denominator)
                       * sp.prod(variable**power for variable, power in zip(parameters, monomial))
                       for monomial, coefficient in terms)
        return polynomial(value.numerator) / polynomial(value.denominator)

    def element(value):
        numerator, denominator = sp.cancel(value).as_numer_denom()
        def terms(polynomial):
            if not parameters:
                return {(): Fraction(polynomial)}
            return {monomial: Fraction(coefficient)
                    for monomial, coefficient in sp.Poly(polynomial, *parameters).terms()}
        return FieldElement(matrix.field, terms(numerator), terms(denominator))

    original = sp.zeros(matrix.nrows, matrix.ncols)
    for row, column, value in matrix.entries:
        original[row, column] = expression(value)
    reduced, pivots = original.rref()
    free = tuple(index for index in range(matrix.ncols) if index not in pivots)
    entries = [(row, column, element(-reduced[index, free_column]))
               for index, row in enumerate(pivots) for column, free_column in enumerate(free)]
    entries.extend((row, column, matrix.field.one) for column, row in enumerate(free))
    kernel = ExactMatrix(matrix.field, matrix.ncols, len(free), tuple(entries))
    rows = original[:, list(pivots)].T.rref()[1] if pivots else ()
    prime, point = 101, tuple(11 + 2 * index for index in range(len(parameters)))
    values = {(i, j): _evaluate(value, point, prime) for i, j, value in matrix.entries}
    minor = sp.Matrix([[values.get((i, j), 0) for j in pivots] for i in rows])
    inverse = (tuple(tuple(int(value) for value in row) for row in minor.inv_mod(prime).tolist())
               if pivots else ())
    result = KernelCertificate(matrix.content_hash(), kernel, free, rows, pivots, prime, point, inverse)
    return verify_kernel(matrix, result)


@pytest.fixture
def oracle(monkeypatch):
    pytest.importorskip('sympy')
    calls = []
    def reconstruct(matrix, **options):
        calls.append((matrix, options))
        return exact_kernel(matrix, **options)
    monkeypatch.setattr(invariants, 'reconstruct_kernel', reconstruct)
    return calls


def test_element_nilpotence_has_exact_degree_four_and_geometric_prefixes(nilpotent, oracle):
    algebra, x = nilpotent
    certificate = minimal_polynomial(algebra, x)
    assert certificate.coefficients == (algebra.field.zero,) * 4 + (algebra.field.one,)
    assert certificate.kernel.rank == 4 and certificate.kernel.nullity == 1
    assert [matrix.ncols for matrix, _ in oracle] == [2, 3, 5]
    assert certificate.element == x.coordinates
    assert certificate.seed == algebra.one.coordinates


def test_relative_orbit_truncates_to_first_dependency_and_differs_from_element(nilpotent, oracle):
    algebra, x = nilpotent
    relative = minimal_polynomial(algebra, x, seed=x)
    assert relative.coefficients == (algebra.field.zero,) * 3 + (algebra.field.one,)
    assert relative.kernel.rank == 3 and relative.kernel.nullity == 1
    assert [matrix.ncols for matrix, _ in oracle] == [2, 3, 5, 4]
    # Degree three annihilates x, although the same polynomial does not annihilate 1.
    assert algebra.multiply(x, algebra.multiply(x, x)) != algebra.zero


def test_noncommuting_matrix_unit_control_has_mixed_polynomial(oracle):
    pres, witness, x, y = matrix_algebra_certificate()
    algebra = verify_word_basis(pres, witness)
    a = algebra.element(algebra.coordinates(x + y))
    result = minimal_polynomial(algebra, a)
    assert result.coefficients == tuple(algebra.field.scalar(value) for value in (-1, 0, 1))
    assert algebra.multiply(algebra.element(algebra.coordinates(x)),
                            algebra.element(algebra.coordinates(y))) != algebra.multiply(
                                algebra.element(algebra.coordinates(y)),
                                algebra.element(algebra.coordinates(x)))


def test_generic_parameterized_answer_and_native_controls_are_preserved(oracle):
    pres, witness, x, s = quadratic_certificate()
    algebra = verify_word_basis(pres, witness)
    element = algebra.element(algebra.coordinates(x))
    result = minimal_polynomial(algebra, element, max_degree=9, max_operations=1234,
        max_products=5678, max_primes=3, reconstruction_degree=8, threads=2,
        start_prime=1, shift=[7])
    assert result.coefficients == (-s, algebra.field.zero, algebra.field.one)
    assert result.coefficients[-1].to_json() == algebra.field.one.to_json()
    assert [matrix.ncols for matrix, _ in oracle] == [2, 3]
    assert all(options == {'max_products': 5678, 'max_primes': 3, 'max_degree': 8,
                           'threads': 2, 'start_prime': 1, 'shift': (7,)}
               for _, options in oracle)


@pytest.mark.parametrize('degree', [0, 1, 2, 3])
def test_exhausted_degree_budget_has_no_certificate(nilpotent, oracle, degree):
    algebra, x = nilpotent
    with pytest.raises(VerificationLimit, match='degree budget'):
        minimal_polynomial(algebra, x, max_degree=degree)
    assert all(matrix.ncols <= degree + 1 for matrix, _ in oracle)


@pytest.mark.parametrize('characteristic', [0, 2, 7])
def test_zero_seed_needs_no_native_backend(monkeypatch, characteristic):
    algebra, x = truncated_algebra(characteristic=characteristic)
    def forbidden(*args, **kwargs):
        raise AssertionError('zero-orbit proof must not call the reconstruction backend')
    monkeypatch.setattr(invariants, 'reconstruct_kernel', forbidden)
    result = minimal_polynomial(algebra, x, seed=algebra.zero, max_degree=0, max_products=0)
    assert result.coefficients == (algebra.field.one,)
    assert result.kernel.rank == 0 and result.kernel.nullity == 1
    assert result.kernel.prime == (characteristic or 2147483647)


def test_zero_algebra_has_minimal_polynomial_one_without_native_backend(monkeypatch):
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    pres = AlgebraPresentation(parent, (parent.one,), ('unit',))
    graph = DerivationGraph(pres.content_hash(), (parent.one,),
                            ((DerivationTerm(field.one, (), 'input', 0, ()),),))
    cert = WordBasisCertificate(pres.content_hash(), (), (GeneratorAction(field, ()),),
                                graph, ((),), unit_node=0)
    algebra = verify_word_basis(pres, cert)
    monkeypatch.setattr(invariants, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    result = minimal_polynomial(algebra, algebra.one)
    assert result.element == result.seed == ()
    assert result.coefficients == (field.one,)
    assert result.kernel.kernel.nrows == 1


def test_nonzero_orbit_in_positive_characteristic_is_explicitly_unsupported(nilpotent):
    algebra, x = truncated_algebra(characteristic=2)
    with pytest.raises(ValueError, match='characteristic zero'):
        minimal_polynomial(algebra, x)


@pytest.mark.parametrize('kwargs', [
    {'max_degree': -1}, {'max_degree': True}, {'max_degree': 1.5},
    {'max_operations': 0}, {'max_operations': True},
    {'max_products': -1}, {'max_products': True},
    {'max_primes': 0}, {'max_primes': True}, {'max_primes': 2048},
    {'reconstruction_degree': 0}, {'reconstruction_degree': 65536},
    {'threads': 0}, {'threads': 2**32}, {'start_prime': -1},
    {'start_prime': 2047, 'max_primes': 1}, {'shift': [1]},
    {'shift': [-1]}, {'shift': [2**64]},
])
def test_invalid_options_fail_before_native_calls_even_for_zero_seed(nilpotent, monkeypatch, kwargs):
    algebra, x = nilpotent
    monkeypatch.setattr(invariants, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises(ValueError):
        minimal_polynomial(algebra, x, seed=algebra.zero, **kwargs)


def test_wrong_parent_and_unverified_inputs_are_rejected(nilpotent):
    algebra, x = nilpotent
    other, other_x = truncated_algebra()
    for invalid_algebra, invalid_element, seed in [
        (None, x, None), (algebra, 1, None), (algebra, other_x, None),
        (algebra, x, other.one), (algebra, x, 0),
    ]:
        with pytest.raises(TypeError):
            minimal_polynomial(invalid_algebra, invalid_element, seed=seed)


def test_backend_failure_remains_inconclusive_and_input_is_unchanged(nilpotent, monkeypatch):
    algebra, x = nilpotent
    before = x.coordinates
    def unavailable(*args, **kwargs):
        raise RuntimeError('bounded native search did not reconstruct')
    monkeypatch.setattr(invariants, 'reconstruct_kernel', unavailable)
    with pytest.raises(RuntimeError, match='bounded native'):
        minimal_polynomial(algebra, x)
    assert x.coordinates == before


def test_real_multiplication_operation_limit_is_preserved(nilpotent, monkeypatch):
    algebra, x = nilpotent
    monkeypatch.setattr(invariants, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises(VerificationLimit, match='operation budget'):
        minimal_polynomial(algebra, x, max_operations=1)


def test_wrong_reconstructed_coefficients_cannot_escape_independent_verifier(nilpotent, oracle, monkeypatch):
    algebra, x = nilpotent
    def corrupted(matrix, **options):
        result = exact_kernel(matrix, **options)
        if result.nullity:
            values = list(result.kernel.entries) + [(0, 0, algebra.field.one)]
            return replace(result, kernel=ExactMatrix(matrix.field, matrix.ncols, result.nullity,
                                                       tuple(values)))
        return result
    monkeypatch.setattr(invariants, 'reconstruct_kernel', corrupted)
    with pytest.raises(ValueError):
        minimal_polynomial(algebra, x)


def test_stale_zero_system_answer_cannot_certify_a_nonzero_orbit(nilpotent, oracle, monkeypatch):
    algebra, x = nilpotent
    answers = 0
    def stale_first_answer(matrix, **options):
        nonlocal answers
        answers += 1
        # This is an honestly verified kernel for a different matrix. Inject
        # reuse of that old answer across the producer boundary, without
        # manufacturing malformed certificate internals.
        if answers == 1:
            return exact_kernel(ExactMatrix(matrix.field, matrix.nrows, matrix.ncols))
        return exact_kernel(matrix)
    monkeypatch.setattr(invariants, 'reconstruct_kernel', stale_first_answer)
    with pytest.raises(RuntimeError, match='unique orbit relation'):
        minimal_polynomial(algebra, x)
    assert answers == 2


def test_stale_relation_with_zero_last_coordinate_cannot_be_normalized(nilpotent, oracle, monkeypatch):
    algebra, x = nilpotent
    other = ExactMatrix.from_rows(algebra.field, [[0, 1]])
    old_answer = exact_kernel(other)
    assert verify_kernel(other, old_answer) is old_answer
    monkeypatch.setattr(invariants, 'reconstruct_kernel', lambda *a, **kw: old_answer)
    with pytest.raises(RuntimeError, match='monic orbit relation'):
        minimal_polynomial(algebra, x)


@pytest.mark.parametrize('shift', [(True,), (-1,), (2**64,), (), (1, 2)])
def test_parameter_shift_validation_cannot_be_skipped_for_zero_seed(monkeypatch, shift):
    pytest.importorskip('sympy')
    pres, witness, x, _ = quadratic_certificate()
    algebra = verify_word_basis(pres, witness)
    monkeypatch.setattr(invariants, 'reconstruct_kernel', lambda *a, **kw: pytest.fail('backend'))
    with pytest.raises(ValueError, match='shift'):
        minimal_polynomial(algebra, algebra.element(algebra.coordinates(x)),
                           seed=algebra.zero, shift=shift)
