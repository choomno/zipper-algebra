"""Small actual FiniteFlow reconstructions with independent exact verification."""
from fractions import Fraction

import pytest

from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.matrices import ExactMatrix
from zipper_algebra.reconstruction import reconstruct_kernel
from zipper_algebra.verification.linear import verify_kernel


@pytest.fixture
def native_finiteflow():
    pytest.importorskip('sympy', reason='reconstruction uses the optional symbolic extra')
    return pytest.importorskip('fflow', reason='FiniteFlow native Python bindings are optional')


def assert_exact_nullspace(original, certificate, rank):
    assert verify_kernel(original, certificate) is certificate
    assert len(certificate.pivot_columns) == rank
    assert certificate.kernel.nrows == original.ncols
    assert certificate.kernel.ncols == original.ncols - rank
    # This deliberately uses public literal coefficient arithmetic, separately
    # from the verifier's normalized working coefficients and native solver.
    a = {(i, j): value for i, j, value in original.entries}
    n = {(i, j): value for i, j, value in certificate.kernel.entries}
    for row in range(original.nrows):
        for column in range(certificate.kernel.ncols):
            residual = sum((a.get((row, k), original.field.zero)
                            * n.get((k, column), original.field.zero)
                            for k in range(original.ncols)), original.field.zero)
            assert residual.is_zero()


@pytest.mark.parametrize('threads', [1, 2])
def test_native_parameterized_sparse_reconstruction_and_repeat_cleanup(native_finiteflow, threads):
    field = RationalFunctionField(('q', 't'))
    q, t = (field.parameter(name) for name in field.parameters)
    original = ExactMatrix.from_rows(field, [[1, q / (t + 1), t],
                                            [2, 2 * q / (t + 1), 2 * t]])
    before = original.to_json()
    for _ in range(2):
        certificate = reconstruct_kernel(original, max_primes=3, max_degree=8, threads=threads)
        assert_exact_nullspace(original, certificate, 1)
    assert original.to_json() == before


def test_native_parameter_names_are_data_and_keep_order(native_finiteflow):
    field = RationalFunctionField(('lambda', 'I'))
    lam, imaginary_name = (field.parameter(name) for name in field.parameters)
    original = ExactMatrix.from_rows(field, [[lam + 1, imaginary_name / 2, Fraction(2, 3)]])
    certificate = reconstruct_kernel(original, max_primes=3, max_degree=6, threads=1)
    assert_exact_nullspace(original, certificate, 1)
    assert certificate.kernel.field.parameters == ('lambda', 'I')


@pytest.mark.parametrize('rows,rank', [([[1, 2], [3, 4]], 2),
                                      ([[0, 0], [0, 0]], 0),
                                      ([[1, 2, 3], [2, 4, 6]], 1)])
def test_native_constant_systems_full_zero_and_deficient_rank(native_finiteflow, rows, rank):
    original = ExactMatrix.from_rows(RationalFunctionField(), rows)
    certificate = reconstruct_kernel(original, max_primes=3, max_degree=4, threads=1)
    assert_exact_nullspace(original, certificate, rank)


def test_native_custom_shift_and_start_prime_recover_after_invalid_shift(native_finiteflow):
    field = RationalFunctionField(('q', 't'))
    q, t = (field.parameter(name) for name in field.parameters)
    original = ExactMatrix.from_rows(field, [[1, (q + 1) / (t + 1), q * t]])
    smallest = native_finiteflow.PrimeNo(native_finiteflow.NAvailablePrimes() - 1)
    with pytest.raises(ValueError, match='smallest'):
        reconstruct_kernel(original, shift=(smallest, 3), max_primes=3, max_degree=8)
    proof = reconstruct_kernel(original, shift=[2, 3], start_prime=1,
                               max_primes=3, max_degree=8)
    assert_exact_nullspace(original, proof, 1)
    # Graph-local controls must not change a subsequent call's defaults.
    default_proof = reconstruct_kernel(original, max_primes=3, max_degree=8)
    assert_exact_nullspace(original, default_proof, 1)


def test_native_explicit_start_prime_avoids_rank_drop_at_first_learning_prime(native_finiteflow):
    # The generic determinant is the first native prime. Learning in that prime
    # would erase the first column; the explicit next index restores full rank.
    first_prime = native_finiteflow.PrimeNo(0)
    original = ExactMatrix.from_rows(RationalFunctionField(), [[first_prime, 0], [0, 1]])
    proof = reconstruct_kernel(original, start_prime=1, max_primes=3, max_degree=4)
    assert_exact_nullspace(original, proof, 2)


def test_native_coefficient_denominator_bad_prime_then_numeric_recovery(native_finiteflow):
    first_prime = native_finiteflow.PrimeNo(0)
    original = ExactMatrix.from_rows(RationalFunctionField(), [[Fraction(1, first_prime), 1]])
    with pytest.raises(RuntimeError, match="FiniteFlow reconstruction failed"):
        reconstruct_kernel(original, start_prime=0, max_primes=3, max_degree=4)
    proof = reconstruct_kernel(original, start_prime=1, max_primes=3, max_degree=4)
    assert_exact_nullspace(original, proof, 1)


def test_native_shift_at_parameter_pole_and_alternate_shift(native_finiteflow, record_property):
    field = RationalFunctionField(('q',))
    q = field.parameter('q')
    original = ExactMatrix.from_rows(field, [[1, 1 / (q - 2)]])
    failures = (native_finiteflow.FFlowError, native_finiteflow.Failed,
                native_finiteflow.MissingPoints, native_finiteflow.MissingPrimes)
    try:
        pole_result = reconstruct_kernel(original, shift=(2,), max_primes=3, max_degree=6)
    except failures as error:
        record_property('pole_shift_outcome', type(error).__name__)
    except RuntimeError as error:
        if (str(error) != 'FiniteFlow did not reconstruct the complete kernel coefficients'
                and not str(error).startswith('FiniteFlow reconstruction failed:')):
            raise
        record_property('pole_shift_outcome', str(error))
    else:
        assert_exact_nullspace(original, pole_result, 1)
        record_property('pole_shift_outcome', 'exactly verified')
    # A rejected graph or tolerated pole must leave a fresh graph usable.
    proof = reconstruct_kernel(original, shift=(3,), max_primes=3, max_degree=6)
    assert_exact_nullspace(original, proof, 1)


def test_two_python_callers_get_independently_verified_native_results(native_finiteflow):
    from concurrent.futures import ThreadPoolExecutor
    field = RationalFunctionField(('q',))
    q = field.parameter('q')
    originals = [ExactMatrix.from_rows(field, [[1, q + 1, 2]]),
                 ExactMatrix.from_rows(field, [[1, 1 / (q + 2), 3]])]
    with ThreadPoolExecutor(max_workers=2) as callers:
        futures = [callers.submit(reconstruct_kernel, original, max_primes=3,
                                   max_degree=6, threads=1) for original in originals]
        proofs = [future.result(timeout=20) for future in futures]
    for original, proof in zip(originals, proofs):
        assert_exact_nullspace(original, proof, 1)


def test_native_signed_ground_rationals_with_explicit_empty_shift(native_finiteflow):
    original = ExactMatrix.from_rows(RationalFunctionField(),
                                     [[Fraction(-1, 2), Fraction(3, 7)],
                                      [-1, Fraction(6, 7)]])
    proof = reconstruct_kernel(original, shift=(), max_primes=3, max_degree=4)
    assert_exact_nullspace(original, proof, 1)
