"""Independent exact controls for reconstructed linear-system certificates."""
from copy import deepcopy
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction

import pytest

from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.matrices import ExactMatrix
from zipper_algebra.verification.linear import KernelCertificate, verify_kernel


def matrix(field, rows, ncols=None):
    return ExactMatrix.from_rows(field, rows, ncols=ncols)


def rank_two():
    field = RationalFunctionField()
    original = matrix(field, [[1, 2, 3, 0], [0, 1, 1, 1], [1, 3, 4, 1]])
    nullspace = matrix(field, [[-1, 2], [-1, -1], [1, 0], [0, 1]])
    proof = KernelCertificate(original.content_hash(), nullspace, (2, 3), (0, 1),
                              (0, 1), 5, (), ((1, 3), (0, 1)))
    return original, proof


def reject(original, proof, **updates):
    with pytest.raises((TypeError, ValueError, ZeroDivisionError)):
        verify_kernel(original, replace(proof, **updates))


def test_hand_rank_two_witness_proves_both_rank_bounds_and_complete_kernel():
    original, proof = rank_two()
    assert verify_kernel(original, proof) is proof
    assert proof.kernel.nrows == 4 and proof.kernel.ncols == 2
    # The free-coordinate identity block proves these two null vectors independent;
    # the invertible leading minor proves no third independent null vector exists.
    assert proof.identity_rows == (2, 3)
    assert proof.pivot_columns == (0, 1)


def test_nonzero_residual_rejected_despite_correct_free_identity_block():
    original, proof = rank_two()
    bad = matrix(original.field, [[-2, 2], [-1, -1], [1, 0], [0, 1]])
    reject(original, proof, kernel=bad)


def test_dependent_incomplete_and_wrong_shape_kernels_are_rejected():
    original, proof = rank_two()
    field = original.field
    for bad in [matrix(field, [[-1, -1], [-1, -1], [1, 1], [0, 0]]),
                matrix(field, [[-1], [-1], [1], [0]]),
                matrix(field, [[-1, 2], [-1, -1], [1, 0]]),
                matrix(field, [[-1, 2, 0], [-1, -1, 0], [1, 0, 0], [0, 1, 0]])]:
        reject(original, proof, kernel=bad)
    reject(original, proof, identity_rows=(3, 2))


def test_false_full_rank_is_rejected_even_when_upper_bound_is_vacuous():
    field = RationalFunctionField()
    original = matrix(field, [[1, 2], [2, 4]])
    proof = KernelCertificate(original.content_hash(), ExactMatrix(field, 2, 0), (),
                              (0, 1), (0, 1), 5, (), ((1, 0), (0, 1)))
    with pytest.raises(ValueError):
        verify_kernel(original, proof)


def test_undercounting_rank_is_rejected_by_exact_kernel_identity():
    field = RationalFunctionField()
    original = matrix(field, [[1, 0], [0, 1]])
    proof = KernelCertificate(original.content_hash(), matrix(field, [[0], [1]]), (1,),
                              (0,), (0,), 5, (), ((1,),))
    with pytest.raises(ValueError):
        verify_kernel(original, proof)


@pytest.mark.parametrize('updates', [
    {'matrix_hash': '0' * 64}, {'minor_inverse': ((2, 3), (0, 1))},
    {'prime': 9}, {'prime': True}, {'prime': 0}, {'prime': -5},
    {'point': (1,)}, {'pivot_rows': (0,)}, {'pivot_columns': (0,)},
    {'pivot_rows': (0, 0)}, {'pivot_columns': (0, 0)},
    {'pivot_rows': (0, 3)}, {'pivot_columns': (0, 4)},
    {'identity_rows': (2, 2)}, {'identity_rows': (2, 4)},
    {'minor_inverse': ((1,), (0,))}, {'minor_inverse': ((1, 3),)},
    {'minor_inverse': ((True, 3), (0, 1))},
])
def test_malformed_or_false_minor_witness_is_rejected(updates):
    original, proof = rank_two()
    reject(original, proof, **updates)


def test_unsorted_minor_and_kernel_identity_rows_remain_mathematically_valid():
    original, proof = rank_two()
    proof = replace(proof, pivot_rows=(1, 0), pivot_columns=(1, 0),
                    minor_inverse=((1, 0), (3, 1)))
    assert verify_kernel(original, proof) is proof
    switched = matrix(original.field, [[2, -1], [-1, -1], [0, 1], [1, 0]])
    proof = replace(proof, kernel=switched, identity_rows=(3, 2))
    assert verify_kernel(original, proof) is proof


@pytest.mark.parametrize('shape', [(0, 0), (0, 3), (2, 0), (3, 3)])
def test_zero_matrices_and_empty_dimensions_have_exact_complete_kernels(shape):
    field = RationalFunctionField()
    rows, columns = shape
    original = ExactMatrix(field, rows, columns)
    nullspace = ExactMatrix(field, columns, columns,
                            tuple((i, i, 1) for i in range(columns)))
    proof = KernelCertificate(original.content_hash(), nullspace, tuple(range(columns)),
                              (), (), 5, (), ())
    assert verify_kernel(original, proof) is proof


def test_full_rank_with_no_nullspace_columns():
    field = RationalFunctionField()
    original = matrix(field, [[2, 0], [0, 3], [2, 3]])
    proof = KernelCertificate(original.content_hash(), ExactMatrix(field, 2, 0), (),
                              (0, 1), (0, 1), 5, (), ((3, 0), (0, 2)))
    assert verify_kernel(original, proof) is proof


@pytest.mark.parametrize('characteristic', [0, 5])
def test_parameterized_rational_nullspace_and_safe_modular_minor(characteristic):
    pytest.importorskip('sympy')
    field = RationalFunctionField(('t', 'q'), characteristic)
    t, q = (field.parameter(name) for name in field.parameters)
    original = matrix(field, [[1 / (q + t), q / (q + t)]])
    nullspace = matrix(field, [[-q], [1]])
    proof = KernelCertificate(original.content_hash(), nullspace, (1,), (0,), (0,),
                              5, (1, 2), ((3,),))
    assert verify_kernel(original, proof) is proof
    reject(original, proof, point=(1, 4))
    reject(original, proof, point=(Fraction(1, 5), 2))
    reject(original, proof, point=(True, 2))
    reject(original, proof, point=(1.0, 2))
    if characteristic:
        reject(original, proof, prime=7)


def test_every_stored_input_denominator_is_checked_not_only_selected_minor():
    pytest.importorskip('sympy')
    field = RationalFunctionField(('q',))
    q = field.parameter('q')
    # The selected first column is constant, but the other column has a pole.
    original = matrix(field, [[1, 1 / (q - 1)]])
    proof = KernelCertificate(original.content_hash(), matrix(field, [[-1 / (q - 1)], [1]]),
                              (1,), (0,), (0,), 5, (1,), ((1,),))
    with pytest.raises((ValueError, ZeroDivisionError)):
        verify_kernel(original, proof)


def test_scalar_denominator_vanishing_modulo_prime_is_rejected():
    field = RationalFunctionField()
    original = matrix(field, [[1, Fraction(1, 5)]])
    proof = KernelCertificate(original.content_hash(), matrix(field, [[Fraction(-1, 5)], [1]]),
                              (1,), (0,), (0,), 5, (), ((1,),))
    with pytest.raises((ValueError, ZeroDivisionError)):
        verify_kernel(original, proof)


def test_kernel_parent_is_checked_even_if_all_entries_are_zero():
    original, proof = rank_two()
    for foreign in [RationalFunctionField((), 5), RationalFunctionField(('q',))]:
        reject(original, proof, kernel=ExactMatrix(foreign, 4, 2))
    field = RationalFunctionField(('q', 't'))
    a = ExactMatrix(field, 1, 1)
    k = ExactMatrix(RationalFunctionField(('t', 'q')), 1, 1, ((0, 0, 1),))
    candidate = KernelCertificate(a.content_hash(), k, (0,), (), (), 5, (1, 2), ())
    with pytest.raises((ValueError, TypeError)):
        verify_kernel(a, candidate)


def test_operation_budget_exhaustion_does_not_accept_the_certificate():
    original, proof = rank_two()
    with pytest.raises(RuntimeError, match='budget'):
        verify_kernel(original, proof, max_products=1)


def test_sparse_matrix_canonicalization_snapshots_and_serialization():
    field = RationalFunctionField(('q',))
    q = field.parameter('q')
    original = ExactMatrix(field, 2, 3, ((1, 2, q), (0, 1, Fraction(1, 2)), (0, 0, 0)))
    assert original.entries == ((0, 1, field.scalar(Fraction(1, 2))), (1, 2, q))
    equivalent = matrix(field, [[0, Fraction(1, 2), 0], [0, 0, q]])
    assert original == equivalent and original.content_hash() == equivalent.content_hash()
    restored = ExactMatrix.from_json(original.to_json())
    assert restored == original and restored.content_hash() == original.content_hash()
    encoded = original.to_json()
    encoded.clear()
    assert len(original.entries) == 2
    with pytest.raises((FrozenInstanceError, AttributeError)):
        original.nrows = 10
    changed = ExactMatrix(field, 2, 3, ((1, 2, q + 1), (0, 1, Fraction(1, 2))))
    assert changed.content_hash() != original.content_hash()
    assert ExactMatrix(field, 0, 3).content_hash() != ExactMatrix(field, 0, 2).content_hash()


def test_dense_matrix_input_does_not_alias_mutable_rows():
    field = RationalFunctionField()
    rows = [[1, 2], [3, 4]]
    value = ExactMatrix.from_rows(field, rows)
    rows[0][0] = 7
    assert value.entries[0] == (0, 0, field.one)
    assert ExactMatrix.from_rows(field, [], ncols=3).ncols == 3
    assert ExactMatrix.from_rows(field, []).ncols == 0
    with pytest.raises(ValueError):
        ExactMatrix.from_rows(field, [[1], [1, 2]])
    with pytest.raises(ValueError):
        ExactMatrix.from_rows(field, [[1]], ncols=2)


@pytest.mark.parametrize('shape', [(-1, 2), (2, -1), (True, 2), (2, 2.0)])
def test_matrix_dimensions_are_exact_nonnegative_integers(shape):
    with pytest.raises((TypeError, ValueError)):
        ExactMatrix(RationalFunctionField(), *shape)


@pytest.mark.parametrize('entries', [
    ((0, 0, 0), (0, 0, 1)), ((0, 0, 1), (0, 0, 0)),
    ((-1, 0, 1),), ((0, 2, 1),), ((True, 0, 1),),
    ((0, 0),), ((0, 0, 1, 2),), ((0, 0, True),),
    ((0, 0, 0.5),), ((0, 0, '1'),),
])
def test_invalid_sparse_locations_and_inexact_coefficients_are_rejected(entries):
    with pytest.raises((TypeError, ValueError)):
        ExactMatrix(RationalFunctionField(), 2, 2, entries)


def test_matrix_rejects_undeclared_field_and_foreign_zero():
    with pytest.raises(TypeError):
        ExactMatrix('Q', 1, 1)
    with pytest.raises(TypeError):
        ExactMatrix(RationalFunctionField(), 1, 1, ((0, 0, RationalFunctionField((), 5).zero),))


def test_matrix_json_rejects_unknown_missing_and_ambiguous_data():
    source = ExactMatrix(RationalFunctionField(), 1, 1, ((0, 0, 1),)).to_json()
    bad_extra = dict(source, unexpected='not accepted')
    bad_missing = dict(source)
    bad_missing.pop('field')
    bad_schema = dict(source, schema='future-or-unrelated-schema')
    for bad in (None, [], bad_extra, bad_missing, bad_schema):
        with pytest.raises((TypeError, ValueError)):
            ExactMatrix.from_json(bad)
    duplicate = deepcopy(source)
    duplicate['entries'].append(deepcopy(duplicate['entries'][0]))
    with pytest.raises((TypeError, ValueError)):
        ExactMatrix.from_json(duplicate)


@pytest.mark.parametrize('budget', [-1, True, 1.5, None])
def test_verification_budget_must_be_an_exact_nonnegative_integer(budget):
    original, proof = rank_two()
    with pytest.raises(ValueError):
        verify_kernel(original, proof, max_products=budget)


def test_verifier_rejects_objects_without_declared_matrix_or_certificate_types():
    original, proof = rank_two()
    with pytest.raises(TypeError):
        verify_kernel(original.to_json(), proof)
    with pytest.raises(TypeError):
        verify_kernel(original, proof.to_json())


def test_kernel_certificate_snapshot_json_and_immutable_state():
    original, proof = rank_two()
    mutable_inverse = [[1, 3], [0, 1]]
    snapshot = replace(proof, identity_rows=[2, 3], point=[], minor_inverse=mutable_inverse)
    mutable_inverse[0][0] = 2
    assert snapshot.minor_inverse == ((1, 3), (0, 1))
    assert snapshot.rank == 2 and snapshot.nullity == 2
    restored = KernelCertificate.from_json(snapshot.to_json())
    assert restored == snapshot
    assert verify_kernel(original, restored) is restored
    with pytest.raises((FrozenInstanceError, AttributeError)):
        proof.prime = 7
    with pytest.raises(TypeError):
        replace(proof, kernel={})
    for digest in ['short', 'G' * 64, 0]:
        with pytest.raises(ValueError):
            replace(proof, matrix_hash=digest)


@pytest.mark.parametrize('mutate', [
    lambda d: d.update(schema='future'),
    lambda d: d.update(unknown=1),
    lambda d: d.pop('kernel'),
    lambda d: d.update(identity_rows=()),
    lambda d: d.update(point=[1]),
    lambda d: d.update(point=['unparsed']),
    lambda d: d.update(minor_inverse=[(1, 3), (0, 1)]),
])
def test_kernel_certificate_json_rejects_unrecognized_or_inexact_data(mutate):
    _, proof = rank_two()
    encoded = proof.to_json()
    mutate(encoded)
    with pytest.raises((TypeError, ValueError)):
        KernelCertificate.from_json(encoded)


@pytest.mark.parametrize('mutate', [
    lambda d: d.update(orientation='row-vectors'),
    lambda d: d.update(shape=()),
    lambda d: d.update(shape=[1]),
    lambda d: d.update(entries=()),
    lambda d: d.update(entries=[[0, 0]]),
    lambda d: d.update(entries=[(0, 0, {})]),
])
def test_matrix_json_rejects_ambiguous_shape_orientation_or_entries(mutate):
    encoded = ExactMatrix(RationalFunctionField(), 1, 1).to_json()
    mutate(encoded)
    with pytest.raises((TypeError, ValueError)):
        ExactMatrix.from_json(encoded)


def test_modular_minor_evaluation_handles_huge_sparse_exponents_without_integer_growth():
    pytest.importorskip('sympy')
    field = RationalFunctionField(('q',))
    value = field.polynomial({(10**9,): 1})
    original = matrix(field, [[value]])
    proof = KernelCertificate(original.content_hash(), ExactMatrix(field, 1, 0), (),
                              (0,), (0,), 5, (2,), ((1,),))
    assert verify_kernel(original, proof) is proof


@pytest.mark.parametrize('kwargs', [
    {'max_primes': 0}, {'max_primes': True}, {'max_degree': -1}, {'max_degree': 1.0},
    {'threads': 0}, {'threads': '1'}, {'max_products': -1}, {'max_products': False},
])
def test_producer_rejects_invalid_limits_before_loading_native_bindings(kwargs):
    from zipper_algebra.reconstruction import reconstruct_kernel
    with pytest.raises(ValueError):
        reconstruct_kernel(ExactMatrix(RationalFunctionField(), 0, 1), **kwargs)


def test_producer_rejects_wrong_parent_characteristic_and_native_exponent_overflow():
    from zipper_algebra.reconstruction import reconstruct_kernel
    with pytest.raises(TypeError):
        reconstruct_kernel([[1]])
    with pytest.raises(ValueError, match='characteristic'):
        reconstruct_kernel(ExactMatrix(RationalFunctionField((), 5), 0, 0))
    field = RationalFunctionField(('q',))
    for value in [field.polynomial({(2**16,): 1}), 1 / field.polynomial({(2**16,): 1})]:
        with pytest.raises(ValueError, match='exponents'):
            reconstruct_kernel(matrix(field, [[value]]))


def test_optional_binding_absence_has_an_actionable_error(monkeypatch):
    from types import SimpleNamespace
    from zipper_algebra import reconstruction
    def missing(name):
        raise ImportError('native bindings unavailable')
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=missing))
    with pytest.raises(ImportError, match='FiniteFlow Python bindings'):
        reconstruction.reconstruct_kernel(matrix(RationalFunctionField(), [[1]]))
    # Empty input has a direct exact witness and does not require loading native code.
    original = ExactMatrix(RationalFunctionField(), 0, 2)
    proof = reconstruction.reconstruct_kernel(original)
    assert verify_kernel(original, proof) is proof and proof.nullity == 2


def test_minor_witness_uses_actual_library_inverse_and_retries_poles_and_zero_minors():
    import random
    pytest.importorskip('sympy')
    from zipper_algebra.reconstruction import _lower_inverse
    field = RationalFunctionField(('q',))
    q = field.parameter('q')
    prime = 2147483647
    first = random.Random(0).randrange(prime)
    for coefficient in [1 / (q - first), q - first]:
        used_prime, point, inverse = _lower_inverse(matrix(field, [[coefficient]]), (0,), (0,))
        assert used_prime == prime and point != (first,)
        scalar = coefficient.specialize({'q': point[0]}, characteristic=prime)
        assert scalar * inverse[0][0] == 1
    used_prime, point, inverse = _lower_inverse(ExactMatrix(field, 0, 2), (), ())
    assert used_prime == prime and len(point) == 1 and inverse == ()


def test_minor_schedule_has_independent_coordinates_not_a_fixed_affine_line():
    pytest.importorskip('sympy')
    from zipper_algebra.reconstruction import _lower_inverse
    field = RationalFunctionField(('q', 't'))
    q, t = (field.parameter(name) for name in field.parameters)
    coefficient = q - t + 1  # Identically zero on the former schedule (a, a+1).
    prime, point, inverse = _lower_inverse(matrix(field, [[coefficient]]), (0,), (0,))
    assert (point[0] - point[1] + 1) % prime != 0
    scalar = coefficient.specialize(dict(zip(field.parameters, point)), characteristic=prime)
    assert scalar * inverse[0][0] == 1


@pytest.mark.parametrize('coefficient', [2147483647, Fraction(1, 2147483647)])
def test_minor_search_changes_prime_when_first_prime_erases_determinant_or_denominator(coefficient):
    pytest.importorskip('sympy')
    from zipper_algebra.reconstruction import _lower_inverse
    field = RationalFunctionField()
    prime, point, inverse = _lower_inverse(matrix(field, [[coefficient]]), (0,), (0,))
    assert prime == 2147483629 and point == ()
    assert field.scalar(coefficient).specialize({}, characteristic=prime) * inverse[0][0] == 1


def test_minor_search_bounded_failure_is_explicit_after_all_witness_primes():
    pytest.importorskip('sympy')
    from zipper_algebra.reconstruction import _lower_inverse
    erased_everywhere = 2147483647 * 2147483629 * 2147483587
    with pytest.raises(RuntimeError, match='bounded witness search'):
        _lower_inverse(matrix(RationalFunctionField(), [[erased_everywhere]]), (0,), (0,))


class NativeOutput:
    """Explicit native-wire result for fault-injection tests; not a solver."""
    def __init__(self, count=1, nvars=1, terms=None):
        self.count, self.variables = count, nvars
        self.terms = terms if terms is not None else ([('-1', (1,))], [('1', (0,))])

    def __len__(self):
        return self.count

    def nvars(self):
        return self.variables

    def monomials(self, index):
        return self.terms


class NativeProtocol:
    """Only the small graph protocol, used to reject fabricated native faults."""
    RatFunList = NativeOutput

    def __init__(self, dependent=(0,), free=(1,), rows=(0,), output=None, fail=False):
        self.dependent, self.free, self.rows = dependent, free, rows
        self.output = NativeOutput() if output is None else output
        self.deleted = []
        self.calls = []
        self.fail = fail

    def NAvailablePrimes(self):
        return 2048

    def PrimeNo(self, index):
        return 2147483647

    def SetReconstructionShift(self, graph, shift):
        self.calls.append(('shift', tuple(shift)))

    def NewGraphWithInput(self, nvars):
        self.calls.append(('graph', nvars))
        return 17, 23

    def NewIdxRatFunList(self, nvars, coefficients):
        return coefficients

    def AlgAnalyticSparseLSolve(self, graph, input_node, ncols, sparse_rows, functions):
        return 29

    def AlgNumericSparseLSolve(self, graph, ncols, sparse_rows, coefficients):
        assert all(type(value) is str for value in coefficients)
        assert sum(map(len, sparse_rows)) == len(coefficients)
        self.numeric_constructor = (ncols, sparse_rows, coefficients)
        return 29

    def LSolveValidateFullSystem(self, graph, node, enabled):
        assert enabled is True

    def LSolveOnlyHomogeneous(self, graph, node):
        pass

    def LSolveSparseOutput(self, graph, node, sparse):
        pass

    def SetOutputNode(self, graph, node):
        pass

    def LearnEx(self, graph, *, prime_no):
        self.calls.append(('learn', prime_no))
        if self.fail:
            raise ArithmeticError('native learning failed')

    def LSolveNDepVars(self, graph, node):
        return len(self.dependent)

    def LSolveNIndepVars(self, graph, node):
        return len(self.free)

    def LSolveNIndepEqs(self, graph, node):
        return len(self.rows)

    def LSolveDepVars(self, graph, node):
        return self.dependent

    def LSolveIndepVars(self, graph, node):
        return self.free

    def LSolveIndepEqs(self, graph, node):
        return self.rows

    def ReconstructFunction(self, graph, **kwargs):
        self.calls.append(('function', dict(kwargs)))
        return self.output

    def DeleteGraph(self, graph):
        self.deleted.append(graph)


def invoke_protocol(monkeypatch, native, original=None, **options):
    from types import SimpleNamespace
    from zipper_algebra import reconstruction
    field = RationalFunctionField(('q',))
    if original is None:
        original = matrix(field, [[1, field.parameter('q')]])
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=lambda _: native))
    return reconstruction.reconstruct_kernel(original, **options)


@pytest.mark.parametrize('metadata', [
    {'dependent': (0,), 'free': (0,)},
    {'dependent': (0,), 'free': (), 'rows': (0,)},
    {'dependent': (True,), 'free': (1,)},
    {'dependent': (0,), 'free': (1,), 'rows': ()},
    {'dependent': (0,), 'free': (1,), 'rows': (1,)},
    {'dependent': (0, 1), 'free': (), 'rows': (0, 0)},
])
def test_native_inconsistent_metadata_is_rejected_and_graph_released(monkeypatch, metadata):
    native = NativeProtocol(**metadata)
    with pytest.raises(RuntimeError, match='indices'):
        invoke_protocol(monkeypatch, native)
    assert native.deleted == [17]


@pytest.mark.parametrize('output', [object(), NativeOutput(count=0), NativeOutput(nvars=2),
                                    NativeOutput(terms=([('-1', (1,)), ('-1', (1,))], [('1', (0,))])),
                                    NativeOutput(terms=([('-1', (1,))], [('1', (0,)), ('1', (0,))]))])
def test_native_incomplete_output_and_repeated_monomials_are_rejected(monkeypatch, output):
    native = NativeProtocol(output=output)
    with pytest.raises(RuntimeError):
        invoke_protocol(monkeypatch, native)
    assert native.deleted == [17]


def test_native_failure_propagates_and_releases_graph(monkeypatch):
    native = NativeProtocol(fail=True)
    with pytest.raises(ArithmeticError, match='native learning failed'):
        invoke_protocol(monkeypatch, native)
    assert native.deleted == [17]


def test_false_reconstructed_formula_is_rejected_by_independent_exact_verifier(monkeypatch):
    pytest.importorskip('sympy')
    native = NativeProtocol(output=NativeOutput(terms=([('-2', (1,))], [('1', (0,))])))
    with pytest.raises(RuntimeError, match=r'independent exact verification:.*A\*N'):
        invoke_protocol(monkeypatch, native)
    assert native.deleted == [17]


def test_native_rank_zero_claim_on_nonzero_system_is_rejected(monkeypatch):
    pytest.importorskip('sympy')
    native = NativeProtocol(dependent=(), free=(0, 1), rows=())
    with pytest.raises(RuntimeError, match=r'independent exact verification:.*A\*N'):
        invoke_protocol(monkeypatch, native)
    assert native.deleted == [17]


def test_upstream_binding_without_full_system_validation_is_rejected_before_graph_creation(monkeypatch):
    native = NativeProtocol()
    monkeypatch.delattr(NativeProtocol, 'LSolveValidateFullSystem')
    with pytest.raises(ImportError, match='fork'):
        invoke_protocol(monkeypatch, native)
    assert native.deleted == []


@pytest.mark.parametrize('point', ['1.5', '1e3', 'NaN', '1/2/3', ' 1 ', 'q'])
def test_kernel_certificate_json_accepts_only_exact_rational_notation(point):
    _, proof = rank_two()
    data = proof.to_json()
    data['point'] = [point]
    with pytest.raises(ValueError):
        KernelCertificate.from_json(data)


@pytest.mark.parametrize('kwargs', [{'max_degree': 0}, {'max_degree': 2**16},
                                    {'max_primes': 2048}, {'threads': 2**32}])
def test_native_option_boundary_rejects_truncation_and_silent_default_substitution(kwargs):
    from zipper_algebra.reconstruction import reconstruct_kernel
    # A zero-size problem exercises validation without constructing a native graph.
    with pytest.raises(ValueError):
        reconstruct_kernel(ExactMatrix(RationalFunctionField(), 0, 0), **kwargs)


@pytest.mark.parametrize('shape', [(2**32, 0), (0, 2**32), (2**64, 2)])
def test_native_dimension_overflow_is_rejected_without_dense_allocation(shape):
    from zipper_algebra.reconstruction import reconstruct_kernel
    original = ExactMatrix(RationalFunctionField(), *shape)
    with pytest.raises(ValueError, match='dimensions'):
        reconstruct_kernel(original)


def test_maximum_native_options_and_large_row_dimension_are_valid_without_allocation():
    from zipper_algebra.reconstruction import reconstruct_kernel
    original = ExactMatrix(RationalFunctionField(), 2**32 - 1, 0)
    proof = reconstruct_kernel(original, max_primes=2047, max_degree=2**16 - 1,
                               threads=2**32 - 1, max_products=0)
    assert verify_kernel(original, proof, max_products=0) is proof
    assert proof.rank == proof.nullity == 0


@pytest.mark.parametrize('bad_row', [0, 9])
def test_streamed_kernel_residual_rejects_early_and_late_rows(bad_row):
    field = RationalFunctionField()
    rows = [[1, 2] for _ in range(10)]
    rows[bad_row] = [1, 3]
    original = matrix(field, rows)
    proof = KernelCertificate(original.content_hash(), matrix(field, [[-2], [1]]),
                              (1,), (int(bad_row == 0),), (0,), 5, (), ((1,),))
    with pytest.raises(ValueError, match=r'A\*N'):
        verify_kernel(original, proof)


@pytest.mark.parametrize('output', [None, (), ['-2', '0'], [], object()])
def test_numeric_native_output_wrong_type_or_cardinality_is_rejected_with_cleanup(monkeypatch, output):
    native = NativeProtocol()
    native.ReconstructNumeric = lambda *args, **kwargs: output
    original = matrix(RationalFunctionField(), [[1, 2]])
    with pytest.raises(RuntimeError, match='numeric kernel'):
        invoke_protocol(monkeypatch, native, original)
    assert native.deleted == [17]


@pytest.mark.parametrize('output', [['not-a-rational'], ['1/0']])
def test_numeric_native_malformed_rational_is_rejected_with_cleanup(monkeypatch, output):
    native = NativeProtocol()
    native.ReconstructNumeric = lambda *args, **kwargs: output
    with pytest.raises(RuntimeError, match='FiniteFlow reconstruction failed'):
        invoke_protocol(monkeypatch, native, matrix(RationalFunctionField(), [[1, 2]]))
    assert native.deleted == [17]


def test_numeric_native_failure_propagates_and_releases_graph(monkeypatch):
    native = NativeProtocol()
    def fail(*args, **kwargs):
        raise ArithmeticError('native numeric reconstruction failed')
    native.ReconstructNumeric = fail
    with pytest.raises(ArithmeticError, match='numeric reconstruction failed'):
        invoke_protocol(monkeypatch, native, matrix(RationalFunctionField(), [[1, 2]]))
    assert native.deleted == [17]


def test_numeric_native_false_rational_formula_is_exactly_rejected(monkeypatch):
    pytest.importorskip('sympy')
    native = NativeProtocol()
    native.ReconstructNumeric = lambda *args, **kwargs: ['-3']
    with pytest.raises(RuntimeError, match=r'independent exact verification:.*A\*N'):
        invoke_protocol(monkeypatch, native, matrix(RationalFunctionField(), [[1, 2]]))
    assert native.deleted == [17]


@pytest.mark.parametrize('start', [-1, True, 1.5, 2047])
def test_start_prime_requires_a_valid_table_index_with_unused_checking_prime(start):
    from zipper_algebra.reconstruction import reconstruct_kernel
    with pytest.raises(ValueError):
        reconstruct_kernel(ExactMatrix(RationalFunctionField(), 0, 0),
                           start_prime=start, max_primes=1)


@pytest.mark.parametrize('shift', [(), (1, 2), (True,), (-1,), (2**64,), (1.5,), 1])
def test_shift_arity_and_exact_uint64_coordinates_are_checked_before_native_import(shift):
    from zipper_algebra.reconstruction import reconstruct_kernel
    field = RationalFunctionField(('q',))
    with pytest.raises((TypeError, ValueError)):
        reconstruct_kernel(matrix(field, [[1]]), shift=shift)


def test_actual_native_prime_table_budget_is_checked_before_graph_allocation(monkeypatch):
    native = NativeProtocol()
    native.NAvailablePrimes = lambda: 6
    with pytest.raises(ValueError, match='checking prime'):
        invoke_protocol(monkeypatch, native, start_prime=1, max_primes=5)
    assert native.calls == [] and native.deleted == []


@pytest.mark.parametrize('shift', [(7,), (8,)])
def test_shift_must_be_strictly_below_actual_smallest_native_prime(monkeypatch, shift):
    native = NativeProtocol()
    native.PrimeNo = lambda index: 7
    with pytest.raises(ValueError, match='smallest'):
        invoke_protocol(monkeypatch, native, shift=shift)
    assert native.calls == [] and native.deleted == []


def test_custom_shift_requires_fork_shift_setter_before_graph_allocation(monkeypatch):
    native = NativeProtocol()
    monkeypatch.delattr(NativeProtocol, 'SetReconstructionShift')
    with pytest.raises(ImportError, match='SetReconstructionShift'):
        invoke_protocol(monkeypatch, native, shift=(3,))
    assert native.calls == [] and native.deleted == []


def test_custom_shift_and_start_prime_are_forwarded_before_parameter_learning(monkeypatch):
    pytest.importorskip('sympy')
    native = NativeProtocol(output=NativeOutput(terms=([('-2', (1,))], [('1', (0,))])))
    shift = [3]
    # A deliberately false output keeps this a fault-injection/configuration test;
    # actual successful reconstruction is covered by the native integration tests.
    with pytest.raises(RuntimeError, match=r'independent exact verification:.*A\*N'):
        invoke_protocol(monkeypatch, native, start_prime=2, max_primes=3, shift=shift)
    shift[0] = 4
    assert native.calls[:3] == [('graph', 1), ('shift', (3,)), ('learn', 2)]
    options = native.calls[3][1]
    assert native.calls[3][0] == 'function'
    assert options['start_mod'] == 2 and options['max_primes'] == 3
    assert native.deleted == [17]


def test_start_prime_is_also_forwarded_to_numeric_reconstruction(monkeypatch):
    pytest.importorskip('sympy')
    native = NativeProtocol()
    def numeric(graph, **kwargs):
        native.calls.append(('numeric', kwargs))
        return ['-3']
    native.ReconstructNumeric = numeric
    with pytest.raises(RuntimeError, match=r'independent exact verification:.*A\*N'):
        invoke_protocol(monkeypatch, native, matrix(RationalFunctionField(), [[1, 2]]),
                        start_prime=1, max_primes=3)
    assert native.calls[1] == ('learn', 1)
    assert native.calls[2][0] == 'numeric'
    assert native.calls[2][1]['start_mod'] == 1
    assert native.deleted == [17]


def test_ground_verification_and_trivial_kernel_reconstruction_need_no_symbolic_dependency(monkeypatch):
    import builtins
    from zipper_algebra.reconstruction import reconstruct_kernel
    actual_import = builtins.__import__
    def without_sympy(name, *args, **kwargs):
        if name == 'sympy' or name.startswith('sympy.'):
            raise ImportError('symbolic dependency deliberately unavailable')
        return actual_import(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', without_sympy)
    original, proof = rank_two()
    assert verify_kernel(original, proof) is proof
    reject(original, proof, minor_inverse=((2, 3), (0, 1)))
    zero = ExactMatrix(RationalFunctionField(), 0, 2)
    trivial = reconstruct_kernel(zero)
    assert verify_kernel(zero, trivial) is trivial
    assert trivial.rank == 0 and trivial.nullity == 2
    large = replace(proof, prime=9223372036854775783)
    with pytest.raises(ImportError, match='symbolic extra'):
        verify_kernel(original, large)


def test_large_prime_witness_uses_definitive_optional_primality_check():
    pytest.importorskip('sympy')
    field = RationalFunctionField()
    original = matrix(field, [[1]])
    proof = KernelCertificate(original.content_hash(), ExactMatrix(field, 1, 0), (),
                              (0,), (0,), 9223372036854775783, (), ((1,),))
    assert verify_kernel(original, proof) is proof
    reject(original, proof, prime=2**63)
    reject(original, proof, prime=2**64)


def test_python_callers_serialize_native_graphs_and_release_lock_after_failure(monkeypatch):
    import threading
    from types import SimpleNamespace
    from zipper_algebra import reconstruction
    first_learning = threading.Event()
    second_lock_attempt = threading.Event()
    release_first = threading.Event()
    overlap = threading.Event()
    state_lock = threading.Lock()
    original_lock = reconstruction._NATIVE_LOCK

    class ObservedLock:
        attempts = 0

        def __enter__(self):
            with state_lock:
                self.attempts += 1
                if self.attempts == 2:
                    second_lock_attempt.set()
            return original_lock.__enter__()

        def __exit__(self, *exception):
            return original_lock.__exit__(*exception)

    class OverlapProtocol(NativeProtocol):
        active_graphs = 0
        graphs_created = 0

        def NewGraphWithInput(self, nvars):
            with state_lock:
                self.active_graphs += 1
                self.graphs_created += 1
                if self.active_graphs > 1:
                    overlap.set()
            return super().NewGraphWithInput(nvars)

        def LearnEx(self, graph, *, prime_no):
            if self.graphs_created == 1:
                first_learning.set()
                if not release_first.wait(3):
                    raise AssertionError('test did not release first native graph')
            raise ArithmeticError('deliberate per-call native failure')

        def DeleteGraph(self, graph):
            with state_lock:
                self.active_graphs -= 1
            super().DeleteGraph(graph)

    native = OverlapProtocol()
    monkeypatch.setattr(reconstruction, '_NATIVE_LOCK', ObservedLock())
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=lambda _: native))
    original = matrix(RationalFunctionField(), [[1, 2]])
    outcomes = []

    def worker():
        try:
            reconstruction.reconstruct_kernel(original)
        except Exception as error:
            outcomes.append(error)

    workers = [threading.Thread(target=worker, daemon=True) for _ in range(2)]
    workers[0].start()
    try:
        assert first_learning.wait(3)
        workers[1].start()
        assert second_lock_attempt.wait(3)
        assert native.graphs_created == 1 and native.active_graphs == 1
        assert not overlap.is_set()
    finally:
        release_first.set()
        for thread in workers:
            if thread.ident is not None:
                thread.join(3)
    assert not any(thread.is_alive() for thread in workers)
    assert len(outcomes) == 2 and all(isinstance(error, ArithmeticError) for error in outcomes)
    assert native.graphs_created == 2 and native.active_graphs == 0
    assert native.deleted == [17, 17] and not overlap.is_set()
    assert not original_lock.locked()


def test_ground_coefficients_reach_dedicated_numeric_constructor_as_exact_rationals(monkeypatch):
    pytest.importorskip('sympy')
    native = NativeProtocol()
    native.ReconstructNumeric = lambda *args, **kwargs: ['-1']
    original = matrix(RationalFunctionField(), [[Fraction(2, 3), 1]])
    with pytest.raises(RuntimeError, match=r'independent exact verification:.*A\*N'):
        invoke_protocol(monkeypatch, native, original)
    assert native.numeric_constructor == (2, [[0, 1]], ['2/3', '1'])
    assert native.deleted == [17]


def test_empty_parameterized_reconstruction_has_complete_proof_without_native_graph(monkeypatch):
    from types import SimpleNamespace
    pytest.importorskip('sympy')
    from zipper_algebra import reconstruction
    def no_native_graph(name):
        raise AssertionError('empty matrix must not load FiniteFlow')
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=no_native_graph))
    field = RationalFunctionField(('q',))
    original = ExactMatrix(field, 2, 3)
    proof = reconstruction.reconstruct_kernel(original, max_products=0)
    assert verify_kernel(original, proof, max_products=0) is proof
    assert proof.kernel == ExactMatrix(field, 3, 3, tuple((i, i, 1) for i in range(3)))
    assert proof.rank == 0 and proof.nullity == 3 and proof.point == (2,)


def test_characteristic_mismatch_precedes_optional_large_prime_dependency(monkeypatch):
    import builtins
    field = RationalFunctionField((), 5)
    original = matrix(field, [[1]])
    proof = KernelCertificate(original.content_hash(), ExactMatrix(field, 1, 0), (),
                              (0,), (0,), 9223372036854775783, (), ((1,),))
    actual_import = builtins.__import__
    def without_sympy(name, *args, **kwargs):
        if name == 'sympy' or name.startswith('sympy.'):
            raise ImportError('symbolic dependency deliberately unavailable')
        return actual_import(name, *args, **kwargs)
    monkeypatch.setattr(builtins, '__import__', without_sympy)
    with pytest.raises(ValueError, match='characteristic'):
        verify_kernel(original, proof)


@pytest.mark.parametrize('status', ['FFlowError', 'Failed', 'MissingPoints', 'MissingPrimes'])
def test_documented_native_failures_are_normalized_and_graph_released(monkeypatch, status):
    native = NativeProtocol()
    class NativeFailure(Exception):
        pass
    setattr(native, status, NativeFailure)
    def failure(*args, **kwargs):
        raise NativeFailure('bounded native failure')
    native.LearnEx = failure
    with pytest.raises(RuntimeError, match='FiniteFlow reconstruction failed') as caught:
        invoke_protocol(monkeypatch, native)
    assert isinstance(caught.value.__cause__, NativeFailure)
    assert native.deleted == [17]


def test_zero_denominator_in_serialized_witness_point_is_invalid_input():
    original, proof = rank_two()
    data = proof.to_json()
    data['point'] = ['1/0']
    with pytest.raises(ValueError, match='zero scalar denominator'):
        KernelCertificate.from_json(data)


def test_kernel_product_exhaustion_is_a_typed_verification_limit():
    from zipper_algebra.verification.derivations import VerificationLimit
    original, proof = rank_two()
    with pytest.raises(VerificationLimit, match='scalar-product budget'):
        verify_kernel(original, proof, max_products=0)


def test_algebraic_matrices_and_kernel_certificates_are_explicitly_unhashable():
    original, proof = rank_two()
    empty = ExactMatrix(RationalFunctionField(), 0, 0)
    assert ExactMatrix.__hash__ is None
    assert KernelCertificate.__hash__ is None
    for value in (original, empty, proof):
        with pytest.raises(TypeError, match='unhashable'):
            hash(value)
