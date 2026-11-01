"""Bounded original-input support and independently checked native-wire candidates."""
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
from types import SimpleNamespace

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.cofactor import reconstruction
from zipper_algebra.cofactor._contexts import _collect_contexts, _validate_graph
from zipper_algebra.cofactor.general import verify_cofactor
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.derivations import DerivationGraph, DerivationTerm, VerificationLimit, _Budget


def presentation(parent, relations):
    return AlgebraPresentation(parent, tuple(relations), tuple(f'r{i}' for i in range(len(relations))))


def term(field, index, *, left=(), right=(), value=1, kind='input'):
    return DerivationTerm(field.scalar(value), left, kind, index, right)


def fixture(parameters=()):
    field = RationalFunctionField(parameters)
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = map(parent.generator, parent.generators)
    p = presentation(parent, (x, y))
    graph = DerivationGraph(p.content_hash(), (x,), ((term(field, 0),),))
    return p, graph, x, y


def support(p, graph, node, *, operations=100000, contexts=100000, nodes=100000):
    budget = _Budget(operations)
    order = _validate_graph(p, graph, budget, nodes)
    result, statistics = _collect_contexts(p, graph, node, order, budget, contexts)
    return result, statistics


def test_structural_support_contains_canceling_contexts_without_coefficient_propagation():
    p, _, x, y = fixture()
    field = p.algebra.field
    graph = DerivationGraph(p.content_hash(), (x,), ((
        term(field, 0), term(field, 1, left=(1,), value=3),
        term(field, 1, left=(1,), value=-3),
        term(field, 1, left=(2,), value=0)),))
    contexts, _ = support(p, graph, 0)
    assert set(contexts) == {((), 0, ()), ((1,), 1, ())}
    assert graph.nodes[0] == x and not (x+y).is_zero()


def test_forward_shared_contexts_preserve_noncommuting_prefix_and_suffix_order():
    p, _, x, y = fixture()
    field = p.algebra.field
    graph = DerivationGraph(p.content_hash(), (y*x*y*x, x), (
        (term(field, 1, left=(2,), right=(2, 1), kind='node'),),
        (term(field, 0),)))
    contexts, statistics = support(p, graph, 0)
    assert contexts == (((2,), 0, (2, 1)),)
    assert statistics['needed_nodes'] == 2


def test_last_user_release_bounds_live_sets_in_a_forward_chain():
    p, _, x, _ = fixture()
    field = p.algebra.field
    count = 25
    rows = tuple((term(field, i+1, kind='node'),) if i<count-1
                 else (term(field, 0),) for i in range(count))
    graph = DerivationGraph(p.content_hash(), tuple(x for _ in rows), rows)
    contexts, statistics = support(p, graph, 0, contexts=2)
    assert contexts == (((), 0, ()),)
    assert statistics['expanded_nodes'] == count
    assert statistics['peak_contexts'] == 2
    with pytest.raises(VerificationLimit):
        support(p, graph, 0, contexts=1)


def test_zero_inputs_and_zero_edges_do_not_create_structural_columns():
    p, _, x, _ = fixture()
    p = presentation(p.algebra, (p.algebra.zero, x))
    graph = DerivationGraph(p.content_hash(), (x, x), (
        (term(p.algebra.field, 0), term(p.algebra.field, 1),
         term(p.algebra.field, 1, kind='node', value=0)), ()))
    contexts, statistics = support(p, graph, 0)
    assert contexts == (((), 1, ()),)
    assert statistics['needed_nodes'] == 1


def test_unit_input_contexts_include_empty_word_and_all_ordered_splits():
    p, _, x, y = fixture()
    p = presentation(p.algebra, (p.algebra.one,))
    f = p.algebra.field
    graph = DerivationGraph(p.content_hash(), (x*y,), ((
        term(f, 0, left=(), right=(1, 2)),
        term(f, 0, left=(1,), right=(2,)),
        term(f, 0, left=(1, 2), right=())),))
    assert set(support(p, graph, 0)[0]) == {
        ((), 0, (1, 2)), ((1,), 0, (2,)), ((1, 2), 0, ())}


@pytest.mark.parametrize('bad', ['index', 'zero-index', 'zero-cycle', 'unrelated-cycle'])
def test_every_raw_reference_is_validated_before_context_pruning(bad):
    p, _, x, _ = fixture()
    f = p.algebra.field
    if bad in ('index', 'zero-index'):
        rows = ((term(f, 0),), (term(f, 99, value=0 if bad=='zero-index' else 1),))
    elif bad=='zero-cycle':
        rows = ((term(f, 1, kind='node'),), (term(f, 0, kind='node', value=0),))
    else:
        rows = ((term(f, 0),), (term(f, 1, kind='node', value=0),))
    graph = DerivationGraph(p.content_hash(), (x, x), rows)
    with pytest.raises(ValueError):
        support(p, graph, 0)


def test_node_hash_parent_types_and_graph_budget_are_checked():
    p, graph, _, _ = fixture()
    with pytest.raises(ValueError):
        support(p, replace(graph, presentation_hash='0'*64), 0)
    other = FreeAlgebra(p.algebra.field, ('a', 'b'))
    foreign = DerivationGraph(p.content_hash(), (other.zero,), ((),))
    with pytest.raises(TypeError):
        support(p, foreign, 0)
    with pytest.raises(VerificationLimit):
        support(p, graph, 0, operations=1)
    doubled = DerivationGraph(p.content_hash(), graph.nodes*2, graph.derivations*2)
    with pytest.raises(VerificationLimit):
        support(p, doubled, 0, nodes=1)


class AnalyticOutput:
    """Actual RatFunList wire methods, with explicit adversarial values."""
    def __init__(self, values=None, arity=1):
        self.values = values if values is not None else [([('1', (0,))], [('1', (0,))])]
        self.arity = arity

    def __len__(self):
        return len(self.values)

    def nvars(self):
        return self.arity

    def monomials(self, index):
        return self.values[index]


class NativeWire:
    """Faithful graph-wire facade, never an implementation of linear solving."""
    RatFunList = AnalyticOutput

    class Failed(Exception):
        pass

    def __init__(self, numeric=None, analytic=None, dependent=None):
        self.numeric = ['1'] if numeric is None else numeric
        self.analytic = AnalyticOutput() if analytic is None else analytic
        self.dependent = dependent
        self.calls = []
        self.deleted = []

    def NAvailablePrimes(self):
        return 2048

    def PrimeNo(self, index):
        return 2147483647

    def NewGraphWithInput(self, parameters):
        self.calls.append(('new', parameters))
        return 17, 23

    def NewIdxRatFunList(self, parameters, functions, indexes=None):
        self.analytic_input = parameters, functions, indexes
        return functions

    def SetReconstructionShift(self, graph, coordinates):
        self.calls.append(('shift', coordinates))

    def AlgNumericSparseLSolve(self, graph, variables, rows, coefficients):
        assert sum(map(len, rows)) == len(coefficients)
        self.shape = len(rows), variables
        self.numeric_input = rows, coefficients
        return 29

    def AlgAnalyticSparseLSolve(self, graph, input_node, variables, rows, coefficients):
        self.shape = len(rows), variables
        return 29

    def LSolveValidateFullSystem(self, graph, node, enabled):
        self.calls.append(('full', enabled))

    def LSolveOnlyNonHomogeneous(self, graph, node):
        self.calls.append(('particular',))

    def LSolveSparseOutput(self, graph, node, enabled):
        self.calls.append(('sparse', enabled))

    def LSolveEqWeight(self, graph, node, weights):
        self.calls.append(('weights', weights))

    def SetOutputNode(self, graph, node):
        pass

    def LearnEx(self, graph, *, prime_no):
        self.calls.append(('learn', prime_no))

    def LSolveIsImpossible(self, graph, node):
        return False

    def LSolveNEqsNVars(self, graph, node):
        return self.shape

    def LSolveNDepVars(self, graph, node):
        return self.shape[1]

    def LSolveDepVars(self, graph, node):
        return list(range(self.shape[1])) if self.dependent is None else self.dependent

    def LSolveNIndepVars(self, graph, node):
        return 0

    def GraphNParsOut(self, graph):
        return self.shape[1]

    def LSolveOutputIsSparse(self, graph, node):
        return False

    def ReconstructNumeric(self, graph, **options):
        self.calls.append(('numeric', options))
        return self.numeric

    def ReconstructFunction(self, graph, **options):
        self.calls.append(('analytic', options))
        return self.analytic

    def DeleteGraph(self, graph):
        self.deleted.append(graph)


def invoke(monkeypatch, native=None, *, p=None, graph=None, **options):
    if p is None:
        p, graph, _, _ = fixture()
    native = NativeWire() if native is None else native
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=lambda _: native))
    return reconstruction.reconstruct_derivation(p, graph, 0, **options)


def test_numeric_wire_preserves_exact_rational_rhs_and_graph_mode_is_local(monkeypatch):
    p, _, x, _ = fixture()
    p = presentation(p.algebra, (Fraction(2, 3)*x,))
    graph = DerivationGraph(p.content_hash(), (x,), ((term(p.algebra.field, 0, value=9),),))
    native = NativeWire(numeric=['3/2'])
    result = invoke(monkeypatch, native, p=p, graph=graph, start_prime=1, threads=2)
    assert verify_cofactor(p, result.identity)
    assert native.numeric_input == ([[0, 1]], ['2/3', '1'])
    assert ('full', False) in native.calls and ('particular',) in native.calls
    assert ('learn', 1) in native.calls and native.deleted == [17]
    assert next(call[1] for call in native.calls if call[0]=='numeric') == {
        'n_threads': 2, 'max_primes': 5, 'start_mod': 1}
    stats = result.statistics
    stats.clear()
    assert result.statistics
    assert all(type(value) is int for value in result.statistics.values())
    with pytest.raises((FrozenInstanceError, AttributeError)):
        result.identity = None


@pytest.mark.parametrize('output', [[], ['1', '0'], (), object(), ['not-a-rational'], ['1/0'], ['2']])
def test_incomplete_malformed_or_false_numeric_candidates_are_rejected(monkeypatch, output):
    native = NativeWire()
    native.numeric = output
    with pytest.raises(RuntimeError):
        invoke(monkeypatch, native)
    assert native.deleted == [17]


@pytest.mark.parametrize('method,value', [
    ('LSolveIsImpossible', True), ('LSolveNEqsNVars', (2, 1)),
    ('LSolveNDepVars', 0), ('LSolveNDepVars', True),
    ('LSolveNIndepVars', 1), ('LSolveDepVars', [True]),
    ('LSolveDepVars', [0, 0]), ('LSolveDepVars', [1]),
    ('GraphNParsOut', 2), ('LSolveOutputIsSparse', True),
])
def test_native_layout_faults_fail_closed_and_release_graph(monkeypatch, method, value):
    native = NativeWire()
    setattr(native, method, lambda *args: value)
    with pytest.raises(RuntimeError):
        invoke(monkeypatch, native)
    assert native.deleted == [17]


def test_unsorted_native_ids_map_values_to_original_inputs_and_drop_forced_zero(monkeypatch):
    p, _, x, y = fixture()
    p = presentation(p.algebra, (x, y, x+y))
    f = p.algebra.field
    graph = DerivationGraph(p.content_hash(), (x+y,),
        ((term(f, 0), term(f, 1), term(f, 2)),))
    native = NativeWire(numeric=['0', '1', '1'], dependent=[2, 1, 0])
    result = invoke(monkeypatch, native, p=p, graph=graph)
    assert verify_cofactor(p, result.identity)
    assert [index for _, index, _ in result.identity.terms] == [0, 1]
    assert result.statistics['candidate_contexts'] == 3
    assert result.statistics['contexts'] == result.statistics['cofactor_terms'] == 2


def test_parameter_wire_decodes_a_generic_nonconstant_cofactor_and_custom_shift(monkeypatch):
    pytest.importorskip('sympy')
    p, _, x, _ = fixture(('q',))
    q = p.algebra.field.parameter('q')
    graph = DerivationGraph(p.content_hash(), (q*x,), ((term(p.algebra.field, 0),),))
    output = AnalyticOutput([([('1', (1,))], [('1', (0,))])])
    native = NativeWire(analytic=output)
    result = invoke(monkeypatch, native, p=p, graph=graph, shift=(3,), max_degree=8)
    assert verify_cofactor(p, result.identity)
    assert ('shift', [3]) in native.calls
    assert native.deleted == [17]
    assert result.identity.terms[0][0] == p.algebra.scalar(q)


@pytest.mark.parametrize('output', [
    object(), AnalyticOutput([], 1), AnalyticOutput(arity=2),
    AnalyticOutput([([('1', ())], [('1', (0,))])]),
    AnalyticOutput([([('1', (0,)), ('1', (0,))], [('1', (0,))])]),
    AnalyticOutput([([('1', (0,))], [('1', (0,)), ('1', (0,))])]),
    AnalyticOutput([([('1', (0,))], [])]),
    AnalyticOutput([([('1/0', (0,))], [('1', (0,))])]),
    AnalyticOutput([([('2', (0,))], [('1', (0,))])]),
])
def test_bad_analytic_wire_data_and_false_generic_formula_are_rejected(monkeypatch, output):
    pytest.importorskip('sympy')
    p, graph, _, _ = fixture(('q',))
    native = NativeWire(analytic=output)
    with pytest.raises(RuntimeError):
        invoke(monkeypatch, native, p=p, graph=graph)
    assert native.deleted == [17]


@pytest.mark.parametrize('which', ['import', 'fork', 'shift'])
def test_missing_optional_dependency_or_required_capability_is_explicit(monkeypatch, which):
    if which=='import':
        def missing(name):
            raise ImportError('bindings absent')
        monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=missing))
        p, graph, _, _ = fixture()
        with pytest.raises(ImportError):
            reconstruction.reconstruct_derivation(p, graph, 0)
    else:
        native = NativeWire()
        missing_name = 'LSolveValidateFullSystem' if which=='fork' else 'SetReconstructionShift'
        # A namespace contains exactly the callable surface exposed by this binding.
        native = SimpleNamespace(**{name:getattr(native, name) for name in dir(native)
                                    if not name.startswith('_') and name!=missing_name})
        with pytest.raises(ImportError):
            invoke(monkeypatch, native, **({'shift':()} if which=='shift' else {}))


def test_native_exception_does_not_leak_a_graph(monkeypatch):
    native = NativeWire()
    def fail(*args, **kwargs):
        raise native.Failed('learning failure')
    native.LearnEx = fail
    with pytest.raises(RuntimeError):
        invoke(monkeypatch, native)
    assert native.deleted == [17]


def test_zero_target_needs_no_optional_native_dependency(monkeypatch):
    p, _, _, _ = fixture()
    graph = DerivationGraph(p.content_hash(), (p.algebra.zero,), ((),))
    def unexpected(name):
        raise AssertionError('zero target must need no native graph')
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=unexpected))
    result = reconstruction.reconstruct_derivation(p, graph, 0)
    assert result.identity.terms == () and verify_cofactor(p, result.identity)
    assert result.statistics['candidate_contexts'] == result.statistics['contexts'] == 0
    assert result.statistics['cofactor_terms'] == 0


@pytest.mark.parametrize('which', ['presentation', 'graph'])
def test_public_types_are_checked_before_native_import(which):
    p, graph, _, _ = fixture()
    with pytest.raises(TypeError):
        reconstruction.reconstruct_derivation(None if which=='presentation' else p,
            None if which=='graph' else graph, 0)


@pytest.mark.parametrize('node', [-1, 1, True, 0.5, None])
def test_public_root_must_be_an_existing_integer(node):
    p, graph, _, _ = fixture()
    with pytest.raises((ValueError, TypeError)):
        reconstruction.reconstruct_derivation(p, graph, node)


@pytest.mark.parametrize('option', ['max_operations', 'max_contexts', 'max_nodes',
    'max_matrix_entries', 'max_primes', 'max_degree', 'threads'])
@pytest.mark.parametrize('value', [0, True, -1])
def test_public_bounds_are_positive_exact_integers(option, value):
    p, graph, _, _ = fixture()
    with pytest.raises(ValueError):
        reconstruction.reconstruct_derivation(p, graph, 0, **{option:value})


@pytest.mark.parametrize('options', [
    {'max_primes':2048}, {'max_degree':2**16}, {'threads':2**32},
    {'start_prime':-1}, {'start_prime':True}, {'start_prime':1.5},
    {'start_prime':2047, 'max_primes':1},
    {'shift':(1,)}, {'shift':(True,)}, {'shift':(-1,)}, {'shift':(2**64,)}, {'shift':1},
])
def test_native_options_reject_out_of_range_or_inexact_controls(options):
    p, graph, _, _ = fixture()
    with pytest.raises((ValueError, TypeError)):
        reconstruction.reconstruct_derivation(p, graph, 0, **options)


@pytest.mark.parametrize('limit', ['max_operations', 'max_contexts', 'max_nodes', 'max_matrix_entries'])
def test_real_production_resource_caps_are_explicit_before_native_import(monkeypatch, limit):
    p, _, x, y = fixture()
    field = p.algebra.field
    graph = DerivationGraph(p.content_hash(), (x*y, x),
        ((term(field, 1, right=(2,), kind='node'),), (term(field, 0),)))
    def unexpected(name):
        raise AssertionError('resource cap should precede native import')
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=unexpected))
    with pytest.raises(VerificationLimit):
        reconstruction.reconstruct_derivation(p, graph, 0, **{limit:1})


def test_positive_characteristic_is_not_silently_sent_to_rational_native_backend():
    field = RationalFunctionField((), 5)
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    p = presentation(parent, (x,))
    graph = DerivationGraph(p.content_hash(), (x,), ((term(field, 0),),))
    with pytest.raises(ValueError, match='characteristic'):
        reconstruction.reconstruct_derivation(p, graph, 0)


def test_true_zero_edges_cannot_hide_a_nonzero_target_from_the_exact_gate(monkeypatch):
    p, _, x, y = fixture()
    graph = DerivationGraph(p.content_hash(), (y,), ((term(p.algebra.field, 0),),))
    native = NativeWire(numeric=['1'])
    with pytest.raises(RuntimeError, match='exact|identity|verification'):
        invoke(monkeypatch, native, p=p, graph=graph)
    assert native.deleted == [17]


def test_nonzero_target_with_only_zero_source_support_is_rejected_before_native_import(monkeypatch):
    p, _, x, _ = fixture()
    p = presentation(p.algebra, (p.algebra.zero,))
    graph = DerivationGraph(p.content_hash(), (x,), ((term(p.algebra.field, 0),),))
    def unexpected(name):
        raise AssertionError('no nonzero context support must not allocate native state')
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=unexpected))
    with pytest.raises(RuntimeError):
        reconstruction.reconstruct_derivation(p, graph, 0)


def test_context_duplicates_are_unioned_and_known_support_is_not_capped_at_4000(monkeypatch):
    p, _, x, _ = fixture()
    count = 4001
    p = presentation(p.algebra, (x,)*count)
    f = p.algebra.field
    graph = DerivationGraph(p.content_hash(), (x,),
        (tuple(term(f, i) for i in range(count)) + (term(f, 0),),))
    native = NativeWire(numeric=['1']+['0']*(count-1))
    result = invoke(monkeypatch, native, p=p, graph=graph, max_contexts=5000)
    assert verify_cofactor(p, result.identity)
    assert len(result.identity.terms) == 1
    assert native.shape == (1, count)


def test_analytic_coefficient_owner_survives_native_graph_deletion(monkeypatch):
    import weakref
    pytest.importorskip('sympy')
    p, graph, _, _ = fixture(('q',))
    native = NativeWire()
    class Owner:
        pass
    owners = []
    def functions(*args):
        owner = Owner()
        owners.append(weakref.ref(owner))
        return owner
    def delete(graph):
        assert owners[-1]() is not None
        native.deleted.append(graph)
    native.NewIdxRatFunList = functions
    native.DeleteGraph = delete
    invoke(monkeypatch, native, p=p, graph=graph)
    assert owners[-1]() is None and native.deleted == [17]


@pytest.mark.parametrize('node', [-1, 1, True])
def test_structural_collection_rejects_an_invalid_root_even_after_valid_graph_order(node):
    p, graph, _, _ = fixture()
    with pytest.raises(ValueError):
        support(p, graph, node)


def test_shared_diamond_is_unioned_after_both_last_consumers():
    p, _, x, _ = fixture()
    field = p.algebra.field
    graph = DerivationGraph(p.content_hash(), (2*x, x, x, x), (
        (term(field, 1, kind='node'), term(field, 2, kind='node')),
        (term(field, 3, kind='node'),), (term(field, 3, kind='node'),),
        (term(field, 0),)))
    contexts, stats = support(p, graph, 0, contexts=3)
    assert contexts == (((), 0, ()),)
    assert stats['needed_nodes'] == stats['expanded_nodes'] == 4
    assert stats['peak_contexts'] == 3
    with pytest.raises(VerificationLimit):
        support(p, graph, 0, contexts=2)


@pytest.mark.parametrize('value', [True, 0.5])
def test_native_numeric_scalars_must_be_exact_and_not_booleans(monkeypatch, value):
    native = NativeWire(numeric=[value])
    with pytest.raises(RuntimeError, match='exact rational'):
        invoke(monkeypatch, native)
    assert native.deleted == [17]


@pytest.mark.parametrize('available', [0, True])
def test_invalid_native_prime_table_size_is_rejected_before_graph_creation(monkeypatch, available):
    native = NativeWire()
    native.NAvailablePrimes = lambda: available
    with pytest.raises(RuntimeError, match='prime-table'):
        invoke(monkeypatch, native)
    assert native.deleted == [] and not native.calls


def test_native_prime_budget_leaves_an_unused_checking_prime(monkeypatch):
    native = NativeWire()
    native.NAvailablePrimes = lambda: 6
    with pytest.raises(ValueError, match='checking prime'):
        invoke(monkeypatch, native, start_prime=1, max_primes=5)
    assert native.deleted == []


def test_shift_coordinates_respect_actual_native_prime_range(monkeypatch):
    pytest.importorskip('sympy')
    p, graph, _, _ = fixture(('q',))
    native = NativeWire()
    with pytest.raises(ValueError, match='smallest FiniteFlow prime'):
        invoke(monkeypatch, native, p=p, graph=graph, shift=(2147483647,))
    assert native.deleted == []


@pytest.mark.parametrize('shift', [(True,), (-1,), (2**64,), (1.5,)])
def test_parameter_shift_checks_each_exact_unsigned_coordinate(shift):
    p, graph, _, _ = fixture(('q',))
    with pytest.raises(ValueError, match='unsigned exact'):
        reconstruction.reconstruct_derivation(p, graph, 0, shift=shift)


def test_native_unsupported_exponent_is_invalid_input_before_graph_creation(monkeypatch):
    p, _, x, _ = fixture(('q',))
    coefficient = p.algebra.field.polynomial({(2**16,):1})
    p = presentation(p.algebra, (p.algebra.element({(1,):coefficient}),))
    graph = DerivationGraph(p.content_hash(), (x,), ((term(p.algebra.field, 0),),))
    native = NativeWire()
    with pytest.raises(ValueError, match='monomial exponents'):
        invoke(monkeypatch, native, p=p, graph=graph)
    assert native.deleted == [] and not native.calls


def test_bad_native_ratfun_structure_is_a_reconstruction_error_with_cleanup(monkeypatch):
    pytest.importorskip('sympy')
    p, graph, _, _ = fixture(('q',))
    output = AnalyticOutput([None])
    native = NativeWire(analytic=output)
    with pytest.raises(RuntimeError, match='malformed'):
        invoke(monkeypatch, native, p=p, graph=graph)
    assert native.deleted == [17]


def test_analytic_zero_coordinate_is_omitted_without_changing_exact_identity(monkeypatch):
    pytest.importorskip('sympy')
    p, _, x, _ = fixture(('q',))
    field = p.algebra.field
    graph = DerivationGraph(p.content_hash(), (x,), ((term(field, 0), term(field, 1)),))
    native = NativeWire(analytic=AnalyticOutput([
        ([('1', (0,))], [('1', (0,))]), ([], [('1', (0,))])]))
    result = invoke(monkeypatch, native, p=p, graph=graph)
    assert verify_cofactor(p, result.identity)
    assert len(result.identity.terms) == 1 and result.identity.terms[0][1] == 0


def test_rhs_entries_alone_can_exhaust_matrix_budget_before_native_allocation(monkeypatch):
    p, _, x, y = fixture()
    graph = DerivationGraph(p.content_hash(), (x+y,), ((term(p.algebra.field, 0),),))
    def unexpected(name):
        raise AssertionError('RHS storage cap must precede native import')
    monkeypatch.setattr(reconstruction, 'importlib', SimpleNamespace(import_module=unexpected))
    with pytest.raises(VerificationLimit, match='matrix-entry'):
        reconstruction.reconstruct_derivation(p, graph, 0, max_matrix_entries=1)


def test_scalar_unit_identity_requires_no_nonempty_alphabet(monkeypatch):
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ())
    p = presentation(parent, (parent.one,))
    graph = DerivationGraph(p.content_hash(), (parent.scalar(3),), ((term(field, 0),),))
    native = NativeWire(numeric=['3'])
    result = invoke(monkeypatch, native, p=p, graph=graph)
    assert verify_cofactor(p, result.identity)
    assert native.numeric_input == ([[0, 1]], ['1', '3'])


@pytest.mark.parametrize('shape', [(0, 0), (0, 5), (7, 0), (2**32-1, 2**32-1)])
def test_native_dimension_counts_accept_exact_uint32_boundaries_without_allocations(shape):
    assert reconstruction._native_dimensions(*shape) is None


@pytest.mark.parametrize('shape', [(-1, 1), (1, -1), (True, 1), (1, True),
    (2**32, 1), (1, 2**32), (1.0, 0), (0, None)])
def test_native_dimension_counts_reject_unrepresentable_or_inexact_shapes(shape):
    with pytest.raises(ValueError):
        reconstruction._native_dimensions(*shape)


def test_inherited_word_copying_consumes_the_structural_operation_budget():
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    length = 1000
    left = (1,)*length
    polynomial = parent.element({left+(1,):1})
    p = presentation(parent, (x,))
    graph = DerivationGraph(p.content_hash(), (polynomial, polynomial), (
        (term(field, 0, left=left),), (term(field, 0, kind='node'),)))
    # Both graph identities are literal equalities. Copying the inherited long
    # prefix at node 1 is real work even though that edge has empty contexts.
    with pytest.raises(VerificationLimit, match='operation budget'):
        support(p, graph, 1, operations=1500)
    contexts, stats = support(p, graph, 1, operations=3000)
    assert contexts == ((left, 0, ()),)
    assert stats['structural_operations'] >= 2*(length+1)
