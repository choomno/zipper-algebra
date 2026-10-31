"""Proof contexts are hints; exact word equations decide acceptance."""
from dataclasses import replace
from fractions import Fraction

import pytest

from zipper_algebra import Presentation
from zipper_algebra.coefficients import NumericRing
from zipper_algebra.cofactor import discovery as D
from zipper_algebra.cofactor.proof import Basis, Proof
from zipper_algebra.cofactor.skeleton import Skeleton
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.matrices import ExactMatrix
from zipper_algebra.untwisting import UntwistingPolynomial
from zipper_algebra.verification.linear import KernelCertificate


@pytest.fixture
def graph():
    pres = Presentation.from_untwistings(3, {3: UntwistingPolynomial([1])})
    polynomial = pres.target(3, 1, NumericRing(2, 3, 0))
    basis = Basis([polynomial], 0, pres.var_names())
    index = pres.labels().index('T_3')
    proof = Proof([[(Fraction(1), (), 'input', index, ())]])
    skeleton = Skeleton(pres.slug(), pres.content_hash(), 'X_3', (2, 3, 0),
                        [{'left': (), 'basis_key': basis.key[0], 'right': (), 'mult': Fraction(1)}],
                        [basis.key[0]])
    return pres, basis, proof, skeleton, index


def build(graph, **kwargs):
    pres, basis, proof, skeleton, _ = graph
    return D.build_system(pres, 'X_3', skeleton, basis, proof, **kwargs)


def test_word_equations_have_independent_expected_sign_and_full_generic_input(graph):
    system = build(graph)
    assert system.contexts == (((), graph[-1], ()),)
    values = {(i, j): value for i, j, value in system.matrix.entries}
    assert system.matrix.ncols == 2
    assert system.matrix.nrows == len(system.words)
    assert all(values[i, 0] == -values[i, 1] for i in range(system.matrix.nrows))
    target = graph[0].symbolic_ncpolys()[0][graph[-1]]
    assert {word: dict(values[i, 0].numerator) for i, word in enumerate(system.words)} == target
    saved = system.to_json()
    assert saved['target_column'] == 1
    assert saved['matrix_hash'] == system.matrix.content_hash()
    assert saved['contexts'] == [[[], graph[-1], []]]


def test_union_expansion_preserves_zero_paths_and_noncommutative_context_order(graph):
    pres, basis, proof, skeleton, index = graph
    # A forward edge and two cancelling numeric paths still contribute their
    # syntactic contexts. Coefficients are deliberately irrelevant here.
    basis = Basis([basis.polys[0], basis.polys[0]], 0, pres.var_names())
    proof = Proof([[(0, (1,), 'basis', 1, (2,))],
                   [(1, (3,), 'input', index, (4,)),
                    (-1, (3,), 'input', index, (4,)),
                    (0, (), 'input', index, ())]])
    skeleton = replace(skeleton, steps=[{'left': (2,), 'basis_key': basis.key[0],
                                         'right': (1,), 'mult': 0}] * 2)
    system = D.build_system(pres, 'X_3', skeleton, basis, proof)
    assert set(system.contexts) == {((2, 1, 3), index, (4, 2, 1)),
                                    ((2, 1), index, (2, 1))}
    assert system.statistics['needed_basis_nodes'] == 2
    assert system.statistics['distinct_contexts'] == 2


@pytest.mark.parametrize('name,value', [('max_contexts', 0), ('max_contexts', True),
                                        ('max_products', -1), ('max_products', 1.5)])
def test_invalid_budget(graph, name, value):
    with pytest.raises(ValueError, match='budgets'):
        build(graph, **{name: value})


@pytest.mark.parametrize('limit', [1, 2])
def test_bounded_context_substitution_is_inconclusive(graph, limit):
    with pytest.raises(RuntimeError, match='budget exhausted'):
        build(graph, max_products=limit)
    if limit == 1:
        with pytest.raises(RuntimeError, match='storage budget'):
            build(graph, max_contexts=limit)


@pytest.mark.parametrize('change', [{'presentation_hash': 'wrong'}, {'target': 'X_2'},
                                    {'remainder_terms': 1}])
def test_skeleton_binding_and_nonzero_remainder(graph, change):
    pres, basis, proof, skeleton, _ = graph
    with pytest.raises(ValueError, match='zero-remainder skeleton'):
        D.build_system(pres, 'X_3', replace(skeleton, **change), basis, proof)


def test_basis_domain_shape_names_and_missing_key(graph):
    pres, basis, proof, skeleton, _ = graph
    for invalid in [Basis(basis.polys, 7, basis.var_names),
                    Basis(basis.polys, 0, ['a', 'b', 'c', 'd']), Basis([], 0, basis.var_names)]:
        with pytest.raises(ValueError, match='basis must match'):
            D.build_system(pres, 'X_3', skeleton, invalid, proof)
    with pytest.raises(ValueError, match='existing basis keys'):
        D.build_system(pres, 'X_3', replace(skeleton, steps=[{'left': (),
                       'basis_key': 'absent', 'right': ()}]), basis, proof)


@pytest.mark.parametrize('term', [(1, (), 'basis', 8, ()), (1, (), 'input', -1, ()),
                                  (1, (), 'input', True, ()), (1, (), 'bad', 0, ()),
                                  (1, (99,), 'input', 0, ()), (1, (), 'input', 0, (False,))])
def test_malformed_proof_reference(graph, term):
    pres, basis, _, skeleton, _ = graph
    with pytest.raises(ValueError):
        D.build_system(pres, 'X_3', skeleton, basis, Proof([[term]]))


def test_cycle_is_rejected_without_recursion(graph):
    pres, basis, _, skeleton, _ = graph
    with pytest.raises(ValueError, match='cyclic proof graph'):
        D.build_system(pres, 'X_3', skeleton, basis, Proof([[(1, (), 'basis', 0, ())]]))


def test_empty_roots_and_empty_proof_line_are_valid_finite_searches(graph):
    pres, basis, proof, skeleton, _ = graph
    for steps, lines in [([], proof.terms), (skeleton.steps, [[]])]:
        system = D.build_system(pres, 'X_3', replace(skeleton, steps=steps), basis, Proof(lines))
        assert not system.contexts
        assert system.matrix.ncols == 1
        assert system.matrix.nrows > 0


def certificate(system):
    f = system.matrix.field
    # A=h; target is a dependent coordinate in this valid alternative layout.
    kernel = ExactMatrix.from_rows(f, [[1], [1]])
    return KernelCertificate(system.matrix.content_hash(), kernel, (0,), (0,), (1,),
                             101, (2, 3), ((100,),))


def test_cofactor_conversion_clears_rational_polynomial_lcm(graph):
    pytest.importorskip('sympy')
    f = RationalFunctionField(('q', 't'))
    q, t = f.parameter('q'), f.parameter('t')
    coefficients = [Fraction(2, 3) / (q + 1), (t + 1) / (q + 1)**2,
                    f.zero, f.scalar(7)]
    pres, _, _, _, index = graph
    contexts = [((), index, ())] * len(coefficients)
    cofactor = D._cofactor(pres, 'X_3', contexts, coefficients, (2, 3, 0))
    assert len(cofactor.atoms) == 1
    assert len(cofactor.terms) == 3
    assert cofactor.atoms.atoms[0].coefs == {(2, 0): 1, (1, 0): 2, (0, 0): 1}
    for term, expected in zip(cofactor.terms, [coefficients[0], coefficients[1], coefficients[3]]):
        _, _, _, numerator, powers = term
        actual = f.polynomial(numerator) / f.polynomial(cofactor.atoms.atoms[0].coefs)
        assert actual == expected
        assert powers == (1,)


def test_cofactor_conversion_integer_denominator_and_empty_result(graph):
    pytest.importorskip('sympy')
    f = RationalFunctionField(('q', 't'))
    pres, _, _, _, index = graph
    for coefficients in [[f.scalar(Fraction(2, 3))], [f.zero], []]:
        cofactor = D._cofactor(pres, 'X_3', [((), index, ())] * len(coefficients),
                              coefficients, None)
        assert not len(cofactor.atoms)
        assert len(cofactor.terms) == int(bool(coefficients and not coefficients[0].is_zero()))
        if cofactor.terms:
            assert cofactor.terms[0][3] == {(0, 0): Fraction(2, 3)}
            assert cofactor.terms[0][4] == ()


def test_reconstruction_keeps_independent_kernel_and_free_algebra_checks(graph, monkeypatch):
    pytest.importorskip('sympy')
    system = build(graph)
    result = certificate(system)
    calls = []
    monkeypatch.setattr(D, 'reconstruct_kernel', lambda matrix, **kw: calls.append(('produce', matrix, kw)) or result)
    original_verify = D.verify_kernel
    def checked(matrix, cert, **kw):
        calls.append(('check', matrix, cert, kw))
        return original_verify(matrix, cert, **kw)
    monkeypatch.setattr(D, 'verify_kernel', checked)
    discovered = D.reconstruct_cofactor(graph[0], 'X_3', system, anchor=(2, 3, 0),
                                       max_primes=3, max_degree=8, threads=2,
                                       start_prime=1, shift=(5, 7), max_products=1234)
    assert discovered.system is system
    assert discovered.kernel is result
    assert discovered.exact['pass']
    assert discovered.cofactor.n_terms == 1
    assert discovered.cofactor.anchor == (2, 3, 0)
    assert discovered.cofactor.terms[0][3] == {(0, 0): 1}
    assert [entry[0] for entry in calls] == ['produce', 'check']
    assert calls[0][2] == dict(max_primes=3, max_degree=8, threads=2, start_prime=1,
                               shift=(5, 7), max_products=1234)


@pytest.mark.parametrize('change', [{'presentation_hash': 'wrong'}, {'target': 'X_2'},
                                    {'matrix': ExactMatrix(RationalFunctionField(), 1, 2)},
                                    {'contexts': ()}])
def test_system_binding_before_native_use(graph, change):
    with pytest.raises(ValueError, match='context system does not match'):
        D.reconstruct_cofactor(graph[0], 'X_3', replace(build(graph), **change))


def test_exact_checks_fail_closed_and_missing_target_coordinate_is_inconclusive(graph, monkeypatch):
    pytest.importorskip('sympy')
    system = build(graph)
    result = certificate(system)
    monkeypatch.setattr(D, 'reconstruct_kernel', lambda *a, **kw: result)
    def reject(*args, **kwargs):
        raise ValueError('bad generic kernel')
    monkeypatch.setattr(D, 'verify_kernel', reject)
    with pytest.raises(RuntimeError, match='failed independent exact verification'):
        D.reconstruct_cofactor(graph[0], 'X_3', system)
    monkeypatch.setattr(D, 'verify_kernel', lambda *a, **kw: None)
    monkeypatch.setattr(D, 'check', lambda *a, **kw: {'pass': False})
    with pytest.raises(RuntimeError, match='free-algebra identity'):
        D.reconstruct_cofactor(graph[0], 'X_3', system)
    result = replace(result, kernel=ExactMatrix(system.matrix.field, 2, 0))
    with pytest.raises(RuntimeError, match='finite proof-graph context family'):
        D.reconstruct_cofactor(graph[0], 'X_3', system)


def test_actual_native_cofactor_reconstruction(graph):
    pytest.importorskip('sympy')
    pytest.importorskip('fflow')
    discovered = D.reconstruct_cofactor(graph[0], 'X_3', build(graph), max_primes=3,
                                       max_degree=8, threads=2)
    assert discovered.exact['pass']
    assert discovered.kernel.nullity == 1
    assert discovered.cofactor.terms == [((), 'T_3', (), {(0, 0): Fraction(1)}, ())]


def test_nonunit_dependent_target_coordinate_is_normalized_exactly(graph, monkeypatch):
    pytest.importorskip('sympy')
    _, basis, _, skeleton, index = graph
    proof = Proof([[(Fraction(1, 2), (), 'input', index, ())]])
    pres = Presentation.from_untwistings(3, {3: UntwistingPolynomial([2])})
    skeleton = replace(skeleton, presentation_hash=pres.content_hash())
    system = D.build_system(pres, 'X_3', skeleton, basis, proof)
    result = certificate(system)
    result = replace(result, kernel=ExactMatrix.from_rows(system.matrix.field, [[1], [2]]))
    # This is a valid full kernel with the *first* row the identity block; the
    # target row is dependent, its value is two, and the producer is untrusted.
    assert D.verify_kernel(system.matrix, result) is result
    monkeypatch.setattr(D, 'reconstruct_kernel', lambda *a, **kw: result)
    discovery = D.reconstruct_cofactor(pres, 'X_3', system)
    assert discovery.exact['pass']
    assert discovery.cofactor.terms[0][3] == {(0, 0): Fraction(1, 2)}


def test_fractional_polynomial_denominator_is_cleared_to_integral_atom(graph):
    pytest.importorskip('sympy')
    f = RationalFunctionField(('q', 't'))
    value = f.one / (f.parameter('q') + Fraction(1, 2))
    pres, _, _, _, index = graph
    cofactor = D._cofactor(pres, 'X_3', [((), index, ())], [value], None)
    assert cofactor.atoms.atoms[0].coefs == {(1, 0): 2, (0, 0): 1}
    assert cofactor.terms[0][3] == {(0, 0): Fraction(2)}


def test_deep_forward_proof_graph_does_not_require_python_recursion(graph):
    pres, basis, _, skeleton, index = graph
    length = 1200
    basis = Basis([basis.polys[0]] * length, 0, pres.var_names())
    proof = Proof([[(1, (), 'basis', i + 1, ())] for i in range(length - 1)]
                  + [[(1, (), 'input', index, ())]])
    skeleton = replace(skeleton, steps=[{'left': (), 'basis_key': basis.key[0], 'right': ()}])
    system = D.build_system(pres, 'X_3', skeleton, basis, proof)
    assert system.contexts == (((), index, ()),)
    assert system.statistics['needed_basis_nodes'] == length
    assert system.statistics['stored_contexts'] == length + 1


def test_valid_kernel_for_forged_context_mapping_fails_actual_free_algebra_check(graph, monkeypatch):
    pytest.importorskip('sympy')
    system = build(graph)
    result = certificate(system)
    # The unchanged word equations have a perfectly valid kernel, but this
    # context now names INV_1_p instead of T_3. Only the separate free-algebra
    # accumulation can reject that semantic discrepancy.
    system = replace(system, contexts=(((), 0, ()),))
    monkeypatch.setattr(D, 'reconstruct_kernel', lambda *a, **kw: result)
    with pytest.raises(RuntimeError, match='free-algebra identity'):
        D.reconstruct_cofactor(graph[0], 'X_3', system)


@pytest.mark.parametrize('step', [{}, {'left': (), 'right': (), 'basis_key': []}])
def test_malformed_skeleton_step_has_input_error_classification(graph, step):
    pres, basis, proof, skeleton, _ = graph
    with pytest.raises(ValueError, match='existing basis keys'):
        D.build_system(pres, 'X_3', replace(skeleton, steps=[step]), basis, proof)


def test_unreachable_cycle_does_not_change_finite_context_closure(graph):
    pres, basis, proof, skeleton, index = graph
    basis = Basis([basis.polys[0]] * 3, 0, pres.var_names())
    proof = Proof(proof.terms + [[(1, (), 'basis', 2, ())], [(1, (), 'basis', 1, ())]])
    skeleton = replace(skeleton, steps=[{'left': (), 'basis_key': basis.key[0], 'right': ()}])
    system = D.build_system(pres, 'X_3', skeleton, basis, proof)
    assert system.contexts == (((), index, ()),)
    assert system.statistics['needed_basis_nodes'] == 1


@pytest.mark.parametrize('error_type', [TypeError, ZeroDivisionError])
def test_invalid_producer_verification_failures_are_inconclusive(graph, monkeypatch, error_type):
    system = build(graph)
    monkeypatch.setattr(D, 'reconstruct_kernel', lambda *a, **kw: certificate(system))
    def reject(*args, **kwargs):
        raise error_type('invalid reconstructed candidate')
    monkeypatch.setattr(D, 'verify_kernel', reject)
    with pytest.raises(RuntimeError, match='failed independent exact verification'):
        D.reconstruct_cofactor(graph[0], 'X_3', system)


def test_unreduced_kernel_factors_do_not_inflate_the_clearing_polynomial(graph):
    pytest.importorskip('sympy')
    f = RationalFunctionField(('q', 't'))
    q, t = f.parameter('q'), f.parameter('t')
    value = ((q + 1) * (t + 1)) / ((q + 1) * (t + 2))
    pres, _, _, _, index = graph
    cofactor = D._cofactor(pres, 'X_3', [((), index, ())], [value], None)
    assert cofactor.atoms.atoms[0].coefs == {(0, 1): 1, (0, 0): 2}
    assert cofactor.terms[0][3] == {(0, 1): Fraction(1), (0, 0): Fraction(1)}


def test_multiple_target_kernel_columns_choose_the_first_valid_identity(graph, monkeypatch):
    pytest.importorskip('sympy')
    pres = Presentation.from_untwistings(4, {4: UntwistingPolynomial([1])}, extra=('X_4',))
    polynomial = pres.target(4, 1, NumericRing(2, 3, 0))
    basis = Basis([polynomial], 0, pres.var_names())
    indices = [pres.labels().index(label) for label in ('T_4', 'EXTRA[X_4]')]
    proof = Proof([[(Fraction(1, 2), (), 'input', index, ()) for index in indices]])
    skeleton = Skeleton(pres.slug(), pres.content_hash(), 'X_4', (2, 3, 0),
                        [{'left': (), 'basis_key': basis.key[0], 'right': ()}])
    system = D.build_system(pres, 'X_4', skeleton, basis, proof)
    # Both presentation columns equal h. This alternative full kernel has
    # target coordinate one in *both* columns and an identity on rows 0,1.
    kernel = ExactMatrix.from_rows(system.matrix.field, [[1, 0], [0, 1], [1, 1]])
    point, prime = (2, 3), 101
    minor_values = [(row, value.specialize(dict(zip(('q', 't'), point)), characteristic=prime))
                    for row, column, value in system.matrix.entries if column == 2]
    row, value = next((row, value) for row, value in minor_values if not value.is_zero())
    inverse = pow(int(value.numerator[0][1]), -1, prime)
    result = KernelCertificate(system.matrix.content_hash(), kernel, (0, 1), (row,), (2,),
                               prime, point, ((inverse,),))
    assert D.verify_kernel(system.matrix, result) is result
    monkeypatch.setattr(D, 'reconstruct_kernel', lambda *a, **kw: result)
    discovered = D.reconstruct_cofactor(pres, 'X_4', system)
    assert discovered.exact['pass']
    assert discovered.cofactor.terms == [((), 'T_4', (), {(0, 0): Fraction(1)}, ())]
