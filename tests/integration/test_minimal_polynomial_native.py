"""Tiny installed-native computations with literal independent algebra witnesses."""
from fractions import Fraction

import pytest

from zipper_algebra.algebra import AlgebraPresentation, FreeAlgebra
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.invariants import minimal_polynomial
from zipper_algebra.verification.derivations import DerivationGraph, DerivationTerm
from zipper_algebra.verification.finite import GeneratorAction, WordBasisCertificate, verify_word_basis


@pytest.fixture
def native_finiteflow():
    pytest.importorskip('sympy', reason='generic reconstruction requires symbolic arithmetic')
    return pytest.importorskip('fflow', reason='native reconstruction is optional')


def matrix_algebra():
    """M2(Q), with x=E12, y=E21 and named basis (1,x,y,xy)."""
    field = RationalFunctionField()
    parent = FreeAlgebra(field, ('x', 'y'))
    x, y = parent.generator('x'), parent.generator('y')
    relations = (x*x, y*y, y*x+x*y-1)
    pres = AlgebraPresentation(parent, relations, ('xx', 'yy', 'yx'))
    def term(index, left=(), right=(), value=1):
        return DerivationTerm(field.scalar(value), left, 'input', index, right)
    def action(columns):
        return GeneratorAction(field, tuple(tuple((row, field.scalar(value)) for row, value in column)
                                            for column in columns))
    mx = action((((1, 1),), (), ((3, 1),), ()))
    my = action((((2, 1),), ((0, 1), (3, -1)), (), ((2, 1),)))
    nodes = (parent.zero, x*x, x*x*y, relations[2], y*y, y*x*y-y)
    proof = DerivationGraph(pres.content_hash(), nodes,
                            ((), (term(0),), (term(0, right=(2,)),), (term(2),),
                             (term(1),), (term(2, right=(2,)), term(1, left=(1,), value=-1))))
    cert = WordBasisCertificate(pres.content_hash(), ((), (1,), (2,), (1, 2)),
                                (mx, my), proof, ((0, 1, 0, 2), (0, 3, 4, 5)))
    return verify_word_basis(pres, cert)


def test_native_noncommutative_element_and_relative_minimal_polynomials(native_finiteflow):
    algebra = matrix_algebra()
    field = algebra.field
    x, y = algebra.element((0, 1, 0, 0)), algebra.element((0, 0, 1, 0))
    assert x*y != y*x
    e = x*y
    full = minimal_polynomial(algebra, e, max_primes=3, reconstruction_degree=4, threads=2)
    on_image = minimal_polynomial(algebra, e, seed=e, max_primes=3, reconstruction_degree=4)
    on_complement = minimal_polynomial(algebra, e, seed=algebra.one-e,
                                       max_primes=3, reconstruction_degree=4)
    assert full.coefficients == tuple(field.scalar(value) for value in (0, -1, 1))
    assert on_image.coefficients == tuple(field.scalar(value) for value in (-1, 1))
    assert on_complement.coefficients == tuple(field.scalar(value) for value in (0, 1))


def test_native_generic_quadratic_keeps_parameter_order_and_shift(native_finiteflow):
    field = RationalFunctionField(('c0', 'c1'))
    c0, c1 = field.parameter('c0'), field.parameter('c1')
    parent = FreeAlgebra(field, ('x',))
    x = parent.generator('x')
    relation = x*x+c1*x+c0
    pres = AlgebraPresentation(parent, (relation,), ('quadratic',))
    action = GeneratorAction(field, (((1, field.one),), ((0, -c0), (1, -c1))))
    proof = DerivationGraph(pres.content_hash(), (parent.zero, relation),
                            ((), (DerivationTerm(field.one, (), 'input', 0, ()),)))
    cert = WordBasisCertificate(pres.content_hash(), ((), (1,)), (action,), proof, ((0, 1),))
    algebra = verify_word_basis(pres, cert)
    result = minimal_polynomial(algebra, algebra.element((0, 1)), max_primes=3,
                                 reconstruction_degree=4, threads=2, start_prime=1, shift=(3, 5))
    assert result.coefficients == (c0, c1, field.one)
    assert result.kernel.kernel.field.parameters == ('c0', 'c1')


def test_native_bad_coefficient_prime_then_explicit_recovery(native_finiteflow):
    algebra = matrix_algebra()
    coefficient = Fraction(1, native_finiteflow.PrimeNo(0))
    element = algebra.scalar(coefficient)
    with pytest.raises(RuntimeError, match='FiniteFlow reconstruction failed'):
        minimal_polynomial(algebra, element, start_prime=0, max_primes=3)
    result = minimal_polynomial(algebra, element, start_prime=1, max_primes=3)
    assert result.coefficients == (-algebra.field.scalar(coefficient), algebra.field.one)
