"""Independent arithmetic and rejection controls for exact coefficient fields."""
from copy import deepcopy
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
import json
import operator

import pytest

from zipper_algebra.fields import FieldElement, RationalFunctionField


@pytest.mark.timeout(2)
def test_prime_specialization_uses_modular_power_for_large_sparse_exponents():
    for characteristic in (0, 5):
        field = RationalFunctionField(('s',), characteristic)
        huge = field.polynomial({(10**9,): 1})
        # The exponent is divisible by four, and 2 has order four modulo five.
        assert huge.specialize({'s': 2}, characteristic=5) == 1


def test_domains_are_immutable_and_parameter_order_is_mathematical_data():
    field = RationalFunctionField(('q', 'a', 'b'))
    assert field.parameters == ('q', 'a', 'b')
    assert field.characteristic == 0
    assert field == RationalFunctionField(('q', 'a', 'b'))
    assert field != RationalFunctionField(('a', 'q', 'b'))
    assert hash(field) == hash(RationalFunctionField(('q', 'a', 'b')))
    with pytest.raises(FrozenInstanceError):
        field.characteristic = 3
    assert RationalFunctionField().scalar(F(2, 3)).numerator == (((), F(2, 3)),)


@pytest.mark.parametrize('parameters', [
    ['q'], 'q', (1,), ('',), ('x+y',), ('x-y',), ('1q',), ('α',), ('q', 'q'),
])
def test_bad_parameter_declarations_are_rejected(parameters):
    with pytest.raises((TypeError, ValueError)):
        RationalFunctionField(parameters)


@pytest.mark.parametrize('characteristic', [True, 1, -1, 4, 9, 15, 2.0, '3'])
def test_composite_and_inexact_characteristics_are_rejected(characteristic):
    with pytest.raises((TypeError, ValueError)):
        RationalFunctionField((), characteristic)


def test_valid_prime_fields_use_exact_modular_rationals():
    field = RationalFunctionField(('x',), 7)
    x = field.parameter('x')
    assert field.scalar(F(1, 2)).numerator == (((0,), 4),)
    assert (x + 10).numerator == (((0,), 3), ((1,), 1))
    assert (x * 7).is_zero()
    assert field.polynomial({(2,): 14}).is_zero()
    assert field.scalar(-1) == 6
    assert field.scalar(F(5, 2)) == 6
    assert field.scalar(F(3, 2)) * 2 == 3
    with pytest.raises(ZeroDivisionError):
        field.scalar(F(1, 7))


def test_sparse_input_is_an_immutable_sorted_snapshot():
    field = RationalFunctionField(('x', 'y'))
    source = {(0, 1): 2, (1, 0): F(1, 3), (0, 0): 0}
    value = field.polynomial(source)
    source[(0, 1)] = 99
    assert value.numerator == (((0, 1), F(2)), ((1, 0), F(1, 3)))
    assert value.denominator == (((0, 0), F(1)),)
    assert field.scalar(value) is value
    assert value.field is field
    with pytest.raises(FrozenInstanceError):
        value.numerator = ()
    with pytest.raises(TypeError):
        hash(value)


@pytest.mark.parametrize('bad', [
    [], {(1,): 1}, {(0, -1): 1}, {(0, True): 1}, {'x': 1},
    {(0, 0): True}, {(0, 0): 0.5}, {(0, 0): '1/2'}, {(0, 0): None},
])
def test_sparse_polynomial_input_requires_exact_coefficients_and_exponents(bad):
    with pytest.raises((TypeError, ValueError)):
        RationalFunctionField(('x', 'y')).polynomial(bad)


def test_constructor_and_unknown_parameters_validate_the_domain():
    field = RationalFunctionField(('q',))
    with pytest.raises(TypeError):
        FieldElement(None, {(0,): 1})
    with pytest.raises(ValueError):
        field.parameter('absent')
    with pytest.raises(ZeroDivisionError):
        FieldElement(field, {(0,): 1}, {})
    value = FieldElement(field, {(0,): 2}, {(0,): 4})
    assert value == F(1, 2)
    assert value.denominator == (((0,), F(1)),)
    assert FieldElement(field, {}, {(1,): 1}).denominator == field.one.denominator


def test_independent_three_parameter_polynomial_identity_and_literal_coefficients():
    field = RationalFunctionField(('q', 'a', 'b'))
    q, a, b = (field.parameter(name) for name in field.parameters)
    product = (q + a + b) * (q + a - b)
    # Expansion is specified independently, not by production multiplication.
    assert product.numerator == (
        ((0, 0, 2), F(-1)), ((0, 2, 0), F(1)),
        ((1, 1, 0), F(2)), ((2, 0, 0), F(1)),
    )
    assert product == q**2 + 2*q*a + a**2 - b**2
    fraction = (q + a) / (q - b)
    assert fraction * (q - b) == q + a
    assert fraction + (a - q) / (q - b) == 2*a / (q - b)
    assert fraction - fraction == 0
    assert (q / a) + (a / q) == (q**2 + a**2) / (q*a)
    assert (2 - q) + (q - 2) == 0
    assert q + 0 == q and 0 + q == q
    assert field.zero + q == q
    assert 1*q == q and q*1 == q and q*0 == 0 and 0*q == 0
    assert q / 2 == F(1, 2)*q
    assert 2/q == field.scalar(2) / q
    assert -(-q) == q


def test_prime_field_frobenius_uses_the_declared_characteristic():
    field = RationalFunctionField(('x', 'y'), 2)
    x, y = (field.parameter(name) for name in field.parameters)
    assert (x + y)**2 == x**2 + y**2
    assert (x + y)**2 != x**2 + 2*x*y + y**2 + 1
    assert x - y == x + y
    assert (x / (y + 1))**2 == x**2 / (y**2 + 1)


def test_integer_powers_and_field_inverse_are_exact():
    field = RationalFunctionField(('s',))
    s = field.parameter('s')
    assert s**0 == 1
    assert s**1 == s
    assert (s + 1)**3 == field.polynomial({(0,): 1, (1,): 3, (2,): 3, (3,): 1})
    assert s**-3 * s**3 == 1
    assert (s + 1)**-2 == 1 / (s + 1)**2
    assert field.zero**0 == 1
    assert field.zero / s == 0
    for zero in (0, field.zero):
        with pytest.raises(ZeroDivisionError):
            s / zero
        with pytest.raises(ZeroDivisionError):
            field.zero / zero
    with pytest.raises(ZeroDivisionError):
        field.zero**-1
    for exponent in (True, F(1, 2), 2.0):
        with pytest.raises(TypeError):
            s**exponent


@pytest.mark.parametrize('operation', [operator.add, operator.sub, operator.mul, operator.truediv])
def test_arithmetic_rejects_different_full_domains_even_for_zero_operands(operation):
    fields = [RationalFunctionField(('x', 'y')), RationalFunctionField(('y', 'x')),
              RationalFunctionField(('x',)), RationalFunctionField(('x', 'y'), 3)]
    for other in fields[1:]:
        with pytest.raises(TypeError, match='fields'):
            operation(fields[0].zero, other.one)
        with pytest.raises(TypeError, match='fields'):
            operation(fields[0].one, other.zero)
        assert fields[0].zero != other.zero
        with pytest.raises(TypeError, match='fields'):
            fields[0].scalar(other.zero)


@pytest.mark.parametrize('value', [True, 0.5, '2', None, object()])
def test_unsupported_arithmetic_operands_are_not_guessed(value):
    field = RationalFunctionField(('s',))
    s = field.parameter('s')
    assert s != value
    with pytest.raises(TypeError):
        field.scalar(value)
    for operation in (operator.add, operator.sub, operator.mul, operator.truediv):
        with pytest.raises(TypeError):
            operation(s, value)
        with pytest.raises(TypeError):
            operation(value, s)


def test_cross_multiplication_equality_retains_the_stored_exceptional_locus():
    field = RationalFunctionField(('s',))
    s = field.parameter('s')
    unreduced = (s**2 - 1) / (s - 1)
    assert unreduced == s + 1
    assert unreduced.to_json() != (s + 1).to_json()
    assert unreduced.specialize({'s': 2}) == 3
    assert (s + 1).specialize({'s': 1}) == 2
    with pytest.raises(ZeroDivisionError):
        unreduced.specialize({'s': 1})


def test_zero_cancellation_changes_the_retained_expression_specialization_locus():
    field = RationalFunctionField(('s',))
    summand = 1 / (field.parameter('s') - 1)
    with pytest.raises(ZeroDivisionError):
        summand.specialize({'s': 1})
    cancelled = summand - summand
    assert cancelled.numerator == ()
    assert cancelled.denominator == field.one.denominator
    # Specialization checks the stored zero/one fraction. Arithmetic history
    # does not supply an additional excluded locus for the simplified value.
    assert cancelled.specialize({'s': 1}) == 0


def test_specialization_has_explicit_exact_domain_and_prime_reduction():
    field = RationalFunctionField(('q', 'a', 'b'))
    q, a, b = (field.parameter(name) for name in field.parameters)
    value = (q**2 + a*b + F(1, 2)) / (q - a)
    point = {'b': 5, 'q': 3, 'a': 2}
    assert value.specialize(point) == F(39, 2)
    assert value.specialize(point).field == RationalFunctionField()
    residue = value.specialize(point, characteristic=7)
    assert residue == 2
    assert residue.field == RationalFunctionField((), 7)
    with pytest.raises(ZeroDivisionError):
        value.specialize(point, characteristic=2)
    finite = RationalFunctionField(('q',), 7)
    assert finite.parameter('q').specialize({'q': F(1, 2)}) == 4
    assert RationalFunctionField().scalar(3).specialize({}) == 3
    with pytest.raises(ValueError):
        finite.one.specialize({'q': 1}, characteristic=0)
    with pytest.raises(ValueError):
        finite.one.specialize({'q': 1}, characteristic=3)


@pytest.mark.parametrize('point', [None, [], {}, {'s': 2, 'x': 3}, {'x': 2}, {'s': 0.5}, {'s': True}])
def test_specialization_rejects_incomplete_wrong_or_inexact_points(point):
    with pytest.raises((TypeError, ValueError)):
        RationalFunctionField(('s',)).parameter('s').specialize(point)


def test_field_and_element_json_round_trips_have_exact_literal_data():
    field = RationalFunctionField(('q', 'a'))
    q, a = (field.parameter(name) for name in field.parameters)
    expression = (q + F(2, 3)*a) / (q - 1)
    declaration = {'schema': 'zipper-field-v1', 'characteristic': 0, 'parameters': ['q', 'a']}
    assert field.to_json() == declaration
    assert RationalFunctionField.from_json(json.loads(json.dumps(declaration))) == field
    encoded = {'schema': 'zipper-field-element-v1', 'field': declaration,
               'numerator': [[[0, 1], '2/3'], [[1, 0], 1]],
               'denominator': [[[0, 0], -1], [[1, 0], 1]]}
    assert expression.to_json() == encoded
    loaded = FieldElement.from_json(json.loads(json.dumps(encoded)))
    assert loaded == expression
    assert loaded.to_json() == encoded
    assert FieldElement.from_json(field.zero.to_json()).is_zero()
    prime = RationalFunctionField(('x',), 5)
    assert FieldElement.from_json((prime.parameter('x') / 2).to_json()) == 3*prime.parameter('x')


@pytest.mark.parametrize('value', [
    None, [], {}, {'schema': 'unknown', 'characteristic': 0, 'parameters': []},
    {'schema': 'zipper-field-v1', 'characteristic': 0, 'parameters': ()},
    {'schema': 'zipper-field-v1', 'characteristic': 0, 'parameters': [], 'extra': 1},
    {'schema': 'zipper-field-v1', 'parameters': []},
])
def test_malformed_field_declarations_are_rejected(value):
    with pytest.raises(ValueError):
        RationalFunctionField.from_json(value)


@pytest.mark.parametrize('replacement', [
    None, (), 'x', [1], [[[], 1]], [[[-1], 1]], [[[True], 1]],
    [[[0], 1], [[0], 2]], [[[0], True]], [[[0], 0.5]], [[[0], '1/0']],
    [[[0], '0.5']], [[[0], '1e2']], [[[0], '1/2/3']], [[[0], None]],
    [[(0,), 1]], [[[0], 1, 2]],
])
def test_malformed_sparse_json_is_rejected(replacement):
    encoded = RationalFunctionField(('s',)).one.to_json()
    encoded['numerator'] = replacement
    with pytest.raises((TypeError, ValueError)):
        FieldElement.from_json(encoded)


def test_element_json_rejects_wrong_schema_keys_and_denominators():
    valid = RationalFunctionField(('s',)).one.to_json()
    invalid = [None, [], {}, dict(valid, schema='wrong'), dict(valid, extra=1)]
    missing = deepcopy(valid)
    del missing['numerator']
    invalid.append(missing)
    for value in invalid:
        with pytest.raises(ValueError):
            FieldElement.from_json(value)
    zero_denominator = dict(valid, denominator=[])
    with pytest.raises(ZeroDivisionError):
        FieldElement.from_json(zero_denominator)


def test_text_is_deterministic_and_not_an_implicit_parser():
    field = RationalFunctionField(('q', 's'))
    q, s = (field.parameter(name) for name in field.parameters)
    assert str(field.zero) == '0'
    assert str(field.scalar(F(2, 3))) == '2/3'
    assert str(q**2 + s) == '1*s+1*q^2'
    assert str(q/(s + 1)) == '(1*q)/(1+1*s)'
