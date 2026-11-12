"""Exact normalization for working coefficients, preserving literal inputs."""
import builtins
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
import json
import subprocess
import sys
from threading import Event

import pytest

from zipper_algebra.fields import FieldElement, RationalFunctionField
from zipper_algebra.verification.working_field import WorkingField


def test_rational_working_scalars_use_fractions_without_symbolic_arithmetic():
    field = RationalFunctionField()
    working = WorkingField(field)
    assert type(working.zero) is F and working.zero == 0
    assert type(working.one) is F and working.one == 1
    assert working.from_element(field.scalar(F(2, 3))) == F(2, 3)
    assert working.from_element(field.zero) == 0
    assert working.to_element(F(7, 11)) == F(7, 11)
    assert working.to_element(7) == 7
    assert working.to_element(working.from_element(field.scalar(-3))) == -3
    left = working.from_element(field.scalar(F(1, 3)))
    right = working.from_element(field.scalar(F(1, 6)))
    assert working.to_element(left + right) == F(1, 2)


def test_ground_prime_field_keeps_exact_residue_arithmetic():
    field = RationalFunctionField((), 7)
    working = WorkingField(field)
    assert isinstance(working.zero, FieldElement) and working.zero == 0
    assert isinstance(working.one, FieldElement) and working.one == 1
    scalar = field.scalar(F(1, 2))
    assert working.from_element(scalar) is scalar
    assert working.to_element(scalar) is scalar
    assert working.to_element(working.one / scalar) == 2
    assert working.to_element(F(1, 2)) == 4
    assert working.to_element(9) == 2
    with pytest.raises(ZeroDivisionError):
        working.to_element(F(1, 7))


@pytest.mark.parametrize('field', [None, 'Q', 0, {}])
def test_working_field_requires_a_validated_domain(field):
    with pytest.raises(TypeError, match='declared'):
        WorkingField(field)


@pytest.mark.parametrize('cache_size', [-1, True, 2.0, None, '2'])
def test_working_field_cache_size_requires_an_exact_nonnegative_integer(cache_size):
    with pytest.raises(ValueError, match='cache size'):
        WorkingField(RationalFunctionField(), cache_size=cache_size)


@pytest.mark.parametrize('value', [True, 0.5, '1/2', None, object()])
@pytest.mark.parametrize('characteristic', [0, 7])
def test_native_ground_output_rejects_inexact_or_unsupported_values(value, characteristic):
    with pytest.raises(TypeError):
        WorkingField(RationalFunctionField((), characteristic)).to_element(value)


@pytest.mark.parametrize('value', [F(1, 2), 1, None, 'q', 0.5])
def test_input_requires_an_exact_public_coefficient_object(value):
    with pytest.raises(TypeError, match='different working field'):
        WorkingField(RationalFunctionField()).from_element(value)


def test_ground_working_fields_reject_mixed_characteristics_and_parameter_parents():
    rational = RationalFunctionField()
    prime = RationalFunctionField((), 7)
    parameterized = RationalFunctionField(('s',))
    for working, foreign in [(WorkingField(rational), prime.zero),
                             (WorkingField(prime), rational.zero),
                             (WorkingField(rational), parameterized.zero)]:
        with pytest.raises(TypeError):
            working.from_element(foreign)
        with pytest.raises(TypeError):
            working.to_element(foreign)


@pytest.mark.parametrize('characteristic', [0, 5])
def test_parameterized_normalization_cancels_gcd_but_preserves_literal_input(characteristic):
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s', 't'), characteristic)
    s, t = (field.parameter(name) for name in field.parameters)
    expression = (s**2 - t**2) / (s - t)
    original = json.dumps(expression.to_json(), sort_keys=True)
    working = WorkingField(field)
    native = working.from_element(expression)
    normalized = working.to_element(native)
    assert normalized == s + t
    assert normalized.denominator == field.one.denominator
    assert working.to_element(working.zero) == field.zero
    assert working.to_element(working.one) == field.one
    assert native == working.from_element(s + t)
    assert original == json.dumps(expression.to_json(), sort_keys=True)
    with pytest.raises(ZeroDivisionError):
        expression.specialize({'s': 2, 't': 2})
    assert normalized.specialize({'s': 2, 't': 2}) == field.scalar(4).specialize({'s': 2, 't': 2})
    # Arithmetic stays normalized through repeated exact multiply/divide steps.
    sn = working.from_element(s)
    tn = working.from_element(t)
    for _ in range(8):
        native = (native * sn - native * tn) / (sn - tn)
    assert working.to_element(native) == s + t


@pytest.mark.parametrize('characteristic', [0, 7])
def test_parameterized_rational_coefficients_and_parameter_order_round_trip(characteristic):
    pytest.importorskip('sympy')
    field = RationalFunctionField(('z', 'a', 'q'), characteristic)
    z, a, q = (field.parameter(name) for name in field.parameters)
    working = WorkingField(field)
    value = (F(2, 3)*z**2 + a - 2*q) / (q*a + 1)
    result = working.to_element(working.from_element(value))
    assert result == value
    assert result.field.parameters == ('z', 'a', 'q')
    assert result.specialize({'z': 2, 'a': 3, 'q': 4}) == value.specialize({'z': 2, 'a': 3, 'q': 4})
    assert working.to_element(F(1, 2)) == field.scalar(F(1, 2))


def test_parameter_names_are_symbols_without_expression_parsing():
    pytest.importorskip('sympy')
    names = ('lambda', 'I', 'E', 'pi', 'True', 'q_2')
    field = RationalFunctionField(names)
    working = WorkingField(field)
    coefficients = {name: i + 2 for i, name in enumerate(names)}
    for name in names:
        parameter = field.parameter(name)
        assert working.to_element(working.from_element(parameter)) == parameter
    value = sum((field.parameter(name) for name in names), field.zero)
    assert working.to_element(working.from_element(value)).specialize(coefficients) == sum(coefficients.values())


def test_parameterized_input_rejects_wrong_parameter_order_and_characteristic():
    pytest.importorskip('sympy')
    field = RationalFunctionField(('x', 'y'))
    working = WorkingField(field)
    for foreign_field in (RationalFunctionField(('y', 'x')),
                          RationalFunctionField(('x', 'y'), 5),
                          RationalFunctionField(('x',))):
        with pytest.raises(TypeError):
            working.from_element(foreign_field.zero)
        foreign = WorkingField(foreign_field)
        with pytest.raises(TypeError):
            working.to_element(foreign.from_element(foreign_field.zero))
    with pytest.raises(TypeError):
        working.from_element(RationalFunctionField().zero)
    # Native Fractions explicitly embed as constant coefficients.
    assert working.to_element(F(0)) == field.zero


@pytest.mark.parametrize('value', [True, 0.5, 's', None, object()])
def test_parameterized_output_requires_a_native_exact_value(value):
    pytest.importorskip('sympy')
    with pytest.raises(TypeError):
        WorkingField(RationalFunctionField(('s',))).to_element(value)


def test_native_conversion_cache_uses_literal_tuples_without_hashing_public_values():
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s',))
    working = WorkingField(field)
    value = (field.parameter('s')**2 - 1) / (field.parameter('s') - 1)
    same_literal = FieldElement.from_json(value.to_json())
    first = working.from_element(value)
    second = working.from_element(same_literal)
    assert first is second
    assert working.to_element(first) == field.parameter('s') + 1


def test_bounded_cache_evicts_least_recently_used_values_without_changing_arithmetic():
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s',))
    s = field.parameter('s')
    working = WorkingField(field, cache_size=2)
    first = working.from_element(s + 1)
    second = working.from_element(s + 2)
    assert working.from_element(s + 1) is first
    third = working.from_element(s + 3)
    assert len(working._cache) == 2
    assert working.from_element(s + 1) is first
    assert working.to_element(first + second - third) == s
    reloaded = working.from_element(s + 2)
    assert reloaded is not second
    assert working.to_element(reloaded) == s + 2
    assert len(working._cache) == 2
    for integer in range(10):
        assert working.to_element(working.from_element(s + integer)) == s + integer
        assert len(working._cache) <= 2


def test_conversion_cache_can_be_disabled_without_changing_exact_results():
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s',))
    working = WorkingField(field, cache_size=0)
    s = field.parameter('s')
    native = working.from_element((s**2 - 1)/(s - 1))
    again = working.from_element(s + 1)
    assert working.to_element(native - again) == 0
    assert working._cache == {}


def test_prime_field_nonmonic_native_fractions_compare_by_exact_residual():
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s',), 5)
    working = WorkingField(field)
    s = working.from_element(field.parameter('s'))
    first = s / (2*s + 1)
    second = 3*s / (s + 3)
    # SymPy can retain different scalar multiples of numerator/denominator in
    # a finite field. Subtract before checking zero in the working arithmetic.
    assert first - second == 0
    assert working.to_element(first) == working.to_element(second)
    assert working.to_element(first).specialize({'s': 1}) == 2


def test_parameter_field_constructor_reports_the_optional_dependency(monkeypatch):
    original = builtins.__import__

    def block_sympy(name, *args, **kwargs):
        if name == 'sympy':
            raise ImportError('symbolic dependency unavailable')
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', block_sympy)
    WorkingField(RationalFunctionField())
    WorkingField(RationalFunctionField((), 5))
    with pytest.raises(ImportError, match=r'zipper-algebra\[symbolic\]'):
        WorkingField(RationalFunctionField(('s',)))


def test_module_and_ground_working_fields_do_not_import_sympy(tmp_path):
    result = subprocess.run([sys.executable, '-I', '-c', '''
import sys
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.verification.working_field import WorkingField
WorkingField(RationalFunctionField())
WorkingField(RationalFunctionField((), 5))
assert 'sympy' not in sys.modules
print('exact ground fields available')
'''], cwd=tmp_path, text=True, capture_output=True, check=True, timeout=10)
    assert result.stdout.strip() == 'exact ground fields available'


def test_parameter_cache_hit_survives_concurrent_eviction():
    pytest.importorskip('sympy')
    field = RationalFunctionField(('s',))
    working = WorkingField(field, cache_size=1)
    hit_started = Event()
    eviction_finished = Event()

    class PausingCache(OrderedDict):
        def move_to_end(self, key, last=True):
            hit_started.set()
            # Let a concurrent conversion attempt to evict this hit.
            eviction_finished.wait(timeout=0.5)
            super().move_to_end(key, last)

    working._cache = PausingCache()
    parameter = field.parameter('s')
    expected = working.from_element(parameter)
    with ThreadPoolExecutor(max_workers=2) as pool:
        hit = pool.submit(working.from_element, parameter)
        assert hit_started.wait(timeout=5)
        eviction = pool.submit(working.from_element, field.one)
        eviction.add_done_callback(lambda _: eviction_finished.set())
        assert hit.result(timeout=5) is expected
        assert working.to_element(eviction.result(timeout=5)) == field.one
