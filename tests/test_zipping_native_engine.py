"""Differential checks of the retained C++ state engine against Python replay."""
from copy import deepcopy
from itertools import product

import pytest

from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession, reference_certificate


@pytest.fixture(scope="module")
def engine():
    with NativeSession() as session:
        yield session


def as_arrays(value):
    if isinstance(value, (list, tuple)):
        return [as_arrays(item) for item in value]
    if isinstance(value, dict):
        return {key: as_arrays(item) for key, item in value.items()}
    return value


def first_difference(actual, expected, path="$ "):
    if type(actual) is not type(expected):
        return f"{path}: type {type(actual).__name__} != {type(expected).__name__}"
    if isinstance(actual, dict):
        if actual.keys() != expected.keys():
            return f"{path}: keys {set(actual) ^ set(expected)}"
        for key in actual:
            difference = first_difference(actual[key], expected[key], f"{path}.{key}")
            if difference:
                return difference
    elif isinstance(actual, (list, tuple)):
        if len(actual) != len(expected):
            return f"{path}: lengths {len(actual)} != {len(expected)}"
        for index, (left, right) in enumerate(zip(actual, expected, strict=True)):
            difference = first_difference(left, right, f"{path}[{index}]")
            if difference:
                return difference
    elif actual != expected:
        if isinstance(actual, int) and max(actual.bit_length(), expected.bit_length()) > 1000:
            return f"{path}: different arbitrary-precision integers"
        return f"{path}: {actual!r} != {expected!r}"
    return None


def compare(engine, word, n, **options):
    oracle = CollectedSession()
    oracle._catalog.configure_rule_cache(max_bytes=0, max_entries=0)
    expected = oracle.reduce(word, n=n, **options)
    actual = engine.reduce(word, n=n, **options)
    difference = first_difference(reference_certificate(actual), as_arrays(expected))
    assert difference is None, (n, word, options, difference)
    return actual


@pytest.mark.parametrize("n", [1, 2, 3, 4, 8, 16])
def test_empty_literal_words_and_recursive_complements(engine, n):
    compare(engine, (), n, trace=True, batch_size=4)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_all_signed_literal_words_through_length_two(engine, n):
    alphabet = tuple(range(-(n - 1), 0)) + tuple(range(1, n))
    for length in range(3):
        for word in product(alphabet, repeat=length):
            compare(engine, word, n, batch_size=4)


@pytest.mark.parametrize(
    "n,word",
    [
        (2, (1, -1, 1, -1, -1, 1)),
        (3, (1, -2, 1, -2)),
        (3, (-1, -2, -1, -2)),
        (4, (1, 2, 3)),
        (5, (3, -4)),
        (8, (6, -7)),
    ],
)
@pytest.mark.parametrize("trace", [False, True])
def test_all_phases_mirrors_and_unreduced_words(engine, n, word, trace):
    compare(engine, word, n, trace=trace, batch_size=4)


@pytest.mark.parametrize(
    "word,n,limits",
    [
        ((1,), 2, Limits(max_terms=1)),
        ((1,), 2, Limits(max_coefficient_support=1)),
        ((1, 1), 2, Limits(max_integer_bits=1)),
        ((1,), 2, Limits(max_word_size=1)),
        ((1, 2), 3, Limits(max_word_size=5)),
        ((1, 2), 3, Limits(expansion_budget=1)),
        ((1, 2), 3, Limits(max_steps=2)),
        ((1, -2) * 5, 3, Limits(max_steps=128)),
        ((1, 2, 3, 4, 5, 6, 7), 8, Limits(max_steps=128)),
    ],
)
def test_atomic_limits_preserve_exact_pending_states_and_trace(engine, word, n, limits):
    result = compare(engine, word, n, limits=limits, trace=True)
    assert result["status"] == "incomplete"
    assert result["diagnostic"]["batch_rolled_back"] is True
    assert result["pending"]


@pytest.mark.parametrize("batch_size", [1, 3, 16, 64])
def test_batch_order_collection_and_cancellation(engine, batch_size):
    result = compare(engine, (1, 2, -1), 3, trace=True, batch_size=batch_size)
    assert result["metrics"]["merged_terms"] > 0


def test_repeated_native_calls_are_isolated_and_reproducible(engine):
    first = engine.reduce((1, 2), n=3, trace=True)
    engine.reduce((-1,), n=8)
    second = engine.reduce((1, 2), n=3, trace=True)
    assert first == second
    assert engine.check(first)
    changed = deepcopy(first)
    changed["output"][0]["coefficient"][0][1] += 1
    with pytest.raises(ValueError):
        engine.check(changed)


def test_native_cache_warm_cold_and_clear_preserve_certificate(engine):
    engine.clear_cache()
    assert engine.cache_info()["entries"] == 0
    cold = engine.reduce((1, 2), n=3, trace=True)
    cold_info = engine.cache_info()
    warm = engine.reduce((1, 2), n=3, trace=True)
    assert cold == warm
    assert engine.cache_info()["hits"] > cold_info["hits"]
    assert engine.cache_info()["evaluations"] == cold_info["evaluations"]
    engine.clear_cache()
    assert engine.cache_info()["entries"] == engine.cache_info()["bytes"] == 0
    assert engine.cache_info()["hits"] == engine.cache_info()["evaluations"] == 0
    assert engine.reduce((1, 2), n=3, trace=True) == cold


@pytest.mark.parametrize("configuration", [{"cache_entries": 0}, {"cache_bytes": 0},
                                            {"cache_entries": 1}, {"cache_bytes": 1}])
def test_disabled_evicted_or_oversized_cache_preserves_exact_result(configuration):
    with NativeSession(**configuration) as session:
        compare(session, (1, 2), 3, trace=True)
        info = session.cache_info()
        if not configuration.get("cache_entries", 1) or not configuration.get("cache_bytes", 1):
            assert info["entries"] == info["bytes"] == info["hits"] == 0
        elif configuration.get("cache_entries") == 1:
            assert info["entries"] <= 1 and info["evictions"] > 0
        else:
            assert info["entries"] == 0 and info["oversized"] > 0


def test_native_cache_keys_retain_limits_and_trace_policy(engine):
    engine.clear_cache()
    compare(engine, (1, 2), 3, trace=False)
    compare(engine, (1, 2), 3, trace=True)
    compare(engine, (1, 2), 3, trace=True, limits=Limits(expansion_budget=1))
    compare(engine, (1, 2), 3, trace=True, limits=Limits(max_steps=2))


def test_native_sparse_coefficients_above_machine_width_in_every_phase(engine):
    from dataclasses import asdict
    from zipper_algebra.zipping_algorithm import collected
    from zipper_algebra.zipping_algorithm.native import _decode_coefficients, _encode_coefficients

    oracle = CollectedSession()
    run = oracle.reduce((1, 2), n=3, trace=True)
    states = {}
    for batch in run["transitions"]:
        for state in batch["predecessors"]:
            states.setdefault(state["phase"], state)
    assert set(states) == {"single", "zip_up", "zip_down", "layer", "extract"}
    scalar = (1 << 32770) + (1 << 511) + 19
    for original in states.values():
        work = deepcopy(original)
        collected._set_coefficient(work, [[-17, -scalar], [0, scalar + 1], [23, -3]])
        expected = collected._transition(work, Limits(), True, oracle._catalog)
        request = {"op": "transition", "work": work, "limits": asdict(Limits()), "trace": True}
        cold = _decode_coefficients(engine._request(_encode_coefficients(request)))
        warm = _decode_coefficients(engine._request(_encode_coefficients(request)))
        assert cold == warm == as_arrays(expected), work["phase"]


def test_reduce_does_not_call_python_transition_or_collection(engine, monkeypatch):
    from zipper_algebra.zipping_algorithm import algorithm, collected

    expected = CollectedSession().reduce((1, 2), n=3, trace=True)
    def forbidden(*args, **kwargs):
        raise AssertionError("native reduction called the Python state engine")
    monkeypatch.setattr(collected, "_transition", forbidden)
    monkeypatch.setattr(collected, "_stage", forbidden)
    for name in ("_single_step", "_zip_step", "_layer_step", "_extract_step"):
        monkeypatch.setattr(algorithm, name, forbidden)
    actual = engine.reduce((1, 2), n=3, trace=True)
    assert first_difference(reference_certificate(actual), as_arrays(expected)) is None


@pytest.mark.parametrize("coordinate", ["q_exponent", "depth"])
def test_native_coordinate_overflow_is_a_bounded_refusal(engine, coordinate):
    from dataclasses import asdict

    reference = CollectedSession().reduce((1,), n=2, trace=True)
    work = deepcopy(reference["transitions"][0]["predecessors"][0])
    if coordinate == "q_exponent":
        work["coefficient"] = [[2**63 - 1, 1]]
    else:
        work["depth"] = 2**63 - 1
    result = engine._request({"op": "transition", "work": work,
                              "limits": asdict(Limits()), "trace": False})
    assert result == {"error": "incomplete", "message": "signed integer coordinate overflow"}
