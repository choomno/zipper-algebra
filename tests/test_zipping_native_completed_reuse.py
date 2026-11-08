"""Completed-subcomputation reuse retains exact full mathematical outputs."""
from itertools import product

import pytest

from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeError, NativeSession
from zipper_algebra.zipping_algorithm.output_normal_form import canonical_output

pytestmark = pytest.mark.timeout(300)


@pytest.fixture(scope="module")
def engine():
    with NativeSession() as session:
        yield session


def output(result):
    assert result["status"] == "complete", result["diagnostic"]
    return canonical_output(result["output"], n=result["n"], complete=True)


@pytest.mark.parametrize("n,length", [(2, 4), (3, 3), (4, 2)])
def test_all_short_words_match_standalone(engine, n, length):
    alphabet = tuple(range(-(n - 1), 0)) + tuple(range(1, n))
    for size in range(length + 1):
        for word in product(alphabet, repeat=size):
            reused = engine.reduce_reused(word, n=n)
            assert reused["schema"] == "zipper-native-reuse-v1"
            assert reused["trace_enabled"] is False
            assert output(reused) == output(engine.reduce(word, n=n))


@pytest.mark.parametrize("n,word", [
    (1, ()), (8, (6, -7, 6)), (5, (1, -2, 3, -4)),
    (3, (1, 2, -1)), (4, (3, -2, 1, 2)),
])
def test_recursive_contexts_labels_frames_and_full_trace(engine, n, word):
    before = engine.reduce(word, n=n, trace=True)
    first = engine.reduce_reused(word, n=n)
    second = engine.reduce_reused(word, n=n)
    after = engine.reduce(word, n=n, trace=True)
    assert before == after
    assert output(first) == output(second) == output(before)
    assert second["metrics"]["transition_evaluations"] == 0
    assert second["metrics"]["completed_cache_hits"] == 1


@pytest.mark.parametrize("entries,bytes_", [(0, 0), (2, 20000), (1, 1)])
def test_disabled_evicted_and_oversized_completed_caches(entries, bytes_):
    with NativeSession(cache_entries=entries, cache_bytes=bytes_) as session:
        for word in [(1, 2), (2, -1), (1, 2)]:
            actual = session.reduce_reused(word, n=3)
            assert output(actual) == output(session.reduce(word, n=3))
        info = session.cache_info()["completed"]
        assert info["entries"] <= entries
        assert info["bytes"] <= bytes_
        if entries == 2:
            assert info["evictions"] > 0
        if entries == 1:
            assert info["oversized"] > 0
        session.clear_cache()
        assert session.cache_info()["completed"]["entries"] == 0
        assert session.cache_info()["completed"]["hits"] == 0


@pytest.mark.parametrize("limits", [
    Limits(max_steps=1), Limits(max_terms=1), Limits(max_word_size=1),
    Limits(max_coefficient_support=1), Limits(expansion_budget=1),
    Limits(max_integer_bits=1),
])
def test_partial_dag_never_claims_complete_output(engine, limits):
    word = (1, 2, -1)
    engine.reduce_reused(word, n=3)  # A complete result for different limits cannot satisfy this request.
    result = engine.reduce_reused(word, n=3, limits=limits)
    assert result["status"] == "incomplete"
    assert result["output"] == []
    assert result["pending"] == []
    assert result["diagnostic"]["scope"] == "reuse-evaluator"


@pytest.mark.parametrize("word,n,limits", [
    ((), True, None), ((), 0, None), ((True,), 2, None), ((0,), 2, None),
    ((2,), 2, None), (iter([1]), 2, None), ((1,), 2, {}),
    ((), 2**63, None), ((), 2, Limits(max_steps=2**63)),
])
def test_reused_wrapper_rejects_invalid_inputs(engine, word, n, limits):
    with pytest.raises(NativeError):
        engine.reduce_reused(word, n=n, limits=limits)


@pytest.mark.parametrize('entries,bytes_', [(1, 100000), (7, 16000), (10000, 16000000)])
def test_shared_payload_accounting_survives_last_alias_eviction(entries, bytes_):
    """Eviction can remove a result's final alias just before its parent retains it."""
    with NativeSession(cache_entries=entries, cache_bytes=bytes_) as session:
        for word in [(1, 2, -1), (2, -1), (1, 2, -1)]:
            result = session.reduce_reused(word, n=3)
            assert output(result) == output(session.reduce(word, n=3))
            info = session.cache_info()['completed']
            assert info['bytes'] == info['key_bytes'] + info['unique_payload_bytes']
            assert info['bytes'] <= bytes_
            assert 0 < info['unique_payloads'] <= info['entries'] <= entries
            if entries == 1:
                assert info['unique_payloads'] == info['entries'] == 1
            if entries == 10000:
                assert info['unique_payloads'] < info['entries']
        session.clear_cache()
        info = session.cache_info()['completed']
        assert all(info[k] == 0 for k in ('bytes', 'key_bytes', 'unique_payload_bytes',
                                        'unique_payloads', 'entries'))
        assert output(session.reduce_reused((1, 2), n=3)) == output(session.reduce((1, 2), n=3))
