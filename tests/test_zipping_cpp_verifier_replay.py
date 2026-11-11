"""Full deterministic replay parity, including bounded atomic rollback and cache eviction."""
from dataclasses import asdict

import pytest

from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier, _same, _wire

pytestmark = pytest.mark.timeout(300)


@pytest.fixture(scope="module")
def producer_metadata():
    with NativeSession() as session:
        yield session.metadata


def assert_replay(verifier, reference, word, n, limits, batch, trace):
    expected = _wire(reference.reduce(word, n=n, limits=limits, batch_size=batch, trace=trace))
    actual = verifier.diagnostic("replay", n=n, input_word=word, limits=asdict(limits),
                                 batch_size=batch, trace_enabled=trace)
    assert _same(actual, expected)
    return actual


@pytest.mark.parametrize("trace", [False, True])
def test_all_ranks_signed_unreduced_and_batch_sizes(producer_metadata, trace):
    reference = CollectedSession()
    with CppVerifier(producer_metadata=producer_metadata) as verifier:
        for n in range(2, 9):
            for word in ([], [1, -1], [n - 1, 1, -(n - 1), -1]):
                for batch in (1, 4, 64):
                    assert_replay(verifier, reference, word, n, Limits(max_steps=40), batch, trace)


@pytest.mark.parametrize("field", list(asdict(Limits())))
@pytest.mark.parametrize("trace", [False, True])
def test_each_resource_bound_preserves_atomic_rollback(producer_metadata, field, trace):
    reference = CollectedSession()
    values = asdict(Limits(max_steps=40))
    values[field] = 1
    limits = Limits(**values)
    with CppVerifier(producer_metadata=producer_metadata) as verifier:
        result = assert_replay(verifier, reference, [1, -2, 1, 2, -1], 3, limits, 4, trace)
        assert result["status"] == "incomplete"
        assert result["diagnostic"]["batch_rolled_back"] is True
        assert result["pending"]
        if field == "max_steps":
            assert result["metrics"]["steps"] == 1
            assert result["metrics"]["attempted_steps"] == 1


@pytest.mark.parametrize("cache_entries,cache_bytes", [(0, 0), (1, 1000000), (10000, 10000000)])
def test_cold_warm_and_evicted_cache_match_full_reference(producer_metadata, cache_entries, cache_bytes):
    reference = CollectedSession()
    with CppVerifier(producer_metadata=producer_metadata, cache_entries=cache_entries,
                     cache_bytes=cache_bytes) as verifier:
        for _ in range(2):
            for trace in (False, True):
                for word in ([1, -2, 1], [1, -1], [-2, 1, -2]):
                    assert_replay(verifier, reference, word, 3, Limits(max_steps=100), 4, trace)
        info = verifier.cache_info()
        if cache_entries == 0:
            assert info["entries"] == info["hits"] == 0
        elif cache_entries == 1:
            assert info["entries"] == 1
            assert info["evictions"] > 0
        else:
            assert info["cross_input_hits"] > 0
        verifier.clear_cache()
        assert verifier.cache_info()["entries"] == 0
        assert_replay(verifier, reference, [1, -2, 1], 3, Limits(max_steps=100), 4, True)
