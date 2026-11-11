"""Intrinsic transition reuse preserves contextual collected calculations."""

from copy import deepcopy
from itertools import product
from pathlib import Path

import pytest

from zipper_algebra.zipping_algorithm import collected, memoized, framed, algorithm
from tests.test_zipping_algorithm import _root_common
from zipper_algebra.zipping_algorithm.interpreter import Limits


def _mathematical(certificate):
    value = deepcopy(certificate)
    for name in ("schema", "scope", "source_sha256"):
        value.pop(name)
    value["execution"].pop("memoization", None)
    return value


def _states(word=(1, 2), n=3):
    run = collected.reduce_word(word, n=n, trace=True)
    return [state for batch in run["transitions"] for state in batch["predecessors"]]


def _worker_snapshot():
    assert collected._CATALOG is not None
    return collected._CATALOG.rule_cache_info(), collected.__file__


def test_all_three_strand_literal_inputs_through_length_two_match_collected_oracle():
    oracle = collected.CollectedSession()
    session = memoized.MemoizedSession()
    for length in range(3):
        for word in product((-2, -1, 1, 2), repeat=length):
            expected = oracle.reduce(word, n=3)
            actual = session.reduce(word, n=3)
            assert _mathematical(actual) == _mathematical(expected), word
    assert session.cache_info()["cross_input_hits"] > 0
    assert session.cache_info()["complement_hits"] > 0
    info = session.cache_info()
    assert sum(info["actual_unit_evaluations"].values()) == info["evaluations"]
    assert sum(info["unit_lru_hits"].values()) == info["hits"]
    assert sum(info["unit_batch_reuses"].values()) == info["batch_reuses"]


@pytest.mark.parametrize(
    "word,n", [((1, 2), 3), ((1, -2, 1, -2), 3), ((1, 2, 3), 4), ((3, -4), 5)]
)
@pytest.mark.parametrize("trace", [False, True])
def test_cold_warm_disabled_and_parallel_cache_return_identical_exact_certificates(
    word, n, trace
):
    oracle = collected.reduce_word(word, n=n, trace=trace, batch_size=4)
    session = memoized.MemoizedSession()
    cold = session.reduce(word, n=n, trace=trace, batch_size=4)
    before = session.cache_info()
    warm = session.reduce(word, n=n, trace=trace, batch_size=4)
    assert cold == warm
    assert session.cache_info()["hits"] > before["hits"]
    assert session.cache_info()["evaluations"] == before["evaluations"]
    disabled = memoized.MemoizedSession(cache_bytes=0, cache_entries=0).reduce(
        word, n=n, trace=trace, batch_size=4
    )
    assert cold == disabled
    assert _mathematical(cold) == _mathematical(oracle)


def test_spawned_miss_workers_and_serial_execution_are_identical():
    a = memoized.reduce_word((1, 2, 1), n=3, trace=True, batch_size=4)
    b = memoized.reduce_word((1, 2, 1), n=3, trace=True, batch_size=4, workers=2)
    assert a == b


@pytest.mark.parametrize("rule_bytes,rule_entries", [(0, 0), (4096, 2)])
def test_spawned_workers_use_worktree_sources_and_requested_family_cache_caps(
    rule_bytes, rule_entries
):
    session = memoized.MemoizedSession(
        rule_cache_bytes=rule_bytes, rule_cache_entries=rule_entries
    )
    with session._new_pool(2, 4) as pool:
        info, source = pool.submit(_worker_snapshot).result(timeout=20)
    assert info["max_bytes"] == rule_bytes and info["max_entries"] == rule_entries
    assert Path(source).resolve().is_relative_to(Path(__file__).resolve().parents[1])


def test_both_cache_layers_can_be_disabled_without_changing_certificate():
    session = memoized.MemoizedSession(
        cache_bytes=0, cache_entries=0, rule_cache_bytes=0, rule_cache_entries=0
    )
    value = session.reduce((1, 2), n=3, trace=True)
    assert value == memoized.reduce_word((1, 2), n=3, trace=True)
    info = session.cache_info()
    assert info["entries"] == info["bytes"] == 0
    assert info["rule_cache"]["entries"] == info["rule_cache"]["retained_bytes"] == 0


def test_shallow_key_and_intrinsic_views_never_mutate_caller_state():
    state = next(x for x in _states() if "branch" in x)
    before = deepcopy(state)
    collected.continuation_key(state)
    active = memoized._intrinsic(state)
    assert state == before
    assert active["branch"] is not state["branch"]
    memoized.MemoizedSession()._execute_batch(
        [(collected.continuation_key(state), state)], Limits(), True, 1, 64, None
    )
    assert state == before


@pytest.mark.parametrize("trace", [False, True])
def test_complement_templates_rebase_nonconsecutive_labels_contexts_and_signed_scalars(
    trace,
):
    states = _states((3, -4), n=5) + _states((-1, -2), n=3)
    # Include both extraction-to-complement and final extraction; at least one
    # caller has endpoint values outside its current local rank.
    chosen = [
        work for work in states if work["contexts"] and work["phase"] == "extract"
    ]
    assert any(max(work["top_labels"]) > work["rank"] for work in chosen)
    assert any(work["rank"] == 1 for work in chosen)
    assert any(work["rank"] > 1 for work in chosen)
    session = memoized.MemoizedSession()
    for original in chosen:
        work = deepcopy(original)
        collected._set_coefficient(work, [[-3, -2], [0, 3], [4, -1]])
        expected = collected._transition(work, Limits(), trace, session._catalog)
        selected = [(collected.continuation_key(work), work)]
        cold, pool = session._execute_batch(selected, Limits(), trace, 1, 64, None)
        warm, pool = session._execute_batch(selected, Limits(), trace, 1, 64, pool)
        assert pool is None
        assert cold == warm == [expected]
        if cold[0]["output"]:
            assert (
                cold[0]["output"][0]["components"][: len(work["contexts"])]
                == work["contexts"]
            )
    assert session.cache_info()["complement_hits"] > 0


def test_different_outer_contexts_share_intrinsic_complement_but_never_accumulators():
    states = _states()
    candidates = [x for x in states if x["contexts"] and x["phase"] == "single"]
    groups = {}
    pair = None
    for work in candidates:
        key = collected.continuation_key(memoized._intrinsic(work))
        if key in groups and groups[key]["contexts"] != work["contexts"]:
            pair = groups[key], work
            break
        groups[key] = work
    assert pair is not None
    session = memoized.MemoizedSession()
    for work in pair:
        expected = collected._transition(work, Limits(), True, session._catalog)
        result, _ = session._execute_batch(
            [(collected.continuation_key(work), work)], Limits(), True, 1, 64, None
        )
        assert result == [expected]
    info = session.cache_info()
    assert info["misses"] == info["evaluations"] == 1
    assert info["hits"] == info["complement_hits"] == 1
    # A returned result may be edited without changing the immutable cache.
    result[0]["children"][0]["contexts"].clear()
    work = pair[1]
    again, _ = session._execute_batch(
        [(collected.continuation_key(work), work)], Limits(), True, 1, 64, None
    )
    assert again == [collected._transition(work, Limits(), True, session._catalog)]


@pytest.mark.parametrize(
    "configuration",
    [
        {"cache_entries": 1},
        {"cache_bytes": 1},
        {"cache_entries": 0},
        {"cache_bytes": 0},
    ],
)
def test_eviction_oversized_and_disabled_cache_change_only_runtime_work(configuration):
    session = memoized.MemoizedSession(**configuration)
    actual = session.reduce((1, 2), n=3, trace=True)
    expected = collected.reduce_word((1, 2), n=3, trace=True)
    assert _mathematical(actual) == _mathematical(expected)
    info = session.cache_info()
    assert info["entries"] <= info["max_entries"]
    assert info["bytes"] <= info["max_bytes"]
    if configuration == {"cache_entries": 1}:
        assert info["evictions"] > 0
    if configuration == {"cache_bytes": 1}:
        assert info["oversized"] > 0 and info["entries"] == 0
    session.clear_cache()
    assert session.cache_info()["bytes"] == session.cache_info()["entries"] == 0
    assert session.cache_info()["hits"] == session.cache_info()["misses"] == 0


@pytest.mark.parametrize(
    "limits",
    [
        Limits(max_steps=15),
        Limits(max_word_size=5),
        Limits(max_terms=3),
        Limits(max_coefficient_support=1),
        Limits(max_integer_bits=1),
        Limits(expansion_budget=1),
    ],
)
def test_limits_and_refusal_checkpoints_match_even_with_warm_templates(limits):
    session = memoized.MemoizedSession()
    session.reduce((1, 2), n=3, trace=True)
    actual = session.reduce((1, 2), n=3, trace=True, limits=limits)
    expected = collected.reduce_word((1, 2), n=3, trace=True, limits=limits)
    assert _mathematical(actual) == _mathematical(expected)
    assert actual["status"] == "incomplete"
    assert session.reduce((1, 2), n=3, trace=True, limits=limits) == actual


def test_local_failures_are_not_retained_and_trace_policies_do_not_cross():
    session = memoized.MemoizedSession()
    limits = Limits(expansion_budget=1)
    a = session.reduce((1, 2), n=3, limits=limits)
    b = session.reduce((1, 2), n=3, limits=limits)
    assert a == b and session.cache_info()["entries"] == 0
    assert session.cache_info()["misses"] == session.cache_info()["evaluations"] == 2
    session.reduce((1, 2), n=3, trace=False)
    before = session.cache_info()["misses"]
    session.reduce((1, 2), n=3, trace=True)
    assert session.cache_info()["misses"] > before


def test_replay_is_cold_and_resume_and_tamper_validation_remain_exact():
    session = memoized.MemoizedSession()
    stopped = session.reduce((1, 2), n=3, trace=True, limits=Limits(max_steps=15))
    before = session.cache_info()
    assert session.check(stopped)
    assert session.cache_info() == before
    assert memoized.check_certificate(stopped)
    assert memoized.resume(stopped) == stopped
    assert memoized.resume(stopped, limits=Limits()) == memoized.reduce_word(
        (1, 2), n=3, trace=True
    )
    with pytest.raises(collected.CollectedError, match="must not decrease"):
        memoized.resume(stopped, limits=Limits(max_steps=14))
    altered = deepcopy(stopped)
    altered["pending"][0]["coefficient"] = [[0, 99]]
    with pytest.raises(collected.CollectedError, match="differs"):
        session.check(altered)


def test_repeated_checks_reuse_only_parsed_catalog_not_results(monkeypatch):
    session = memoized.MemoizedSession()
    value = session.reduce((1, 2), n=3, trace=True)
    loads = 0
    original = collected.Catalog.bundled

    def counted_catalog(cls):
        nonlocal loads
        loads += 1
        return original()

    monkeypatch.setattr(collected.Catalog, "bundled", classmethod(counted_catalog))
    assert session.check(value)
    verifier = session._checker
    for _ in range(3):
        assert session.check(value)
    assert loads == 1 and session._checker is verifier
    info = verifier.cache_info()
    assert info["max_bytes"] == info["max_entries"] == info["entries"] == 0
    assert info["rule_cache"]["max_bytes"] == info["rule_cache"]["max_entries"] == 0
    assert info["rule_cache"]["entries"] == 0
    # Replaying altered data still executes the full unfactored verifier.
    altered = deepcopy(value)
    altered["output"][0]["coefficient"][0][1] += 7
    with pytest.raises(collected.CollectedError, match="differs"):
        session.check(altered)


@pytest.mark.parametrize(
    "configuration",
    [
        {"cache_bytes": -1},
        {"cache_entries": -1},
        {"cache_bytes": True},
        {"cache_entries": 0.5},
    ],
)
def test_invalid_cache_budgets_rejected(configuration):
    with pytest.raises(collected.CollectedError, match="nonnegative integers"):
        memoized.MemoizedSession(**configuration)


@pytest.mark.parametrize("size", [2, 3, 5, 8, 9])
@pytest.mark.parametrize("mirrored", [False, True])
def test_parameterized_root_child_flypes_reuse_beyond_size_five(size, mirrored):
    n = size + 1
    incoming = (
        tuple(range(1, size)) if mirrored else tuple(-x for x in range(size - 1, 0, -1))
    )
    if size == 2 and mirrored:
        with pytest.raises(framed.ControlError, match="normal X2 convention"):
            framed.initialize(n, incoming, size=size, mirrored=True)
        return
    state = framed.initialize(
        n, incoming, size=size, start=1, mirrored=mirrored, mirrored_policy="bar_mirror"
    )
    work = algorithm._pass_work(_root_common(n), "zip_up", state.branches[0])
    collected._set_coefficient(work, [[-2, -3], [1, 2]])
    original = deepcopy(work)
    session = memoized.MemoizedSession()
    expected = collected._transition(work, Limits(), True, session._catalog)
    selected = [(collected.continuation_key(work), work)]
    cold, _ = session._execute_batch(selected, Limits(), True, 1, 64, None)
    warm, _ = session._execute_batch(selected, Limits(), True, 1, 64, None)
    assert cold == warm == [expected]
    assert work == original
    assert session.cache_info()["hits"] == 1


@pytest.mark.parametrize("size", [2, 3, 5, 8, 9])
def test_growth_to_all_ambient_slider_sizes_is_uncapped(size):
    n = size + 1
    work = dict(
        _root_common(n),
        phase="layer",
        payload_length=0,
        completed_sides=0,
        word=list(range(size, 0, -1)) + [["slider", size, 2, False, False, False]],
    )
    session = memoized.MemoizedSession()
    expected = collected._transition(work, Limits(), True, session._catalog)
    selected = [(collected.continuation_key(work), work)]
    cold, _ = session._execute_batch(selected, Limits(), True, 1, 64, None)
    warm, _ = session._execute_batch(selected, Limits(), True, 1, 64, None)
    assert cold == warm == [expected]
    assert (
        max(x["branch"]["control"]["i"] for x in cold[0]["children"] if "branch" in x)
        == n
    )
