"""Exact rank-two schedule, checkpoint, branch, and coefficient regressions."""
from copy import deepcopy
from itertools import product
import json

import pytest

from zipper_algebra.zipping_algorithm.two_strand import (
    TwoStrandSession, check_two_strand, reduce_two_strand, resume_two_strand,
)


def output_map(run):
    return {(term["kind"], term["sigma_power"]):
            {c["q_power"]: c["integer"] for c in term["coefficient"]}
            for term in run["output"]}


def test_two_strand_definition_and_inverse_foe():
    positive = reduce_two_strand([1])
    assert output_map(positive) == {
        ("slider", 0): {0: -1},
        ("singletons", 0): {0: 1, 1: -1},
        ("singletons", -1): {1: 1},
    }
    negative = reduce_two_strand([-1, -1], trace=True)
    assert output_map(negative) == {
        ("slider", -1): {-1: 1},
        ("singletons", -1): {-1: -1, 0: 1},
        ("singletons", 0): {-1: 1},
    }
    assert negative["transitions"][1]["rule_id"] == "single.foe.right.negative.I"
    assert check_two_strand(positive)["pass"]
    assert check_two_strand(negative)["pass"]


def test_retains_processed_tail_and_branch_isolation():
    run = reduce_two_strand([-1, -1], trace=True)
    step = run["transitions"][1]
    assert step["before"]["processed"] == [-1]
    assert step["evidence"]["extracted_tau"] == -1
    assert all(s["processed"] == [-1] for s in step["pending_after"])
    checkpoint = deepcopy(run)
    step["pending_after"][0]["processed"].append(1)
    assert step["pending_after"][1]["processed"] == [-1]
    assert checkpoint["transitions"][1]["pending_after"][0]["processed"] == [-1]
    with pytest.raises(ValueError, match="replay mismatch"):
        check_two_strand(run)


def test_inverse_cancellation_and_exact_collection():
    assert output_map(reduce_two_strand([1, -1])) == {("singletons", 0): {0: 1}}
    assert output_map(reduce_two_strand([-1, 1])) == {("singletons", 0): {0: 1}}
    run = reduce_two_strand([1, -1, 1, 1])
    assert output_map(run) == output_map(reduce_two_strand([1, 1]))
    assert run["metrics"]["uncollected_output_terms"] > run["metrics"]["collected_output_terms"]


@pytest.mark.parametrize("word", [
    word for length in range(6) for word in product((-1, 1), repeat=length)
])
def test_exhaustive_small_words_independently_expand_in_group_algebra(word):
    length = len(word)
    run = reduce_two_strand(word)
    assert run["status"] == "complete"
    assert not run["pending"]
    assert check_two_strand(json.loads(json.dumps(run)))["pass"]
    assert all(t["sigma_power"] in (-1, 0) for t in run["output"] if t["kind"] == "singletons")
    assert run["metrics"]["steps"] <= 3*2**length-2
    assert run["metrics"]["max_stored_braid_length"] <= length


def test_step_and_term_limits_resume_exact_queue():
    word = [1, 1, -1, 1]
    checkpoint = reduce_two_strand(word, max_steps=0)
    assert checkpoint["status"] == "incomplete"
    assert checkpoint["diagnostic"] == "step_limit"
    assert check_two_strand(checkpoint)["pass"]
    limited = reduce_two_strand(word, max_terms=1)
    assert limited["diagnostic"] == "term_limit"
    assert check_two_strand(limited)["pass"]
    complete = reduce_two_strand(word)
    assert resume_two_strand(checkpoint) == complete
    assert resume_two_strand(limited) == complete
    assert resume_two_strand(complete) == complete
    partial = resume_two_strand(checkpoint, max_steps=2)
    assert partial["metrics"]["steps"] == 2
    assert partial["transitions"] == []
    assert check_two_strand(partial)["pass"]
    with pytest.raises(ValueError, match="precedes"):
        resume_two_strand(partial, max_steps=1)


def test_rank_one_and_empty_rank_two():
    for n in (1, 2):
        run = reduce_two_strand([], n=n)
        assert output_map(run) == {("singletons", 0): {0: 1}}
        assert check_two_strand(run)["pass"]


@pytest.mark.parametrize("word,n", [([1], 1), ([], 0), ([], 3), ([], True), ([0], 2), ([2], 2), ([True], 2), ("1", 2)])
def test_reject_invalid_input(word, n):
    with pytest.raises(ValueError):
        reduce_two_strand(word, n=n)


@pytest.mark.parametrize("kwargs", [{"max_steps": -1}, {"max_steps": True}, {"max_terms": 0}, {"max_terms": False}])
def test_reject_invalid_limits(kwargs):
    with pytest.raises(ValueError):
        reduce_two_strand([], **kwargs)
    with pytest.raises(ValueError):
        resume_two_strand(reduce_two_strand([]), **kwargs)


def test_certificate_strict_replay_and_malformed_records():
    with pytest.raises(ValueError, match="schema"):
        check_two_strand([])
    with pytest.raises(ValueError, match="schema"):
        check_two_strand({"schema": "other"})
    for field in ("input_word", "limits", "n"):
        run = reduce_two_strand([])
        run.pop(field)
        with pytest.raises(ValueError, match="Malformed"):
            check_two_strand(run)
    run = reduce_two_strand([])
    run["limits"] = None
    with pytest.raises(ValueError, match="Malformed"):
        check_two_strand(run)
    for field, replacement in (("catalog_sha256", "0"*64), ("status", "incomplete"), ("source_sha256", {})):
        run = reduce_two_strand([1])
        run[field] = replacement
        with pytest.raises(ValueError, match="replay mismatch"):
            check_two_strand(run)


def test_independent_checker_rejects_false_output_identity():
    run = reduce_two_strand([1])
    run["output"][0]["coefficient"][0]["integer"] = 2
    with pytest.raises(ValueError, match="Laurent-group identity failed"):
        check_two_strand(run)


def test_certificate_rejects_boolean_integer_alias():
    run = reduce_two_strand([], trace=True)
    run["transitions"][0]["evidence"]["largest_bottom_label"] = 2.0
    with pytest.raises(ValueError, match="replay mismatch"):
        check_two_strand(run)


def test_checkpoint_queue_does_not_alias_saved_transition():
    run = reduce_two_strand([1], max_steps=1, trace=True)
    run["pending"][0]["coefficient"][0]["integer"] = 999
    assert run["transitions"][0]["pending_after"][0]["coefficient"][0]["integer"] != 999


def test_resume_rejects_term_bound_below_saved_peak():
    complete = reduce_two_strand([1, 1])
    with pytest.raises(ValueError, match="recorded peak"):
        resume_two_strand(complete, max_terms=1)
    changed = resume_two_strand(complete, max_terms=complete["metrics"]["max_total_terms"])
    assert check_two_strand(changed)["pass"]


def test_optional_trace_has_same_calculation_and_streaming_metrics():
    word = [1, 1, -1, -1, 1, -1, 1]
    compact = reduce_two_strand(word)
    detailed = reduce_two_strand(word, trace=True)
    assert compact["input_word"] == detailed["input_word"] == word
    assert compact["transitions"] == []
    assert compact["trace_enabled"] is False
    assert detailed["trace_enabled"] is True
    assert detailed["transitions"]
    assert compact["output"] == detailed["output"]
    assert compact["pending"] == detailed["pending"]
    assert compact["metrics"] == detailed["metrics"]
    assert check_two_strand(compact)["pass"]
    assert check_two_strand(detailed)["pass"]
    def reject_history(value):
        if isinstance(value, dict):
            assert not {"before", "pending_after", "output_added", "evidence"} & set(value)
            for child in value.values():
                reject_history(child)
        elif isinstance(value, list):
            for child in value:
                reject_history(child)
    reject_history(compact)


def test_compact_and_detailed_resume_preserve_original_modes():
    word = [-1, -1, 1, 1]
    for trace in (False, True):
        checkpoint = reduce_two_strand(word, max_steps=2, trace=trace)
        complete = resume_two_strand(checkpoint, trace=trace)
        assert complete == reduce_two_strand(word, trace=trace)
        assert check_two_strand(complete)["pass"]
        with pytest.raises(ValueError, match="trace mode"):
            resume_two_strand(checkpoint, trace=not trace)
        with pytest.raises(ValueError, match="trace mode"):
            resume_two_strand(checkpoint, trace=1)


def test_trace_is_strict_boolean_and_certificate_binds_mode():
    for value in (None, 1, 0, "true"):
        with pytest.raises(ValueError, match="trace must be a boolean"):
            reduce_two_strand([], trace=value)
    run = reduce_two_strand([])
    run["trace_enabled"] = 0
    with pytest.raises(ValueError, match="trace must be a boolean"):
        check_two_strand(run)
    run = reduce_two_strand([1])
    run["transitions"] = [{"kind": "fake"}]
    with pytest.raises(ValueError, match="replay mismatch"):
        check_two_strand(run)


def test_literal_inverse_input_is_preserved_with_optional_history():
    word = [1, -1] * 40
    compact = reduce_two_strand(word)
    detailed = reduce_two_strand(word, trace=True)
    assert compact["input_word"] == detailed["input_word"] == word
    assert len(compact["input_word"]) == 80
    assert compact["metrics"] == detailed["metrics"]
    assert compact["metrics"]["steps"] > 0
    assert compact["transitions"] == []
    assert len(detailed["transitions"]) == detailed["metrics"]["steps"]


def test_session_reuses_catalog_and_matches_fresh_operations(monkeypatch):
    from zipper_algebra.zipping_algorithm.catalog import Catalog
    bundled = Catalog.bundled
    loads = []
    def counted_load():
        loads.append(True)
        return bundled()
    monkeypatch.setattr(Catalog, "bundled", counted_load)
    session = TwoStrandSession()
    assert len(loads) == 1
    for trace in (False, True):
        first = session.reduce([1, -1, 1], trace=trace)
        checkpoint = session.reduce([-1, -1], max_steps=1, trace=trace)
        assert session.check(first)["pass"]
        resumed = session.resume(checkpoint)
        assert session.check(resumed)["pass"]
        assert len(loads) == (1 if not trace else 5)
        assert first == reduce_two_strand([1, -1, 1], trace=trace)
        assert resumed == reduce_two_strand([-1, -1], trace=trace)
        assert session.check(first) == check_two_strand(first)
        assert resumed == resume_two_strand(checkpoint)
    assert len(loads) == 9


def test_batch_runs_do_not_share_algorithm_state_or_modify_inputs():
    session = TwoStrandSession()
    word = [1, 1]
    reference = session.reduce(word, trace=True)
    saved = deepcopy(reference)
    word.append(-1)
    reference["transitions"][0]["before"]["coefficient"][0]["integer"] = 999
    reference["output"].clear()
    reference["metrics"]["steps"] = -1
    again = session.reduce([1, 1], trace=True)
    assert again == saved
    assert session.check(again)["pass"]
    checkpoint = session.reduce(word, max_steps=1)
    checkpoint_saved = deepcopy(checkpoint)
    session.resume(checkpoint)
    assert checkpoint == checkpoint_saved
    failed = deepcopy(again)
    failed["output"].clear()
    with pytest.raises(ValueError):
        session.check(failed)
    assert session.reduce([1, 1], trace=True) == saved


def test_batch_and_fresh_match_all_small_literal_words():
    session = TwoStrandSession()
    for length in range(4):
        for word in product((-1, 1), repeat=length):
            for trace in (False, True):
                actual = session.reduce(word, trace=trace)
                assert actual == reduce_two_strand(word, trace=trace)
                assert session.check(actual)["pass"]
