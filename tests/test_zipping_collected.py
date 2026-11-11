"""Exact collection, weighted local edges, and deterministic process batches."""

from copy import deepcopy
from itertools import product
import json
import os

import pytest

from zipper_algebra.zipping_algorithm import algorithm as reference
from zipper_algebra.zipping_algorithm import collected
from zipper_algebra.zipping_algorithm.interpreter import Limits


def _outputs(result):
    return {
        collected.continuation_key(item, pending=False): item["coefficient"]
        for item in collected.collect_terms(result["output"], pending=False)
    }


def _scale(item, coefficient):
    after = deepcopy(item)
    after["coefficient"] = deepcopy(coefficient)
    if "branch" in after:
        after["branch"]["coefficient"] = deepcopy(coefficient)
    return after


def _zip_state():
    result = reference.reduce_word((1,), n=2, limits=Limits(max_steps=1))
    return next(item for item in result["pending"] if item["phase"] == "zip_up")


def _crash_worker(_job):
    os._exit(23)


def test_every_three_strand_literal_through_length_two_matches_reference_sum():
    session = collected.CollectedSession()
    baseline = reference.AlgorithmSession()
    for length in range(3):
        for word in product((-2, -1, 1, 2), repeat=length):
            actual = session.reduce(word, n=3)
            expected = baseline.reduce(word, n=3)
            assert actual["status"] == expected["status"] == "complete"
            assert _outputs(actual) == _outputs(expected), word
            assert actual["input_word"] == list(word)
    assert session.reduce((), n=1)["status"] == "complete"


def test_collection_synchronizes_scalar_copies_depth_and_removes_zeroes():
    state = _zip_state()
    original = deepcopy(state)
    a = _scale(state, [[-2, 2], [3, 1]])
    b = _scale(state, [[-2, -2], [1, 4]])
    a["depth"], b["depth"] = 3, 8
    terms = collected.collect_terms([a, b])
    assert len(terms) == 1
    assert terms[0]["coefficient"] == [[1, 4], [3, 1]]
    assert terms[0]["branch"]["coefficient"] == terms[0]["coefficient"]
    assert terms[0]["depth"] == 8
    opposite = _scale(terms[0], [[1, -4], [3, -1]])
    assert collected.collect_terms([terms[0], opposite]) == []
    assert collected.collect_terms([terms[0], opposite, state]) == [state]
    assert state == original


@pytest.mark.parametrize(
    "field,value",
    [
        ("world_phi", True),
        ("bottom_labels", [2, 1]),
        ("top_labels", [2, 1]),
        ("boundary_left", [1]),
        ("boundary_right", [-1]),
        ("boundary_left_factors", [{"growth_size": 2, "word": [1]}]),
        ("boundary_right_factors", [{"growth_size": 2, "word": [-1]}]),
        ("contexts", [{"rank": 3, "size": 1}]),
    ],
)
def test_full_continuation_key_retains_every_boundary_and_label(field, value):
    state = _zip_state()
    different = deepcopy(state)
    different[field] = value
    assert len(collected.collect_terms([state, different])) == 2


def test_branch_control_and_completed_component_records_remain_distinct():
    state = _zip_state()
    different = deepcopy(state)
    different["branch"]["control"]["k"] += 1
    assert collected.continuation_key(state) != collected.continuation_key(different)
    term = collected.reduce_word((1,), n=2)["output"][0]
    other = deepcopy(term)
    other["components"][0]["record"]["orientation"] = "mirrored"
    assert collected.continuation_key(
        term, pending=False
    ) != collected.continuation_key(other, pending=False)
    assert len(collected.collect_terms([term, other], pending=False)) == 2
    doubled = collected.collect_terms([term, term], pending=False)
    assert doubled[0]["coefficient"] == reference._multiply(
        term["coefficient"], [[0, 2]]
    )


def test_real_branching_collects_and_cancels_during_the_calculation():
    word = (1, -2, 1, -2)
    actual = collected.reduce_word(word, n=3, trace=True)
    baseline = reference.reduce_word(word, n=3)
    assert actual["status"] == "complete"
    assert actual["metrics"]["merged_terms"] > 0
    assert actual["metrics"]["zero_terms"] > 0
    assert actual["metrics"]["steps"] < baseline["metrics"]["steps"]
    assert _outputs(actual) == _outputs(baseline)
    assert any(
        value is None
        for batch in actual["transitions"]
        for value in batch["pending_updates"].values()
    )
    compact = collected.reduce_word(word, n=3)
    assert compact["output"] == actual["output"]
    assert compact["metrics"] == actual["metrics"]
    assert not compact["transitions"]
    assert collected.check_certificate(json.loads(json.dumps(actual)))


def test_collected_cap_can_complete_when_uncollected_fanout_exceeds_it():
    word = (1, -2, 1, -2)
    limits = Limits(max_terms=26)
    actual = collected.reduce_word(word, n=3, limits=limits)
    baseline = reference.reduce_word(word, n=3, limits=limits)
    assert actual["status"] == "complete"
    assert actual["metrics"]["max_terms"] == 26
    assert baseline["status"] == "incomplete"
    assert "term limit" in baseline["diagnostic"]["message"]


def test_incremental_batches_merge_into_unselected_states_and_keep_finished_terms():
    result = collected.reduce_word((1, 2, 1, 2), n=3, batch_size=4, trace=True)
    initial = result["transitions"][0]["predecessors"][0]
    frontier = {collected.continuation_key(initial): initial}
    merged_into_unselected = False
    for batch in result["transitions"]:
        selected = {collected.continuation_key(x) for x in batch["predecessors"]}
        incoming = {
            collected.continuation_key(x)
            for edge in batch["rewrites"]
            for x in edge["children"]
        }
        merged_into_unselected |= bool((frontier.keys() - selected) & incoming)
        for key in selected:
            del frontier[key]
        for key, value in batch["pending_updates"].items():
            if value is None:
                frontier.pop(key, None)
            else:
                frontier[key] = value
    assert merged_into_unselected and not frontier
    stopped = collected.reduce_word(
        (1, -2, 1, -2), n=3, batch_size=4, limits=Limits(max_steps=120)
    )
    assert stopped["status"] == "incomplete" and stopped["output"]
    assert collected.resume(stopped, limits=Limits()) == collected.reduce_word(
        (1, -2, 1, -2), n=3, batch_size=4
    )


def test_two_processes_produce_identical_certificates_and_shut_down_on_refusal():
    word = (1, 2, 1)
    serial = collected.reduce_word(word, n=3, trace=True, batch_size=4)
    parallel = collected.reduce_word(word, n=3, trace=True, batch_size=4, workers=2)
    assert parallel == serial
    assert collected.check_certificate(parallel)
    bounded = Limits(max_word_size=5)
    serial_stop = collected.reduce_word((1, 2), n=3, limits=bounded, batch_size=4)
    parallel_stop = collected.reduce_word(
        (1, 2), n=3, limits=bounded, batch_size=4, workers=2
    )
    assert serial_stop == parallel_stop
    assert parallel_stop["status"] == "incomplete"
    assert parallel_stop["diagnostic"]["batch_rolled_back"]


def test_real_worker_death_reports_committed_prefix_without_a_replay_certificate(
    monkeypatch,
):
    # The initial singleton transition commits in the coordinator. Its three
    # genuine successors then reach spawned workers, one of which exits.
    monkeypatch.setattr(collected, "_worker", _crash_worker)
    with pytest.raises(
        collected.CollectedError,
        match="1 committed batches and 1 committed term transitions",
    ):
        collected.reduce_word((1, 2), n=3, workers=2, batch_size=4)


def test_unit_worker_edges_are_linear_for_every_real_phase():
    run = reference.reduce_word((1, 2), n=3, trace=True)
    states = {}
    for event in run["transitions"]:
        state = event["predecessor"]
        states.setdefault(state["phase"], state)
    assert {"single", "zip_up", "zip_down", "layer", "extract"} <= states.keys()
    collected._initialize_worker()
    weight = [[-2, 3], [1, -2]]
    for state in states.values():
        unit = _scale(state, [[0, 1]])
        base = collected._worker((unit, Limits(), False))
        weighted = collected._worker((_scale(state, weight), Limits(), False))
        for kind in ("children", "output"):
            assert len(base[kind]) == len(weighted[kind])
            for a, b in zip(base[kind], weighted[kind], strict=True):
                assert b["coefficient"] == reference._multiply(weight, a["coefficient"])
                assert collected.continuation_key(
                    a, pending=kind == "children"
                ) == collected.continuation_key(b, pending=kind == "children")


@pytest.mark.parametrize(
    "word,n,limits,message",
    [
        ((1,), 2, Limits(max_terms=1), "collected term"),
        ((1,), 2, Limits(max_coefficient_support=1), "coefficient support"),
        ((1, 1), 2, Limits(max_integer_bits=1), "integer bit"),
        ((1,), 2, Limits(max_word_size=1), "word size"),
        ((1, 2), 3, Limits(max_word_size=5), "word size"),
        ((1, 2), 3, Limits(expansion_budget=1), "budget"),
        ((1, 2), 3, Limits(max_steps=2), "atomic batch"),
    ],
)
def test_real_resource_refusals_keep_exact_predecessor_batch(word, n, limits, message):
    result = collected.reduce_word(word, n=n, limits=limits)
    assert result["status"] == "incomplete"
    assert message in result["diagnostic"]["message"]
    assert result["diagnostic"]["batch_rolled_back"]
    assert result["pending"]
    assert collected.check_certificate(result)


def test_coefficient_and_term_caps_are_after_collection_and_rollback_is_atomic():
    # A local 3-branch rewrite is allowed to finish before its frontier cap.
    term_stop = collected.reduce_word((1,), n=2, limits=Limits(max_terms=1))
    assert term_stop["metrics"]["steps"] == 0
    assert term_stop["metrics"]["attempted_steps"] == 1
    assert term_stop["diagnostic"]["phase"] == "batch"
    assert "batch_representative" in term_stop["diagnostic"]
    assert "state" not in term_stop["diagnostic"]
    # This later refusal follows three simultaneous transitions. None commits.
    bit_stop = collected.reduce_word((1, 1), n=2, limits=Limits(max_integer_bits=1))
    step_stop = collected.reduce_word((1, 1), n=2, limits=Limits(max_steps=1))
    assert bit_stop["metrics"]["steps"] == 1
    assert bit_stop["metrics"]["attempted_steps"] == 4
    assert bit_stop["pending"] == step_stop["pending"]
    assert bit_stop["output"] == step_stop["output"]


def test_resume_replays_atomic_prefix_and_workers_are_execution_only():
    stopped = collected.reduce_word(
        (1, 2), n=3, limits=Limits(max_steps=15), trace=True
    )
    assert stopped["metrics"]["steps"] == 10
    continued = collected.resume(stopped, limits=Limits(), workers=2)
    assert continued == collected.reduce_word((1, 2), n=3, trace=True)
    assert collected.resume(stopped) == stopped
    with pytest.raises(collected.CollectedError, match="must not decrease"):
        collected.resume(stopped, limits=Limits(max_steps=14))


@pytest.mark.parametrize(
    "word,options",
    [
        ((1,), {"n": True}),
        ((1,), {"n": 0}),
        ((1,), {"n": 2, "trace": 1}),
        ((1,), {"n": 2, "workers": True}),
        ((1,), {"n": 2, "workers": 0}),
        ((1,), {"n": 2, "batch_size": 0}),
        ((1,), {"n": 2, "batch_size": False}),
        ((True,), {"n": 2}),
        ((0,), {"n": 2}),
        ((2,), {"n": 2}),
        (iter((1,)), {"n": 2}),
    ],
)
def test_public_input_validation(word, options):
    with pytest.raises(collected.CollectedError):
        collected.reduce_word(word, **options)


def test_replay_rejects_missing_fields_and_altered_sums_and_sessions_are_isolated():
    session = collected.CollectedSession()
    first = session.reduce((1, 2), n=3)
    other = session.reduce((-1,), n=2)
    assert first["input_word"] == [1, 2] and other["input_word"] == [-1]
    with pytest.raises(collected.CollectedError, match="malformed"):
        session.check({})
    altered = deepcopy(first)
    altered["output"][0]["coefficient"][0][1] += 1
    with pytest.raises(collected.CollectedError, match="differs"):
        session.check(altered)


@pytest.mark.parametrize("invalid", [0, -1, True, "ten"])
def test_replay_wraps_invalid_certificate_limits(invalid):
    result = collected.reduce_word((), n=1)
    result["limits"]["max_steps"] = invalid
    with pytest.raises(collected.CollectedError, match="malformed"):
        collected.check_certificate(result)
