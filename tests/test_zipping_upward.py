"""Exact contexts, exceptional branches, frame refusals and replay checks."""
from copy import deepcopy
import json

import pytest

from zipper_algebra.zipping_algorithm.catalog import Catalog, CatalogError
from zipper_algebra.zipping_algorithm.interpreter import Branch, Limits, ReplayError, RunState
from zipper_algebra.zipping_algorithm.upward import (
    ControlError, check_certificate, initialize, resume, run, step,
)


@pytest.fixture(scope="module")
def catalog():
    return Catalog.bundled()


def slider(i=2, j=1, mirrored=False):
    return ("slider", i, j, mirrored, False, False)


def test_far_pass_retains_full_processed_tail_and_inverse_pair(catalog):
    initial = initialize(5, (4, -4), processed=(1,), tail=(2,))
    certificate = run(initial, catalog=catalog, trace=True)
    assert certificate["status"] == "upper_pass_complete"
    final = RunState.from_json(certificate["final"])
    assert final.branches[0].word == (slider(), 4, -4, 1, 2)
    assert final.branches[0].control["processed_length"] == 3
    assert final.branches[0].control["tail_length"] == 1
    assert initial.branches[0].word == (4, -4, slider(), 1, 2)
    assert [x["catalog_transition"]["request"]["rule_id"] for x in certificate["transitions"]] == [
        "zip.far.negative.I", "zip.far.positive.I"]
    assert check_certificate(json.loads(json.dumps(certificate)), catalog=catalog)


def test_exceptional_nearest_branch_reprocesses_its_prefix(catalog):
    initial = initialize(3, (-1,), friends=1, processed=(2,), tail=(-2,))
    following, event = step(initial, catalog=catalog)
    assert len(following.branches) == 4
    assert [b.coefficient for b in following.branches] == [
        ((1, -1),), ((0, 1),), ((1, 1),), ((0, 1), (1, -1))]
    fourth = following.branches[3]
    assert fourth.word == (-1, slider(), 2, -2)
    assert fourth.control["left_length"] == 1
    assert fourth.control["k"] == 0
    assert fourth.control["processed_length"] == 1
    assert fourth.control["tail_length"] == 1
    # Sibling control dictionaries are independent values.
    changed = fourth.control
    changed["processed_length"] = 999
    assert following.branches[0].control["processed_length"] == 5
    certificate = run(initial, catalog=catalog, trace=True)
    assert certificate["status"] == "upper_pass_complete"
    assert len(certificate["transitions"]) == 2
    assert certificate["transitions"][1]["catalog_transition"]["request"]["rule_id"] == "zip.internal.i2.negative.I"
    assert event["normalizations"] == []


def test_farthest_positive_keeps_four_branches_and_records_mirror_identity(catalog):
    after, event = step(initialize(3, (2,), friends=1), catalog=catalog)
    assert [b.word for b in after.branches] == [
        (slider(j=2), 2, -1, -2), (slider(j=2), 2, 1, 2),
        (slider(),), (2, slider())]
    assert [b.coefficient for b in after.branches] == [
        ((1, -1),), ((0, 1),), ((1, 1),), ((0, 1), (1, -1))]
    assert event["normalizations"] == [
        {"branch": 0, "identity": "mirrored_X2_equals_X2"},
        {"branch": 1, "identity": "mirrored_X2_equals_X2"}]
    assert after.branches[-1].control["k"] == 1


@pytest.mark.parametrize("word,expected", [
    ((2, 1), (-1, slider(3, mirrored=True), 2, 1, 2)),
    ((-2, -1), (1, slider(3, mirrored=True), -2, -1, -2)),
])
def test_higher_flypes_retain_both_half_twists(catalog, word, expected):
    certificate = run(initialize(3, word, size=3), catalog=catalog, trace=True)
    assert certificate["status"] == "unsupported_frame"
    assert RunState.from_json(certificate["final"]).branches[0].word == expected
    assert certificate["transitions"][0]["normalizations"] == []
    assert check_certificate(certificate, catalog=catalog)
    continued = resume(certificate, catalog=catalog)
    assert continued["status"] == "unsupported_frame"
    assert continued["transitions"] == []
    assert continued["final"] == certificate["final"]


def test_higher_flype_keeps_unprocessed_exterior(catalog):
    certificate = run(initialize(4, (3, 2, 1), size=3, processed=(-3,), tail=(3,)), catalog=catalog, trace=True)
    assert certificate["status"] == "unsupported_frame"
    final = RunState.from_json(certificate["final"])
    assert final.branches[0].word == (3, -1, slider(3, mirrored=True), 2, 1, 2, -3, 3)


@pytest.mark.parametrize("n,left,i,j,k,rule", [
    (5, (4,), 2, 1, 0, "zip.far.positive.I"),
    (5, (-4,), 2, 1, 0, "zip.far.negative.I"),
    (5, (2,), 2, 1, 3, "zip.close.positive.I"),
    (5, (-2,), 2, 1, 3, "zip.close.negative.I"),
    (2, (1,), 2, 1, 0, "zip.internal.i2.positive.I"),
    (2, (-1,), 2, 1, 0, "zip.internal.i2.negative.I"),
    (3, (1,), 2, 1, 1, "zip.nearest.positive.I"),
    (3, (-1,), 2, 1, 1, "zip.nearest.negative.i2.I"),
    (3, (1,), 2, 2, 0, "zip.left.positive.i2.I"),
    (3, (-1,), 2, 2, 0, "zip.left.negative.I"),
    (3, (-2,), 2, 1, 1, "zip.farthest.negative.I"),
    (3, (2,), 2, 1, 1, "zip.farthest.positive.i2.I"),
    (3, (2,), 2, 1, 0, "zip.outer.positive.I"),
    (3, (-2,), 2, 1, 0, "zip.outer.negative.i2.I"),
    (4, (-1,), 3, 1, 1, "zip.nearest.negative.i_ge_3.I"),
    (4, (1,), 3, 2, 0, "zip.left.positive.i_ge_3.I"),
    (4, (3, 2), 3, 1, 1, "zip.farthest.positive.i_ge_3.I"),
    (5, (-3, -4), 3, 1, 1, "zip.outer.negative.i_ge_3.I"),
])
def test_dispatch_families_and_strict_i2_progress(catalog, n, left, i, j, k, rule):
    initial = initialize(n, left, size=i, start=j, friends=k)
    after, event = step(initial, catalog=catalog)
    assert event["catalog_transition"]["request"]["rule_id"] == rule
    if i == 2:
        potential = (n + 1) * len(left) + n * k
        assert all((n + 1) * b.control["left_length"] + n * b.control["k"] < potential
                   for b in after.branches)


def test_exceptional_long_friend_prefix_is_not_misclassified_as_friends(catalog):
    after, _ = step(initialize(6, (-1,), friends=4), catalog=catalog)
    fourth = after.branches[3]
    assert fourth.word == (-1, 5, 4, 3, slider())
    assert fourth.control["left_length"] == 4
    assert fourth.control["k"] == 0
    certificate = run(initialize(6, (-1,), friends=4), catalog=catalog, trace=True)
    assert certificate["status"] == "upper_pass_complete"


def test_missing_full_tape_retains_predecessor(catalog):
    initial = initialize(3, (1,), size=3)
    certificate = run(initial, catalog=catalog, trace=True)
    assert certificate["status"] == "missing_rule"
    assert certificate["final"] == initial.to_json()
    assert certificate["transitions"] == []
    assert certificate["diagnostic"]["message"] == "missing full-tape rule"
    assert check_certificate(certificate, catalog=catalog)


def test_competing_rules_fail_explicitly(catalog):
    records = list(catalog.rules)
    duplicate = catalog.get("zip.far.positive.I")
    duplicate["id"] = "zip.duplicate.positive.I"
    records.append(duplicate)
    ambiguous = Catalog({"schema_version": "1.0", "rules": records})
    certificate = run(initialize(4, (3,)), catalog=ambiguous, trace=True)
    assert certificate["status"] == "missing_rule"
    assert "competing upward rules" in certificate["diagnostic"]["message"]
    assert "zip.duplicate.positive.I" in certificate["diagnostic"]["message"]


def test_limits_preserve_atomic_checkpoint_and_resume(catalog):
    initial = initialize(3, (-1,), friends=1)
    too_many = run(initial, catalog=catalog, trace=True, limits=Limits(max_terms=3))
    assert too_many["status"] == "incomplete"
    assert too_many["final"] == initial.to_json()
    assert too_many["transitions"] == []
    limited = run(initial, catalog=catalog, trace=True, limits=Limits(max_steps=1))
    assert limited["status"] == "incomplete"
    assert len(limited["transitions"]) == 1
    continued = resume(json.loads(json.dumps(limited)), catalog=catalog)
    assert continued["status"] == "upper_pass_complete"
    assert check_certificate(continued, catalog=catalog)


def test_initial_limit_even_when_upward_word_empty(catalog):
    state = initialize(2, (), processed=(1, 1, 1))
    result = run(state, catalog=catalog, trace=True, limits=Limits(max_word_size=3))
    assert result["status"] == "incomplete"
    assert result["final"] == state.to_json()
    assert check_certificate(result, catalog=catalog)
    assert run(initialize(2, ()), catalog=catalog, trace=True)["status"] == "upper_pass_complete"
    assert run(RunState(2, ()), catalog=catalog, trace=True)["status"] == "upper_pass_complete"


def test_expansion_limit_is_an_incomplete_resumable_checkpoint(catalog):
    state = initialize(2, (1,))
    limited = run(state, catalog=catalog, trace=True, limits=Limits(expansion_budget=1))
    assert limited["status"] == "incomplete"
    assert limited["final"] == state.to_json()
    assert limited["diagnostic"]["message"] == "Expansion word budget exceeded"
    assert check_certificate(limited, catalog=catalog)
    assert resume(limited, catalog=catalog)["status"] == "upper_pass_complete"


def test_certificate_tamper_rejected(catalog):
    certificate = run(initialize(3, (2,), friends=1), catalog=catalog, trace=True)
    mutations = [
        lambda c: c.update(extra=True),
        lambda c: c.update(source_sha256="0" * 64),
        lambda c: c.update(catalog_sha256="0" * 64),
        lambda c: c.update(status="complete"),
        lambda c: c["transitions"][0].update(normalizations=[]),
        lambda c: c["final"]["branches"][0]["control"].update(j=1),
        lambda c: c["final"]["branches"][0]["word"].append(2),
        lambda c: c["final"]["branches"][0]["coefficient"][0].__setitem__(1, 10),
    ]
    for mutate in mutations:
        altered = deepcopy(certificate)
        mutate(altered)
        with pytest.raises(ControlError, match="differs"):
            check_certificate(altered, catalog=catalog)


@pytest.mark.parametrize("kwargs", [
    {"n": True}, {"size": 1}, {"size": True}, {"start": 0}, {"friends": -1},
    {"friends": True}, {"size": 4}, {"start": 3}, {"friends": 2},
    {"left": (True,)}, {"processed": (slider(),)}, {"tail": (False,)},
    {"left": (0,)},
])
def test_invalid_initialization(kwargs):
    args = dict(n=3, left=())
    args.update(kwargs)
    with pytest.raises(ReplayError):
        initialize(**args)


@pytest.mark.parametrize("branch", [-1, 1, True])
def test_invalid_branch_selection(catalog, branch):
    with pytest.raises(ControlError, match="does not exist"):
        step(initialize(2, (1,)), branch, catalog=catalog)


def test_complete_pass_cannot_be_stepped(catalog):
    with pytest.raises(ControlError, match="already complete"):
        step(initialize(2, ()), catalog=catalog)


@pytest.mark.parametrize("mutation", [
    lambda c: c.update(extra=True),
    lambda c: c.update(phase="other"),
    lambda c: c.update(i=True),
    lambda c: c.update(i=1),
    lambda c: c.update(j=0),
    lambda c: c.update(k=-1),
    lambda c: c.update(j=3),
    lambda c: c.update(left_length=-1),
    lambda c: c.update(processed_length=1),
    lambda c: c.update(phase="upper_pass_complete"),
    lambda c: c.update(k=1),
])
def test_malformed_saved_boundaries_rejected(catalog, mutation):
    initial = initialize(3, (1,))
    c = initial.branches[0].control
    mutation(c)
    state = RunState(3, (Branch(initial.branches[0].word, control_json=json.dumps(c)),))
    with pytest.raises(ControlError):
        run(state, catalog=catalog, trace=True)


def test_slider_in_exterior_rejected(catalog):
    initial = initialize(3, (1,))
    state = RunState(3, (Branch((slider(), slider()), control_json=initial.branches[0].control_json),))
    with pytest.raises(ControlError, match="pure braid"):
        run(state, catalog=catalog, trace=True)


def test_retained_frame_metadata_is_strict(catalog):
    state = RunState(3, (Branch((slider(3, mirrored=True),), control_json='{"phase":"unsupported_frame"}'),))
    with pytest.raises(ControlError, match="frame record"):
        run(state, catalog=catalog, trace=True)


def test_optional_trace_preserves_arithmetic_and_streaming_metrics(catalog):
    state = initialize(3, (-1,), friends=1)
    compact = run(state, catalog=catalog)
    detailed = run(state, catalog=catalog, trace=True)
    assert compact["trace_enabled"] is False
    assert detailed["trace_enabled"] is True
    assert compact["final"] == detailed["final"]
    assert compact["summary"] == detailed["summary"]
    assert compact["summary"]["steps"] == 2
    assert compact["summary"]["rule_counts"] == {
        "zip.nearest.negative.i2.I": 1, "zip.internal.i2.negative.I": 1}
    assert compact["summary"]["peak_uncollected_terms"] == 4
    assert compact["transitions"] == []
    assert len(compact["metrics"]) == 2
    assert len(detailed["metrics"]) == 3
    assert compact["metrics"] == [detailed["metrics"][0], detailed["metrics"][-1]]
    assert check_certificate(compact, catalog=catalog)
    assert check_certificate(detailed, catalog=catalog)
    resumed = resume(compact, catalog=catalog, trace=True)
    assert resumed["trace_enabled"] is True
    assert resumed["final"] == compact["final"]


def test_compact_mode_serializes_only_input_and_final_state(catalog, monkeypatch):
    state = initialize(3, (-1,), friends=1)
    original = RunState.to_json
    calls = []

    def record(self):
        calls.append(self)
        return original(self)

    monkeypatch.setattr(RunState, "to_json", record)
    result = run(state, catalog=catalog)
    assert result["summary"]["steps"] == 2
    assert len(calls) == 2


@pytest.mark.parametrize("value", [1, None, "yes"])
def test_trace_flag_requires_boolean(catalog, value):
    state = initialize(2, (1,))
    with pytest.raises(ControlError, match="boolean"):
        run(state, catalog=catalog, trace=value)
    with pytest.raises(ControlError, match="boolean"):
        step(state, catalog=catalog, trace=value)


def test_public_catalog_unbound_context_is_an_input_error(catalog):
    record = catalog.get("zip.far.positive.I")
    word = {"type": "concat", "items": [
        {"type": "word_ref", "name": "exterior_braid"},
        {"type": "slider", "size": "i", "start": "j"}]}
    # P X = P X is a valid identity schema. This controller requires a literal
    # braid context and cannot execute an unspecified named word P.
    record["lhs"] = record["rhs"] = [
        {"coefficient": [{"integer": 1, "q_power": 0}], "word": word}]
    custom = Catalog({"schema_version": "1.0", "rules": [record]})
    state = initialize(4, (3,))
    original = state.to_json()
    for operation in (step, run):
        with pytest.raises(CatalogError, match="Unbound word reference"):
            operation(state, catalog=custom)
        assert state.to_json() == original
