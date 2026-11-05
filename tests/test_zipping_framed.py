"""Retained flype frames, whole-tape contexts and exact contextual defects."""
from collections import defaultdict
from copy import deepcopy
import json

import pytest

from zipper_algebra.zipping_algorithm import framed
from zipper_algebra.zipping_algorithm.catalog import Catalog, CatalogError
from zipper_algebra.zipping_algorithm.interpreter import Branch, Limits, ReplayError, RunState


@pytest.fixture(scope="module")
def catalog():
    return Catalog.bundled()


def slider(i=3, j=1, mirrored=False):
    return ("slider", i, j, mirrored, False, False)


def test_positive_flype_records_full_frames_and_continues_mirrored(catalog):
    state = framed.initialize(3, (1, 2, 2, 1), processed=(-1,), tail=(2,))
    result = framed.run(state, catalog=catalog, trace=True)
    assert result["status"] == "framed_upper_pass_complete"
    assert result["summary"]["steps"] == 2
    assert result["summary"]["rule_counts"] == {"zip.flype.positive.I": 1, "zip.flype.positive.M": 1}
    branch, = RunState.from_json(result["final"]).branches
    assert branch.control["output_power"] == -2
    assert branch.control["input_power"] == 2
    assert branch.control["mirrored"] is False
    assert branch.word == (-2, -2, slider(), 1, 1, 2, 1, 1, 2, -1, 2)
    assert branch.control["processed_length"] == 5
    assert branch.control["tail_length"] == 1
    assert framed.check_certificate(json.loads(json.dumps(result)), catalog=catalog)


CASES = [
    (5, (4,), 1, 0, "zip.far.positive"),
    (5, (-4,), 1, 0, "zip.far.negative"),
    (6, (2,), 1, 3, "zip.close.positive"),
    (6, (-2,), 1, 3, "zip.close.negative"),
    (3, (2, 1), 1, 0, "zip.flype.positive"),
    (3, (-2, -1), 1, 0, "zip.flype.negative"),
    (4, (1,), 1, 1, "zip.nearest.positive"),
    (4, (-1,), 1, 1, "zip.nearest.negative.i_ge_3"),
    (4, (-1,), 2, 0, "zip.left.negative"),
    (4, (1,), 2, 0, "zip.left.positive.i_ge_3"),
    (4, (-3, -2), 1, 1, "zip.farthest.negative"),
    (4, (3, 2), 1, 1, "zip.farthest.positive.i_ge_3"),
    (4, (2, 3), 1, 0, "zip.outer.positive"),
    (4, (-2, -3), 1, 0, "zip.outer.negative.i_ge_3"),
]


@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("mirror", [False, True])
@pytest.mark.parametrize("a,b", [(-1, 2), (2, -1)])
def test_all_higher_rule_families_preserve_frame_coordinates(catalog, case, mirror, a, b):
    n, word, j, k, name = case
    if mirror:
        word = tuple((1 if x > 0 else -1) * (n - abs(x)) for x in word)
        j = n - 3 - j + 2
    state = framed.initialize(n, word, start=j, friends=k, mirrored=mirror,
                              output_power=a, input_power=b)
    following, event = framed.step(state, catalog=catalog)
    assert event["rule_id"] == name + (".M" if mirror else ".I")
    assert all(x.control["left_length"] == 0 for x in following.branches)
    assert all(x.control["output_power"] + x.control["input_power"] == a + b
               for x in following.branches)
    assert all(x.control["phase"] == "framed_upper_pass_complete" for x in following.branches)
    assert len(event["branch_frames"]) == len(following.branches)
    assert all(x["tail_transport"]["transformed"] == x["tail_transport"]["word"]
               for x in event["branch_frames"])


def _inverse(word):
    return tuple(-x for x in reversed(word))


def _substitute(word, images):
    stack = []
    for letter in word:
        image = images[letter - 1] if letter > 0 else _inverse(images[-letter - 1])
        for x in image:
            if stack and stack[-1] == -x:
                stack.pop()
            else:
                stack.append(x)
    return tuple(stack)


def _braid_key(word, n):
    images = tuple((x,) for x in range(1, n + 1))
    for letter in word:
        j = abs(letter)
        generators = [(x,) for x in range(1, n + 1)]
        generators[j - 1:j + 1] = ([(j, j + 1, -j), (j,)] if letter > 0
                                    else [(j + 1,), (-j - 1, j, j + 1)])
        images = tuple(_substitute(w, images) for w in generators)
    return images


def _slider_polynomial(size, start, mirrored):
    result = {(1, (-start,)): 1, (0, ()): 1, (1, ()): -1, (0, (start,)): -1}
    for i in range(3, size + 1):
        out = defaultdict(int)
        w = tuple(range(start, start + i - 1))
        for (q, word), c in result.items():
            out[(q + i - 1, tuple(-x for x in w) + word)] += c
            out[(q, w + word)] -= c
        result = dict(out)
    if mirrored:
        result = {(q, tuple((1 if x > 0 else -1) * (2 * start + size - 2 - abs(x))
                             for x in word)): c for (q, word), c in result.items()}
    return result


def _expanded(branches, n):
    result = defaultdict(int)
    for branch in branches:
        token_index = next(idx for idx, x in enumerate(branch.word) if isinstance(x, tuple))
        token = branch.word[token_index]
        for (q, middle), scalar in _slider_polynomial(token[1], token[2], token[3]).items():
            braid = branch.word[:token_index] + middle + branch.word[token_index + 1:]
            key = _braid_key(braid, n)
            for power, c in branch.coefficient:
                result[(power + q, key)] += scalar * c
    return {key: value for key, value in result.items() if value}


def _difference(a, b):
    result = defaultdict(int, a)
    for key, value in b.items():
        result[key] -= value
    return {key: value for key, value in result.items() if value}


@pytest.mark.parametrize("case", [CASES[4], CASES[5], CASES[7], CASES[9], CASES[11], CASES[13]])
@pytest.mark.parametrize("mirror", [False, True])
def test_full_contextual_defect_equals_catalog_defect_in_faithful_braid_algebra(catalog, case, mirror):
    n, word, j, k, _ = case
    if mirror:
        word = tuple((1 if x > 0 else -1) * (n - abs(x)) for x in word)
        j = n - 3 - j + 2
    state = framed.initialize(n, word, start=j, friends=k, mirrored=mirror,
                              output_power=-1, input_power=1, processed=(1,), tail=(-1,))
    following, event = framed.step(state, catalog=catalog)
    equality = event["catalog_transition"]
    lhs = (Branch.from_json(equality["predecessor"]),)
    rhs = tuple(Branch.from_json(x) for x in equality["successors"])
    # The free braid-group-algebra defects agree exactly. Thus the only
    # quotient identity used is the displayed, contextual catalog identity.
    actual = _difference(_expanded(state.branches, n), _expanded(following.branches, n))
    expected = _difference(_expanded(lhs, n), _expanded(rhs, n))
    assert actual == expected


def test_optional_trace_has_same_arithmetic_and_only_two_compact_snapshots(catalog, monkeypatch):
    state = framed.initialize(3, (1, 2, 2, 1))
    original = RunState.to_json
    calls = []

    def record(self):
        calls.append(self)
        return original(self)

    monkeypatch.setattr(RunState, "to_json", record)
    compact = framed.run(state, catalog=catalog)
    assert len(calls) == 2
    detailed = framed.run(state, catalog=catalog, trace=True)
    assert compact["final"] == detailed["final"]
    assert compact["summary"] == detailed["summary"]
    assert compact["transitions"] == []
    assert compact["metrics"] == [detailed["metrics"][0], detailed["metrics"][-1]]
    assert framed.check_certificate(compact, catalog=catalog)
    assert framed.check_certificate(detailed, catalog=catalog)


def test_internal_child_crossings_are_absorbed_without_algebraic_change(catalog):
    initial = framed.initialize(3, (2, -2), output_payload=(2,))
    following, event = framed.step(initial, catalog=catalog)
    assert following.branches[0].word == initial.branches[0].word
    assert following.branches[0].control["output_payload"] == [-2, 2]
    assert following.branches[0].control["left_length"] == 1
    assert following.branches[0].control["processed_events"] == []
    assert event["rule_id"] == "control.absorb_internal_tape"
    assert event["equality_reference"] == "identical_literal_word"
    result = framed.run(initial, catalog=catalog)
    assert result["status"] == "framed_upper_pass_complete"
    assert result["final"]["branches"][0]["control"]["output_payload"] == [2, -2, 2]
    assert result["summary"]["rule_counts"] == {"control.absorb_internal_tape": 2}
    assert framed.check_certificate(result, catalog=catalog)
    # An internal crossing above a flyped slider remains explicit in its child.
    result = framed.run(framed.initialize(3, (1, 2, 1)), catalog=catalog, trace=True)
    assert result["status"] == "framed_upper_pass_complete"
    assert result["summary"]["rule_counts"] == {"zip.flype.positive.I": 1, "control.absorb_internal_tape": 1}
    assert result["final"]["branches"][0]["control"]["output_payload"] == [1]


@pytest.mark.parametrize("mirror", [False, True])
def test_whole_tape_moves_transport_the_complete_child_payload(catalog, mirror):
    n, j, word, payload = (5, 3, (3, 2), (3, -3)) if mirror else (5, 1, (2, 3), (2, -2))
    initial = framed.initialize(n, word, start=j, mirrored=mirror, output_payload=payload)
    following, event = framed.step(initial, catalog=catalog)
    expected = [2, -2] if mirror else [3, -3]
    assert following.branches[0].control["output_payload"] == expected
    assert event["payload_transport"]["transformed"] == expected
    raw = event["catalog_transition"]
    assert _difference(_expanded(initial.branches, n), _expanded(following.branches, n)) == _difference(
        _expanded((Branch.from_json(raw["predecessor"]),), n),
        _expanded(tuple(Branch.from_json(b) for b in raw["successors"]), n))


def test_noncommuting_three_strand_child_braid_is_retained(catalog):
    initial = framed.initialize(5, (3, 2, 3, 4), size=4, output_payload=(2, 3, -2))
    result = framed.run(initial, catalog=catalog)
    assert result["status"] == "framed_upper_pass_complete"
    assert result["summary"]["steps"] == 2
    assert result["final"]["branches"][0]["control"]["output_payload"] == [3, 3, 4, -3]
    assert framed.check_certificate(result, catalog=catalog)


def test_limits_and_missing_tape_preserve_checkpoint(catalog):
    initial = framed.initialize(3, (1, 2, 2, 1))
    result = framed.run(initial, catalog=catalog, limits=Limits(max_steps=1), trace=True)
    assert result["status"] == "incomplete"
    assert result["summary"]["steps"] == 1
    assert framed.check_certificate(result, catalog=catalog)
    resumed = framed.resume(result, catalog=catalog, trace=False)
    assert resumed["status"] == "framed_upper_pass_complete"
    assert resumed["trace_enabled"] is False
    missing = framed.run(framed.initialize(3, (1,)), catalog=catalog)
    assert missing["status"] == "missing_rule"
    assert missing["summary"]["steps"] == 0
    assert missing["diagnostic"]["message"] == "missing framed full-tape rule"
    assert framed.check_certificate(missing, catalog=catalog)
    continued = framed.resume(missing, catalog=catalog)
    assert continued["status"] == "missing_rule"


def test_branch_and_size_limits_are_atomic(catalog):
    state = framed.initialize(4, (-1,), friends=1)
    for limits in [Limits(max_terms=2), Limits(max_word_size=3), Limits(expansion_budget=1)]:
        result = framed.run(state, catalog=catalog, limits=limits)
        assert result["status"] == "incomplete"
        assert result["final"] == state.to_json()
        assert result["transitions"] == []
        assert framed.check_certificate(result, catalog=catalog)


def test_frame_expansion_budget_precedes_allocation():
    with pytest.raises(ReplayError, match="word limit"):
        framed.initialize(3, (), output_power=10**9, limits=Limits(max_word_size=10))
    with pytest.raises(ReplayError, match="frame word limit"):
        framed._power(4, 1, 20, 10)


def test_two_strand_delegate_and_empty_polynomial(catalog):
    state = framed.initialize(2, (1,), size=2, output_power=999, input_power=-999)
    result = framed.run(state, catalog=catalog)
    assert result["schema"] == "zipper-normal-upward-v1"
    assert framed.check_certificate(result, catalog=catalog)
    assert framed.resume(result, catalog=catalog)["status"] == "upper_pass_complete"
    assert framed.run(RunState(3, ()), catalog=catalog)["status"] == "framed_upper_pass_complete"
    complete = framed.initialize(3, ())
    assert framed.run(complete, catalog=catalog)["status"] == "framed_upper_pass_complete"
    with pytest.raises(ReplayError, match="already complete"):
        framed.step(complete, catalog=catalog)


@pytest.mark.parametrize("kwargs", [
    {"n": True}, {"size": True}, {"size": 1}, {"size": 5}, {"start": 0},
    {"start": 4}, {"friends": -1}, {"friends": 3}, {"mirrored": 1},
    {"mirrored": True, "friends": 1}, {"output_power": True}, {"input_power": 0.5},
    {"left": (False,)}, {"processed": (slider(),)}, {"tail": (True,)},
    {"left": (0,)}, {"size": 2, "mirrored": True},
    {"output_payload": (1,)}, {"output_payload": (True,)},
    {"size": 2, "output_payload": (2,)},
])
def test_invalid_initialization(kwargs):
    args = dict(n=4, left=())
    args.update(kwargs)
    with pytest.raises(ReplayError):
        framed.initialize(**args)


@pytest.mark.parametrize("branch", [-1, 1, True])
def test_invalid_branch(catalog, branch):
    with pytest.raises(ReplayError, match="does not exist"):
        framed.step(framed.initialize(3, (2, 1)), branch, catalog=catalog)


def test_invalid_trace_flag(catalog):
    state = framed.initialize(3, (2, 1))
    with pytest.raises(ReplayError, match="boolean"):
        framed.run(state, catalog=catalog, trace=1)
    with pytest.raises(ReplayError, match="boolean"):
        framed.step(state, catalog=catalog, trace=1)


def test_certificate_tampering_rejected(catalog):
    certificate = framed.run(framed.initialize(3, (2, 1)), catalog=catalog, trace=True)
    for key, value in [("source_sha256", "0" * 64), ("status", "complete"), ("extra", True)]:
        altered = deepcopy(certificate)
        altered[key] = value
        with pytest.raises(ReplayError, match="differs"):
            framed.check_certificate(altered, catalog=catalog)
    altered = deepcopy(certificate)
    altered["transitions"][0]["branch_frames"][0]["emitted_output_power"] = 0
    with pytest.raises(ReplayError, match="differs"):
        framed.check_certificate(altered, catalog=catalog)


def test_invalid_boundaries_rejected(catalog):
    state = framed.initialize(3, (2, 1))
    for changes in [{"extra": 1}, {"mirrored": 1}, {"i": True}, {"i": 2}, {"j": 0},
                    {"k": -1}, {"j": 2}, {"k": 1}, {"left_length": -1}, {"phase": "other"},
                    {"output_power": 1}]:
        c = state.branches[0].control
        c.update(changes)
        malformed = RunState(3, (Branch(state.branches[0].word, control_json=json.dumps(c)),))
        with pytest.raises(ReplayError):
            framed.run(malformed, catalog=catalog)
    c = framed.initialize(3, (), mirrored=True).branches[0].control
    c["k"] = 1
    with pytest.raises(ReplayError):
        framed.run(RunState(3, (Branch((slider(mirrored=True),), control_json=json.dumps(c)),)), catalog=catalog)
    c = state.branches[0].control
    malformed = RunState(3, (Branch((slider(), 1, slider()), control_json=json.dumps(c)),))
    with pytest.raises(ReplayError, match="pure braid"):
        framed.run(malformed, catalog=catalog)


@pytest.mark.parametrize("field,value", [
    ("output_payload", None), ("output_payload", [1]), ("output_payload", [True]),
    ("initial_slider", None), ("initial_slider", {}),
    ("initial_slider", {"j": 1, "k": 0, "mirrored": 1, "output_power": 0, "input_power": 0, "processed_length": 0}),
    ("initial_slider", {"j": True, "k": 0, "mirrored": False, "output_power": 0, "input_power": 0, "processed_length": 0}),
    ("processed_events", None), ("processed_events", [None]), ("processed_events", [{}]),
    ("processed_events", [{"rule_id": 1, "parameters": {}, "rhs_index": 0}]),
    ("processed_events", [{"rule_id": "zip.flype.positive.I", "parameters": None, "rhs_index": 0}]),
    ("processed_events", [{"rule_id": "zip.flype.positive.I", "parameters": {"n": True}, "rhs_index": 0}]),
    ("processed_events", [{"rule_id": "zip.flype.positive.I", "parameters": {}, "rhs_index": -1}]),
    ("processed_events", [{"rule_id": "zip.flype.positive.I", "parameters": {}, "rhs_index": True}]),
])
def test_saved_band_records_reject_malformed_inputs(catalog, field, value):
    state = framed.initialize(3, (2, 1))
    control = state.branches[0].control
    control[field] = value
    malformed = RunState(3, (Branch(state.branches[0].word, control_json=json.dumps(control)),))
    with pytest.raises(ReplayError):
        framed.run(malformed, catalog=catalog)


@pytest.mark.parametrize("case", CASES)
def test_semilinear_mirrored_policy_all_families(catalog, case):
    n, word, j, k, name = case
    word = tuple(-(1 if x > 0 else -1) * (n - abs(x)) for x in word)
    state = framed.initialize(n, word, start=n - 3 - j + 2, friends=k, mirrored=True,
                              output_power=-1, input_power=1, mirrored_policy="bar_mirror")
    following, event = framed.step(state, catalog=catalog)
    assert event["rule_id"] == name + ".BM"
    assert all(b.control["mirror_policy"] == "bar_mirror" for b in following.branches)
    assert all(b.control["left_length"] == 0 for b in following.branches)
    for branch in following.branches:
        c = branch.control
        assert c["processed_events"][0]["mirror_policy"] == "bar_mirror"
        if c["mirrored"]:
            assert all(x < 0 for x in framed.friend_word(n, c))
        else:
            assert c["k"] == 0
        restored = framed.rule_for_event(catalog, c["processed_events"][0])
        assert restored["id"] == event["rule_id"]
        assert restored["source_record_id"] == name + ".M"
    equality = event["catalog_transition"]
    assert _difference(_expanded(state.branches, n), _expanded(following.branches, n)) == _difference(
        _expanded((Branch.from_json(equality["predecessor"]),), n),
        _expanded(tuple(Branch.from_json(x) for x in equality["successors"]), n))


def test_semilinear_policy_inverts_q_and_keeps_negative_friends(catalog):
    state = framed.initialize(4, (3,), start=2, friends=1, mirrored=True,
                              mirrored_policy="bar_mirror")
    following, event = framed.step(state, catalog=catalog)
    assert event["rule_id"] == "zip.nearest.negative.i_ge_3.BM"
    assert [b.coefficient for b in following.branches] == [((-2, -1),), ((0, 1),), ((-2, 1),)]
    assert framed.friend_word(4, state.branches[0].control) == (-2, -1)
    certificate = framed.run(state, catalog=catalog, trace=True)
    assert framed.check_certificate(json.loads(json.dumps(certificate)), catalog=catalog)
    compact = framed.run(state, catalog=catalog)
    assert compact["final"] == certificate["final"]
    assert compact["summary"] == certificate["summary"]
    assert framed.check_certificate(compact, catalog=catalog)


def test_normal_to_semilinear_mirror_arrives_without_friends(catalog):
    state = framed.initialize(4, (-1, -2, 2, 1), mirrored_policy="bar_mirror")
    following, event = framed.step(state, catalog=catalog)
    assert event["rule_id"] == "zip.flype.positive.I"
    c = following.branches[0].control
    assert c["mirrored"] is True and c["k"] == 0 and c["mirror_policy"] == "bar_mirror"
    # The selected policy is retained across orientation changes and checkpoints.
    certificate = framed.run(following, catalog=catalog, limits=Limits(max_steps=1), trace=True)
    assert certificate["transitions"][0]["rule_id"].endswith(".BM")
    assert framed.check_certificate(certificate, catalog=catalog)


def test_invalid_mirrored_policy_and_saved_policy_mismatch(catalog):
    with pytest.raises(ReplayError):
        framed.initialize(3, (), mirrored_policy="bar_with_q_fixed")
    state = framed.initialize(3, (2, 1))
    control = state.branches[0].control
    control["mirror_policy"] = "bar_mirror"
    malformed = RunState(3, (Branch(state.branches[0].word, control_json=json.dumps(control)),))
    with pytest.raises(ReplayError, match="input-band record"):
        framed.run(malformed, catalog=catalog)


def _single_record(record):
    return Catalog({"schema_version": "1.0", "rules": [record]})


def _unit_pair(index=1):
    return [{"type": "sigma", "index": index, "sign": 1},
            {"type": "sigma", "index": index, "sign": -1}]


def test_original_event_rule_is_recovered_without_semilinear_transformation(catalog):
    event = {"rule_id": "zip.flype.positive.I", "parameters": {"n": 3, "i": 3, "j": 1, "k": 0}}
    assert framed.rule_for_event(catalog, event) == catalog.get(event["rule_id"])


def test_public_catalog_rejects_an_unbound_named_context(catalog):
    record = catalog.get("zip.far.positive.I")
    word = {"type": "concat", "items": [
        {"type": "word_ref", "name": "external_braid"},
        {"type": "slider", "size": "i", "start": "j"}]}
    # This is the literal identity P X = P X. The framed entry point requires
    # its braid context explicitly and cannot bind an unspecified named P.
    record["lhs"] = record["rhs"] = [{"coefficient": [{"integer": 1, "q_power": 0}], "word": word}]
    state = framed.initialize(5, (4,))
    before = state.to_json()
    with pytest.raises(CatalogError, match="Unbound word reference"):
        framed.step(state, catalog=_single_record(record))
    assert state.to_json() == before


def test_public_catalog_rejects_identity_with_a_different_active_core(catalog):
    record = catalog.get("zip.far.positive.I")
    word = {"type": "concat", "items": [
        {"type": "sigma", "index": "l", "sign": 1},
        {"type": "slider", "size": "i", "start": {"op": "add", "args": ["j", 1]}}]}
    # An exact identity with another slider anchor cannot replace the active
    # marked core merely because its numerical guard is satisfied.
    record["lhs"] = record["rhs"] = [{"coefficient": [{"integer": 1, "q_power": 0}], "word": word}]
    with pytest.raises(ReplayError, match="left side does not retain the framed core"):
        framed.step(framed.initialize(5, (4,)), catalog=_single_record(record))


def test_public_catalog_identity_must_preserve_the_output_tape(catalog):
    record = catalog.get("zip.outer.positive.I")
    # The identity W X = W X is algebraically valid, but it does not perform
    # the advertised tape crossing or advance the marked output interval.
    record["rhs"] = deepcopy(record["lhs"])
    with pytest.raises(ReplayError, match="output tape disagrees"):
        framed.step(framed.initialize(4, (2, 3)), catalog=_single_record(record))


def test_public_catalog_requires_the_explicit_flype_frame_shape(catalog):
    record = catalog.get("zip.flype.positive.I")
    # Inserting sigma_1 sigma_1^-1 preserves this genuine braid identity.
    # Its additional exterior letters require a separate control witness.
    record["rhs"][0]["word"] = {"type": "concat", "items": _unit_pair() + [record["rhs"][0]["word"]]}
    with pytest.raises(ReplayError, match="missing its retained flype frames"):
        framed.step(framed.initialize(3, (2, 1)), catalog=_single_record(record))


def test_public_catalog_rejects_an_unrecorded_internal_tail_conjugation(catalog):
    record = catalog.get("zip.far.positive.I")
    # Appending an inverse pair preserves far commutativity. Transport across
    # the odd three-strand input frame reflects this internal pair, so the
    # declared literal tail is not a supported marked-tape tail.
    record["rhs"][0]["word"] = {"type": "concat", "items": [record["rhs"][0]["word"]] + _unit_pair()}
    with pytest.raises(ReplayError, match="tail does not preserve the marked input tape"):
        framed.step(framed.initialize(6, (5,), size=4, input_power=1), catalog=_single_record(record))
