"""Turnaround preserves each whole generic-q branch, including flype frames."""

from collections import defaultdict
from copy import deepcopy
from dataclasses import replace
import json

import pytest

from zipper_algebra.zipping_algorithm import framed, turnaround
from zipper_algebra.zipping_algorithm.catalog import Catalog
from zipper_algebra.zipping_algorithm.interpreter import Branch, Limits, RunState


def _inverse(word):
    return tuple(-x for x in reversed(word))


def _substitute(word, images):
    stack = []
    for x in word:
        for y in images[x - 1] if x > 0 else _inverse(images[-x - 1]):
            if stack and stack[-1] == -y:
                stack.pop()
            else:
                stack.append(y)
    return tuple(stack)


def _key(word, n):
    images = tuple((i,) for i in range(1, n + 1))
    for letter in word:
        j = abs(letter)
        elementary = [(i,) for i in range(1, n + 1)]
        elementary[j - 1 : j + 1] = (
            [(j, j + 1, -j), (j,)] if letter > 0 else [(j + 1,), (-j - 1, j, j + 1)]
        )
        images = tuple(_substitute(w, images) for w in elementary)
    return images


def _slider(size, start, mirrored):
    """Independent defining expansion with no catalog or turnaround helper."""
    out = {(1, (-start,)): 1, (0, ()): 1, (1, ()): -1, (0, (start,)): -1}
    for level in range(3, size + 1):
        nxt = defaultdict(int)
        positive = tuple(range(start, start + level - 1))
        for (q, word), coefficient in out.items():
            nxt[(q + level - 1, tuple(-x for x in positive) + word)] += coefficient
            nxt[(q, positive + word)] -= coefficient
        out = dict(nxt)
    if mirrored:
        out = {
            (
                q,
                tuple(
                    (1 if x > 0 else -1) * (2 * start + size - 2 - abs(x)) for x in word
                ),
            ): c
            for (q, word), c in out.items()
        }
    return out


def _expanded(branch, n):
    at = next(
        index for index, token in enumerate(branch.word) if isinstance(token, tuple)
    )
    _, size, start, mirror, reverse, bar = branch.word[at]
    assert not reverse and not bar
    out = defaultdict(int)
    for (q, word), c in _slider(size, start, mirror).items():
        braid = branch.word[:at] + word + branch.word[at + 1 :]
        key = _key(braid, n)
        for p, d in branch.coefficient:
            out[(q + p, key)] += c * d
    return {key: value for key, value in out.items() if value}


def _complete(n, left, **options):
    state = framed.initialize(n, left, **options)
    report = framed.run(state, trace=False)
    assert report["status"] in ("framed_upper_pass_complete", "upper_pass_complete")
    return RunState.from_json(report["final"]).branches


@pytest.mark.parametrize("mirror", [False, True])
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("powers", [(-1, 2), (0, 0), (2, -1)])
def test_single_flype_exact_braid_group_algebra_with_frames_and_exterior(
    mirror, sign, powers
):
    a, b = powers
    word = tuple(sign * x for x in ((1, 2) if mirror else (2, 1)))
    (branch,) = _complete(
        3,
        word,
        size=3,
        start=1,
        mirrored=mirror,
        output_power=a,
        input_power=b,
        processed=(-1, 1),
        tail=(1, -2),
    )
    following, cert = turnaround.turnaround(branch, n=3, trace=True)
    assert _expanded(branch, 3) == _expanded(following, 3)
    control, exterior, processed, tail = turnaround.parts(following, n=3)
    assert exterior == word
    assert processed == (-1, 1) and tail == (1, -2)
    assert (
        control["j"],
        control["mirrored"],
        control["output_power"],
        control["input_power"],
    ) == (1, mirror, a, b)
    assert turnaround.check_certificate(json.loads(json.dumps(cert)))


def test_multiple_flypes_are_unwrapped_newest_first_without_reversing_branching():
    (branch,) = _complete(3, (2, 1, 1, 2, 2, 1), size=3, start=1, tail=(-1, 2))
    following, cert = turnaround.turnaround(branch, n=3)
    assert len(branch.control["processed_events"]) == 3
    assert following.control["output_power"] == following.control["input_power"] == 0
    assert cert["events"] == [] and cert["summary"]["events"] == 3
    assert _expanded(branch, 3) == _expanded(following, 3)
    assert turnaround.check_certificate(cert)
    with pytest.raises(ValueError, match="event limit"):
        turnaround.turnaround(branch, n=3, limits=Limits(max_steps=2))


def test_higher_nested_child_payload_with_repeated_flypes():
    word = (3, 2, 1, 1, 2, 3, 3, 2, 1, 2)
    (branch,) = _complete(4, word, size=4, start=1)
    following, cert = turnaround.turnaround(branch, n=4)
    assert len(branch.control["processed_events"]) == 3
    assert _expanded(branch, 4) == _expanded(following, 4)
    assert turnaround.check_certificate(cert)


def test_two_strand_internal_tails_and_exceptional_branches():
    for left, friends in [((1, -1), 0), ((-2,), 0), ((-1,), 1), ((2,), 1), ((1,), 1)]:
        branches = _complete(3, left, size=2, start=1, friends=friends, tail=(2, -1))
        for branch in branches:
            following, cert = turnaround.turnaround(branch, n=3, trace=True)
            assert _expanded(branch, 3) == _expanded(following, 3)
            assert turnaround.parts(following, n=3)[3] == (2, -1)
            assert turnaround.check_certificate(cert)
    unfinished = framed.initialize(3, (1,), size=2)
    with pytest.raises(ValueError, match="completed"):
        turnaround.turnaround(unfinished.branches[0], n=3)


def test_farthest_positive_branches_and_full_band_tails():
    for mirror in (False, True):
        word = (3, 2) if not mirror else (2, 3)
        start = 1 if not mirror else 3
        branches = _complete(
            5,
            word,
            size=3,
            start=start,
            friends=1,
            mirrored=mirror,
            output_power=1,
            input_power=-1,
            tail=(-4, 4),
        )
        assert len(branches) == 3
        for branch in branches:
            following, cert = turnaround.turnaround(branch, n=5, trace=True)
            assert _expanded(branch, 5) == _expanded(following, 5)
            assert turnaround.parts(following, n=5)[3] == (-4, 4)
            assert turnaround.check_certificate(cert)


def test_zero_tail_events_and_friend_exterior_are_retained():
    for left in ((-3, -2), (1,)):
        (branch,) = _complete(5, left, size=3, start=1, friends=1)
        following, _ = turnaround.turnaround(branch, n=5)
        assert _expanded(branch, 5) == _expanded(following, 5)
    (branch,) = _complete(5, (), size=3, start=1, friends=1)
    following, _ = turnaround.turnaround(branch, n=5)
    assert turnaround.parts(following, n=5)[1] == (2, 3)
    assert (
        following.control["payload_length"],
        following.control["routing_length"],
        following.control["transported_length"],
    ) == (0, 2, 0)


def test_gathered_child_payload_remains_in_the_exterior():
    for word in ((2,), (2, 1, 2)):
        (branch,) = _complete(3, word, size=3, start=1)
        assert branch.control["output_payload"]
        following, cert = turnaround.turnaround(branch, n=3)
        assert _expanded(branch, 3) == _expanded(following, 3)
        assert turnaround.parts(following, n=3)[1][:1] == tuple(
            branch.control["output_payload"]
        )
        assert following.control["payload_length"] == len(
            branch.control["output_payload"]
        )
        assert following.control["routing_length"] == 0
        assert turnaround.check_certificate(cert)


def test_bar_mirror_policy_retains_negative_friends_and_unwinds_actual_packets():
    cases = [
        (4, (3, 2, 1), dict(size=3, start=1)),
        (5, (3,), dict(size=3, start=2, friends=1, mirrored=True)),
        (5, (-2, -3), dict(size=3, start=3, friends=1, mirrored=True)),
        (3, (2, 1, 1, 2, 2, 1), dict(size=3, start=1)),
    ]
    for n, left, options in cases:
        branches = _complete(
            n,
            left,
            mirrored_policy="bar_mirror",
            output_power=-1,
            input_power=2,
            tail=(1, -2),
            **options,
        )
        for branch in branches:
            following, cert = turnaround.turnaround(branch, n=n, trace=True)
            assert _expanded(branch, n) == _expanded(following, n)
            assert following.control["mirror_policy"] == "bar_mirror"
            assert turnaround.check_certificate(cert)
    (branch,) = _complete(4, (3, 2, 1), size=3, start=1, mirrored_policy="bar_mirror")
    following, _ = turnaround.turnaround(branch, n=4)
    control, exterior, _, _ = turnaround.parts(following, n=4)
    assert control["routing_length"] == 2
    assert exterior[:2] == (-2, -1)


def _with_control(branch, change):
    control = branch.control
    change(control)
    return replace(branch, control_json=json.dumps(control))


def test_rejects_forged_processed_events_and_modified_words():
    (branch,) = _complete(3, (2, 1), size=3, start=1)
    changes = [
        lambda c: c.update(processed_events=[{}]),
        lambda c: c["processed_events"][0].update(rhs_index=True),
        lambda c: c["processed_events"][0]["parameters"].update(i=2),
        lambda c: c["processed_events"][0]["parameters"].update(j=0),
        lambda c: c["processed_events"][0].update(rhs_index=9),
        lambda c: c["processed_events"][0].update(rule_id="layer.positive.I"),
        lambda c: c["initial_slider"].update(output_power=100),
        lambda c: c.update(processed_events=[]),
        lambda c: c["processed_events"][0].update(rule_id="zip.flype.negative.I"),
    ]
    for change in changes:
        with pytest.raises(ValueError):
            turnaround.turnaround(_with_control(branch, change), n=3)
    unfinished = framed.initialize(3, (2, 1), size=3)
    with pytest.raises(ValueError, match="completed"):
        turnaround.turnaround(unfinished.branches[0], n=3)
    with pytest.raises(ValueError, match="boolean"):
        turnaround.turnaround(branch, n=3, trace=1)
    with pytest.raises(ValueError):
        turnaround.turnaround(branch, n=3, limits=Limits(max_word_size=1))


def test_saved_certificate_tampering_and_output_boundary_validation():
    (branch,) = _complete(3, (2, 1), size=3)
    following, cert = turnaround.turnaround(branch, n=3)
    cert["schema"] = "wrong"
    with pytest.raises(ValueError, match="schema"):
        turnaround.check_certificate(cert)
    cert["schema"] = "zipper-turnaround-v1"
    cert["source_sha256"] = "changed"
    with pytest.raises(ValueError, match="replay"):
        turnaround.check_certificate(cert)
    for change in (
        lambda c: c.update(exterior_length=True),
        lambda c: c.update(exterior_length=0),
    ):
        with pytest.raises(ValueError):
            turnaround.parts(_with_control(following, change), n=3)
    with pytest.raises(ValueError):
        turnaround.parts(Branch((1,)), n=3)


def _changed_catalog(name, change):
    """A public, schema-valid custom catalog with explicit noncanonical input."""
    records = list(Catalog.bundled().rules)
    record = next(row for row in records if row["id"] == name)
    change(record)
    return Catalog({"schema_version": "1.0", "rules": records})


def _packet_branch():
    # The tail references a flype packet from another accepted local rule,
    # so changing that packet tests the standalone checker independently
    # of the producer that emitted the original branch.
    return _complete(5, (3, 2), size=3, start=1, friends=1)[0]


def _cancellation():
    return [
        {"type": "sigma", "index": 1, "sign": 1},
        {"type": "sigma", "index": 1, "sign": -1},
    ]


def test_public_custom_packet_refuses_nonmonomial_or_scaled_identity():
    branch = _packet_branch()

    def add_same_unit(record):
        unit = {
            "coefficient": [{"integer": 1, "q_power": 0}],
            "word": {"type": "identity"},
        }
        record["lhs"].append(deepcopy(unit))
        record["rhs"].append(deepcopy(unit))

    catalog = _changed_catalog("zip.flype.positive.I", add_same_unit)
    with pytest.raises(ValueError, match="monomial packet"):
        turnaround.turnaround(branch, n=5, catalog=catalog)

    def scale_identity(record):
        for side in ("lhs", "rhs"):
            for term in record[side]:
                for coefficient in term["coefficient"]:
                    coefficient["integer"] *= 2

    catalog = _changed_catalog("zip.flype.positive.I", scale_identity)
    with pytest.raises(ValueError, match="coefficient one"):
        turnaround.turnaround(branch, n=5, catalog=catalog)


def test_public_custom_packet_refuses_equivalent_noncanonical_braid_spellings():
    branch = _packet_branch()

    def prefix_cancel(record):
        old = record["rhs"][0]["word"]
        record["rhs"][0]["word"] = {"type": "concat", "items": _cancellation() + [old]}

    catalog = _changed_catalog("zip.flype.positive.I", prefix_cancel)
    with pytest.raises(ValueError, match="retained frame"):
        turnaround.turnaround(branch, n=5, catalog=catalog)

    def cancel_before_outgoing_tape(record):
        items = record["rhs"][0]["word"]["items"]
        items[-1:-1] = _cancellation()

    catalog = _changed_catalog("zip.flype.positive.I", cancel_before_outgoing_tape)
    with pytest.raises(ValueError, match="outgoing flype tape"):
        turnaround.turnaround(branch, n=5, catalog=catalog)
    catalog = _changed_catalog("zip.farthest.positive.i_ge_3.I", prefix_cancel)
    with pytest.raises(ValueError, match="catalog frame shape"):
        turnaround.turnaround(branch, n=5, catalog=catalog)


def test_public_custom_packet_requires_the_declared_input_slider_orientation():
    branch = _packet_branch()
    other = Catalog.bundled().get("zip.flype.positive.M")

    def different_valid_identity(record):
        # At this reflected center anchor this is a valid flype identity,
        # but it has the opposite source slider orientation.
        record["lhs"], record["rhs"] = other["lhs"], other["rhs"]

    catalog = _changed_catalog("zip.flype.positive.I", different_valid_identity)
    with pytest.raises(ValueError, match="input slider"):
        turnaround.turnaround(branch, n=5, catalog=catalog)


def test_custom_catalog_acceptance_flag_does_not_authorize_missing_flype_frames():
    branch = _packet_branch()

    def missing_output_frame(record):
        record["rhs"][0]["word"]["items"][0]["items"].pop(0)

    catalog = _changed_catalog("zip.farthest.positive.i_ge_3.I", missing_output_frame)
    with pytest.raises(ValueError, match="frame increments"):
        turnaround.turnaround(branch, n=5, catalog=catalog)

    def missing_both_frames(record):
        items = record["rhs"][0]["word"]["items"][0]["items"]
        items.pop(-1)
        items.pop(0)

    catalog = _changed_catalog("zip.farthest.positive.i_ge_3.I", missing_both_frames)
    with pytest.raises(ValueError, match="explicit flype packet"):
        turnaround.turnaround(branch, n=5, catalog=catalog)


def test_reloaded_checkpoint_validates_complete_event_bindings_and_output_word():
    branch = _packet_branch()
    incomplete_binding = _with_control(
        branch, lambda c: c["processed_events"][0]["parameters"].pop("eps")
    )
    with pytest.raises(ValueError, match="processed-tail event"):
        turnaround.turnaround(incomplete_binding, n=5)
    following, _ = turnaround.turnaround(branch, n=5)
    word = list(following.word)
    at = next(i for i, token in enumerate(word) if isinstance(token, tuple))
    token = word[at]
    word[at] = (*token[:3], not token[3], *token[4:])
    with pytest.raises(ValueError, match="full expression"):
        turnaround.parts(replace(following, word=tuple(word)), n=5)


def test_forged_checkpoint_cannot_change_a_tail_return_anchor():
    branch = _packet_branch()

    def truncate_residual_band(record):
        # A forged accepted-source record retains the packet but omits
        # the full-band return word. Its syntactic validity is insufficient.
        record["rhs"][0]["word"]["items"] = record["rhs"][0]["word"]["items"][:2]

    catalog = _changed_catalog("zip.farthest.positive.i_ge_3.I", truncate_residual_band)
    c = branch.control
    begin = len(branch.word) - c["processed_length"] - c["tail_length"]
    shortened_tail = (2, 3)
    c["processed_length"] = len(shortened_tail)
    forged = replace(
        branch, word=branch.word[:begin] + shortened_tail, control_json=json.dumps(c)
    )
    with pytest.raises(ValueError, match="previous band anchor"):
        turnaround.turnaround(forged, n=5, catalog=catalog)
