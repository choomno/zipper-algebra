"""Literal-word scheduling, exact control identities, and complete rank-three runs."""

from collections import defaultdict
from copy import deepcopy
from itertools import product
import json
from math import comb

import pytest

from zipper_algebra.zipping_algorithm.algorithm import (
    AlgorithmError,
    AlgorithmSession,
    ControlObligation,
    _common,
    _extract_step,
    _layer_step,
    check_certificate,
    _pass_work,
    _zip_step,
    reduce_word,
    resume,
    work_expression,
)
from zipper_algebra.zipping_algorithm import framed
from zipper_algebra.zipping_algorithm.catalog import Catalog
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.extraction import RibbonLimit
from tests.test_zipping_transport import _artin_key, _slider


def _expanded(word, n):
    polynomial = {(0, ()): 1}
    for token in word:
        if type(token) is int:
            factor = {(0, (token,)): 1}
        else:
            assert not token[4] and not token[5]
            factor = _slider(token[1], token[2], token[3])
        out = defaultdict(int)
        for (q, u), a in polynomial.items():
            for (p, v), b in factor.items():
                out[q + p, u + v] += a * b
        polynomial = {key: c for key, c in out.items() if c}
    result = defaultdict(int)
    for (q, w), c in polynomial.items():
        result[q, _artin_key(w, n)] += c
    return {key: c for key, c in result.items() if c}


def _expression_sum(items, n, *, pending):
    out = defaultdict(int)
    for item in items:
        word = (
            work_expression(item)
            if pending
            else tuple(tuple(x) if isinstance(x, list) else x for x in item["word"])
        )
        for (q, key), c in _expanded(word, n).items():
            for p, d in item["coefficient"]:
                out[q + p, key] += c * d
    return {key: c for key, c in out.items() if c}


def _root_common(n):
    return {
        "rank": n,
        "coefficient": [[0, 1]],
        "contexts": [],
        "world_phi": False,
        "boundary_left": [],
        "boundary_right": [],
        "boundary_left_factors": [],
        "boundary_right_factors": [],
        "depth": 0,
        "bottom_labels": list(range(1, n + 1)),
        "top_labels": list(range(1, n + 1)),
    }


def test_every_literal_three_strand_word_through_length_two_completes():
    session = AlgorithmSession()
    count = 0
    for length in range(3):
        for word in product((-2, -1, 1, 2), repeat=length):
            result = session.reduce(word, n=3)
            assert result["status"] == "complete", (word, result["diagnostic"])
            assert result["input_word"] == list(word)
            assert not result["pending"] and not result["transitions"]
            assert all(len(term["components"]) <= 3 for term in result["output"])
            for term in result["output"]:
                for component in term["components"]:
                    size, record = component["size"], component["record"]
                    assert record["orientation"] == "normal"
                    assert all(
                        type(token) is int or token[3:] == [False, False, False]
                        for token in component["component_word"]
                    )
                    # These are the projected internal collars and residual
                    # comb routes, not the ambient endpoint-overcrossing routes.
                    bound = comb(size, 3) + comb(size, 2) - 1 if size > 1 else 0
                    assert len(record["top_route"]) <= bound
                    assert len(record["bottom_route"]) <= bound
                assert sorted(
                    p
                    for c in term["components"]
                    for p in c["record"]["original_bottom_ports"]
                ) == [1, 2, 3]
                assert sorted(
                    p
                    for c in term["components"]
                    for p in c["record"]["original_top_ports"]
                ) == [1, 2, 3]
            count += 1
    assert count == 21


def test_real_three_strand_trace_has_growth_and_exact_reloaded_control_identities():
    result = reduce_word((1, 2), n=3, trace=True)
    assert result["status"] == "complete"
    assert result["metrics"]["rule_applications"]["layer.route.left"] > 0
    assert any(c["size"] == 3 for t in result["output"] for c in t["components"])
    assert check_certificate(json.loads(json.dumps(result)))
    checked = set()
    for transition in result["transitions"]:
        evidence = transition["evidence"]
        name = evidence.get("rule_id", "")
        if name.startswith(("control.", "single.", "layer.")):
            before = _expression_sum([transition["predecessor"]], 3, pending=True)
            after = _expression_sum(transition["successors"], 3, pending=True)
            extra = _expression_sum(transition["output_added"], 3, pending=False)
            combined = defaultdict(int, after)
            for key, value in extra.items():
                combined[key] += value
            assert before == {key: value for key, value in combined.items() if value}, (
                name
            )
            checked.add(name)
    assert {"control.extract", "control.turnaround", "layer.route.left"} <= checked


@pytest.mark.parametrize("mirrored", [False, True])
@pytest.mark.parametrize("side", ["left", "right"])
@pytest.mark.parametrize("size", [2, 3, 4])
def test_four_layer_routing_identities_are_generic_braid_group_equalities(
    size, side, mirrored
):
    n = size + 1
    start = 2 if side == "left" else 1
    sign = 1 if side == "left" else -1
    crossing = tuple(
        sign * x
        for x in (
            range(start + size - 2, start - 2, -1)
            if side == "left"
            else range(start, start + size)
        )
    )
    slider = ("slider", size, start, mirrored, False, False)
    state = dict(
        _root_common(n),
        phase="layer",
        word=list(crossing) + [list(slider)],
        payload_length=0,
        completed_sides=0,
    )
    children, event = _layer_step(state, Limits())
    assert len(children) == 2 and event["offensive"]
    assert set(event["growth_child_type"].values()) == {mirrored}
    assert _expression_sum([state], n, pending=True) == _expression_sum(
        children, n, pending=True
    )
    assert children[0]["phase"] == "zip_up"


def test_growth_retains_the_ordered_old_child_stream_and_transports_only_payload():
    # The gathered sigma_2 is internal to the old three-strand component.
    # After growth it shifts to sigma_3 in the new output child {2,3,4}.
    # The prefix C_3(z) still crosses only that child with the outside z=5.
    # On the opposite side H_2 C_3(y) H_1 contains the adopted strand y;
    # it must remain literal for the next sweep, rather than becoming payload.
    prefix, payload = (2, 3, 4), (2,)
    opposite = (1, 3, 2, 1, 2)
    state = dict(
        _root_common(5),
        phase="layer",
        word=list(prefix + (-1, -2, -3) + payload)
        + [["slider", 3, 1, False, False, False]]
        + list(opposite),
        payload_length=1,
        completed_sides=0,
    )
    children, _ = _layer_step(state, Limits())
    growing, same_size = children
    control = growing["branch"]["control"]
    assert control["i"] == 4 and control["left_length"] == len(prefix)
    assert control["output_payload"] == [3]
    assert growing["branch"]["word"] == list(prefix) + [3] + [
        ["slider", 4, 1, False, False, False]
    ] + list(opposite)
    assert _expression_sum([state], 5, pending=True) == _expression_sum(
        children, 5, pending=True
    )

    # C_3(z) is consumed as a complete old-child packet, gaining a friend.
    # It is not mistaken for a whole-parent C_4(z) crossing or an internal run.
    following, event = _zip_step(growing, Catalog.bundled(), Limits(), True)
    assert event["rule_id"] == "zip.outer.positive.I"
    assert following[0]["branch"]["control"]["i"] == 4
    assert following[0]["branch"]["control"]["k"] == 1
    assert following[0]["branch"]["control"]["output_payload"] == [4]
    assert following[0]["branch"]["word"][-len(opposite):] == list(opposite)
    assert _expression_sum([growing], 5, pending=True) == _expression_sum(
        following, 5, pending=True
    )
    assert same_size["phase"] == "layer"
    continued, _ = _layer_step(same_size, Limits())
    assert all(child["phase"] == "layer" for child in continued)


@pytest.mark.parametrize("flags", [(True, False), (False, True), (True, True)])
@pytest.mark.parametrize("prefix", [(), (-1, -2, -3)])
def test_layering_rejects_unsupported_coupon_flags_without_clearing_them(flags, prefix):
    state = dict(
        _root_common(4),
        phase="layer",
        word=list(prefix) + [["slider", 3, 1, False, *flags]],
        payload_length=0,
        completed_sides=0,
    )
    original = deepcopy(state)
    with pytest.raises(ControlObligation, match="reversed or barred flags"):
        _layer_step(state, Limits())
    assert state == original


@pytest.mark.parametrize("size", [2, 3])
@pytest.mark.parametrize("flags", [(True, False), (False, True), (True, True)])
def test_extraction_rejects_unsupported_coupon_flags_before_normalizing(size, flags):
    state = dict(
        _root_common(size),
        phase="extract",
        word=[["slider", size, 1, False, *flags]],
    )
    original = deepcopy(state)
    with pytest.raises(ControlObligation, match="reversed or barred flags"):
        _extract_step(state, limits=Limits(), trace=False)
    assert state == original


@pytest.mark.parametrize("growth_size,fan", [(2, (1,)), (3, (1, 2, 1))])
def test_scheduler_extraction_rejects_collars_without_the_recorded_fan_bound(
    growth_size, fan
):
    # These are literal exterior permutation braids. They satisfy the public
    # normalizer's broader collar contract, but cannot be scheduler fans at
    # their stated growth sizes: the sole size-two fan has at most one leaf
    # crossing in a final three-strand component.
    state = dict(
        _root_common(3),
        phase="extract",
        boundary_left=list(fan),
        boundary_left_factors=[{"growth_size": growth_size, "word": list(fan)}],
        word=[["slider", 3, 1, False, False, False]],
    )
    original = deepcopy(state)
    with pytest.raises(ControlObligation, match="recorded growth-size bound"):
        _extract_step(state, limits=Limits(), trace=False)
    assert state == original


def test_endpoint_fans_are_saved_once_per_size_and_same_size_layers_do_not_rezip():
    result = reduce_word((1, 2, 3), n=4, trace=True)
    assert result["status"] == "complete", result["diagnostic"]
    seen_fan = False
    for transition in result["transitions"]:
        before = transition["predecessor"]
        for state in [before] + transition["successors"]:
            for side in ("left", "right"):
                factors = state[f"boundary_{side}_factors"]
                sizes = [factor["growth_size"] for factor in factors]
                assert len(sizes) == len(set(sizes))
                assert state[f"boundary_{side}"] == [
                    letter for factor in factors for letter in factor["word"]
                ]
                seen_fan |= bool(factors)
        if transition["evidence"].get("rule_id", "").startswith("layer.route."):
            growing, same_size = transition["successors"]
            assert growing["phase"] == "zip_up"
            assert same_size["phase"] == "layer"
            for side in ("left", "right"):
                assert same_size[f"boundary_{side}_factors"] == before[
                    f"boundary_{side}_factors"
                ]
    assert seen_fan


def test_trace_is_optional_and_word_inverse_pairs_are_preserved():
    session = AlgorithmSession()
    literal = [1, -1, 2, -2]
    compact = session.reduce(literal, n=3)
    traced = session.reduce(literal, n=3, trace=True)
    assert compact["status"] == traced["status"] == "complete"
    assert compact["input_word"] == traced["input_word"] == literal
    assert compact["output"] == traced["output"]
    assert compact["metrics"] == traced["metrics"]
    assert compact["transitions"] == [] and traced["transitions"]
    assert session.check(compact) and session.check(traced)


def test_limits_keep_exact_pending_expression_and_resume():
    stopped = reduce_word((1, 2), n=3, limits=Limits(max_steps=2), trace=True)
    assert stopped["status"] == "incomplete"
    assert stopped["pending"]
    assert check_certificate(json.loads(json.dumps(stopped)))
    following = resume(stopped, limits=Limits(max_steps=1000))
    assert following["status"] == "complete"
    assert following["trace_enabled"] is True
    assert following["output"] == reduce_word((1, 2), n=3)["output"]
    assert resume(stopped) == stopped
    with pytest.raises(AlgorithmError, match="verified prefix"):
        resume(following, limits=Limits(max_steps=1))


@pytest.mark.parametrize(
    "limits,diagnostic",
    [
        (Limits(max_terms=1), "term limit"),
        (Limits(max_word_size=1), "word size"),
        (Limits(max_word_size=3), "word size"),
        (Limits(max_coefficient_support=1), "coefficient support"),
    ],
)
def test_limits_reject_transition_without_losing_its_predecessor(limits, diagnostic):
    result = reduce_word((1, 2), n=3, limits=limits)
    assert result["status"] == "incomplete"
    assert diagnostic in result["diagnostic"]["message"]
    assert result["pending"]


def test_terminal_output_size_limit_is_checked_and_can_resume():
    # Odd negative root counts can expand when written as a positive route
    # followed by inverse full lassos. The predecessor fits the word limit,
    # whereas its exact canonical output can be longer.
    result = reduce_word((1, -2, 1), n=3, limits=Limits(max_word_size=6))
    assert result["status"] == "incomplete"
    assert "word size" in result["diagnostic"]["message"]
    assert result["pending"] and result["diagnostic"]["phase"] == "extract"
    assert all(len(term["word"]) <= 6 for term in result["output"])
    assert check_certificate(result)
    complete = resume(result, limits=Limits(max_word_size=100))
    assert complete["status"] == "complete", complete["diagnostic"]
    assert complete["metrics"]["max_local_word_tokens"] >= max(
        len(term["word"]) for term in complete["output"]
    )


@pytest.mark.parametrize(
    "changes",
    [
        {"n": True},
        {"n": 0},
        {"word": [0]},
        {"word": [3]},
        {"word": [True]},
        {"word": [1.0]},
        {"trace": 1},
        {"word": iter([1])},
    ],
)
def test_invalid_literal_inputs(changes):
    args = dict(word=[1], n=3)
    args.update(changes)
    with pytest.raises(AlgorithmError):
        reduce_word(**args)


def test_saved_certificate_tampering_and_session_isolation():
    session = AlgorithmSession()
    result = session.reduce((1, 2), n=3)
    original = deepcopy(result)
    result["output"][0]["coefficient"][0][1] += 1
    with pytest.raises(AlgorithmError):
        session.check(result)
    assert session.reduce((1, 2), n=3) == original
    with pytest.raises(AlgorithmError):
        check_certificate({})


def test_all_strands_in_a_longer_rank_are_scanned_without_isotopy_pruning():
    result = reduce_word((-3, 1, -1), n=4)
    assert result["status"] == "complete"
    assert result["input_word"] == [-3, 1, -1]
    assert result["metrics"]["rule_applications"]["single.far.positive.I"] > 0
    assert result["output"]


def test_singleton_and_empty_word_extraction():
    result = reduce_word((), n=1, trace=True)
    assert result["status"] == "complete"
    assert result["output"][0]["word"] == []
    assert result["output"][0]["components"][0]["record"]["twist_powers"] == []
    assert _common(result["transitions"][0]["predecessor"])["rank"] == 1


@pytest.mark.parametrize("mirrored", [False, True])
def test_two_strand_extraction_also_emits_the_default_normal_coupon(mirrored):
    state = dict(
        _root_common(2),
        phase="extract",
        word=[-1, -1, ["slider", 2, 1, mirrored, False, False], 1],
    )
    children, output, _ = _extract_step(state, limits=Limits(), trace=False)
    assert children == [] and len(output) == 1
    record = output[0]["components"][0]["record"]
    assert record["orientation"] == "normal" and record["twist_powers"] == [-1]
    assert output[0]["word"] == [-1, ["slider", 2, 1, False, False, False]]
    assert _expression_sum([state], 2, pending=True) == _expression_sum(
        output, 2, pending=False
    )


def test_integer_encoding_limit_is_a_resumable_refusal():
    result = reduce_word((1, 1), n=2, limits=Limits(max_integer_bits=1))
    assert result["status"] == "incomplete"
    assert "integer bit" in result["diagnostic"]["message"]
    assert check_certificate(result)


@pytest.mark.parametrize(
    "rank,word",
    [
        (3, (1, 2, 1)),
        (4, (1, -3, 2, 3)),
        (5, (4, 3, 2, 1, 4)),
    ],
)
def test_repeated_band_interactions_and_the_independent_mirrored_obstruction_complete(
    rank, word
):
    result = reduce_word(word, n=rank)
    assert result["status"] == "complete", result["diagnostic"]
    assert result["input_word"] == list(word)
    assert result["metrics"]["steps"] > len(word)
    assert all(
        c["record"].get("twist_powers") is not None
        for term in result["output"]
        for c in term["components"]
    )


@pytest.mark.parametrize("mirrored", [False, True])
@pytest.mark.parametrize("policy", ["M", "bar_mirror"])
@pytest.mark.parametrize("power", [-1, 1])
@pytest.mark.parametrize("friends", [0, 1])
def test_phase_rotation_canonicalizes_three_strand_half_twists_by_garside_identity(
    mirrored, policy, power, friends
):
    n, size, start = 5, 4, 2 if mirrored and friends else 1
    output_anchor = start - friends if mirrored else start + friends + 1
    state = framed.initialize(
        n,
        (),
        size=size,
        start=start,
        mirrored=mirrored,
        output_power=power,
        input_power=-power,
        output_payload=(output_anchor, -(output_anchor + 1)),
        mirrored_policy=policy,
        friends=friends,
    )
    work = _pass_work(_root_common(n), "zip_up", state.branches[0])
    following, evidence = _zip_step(work, Catalog.bundled(), Limits(), True)
    assert following[0]["phase"] == "zip_down"
    assert "Garside" in evidence["rotation_identity"]
    assert _expression_sum([work], n, pending=True) == _expression_sum(
        following, n, pending=True
    )


def test_actual_five_strand_staircase_terminal_gathering_and_refusal():
    # This exact branch occurs for the literal input (1, 2, 3, 4). Its child
    # decorations represent nested ribbons, although their raw letters are
    # interleaved by the recorded frame and rotation identities.
    upper = (1, -4, -3, -2, -1, -2, -3, -2, -4, -3, -2)
    lower = (1, 2, 3, 1, 2, 1, 1, 2, 3, 3, 2, 1, -3, 2, 3)
    state = dict(
        _root_common(5),
        phase="extract",
        world_phi=True,
        coefficient=[[-3, 1], [-2, -1]],
        depth=42,
        boundary_right=[-4, -3, -2],
        boundary_right_factors=[{"growth_size": 5, "word": [-4, -3, -2]}],
        word=list(upper) + [["slider", 5, 1, False, False, False]] + list(lower),
    )
    children, output, event = _extract_step(state, limits=Limits(), trace=True)
    assert children == [] and len(output) == 1
    assert output[0]["components"][0]["record"]["twist_powers"] == [0, 1, 1, 0]
    assert event["nested_slider"]["upper"]["gathering"] is not None
    with pytest.raises(RibbonLimit):
        _extract_step(state, limits=Limits(expansion_budget=1), trace=False)
