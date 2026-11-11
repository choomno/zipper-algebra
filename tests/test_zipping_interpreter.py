"""Context retention, branch isolation, persisted replay and refusal regressions."""
from __future__ import annotations

import copy
import json

import pytest

from zipper_algebra.zipping_algorithm.interpreter import (
    Branch, Limits, LimitExceeded, ReplayError, Request, RunState,
    apply_rule, check_certificate, execute, metrics, resume,
)


def slider(i=1, j=1):
    return ("slider", i, j, False, False, False)


def far_request(**kwargs):
    return Request(0, 1, "single.far.positive.I", {"n": 4, "i": 1, "j": 1, "k": 0, "l": 3}, **kwargs)


def test_literal_context_retains_both_tails_and_inverses():
    state = RunState(4, (Branch((1, 3, slider(), -1)),))
    after, transition = apply_rule(state, far_request())
    assert after.branches[0].word == (1, slider(), 3, -1)
    assert transition["predecessor"] == state.branches[0].to_json()
    assert transition["successors"] == [after.branches[0].to_json()]
    # Adjacent inverses remain literal tokens throughout this layer.
    state = RunState(4, (Branch((1, -1, 3, slider())),))
    after, _ = apply_rule(state, Request(0, 2, "single.far.positive.I", {"n": 4, "i": 1, "j": 1, "k": 0, "l": 3}))
    assert after.branches[0].word[:2] == (1, -1)


def test_branching_coefficients_and_control_are_isolated():
    original = {"labels": [1, 2], "phase": "single"}
    state = RunState(2, (Branch((1, slider()), ((0, 2), (1, -1)), json.dumps(original)),))
    updates = [{"phase": "up"}, {"phase": "single"}, {"phase": "single"}]
    request = Request(0, 0, "single.foe.right.positive.I", {"n": 2, "i": 1, "j": 1, "k": 0, "l": 1},
                      evidence={"offensive": True}, successor_controls=updates)
    updates[0]["phase"] = "changed"
    after, transition = apply_rule(state, request)
    assert [b.word for b in after.branches] == [(slider(2),), (-1, slider()), (slider(),)]
    assert [b.coefficient for b in after.branches] == [((0, -2), (1, 1)),
                                                      ((1, 2), (2, -1)),
                                                      ((0, 2), (1, -3), (2, 1))]
    assert after.branches[0].control == {"phase": "up"}
    mutable = after.branches[1].control
    mutable["phase"] = "changed"
    assert after.branches[2].control == {"phase": "single"}
    assert state.branches[0].control == original
    assert transition["request"]["evidence"] == {"offensive": True}


def test_saved_reloaded_certificate_and_strict_tamper_rejection():
    state = RunState(4, (Branch((1, 3, slider(), -1)),))
    certificate = json.loads(json.dumps(execute(state, [far_request()])))
    assert check_certificate(certificate)
    assert certificate["status"] == "schedule_exhausted"
    assert certificate["scope"] == "algebraic_catalog_replay"
    mutations = [
        lambda c: c.update(extra=True),
        lambda c: c.update(schema="unknown"),
        lambda c: c.update(scope="normal_form"),
        lambda c: c.update(source_sha256="0" * 64),
        lambda c: c.update(catalog_sha256="0" * 64),
        lambda c: c["transitions"][0].update(rule_sha256="0" * 64),
        lambda c: c["transitions"][0]["successors"][0]["word"].append(2),
        lambda c: c["final"]["branches"][0]["coefficient"][0].__setitem__(1, 2),
        lambda c: c["metrics"][0].update(word_tokens=999),
        lambda c: c.update(status="complete"),
    ]
    for mutate in mutations:
        altered = copy.deepcopy(certificate)
        mutate(altered)
        with pytest.raises(ReplayError):
            check_certificate(altered)


def test_limits_save_predecessor_and_resume_pending_work():
    state = RunState(2, (Branch((1, slider())),))
    request = Request(0, 0, "single.foe.right.positive.I", {"n": 2, "i": 1, "j": 1, "k": 0, "l": 1})
    incomplete = execute(state, [request], limits=Limits(max_terms=2))
    assert incomplete["status"] == "incomplete"
    assert incomplete["final"] == state.to_json()
    assert len(incomplete["pending_requests"]) == 1
    assert incomplete["transitions"] == []
    assert check_certificate(json.loads(json.dumps(incomplete)))
    continued = resume(incomplete, limits=Limits(max_terms=3))
    assert continued["status"] == "schedule_exhausted"
    assert len(continued["final"]["branches"]) == 3
    assert check_certificate(continued)


def test_step_limit_is_atomic():
    state = RunState(4, (Branch((3, 3, slider())),))
    requests = [Request(0, 1, "single.far.positive.I", {"n": 4, "i": 1, "j": 1, "k": 0, "l": 3}),
                Request(0, 0, "single.far.positive.I", {"n": 4, "i": 1, "j": 1, "k": 0, "l": 3})]
    certificate = execute(state, requests, limits=Limits(max_steps=1))
    assert certificate["status"] == "incomplete"
    assert certificate["diagnostic"] == "step limit exceeded"
    assert len(certificate["transitions"]) == 1
    assert check_certificate(certificate)
    final = resume(certificate)
    assert final["final"]["branches"][0]["word"] == [list(slider()), 3, 3]


def test_coefficient_and_word_limits():
    states = [
        (RunState(2, (Branch((1, 1)),)), Limits(max_word_size=1)),
        (RunState(2, (Branch((), ((0, 1), (1, 1))),)), Limits(max_coefficient_support=1)),
        (RunState(2, (Branch((), ((0, 256),)),)), Limits(max_integer_bits=8)),
    ]
    for state, limits in states:
        with pytest.raises(LimitExceeded):
            limits.check(state)


def test_exact_collection_metrics_include_cancellation():
    state = RunState(2, (Branch((1,), ((-2, 9), (1, 1))),
                         Branch((1,), ((-2, -9), (1, -1))), Branch(())))
    assert metrics(state) == {"uncollected_terms": 3, "collected_terms": 1,
                              "word_tokens": 2, "max_word_tokens": 1,
                              "coefficient_support": 5, "q_exponent_min": -2,
                              "q_exponent_max": 1, "max_integer_bits": 4}
    assert metrics(RunState(1, ())) == {"uncollected_terms": 0, "collected_terms": 0,
                                       "word_tokens": 0, "max_word_tokens": 0,
                                       "coefficient_support": 0, "q_exponent_min": 0,
                                       "q_exponent_max": 0, "max_integer_bits": 0}


@pytest.mark.parametrize("word", [(0,), (True,), (4,), (("slider", 2, 4, False, False, False),),
                                  (("slider", 2, 1, 0, False, False),), (("bad",),)])
def test_invalid_tokens_rejected(word):
    with pytest.raises(ReplayError):
        RunState(4, (Branch(word),))


def test_invalid_states_coefficients_and_requests():
    for coefficient in [(), ((0, 0),), ((0, True),), ((1, 1), (0, 2)), ((0, 1), (0, 2))]:
        with pytest.raises(ReplayError):
            Branch((), coefficient)
    with pytest.raises(ReplayError):
        Branch((), control_json="[]")
    for n in [0, True]:
        with pytest.raises(ReplayError):
            RunState(n, ())
    with pytest.raises(ReplayError):
        RunState(2, ("not a branch",))
    for args in [(-1, 0, "rule", {}), (0, True, "rule", {}), (0, 0, "", {}),
                 (0, 0, "rule", {"n": True})]:
        with pytest.raises(ReplayError):
            Request(*args)
    with pytest.raises(ReplayError):
        Request(0, 0, "rule", {}, words={"P": [True]})
    with pytest.raises(ReplayError):
        Limits(max_steps=0)
    for method, record in [(Branch.from_json, {"word": []}),
                           (RunState.from_json, {"n": 2}),
                           (Request.from_json, {"branch": 0})]:
        with pytest.raises(ReplayError):
            method(record)


def test_missing_occurrence_rank_guard_and_control_updates_rejected():
    state = RunState(4, (Branch((1, 3, slider(), -1)),))
    invalid = [Request(1, 1, "single.far.positive.I", {"n": 4, "i": 1, "j": 1, "k": 0, "l": 3}),
               Request(0, 0, "single.far.positive.I", {"n": 4, "i": 1, "j": 1, "k": 0, "l": 3}),
               Request(0, 1, "single.far.positive.I", {"n": 3, "i": 1, "j": 1, "k": 0, "l": 3}),
               Request(0, 1, "single.far.positive.I", {"n": 4, "i": 1, "j": 1, "k": 0, "l": 1}),
               far_request(successor_controls=[])]
    for request in invalid:
        with pytest.raises(ReplayError):
            apply_rule(state, request)


def test_empty_schedule_and_roundtrips():
    branch = Branch((slider(),), ((-1, 2),), '{"a":1}')
    assert Branch.from_json(branch.to_json()) == branch
    state = RunState(1, (branch,))
    assert RunState.from_json(state.to_json()) == state
    request = Request(0, 0, "unused", {"n": 1}, words={"P": []})
    assert Request.from_json(request.to_json()).to_json() == request.to_json()
    certificate = execute(state, [])
    assert check_certificate(certificate)


def test_external_catalog_rejects_unaccepted_or_noncontextual_rules():
    from zipper_algebra.zipping_algorithm.catalog import Catalog
    bundled = Catalog.bundled()
    source = bundled.get("single.far.positive.I")
    variants = []
    unaccepted = copy.deepcopy(source)
    unaccepted["execution"]["algebra_verified"] = False
    variants.append(unaccepted)
    polynomial_lhs = copy.deepcopy(source)
    polynomial_lhs["lhs"].append({"coefficient": [{"integer": 1, "q_power": 0}],
                                    "word": {"type": "identity"}})
    variants.append(polynomial_lhs)
    scalar_lhs = copy.deepcopy(source)
    scalar_lhs["lhs"][0]["coefficient"][0]["integer"] = 2
    variants.append(scalar_lhs)
    for rule in variants:
        catalog = Catalog({"schema_version": "1.0", "rules": [rule]})
        with pytest.raises(ReplayError):
            apply_rule(RunState(4, (Branch((1, 3, slider(), -1)),)), far_request(), catalog=catalog)


def test_initial_refusal_retains_nonempty_schedule():
    state = RunState(4, (Branch((1, 3, slider(), -1)),))
    certificate = execute(state, [far_request()], limits=Limits(max_word_size=3))
    assert certificate["status"] == "incomplete"
    assert certificate["final"] == state.to_json()
    assert certificate["pending_requests"] == [far_request().to_json()]
    assert check_certificate(certificate)


def test_external_context_with_missing_named_word_rejected():
    from zipper_algebra.zipping_algorithm.catalog import Catalog
    rule = Catalog.bundled().get("single.far.positive.I")
    rule["rhs"][0]["word"] = {"type": "word_ref", "name": "P"}
    external = Catalog({"schema_version": "1.0", "rules": [rule]})
    with pytest.raises(ReplayError, match="Unbound word reference"):
        apply_rule(RunState(4, (Branch((1, 3, slider(), -1)),)), far_request(), catalog=external)


def test_optional_transition_trace_preserves_exact_result():
    state = RunState(4, (Branch((1, 3, slider(), -1)),))
    full, event = apply_rule(state, far_request())
    compact, omitted = apply_rule(state, far_request(), trace=False)
    assert compact == full
    assert event["request"]["rule_id"] == "single.far.positive.I"
    assert omitted == {}
    with pytest.raises(ReplayError, match="trace must"):
        apply_rule(state, far_request(), trace=1)
