"""Independent regressions for bounded contextual replay."""
import json

import pytest

from zipper_algebra.zipping_algorithm.interpreter import (
    Branch,
    Limits,
    ReplayError,
    Request,
    RunState,
    check_certificate,
    execute,
    resume,
)


def _far_example():
    slider = ("slider", 1, 1, False, False, False)
    state = RunState(4, (Branch((3, slider)),))
    request = Request(0, 0, "single.far.positive.I",
                      {"n": 4, "i": 1, "j": 1, "k": 0, "l": 3})
    return state, request


def test_expansion_limit_saves_exact_predecessor_and_can_resume():
    state, request = _far_example()
    certificate = execute(state, [request], limits=Limits(expansion_budget=1))
    assert certificate["status"] == "incomplete"
    assert certificate["initial"] == certificate["final"] == state.to_json()
    assert certificate["transitions"] == []
    assert certificate["pending_requests"] == [request.to_json()]
    assert check_certificate(json.loads(json.dumps(certificate)))
    continued = resume(certificate, limits=Limits(expansion_budget=2))
    assert continued["status"] == "schedule_exhausted"
    assert continued["final"]["branches"][0]["word"] == [list(state.branches[0].word[1]), 3]
    assert check_certificate(json.loads(json.dumps(continued)))


def test_empty_schedule_checks_initial_limits():
    state = RunState(2, (Branch((1, 1)),))
    certificate = execute(state, [], limits=Limits(max_word_size=1))
    assert certificate["status"] == "incomplete"
    assert certificate["diagnostic"] == "word size limit exceeded"
    assert certificate["final"] == state.to_json()
    assert check_certificate(json.loads(json.dumps(certificate)))


def test_schedule_iterator_refused_without_consumption():
    state, request = _far_example()

    def unbounded_requests():
        raise AssertionError("execution must reject the iterator before consuming it")
        yield request

    with pytest.raises(ReplayError):
        execute(state, unbounded_requests(), limits=Limits(max_steps=1))
