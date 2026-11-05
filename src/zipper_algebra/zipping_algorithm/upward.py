"""Deterministic normal upward passes with every algebraic branch retained.

The input represents U F_(i,j,k) X_(i,j) P T in algebraic product order.
This module completes a single upward pass when the normal-slider rules
suffice. A higher flype retains its two half twists and stops with an explicit
frame obligation. Completion of a pass does not assert completion of zipping.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from .catalog import Catalog, CatalogError, condition, expand_polynomial, expand_word, integer
from .interpreter import (
    Branch, Limits, LimitExceeded, ReplayError, Request, RunState, apply_rule, metrics,
)


class ControlError(ReplayError):
    """A missing, ambiguous, or unsupported literal upward context."""


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _friends(i: int, j: int, k: int) -> tuple[int, ...]:
    return tuple(x for a in range(k, 0, -1) for x in range(j + a, j + a + i - 1))


def _control(i: int, j: int, k: int, left: int, processed: int, tail: int) -> dict[str, Any]:
    return {"phase": "zip_up" if left else "upper_pass_complete", "i": i, "j": j, "k": k,
            "left_length": left, "processed_length": processed, "tail_length": tail}


def initialize(n: int, left: Iterable[int], *, size: int = 2, start: int = 1,
               friends: int = 0, processed: Iterable[int] = (),
               tail: Iterable[int] = ()) -> RunState:
    """Initialize an explicitly supplied normal upward state.

    The complete processed word and exterior tail remain separate boundaries.
    Admissibility of this state as part of a full zipper run is a separate
    invariant; this entry point accepts the stated local context as its input.
    """
    if (any(type(x) is not int for x in (n, size, start, friends))
            or not 2 <= size <= n or start < 1 or friends < 0
            or size + start + friends - 1 > n):
        raise ControlError("valid upward slider and friend intervals required")
    u, p, t = tuple(left), tuple(processed), tuple(tail)
    if any(type(x) is not int for x in u + p + t):
        raise ControlError("exterior words must contain signed braid integers")
    word = u + _friends(size, start, friends) + (("slider", size, start, False, False, False),) + p + t
    control = _control(size, start, friends, len(u), len(p), len(t))
    return RunState(n, (Branch(word, control_json=_json(control)),))


def _parts(branch: Branch, n: int) -> tuple[dict[str, Any], tuple[Any, ...], tuple[Any, ...], tuple[Any, ...]]:
    c = branch.control
    if set(c) != {"phase", "i", "j", "k", "left_length", "processed_length", "tail_length"}:
        raise ControlError("upward state boundaries are missing or unexpected")
    if (c["phase"] not in ("zip_up", "upper_pass_complete")
            or any(type(c[x]) is not int for x in ("i", "j", "k", "left_length", "processed_length", "tail_length"))
            or not 2 <= c["i"] <= n or c["j"] < 1 or c["k"] < 0
            or c["i"] + c["j"] + c["k"] - 1 > n
            or min(c["left_length"], c["processed_length"], c["tail_length"]) < 0):
        raise ControlError("invalid upward state boundaries")
    a, b, d = c["left_length"], c["processed_length"], c["tail_length"]
    core = _friends(c["i"], c["j"], c["k"]) + (("slider", c["i"], c["j"], False, False, False),)
    if (len(branch.word) != a + len(core) + b + d
            or branch.word[a:a + len(core)] != core
            or c["phase"] != ("zip_up" if a else "upper_pass_complete")):
        raise ControlError("upward boundaries do not describe the full expression")
    u = branch.word[:a]
    p = branch.word[a + len(core):a + len(core) + b]
    t = branch.word[a + len(core) + b:]
    if any(type(x) is not int for x in u + p + t):
        raise ControlError("upward exterior words must be pure braid words")
    return c, u, p, t


def _validate_state(state: RunState) -> None:
    for branch in state.branches:
        if branch.control.get("phase") == "unsupported_frame":
            if set(branch.control) != {"phase", "rule_id", "message"}:
                raise ControlError("invalid retained flype frame record")
        else:
            _parts(branch, state.n)


def _flatten(node: dict[str, Any]) -> list[dict[str, Any]]:
    if node["type"] == "concat":
        return [x for child in node["items"] for x in _flatten(child)]
    return [node]


def _literal(nodes: list[dict[str, Any]], env: dict[str, int], budget: int) -> tuple[Any, ...]:
    poly = expand_word({"type": "concat", "items": nodes}, env, max_terms=budget)
    return next(iter(poly))[1]


def _successor_control(term: dict[str, Any], env: dict[str, int], remaining: int,
                       processed: int, tail: int, budget: int, rule_id: str) -> dict[str, Any]:
    nodes = _flatten(term["word"])
    at = next(idx for idx, x in enumerate(nodes) if x["type"] == "slider")
    slider = nodes[at]
    i, j = integer(slider["size"], env), integer(slider["start"], env)
    if slider.get("orientation", "normal") == "mirrored" and i > 2:
        return {"phase": "unsupported_frame", "rule_id": rule_id,
                "message": "higher mirrored slider requires transport of retained half twists"}
    k, before = 0, at
    if at and nodes[at - 1]["type"] == "block":
        block = nodes[at - 1]
        if (integer(block["length"], env) == i - 1 and integer(block["start"], env) == j + 1
                and integer(block.get("sign", 1), env) == 1 and block["order"] == "descending"):
            k, before = integer(block["count"], env), at - 1
    emitted_left = len(_literal(nodes[:before], env, budget))
    emitted_processed = len(_literal(nodes[at + 1:], env, budget))
    return _control(i, j, k, remaining + emitted_left, processed + emitted_processed, tail)


def step(state: RunState, branch: int = 0, *, catalog: Catalog | None = None,
         limits: Limits | None = None, trace: bool = True) -> tuple[RunState, dict[str, Any]]:
    """Perform the unique matching normal rule, requiring its entire tape.

    For i=2, the internal commuting rule takes precedence over a flype of the
    same internal crossing. Its mirrored slider equals the normal slider by
    the defining polynomial q sigma_j^-1 + 1 - q - sigma_j.
    """
    if type(trace) is not bool:
        raise ControlError("trace must be a boolean")
    try:
        return _step(state, branch, catalog=catalog, limits=limits, trace=trace)
    except CatalogError as exc:
        if "budget" in str(exc):
            raise LimitExceeded(str(exc)) from exc
        raise


def _step(state: RunState, branch: int, *, catalog: Catalog | None,
          limits: Limits | None, trace: bool) -> tuple[RunState, dict[str, Any]]:
    catalog = Catalog.bundled() if catalog is None else catalog
    limits = Limits() if limits is None else limits
    if type(branch) is not int or not 0 <= branch < len(state.branches):
        raise ControlError("selected upward branch does not exist")
    predecessor = state.branches[branch]
    c, u, p, t = _parts(predecessor, state.n)
    if not u:
        raise ControlError("upward pass is already complete")
    env = {"n": state.n, "i": c["i"], "j": c["j"], "k": c["k"],
           "l": abs(u[-1]), "eps": 1 if u[-1] > 0 else -1}
    prefix_end = len(u) + (c["i"] - 1) * c["k"] + 1
    candidates = []
    for rule in catalog.rules:
        rid = rule["id"]
        if (not rid.startswith("zip.") or not rid.endswith(".I")
                or (c["i"] == 2 and rid.startswith("zip.flype."))):
            continue
        if not condition(rule["condition"], env):
            continue
        lhs = expand_polynomial(rule["lhs"], env, max_terms=limits.expansion_budget)
        occurrence = next(iter(lhs))[1]
        offset = prefix_end - len(occurrence)
        if 0 <= offset < len(u) and predecessor.word[offset:prefix_end] == occurrence:
            candidates.append((rule, offset))
    if len(candidates) != 1:
        raise ControlError("missing full-tape rule" if not candidates else
                           "competing upward rules: " + ", ".join(r["id"] for r, _ in candidates))
    rule, offset = candidates[0]
    # Ordered word collection agrees with the contextual interpreter. A word
    # appears once in the controller even when a coefficient is serialized in
    # multiple pieces; its complete literal exterior is retained.
    terms_by_word = {_literal([term["word"]], env, limits.expansion_budget): term
                     for term in rule["rhs"]}
    controls = [_successor_control(term, env, offset, len(p), len(t),
                                   limits.expansion_budget, rule["id"])
                for term in terms_by_word.values()]
    request = Request(branch, offset, rule["id"], env,
                      evidence={"controller": "normal_upward", "full_tape_literal_match": True,
                                "i2_internal_precedes_flype": c["i"] == 2},
                      successor_controls=controls)
    raw, equality = apply_rule(state, request, catalog=catalog, limits=limits, trace=trace)
    normalizations = []
    branches = list(raw.branches)
    for idx in range(branch, branch + len(controls)):
        old = branches[idx]
        word = tuple(("slider", 2, x[2], False, False, False)
                     if isinstance(x, tuple) and x[0] == "slider" and x[1] == 2 and x[3]
                     else x for x in old.word)
        if word != old.word:
            branches[idx] = Branch(word, old.coefficient, old.control_json)
            normalizations.append({"branch": idx, "identity": "mirrored_X2_equals_X2"})
    following = RunState(state.n, tuple(branches))
    successors = following.branches[branch:branch + len(controls)]
    for successor in successors:
        if successor.control["phase"] != "unsupported_frame":
            _parts(successor, state.n)
    powers = [q for successor in successors for q, _ in successor.coefficient]
    delta = {"word_tokens": sum(len(b.word) for b in successors) - len(predecessor.word),
             "coefficient_support": sum(len(b.coefficient) for b in successors) - len(predecessor.coefficient),
             "max_word_tokens": max(len(b.word) for b in successors),
             "q_exponent_min": min(powers), "q_exponent_max": max(powers),
             "max_integer_bits": max(abs(c).bit_length() for b in successors for _, c in b.coefficient)}
    event = {"rule_id": rule["id"], "metric_delta": delta}
    if trace:
        event.update(catalog_transition=equality, normalizations=normalizations,
                     successor_state=following.to_json())
    return following, event


def run(state: RunState, *, catalog: Catalog | None = None,
        limits: Limits | None = None, trace: bool = False) -> dict[str, Any]:
    """Run upward passes with an optional detailed execution trace.

    Compact mode retains the exact input, final checkpoint and streaming size
    statistics. Detailed mode additionally serializes every transition and
    every state metric. Both modes use the same rule choices and arithmetic.
    """
    if type(trace) is not bool:
        raise ControlError("trace must be a boolean")
    catalog = Catalog.bundled() if catalog is None else catalog
    limits = Limits() if limits is None else limits
    initial = state
    _validate_state(state)
    transitions: list[dict[str, Any]] = []
    first_metrics = metrics(state)
    history = [first_metrics]
    steps = 0
    counts: Counter[str] = Counter()
    summary = {"peak_uncollected_terms": first_metrics["uncollected_terms"],
               "peak_word_tokens": first_metrics["word_tokens"],
               "peak_coefficient_support": first_metrics["coefficient_support"],
               "max_word_tokens": first_metrics["max_word_tokens"],
               "q_exponent_min": first_metrics["q_exponent_min"],
               "q_exponent_max": first_metrics["q_exponent_max"],
               "max_integer_bits": first_metrics["max_integer_bits"]}
    word_tokens, support = first_metrics["word_tokens"], first_metrics["coefficient_support"]
    diagnostic = None
    while True:
        selected = next((idx for idx, b in enumerate(state.branches) if b.control.get("phase") == "zip_up"), None)
        try:
            limits.check(state)
        except LimitExceeded as exc:
            status, diagnostic = "incomplete", {"branch": selected, "message": str(exc), "state": state.to_json()}
            break
        if selected is None:
            status = ("unsupported_frame" if any(b.control.get("phase") == "unsupported_frame"
                                                for b in state.branches) else "upper_pass_complete")
            break
        try:
            if steps >= limits.max_steps:
                raise LimitExceeded("step limit exceeded")
            following, event = step(state, selected, catalog=catalog, limits=limits, trace=trace)
        except (ControlError, LimitExceeded) as exc:
            status = "incomplete" if isinstance(exc, LimitExceeded) else "missing_rule"
            diagnostic = {"branch": selected, "message": str(exc), "state": state.branches[selected].to_json()}
            break
        state = following
        steps += 1
        counts[event["rule_id"]] += 1
        delta = event["metric_delta"]
        word_tokens += delta["word_tokens"]
        support += delta["coefficient_support"]
        summary["peak_uncollected_terms"] = max(summary["peak_uncollected_terms"], len(state.branches))
        summary["peak_word_tokens"] = max(summary["peak_word_tokens"], word_tokens)
        summary["peak_coefficient_support"] = max(summary["peak_coefficient_support"], support)
        for name in ("max_word_tokens", "q_exponent_max", "max_integer_bits"):
            summary[name] = max(summary[name], delta[name])
        summary["q_exponent_min"] = min(summary["q_exponent_min"], delta["q_exponent_min"])
        if trace:
            transitions.append(event)
            history.append(metrics(state))
    if not trace:
        history.append(metrics(state))
    return {"schema": "zipper-normal-upward-v1", "scope": "normal_upward_pass_only",
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "dependency_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                  for name in ("catalog.py", "interpreter.py")},
            "catalog_sha256": catalog.sha256, "initial": initial.to_json(), "final": state.to_json(),
            "status": status, "diagnostic": diagnostic, "transitions": transitions,
            "limits": asdict(limits), "trace_enabled": trace, "metrics": history,
            "summary": dict(summary, steps=steps, rule_counts=dict(counts))}


def check_certificate(certificate: Mapping[str, Any], *, catalog: Catalog | None = None) -> bool:
    """Reload and replay all rule choices, boundary updates and coefficients."""
    expected = run(RunState.from_json(certificate["initial"]), catalog=catalog,
                   limits=Limits(**certificate["limits"]), trace=certificate["trace_enabled"])
    if _json(expected) != _json(certificate):
        raise ControlError("upward certificate differs from deterministic exact replay")
    return True


def resume(certificate: Mapping[str, Any], *, catalog: Catalog | None = None,
           limits: Limits | None = None, trace: bool | None = None) -> dict[str, Any]:
    """Verify and continue a checkpoint, preserving its trace mode by default."""
    check_certificate(certificate, catalog=catalog)
    return run(RunState.from_json(certificate["final"]), catalog=catalog, limits=limits,
               trace=certificate["trace_enabled"] if trace is None else trace)
