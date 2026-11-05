"""Upward passes retaining both flype frames in exact algebraic order.

A state denotes U D Delta_out^a F X Delta_in^b P T. The protected braid D
is supported on the output tape. The two frame powers
are recorded as integers; its actual word stores every half-twist crossing.
Normal and mirrored catalog rules act only after a checked marked-tape
transport. This controller establishes completion of upward passes in its
supported marked-tape contexts, independently of later layering/extraction.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from . import upward
from .catalog import Catalog, CatalogError, condition, expand_polynomial, integer
from .frames import FrameError, check_frame, transport_frame
from .interpreter import Branch, Limits, LimitExceeded, Request, RunState, apply_rule, metrics
from .semilinear import semilinear_zip_rule

ControlError = upward.ControlError


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _anchors(i: int, j: int, k: int, mirrored: bool) -> tuple[int, int]:
    return (j - k, j + 1) if mirrored else (j + k + 1, j)


def _core(n: int, i: int, j: int, k: int, mirrored: bool, policy: str = "M") -> tuple[Any, ...]:
    base_j = n - i - j + 2 if mirrored else j
    friends = upward._friends(i, base_j, k)
    if mirrored:
        friends = tuple(n - x for x in friends)
        if policy == "bar_mirror":
            friends = tuple(-x for x in friends)
    return friends + (("slider", i, j, mirrored, False, False),)


def friend_word(n: int, control: Mapping[str, Any]) -> tuple[int, ...]:
    """The actual signed friend braid in the selected mirrored policy."""
    return _core(n, control["i"], control["j"], control["k"], control["mirrored"],
                 control["mirror_policy"])[:-1]


def rule_for_event(catalog: Catalog, event: Mapping[str, Any]) -> dict[str, Any]:
    """Recover an original or semilinearly derived local identity."""
    rid = event["rule_id"]
    if rid.endswith(".BM"):
        return semilinear_zip_rule(catalog.get(rid[:-3] + ".M"), dict(event["parameters"]))
    return catalog.get(rid)


class _ExecutionCatalog(Catalog):
    """Read-only derived-rule view; the frozen catalog remains unchanged."""

    def __init__(self, base: Catalog, env: dict[str, int]):
        self.base, self.env, self.sha256 = base, env, base.sha256

    def get(self, rule_id: str) -> dict[str, Any]:
        return rule_for_event(self.base, {"rule_id": rule_id, "parameters": self.env})


def _power(size: int, start: int, exponent: int, budget: int) -> tuple[int, ...]:
    if size * (size - 1) // 2 * abs(exponent) > budget:
        raise LimitExceeded("expanded frame word limit exceeded")
    half = tuple(x for length in range(size - 1, 0, -1) for x in range(start, start + length))
    if exponent < 0:
        half = tuple(-x for x in reversed(half))
    return half * abs(exponent)


def _control(i: int, j: int, k: int, mirrored: bool, a: int, b: int,
             left: int, processed: int, tail: int) -> dict[str, Any]:
    return {"phase": "framed_zip_up" if left else "framed_upper_pass_complete",
            "i": i, "j": j, "k": k, "mirrored": mirrored,
            "output_power": a, "input_power": b, "left_length": left,
            "processed_length": processed, "tail_length": tail, "output_payload": []}


def _expression(n: int, c: dict[str, Any], u: tuple[Any, ...], p: tuple[Any, ...],
                tail: tuple[Any, ...], budget: int) -> tuple[Any, ...]:
    i, j, k = c["i"], c["j"], c["k"]
    payload = tuple(c["output_payload"])
    total = (len(u) + len(payload) + len(p) + len(tail) + (i - 1) * k + 1
             + (i - 1) * (i - 2) // 2 * (abs(c["output_power"]) + abs(c["input_power"])))
    if total > budget:
        raise LimitExceeded("expanded framed state word limit exceeded")
    out, inp = _anchors(i, j, k, c["mirrored"])
    return (u + payload + _power(i - 1, out, c["output_power"], budget)
            + _core(n, i, j, k, c["mirrored"], c["mirror_policy"])
            + _power(i - 1, inp, c["input_power"], budget) + p + tail)


def initialize(n: int, left: Iterable[int], *, size: int = 3, start: int = 1,
               friends: int = 0, mirrored: bool = False, output_power: int = 0,
               input_power: int = 0, processed: Iterable[int] = (),
               tail: Iterable[int] = (), output_payload: Iterable[int] = (),
               mirrored_policy: str = "M",
               limits: Limits | None = None) -> RunState:
    """Initialize a stated marked-tape context, including its entire exterior.

    ``bar_mirror`` derives mirrored rules by simultaneously inverting crossing
    signs and q; their common same-rank slider unit cancels. Its mirrored
    friends have negative crossings. ``M`` retains the reflection reference.
    """
    limits = Limits() if limits is None else limits
    if (any(type(x) is not int for x in (n, size, start, friends, output_power, input_power))
            or type(mirrored) is not bool or mirrored_policy not in ("M", "bar_mirror")
            or not 2 <= size <= n or start < 1
            or friends < 0 or start + size - 1 > n
            or (start - friends < 1 if mirrored else start + size + friends - 1 > n)):
        raise ControlError("valid framed slider, friends and integer powers required")
    if size == 2:
        if mirrored or tuple(output_payload):
            raise ControlError("use the normal X2 convention for the two-strand delegate")
        return upward.initialize(n, left, size=size, start=start, friends=friends,
                                 processed=processed, tail=tail)
    u, p, t = tuple(left), tuple(processed), tuple(tail)
    if any(type(x) is not int for x in u + p + t):
        raise ControlError("framed exterior must contain signed braid integers")
    c = _control(size, start, friends, mirrored, output_power, input_power, len(u), len(p), len(t))
    c["mirror_policy"] = mirrored_policy
    payload = tuple(output_payload)
    out, _ = _anchors(size, start, friends, mirrored)
    if any(type(x) is not int or not out <= abs(x) < out + size - 2 for x in payload):
        raise ControlError("protected output payload must be supported on the marked tape")
    c["output_payload"] = list(payload)
    c["initial_slider"] = {"j": start, "k": friends, "mirrored": mirrored,
                           "output_power": output_power, "input_power": input_power,
                           "processed_length": len(p), "mirror_policy": mirrored_policy}
    c["processed_events"] = []
    state = RunState(n, (Branch(_expression(n, c, u, p, t, limits.max_word_size), control_json=_json(c)),))
    limits.check(state)
    return state


def _parts(branch: Branch, n: int, budget: int) -> tuple[dict[str, Any], tuple[Any, ...], tuple[Any, ...], tuple[Any, ...]]:
    c = branch.control
    keys = {"phase", "i", "j", "k", "mirrored", "output_power", "input_power",
            "left_length", "processed_length", "tail_length", "initial_slider", "processed_events", "output_payload", "mirror_policy"}
    if (set(c) != keys or type(c["mirrored"]) is not bool
            or c["mirror_policy"] not in ("M", "bar_mirror")
            or any(type(c[x]) is not int for x in keys - {"phase", "mirrored", "initial_slider", "processed_events", "output_payload", "mirror_policy"})
            or not 3 <= c["i"] <= n or c["j"] < 1 or c["k"] < 0
            or c["j"] + c["i"] - 1 > n
            or (c["j"] - c["k"] < 1 if c["mirrored"] else c["j"] + c["i"] + c["k"] - 1 > n)
            or min(c["left_length"], c["processed_length"], c["tail_length"]) < 0):
        raise ControlError("invalid framed state boundaries")
    out, _ = _anchors(c["i"], c["j"], c["k"], c["mirrored"])
    if (not isinstance(c["output_payload"], list)
            or any(type(x) is not int or not out <= abs(x) < out + c["i"] - 2 for x in c["output_payload"])):
        raise ControlError("invalid protected output-tape payload")
    initial = c["initial_slider"]
    if (not isinstance(initial, dict)
            or set(initial) != {"j", "k", "mirrored", "output_power", "input_power", "processed_length", "mirror_policy"}
            or type(initial["mirrored"]) is not bool
            or initial["mirror_policy"] != c["mirror_policy"]
            or any(type(initial[x]) is not int for x in initial if x not in ("mirrored", "mirror_policy"))
            or not isinstance(c["processed_events"], list)):
        raise ControlError("invalid framed input-band record")
    for event in c["processed_events"]:
        if (not isinstance(event, dict) or set(event) != {"rule_id", "parameters", "rhs_index", "mirror_policy"}
                or not isinstance(event["rule_id"], str) or not isinstance(event["parameters"], dict)
                or any(type(x) is not int for x in event["parameters"].values())
                or event["mirror_policy"] != c["mirror_policy"]
                or type(event["rhs_index"]) is not int or event["rhs_index"] < 0):
            raise ControlError("invalid processed-band event")
    a, b, d = c["left_length"], c["processed_length"], c["tail_length"]
    if c["phase"] != ("framed_zip_up" if a else "framed_upper_pass_complete"):
        raise ControlError("invalid framed phase")
    p_start = len(branch.word) - b - d
    u, p, t = branch.word[:a], branch.word[p_start:p_start + b], branch.word[p_start + b:]
    if any(type(x) is not int for x in u + p + t):
        raise ControlError("framed exterior must be a pure braid word")
    if branch.word != _expression(n, c, u, p, t, budget):
        raise ControlError("framed boundaries do not describe the complete expression")
    return c, u, p, t


def _branch_shape(base_term: dict[str, Any], env: dict[str, int], mirrored: bool) -> tuple[int, int, bool, int, int]:
    nodes = upward._flatten(base_term["word"])
    at = next(idx for idx, x in enumerate(nodes) if x["type"] == "slider")
    slider = nodes[at]
    j = integer(slider["start"], env)
    base_mirror = slider.get("orientation", "normal") == "mirrored"
    k = 0
    if at and nodes[at - 1]["type"] == "block":
        k = integer(nodes[at - 1]["count"], env)

    def exponent(node):
        if node["type"] == "half_twist":
            return 1
        if node["type"] == "inverse" and node["word"]["type"] == "half_twist":
            return -1
        return 0

    aa = exponent(nodes[0]) if at else 0
    bb = exponent(nodes[at + 1]) if at + 1 < len(nodes) else 0
    return (env["n"] - env["i"] - j + 2 if mirrored else j,
            k, base_mirror != mirrored, aa, bb)


def step(state: RunState, branch: int = 0, *, catalog: Catalog | None = None,
         limits: Limits | None = None, trace: bool = True) -> tuple[RunState, dict[str, Any]]:
    """Transport the leading frame, apply one rule, and transport its tail."""
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
    limits.check(state)
    if type(branch) is not int or not 0 <= branch < len(state.branches):
        raise ControlError("selected framed branch does not exist")
    old = state.branches[branch]
    c, u, p, tail = _parts(old, state.n, limits.max_word_size)
    if not u:
        raise ControlError("framed upward pass is already complete")
    i, j, k, mirror = c["i"], c["j"], c["k"], c["mirrored"]
    out, inp = _anchors(i, j, k, mirror)
    policy = c["mirror_policy"]
    bar_mirror = mirror and policy == "bar_mirror"
    if out <= abs(u[-1]) < out + i - 2:
        # This changes only the boundary between U and the protected child
        # braid. Every crossing stays in exactly the same literal position.
        cc = dict(c)
        cc["left_length"] -= 1
        cc["phase"] = "framed_zip_up" if cc["left_length"] else "framed_upper_pass_complete"
        cc["output_payload"] = [u[-1]] + c["output_payload"]
        new = Branch(old.word, old.coefficient, _json(cc))
        following = RunState(state.n, state.branches[:branch] + (new,) + state.branches[branch + 1:])
        _parts(new, state.n, limits.max_word_size)
        powers = [q for q, _ in new.coefficient]
        event: dict[str, Any] = {"rule_id": "control.absorb_internal_tape", "metric_delta": {
            "word_tokens": 0, "coefficient_support": 0, "max_word_tokens": len(new.word),
            "q_exponent_min": min(powers), "q_exponent_max": max(powers),
            "max_integer_bits": max(abs(v).bit_length() for _, v in new.coefficient)}}
        if trace:
            event.update(equality_reference="identical_literal_word", successor_state=following.to_json())
        return following, event
    env = {"n": state.n, "i": i, "j": state.n - i - j + 2 if mirror else j,
           "k": k, "l": state.n - abs(u[-1]) if mirror else abs(u[-1]),
           "eps": (1 if u[-1] > 0 else -1) * (-1 if bar_mirror else 1)}
    core = _core(state.n, i, j, k, mirror, policy)
    variant = ".M" if mirror else ".I"
    candidates = []
    for rule in catalog.rules:
        if not rule["id"].startswith("zip.") or not rule["id"].endswith(variant):
            continue
        if not condition(rule["condition"], env):
            continue
        if bar_mirror:
            rule = semilinear_zip_rule(rule, env)
        lhs = expand_polynomial(rule["lhs"], env, max_terms=limits.expansion_budget)
        occurrence = next(iter(lhs))[1]
        if occurrence[-len(core):] != core:
            raise ControlError("catalog left side does not retain the framed core")
        w = occurrence[:-len(core)]
        if not w or len(w) > len(u):
            continue
        try:
            witness = transport_frame(u[-len(w):], n=state.n, size=i - 1, start=out,
                                      power=c["output_power"], side="left")
            check_frame(witness)
        except FrameError:
            continue
        if witness.transformed == w:
            candidates.append((rule, w, witness))
    if len(candidates) != 1:
        raise ControlError("missing framed full-tape rule" if not candidates else "competing framed rules")
    rule, w, before = candidates[0]
    urest = u[:-len(w)]
    payload = tuple((1 if x > 0 else -1) * (abs(x) + before.end - out) for x in c["output_payload"])
    leading = _power(i - 1, before.end, c["output_power"], limits.max_word_size)
    middle_word = urest + payload + leading + w + core + _power(i - 1, inp, c["input_power"], limits.max_word_size) + p + tail
    middle = RunState(state.n, state.branches[:branch] + (Branch(middle_word, old.coefficient, old.control_json),)
                      + state.branches[branch + 1:])
    request = Request(branch, len(urest) + len(payload) + len(leading), rule["id"], env,
                      evidence={"controller": "framed_upward", "marked_tape_transport": True})
    execution_catalog = _ExecutionCatalog(catalog, env) if bar_mirror else catalog
    raw, equality = apply_rule(middle, request, catalog=execution_catalog, limits=limits, trace=trace)
    base = catalog.get(rule["base_id"] + ".I")
    new_branches = []
    witnesses = []
    for idx, (term, base_term) in enumerate(zip(rule["rhs"], base["rhs"], strict=True)):
        jj, kk, mm, aa, bb = _branch_shape(base_term, env, mirror)
        if bar_mirror:
            aa, bb = -aa, -bb
        out2, inp2 = _anchors(i, jj, kk, mm)
        if before.end != out2:
            raise ControlError("catalog output tape disagrees with frame transport")
        local = upward._literal([term["word"]], env, limits.expansion_budget)
        prefix = _power(i - 1, out2, aa, limits.max_word_size) + _core(state.n, i, jj, kk, mm, policy)
        right = _power(i - 1, inp2, bb, limits.max_word_size)
        if local[:len(prefix) + len(right)] != prefix + right:
            raise ControlError("catalog branch is missing its retained flype frames")
        residual = local[len(prefix) + len(right):]
        after = transport_frame(residual, n=state.n, size=i - 1, start=inp,
                                power=c["input_power"], side="left")
        check_frame(after)
        if after.end != inp2 or after.transformed != residual:
            raise ControlError("catalog tail does not preserve the marked input tape")
        cc = _control(i, jj, kk, mm, c["output_power"] + aa, c["input_power"] + bb,
                      len(urest), len(residual) + len(p), len(tail))
        cc["initial_slider"] = c["initial_slider"]
        cc["mirror_policy"] = policy
        cc["processed_events"] = [{"rule_id": rule["id"], "parameters": env, "rhs_index": idx,
                                   "mirror_policy": policy}] + c["processed_events"]
        cc["output_payload"] = list(payload)
        emitted = raw.branches[branch + idx]
        word = _expression(state.n, cc, urest, residual + p, tail, limits.max_word_size)
        new_branches.append(Branch(word, emitted.coefficient, _json(cc)))
        if trace:
            witnesses.append({"tail_transport": asdict(after), "emitted_output_power": aa,
                              "emitted_input_power": bb})
    following = RunState(state.n, raw.branches[:branch] + tuple(new_branches)
                         + raw.branches[branch + len(new_branches):])
    limits.check(following)
    for b in new_branches:
        _parts(b, state.n, limits.max_word_size)
    powers = [q for b in new_branches for q, _ in b.coefficient]
    event = {"rule_id": rule["id"], "metric_delta": {
        "word_tokens": sum(len(b.word) for b in new_branches) - len(old.word),
        "coefficient_support": sum(len(b.coefficient) for b in new_branches) - len(old.coefficient),
        "max_word_tokens": max(len(b.word) for b in new_branches),
        "q_exponent_min": min(powers), "q_exponent_max": max(powers),
        "max_integer_bits": max(abs(c).bit_length() for b in new_branches for _, c in b.coefficient)}}
    if trace:
        event.update(incoming_frame_transport=asdict(before), catalog_transition=equality,
                     payload_transport={"start": out, "end": before.end,
                                        "word": c["output_payload"], "transformed": list(payload)},
                     branch_frames=witnesses, successor_state=following.to_json())
        if bar_mirror:
            event["semilinear_derivation"] = rule["semilinear"]
    return following, event


def run(state: RunState, *, catalog: Catalog | None = None,
        limits: Limits | None = None, trace: bool = False) -> dict[str, Any]:
    """Run supported framed passes; stopped contexts remain exact checkpoints."""
    if state.branches and state.branches[0].control.get("i") == 2:
        return upward.run(state, catalog=catalog, limits=limits, trace=trace)
    if type(trace) is not bool:
        raise ControlError("trace must be a boolean")
    catalog = Catalog.bundled() if catalog is None else catalog
    limits = Limits() if limits is None else limits
    initial = state
    for b in state.branches:
        _parts(b, state.n, max(limits.max_word_size, len(b.word)))
    history = [metrics(state)]
    first = history[0]
    summary = {"peak_uncollected_terms": len(state.branches), "peak_word_tokens": first["word_tokens"],
               "peak_coefficient_support": first["coefficient_support"],
               **{name: first[name] for name in ("max_word_tokens", "q_exponent_min", "q_exponent_max", "max_integer_bits")}}
    word_tokens, support = first["word_tokens"], first["coefficient_support"]
    counts: Counter[str] = Counter()
    events: list[dict[str, Any]] = []
    steps, diagnostic = 0, None
    while True:
        selected = next((idx for idx, b in enumerate(state.branches) if b.control["phase"] == "framed_zip_up"), None)
        try:
            limits.check(state)
            if selected is None:
                status = "framed_upper_pass_complete"
                break
            if steps >= limits.max_steps:
                raise LimitExceeded("step limit exceeded")
            following, event = step(state, selected, catalog=catalog, limits=limits, trace=trace)
        except (ControlError, LimitExceeded, FrameError) as exc:
            status = "incomplete" if isinstance(exc, LimitExceeded) else "missing_rule"
            diagnostic = {"branch": selected, "message": str(exc), "state": state.to_json()}
            break
        state = following
        steps += 1
        counts[event["rule_id"]] += 1
        d = event["metric_delta"]
        word_tokens, support = word_tokens + d["word_tokens"], support + d["coefficient_support"]
        summary["peak_uncollected_terms"] = max(summary["peak_uncollected_terms"], len(state.branches))
        summary["peak_word_tokens"] = max(summary["peak_word_tokens"], word_tokens)
        summary["peak_coefficient_support"] = max(summary["peak_coefficient_support"], support)
        for name in ("max_word_tokens", "q_exponent_max", "max_integer_bits"):
            summary[name] = max(summary[name], d[name])
        summary["q_exponent_min"] = min(summary["q_exponent_min"], d["q_exponent_min"])
        if trace:
            events.append(event)
            history.append(metrics(state))
    if not trace:
        history.append(metrics(state))
    return {"schema": "zipper-framed-upward-v1", "scope": "framed_upward_pass_only", "trace_enabled": trace,
            "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "dependency_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                  for name in ("frames.py", "catalog.py", "interpreter.py", "upward.py", "semilinear.py")},
            "catalog_sha256": catalog.sha256, "initial": initial.to_json(), "final": state.to_json(),
            "status": status, "diagnostic": diagnostic, "transitions": events, "metrics": history,
            "summary": dict(summary, steps=steps, rule_counts=dict(counts)), "limits": asdict(limits)}


def check_certificate(certificate: Mapping[str, Any], *, catalog: Catalog | None = None) -> bool:
    if certificate.get("schema") == "zipper-normal-upward-v1":
        return upward.check_certificate(certificate, catalog=catalog)
    expected = run(RunState.from_json(certificate["initial"]), catalog=catalog,
                   limits=Limits(**certificate["limits"]), trace=certificate["trace_enabled"])
    if _json(expected) != _json(certificate):
        raise ControlError("framed certificate differs from exact deterministic replay")
    return True


def resume(certificate: Mapping[str, Any], *, catalog: Catalog | None = None,
           limits: Limits | None = None, trace: bool | None = None) -> dict[str, Any]:
    check_certificate(certificate, catalog=catalog)
    return run(RunState.from_json(certificate["final"]), catalog=catalog, limits=limits,
               trace=certificate["trace_enabled"] if trace is None else trace)
