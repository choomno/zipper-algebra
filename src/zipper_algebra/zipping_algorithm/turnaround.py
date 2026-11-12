"""Transport a completed framed zipper pass through its recorded tails.

Only single-branch braid identities are used here. An emitted whole-band
tail transports every internal operator, including both retained frames.
An emitted flype packet is unwrapped in reverse chronological order; its
other factors stay in the exterior braid. No branching zipper identity
is inverted, and no unrestricted fixed-position slider commute is used.
"""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .catalog import Catalog, condition, expand_polynomial
from .interpreter import Branch, Limits, LimitExceeded, RunState
from .transport import check_transport, transport_band


class TurnaroundError(ValueError):
    """The saved typed tails do not describe the stated completed pass."""


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _word(
    terms: list[dict[str, Any]], env: dict[str, int], budget: int
) -> tuple[Any, ...]:
    polynomial = expand_polynomial(terms, env, max_terms=budget)
    if len(polynomial) != 1:
        raise TurnaroundError("a monomial packet word is required")
    (power, word), scalar = next(iter(polynomial.items()))
    if power != 0 or scalar != 1:
        raise TurnaroundError("packet word must have coefficient one")
    return word


def _event(
    event: Any, n: int, i: int, catalog: Catalog, limits: Limits
) -> dict[str, Any]:
    """Reconstruct a tail's type from its original accepted catalog branch."""
    from .framed import _anchors, _branch_shape, _core, _power, rule_for_event

    if (
        not isinstance(event, dict)
        or set(event) != {"rule_id", "parameters", "rhs_index", "mirror_policy"}
        or not isinstance(event["rule_id"], str)
        or not isinstance(event["parameters"], dict)
        or set(event["parameters"]) != {"n", "i", "j", "k", "l", "eps"}
        or any(type(v) is not int for v in event["parameters"].values())
        or type(event["rhs_index"]) is not int
        or event["mirror_policy"] not in ("M", "bar_mirror")
    ):
        raise TurnaroundError("invalid processed-tail event")
    env = event["parameters"]
    name = event["rule_id"]
    if (
        env.get("n") != n
        or env.get("i") != i
        or not name.startswith("zip.")
        or not name.endswith((".I", ".M", ".BM"))
        or (name.endswith(".BM") and event["mirror_policy"] != "bar_mirror")
        or (name.endswith(".M") and event["mirror_policy"] != "M")
    ):
        raise TurnaroundError("tail event has incompatible rank or orientation")
    rule = rule_for_event(catalog, event)
    if (
        rule.get("execution", {}).get("algebra_verified") is not True
        or not condition(rule["condition"], env)
        or not 0 <= event["rhs_index"] < len(rule["rhs"])
    ):
        raise TurnaroundError("tail event lacks a valid accepted rule branch")
    mirror = name.endswith((".M", ".BM"))
    barred = name.endswith(".BM")
    old_j = n - i - env["j"] + 2 if mirror else env["j"]
    base = catalog.get(rule["base_id"] + ".I")
    index = event["rhs_index"]
    jj, kk, mm, aa, bb = _branch_shape(base["rhs"][index], env, mirror)
    if barred:
        aa, bb = -aa, -bb
    output, inp = _anchors(i, jj, kk, mm)
    term = rule["rhs"][index]
    local = _word(
        [{"coefficient": [{"integer": 1, "q_power": 0}], "word": term["word"]}],
        env,
        limits.expansion_budget,
    )
    prefix = (
        _power(i - 1, output, aa, limits.max_word_size)
        + _core(n, i, jj, kk, mm, policy=event["mirror_policy"])
        + _power(i - 1, inp, bb, limits.max_word_size)
    )
    if local[: len(prefix)] != prefix:
        raise TurnaroundError("tail event does not have the catalog frame shape")
    return {
        "env": dict(env),
        "pre_j": old_j,
        "pre_mirrored": mirror,
        "post_j": jj,
        "post_mirrored": mm,
        "output_delta": aa,
        "input_delta": bb,
        "tail": local[len(prefix) :],
        "mirror_policy": event["mirror_policy"],
    }


def parts(
    branch: Branch, *, n: int, budget: int = 100000
) -> tuple[dict[str, Any], tuple[Any, ...], tuple[Any, ...], tuple[Any, ...]]:
    """Read an actual turnaround output, checking all its word boundaries."""
    from .framed import _anchors, _power

    RunState(n, (branch,))
    c = branch.control
    keys = {
        "phase",
        "i",
        "j",
        "mirrored",
        "mirror_policy",
        "output_power",
        "input_power",
        "exterior_length",
        "payload_length",
        "routing_length",
        "transported_length",
        "processed_length",
        "tail_length",
    }
    if (
        set(c) != keys
        or c["phase"] != "turnaround_complete"
        or type(c["mirrored"]) is not bool
        or c["mirror_policy"] not in ("M", "bar_mirror")
        or any(
            type(c[k]) is not int for k in keys - {"phase", "mirrored", "mirror_policy"}
        )
        or not 2 <= c["i"] <= n
        or not 1 <= c["j"] <= n - c["i"] + 1
        or min(
            c["exterior_length"],
            c["payload_length"],
            c["routing_length"],
            c["transported_length"],
            c["processed_length"],
            c["tail_length"],
        )
        < 0
        or c["exterior_length"]
        != c["payload_length"] + c["routing_length"] + c["transported_length"]
    ):
        raise TurnaroundError("invalid turnaround output boundaries")
    left = c["exterior_length"]
    end = len(branch.word) - c["tail_length"] - c["processed_length"]
    exterior, processed, tail = (
        branch.word[:left],
        branch.word[end : end + c["processed_length"]],
        branch.word[end + c["processed_length"] :],
    )
    output, inp = _anchors(c["i"], c["j"], 0, c["mirrored"])
    core = (
        _power(c["i"] - 1, output, c["output_power"], budget)
        + (("slider", c["i"], c["j"], c["mirrored"], False, False),)
        + _power(c["i"] - 1, inp, c["input_power"], budget)
    )
    if branch.word != exterior + core + processed + tail or any(
        type(x) is not int for x in exterior + processed + tail
    ):
        raise TurnaroundError("turnaround output does not match its full expression")
    return c, exterior, processed, tail


def turnaround(
    branch: Branch,
    *,
    n: int,
    catalog: Catalog | None = None,
    limits: Limits | None = None,
    trace: bool = False,
    _want_certificate: bool = True,
) -> tuple[Branch, dict[str, Any]]:
    """Unwind mandatory processed events of one completed framed pass.

    The original protected processed word and exterior tail are retained.
    The resulting exterior includes the final friends and every transformed
    processed tail. Arbitrary initial frame powers are restored exactly.
    """
    from .framed import _anchors, _core, _parts, _power, rule_for_event

    if type(trace) is not bool:
        raise TurnaroundError("trace must be a boolean")
    catalog = Catalog.bundled() if catalog is None else catalog
    limits = Limits() if limits is None else limits
    limits.check(RunState(n, (branch,)))
    if branch.control.get("i") == 2:
        from .upward import _friends, _parts as two_parts

        c, u, p, tail = two_parts(branch, n)
        if u:
            raise TurnaroundError("completed upward pass required")
        band = transport_band(p, n=n, size=2, start=c["j"])
        check_transport(band)
        exterior = _friends(2, c["j"], c["k"]) + p
        cc = {
            "phase": "turnaround_complete",
            "i": 2,
            "j": band.end,
            "mirrored": False,
            "mirror_policy": "M",
            "output_power": 0,
            "input_power": 0,
            "exterior_length": len(exterior),
            "payload_length": 0,
            "routing_length": c["k"],
            "transported_length": len(p),
            "processed_length": 0,
            "tail_length": len(tail),
        }
        following = Branch(
            exterior + (("slider", 2, band.end, False, False, False),) + tail,
            branch.coefficient,
            _json(cc),
        )
        evidence: list[dict[str, Any]] = (
            [{"band_transport": asdict(band)}] if trace else []
        )
        if not (trace or _want_certificate):
            return following, {}
        return following, _certificate(
            branch, following, n, catalog, limits, trace, evidence, len(band.steps)
        )
    c, u, p, tail = _parts(branch, n, limits.max_word_size)
    if u or c["phase"] != "framed_upper_pass_complete":
        raise TurnaroundError("completed framed upward pass required")
    i, j, mirror = c["i"], c["j"], c["mirrored"]
    policy = c["mirror_policy"]
    a, b = c["output_power"], c["input_power"]
    # F carries the output tape; naturality moves its frame through F
    # toward X, retaining F as exterior braid and preserving its exponent.
    exterior = (
        tuple(c.get("output_payload", ()))
        + _core(n, i, j, c["k"], mirror, policy=policy)[:-1]
    )
    payload_length = len(c.get("output_payload", ()))
    routing_length = (i - 1) * c["k"]
    events = c["processed_events"]
    if len(events) > limits.max_steps:
        raise LimitExceeded("turnaround event limit exceeded")
    evidence = []
    for event in events:
        typed = _event(event, n, i, catalog, limits)
        if (
            typed["post_j"] != j
            or typed["post_mirrored"] != mirror
            or typed["mirror_policy"] != policy
            or p[: len(typed["tail"])] != typed["tail"]
        ):
            raise TurnaroundError("tail event disagrees with current slider or word")
        T = typed["tail"]
        p = p[len(T) :]
        aa, bb = typed["output_delta"], typed["input_delta"]
        packet_id = None
        if aa or bb:
            if aa not in (-1, 1) or bb != -aa or typed["pre_mirrored"] == mirror:
                raise TurnaroundError("invalid emitted flype frame increments")
            barred = typed["pre_mirrored"] and policy == "bar_mirror"
            sign = "positive" if aa == (1 if barred else -1) else "negative"
            variant = ".BM" if barred else ".M" if typed["pre_mirrored"] else ".I"
            packet_id = "zip.flype." + sign + variant
            env = dict(
                typed["env"], j=n - i - j + 2 if typed["pre_mirrored"] else j, k=0
            )
            packet = rule_for_event(
                catalog,
                {
                    "rule_id": packet_id,
                    "parameters": env,
                    "rhs_index": 0,
                    "mirror_policy": policy,
                },
            )
            old_slider = ("slider", i, j, typed["pre_mirrored"], False, False)
            incoming = _word(packet["lhs"], env, limits.expansion_budget)
            if incoming[-1:] != (old_slider,):
                raise TurnaroundError("flype packet has an invalid input slider")
            W = incoming[:-1]
            packet_word = _word(packet["rhs"], env, limits.expansion_budget)
            output, inp = _anchors(i, j, 0, mirror)
            packet_core = (
                _power(i - 1, output, aa, limits.max_word_size)
                + (("slider", i, j, mirror, False, False),)
                + _power(i - 1, inp, bb, limits.max_word_size)
            )
            if packet_word[: len(packet_core)] != packet_core:
                raise TurnaroundError("flype packet has an invalid retained frame")
            C = packet_word[len(packet_core) :]
            if len(C) != i - 1 or T[: len(C)] != C:
                raise TurnaroundError(
                    "tail is missing its complete outgoing flype tape"
                )
            R = T[len(C) :]
            a, b = a - aa, b - bb
            mirror = typed["pre_mirrored"]
        else:
            if typed["pre_mirrored"] != mirror:
                raise TurnaroundError(
                    "orientation change requires an explicit flype packet"
                )
            W, R = (), T
        band = transport_band(R, n=n, size=i, start=j, mirrored=mirror)
        check_transport(band)
        if band.end != typed["pre_j"]:
            raise TurnaroundError(
                "residual tail does not return to the previous band anchor"
            )
        exterior += W + R
        j = band.end
        if trace:
            evidence.append(
                {
                    "event": event,
                    "packet_rule": packet_id,
                    "band_transport": asdict(band),
                    "output_power": a,
                    "input_power": b,
                    "j": j,
                    "mirrored": mirror,
                }
            )
    initial = c["initial_slider"]
    if (
        j != initial["j"]
        or mirror != initial["mirrored"]
        or a != initial["output_power"]
        or b != initial["input_power"]
        or len(p) != initial["processed_length"]
        or policy != initial["mirror_policy"]
    ):
        raise TurnaroundError("typed tails do not restore the initial slider data")
    output, inp = _anchors(i, j, 0, mirror)
    word = (
        exterior
        + _power(i - 1, output, a, limits.max_word_size)
        + (("slider", i, j, mirror, False, False),)
        + _power(i - 1, inp, b, limits.max_word_size)
        + p
        + tail
    )
    cc = {
        "phase": "turnaround_complete",
        "i": i,
        "j": j,
        "mirrored": mirror,
        "mirror_policy": policy,
        "output_power": a,
        "input_power": b,
        "exterior_length": len(exterior),
        "payload_length": payload_length,
        "routing_length": routing_length,
        "transported_length": len(exterior) - payload_length - routing_length,
        "processed_length": len(p),
        "tail_length": len(tail),
    }
    following = Branch(word, branch.coefficient, _json(cc))
    limits.check(RunState(n, (following,)))
    parts(following, n=n, budget=limits.max_word_size)
    if not (trace or _want_certificate):
        return following, {}
    return following, _certificate(
        branch, following, n, catalog, limits, trace, evidence, len(events)
    )


def _certificate(
    branch: Branch,
    following: Branch,
    n: int,
    catalog: Catalog,
    limits: Limits,
    trace: bool,
    evidence: list[dict[str, Any]],
    events: int,
) -> dict[str, Any]:
    return {
        "schema": "zipper-turnaround-v1",
        "n": n,
        "trace_enabled": trace,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "dependency_sha256": {
            name: hashlib.sha256(
                Path(__file__).with_name(name).read_bytes()
            ).hexdigest()
            for name in (
                "framed.py",
                "frames.py",
                "transport.py",
                "catalog.py",
                "upward.py",
                "semilinear.py",
            )
        },
        "catalog_sha256": catalog.sha256,
        "initial": branch.to_json(),
        "final": following.to_json(),
        "events": evidence,
        "summary": {
            "events": events,
            "word_tokens_before": len(branch.word),
            "word_tokens_after": len(following.word),
        },
        "limits": asdict(limits),
    }


def check_certificate(
    certificate: Mapping[str, Any], *, catalog: Catalog | None = None
) -> bool:
    """Reload and recheck every typed tail, frame update and whole-band move."""
    if certificate.get("schema") != "zipper-turnaround-v1":
        raise TurnaroundError("unsupported turnaround certificate schema")
    _, expected = turnaround(
        Branch.from_json(certificate["initial"]),
        n=certificate["n"],
        catalog=catalog,
        limits=Limits(**certificate["limits"]),
        trace=certificate["trace_enabled"],
    )
    if _json(expected) != _json(certificate):
        raise TurnaroundError(
            "turnaround certificate differs from exact typed-tail replay"
        )
    return True
