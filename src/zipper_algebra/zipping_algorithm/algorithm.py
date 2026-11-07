"""Deterministic literal-word zipper scheduling with exact stopped states.

The scheduler scans the largest remaining bottom label, completes both
zipping directions, removes geometrically offensive whole-band crossings,
and extracts components using explicit boundary and nested-ribbon witnesses.
Every input generator is retained initially, including inverse pairs.

An unsupported ribbon context is a control obligation, never a completed
spanning result. Optional traces retain the actual calculation; compact
results retain the exact queue and can be checked by deterministic replay.

The sweep invariant is syntactic at the old-child boundary. Turnaround emits
literal full-band shuffle packets and complete internal tokens. Growth keeps
their ordered stream: internal tokens become child tokens, crossings with
the adopted strand become root packets, and crossings with another outside
strand remain incoming old-child packets. Only the gathered internal payload
is installed in the child slot; the unconsumed prefix and opposite word stay
literal. Endpoint friend fans stay outside the next terminator and are
recorded separately on both sides.
Rotation preserves the literal packet forms. Canonicalization modulo braid
isotopy is needed only for the terminal internal nested-ribbon words.
"""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
from math import comb
from pathlib import Path
from typing import Any, Iterable

from . import framed, upward
from .catalog import Catalog, CatalogError
from .extraction import (
    BoundaryRoute,
    RibbonError,
    RibbonLimit,
    check_split,
    split_layered,
    normalize_canonical_slider,
)
from .interpreter import Branch, Limits, LimitExceeded, RunState
from .frames import FrameError
from .transport import TransportError
from .turnaround import TurnaroundError
from .turnaround import parts as turnaround_parts
from .turnaround import turnaround


class AlgorithmError(ValueError):
    """Invalid input or a failed exact run certificate."""


class ControlObligation(AlgorithmError):
    """A reachable context has no established control transition."""


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _word(value: Iterable[Any]) -> tuple[Any, ...]:
    return tuple(tuple(x) if isinstance(x, list) else x for x in value)


def _encoded(word: Iterable[Any]) -> list[Any]:
    return [list(x) if isinstance(x, tuple) else x for x in word]


def _phi(word: Iterable[Any], n: int) -> tuple[Any, ...]:
    result: list[Any] = []
    for token in reversed(tuple(word)):
        if type(token) is int:
            result.append((1 if token > 0 else -1) * (n - abs(token)))
        else:
            _, size, start, mirrored, reversed_, barred = token
            result.append(
                ("slider", size, n - size - start + 2, mirrored, reversed_, barred)
            )
    return tuple(result)


def _multiply(
    a: Iterable[Iterable[int]], b: Iterable[Iterable[int]]
) -> list[list[int]]:
    values: Counter[int] = Counter()
    for q, c in a:
        for p, d in b:
            values[q + p] += c * d
    return [[q, c] for q, c in sorted(values.items()) if c]


def _common(work: dict[str, Any]) -> dict[str, Any]:
    return {
        key: deepcopy(work[key])
        for key in (
            "rank",
            "coefficient",
            "contexts",
            "world_phi",
            "boundary_left",
            "boundary_right",
            "boundary_left_factors",
            "boundary_right_factors",
            "depth",
            "bottom_labels",
            "top_labels",
        )
    }


def _single(common: dict[str, Any], word: Iterable[int]) -> dict[str, Any]:
    return dict(
        deepcopy(common),
        phase="single",
        left=list(word),
        processed=[],
        position=common["rank"],
    )


def _pass_work(common: dict[str, Any], phase: str, branch: Branch) -> dict[str, Any]:
    return dict(
        deepcopy(common),
        phase=phase,
        branch=branch.to_json(),
        coefficient=[list(x) for x in branch.coefficient],
    )


def _branch_with_coefficient(
    state: RunState, coefficient: Iterable[Iterable[int]]
) -> Branch:
    branch = state.branches[0]
    return Branch(
        branch.word, tuple((q, c) for q, c in coefficient), branch.control_json
    )


def _new_zip(
    common: dict[str, Any],
    *,
    left: Iterable[int],
    right: Iterable[int],
    size: int,
    start: int,
    mirrored: bool = False,
    output_power: int = 0,
    input_power: int = 0,
    output_payload: Iterable[int] = (),
    limits: Limits,
) -> dict[str, Any]:
    if size == 2:
        if output_power or input_power:
            raise ControlObligation("two-strand framing must be the identity")
        # The two mirrored defining polynomials are literally equal.
        state = upward.initialize(
            common["rank"],
            tuple(left) + tuple(output_payload),
            size=2,
            start=start,
            tail=right,
        )
    else:
        state = framed.initialize(
            common["rank"],
            left,
            size=size,
            start=start,
            mirrored=mirrored,
            output_power=output_power,
            input_power=input_power,
            output_payload=output_payload,
            tail=right,
            limits=limits,
            mirrored_policy="bar_mirror",
        )
    return _pass_work(
        common, "zip_up", _branch_with_coefficient(state, common["coefficient"])
    )


def _single_step(
    work: dict[str, Any], catalog: Catalog, limits: Limits
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    n, left, processed, position = (
        work["rank"],
        work["left"],
        work["processed"],
        work["position"],
    )
    if not left:
        slider = ("slider", 1, position, False, False, False)
        local = (slider,) + tuple(processed)
        return [dict(_common(work), phase="extract", word=_encoded(local))], {
            "rule_id": "control.single_scan_complete"
        }
    letter = left[-1]
    index = abs(letter)
    adjacent = position in (index, index + 1)
    active_over = adjacent and ((position == index) == (letter > 0))
    # The active label is the largest remaining bottom label. Every adjacent
    # strand is smaller, so precisely its undercrossings are offensive.
    offensive = adjacent and not active_over
    sign = "positive" if letter > 0 else "negative"
    side = "left" if position == index + 1 else "right"
    family = "foe" if offensive else "close" if adjacent else "far"
    rule_id = f"single.{family}." + (side + "." if adjacent else "") + sign + ".I"
    env = {"n": n, "i": 1, "j": position, "k": 0, "l": index}
    instance = catalog.instantiate(rule_id, env, max_terms=limits.expansion_budget)
    evidence = {
        "rule_id": rule_id,
        "parameters": env,
        "active_label": n,
        "active_position": position,
        "offensive": offensive,
        "active_over": active_over if adjacent else None,
    }
    if not offensive:
        after = deepcopy(work)
        after["left"] = left[:-1]
        after["processed"] = [letter] + processed
        if adjacent:
            after["position"] = index + 1 if position == index else index
        return [after], evidence
    grouped: dict[tuple[Any, ...], Counter[int]] = {}
    for term in instance["rhs"]:
        word = _word(term["word"])
        grouped.setdefault(word, Counter())[term["q_power"]] += term["integer"]
    children = []
    for word, scalar in grouped.items():
        common = _common(work)
        common["coefficient"] = _multiply(common["coefficient"], scalar.items())
        at = next(idx for idx, token in enumerate(word) if isinstance(token, tuple))
        slider = word[at]
        if slider[1] == 2:
            children.append(
                _new_zip(
                    common,
                    left=left[:-1] + list(word[:at]),
                    right=list(word[at + 1 :]) + processed,
                    size=2,
                    start=slider[2],
                    limits=limits,
                )
            )
        else:
            children.append(
                dict(
                    common,
                    phase="single",
                    left=left[:-1] + list(word[:at]),
                    processed=list(word[at + 1 :]) + processed,
                    position=position,
                )
            )
    return children, evidence


def _zip_step(
    work: dict[str, Any], catalog: Catalog, limits: Limits, trace: bool
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    n = work["rank"]
    branch = Branch.from_json(work["branch"])
    c = branch.control
    complete = c["phase"] in ("upper_pass_complete", "framed_upper_pass_complete")
    if not complete:
        engine = upward if c["i"] == 2 else framed
        following, evidence = engine.step(
            RunState(n, (branch,)), catalog=catalog, limits=limits, trace=trace
        )
        return [
            _pass_work(_common(work), work["phase"], b) for b in following.branches
        ], evidence
    turned, certificate = turnaround(
        branch, n=n, catalog=catalog, limits=limits, trace=trace, _want_certificate=trace
    )
    cc, exterior, protected, tail = turnaround_parts(
        turned, n=n, budget=limits.max_word_size
    )
    if protected:
        raise ControlObligation(
            "a zipper traversal must begin with an empty processed word"
        )
    evidence = {
        "rule_id": "control.turnaround",
        "direction": work["phase"],
        "rotation_identity": "phi fixes slider orientation and sends each half twist to the canonical reflected half twist by Garside symmetry",
        "old_child_mirrored": bool(cc["mirrored"] ^ (cc["output_power"] % 2)),
    }
    if trace:
        evidence["turnaround"] = certificate
    # The friend fan opens endpoints at the boundary of the completed child.
    # It remains outside the next active terminator. A protected child braid
    # before that fan moves through it by whole-child naturality, without
    # moving through the slider or becoming a boundary permutation itself.
    payload = tuple(c.get("output_payload", ()))
    friend = (
        upward._friends(c["i"], c["j"], c["k"])
        if c["i"] == 2
        else framed.friend_word(n, c)
    )
    if exterior[: len(payload) + len(friend)] != payload + friend:
        raise ControlObligation(
            "turnaround exterior does not preserve its endpoint friend fan"
        )
    amount = c["k"] if c.get("mirrored", False) else -c["k"]
    payload_shifted = tuple((1 if x > 0 else -1) * (abs(x) + amount) for x in payload)
    body = payload_shifted + exterior[len(payload) + len(friend) :]
    common = _common(work)
    if work["phase"] == "zip_up":
        common["boundary_left"] += list(friend)
        if friend:
            common["boundary_left_factors"].append(
                {"growth_size": cc["i"] + 1, "word": list(friend)}
            )
        jj = n - cc["i"] - cc["j"] + 2
        initial = framed.initialize(
            n,
            _phi(tail, n),
            size=cc["i"],
            start=jj,
            mirrored=cc["mirrored"],
            output_power=cc["input_power"],
            input_power=cc["output_power"],
            tail=_phi(body, n),
            limits=limits,
            mirrored_policy=c.get("mirror_policy", "bar_mirror"),
        )
        return [
            _pass_work(
                common,
                "zip_down",
                _branch_with_coefficient(initial, work["coefficient"]),
            )
        ], evidence
    common["boundary_right"] = list(_phi(friend, n)) + common["boundary_right"]
    if friend:
        common["boundary_right_factors"] = [
            {"growth_size": cc["i"] + 1, "word": list(_phi(friend, n))}
        ] + common["boundary_right_factors"]
    # Replace the down-coordinate exterior by its body and leave the endpoint
    # fan in the protected physical bottom route.
    local = _phi(body + turned.word[len(exterior) :], n)
    return [
        dict(
            common,
            phase="layer",
            word=_encoded(local),
            payload_length=0,
            completed_sides=0,
        )
    ], evidence


def _slider_parts(
    word: tuple[Any, ...],
) -> tuple[tuple[int, ...], tuple[Any, ...], tuple[int, ...]]:
    indices = [idx for idx, token in enumerate(word) if type(token) is not int]
    if len(indices) != 1:
        raise ControlObligation("the active local diagram must have exactly one slider")
    at = indices[0]
    return word[:at], word[at], word[at + 1 :]


def _outside_over(word: Iterable[int], *, n: int, size: int, start: int) -> bool:
    member = [start <= idx < start + size for idx in range(1, n + 1)]
    for letter in reversed(tuple(word)):
        at = abs(letter) - 1
        a, b = member[at : at + 2]
        if a != b and not (a if letter > 0 else b):
            return True
        member[at], member[at + 1] = b, a
    return False


def _layer_step(
    work: dict[str, Any], limits: Limits
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    n = work["rank"]
    word = _word(work["word"])
    left, slider, right = _slider_parts(word)
    _, size, start, mirrored, reversed_, barred = slider
    if reversed_ is not False or barred is not False:
        raise ControlObligation(
            "layering requires an explicit normal or mirrored coupon without reversed or barred flags"
        )
    length = work["payload_length"]
    unprocessed, payload = (left[:-length], left[-length:]) if length else (left, ())
    if not _outside_over(unprocessed, n=n, size=size, start=start):
        if work["completed_sides"] == 0:
            following = dict(
                _common(work),
                phase="layer",
                word=_encoded(_phi(word, n)),
                payload_length=0,
                completed_sides=1,
            )
            following["world_phi"] = not work["world_phi"]
            following["boundary_left"] = list(_phi(work["boundary_right"], n))
            following["boundary_right"] = list(_phi(work["boundary_left"], n))
            following["boundary_left_factors"] = [
                dict(f, word=list(_phi(f["word"], n)))
                for f in reversed(work["boundary_right_factors"])
            ]
            following["boundary_right_factors"] = [
                dict(f, word=list(_phi(f["word"], n)))
                for f in reversed(work["boundary_left_factors"])
            ]
            return [following], {"rule_id": "control.layer_change_side"}
        return [dict(_common(work), phase="extract", word=work["word"])], {
            "rule_id": "control.layer_complete"
        }
    letter = unprocessed[-1]
    index, sign = abs(letter), 1 if letter > 0 else -1
    if start <= index < start + size - 1:
        after = deepcopy(work)
        after["payload_length"] += 1
        return [after], {"rule_id": "control.layer_absorb_internal"}
    if index < start - 1 or index >= start + size:
        remaining = unprocessed[:-1]
        local = remaining + payload + (slider, letter) + right
        return [
            dict(
                _common(work),
                phase="layer",
                word=_encoded(local),
                payload_length=length,
                completed_sides=work["completed_sides"],
            )
        ], {"rule_id": "control.layer_far"}
    move = 1 if index == start + size - 1 else -1
    crossing = tuple(
        sign * x
        for x in (
            range(start, start + size)
            if move == 1
            else range(start + size - 2, start - 2, -1)
        )
    )
    if len(unprocessed) < size or unprocessed[-size:] != crossing:
        raise ControlObligation("layering requires the complete active-band crossing")
    remaining = unprocessed[:-size]
    shifted = tuple((1 if x > 0 else -1) * (abs(x) + move) for x in payload)
    offensive = (move == 1 and sign < 0) or (move == -1 and sign > 0)
    current = dict(
        _common(work),
        phase="layer",
        payload_length=len(shifted),
        completed_sides=work["completed_sides"],
    )
    good = tuple(-x for x in crossing) if offensive else crossing
    moved_slider = ("slider", size, start + move, mirrored, False, False)
    current["word"] = _encoded(remaining + shifted + (moved_slider,) + good + right)
    event = {
        "rule_id": "control.layer_whole_band",
        "crossing": list(crossing),
        "input_anchor": start,
        "output_anchor": start + move,
        "offensive": offensive,
        "orientation": "mirrored" if mirrored else "normal",
    }
    if not offensive:
        return [current], event
    if move == 1:
        scalar, growth_scalar = [[-size, 1]], [[-size, 1]]
        new_start, new_mirror = start, False
        a, b = (1, -1) if mirrored else (0, 0)
    else:
        scalar, growth_scalar = [[size, 1]], [[0, -1]]
        new_start, new_mirror = start - 1, True
        a, b = (0, 0) if mirrored else (-1, 1)
    current["coefficient"] = _multiply(work["coefficient"], scalar)
    common = _common(work)
    common["coefficient"] = _multiply(work["coefficient"], growth_scalar)
    # Only the gathered internal run is a child-supported payload. Keep the
    # remaining prefix and opposite side in their literal order: a previous
    # whole-band crossing with the adopted strand is now a root packet, while
    # one with a still outside strand is an incoming old-child crossing.
    growth = _new_zip(
        common,
        left=remaining,
        right=right,
        size=size + 1,
        start=new_start,
        mirrored=new_mirror,
        output_power=a,
        input_power=b,
        output_payload=shifted,
        limits=limits,
    )
    event["rule_id"] = "layer.route." + ("right" if move == 1 else "left")
    event["growth_frames"] = {"output_power": a, "input_power": b}
    event["growth_child_type"] = {
        "old_child_mirrored": mirrored,
        "output_payload_mirrored": bool(new_mirror ^ (a % 2)),
        "input_payload_mirrored": bool(new_mirror ^ (b % 2)),
    }
    event["equality_reference"] = {
        "catalog_rule": "layer.positive.I"
        if move == 1
        else "layer.opposite_orientation.M",
        "full_band_garside_conjugation": bool(mirrored if move == 1 else not mirrored),
        "scope": "braid-group-algebra slider recursion and whole-band naturality",
    }
    return [growth, current], event


def _bottom_ports(word: tuple[Any, ...], n: int) -> tuple[int, ...]:
    _, slider, right = _slider_parts(word)
    member = [slider[2] <= idx < slider[2] + slider[1] for idx in range(1, n + 1)]
    for letter in right:
        at = abs(letter) - 1
        member[at], member[at + 1] = member[at + 1], member[at]
    return tuple(idx for idx, value in enumerate(member, 1) if value)


def _shift(word: Iterable[Any], amount: int) -> tuple[Any, ...]:
    result: list[Any] = []
    for token in word:
        if type(token) is int:
            result.append((1 if token > 0 else -1) * (abs(token) + amount))
        else:
            result.append((token[0], token[1], token[2] + amount, *token[3:]))
    return tuple(result)


def _lift(word: tuple[Any, ...], contexts: list[dict[str, Any]]) -> tuple[Any, ...]:
    for context in reversed(contexts):
        word = (
            tuple(context["top_route"])
            + _word(context["component_word"])
            + _shift(word, context["size"])
            + tuple(context["bottom_route"])
        )
    return word


def work_expression(work: dict[str, Any]) -> tuple[Any, ...]:
    """Evaluate all saved boundaries and extraction contexts into the full word.

    This is the expression carried by a pending branch, before its Laurent
    coefficient. It performs no algebraic reduction or slider permutation.
    """
    n = work["rank"]
    if work["phase"] == "single":
        local = (
            tuple(work["left"])
            + (("slider", 1, work["position"], False, False, False),)
            + tuple(work["processed"])
        )
    elif work["phase"] in ("zip_up", "zip_down"):
        local = _word(work["branch"]["word"])
        if work["phase"] == "zip_down":
            local = _phi(local, n)
    else:
        local = _word(work["word"])
    local = tuple(work["boundary_left"]) + local + tuple(work["boundary_right"])
    if work["world_phi"]:
        local = _phi(local, n)
    return _lift(local, work["contexts"])


def _project_factors(
    factors: list[dict[str, Any]], member: list[bool]
) -> tuple[list[BoundaryRoute], list[bool]]:
    """Project inherited fan factors to the selected component, in literal order."""
    result = []
    for factor in reversed(factors):
        chronological = []
        for letter in reversed(factor["word"]):
            at = abs(letter) - 1
            a, b = member[at : at + 2]
            if a and b:
                chronological.append((1 if letter > 0 else -1) * (sum(member[:at]) + 1))
            member[at], member[at + 1] = b, a
        projected = tuple(reversed(chronological))
        if projected:
            result.append(BoundaryRoute(factor["growth_size"], projected))
    return list(reversed(result)), member


def _check_generated_fans(factors: list[BoundaryRoute], size: int) -> None:
    """Check sharp fan bounds for scheduler-origin collars after projection."""
    tags: set[int] = set()
    for factor in factors:
        growth = factor.growth_size
        # A saved active-s fan has tag s+1. Nonempty projection requires
        # adopting a friend later, so 3 <= tag <= final component size.
        if (type(growth) is not int or not 3 <= growth <= size or growth in tags
                or len(factor.word) > (growth - 2) * (size - growth + 1)):
            raise ControlObligation(
                "projected endpoint fan exceeds its recorded growth-size bound"
            )
        tags.add(growth)
    if sum(len(factor.word) for factor in factors) > comb(size, 3):
        raise ControlObligation("projected endpoint collar exceeds its size bound")


def _extract_step(
    work: dict[str, Any], *, limits: Limits, trace: bool
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    n, body = work["rank"], _word(work["word"])
    local = tuple(work["boundary_left"]) + body + tuple(work["boundary_right"])
    top_factors, bottom_factors = (
        work["boundary_left_factors"],
        work["boundary_right_factors"],
    )
    if work["world_phi"]:
        local = _phi(local, n)
        body = _phi(body, n)
        top_factors, bottom_factors = (
            [dict(f, word=list(_phi(f["word"], n))) for f in reversed(bottom_factors)],
            [dict(f, word=list(_phi(f["word"], n))) for f in reversed(top_factors)],
        )
    _, source_slider, _ = _slider_parts(local)
    if source_slider[4] is not False or source_slider[5] is not False:
        raise ControlObligation(
            "component extraction requires a coupon without reversed or barred flags"
        )
    selected = _bottom_ports(local, n)
    split = split_layered(local, n=n, bottom_ports=selected)
    check_split(split)
    component = split.component_word
    size = len(selected)
    record: dict[str, Any] = {
        "bottom_ports": list(split.bottom_ports),
        "top_ports": list(split.top_ports),
        "original_bottom_ports": [
            work["bottom_labels"][p - 1] for p in split.bottom_ports
        ],
        "original_top_ports": [work["top_labels"][p - 1] for p in split.top_ports],
    }
    canonical = None
    if size == 1:
        component = ()
        record["twist_powers"] = []
        record.update(
            orientation="normal",
            top_route=[],
            bottom_route=[],
            canonical_top_route=[],
            canonical_bottom_route=[],
            canonical_parities={"upper": [], "lower_inverse": []},
        )
    elif size == 2:
        # The two reflected X2 polynomials coincide, so the base case also
        # emits the default normal coupon rather than a redundant mirror bit.
        slider = ("slider", 2, 1, False, False, False)
        power = sum(x for x in component if type(x) is int)
        component = ((1,) * power if power >= 0 else (-1,) * -power) + (slider,)
        record["twist_powers"] = [power]
        record.update(
            orientation="normal",
            top_route=[],
            bottom_route=[],
            canonical_top_route=[],
            canonical_bottom_route=[],
            canonical_parities={"upper": [], "lower_inverse": []},
        )
    else:
        inner_split = split_layered(body, n=n, bottom_ports=_bottom_ports(body, n))
        inner = inner_split.component_word
        whole_left, whole_slider, whole_right = _slider_parts(component)
        inner_left, inner_slider, inner_right = _slider_parts(inner)
        if (
            whole_slider != inner_slider
            or whole_left[len(whole_left) - len(inner_left) :] != inner_left
            or whole_right[: len(inner_right)] != inner_right
        ):
            raise ControlObligation(
                "inherited endpoint routes do not delimit the component body"
            )
        route_left = whole_left[: len(whole_left) - len(inner_left)]
        route_right = whole_right[len(inner_right) :]
        member = [idx in selected for idx in range(1, n + 1)]
        inherited_bottom, member = _project_factors(bottom_factors, member)
        for token in reversed(body):
            if type(token) is int:
                at = abs(token) - 1
                member[at], member[at + 1] = member[at + 1], member[at]
        inherited_top, _ = _project_factors(top_factors, member)
        if (
            tuple(x for factor in inherited_top for x in factor.word) != route_left
            or tuple(x for factor in inherited_bottom for x in factor.word)
            != route_right
        ):
            raise ControlObligation(
                "projected endpoint fan factors disagree with the full component"
            )
        _check_generated_fans(inherited_top, size)
        _check_generated_fans(inherited_bottom, size)
        witness = normalize_canonical_slider(
            component,
            top_routes=tuple(inherited_top),
            bottom_routes=tuple(inherited_bottom),
            max_image_letters=limits.expansion_budget,
            max_operations=100 * limits.expansion_budget,
        )
        route_bound = comb(size, 3) + comb(size, 2) - 1
        if max(len(witness.top_route), len(witness.bottom_route)) > route_bound:
            raise ControlObligation("terminal finite route exceeds its size bound")
        canonical = witness
        component = witness.normalized_word
        record.update(
            twist_powers=list(witness.twist_powers),
            inherited_top_route=list(route_left),
            inherited_bottom_route=list(route_right),
            orientation="mirrored" if witness.slider[3] else "normal",
            top_route=list(witness.top_route),
            bottom_route=list(witness.bottom_route),
            canonical_top_route=list(witness.top_route[len(route_left) :]),
            canonical_bottom_route=list(
                witness.bottom_route[: len(witness.bottom_route) - len(route_right)]
            ),
            canonical_parities={
                "upper": list(witness.upper.boundary_parities),
                "lower_inverse": list(witness.lower_inverse.boundary_parities),
            },
            output_hierarchy=list(witness.upper.labels),
            input_hierarchy=list(witness.lower_inverse.labels),
        )
    context = {
        "rank": n,
        "size": size,
        "component_word": _encoded(component),
        "top_route": list(split.top_route),
        "bottom_route": list(split.bottom_route),
        "record": record,
    }
    contexts = deepcopy(work["contexts"]) + [context]
    remaining = n - size
    event = {
        "rule_id": "control.extract",
        "component": record,
        "rank_decrease": [n, remaining],
    }
    if trace:
        event["boundary_split"] = asdict(split)
        if canonical is not None:
            event["nested_slider"] = asdict(canonical)
    if remaining:
        common = {
            "rank": remaining,
            "coefficient": work["coefficient"],
            "contexts": contexts,
            "world_phi": False,
            "boundary_left": [],
            "boundary_right": [],
            "boundary_left_factors": [],
            "boundary_right_factors": [],
            "depth": work["depth"],
            "bottom_labels": [
                label
                for p, label in enumerate(work["bottom_labels"], 1)
                if p not in split.bottom_ports
            ],
            "top_labels": [
                label
                for p, label in enumerate(work["top_labels"], 1)
                if p not in split.top_ports
            ],
        }
        return [_single(common, split.complement_word)], [], event
    output = {
        "coefficient": deepcopy(work["coefficient"]),
        "word": _encoded(_lift((), contexts)),
        "components": contexts,
    }
    return [], [output], event


def _size(work: dict[str, Any]) -> int:
    return len(work_expression(work))


def _source() -> dict[str, str]:
    names = (
        "algorithm.py",
        "catalog.py",
        "upward.py",
        "framed.py",
        "frames.py",
        "turnaround.py",
        "transport.py",
        "interpreter.py",
        "extraction.py",
        "semilinear.py",
    )
    return {
        name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in names
    }


class AlgorithmSession:
    """Reuse the frozen local identities, while isolating every input queue."""

    def __init__(self) -> None:
        self._catalog = Catalog.bundled()

    def reduce(
        self,
        word: Iterable[int],
        *,
        n: int,
        limits: Limits | None = None,
        trace: bool = False,
    ) -> dict[str, Any]:
        if type(n) is not int or n < 1 or type(trace) is not bool:
            raise AlgorithmError(
                "positive integer rank and boolean trace flag required"
            )
        if not isinstance(word, (list, tuple)):
            raise AlgorithmError("a finite literal word list or tuple is required")
        literal = tuple(word)
        if any(type(x) is not int or not 1 <= abs(x) < n for x in literal):
            raise AlgorithmError(
                "literal input must contain signed generators in the stated rank"
            )
        limits = Limits() if limits is None else limits
        queue = [
            _single(
                {
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
                },
                literal,
            )
        ]
        output: list[dict[str, Any]] = []
        transitions = []
        counts: Counter[str] = Counter()
        steps, diagnostic, status = 0, None, "complete"
        max_terms, max_word = 1, len(literal) + 1
        peak_support, integer_bits, qmin, qmax, depth, branch_width = 1, 1, 0, 0, 0, 0
        while queue:
            before = queue[0]
            added: list[dict[str, Any]]
            try:
                if steps >= limits.max_steps:
                    raise LimitExceeded("step limit exceeded")
                if len(queue) + len(output) > limits.max_terms:
                    raise LimitExceeded("term limit exceeded")
                if max(_size(x) for x in queue) > limits.max_word_size:
                    raise LimitExceeded("word size limit exceeded")
                for item in queue:
                    if len(item["coefficient"]) > limits.max_coefficient_support:
                        raise LimitExceeded("coefficient support limit exceeded")
                    if any(
                        abs(c).bit_length() > limits.max_integer_bits
                        for _, c in item["coefficient"]
                    ):
                        raise LimitExceeded("coefficient integer bit limit exceeded")
                if before["phase"] == "single":
                    children, event = _single_step(before, self._catalog, limits)
                    added = []
                elif before["phase"] in ("zip_up", "zip_down"):
                    children, event = _zip_step(before, self._catalog, limits, trace)
                    added = []
                elif before["phase"] == "layer":
                    children, event = _layer_step(before, limits)
                    added = []
                else:
                    children, added, event = _extract_step(
                        before, limits=limits, trace=trace
                    )
                if (
                    len(queue) - 1 + len(children) + len(output) + len(added)
                    > limits.max_terms
                ):
                    raise LimitExceeded("term limit exceeded")
                if any(_size(child) > limits.max_word_size for child in children):
                    raise LimitExceeded("word size limit exceeded")
                if any(len(item["word"]) > limits.max_word_size for item in added):
                    raise LimitExceeded("output word size limit exceeded")
                for item in children + added:
                    if len(item["coefficient"]) > limits.max_coefficient_support:
                        raise LimitExceeded("coefficient support limit exceeded")
                    if any(
                        abs(c).bit_length() > limits.max_integer_bits
                        for _, c in item["coefficient"]
                    ):
                        raise LimitExceeded("coefficient integer bit limit exceeded")
            except (
                LimitExceeded,
                ControlObligation,
                upward.ControlError,
                RibbonError,
                FrameError,
                TurnaroundError,
                TransportError,
                CatalogError,
            ) as exc:
                status = (
                    "incomplete"
                    if isinstance(exc, (LimitExceeded, RibbonLimit))
                    or isinstance(exc, CatalogError)
                    and "budget" in str(exc)
                    else "control_obligation"
                )
                diagnostic = {
                    "message": str(exc),
                    "phase": before["phase"],
                    "state": deepcopy(before),
                }
                break
            for child in children:
                child["depth"] = before["depth"] + 1
            queue = children + queue[1:]
            output.extend(added)
            steps += 1
            name = event.get(
                "rule_id",
                event.get("catalog_transition", {}).get("request", {}).get("rule_id"),
            )
            counts[str(name)] += 1
            max_terms = max(max_terms, len(queue) + len(output))
            max_word = max(
                max_word,
                max((_size(x) for x in children), default=0),
                max((len(item["word"]) for item in added), default=0),
            )
            depth = max(depth, before["depth"] + 1)
            branch_width = max(branch_width, len(children) + len(added))
            for item in children + added:
                peak_support = max(peak_support, len(item["coefficient"]))
                for q, c in item["coefficient"]:
                    qmin, qmax = min(qmin, q), max(qmax, q)
                    integer_bits = max(integer_bits, abs(c).bit_length())
            if trace:
                transitions.append(
                    {
                        "predecessor": deepcopy(before),
                        "successors": deepcopy(children),
                        "output_added": deepcopy(added),
                        "evidence": event,
                    }
                )
        return {
            "schema": "zipper-algorithm-run-v1",
            "scope": "literal_word_spanning_procedure",
            "n": n,
            "input_word": list(literal),
            "coefficient_ring": "Z[q,q^-1]",
            "source_sha256": _source(),
            "catalog_sha256": self._catalog.sha256,
            "trace_enabled": trace,
            "limits": asdict(limits),
            "status": status,
            "diagnostic": diagnostic,
            "pending": queue,
            "output": output,
            "transitions": transitions,
            "metrics": {
                "steps": steps,
                "rule_applications": dict(counts),
                "max_terms": max_terms,
                "max_local_word_tokens": max_word,
                "output_terms": len(output),
                "pending_terms": len(queue),
                "max_coefficient_support": peak_support,
                "max_integer_bits": integer_bits,
                "q_exponent_min": qmin,
                "q_exponent_max": qmax,
                "maximum_depth": depth,
                "max_branch_count": branch_width,
            },
        }

    def check(self, certificate: dict[str, Any]) -> bool:
        try:
            expected = self.reduce(
                certificate["input_word"],
                n=certificate["n"],
                limits=Limits(**certificate["limits"]),
                trace=certificate["trace_enabled"],
            )
        except (KeyError, TypeError) as exc:
            raise AlgorithmError("malformed zipper algorithm certificate") from exc
        if _json(certificate) != _json(expected):
            raise AlgorithmError(
                "zipper algorithm certificate differs from deterministic replay"
            )
        return True


def reduce_word(
    word: Iterable[int], *, n: int, limits: Limits | None = None, trace: bool = False
) -> dict[str, Any]:
    return AlgorithmSession().reduce(word, n=n, limits=limits, trace=trace)


def check_certificate(certificate: dict[str, Any]) -> bool:
    return AlgorithmSession().check(certificate)


def resume(
    certificate: dict[str, Any], *, limits: Limits | None = None
) -> dict[str, Any]:
    """Check a stopped prefix, then regenerate it under the enlarged limits.

    The exact queue is saved for inspection; this reference continuation checks
    and replays its deterministic prefix before extending it.
    """
    session = AlgorithmSession()
    session.check(certificate)
    following = Limits(**certificate["limits"]) if limits is None else limits
    metrics = certificate["metrics"]
    requirements = {
        "max_steps": metrics["steps"],
        "max_terms": metrics["max_terms"],
        "max_word_size": metrics["max_local_word_tokens"],
        "max_coefficient_support": metrics["max_coefficient_support"],
        "max_integer_bits": metrics["max_integer_bits"],
        "expansion_budget": certificate["limits"]["expansion_budget"],
    }
    if any(getattr(following, name) < used for name, used in requirements.items()):
        raise AlgorithmError("continuation limits must retain the verified prefix")
    return session.reduce(
        certificate["input_word"],
        n=certificate["n"],
        limits=following,
        trace=certificate["trace_enabled"],
    )
