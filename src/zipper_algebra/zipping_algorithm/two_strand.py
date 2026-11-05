"""Complete singleton/zipping/extraction schedule on at most two strands.

The coefficient ring is Z[q,q^-1]. Output terms are 1, sigma^-1, and
sigma^k X_2 (k arbitrary). The full input word is retained in the certificate.
The checker independently substitutes X_2=q sigma^-1+1-q-sigma in the
Laurent group algebra of B_2; no untwisting relation is used.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
from typing import Any

from .catalog import Catalog

_SCHEMA = "zipper-two-strand-run-v2"


def _coefficient(values: dict[int, int]) -> list[dict[str, int]]:
    return [{"q_power": q, "integer": c} for q, c in sorted(values.items()) if c]


def _multiply(a: list[dict[str, int]], b: list[dict[str, int]]) -> list[dict[str, int]]:
    values: dict[int, int] = {}
    for x in a:
        for y in b:
            q = x["q_power"] + y["q_power"]
            values[q] = values.get(q, 0) + x["integer"] * y["integer"]
    return _coefficient(values)


def _collect(output: list[dict[str, Any]], addition: dict[str, Any]) -> None:
    for index, term in enumerate(output):
        if (term["kind"], term["sigma_power"]) == (addition["kind"], addition["sigma_power"]):
            values = {x["q_power"]: x["integer"] for x in term["coefficient"]}
            for x in addition["coefficient"]:
                values[x["q_power"]] = values.get(x["q_power"], 0) + x["integer"]
            coefficient = _coefficient(values)
            if coefficient:
                output[index] = dict(term, coefficient=coefficient)
            else:
                output.pop(index)
            return
    output.append(deepcopy(addition))


def _transition(state: dict[str, Any], n: int, catalog: Catalog) -> dict[str, Any]:
    left, processed, position = state["unprocessed"], state["processed"], state["position"]
    coefficient = state["coefficient"]
    if not left:
        power = sum(processed)
        return {"kind": "extract_singletons", "before": deepcopy(state), "pending_after": [],
                "output_added": [{"kind": "singletons", "sigma_power": power,
                                  "coefficient": deepcopy(coefficient)}],
                "evidence": {"largest_bottom_label": n, "top_position": position,
                             "no_outside_overcrossing": True}}
    letter = left[-1]
    # Positive generators have their left bottom port over the right one.
    # The active bottom label is 2 and its current position is transported.
    offensive = (letter > 0 and position == 2) or (letter < 0 and position == 1)
    side = "left" if position == 2 else "right"
    sign = "positive" if letter > 0 else "negative"
    family = "foe" if offensive else "close"
    rule_id = f"single.{family}.{side}.{sign}.I"
    parameters = {"n": n, "i": 1, "j": position, "l": 1, "k": 0}
    instance = catalog.instantiate(rule_id, parameters)
    remaining = left[:-1]
    if not offensive:
        children = [{"unprocessed": remaining, "processed": [letter] + processed,
                     "position": 3-position, "coefficient": deepcopy(coefficient), "depth": state["depth"]+1}]
        return {"kind": "consume", "rule_id": rule_id, "parameters": parameters,
                "before": deepcopy(state), "pending_after": children, "output_added": [],
                "evidence": {"active_label": 2, "other_label": 1, "active_over": True}}
    grouped: dict[tuple[Any, ...], dict[int, int]] = {}
    for term in instance["rhs"]:
        word = tuple(tuple(x) if isinstance(x, list) else x for x in term["word"])
        values = grouped.setdefault(word, {})
        q = term["q_power"]
        values[q] = values.get(q, 0) + term["integer"]
    pending, output = [], []
    for word, values in grouped.items():
        factor = _multiply(coefficient, _coefficient(values))
        slider = next(x for x in word if isinstance(x, tuple))
        if slider[1] == 2:
            output.append({"kind": "slider", "sigma_power": sum(remaining)+sum(processed),
                           "coefficient": factor})
        else:
            prefix = [x for x in word if type(x) is int]
            pending.append({"unprocessed": remaining + prefix, "processed": list(processed),
                            "position": position, "coefficient": factor, "depth": state["depth"]+1})
    return {"kind": "foe", "rule_id": rule_id, "parameters": parameters,
            "before": deepcopy(state), "pending_after": pending, "output_added": output,
            "evidence": {"active_label": 2, "other_label": 1, "active_over": False,
                         "slider_transport": "X_2 commutes with sigma_1 and its inverse",
                         "extracted_tau": sum(remaining)+sum(processed)}}


def _limits(max_steps: int, max_terms: int) -> None:
    if type(max_steps) is not int or max_steps < 0:
        raise ValueError("max_steps must be a nonnegative integer")
    if type(max_terms) is not int or max_terms < 1:
        raise ValueError("max_terms must be a positive integer")


def _observe(metrics: dict[str, Any], states: list[dict[str, Any]],
             terms: list[dict[str, Any]]) -> None:
    """Accumulate bounded statistics without retaining previous states."""
    for state in states:
        metrics["max_stored_braid_length"] = max(metrics["max_stored_braid_length"],
                                                   len(state["unprocessed"])+len(state["processed"]))
    for item in states + terms:
        coefficient = item["coefficient"]
        metrics["max_coefficient_support"] = max(metrics["max_coefficient_support"], len(coefficient))
        for entry in coefficient:
            metrics["max_integer_bits"] = max(metrics["max_integer_bits"], abs(entry["integer"]).bit_length())
            metrics["q_exponent_min"] = min(metrics["q_exponent_min"], entry["q_power"])
            metrics["q_exponent_max"] = max(metrics["q_exponent_max"], entry["q_power"])


def _advance(run: dict[str, Any], catalog: Catalog, max_steps: int, max_terms: int) -> dict[str, Any]:
    _limits(max_steps, max_terms)
    run["limits"] = {"max_steps": max_steps, "max_terms": max_terms}
    run["status"] = "incomplete"
    run["diagnostic"] = None
    metrics = run["metrics"]
    while run["pending"]:
        if metrics["steps"] >= max_steps:
            run["diagnostic"] = "step_limit"
            break
        transition = _transition(run["pending"][0], run["n"], catalog)
        output = deepcopy(run["output"])
        for term in transition["output_added"]:
            _collect(output, term)
        pending = transition["pending_after"] + run["pending"][1:]
        if len(pending) + len(output) > max_terms:
            run["diagnostic"] = "term_limit"
            break
        run["pending"] = deepcopy(pending)
        run["output"] = output
        if run["trace_enabled"]:
            run["transitions"].append(transition)
        name = transition.get("rule_id", transition["kind"])
        counts = metrics["rule_applications"]
        counts[name] = counts.get(name, 0)+1
        metrics["steps"] += 1
        metrics["uncollected_output_terms"] += len(transition["output_added"])
        metrics["collected_output_terms"] = len(output)
        metrics["pending_terms"] = len(pending)
        metrics["max_pending_terms"] = max(metrics["max_pending_terms"], len(pending))
        metrics["max_total_terms"] = max(metrics["max_total_terms"], len(pending)+len(output))
        metrics["max_branch_count"] = max(metrics["max_branch_count"],
                                            len(transition["pending_after"])+len(transition["output_added"]))
        metrics["maximum_depth"] = max(metrics["maximum_depth"], transition["before"]["depth"]+1)
        _observe(metrics, transition["pending_after"], transition["output_added"] + output)
    if not run["pending"]:
        run["status"] = "complete"
    return run


def _reduce(word: list[int] | tuple[int, ...], *, catalog: Catalog, n: int = 2,
            max_steps: int = 100000, max_terms: int = 100000,
            trace: bool = False) -> dict[str, Any]:
    """Run the complete rank-one/two schedule, retaining a resumable certificate.

    A limit hit leaves the exact remaining R-linear combination in ``pending``.
    ``max_steps`` counts singleton consumption, foe rewrites and extraction.
    Set ``trace=True`` to retain every calculation step. The default stores
    current states and streaming statistics without calculation history.
    """
    if type(trace) is not bool:
        raise ValueError("trace must be a boolean")
    if type(n) is not int or n not in (1, 2):
        raise ValueError("This schedule requires n=1 or n=2")
    if not isinstance(word, (list, tuple)) or any(type(x) is not int or abs(x) != 1 for x in word):
        raise ValueError("Two-strand words require signed generator 1")
    if n == 1 and word:
        raise ValueError("B_1 has no braid generator")
    _limits(max_steps, max_terms)
    run = {"schema": _SCHEMA, "n": n, "input_word": list(word),
           "catalog_sha256": catalog.sha256, "coefficient_ring": "Z[q,q^-1]",
           "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                             for name in ("two_strand.py", "catalog.py")},
           "pending": [{"unprocessed": list(word), "processed": [], "position": n,
                        "coefficient": [{"q_power": 0, "integer": 1}], "depth": 0}],
           "output": [], "transitions": [], "trace_enabled": trace,
           "metrics": {
               "rule_applications": {}, "steps": 0, "uncollected_output_terms": 0,
               "collected_output_terms": 0, "pending_terms": 1,
               "max_pending_terms": 1, "max_total_terms": 1, "max_branch_count": 0,
               "maximum_depth": 0, "max_stored_braid_length": len(word),
               "max_coefficient_support": 1, "max_integer_bits": 1,
               "q_exponent_min": 0, "q_exponent_max": 0,
           }}
    return _advance(run, catalog, max_steps, max_terms)


def _resume(certificate: dict[str, Any], *, catalog: Catalog, max_steps: int = 100000,
            max_terms: int = 100000, trace: bool | None = None) -> dict[str, Any]:
    """Continue the stored queue with enlarged limits and its original trace mode.

    ``trace=None`` preserves the checkpoint mode. An explicit boolean must agree
    with that mode, since an earlier history cannot be recovered from a queue.
    """
    _limits(max_steps, max_terms)
    _check(certificate, catalog=catalog)
    if trace is not None and (type(trace) is not bool or trace != certificate["trace_enabled"]):
        raise ValueError("Resume trace mode must match the checkpoint")
    if max_steps < certificate["metrics"]["steps"]:
        raise ValueError("Resume step limit precedes the existing checkpoint")
    if max_terms < certificate["metrics"]["max_total_terms"]:
        raise ValueError("Resume term limit is below the recorded peak")
    return _advance(deepcopy(certificate), catalog, max_steps, max_terms)


def _expanded(run: dict[str, Any]) -> dict[tuple[int, int], int]:
    values: dict[tuple[int, int], int] = {}
    def add(q: int, power: int, integer: int) -> None:
        key = (q, power)
        values[key] = values.get(key, 0)+integer
    for state in run["pending"]:
        for coefficient in state["coefficient"]:
            add(coefficient["q_power"], sum(state["unprocessed"])+sum(state["processed"]), coefficient["integer"])
    for term in run["output"]:
        for coefficient in term["coefficient"]:
            q, power, integer = coefficient["q_power"], term["sigma_power"], coefficient["integer"]
            if term["kind"] == "singletons":
                add(q, power, integer)
            else:
                for dq, ds, scale in ((1, -1, 1), (0, 0, 1), (1, 0, -1), (0, 1, -1)):
                    add(q+dq, power+ds, integer*scale)
    return {key: integer for key, integer in values.items() if integer}


def _check(certificate: dict[str, Any], *, catalog: Catalog) -> dict[str, Any]:
    """Check exact deterministic replay and independent Laurent-group equality.

    The equality check also applies to incomplete checkpoints and includes their
    pending terms. Equality holds already in Z[q,q^-1][B_2].
    """
    if not isinstance(certificate, dict) or certificate.get("schema") != _SCHEMA:
        raise ValueError("Invalid two-strand certificate schema")
    try:
        expected = _reduce(certificate["input_word"], catalog=catalog, n=certificate["n"],
                                    trace=certificate["trace_enabled"], **certificate["limits"])
        if _expanded(certificate) != {(0, sum(certificate["input_word"])): 1}:
            raise ValueError("Two-strand Laurent-group identity failed")
        if json.dumps(certificate, sort_keys=True, separators=(",", ":")) != json.dumps(expected, sort_keys=True, separators=(",", ":")):
            raise ValueError("Two-strand certificate replay mismatch")
    except (KeyError, TypeError) as exc:
        raise ValueError("Malformed two-strand certificate") from exc
    return {"pass": True, "status": certificate["status"],
            "scope": "Exact deterministic replay and equality in Z[q,q^-1][B_2]"}


class TwoStrandSession:
    """Reuse one validated frozen catalog across independent rank-two runs.

    Only the catalog is shared. Each reduction constructs a fresh queue, output,
    counters and optional trace; no algorithm states or results are memoized.
    The returned certificates equal those from the fresh-loading functions.
    """

    def __init__(self) -> None:
        self._catalog = Catalog.bundled()

    def reduce(self, word: list[int] | tuple[int, ...], *, n: int = 2,
               max_steps: int = 100000, max_terms: int = 100000,
               trace: bool = False) -> dict[str, Any]:
        """Reduce one literal word, with optional detailed calculation history."""
        return _reduce(word, catalog=self._catalog, n=n, max_steps=max_steps,
                       max_terms=max_terms, trace=trace)

    def check(self, certificate: dict[str, Any]) -> dict[str, Any]:
        """Replay one run and independently check its Laurent-group identity."""
        return _check(certificate, catalog=self._catalog)

    def resume(self, certificate: dict[str, Any], *, max_steps: int = 100000,
               max_terms: int = 100000, trace: bool | None = None) -> dict[str, Any]:
        """Validate and continue one checkpoint without altering its trace mode."""
        return _resume(certificate, catalog=self._catalog, max_steps=max_steps,
                       max_terms=max_terms, trace=trace)


def reduce_two_strand(word: list[int] | tuple[int, ...], *, n: int = 2,
                      max_steps: int = 100000, max_terms: int = 100000,
                      trace: bool = False) -> dict[str, Any]:
    """Reduce one literal word using a freshly loaded and validated catalog.

    A limit hit retains an exact resumable pending combination. Set ``trace``
    to retain detailed steps; the default stores streaming statistics only.
    Use ``TwoStrandSession`` when reducing many independent words.
    """
    return TwoStrandSession().reduce(word, n=n, max_steps=max_steps,
                                    max_terms=max_terms, trace=trace)


def resume_two_strand(certificate: dict[str, Any], *, max_steps: int = 100000,
                      max_terms: int = 100000, trace: bool | None = None) -> dict[str, Any]:
    """Validate and resume one checkpoint with a freshly validated catalog."""
    return TwoStrandSession().resume(certificate, max_steps=max_steps,
                                    max_terms=max_terms, trace=trace)


def check_two_strand(certificate: dict[str, Any]) -> dict[str, Any]:
    """Replay and independently check one run with a freshly validated catalog."""
    return TwoStrandSession().check(certificate)
