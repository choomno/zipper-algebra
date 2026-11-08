"""Versioned exact formal-output interchange for the current zipper procedure.

This normal form collects identical full terms and sorts them. It imposes no
braid identities, no free reduction, and no permutation of components. Every
term/record extension field is retained in its identity. Coefficient pairs are
canonical decimal strings, so arbitrary-precision coefficients survive JSON
tooling. Other fields and extension integers retain ordinary JSON semantics.

Structural validation is NOT a proof that an input was computed correctly.
Keep the original certificate and producer identity separately and replay them
with an independent verifier. A partial output sum is labeled partial and has
a different identity from the same sum claimed complete. It is not a resumable
checkpoint: pending continuation states remain in the original certificate.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from typing import Any

SCHEMA = "zipper-formal-output-normal-form-v1"
RING = "Z[q,q^-1]"


def _integer(value: Any, *, positive: bool = False) -> int:
    if type(value) is not int or (positive and value < 1):
        raise ValueError("integer required" if not positive else "positive integer required")
    return value


def _decimal(value: int) -> str:
    if not value:
        return "0"
    negative, value = value < 0, abs(value)
    chunks = []
    while value:
        value, tail = divmod(value, 10**9)
        chunks.append(tail)
    return ("-" if negative else "") + str(chunks[-1]) + "".join(
        f"{chunk:09d}" for chunk in reversed(chunks[:-1]))


def _exact(value: Any) -> int:
    if type(value) is int:
        return value
    if type(value) is not str or not value:
        raise ValueError("exact integer or canonical decimal string required")
    digits = value[1:] if value.startswith("-") else value
    if (not digits or not digits.isascii() or not digits.isdigit()
            or (len(digits) > 1 and digits[0] == "0") or value == "-0"):
        raise ValueError("noncanonical decimal integer")
    result = 0
    for at in range(0, len(digits), 9):
        chunk = digits[at:at + 9]
        result = result * 10**len(chunk) + int(chunk)
    return -result if value.startswith("-") else result


def _json_value(value: Any) -> Any:
    if value is None or type(value) in (str, int, bool):
        return value
    if type(value) in (list, tuple):
        return [_json_value(x) for x in value]
    if type(value) is dict and all(type(k) is str for k in value):
        return {k: _json_value(v) for k, v in value.items()}
    raise ValueError("exact JSON values required; floats and non-string keys are invalid")


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _word(value: Any, n: int, *, generators_only: bool = False) -> list[Any]:
    if type(value) is not list:
        raise ValueError("word must be an array")
    for token in value:
        if type(token) is int:
            if not 0 < abs(token) < n:
                raise ValueError("generator outside word rank")
        elif (not generators_only and type(token) is list and len(token) == 6
              and token[0] == "slider" and type(token[1]) is int
              and type(token[2]) is int and 1 <= token[1] <= n
              and 1 <= token[2] <= n - token[1] + 1
              and all(type(flag) is bool for flag in token[3:])):
            pass
        else:
            raise ValueError("invalid word token")
    return value


def _ports(value: Any, *, maximum: int, size: int) -> list[int]:
    if (type(value) is not list or len(value) != size
            or any(type(x) is not int or not 1 <= x <= maximum for x in value)
            or value != sorted(set(value))):
        raise ValueError("ports must be an increasing set of the component size")
    return value


def _expected_over_route(selected: list[int]) -> list[int]:
    """Stable-partition route directly from ports, without simulating a braid.

    Selected port p in position d moves across indices p-1,...,d. Earlier
    selected ports are smaller and their moves leave p's current position
    unchanged. Reverse the chronological list to obtain algebraic word order.
    """
    return [-index for destination, port in reversed(list(enumerate(selected, 1)))
            for index in range(destination, port)]


def _shift(word: list[Any], size: int) -> list[Any]:
    return [(1 if t > 0 else -1) * (abs(t) + size) if type(t) is int
            else [t[0], t[1], t[2] + size, *t[3:]] for t in word]


def _term(value: Any, n: int) -> tuple[dict[str, Any], dict[int, int]]:
    term = _json_value(value)
    if type(term) is not dict or not {"coefficient", "word", "components"} <= term.keys():
        raise ValueError("term requires coefficient, word, and components")
    coefficient = term.pop("coefficient")
    if type(coefficient) is not list:
        raise ValueError("coefficient must be an array of exponent/scalar pairs")
    scalars: defaultdict[int, int] = defaultdict(int)
    for pair in coefficient:
        if type(pair) is not list or len(pair) != 2:
            raise ValueError("coefficient requires exponent/scalar pairs")
        scalars[_exact(pair[0])] += _exact(pair[1])
    _word(term["word"], n)
    components = term["components"]
    if type(components) is not list:
        raise ValueError("components must be an ordered array")
    top = bottom = list(range(1, n + 1))
    for component in components:
        required = {"rank", "size", "component_word", "top_route", "bottom_route", "record"}
        if type(component) is not dict or not required <= component.keys():
            raise ValueError("incomplete component")
        rank = _integer(component["rank"], positive=True)
        size = _integer(component["size"], positive=True)
        if rank != len(bottom) or size > rank:
            raise ValueError("component ranks must follow extraction order")
        inner = _word(component["component_word"], size)
        sliders = [t for t in inner if type(t) is list]
        if (size == 1 and inner) or (size > 1 and (
                len(sliders) != 1 or sliders[0][1:3] != [size, 1])):
            raise ValueError("component requires its single full-size slider")
        for side in ("top", "bottom"):
            _word(component[side + "_route"], rank, generators_only=True)
        record = component["record"]
        required_record = {"bottom_ports", "top_ports", "original_bottom_ports", "original_top_ports",
                           "twist_powers", "orientation", "top_route", "bottom_route",
                           "canonical_top_route", "canonical_bottom_route", "canonical_parities"}
        if type(record) is not dict or not required_record <= record.keys():
            raise ValueError("incomplete component record")
        local_top = _ports(record["top_ports"], maximum=rank, size=size)
        local_bottom = _ports(record["bottom_ports"], maximum=rank, size=size)
        expected_bottom = _expected_over_route(local_bottom)
        expected_top = [-x for x in reversed(_expected_over_route(local_top))]
        if (component["bottom_route"] != expected_bottom
                or component["top_route"] != expected_top):
            raise ValueError("ambient component routes disagree with canonical port routes")
        for side, local, labels in (("top", local_top, top), ("bottom", local_bottom, bottom)):
            original = _ports(record["original_" + side + "_ports"], maximum=n, size=size)
            if original != [labels[p - 1] for p in local]:
                raise ValueError("original ports disagree with inherited endpoint labels")
        top = [x for p, x in enumerate(top, 1) if p not in local_top]
        bottom = [x for p, x in enumerate(bottom, 1) if p not in local_bottom]
        powers = record["twist_powers"]
        if type(powers) is not list or len(powers) != size - 1 or any(type(x) is not int for x in powers):
            raise ValueError("one integral twist power per nontrivial nested size required")
        if record["orientation"] not in ("normal", "mirrored"):
            raise ValueError("invalid orientation")
        if record["orientation"] != ("mirrored" if sliders and sliders[0][3] else "normal"):
            raise ValueError("record orientation disagrees with component slider")
        for route in ("top_route", "bottom_route", "canonical_top_route", "canonical_bottom_route"):
            _word(record[route], size, generators_only=True)
        parities = record["canonical_parities"]
        if type(parities) is not dict or not {"upper", "lower_inverse"} <= parities.keys():
            raise ValueError("canonical parities required")
        for name in ("upper", "lower_inverse"):
            values = parities[name]
            if (type(values) is not list or len(values) != max(0, size - 2)
                    or any(type(x) is not int or x not in (0, 1) for x in values)):
                raise ValueError("invalid canonical parities")
        if size > 2:
            for side in ("top", "bottom"):
                _word(record.get("inherited_" + side + "_route"), size, generators_only=True)
            for name in ("input_hierarchy", "output_hierarchy"):
                hierarchy = record.get(name)
                if (type(hierarchy) is not list or any(type(x) is not int for x in hierarchy)
                        or sorted(hierarchy) != list(range(1, size + 1))):
                    raise ValueError("component hierarchy must permute its local strands")
    if bottom:
        raise ValueError("terminal components must partition all endpoints")
    reconstructed: list[Any] = []
    for component in reversed(components):
        reconstructed = (component["top_route"] + component["component_word"]
                         + _shift(reconstructed, component["size"]) + component["bottom_route"])
    if reconstructed != term["word"]:
        raise ValueError("terminal word disagrees with ordered component reconstruction")
    return term, dict(scalars)


def canonical_output(output: Any, *, n: int, complete: bool) -> dict[str, Any]:
    """Collect an explicitly complete or partial sum without algebraic rewrites."""
    _integer(n, positive=True)
    if type(complete) is not bool or type(output) not in (list, tuple):
        raise ValueError("explicit boolean completeness and output array required")
    terms: dict[str, dict[str, Any]] = {}
    coefficients: dict[str, defaultdict[int, int]] = {}
    for value in output:
        term, scalar = _term(value, n)
        key = _json(term)
        terms[key] = term
        total = coefficients.setdefault(key, defaultdict(int))
        for exponent, integer in scalar.items():
            total[exponent] += integer
    result = []
    for key in sorted(terms):
        coefficient = [[_decimal(q), _decimal(c)] for q, c in sorted(coefficients[key].items()) if c]
        if coefficient:
            result.append(dict(terms[key], coefficient=coefficient))
    return {"schema": SCHEMA, "coefficient_ring": RING, "n": n,
            "completion": "complete" if complete else "partial", "terms": result}


def from_certificate(certificate: Any) -> dict[str, Any]:
    """Import current run formats structurally; this does not replay the run."""
    if (type(certificate) is not dict or certificate.get("schema") not in (
            "zipper-algorithm-run-v1", "zipper-collected-run-v1", "zipper-native-run-v1")
            or certificate.get("coefficient_ring") != RING
            or certificate.get("status") not in ("complete", "incomplete", "control_obligation")
            or type(certificate.get("pending")) is not list):
        raise ValueError("unsupported or malformed current run certificate")
    complete = certificate["status"] == "complete"
    if complete and certificate["pending"]:
        raise ValueError("complete certificate has a pending frontier")
    return canonical_output(certificate.get("output"), n=_integer(certificate.get("n"), positive=True), complete=complete)


def canonical_bytes(normal_form: Any) -> bytes:
    """Validate the envelope, re-normalize its sum, and encode deterministically."""
    if (type(normal_form) is not dict or set(normal_form) != {
            "schema", "coefficient_ring", "n", "completion", "terms"}
            or normal_form["schema"] != SCHEMA or normal_form["coefficient_ring"] != RING
            or normal_form["completion"] not in ("complete", "partial")):
        raise ValueError("unsupported output normal form")
    normalized = canonical_output(normal_form["terms"], n=normal_form["n"],
                                  complete=normal_form["completion"] == "complete")
    return _json(normalized).encode("ascii")


def output_fingerprint(normal_form: Any) -> str:
    """SHA256 lookup key of the exact canonical payload, including completion."""
    return hashlib.sha256(canonical_bytes(normal_form)).hexdigest()


def engine_output(normal_form: Any) -> list[dict[str, Any]]:
    """Decode the canonical sum to engine integer pairs, preserving all fields.

    This returns only terminal terms, not a certificate or checkpoint. Callers
    needing a complete result must also require completion == 'complete'.
    """
    normalized = json.loads(canonical_bytes(normal_form))
    for term in normalized["terms"]:
        term["coefficient"] = [[_exact(q), _exact(c)] for q, c in term["coefficient"]]
    return normalized["terms"]
