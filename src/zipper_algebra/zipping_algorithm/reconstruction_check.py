"""Independent exact assembly of a checked input tree and verified terminal leaves.

This checker follows the witness tree directly, without the producer's flattened
plan, route simulation, polynomial convolution, or reconstruction routine. Leaf
rewrite verification is the caller's responsibility.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import hashlib
import itertools
import json
from typing import Any, Mapping

from .input_normal_form import POLICY, check_input_normal_form
from .normal_form_reconstruction import completed_output, execution_identity
from .output_normal_form import canonical_bytes, canonical_output, engine_output


def _key(n: int, word: list[int]) -> str:
    value = {"policy": POLICY, "rank": n, "word": word}
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _shift(word: list[Any], offset: int) -> list[Any]:
    result = deepcopy(word)
    for i, token in enumerate(result):
        if isinstance(token, list):
            token[2] += offset
        else:
            result[i] = token + offset if token > 0 else token - offset
    return result


def _singleton(bottom: int, top: int) -> dict[str, Any]:
    record: dict[str, Any] = {name: [] for name in ("twist_powers", "top_route", "bottom_route",
                                  "canonical_top_route", "canonical_bottom_route")}
    record.update(original_bottom_ports=[bottom], original_top_ports=[top],
                  orientation="normal", canonical_parities={"upper": [], "lower_inverse": []})
    return {"size": 1, "component_word": [], "record": record}


def identity_leaf_request(form: Any, requests: Any) -> str | None:
    """Certify that reconstruction is literally the verified leaf itself.

    Check the original witness, literal request, and both endpoint maps. This
    shortcut has no output assembly to verify: its return value identifies the
    already independently verified complete result to return unchanged.
    """
    if not check_input_normal_form(form) or type(requests) is not list or len(requests) != 1:
        return None
    root = form["root"]
    if root["kind"] != "leaf":
        return None
    # The independent tree checker above requires both root endpoint maps
    # to be the ordered original labels 1,...,n.
    expected = {"key": _key(root["rank"], root["word"]),
                "n": root["rank"], "word": root["word"]}
    if json.dumps(requests[0], sort_keys=True) != json.dumps(expected, sort_keys=True):
        return None
    return expected["key"]


def checked_reconstruction(form: Any, leaves: Mapping[str, Any], *,
                           source_identity: Mapping[str, Any]) -> dict[str, Any]:
    """Derive the exact output independently of ``reconstruct_output``."""
    if not check_input_normal_form(form):
        raise ValueError("invalid input reconstruction witness")
    used: set[str] = set()

    def visit(node: dict[str, Any]) -> list[list[tuple[dict[int, int], list[dict[str, Any]]]]]:
        kind = node["kind"]
        if kind == "rank_two":
            return visit(node["child"])
        if kind == "split":
            return [part for edge in node["children"] for part in visit(edge["node"])]
        if kind == "identity":
            return [[({0: 1}, [_singleton(b, t)])]
                    for b, t in reversed(list(zip(node["bottom_labels"], node["top_labels"])))]
        if kind == "peel":
            head = _singleton(node["bottom_labels"][node["bottom_port"] - 1],
                              node["top_labels"][node["top_port"] - 1])
            return [[({0: 1}, [head])]] + visit(node["child"])
        key = _key(node["rank"], node["word"])
        used.add(key)
        result = leaves[key]
        if (result["n"] != node["rank"] or result["input_word"] != node["word"]
                or json.dumps(execution_identity(result), sort_keys=True)
                != json.dumps(source_identity, sort_keys=True)):
            raise ValueError("leaf identity mismatch")
        normal = completed_output(result)
        if normal["completion"] != "complete":
            raise ValueError("incomplete leaf")
        choices = []
        for term in engine_output(normal):
            components = deepcopy(term["components"])
            for component in components:
                for side in ("top", "bottom"):
                    field = "original_" + side + "_ports"
                    component["record"][field] = [node[side + "_labels"][p - 1]
                                                   for p in component["record"][field]]
            choices.append((dict(term["coefficient"]), components))
        return [choices]

    pieces = visit(form["root"])
    if set(leaves) != used:
        raise ValueError("unexpected leaf dependencies")
    n = form["root"]["rank"]
    terms = []
    for product in itertools.product(*pieces):
        # Balanced polynomial multiplication bounds intermediate support; this
        # checker does not reuse the producer's left-associated convolution.
        factors = [list(c.items()) for c, _ in product] or [[(0, 1)]]
        while len(factors) > 1:
            following = []
            for i in range(0, len(factors), 2):
                if i + 1 == len(factors):
                    following.append(factors[i])
                    continue
                sums: defaultdict[int, int] = defaultdict(int)
                for p, a in factors[i]:
                    for q, b in factors[i + 1]:
                        sums[p + q] += a * b
                following.append([(p, c) for p, c in sums.items() if c])
            factors = following
        scalar = dict(factors[0])
        components = deepcopy([c for _, group in product for c in group])
        remaining = {side: list(range(1, n + 1)) for side in ("top", "bottom")}
        left: list[Any] = []
        right: list[Any] = []
        offset = 0
        for component in components:
            component["rank"] = n - offset
            for side in ("top", "bottom"):
                ports = component["record"]["original_" + side + "_ports"]
                positions = [remaining[side].index(p) + 1 for p in ports]
                component["record"][side + "_ports"] = positions
                inverse = [-i for d, p in reversed(list(enumerate(positions, 1)))
                           for i in range(d, p)]
                route = [-i for i in reversed(inverse)] if side == "top" else inverse
                component[side + "_route"] = route
                remaining[side] = [p for p in remaining[side] if p not in ports]
            left += _shift(component["top_route"] + component["component_word"], offset)
            right = _shift(component["bottom_route"], offset) + right
            offset += component["size"]
        terms.append({"coefficient": [[p, c] for p, c in scalar.items()],
                      "word": left + right, "components": components})
    return canonical_output(terms, n=n, complete=True)


def check_reconstruction(form: Any, leaves: Mapping[str, Any], output: Any, *,
                         source_identity: Mapping[str, Any]) -> bool:
    """Compare exact canonical sums; reject malformed or mismatched witnesses."""
    try:
        expected = checked_reconstruction(form, leaves, source_identity=source_identity)
        return canonical_bytes(expected) == canonical_bytes(output)
    except (ValueError, KeyError, TypeError, IndexError):
        return False
