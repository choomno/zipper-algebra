"""Certified input interface for the existing zipper procedure.

This is a normal form of a reconstruction request, not of a braid group element.
Only literal support separation, index translation, and certified initial
singleton peeling, and the proved two-strand exponent classification are used.
No cancellation at rank >=3 or connected-block commutation is performed. Primitive keys identify reusable local calculations; the ordered
reconstruction plan additionally preserves all original endpoint embeddings.
It identifies full outputs under the separation/peeling theorems, not traces,
resource-limited statuses, or runtime. Source/version binding remains the
responsibility of the result store using these request keys.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any, Sequence

from .word_coverage import decompose

SCHEMA = "zipper-current-input-v1"
POLICY = "literal-support-peeling-and-rank-two-v1"


def _word(word: Sequence[int], n: int) -> list[int]:
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    if not isinstance(word, (list, tuple)):
        raise ValueError("literal word must be a list or tuple")
    result = list(word)
    if any(type(a) is not int or not 0 < abs(a) < n for a in result):
        raise ValueError("invalid signed generator")
    return result


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _digest(value: Any) -> str:
    return hashlib.sha256(_json(value).encode("ascii")).hexdigest()


def _peel(word: list[int], n: int) -> tuple[bool, list[int], int]:
    """Producer: follow all original labels and literally delete bottom n."""
    labels = list(range(1, n + 1))
    backward = []
    over = True
    for a in reversed(word):
        i = abs(a) - 1
        left, right = labels[i:i + 2]
        if n in (left, right):
            over &= (left if a > 0 else right) == n
        else:
            projected = i + 1 - int(labels.index(n) < i)
            backward.append(projected if a > 0 else -projected)
        labels[i], labels[i + 1] = right, left
    return over, list(reversed(backward)), labels.index(n) + 1


def input_normal_form(word: Sequence[int], *, n: int) -> dict[str, Any]:
    """Build a deterministic checked reconstruction tree retaining the input."""
    original = _word(word, n)

    def build(w: list[int], rank: int, top: list[int], bottom: list[int]) -> dict[str, Any]:
        node: dict[str, Any] = {"rank": rank, "word": w, "top_labels": top,
                                "bottom_labels": bottom}
        if not w:
            return dict(node, kind="identity")
        if rank == 2:
            exponent = sum(w)
            canonical = [1 if exponent > 0 else -1] * abs(exponent)
            if w != canonical:
                return dict(node, kind="rank_two", exponent=exponent,
                            child=build(canonical, rank, list(top), list(bottom)))
        if {abs(a) for a in w} != set(range(1, rank)):
            blocks = {b.offset: b for b in decompose(w, n=rank, commute=False).blocks}
            children = []
            offset = 0
            while offset < rank:
                block = blocks.get(offset)
                size = block.n if block else 1
                child_word = list(block.word) if block else []
                positions = list(block.positions) if block else []
                children.append({"offset": offset, "positions": positions,
                                 "node": build(child_word, size, top[offset:offset + size],
                                               bottom[offset:offset + size])})
                offset += size
            return dict(node, kind="split", children=list(reversed(children)))
        over, core, port = _peel(w, rank)
        if over:
            return dict(node, kind="peel", top_port=port, bottom_port=rank,
                        child=build(core, rank - 1, top[:port - 1] + top[port:], bottom[:-1]))
        return dict(node, kind="leaf")

    result = {"schema": SCHEMA, "policy": POLICY,
              "root": build(original, n, list(range(1, n + 1)), list(range(1, n + 1)))}
    if not check_input_normal_form(result):
        raise AssertionError("invalid generated input reconstruction witness")
    return result


def _checked_peel(w: list[int], rank: int) -> tuple[bool, list[int], int]:
    """Verifier: only track the removed strand's current position."""
    position, under = rank, False
    core: list[int] = []
    for a in reversed(w):
        index = abs(a)
        if position in (index, index + 1):
            under |= (position == index) != (a > 0)
            position = index + 1 if position == index else index
        else:
            projected = index if index < position else index - 1
            core.insert(0, projected if a > 0 else -projected)
    return not under, core, position


def check_input_normal_form(form: Any) -> bool:
    """Check every witness without calling the constructor or producer peel."""
    try:
        if (type(form) is not dict or set(form) != {"schema", "policy", "root"}
                or form["schema"] != SCHEMA or form["policy"] != POLICY):
            return False
        root = form["root"]
        rank = root["rank"]
        _word(root["word"], rank)
        labels = list(range(1, rank + 1))

        def visit(node: Any, n: int, top: list[int], bottom: list[int]) -> bool:
            if type(node) is not dict:
                return False
            base = {"kind", "rank", "word", "top_labels", "bottom_labels"}
            if (type(node.get("rank")) is not int or node["rank"] != n
                    or node.get("top_labels") != top or node.get("bottom_labels") != bottom
                    or type(node.get("top_labels")) is not list
                    or type(node.get("bottom_labels")) is not list
                    or any(type(x) is not int for x in node["top_labels"] + node["bottom_labels"])):
                return False
            w = _word(node["word"], n)
            if type(node["word"]) is not list:
                return False
            if not w:
                return set(node) == base and node["kind"] == "identity"
            if n == 2:
                exponent = sum(w)
                canonical = [1 if exponent > 0 else -1] * abs(exponent)
                if w != canonical:
                    if (set(node) != base | {"exponent", "child"}
                            or node["kind"] != "rank_two"
                            or type(node["exponent"]) is not int
                            or node["exponent"] != exponent):
                        return False
                    child = node["child"]
                    return (type(child) is dict and child.get("word") == canonical
                            and visit(child, n, top, bottom))
            support = {abs(a) for a in w}
            if support != set(range(1, n)):
                if set(node) != base | {"children"} or node["kind"] != "split":
                    return False
                if type(node["children"]) is not list:
                    return False
                # Independent partition: a missing generator ends its strand block.
                cuts = [0] + [i for i in range(1, n) if i not in support] + [n]
                intervals = list(reversed(list(zip(cuts, cuts[1:]))))
                if len(intervals) != len(node["children"]):
                    return False
                for (lo, hi), edge in zip(intervals, node["children"]):
                    if (type(edge) is not dict or set(edge) != {"offset", "positions", "node"}
                            or type(edge["offset"]) is not int or edge["offset"] != lo):
                        return False
                    positions = [i for i, a in enumerate(w) if lo < abs(a) < hi]
                    if (type(edge["positions"]) is not list or edge["positions"] != positions
                            or any(type(x) is not int for x in edge["positions"])):
                        return False
                    child_word = [(abs(w[i]) - lo) * (1 if w[i] > 0 else -1) for i in positions]
                    child = edge["node"]
                    if type(child) is not dict or child.get("word") != child_word:
                        return False
                    if not visit(child, hi - lo, top[lo:hi], bottom[lo:hi]):
                        return False
                return True
            over, core, port = _checked_peel(w, n)
            if over:
                if (set(node) != base | {"child", "top_port", "bottom_port"}
                        or node["kind"] != "peel" or type(node["top_port"]) is not int
                        or type(node["bottom_port"]) is not int
                        or node["top_port"] != port or node["bottom_port"] != n):
                    return False
                child = node["child"]
                return (type(child) is dict and child.get("word") == core
                        and visit(child, n - 1, top[:port - 1] + top[port:], bottom[:-1]))
            return set(node) == base and node["kind"] == "leaf"

        return visit(root, rank, labels, labels)
    except (ValueError, TypeError, KeyError, IndexError, RecursionError):
        return False


def canonical_input_json(form: dict[str, Any]) -> str:
    """Strict canonical serialization; reject unknown policy or malformed data."""
    if not check_input_normal_form(form):
        raise ValueError("invalid current-procedure input witness")
    return _json(form)


def _leaf_key(n: int, word: list[int]) -> str:
    return _digest({"policy": POLICY, "rank": n, "word": word})


def leaf_requests(form: dict[str, Any]) -> list[dict[str, Any]]:
    """Distinct primitive local computations, sorted by rank and literal word."""
    canonical_input_json(form)
    leaves: dict[str, dict[str, Any]] = {}

    def walk(node: dict[str, Any]) -> None:
        if node["kind"] == "leaf":
            key = _leaf_key(node["rank"], node["word"])
            leaves[key] = {"key": key, "n": node["rank"], "word": list(node["word"])}
        elif node["kind"] in ("peel", "rank_two"):
            walk(node["child"])
        elif node["kind"] == "split":
            for edge in node["children"]:
                walk(edge["node"])

    walk(form["root"])
    return sorted(leaves.values(), key=lambda item: (item["n"], item["word"]))


def reconstruction_plan(form: dict[str, Any]) -> dict[str, Any]:
    """Ordered singleton/leaf requests with original top and bottom maps.

    Equal plans give equal reconstructed full outputs under the current
    separation/peeling theorems. They need not have equal calculations. Plans
    are sufficient invariants, not complete classifiers of output equality.
    """
    canonical_input_json(form)
    atoms: list[dict[str, Any]] = []

    def walk(node: dict[str, Any]) -> None:
        kind = node["kind"]
        if kind == "identity":
            for b, t in reversed(list(zip(node["bottom_labels"], node["top_labels"]))):
                atoms.append({"kind": "singleton", "bottom": b, "top": t})
        elif kind == "peel":
            atoms.append({"kind": "singleton", "bottom": node["bottom_labels"][-1],
                          "top": node["top_labels"][node["top_port"] - 1]})
            walk(node["child"])
        elif kind == "rank_two":
            walk(node["child"])
        elif kind == "split":
            for edge in node["children"]:
                walk(edge["node"])
        else:
            atoms.append({"kind": "leaf", "key": _leaf_key(node["rank"], node["word"]),
                          "n": node["rank"], "word": list(node["word"]),
                          "bottom_labels": list(node["bottom_labels"]),
                          "top_labels": list(node["top_labels"])})

    walk(form["root"])
    return {"schema": SCHEMA, "policy": POLICY, "n": form["root"]["rank"], "atoms": atoms}


def reconstruction_key(form: dict[str, Any]) -> str:
    """Hash the canonical reconstruction plan, excluding literal provenance."""
    return _digest(reconstruction_plan(form))
