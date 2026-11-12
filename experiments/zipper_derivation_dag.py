"""Bounded experimental transition DAGs; not the production census certificate.

The independent native diagnostic currently recomputes each entire transition,
including rule selection. This experiment measures that boundary; it is not a
small local-rule proof checker. Graph nodes preserve coefficients and contexts,
normalizing only diagnostic depth. There is no intrinsic-context quotient.
"""
from __future__ import annotations

import argparse
from collections import defaultdict, deque
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from time import perf_counter
from typing import Any, Callable
import zlib

from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import (
    NativeSession, _encode_coefficients, _decode_coefficient_object,
)
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier, _wire

Transition = Callable[[dict[str, Any]], dict[str, Any]]
SCHEMA = "zipper-experimental-transition-dag-v2"
SEMANTICS = {"root": "literal-singleton-v1", "depth": "diagnostic-zero",
             "limits": "per-transition-and-separate-graph-assembly-guards"}


def stable_binding(metadata: dict[str, Any]) -> dict[str, Any]:
    """Bind authenticated implementation/build bytes, not session timings/paths."""
    fields = ("abi", "source_sha256", "binary_sha256", "wrapper_source_sha256",
              "compiler", "compiler_version", "flags", "link_libraries", "libraries",
              "reference_source_sha256")
    required = {"abi", "source_sha256", "binary_sha256"}
    def sha(value: Any) -> bool:
        return (type(value) is str and len(value) == 64
                and all(c in "0123456789abcdef" for c in value))
    if (not required <= metadata.keys() or type(metadata.get("abi")) is not str
            or not metadata["abi"] or not sha(metadata.get("binary_sha256"))
            or type(metadata.get("source_sha256")) is not dict or not metadata["source_sha256"]
            or any(type(name) is not str or not name or not sha(value)
                   for name, value in metadata["source_sha256"].items())):
        raise ValueError("authenticated implementation binding required")
    if "wrapper_source_sha256" in metadata and not sha(metadata["wrapper_source_sha256"]):
        raise ValueError("invalid wrapper source binding")
    if "reference_source_sha256" in metadata and (
            type(metadata["reference_source_sha256"]) is not dict
            or not metadata["reference_source_sha256"]
            or any(type(name) is not str or not name or not sha(value)
                   for name, value in metadata["reference_source_sha256"].items())):
        raise ValueError("invalid reference source binding")
    return {key: deepcopy(metadata[key]) for key in fields if key in metadata}



def encoded(value: Any) -> bytes:
    return json.dumps(_encode_coefficients(value), sort_keys=True,
                      separators=(",", ":"), allow_nan=False).encode()


def digest(value: Any) -> str:
    return hashlib.sha256(encoded(value)).hexdigest()


def initial(n: int, word: list[int]) -> dict[str, Any]:
    if type(n) is not int or n < 1 or type(word) is not list or any(
            type(x) is not int or not 0 < abs(x) < n for x in word):
        raise ValueError("invalid literal input")
    return dict(rank=n, coefficient=[[0, 1]], contexts=[], world_phi=False,
                boundary_left=[], boundary_right=[], boundary_left_factors=[],
                boundary_right_factors=[], depth=0, bottom_labels=list(range(1, n + 1)),
                top_labels=list(range(1, n + 1)), phase="single", left=word[:],
                processed=[], position=n)


def normalized(work: dict[str, Any]) -> dict[str, Any]:
    value = deepcopy(work)
    value["depth"] = 0
    return value


def compact(result: dict[str, Any]) -> dict[str, Any]:
    if set(result) != {"children", "output", "unit_evidence"}:
        raise ValueError("transition did not complete or has invalid shape")
    return dict(children=[digest(normalized(c)) for c in result["children"]],
                output=deepcopy(result["output"]), unit_evidence=deepcopy(result["unit_evidence"]))


def order(nodes: dict[str, Any], root: str) -> list[str]:
    """Root-first topological order; reject missing, cyclic and unreachable nodes."""
    for key, node in nodes.items():
        if (type(key) is not str or type(node) is not dict
                or type(node.get("result")) is not dict
                or type(node["result"].get("children")) is not list
                or any(type(child) is not str for child in node["result"]["children"])):
            raise ValueError("invalid graph node shape")
    reached: set[str] = set()
    stack = [root]
    indegree: dict[str, int] = defaultdict(int)
    while stack:
        key = stack.pop()
        if key in reached:
            continue
        if key not in nodes:
            raise ValueError("missing node")
        reached.add(key)
        for child in nodes[key]["result"]["children"]:
            indegree[child] += 1
            stack.append(child)
    if reached != set(nodes):
        raise ValueError("unreachable node")
    ready = deque(k for k in nodes if indegree[k] == 0)
    result = []
    while ready:
        key = ready.popleft()
        result.append(key)
        for child in nodes[key]["result"]["children"]:
            indegree[child] -= 1
            if indegree[child] == 0:
                ready.append(child)
    if len(result) != len(nodes):
        raise ValueError("cyclic graph")
    return result


def terminal_sum(nodes: dict[str, Any], root: str, *, max_terms: int,
                 max_integer_bits: int) -> dict[str, dict[str, list[list[int]]]]:
    """Independent exact assembly by path multiplicities, without cached outputs."""
    multiplicity: dict[str, int] = defaultdict(int, {root: 1})
    sums: dict[str, dict[int, int]] = {}
    for key in order(nodes, root):
        weight = multiplicity[key]
        if weight.bit_length() > max_integer_bits:
            raise ValueError("path multiplicity limit exceeded")
        for child in nodes[key]["result"]["children"]:
            multiplicity[child] += weight
        for term in nodes[key]["result"]["output"]:
            if set(term) != {"coefficient", "word", "components"}:
                raise ValueError("invalid terminal term")
            identity = encoded({"word": term["word"], "components": term["components"]}).decode()
            target = sums.setdefault(identity, {})
            if len(sums) > max_terms:
                raise ValueError("terminal key limit exceeded")
            for q, c in term["coefficient"]:
                if type(q) is not int or type(c) is not int:
                    raise ValueError("noninteger coefficient")
                target[q] = target.get(q, 0) + weight * c
                if abs(target[q]).bit_length() > max_integer_bits:
                    raise ValueError("coefficient bit limit exceeded")
    return {key: {"coefficient": [[q, c] for q, c in sorted(value.items()) if c]}
            for key, value in sorted(sums.items()) if any(value.values())}


def capture(n: int, word: list[int], limits: Limits, identity: dict[str, Any],
            transition: Transition, *, max_nodes: int = 10000) -> dict[str, Any]:
    start = initial(n, word)
    root = digest(start)
    pending = [start]
    nodes: dict[str, Any] = {}
    while pending:
        work = pending.pop()
        key = digest(work)
        if key in nodes:
            if encoded(nodes[key]["work"]) != encoded(work):
                raise ValueError("state digest collision")
            continue
        if len(nodes) >= max_nodes:
            raise ValueError("node limit exceeded")
        result = transition(work)
        nodes[key] = {"work": work, "result": compact(result)}
        pending.extend(normalized(c) for c in result["children"])
    output = terminal_sum(nodes, root, max_terms=limits.max_terms,
                          max_integer_bits=limits.max_integer_bits)
    return dict(schema=SCHEMA, semantics=deepcopy(SEMANTICS), request={"n": n, "word": word[:]},
                transition_limits=asdict(limits),
                identity=deepcopy(identity), root=root, nodes=nodes, output=output)


def check(certificate: dict[str, Any], *, n: int, word: list[int], limits: Limits,
          identity: dict[str, Any], transition: Transition, max_nodes: int = 10000) -> bool:
    if set(certificate) != {"schema", "semantics", "request", "transition_limits", "identity", "root", "nodes", "output"}:
        raise ValueError("invalid certificate shape")
    if (certificate["schema"] != SCHEMA or encoded(certificate["request"]) != encoded({"n": n, "word": word})
            or encoded(certificate["semantics"]) != encoded(SEMANTICS)
            or encoded(certificate["transition_limits"]) != encoded(asdict(limits))
            or encoded(certificate["identity"]) != encoded(identity)):
        raise ValueError("certificate binding mismatch")
    nodes, root = certificate["nodes"], certificate["root"]
    if type(nodes) is not dict or not 0 < len(nodes) <= max_nodes:
        raise ValueError("node limit exceeded or invalid nodes")
    if root != digest(initial(n, word)):
        raise ValueError("root input mismatch")
    sequence = order(nodes, root)
    for key in sequence:
        node = nodes[key]
        if (set(node) != {"work", "result"} or type(node["work"]) is not dict
                or digest(node["work"]) != key or type(node["work"].get("depth")) is not int
                or node["work"].get("depth") != 0):
            raise ValueError("node identity mismatch")
    for key in sequence:
        node = nodes[key]
        if encoded(compact(transition(node["work"]))) != encoded(node["result"]):
            raise ValueError("independent transition mismatch")
    expected = terminal_sum(nodes, root, max_terms=limits.max_terms,
                            max_integer_bits=limits.max_integer_bits)
    if encoded(expected) != encoded(certificate["output"]):
        raise ValueError("terminal assembly mismatch")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cases", default='[[3,[1,2,1]],[3,[-1,-2,-1,-2]],[4,[1,3,2]]]')
    parser.add_argument("--max-nodes", type=int, default=10000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    limits = Limits(max_steps=250000, max_terms=20000, expansion_budget=1000000)
    rows = []
    with NativeSession() as producer, CppVerifier(producer_metadata=producer.metadata) as verifier:
        identity = {"producer": stable_binding(producer.metadata), "verifier": stable_binding(verifier.metadata),
                    "catalog_sha256": producer._catalog_sha256,
                    "reference_sources": producer._reference_source,
                    "experiment_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        for index, (n, word) in enumerate(json.loads(args.cases)):
            producer._request({"op": "clear_cache"})
            verifier.clear_cache()
            start = perf_counter()
            standard = producer.reduce_reused(word, n=n, limits=limits)
            produced = perf_counter() - start
            start = perf_counter()
            verifier.check_reused(standard)
            replay = perf_counter() - start
            verifier.clear_cache()
            start = perf_counter()
            graph = capture(n, word, limits, identity, lambda work: producer._request(
                {"op": "transition", "work": work, "limits": asdict(limits), "trace": False},
                decode_coefficients=True), max_nodes=args.max_nodes)
            captured = perf_counter() - start
            start = perf_counter()
            payload = zlib.compress(encoded(graph))
            target = args.output / f"case-{index}.json.zlib"
            target.write_bytes(payload)
            loaded = json.loads(zlib.decompress(target.read_bytes()), object_hook=_decode_coefficient_object)
            transfer = perf_counter() - start
            start = perf_counter()
            # One frozen-source guard around the whole graph, just as whole-run
            # replay guards its single request. Rehashing every source per node
            # would measure wrapper audit cost rather than the graph boundary.
            verifier._assert_source()
            check(loaded, n=n, word=word, limits=limits, identity=identity,
                  transition=lambda work: verifier._request({"op": "transition", "work": _wire(work),
                      "limits": asdict(limits), "trace": False}), max_nodes=args.max_nodes)
            verifier._assert_source()
            checked = perf_counter() - start
            synthetic = {"only": {"result": {"children": [], "output": standard["output"]}}}
            expected = terminal_sum(synthetic, "only", max_terms=limits.max_terms,
                                    max_integer_bits=limits.max_integer_bits)
            if encoded(expected) != encoded(loaded["output"]):
                raise ValueError("graph differs from completed native output")
            rows.append(dict(n=n, word=word, nodes=len(graph["nodes"]), terms=len(graph["output"]),
                             json_bytes=len(encoded(graph)), compressed_bytes=len(payload),
                             producer_seconds=produced, replay_seconds=replay,
                             capture_seconds=captured, store_load_seconds=transfer,
                             graph_check_seconds=checked, exact_output_equal=True))
            (args.output / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
