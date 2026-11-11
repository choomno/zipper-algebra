"""Adversarial checks of the experimental graph boundary and exact assembly."""
from copy import deepcopy
import json
import zlib

import pytest

from experiments.zipper_derivation_dag import (
    capture, check, digest, initial, order, terminal_sum, encoded, stable_binding,
)
from zipper_algebra.zipping_algorithm.native import _decode_coefficient_object
from zipper_algebra.zipping_algorithm.interpreter import Limits


LIMITS = Limits()
IDENTITY = {"catalog": "example", "producer": "p", "checker": "v"}


def transition(work):
    word = work["left"]
    successors = [[1], [2]] if word == [1, 2] else [[]] if word else []
    children = []
    for successor in successors:
        child = deepcopy(work)
        child.update(left=successor, depth=1)
        children.append(child)
    return {"children": children, "output": [] if word else [
        {"word": [], "components": [], "coefficient": [[-3, -(2**200)], [7, 3]]}],
        "unit_evidence": {"rule": "test"}}


def graph():
    return capture(3, [1, 2], LIMITS, IDENTITY, transition)


def verify(value, **kwargs):
    options = dict(n=3, word=[1, 2], limits=LIMITS, identity=IDENTITY, transition=transition)
    options.update(kwargs)
    return check(value, **options)


def test_shared_node_counts_both_paths_exactly():
    value = graph()
    assert len(value["nodes"]) == 4
    assert next(iter(value["output"].values())) == {"coefficient": [[-3, -(2**201)], [7, 6]]}
    assert verify(value)


@pytest.mark.parametrize("field,replacement", [
    ("schema", "wrong"), ("request", {"n": 3, "word": [2, 1]}),
    ("identity", {"catalog": "other"}), ("transition_limits", {}), ("root", "x"),
    ("output", {}),
])
def test_mutations_rejected(field, replacement):
    value = graph()
    value[field] = replacement
    with pytest.raises(ValueError):
        verify(value)


def test_transition_evidence_and_state_are_checked():
    value = graph()
    value["nodes"][value["root"]]["result"]["unit_evidence"] = {"rule": "wrong"}
    with pytest.raises(ValueError, match="transition"):
        verify(value)
    value = graph()
    value["nodes"][value["root"]]["work"]["rank"] = 4
    with pytest.raises(ValueError, match="identity"):
        verify(value)


def test_missing_unreachable_and_cyclic_graphs_rejected():
    value = graph()
    del value["nodes"][digest(initial(3, []))]
    with pytest.raises(ValueError, match="missing"):
        verify(value)
    value = graph()
    value["nodes"]["unreachable"] = value["nodes"][value["root"]]
    with pytest.raises(ValueError, match="unreachable"):
        verify(value)
    value = graph()
    value["nodes"][digest(initial(3, []))]["result"]["children"] = [value["root"]]
    with pytest.raises(ValueError, match="cyclic"):
        verify(value)


def test_limits_and_refusals_are_not_complete_certificates():
    with pytest.raises(ValueError, match="node limit"):
        capture(3, [1, 2], LIMITS, IDENTITY, transition, max_nodes=3)
    with pytest.raises(ValueError, match="node limit"):
        verify(graph(), max_nodes=3)
    with pytest.raises(ValueError, match="did not complete"):
        capture(3, [1, 2], LIMITS, IDENTITY, lambda _: {"error": "incomplete"})
    value = graph()
    with pytest.raises(ValueError, match="coefficient bit"):
        terminal_sum(value["nodes"], value["root"], max_terms=10, max_integer_bits=100)
    with pytest.raises(ValueError, match="terminal key"):
        terminal_sum(value["nodes"], value["root"], max_terms=0, max_integer_bits=1000)


def test_duplicate_edges_retain_multiplicity():
    nodes = {"r": {"result": {"children": ["t", "t"], "output": []}},
             "t": {"result": {"children": [], "output": [
                 {"word": [], "components": [], "coefficient": [[0, 1]]}]}}}
    assert order(nodes, "r") == ["r", "t"]
    assert list(terminal_sum(nodes, "r", max_terms=1, max_integer_bits=8).values()) == [{"coefficient": [[0, 2]]}]


@pytest.mark.parametrize("n,word", [(True, []), (0, []), (3, [0]), (3, [3]), (3, [True])])
def test_invalid_input(n, word):
    with pytest.raises(ValueError, match="invalid literal"):
        initial(n, word)


def test_huge_coefficient_roundtrip_uses_exact_decimal_codec():
    huge = 1 << 20000
    def huge_transition(work):
        return {"children": [], "output": [{"word": [], "components": [],
                "coefficient": [[-4, -huge], [8, huge + 1]]}], "unit_evidence": {}}
    limits = Limits(max_integer_bits=30000)
    value = capture(3, [1, 2], limits, IDENTITY, huge_transition)
    loaded = json.loads(zlib.decompress(zlib.compress(encoded(value))), object_hook=_decode_coefficient_object)
    assert loaded == value
    assert verify(loaded, limits=limits, transition=huge_transition)


def test_stable_binding_excludes_build_session_changes():
    first = dict(abi="v1", source_sha256={"a": "b" * 64}, binary_sha256="c" * 64, flags=["-O3"],
                 build_seconds=4.0, compiled_during_initialization=True, binary="/first/path")
    second = dict(first, build_seconds=0.0, compiled_during_initialization=False, binary="/new/path")
    assert stable_binding(first) == stable_binding(second)
    assert stable_binding(first) != stable_binding(dict(first, binary_sha256="d" * 64))
    with pytest.raises(ValueError, match="binding"):
        stable_binding({})


@pytest.mark.parametrize("replacement", [None, {}, {"result": {"children": [None]}},
                                        {"result": {"children": {}}}])
def test_malformed_nodes_rejected(replacement):
    value = graph()
    value["nodes"][value["root"]] = replacement
    with pytest.raises(ValueError, match="shape"):
        verify(value)


def test_matching_different_root_and_nonzero_depth_rejected():
    value = graph()
    other = initial(3, [2, 1])
    key = digest(other)
    value["nodes"][key] = {"work": other, "result": {"children": [], "output": [], "unit_evidence": {}}}
    value["root"] = key
    with pytest.raises(ValueError, match="root input"):
        verify(value)
    value = graph()
    old = digest(initial(3, []))
    node = value["nodes"].pop(old)
    node["work"]["depth"] = 8
    new = digest(node["work"])
    value["nodes"][new] = node
    for item in value["nodes"].values():
        item["result"]["children"] = [new if key == old else key for key in item["result"]["children"]]
    with pytest.raises(ValueError, match="node identity"):
        verify(value)


def test_path_multiplicity_guard_is_distinct_from_output_bound():
    nodes = {"r": {"result": {"children": ["t", "t"], "output": []}},
             "t": {"result": {"children": [], "output": []}}}
    with pytest.raises(ValueError, match="path multiplicity"):
        terminal_sum(nodes, "r", max_terms=1, max_integer_bits=1)


def test_semantics_and_well_shaped_output_mutations_rejected():
    value = graph()
    value["semantics"]["depth"] = "actual-depth"
    with pytest.raises(ValueError, match="binding"):
        verify(value)
    value = graph()
    next(iter(value["output"].values()))["coefficient"][0][1] += 1
    with pytest.raises(ValueError, match="terminal assembly"):
        verify(value)


@pytest.mark.parametrize("change", [
    {"abi": 1}, {"abi": ""}, {"binary_sha256": "x" * 64},
    {"source_sha256": "a" * 64}, {"source_sha256": {}},
    {"source_sha256": {"file": 123}}, {"wrapper_source_sha256": None},
    {"reference_source_sha256": {"file": "short"}},
])
def test_stable_binding_requires_authenticated_digest_types(change):
    value = dict(abi="v1", source_sha256={"file": "a" * 64}, binary_sha256="b" * 64)
    value.update(change)
    with pytest.raises(ValueError, match="binding"):
        stable_binding(value)
