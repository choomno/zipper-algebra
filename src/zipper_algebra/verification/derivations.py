"""Check shared acyclic graphs of exact two-sided ideal identities."""
from __future__ import annotations

from dataclasses import dataclass
from graphlib import CycleError, TopologicalSorter
import re

from ..algebra import AlgebraElement, AlgebraPresentation, Word
from ..fields import FieldElement
from .working_field import WorkingField

__all__ = ["DerivationTerm", "DerivationGraph", "VerificationLimit", "verify_derivations"]


class VerificationLimit(RuntimeError):
    """An explicit verification or computation budget has been exhausted."""


def _positive_limit(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


class _Budget:
    def __init__(self, maximum):
        self.maximum = _positive_limit(maximum, "max_operations")
        self.used = 0

    def consume(self, count=1):
        if count > self.maximum - self.used:
            raise VerificationLimit("finite-algebra computation exhausted its operation budget")
        self.used += count


def _word(value):
    if (not isinstance(value, tuple)
            or any(type(index) is not int or index < 1 for index in value)):
        raise ValueError("derivation contexts must be tuples of positive integer indices")


@dataclass(frozen=True, slots=True)
class DerivationTerm:
    """A central coefficient times a contextual input or earlier proved node."""

    coefficient: FieldElement
    left: Word
    kind: str
    index: int
    right: Word

    def __post_init__(self):
        if not isinstance(self.coefficient, FieldElement):
            raise TypeError("derivation coefficients must be exact field elements")
        _word(self.left)
        _word(self.right)
        if self.kind not in ("input", "node"):
            raise ValueError("derivation references must name input or node polynomials")
        if type(self.index) is not int or self.index < 0:
            raise ValueError("derivation indices must be nonnegative integers")

    def to_json(self):
        return {"coefficient": self.coefficient.to_json(), "left": list(self.left),
                "kind": self.kind, "index": self.index, "right": list(self.right)}

    @classmethod
    def from_json(cls, value):
        if (not isinstance(value, dict)
                or set(value) != {"coefficient", "left", "kind", "index", "right"}
                or not isinstance(value["left"], list) or not isinstance(value["right"], list)):
            raise ValueError("invalid derivation term data")
        return cls(FieldElement.from_json(value["coefficient"]), tuple(value["left"]),
                   value["kind"], value["index"], tuple(value["right"]))


@dataclass(frozen=True, slots=True)
class DerivationGraph:
    """Local ideal identities bound to actual inputs; zero nodes are allowed.

    Forward node references are permitted. Acyclicity and every raw reference
    are checked at verification, including terms with zero coefficients.
    """

    presentation_hash: str
    nodes: tuple[AlgebraElement, ...]
    derivations: tuple[tuple[DerivationTerm, ...], ...]

    def __post_init__(self):
        if (not isinstance(self.presentation_hash, str)
                or re.fullmatch(r"[0-9a-f]{64}", self.presentation_hash) is None):
            raise ValueError("a derivation graph needs the full presentation SHA-256")
        if (not isinstance(self.nodes, tuple)
                or any(not isinstance(node, AlgebraElement) for node in self.nodes)):
            raise TypeError("derivation nodes must be an immutable tuple of free polynomials")
        if (not isinstance(self.derivations, tuple)
                or any(not isinstance(row, tuple)
                       or any(not isinstance(term, DerivationTerm) for term in row)
                       for row in self.derivations)):
            raise TypeError("derivations must be immutable tuples of DerivationTerm values")
        if len(self.nodes) != len(self.derivations):
            raise ValueError("each derivation node needs one local identity")
        if self.nodes:
            algebra = self.nodes[0].algebra
            if any(node.algebra != algebra for node in self.nodes):
                raise TypeError("derivation nodes have different free algebras")
            for row in self.derivations:
                for term in row:
                    if term.coefficient.field != algebra.field:
                        raise TypeError("a derivation coefficient has the wrong field")
                    if any(index > len(algebra.generators) for index in term.left + term.right):
                        raise ValueError("a derivation context uses an undeclared generator")

    def to_json(self):
        return {"schema": "zipper-derivation-graph-v1", "presentation_hash": self.presentation_hash,
                "nodes": [node.to_json() for node in self.nodes],
                "derivations": [[term.to_json() for term in row] for row in self.derivations]}

    @classmethod
    def from_json(cls, value):
        if (not isinstance(value, dict)
                or set(value) != {"schema", "presentation_hash", "nodes", "derivations"}
                or value["schema"] != "zipper-derivation-graph-v1"
                or not isinstance(value["nodes"], list)
                or not isinstance(value["derivations"], list)
                or any(not isinstance(row, list) for row in value["derivations"])):
            raise ValueError("invalid derivation graph data")
        return cls(value["presentation_hash"],
                   tuple(AlgebraElement.from_json(node) for node in value["nodes"]),
                   tuple(tuple(DerivationTerm.from_json(term) for term in row)
                         for row in value["derivations"]))


def _accumulate(polynomial, word, coefficient):
    value = polynomial.get(word, 0) + coefficient
    if value == 0:
        polynomial.pop(word, None)
    else:
        polynomial[word] = value


def _equal_sparse(left, right):
    # In particular, native fractions over finite fields need not have a
    # canonical nonzero denominator. Subtraction detects equality exactly.
    return left.keys() == right.keys() and all(left[key] - right[key] == 0 for key in left)


def _verify_graph(presentation, graph, scalar, budget=None):
    """Shared local checker; the caller selects an exact working arithmetic."""
    # Public entry points validate the presentation and graph types before
    # calling this private checker, including the legacy Groebner adapter.
    if graph.presentation_hash != presentation.content_hash():
        raise ValueError("derivation graph presentation hash mismatch")
    if any(node.algebra != presentation.algebra for node in graph.nodes):
        raise TypeError("derivation graph and presentation have different free algebras")
    inputs = [{word: scalar(coefficient) for word, coefficient in node.terms}
              for node in presentation.relations]
    nodes = [{word: scalar(coefficient) for word, coefficient in node.terms} for node in graph.nodes]
    dependencies: dict[int, set[int]] = {}
    for index, row in enumerate(graph.derivations):
        dependencies[index] = set()
        for term in row:
            if budget is not None:
                budget.consume()
            bound = len(inputs) if term.kind == "input" else len(nodes)
            if term.index >= bound:
                raise ValueError("a raw derivation reference names no input or node polynomial")
            if term.kind == "node":
                dependencies[index].add(term.index)
    try:
        order = tuple(TopologicalSorter(dependencies).static_order())
    except CycleError as error:
        raise ValueError("cyclic derivation graph") from error
    products = 0
    for index in order:
        result: dict = {}
        for term in graph.derivations[index]:
            source = inputs[term.index] if term.kind == "input" else nodes[term.index]
            coefficient = scalar(term.coefficient)
            for word, value in source.items():
                if budget is not None:
                    budget.consume()
                products += 1
                _accumulate(result, term.left + word + term.right, coefficient * value)
        if not _equal_sparse(result, nodes[index]):
            raise ValueError(f"derivation node {index} fails its local input-ideal identity")
    return {"nodes": len(nodes), "terms": sum(map(len, graph.derivations)),
            "coefficient_products": products}


def verify_derivations(presentation: AlgebraPresentation, graph: DerivationGraph, *,
                       max_operations: int = 2_000_000, max_nodes: int = 100_000) -> dict:
    """Verify ideal-membership identities, without any completion assertion."""
    _positive_limit(max_nodes, "max_nodes")
    if not isinstance(graph, DerivationGraph):
        raise TypeError("verification requires a DerivationGraph")
    if len(graph.nodes) > max_nodes:
        raise VerificationLimit("derivation graph exceeds its node bound")
    budget = _Budget(max_operations)
    if not isinstance(presentation, AlgebraPresentation):
        raise TypeError("verification requires an AlgebraPresentation")
    arithmetic = WorkingField(presentation.algebra.field)
    counts = _verify_graph(presentation, graph, arithmetic.from_element, budget)
    return {"verified": True, **counts, "operations": budget.used}
