"""Structural validation and original-input contexts for shared ideal proofs."""
from __future__ import annotations

from graphlib import CycleError, TopologicalSorter
from typing import TypeAlias

from ..algebra import AlgebraPresentation, Word
from ..verification.derivations import (DerivationGraph, VerificationLimit, _Budget,
                                       _positive_limit)

Context: TypeAlias = tuple[Word, int, Word]


def _validate_graph(presentation: AlgebraPresentation, graph: DerivationGraph,
                    budget: _Budget, max_nodes: int) -> tuple[int, ...]:
    """Check every raw reference and cycle, including zero and unrelated edges."""
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(graph, DerivationGraph):
        raise TypeError("cofactor expansion requires an AlgebraPresentation and DerivationGraph")
    _positive_limit(max_nodes, "max_nodes")
    if len(graph.nodes) > max_nodes:
        raise VerificationLimit("cofactor expansion exceeds its proof-node bound")
    if graph.presentation_hash != presentation.content_hash():
        raise ValueError("derivation graph presentation hash mismatch")
    if any(node.algebra != presentation.algebra for node in graph.nodes):
        raise TypeError("derivation graph and presentation have different free algebras")
    dependencies: dict[int, set[int]] = {}
    for index, row in enumerate(graph.derivations):
        budget.consume()
        references: set[int] = set()
        for term in row:
            budget.consume()
            bound = len(presentation.relations) if term.kind == "input" else len(graph.nodes)
            if not 0 <= term.index < bound:
                raise ValueError("a raw derivation reference names no input or node polynomial")
            if term.kind == "node":
                references.add(term.index)
        dependencies[index] = references
    try:
        return tuple(TopologicalSorter(dependencies).static_order())
    except CycleError as error:
        raise ValueError("cyclic derivation graph") from error


def _collect_contexts(presentation: AlgebraPresentation, graph: DerivationGraph,
                      node: int, order: tuple[int, ...], budget: _Budget,
                      max_contexts: int) -> tuple[tuple[Context, ...], dict[str, int]]:
    """Union context labels without multiplying coefficients along proof paths.

    The storage bound counts the sum of live completed context sets and the new
    set under construction, before releasing a child after its last consumer.
    The caller must first validate the complete raw graph with _validate_graph.
    """
    _positive_limit(max_contexts, "max_contexts")
    if type(node) is not int or not 0 <= node < len(graph.nodes):
        raise ValueError("the selected derivation node is out of range")
    before = budget.used
    needed: set[int] = set()
    dependencies: dict[int, set[int]] = {}
    pending = [node]
    while pending:
        index = pending.pop()
        budget.consume()
        if index in needed:
            continue
        needed.add(index)
        references: set[int] = set()
        for term in graph.derivations[index]:
            budget.consume()
            if term.kind == "node" and not term.coefficient.is_zero():
                references.add(term.index)
        dependencies[index] = references
        pending.extend(references)
    users = dict.fromkeys(needed, 0)
    for references in dependencies.values():
        for index in references:
            budget.consume()
            users[index] += 1
    expanded: dict[int, set[Context]] = {}
    retained = peak = expanded_nodes = 0
    for index in order:
        if index not in needed:
            continue
        output: set[Context] = set()
        for term in graph.derivations[index]:
            budget.consume()
            if term.coefficient.is_zero():
                continue
            if term.kind == "input":
                if presentation.relations[term.index].is_zero():
                    continue
                children: tuple[Context, ...] | set[Context] = (((), term.index, ()),)
            else:
                children = expanded[term.index]
            for left, relation, right in children:
                budget.consume(len(term.left) + len(left) + len(right) + len(term.right) + 1)
                context = term.left + left, relation, right + term.right
                if context not in output:
                    if retained + len(output) >= max_contexts:
                        raise VerificationLimit("cofactor reconstruction exceeds its context-storage bound")
                    output.add(context)
                    peak = max(peak, retained + len(output))
        expanded[index] = output
        retained += len(output)
        expanded_nodes += 1
        for reference in dependencies[index]:
            users[reference] -= 1
            if users[reference] == 0:
                retained -= len(expanded.pop(reference))
    contexts = expanded[node]
    budget.consume(len(contexts))
    result = tuple(sorted(contexts, key=lambda item: (len(item[0]) + len(item[2]), item)))
    return result, {"needed_nodes": len(needed), "expanded_nodes": expanded_nodes,
                    "contexts": len(result), "peak_contexts": peak,
                    "structural_operations": budget.used - before}
