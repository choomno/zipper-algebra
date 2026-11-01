"""Produce certified regular actions from ideal-derived rewrite proposals.

Each word rewrite records a shared local identity for w minus its normal form.
Existing Gröbner derivations remain shared graph nodes; no original-input
cofactor expansion or full multiplication table is constructed.
"""
from __future__ import annotations

from collections import deque
from typing import Any

from .algebra import AlgebraPresentation
from .verification.derivations import (DerivationGraph, DerivationTerm, VerificationLimit,
                                      _Budget, _positive_limit)
from .verification.finite import GeneratorAction, WordBasisCertificate, verify_word_basis
from .verification.groebner import GroebnerCertificate
from .verification.working_field import WorkingField

__all__ = ["word_basis_from_rewrites"]


def _prune_certificate(presentation, certificate, budget):
    """Keep the raw dependency closure of the finite certificate's proof roots.

    Every original raw reference and cycle is checked before removing nodes,
    including unused and zero-coefficient references. Arithmetic of discarded
    proposal nodes is not asserted. The caller independently checks every
    retained identity and the complete finite-algebra statement afterward.
    """
    graph = certificate.proof
    size = len(graph.nodes)
    for row in graph.derivations:
        budget.consume()
        for term in row:
            budget.consume()
            bound = len(presentation.relations) if term.kind == "input" else size
            if term.index >= bound:
                raise ValueError("a raw derivation reference names no input or node polynomial")
    # A byte per node and a depth-first stack avoid a second full adjacency map.
    state = bytearray(size)
    for root in range(size):
        if state[root]:
            continue
        state[root] = 1
        stack = [(root, iter(graph.derivations[root]))]
        while stack:
            index, children = stack[-1]
            term = next(children, None)
            if term is None:
                state[index] = 2
                stack.pop()
                continue
            budget.consume()
            if term.kind == "node":
                if state[term.index] == 1:
                    raise ValueError("cyclic derivation graph")
                if not state[term.index]:
                    state[term.index] = 1
                    stack.append((term.index, iter(graph.derivations[term.index])))
    roots = ([certificate.unit_node] if not certificate.basis_words
             else [index for row in certificate.action_nodes for index in row])
    needed = bytearray(size)
    pending = []
    for root in roots:
        budget.consume()
        if root >= size:
            raise ValueError("a finite proof root references a missing derivation node")
        if not needed[root]:
            needed[root] = 1
            pending.append(root)
    while pending:
        for term in graph.derivations[pending.pop()]:
            budget.consume()
            if term.kind == "node" and not needed[term.index]:
                needed[term.index] = 1
                pending.append(term.index)
    if all(needed):
        return certificate
    selected = tuple(index for index, keep in enumerate(needed) if keep)
    remap = {index: position for position, index in enumerate(selected)}
    rows = tuple(tuple(DerivationTerm(term.coefficient, term.left, term.kind,
                                     remap[term.index], term.right)
                       if term.kind == "node" and remap[term.index] != term.index else term
                       for term in graph.derivations[index]) for index in selected)
    proof = DerivationGraph(graph.presentation_hash, tuple(graph.nodes[index] for index in selected), rows)
    return WordBasisCertificate(certificate.presentation_hash, certificate.basis_words,
                                certificate.actions, proof,
                                tuple(tuple(remap[index] for index in row) for row in certificate.action_nodes),
                                remap[certificate.unit_node] if certificate.unit_node is not None else None)


def word_basis_from_rewrites(presentation: AlgebraPresentation, certificate: GroebnerCertificate, *,
                             max_operations=10000000, max_nodes=100000,
                             max_words=100000) -> WordBasisCertificate:
    """Construct a finite word-basis witness and check it exactly once.

    The proposal need not have passed a Gröbner-completion check. Degree-lex
    orientation guarantees terminating reductions. Every irreducible word is
    enumerated, with finite enumeration and construction sharing the operation
    budget. Missing closure or a nonconfluent proposal cannot establish a
    dimension: the final independent checker requires all retained local ideal proofs,
    original relations on every matrix column, and cyclic independence.

    The finite proof retains only the raw dependency closure of the action
    identities (or the unit identity for the zero algebra). All proposal graph
    references and cycles are validated before pruning; no arithmetic assertion
    is made about discarded unused proposal nodes. The separate rewrite
    certificate is unchanged.

    ``max_nodes`` bounds the full temporary graph before pruning as well as
    the final retained graph. ``max_operations`` separately
    bounds construction and the independent finite-basis verification. Stored
    word lengths consume construction budget; ``max_words`` bounds their count.
    These bounds are not process time or RAM limits.
    """
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(certificate, GroebnerCertificate):
        raise TypeError("action production requires an AlgebraPresentation and GroebnerCertificate")
    if certificate.presentation_hash != presentation.content_hash():
        raise ValueError("rewrite proposal presentation hash mismatch")
    if any(polynomial.algebra != presentation.algebra for polynomial in certificate.basis):
        raise TypeError("rewrite proposal and presentation have different free algebras")
    _positive_limit(max_nodes, "max_nodes")
    _positive_limit(max_words, "max_words")
    budget = _Budget(max_operations)
    if len(certificate.basis) > max_nodes:
        raise VerificationLimit("action production exceeds its proof-node bound")
    algebra, field = presentation.algebra, presentation.algebra.field
    arithmetic = WorkingField(field)
    nodes = list(certificate.basis)
    derivations = []
    for row in certificate.derivations:
        converted = []
        for term in row:
            budget.consume()
            bound = len(presentation.relations) if term.kind == "input" else len(nodes)
            if term.index >= bound:
                raise ValueError("a raw rewrite reference names no original input or proposal node")
            converted.append(DerivationTerm(term.coefficient, term.left,
                              "input" if term.kind == "input" else "node", term.index, term.right))
        derivations.append(tuple(converted))
    leading = [max((word for word, _ in polynomial.terms), key=lambda word: (len(word), word))
               for polynomial in certificate.basis]
    selected = sorted(range(len(leading)), key=lambda index: (len(leading[index]), leading[index], index))
    rules: list[tuple[int, tuple[int, ...], Any, tuple]] = []
    for index in selected:
        word = leading[index]
        budget.consume(len(word) + 1)
        if any(any(word[position:position + len(divisor)] == divisor
                   for position in range(len(word) - len(divisor) + 1))
               for _, divisor, _, _ in rules):
            continue
        polynomial = certificate.basis[index]
        budget.consume(len(polynomial.terms))
        native = {monomial: arithmetic.from_element(coefficient)
                  for monomial, coefficient in polynomial.terms}
        inverse = arithmetic.one / native[word]
        tail = tuple((monomial, -coefficient * inverse) for monomial, coefficient in native.items()
                     if monomial != word)
        rules.append((index, word, inverse, tail))
    enumerated: list[tuple[int, ...]] = [] if rules and not rules[0][1] else [()]
    frontier = deque(enumerated)
    while frontier:
        word = frontier.popleft()
        for generator in range(1, len(algebra.generators) + 1):
            # Charge before allocation so a one-generator free algebra cannot
            # retain a quadratic number of word entries before max_words fires.
            budget.consume(len(word) + 1)
            candidate = word + (generator,)
            if any(any(candidate[position:position + len(divisor)] == divisor
                       for position in range(len(candidate) - len(divisor) + 1))
                   for _, divisor, _, _ in rules):
                continue
            enumerated.append(candidate)
            if len(enumerated) > max_words:
                raise VerificationLimit("finite irreducible-word enumeration did not finish")
            frontier.append(candidate)
    words = tuple(enumerated)

    def append(polynomial, terms):
        if len(nodes) >= max_nodes:
            raise VerificationLimit("action production exceeds its proof-node bound")
        budget.consume(len(polynomial.terms) + len(terms) + 1)
        index = len(nodes)
        nodes.append(polynomial)
        derivations.append(tuple(terms))
        return index

    def graph():
        return DerivationGraph(presentation.content_hash(), tuple(nodes), tuple(derivations))

    if not words:
        index, _, inverse, _ = rules[0]
        unit = append(algebra.one, (DerivationTerm(arithmetic.to_element(inverse), (), "node", index, ()),))
        result = WordBasisCertificate(presentation.content_hash(), (),
                                      tuple(GeneratorAction(field, ()) for _ in algebra.generators),
                                      graph(), tuple(() for _ in algebra.generators), unit)
    else:
        zero = append(algebra.zero, ())
        memo: dict = {}
        pending: dict = {}

        def accumulate(destination, word, coefficient):
            budget.consume()
            value = destination.get(word, arithmetic.zero) + coefficient
            if value == 0:
                destination.pop(word, None)
            else:
                destination[word] = value

        def normal_form(word):
            stack = [(word, False)]
            while stack:
                current, ready = stack.pop()
                if current in memo:
                    continue
                if ready:
                    index, left, right, inverse, children = pending.pop(current)
                    normal: dict = {}
                    terms = [DerivationTerm(arithmetic.to_element(inverse), left, "node", index, right)]
                    for child, coefficient in children:
                        child_normal, child_node = memo[child]
                        for monomial, value in child_normal.items():
                            accumulate(normal, monomial, coefficient * value)
                        if child_node != zero:
                            terms.append(DerivationTerm(arithmetic.to_element(coefficient), (),
                                                        "node", child_node, ()))
                    difference = {current: arithmetic.one}
                    for monomial, coefficient in normal.items():
                        accumulate(difference, monomial, -coefficient)
                    polynomial = algebra.element({monomial: arithmetic.to_element(coefficient)
                                                  for monomial, coefficient in difference.items()})
                    memo[current] = normal, append(polynomial, terms)
                    continue
                match = None
                for index, leading, inverse, tail in rules:
                    budget.consume()
                    position = next((position for position in range(len(current) - len(leading) + 1)
                                     if current[position:position + len(leading)] == leading), None)
                    if position is not None:
                        left, right = current[:position], current[position + len(leading):]
                        match = index, left, right, inverse, tuple((left + monomial + right, coefficient)
                                                                 for monomial, coefficient in tail)
                        break
                if match is None:
                    memo[current] = {current: arithmetic.one}, zero
                else:
                    pending[current] = match
                    stack.append((current, True))
                    stack.extend((child, False) for child, _ in match[-1] if child not in memo)
            return memo[word]

        rows = {word: row for row, word in enumerate(words)}
        actions = []
        action_nodes = []
        for generator in range(1, len(algebra.generators) + 1):
            columns, identities = [], []
            for word in words:
                normal, node = normal_form((generator,) + word)
                columns.append(tuple((rows[monomial], arithmetic.to_element(coefficient))
                                     for monomial, coefficient in normal.items()))
                identities.append(node)
            actions.append(GeneratorAction(field, tuple(columns)))
            action_nodes.append(tuple(identities))
        result = WordBasisCertificate(presentation.content_hash(), words, tuple(actions), graph(),
                                      tuple(action_nodes))
    result = _prune_certificate(presentation, result, budget)
    verify_word_basis(presentation, result, max_operations=max_operations,
                      max_nodes=max_nodes, max_dimension=max_words)
    return result
