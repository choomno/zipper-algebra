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


def word_basis_from_rewrites(presentation: AlgebraPresentation, certificate: GroebnerCertificate, *,
                             max_operations=10000000, max_nodes=100000,
                             max_words=100000) -> WordBasisCertificate:
    """Construct a finite word-basis witness and check it exactly once.

    The proposal need not have passed a Gröbner-completion check. Degree-lex
    orientation guarantees terminating reductions. Every irreducible word is
    enumerated, with finite enumeration and construction sharing the operation
    budget. Missing closure or a nonconfluent proposal cannot establish a
    dimension: the final independent checker requires all local ideal proofs,
    original relations on every matrix column, and cyclic independence.

    ``max_nodes`` bounds the retained proof graph. ``max_operations`` separately
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
    derivations = [tuple(DerivationTerm(term.coefficient, term.left,
                                       "input" if term.kind == "input" else "node",
                                       term.index, term.right) for term in row)
                   for row in certificate.derivations]
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
    verify_word_basis(presentation, result, max_operations=max_operations,
                      max_nodes=max_nodes, max_dimension=max_words)
    return result
