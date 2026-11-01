"""Produce independently checked regular actions from complete Gröbner data.

Each word rewrite records a shared local identity for w minus its normal form.
Existing Gröbner derivations remain shared graph nodes; no original-input
cofactor expansion or full multiplication table is constructed.
"""
from __future__ import annotations

from .algebra import AlgebraPresentation
from .verification.derivations import (DerivationGraph, DerivationTerm, VerificationLimit,
                                      _Budget, _positive_limit)
from .verification.finite import GeneratorAction, WordBasisCertificate, verify_word_basis
from .verification.groebner import GroebnerCertificate, verify_groebner
from .verification.working_field import WorkingField

__all__ = ["word_basis_from_groebner"]


def word_basis_from_groebner(presentation: AlgebraPresentation, certificate: GroebnerCertificate, *,
                             max_operations=10000000, max_nodes=100000,
                             max_reductions=2000000, max_words=100000) -> WordBasisCertificate:
    """Return exact regular generator actions with shared spanning proofs.

    Completion is verified afresh. Failure to finish standard-word enumeration
    is inconclusive and raises VerificationLimit. ``max_nodes`` bounds the
    retained proof graph; ``max_operations`` separately bounds construction and
    independent word-basis verification. ``max_reductions`` bounds the initial
    completion check. These bounds are not process time or RAM limits.
    """
    _positive_limit(max_nodes, "max_nodes")
    budget = _Budget(max_operations)
    reduction = verify_groebner(presentation, certificate, max_reductions=max_reductions,
                                max_words=max_words, prune=True)
    words = reduction.standard_words
    if words is None:
        raise VerificationLimit("finite standard-word enumeration did not finish")
    algebra, field = presentation.algebra, presentation.algebra.field
    arithmetic = WorkingField(field)
    nodes = list(certificate.basis)
    derivations = [tuple(DerivationTerm(term.coefficient, term.left,
                                       "input" if term.kind == "input" else "node",
                                       term.index, term.right) for term in row)
                   for row in certificate.derivations]
    if len(nodes) > max_nodes:
        raise VerificationLimit("action production exceeds its proof-node bound")

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
        index, scalar = next((index, polynomial.terms[0][1])
                             for index, polynomial in enumerate(certificate.basis)
                             if len(polynomial.terms) == 1 and polynomial.terms[0][0] == ())
        unit = append(algebra.one, (DerivationTerm(field.one / scalar, (), "node", index, ()),))
        result = WordBasisCertificate(presentation.content_hash(), (),
                                      tuple(GeneratorAction(field, ()) for _ in algebra.generators),
                                      graph(), tuple(() for _ in algebra.generators), unit)
    else:
        zero = append(algebra.zero, ())
        rules = []
        for index, polynomial in enumerate(certificate.basis):
            native = {word: arithmetic.from_element(coefficient) for word, coefficient in polynomial.terms}
            leading = max(native, key=lambda word: (len(word), word))
            inverse = arithmetic.one / native[leading]
            tail = tuple((word, -coefficient * inverse) for word, coefficient in native.items()
                         if word != leading)
            rules.append((index, leading, inverse, tail))
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
