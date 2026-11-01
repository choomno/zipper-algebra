"""Explicit ideal cofactors over arbitrary declared coefficient fields.

Graph expansion retains noncommuting left/right contexts and collects the final
terms against the original ordered inputs. Finite-basis reduction returns both a
representative and an exact cofactor for their difference; no Gröbner completion
or anticipated coefficient support is required.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

from ..algebra import AlgebraElement, AlgebraPresentation, CofactorIdentity
from ..verification.derivations import (DerivationGraph, DerivationTerm, VerificationLimit,
                                       _Budget, _positive_limit)
from ..verification.finite import WordBasisCertificate, verify_word_basis
from ..verification.working_field import WorkingField
from ._contexts import _validate_graph

__all__ = ["CofactorReduction", "expand_derivation", "cofactor_representation", "verify_cofactor"]


def _parents(presentation, identity):
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(identity, CofactorIdentity):
        raise TypeError("cofactor verification requires an AlgebraPresentation and CofactorIdentity")
    if identity.presentation_hash != presentation.content_hash():
        raise ValueError("cofactor presentation hash mismatch")
    if identity.target.algebra != presentation.algebra:
        raise TypeError("cofactor and presentation have different free algebras")


def _add(destination, key, value, zero):
    result = destination.get(key, zero) + value
    if result == 0:
        destination.pop(key, None)
    else:
        destination[key] = result


def verify_cofactor(presentation: AlgebraPresentation, identity: CofactorIdentity, *,
                    max_operations=10000000) -> bool:
    """Check an explicit identity by exact normalized free-polynomial arithmetic.

    Arbitrary polynomial left and right factors are supported. The identity is
    interpreted over the declared field only. Malformed parents/references raise;
    a well-formed false identity returns False. The work bound counts scalar
    conversions and contextual coefficient placements, not RAM or elapsed time.
    """
    _parents(presentation, identity)
    budget = _Budget(max_operations)
    arithmetic = WorkingField(presentation.algebra.field)

    def native(polynomial):
        budget.consume(len(polynomial.terms))
        return tuple((word, arithmetic.from_element(coefficient))
                     for word, coefficient in polynomial.terms)

    expected = dict(native(identity.target))
    inputs = [native(polynomial) for polynomial in presentation.relations]
    result: dict = {}
    for left, index, right in identity.terms:
        budget.consume()
        if index >= len(inputs):
            raise ValueError("cofactor relation index is out of range")
        lhs, rhs = native(left), native(right)
        for u, a in lhs:
            for word, b in inputs[index]:
                for v, c in rhs:
                    budget.consume()
                    _add(result, u + word + v, a*b*c, arithmetic.zero)
    return result.keys() == expected.keys() and all(result[word] - expected[word] == 0
                                                    for word in result)


class _Expansion:
    def __init__(self, presentation, graph, max_operations, max_contexts, max_nodes):
        self.maximum = _positive_limit(max_contexts, "max_contexts")
        self.budget = _Budget(max_operations)
        self.order = _validate_graph(presentation, graph, self.budget, max_nodes)
        self.presentation, self.graph = presentation, graph
        self.arithmetic = WorkingField(presentation.algebra.field)
        self.peak = self.expanded_nodes = 0

    def expand(self, root):
        graph, arithmetic = self.graph, self.arithmetic
        if type(root) is not int or not 0 <= root < len(graph.nodes):
            raise ValueError("the selected derivation node is out of range")
        needed = set()
        pending = [root]
        dependencies = {}
        while pending:
            index = pending.pop()
            if index in needed:
                continue
            needed.add(index)
            references = {term.index for term in graph.derivations[index]
                          if term.kind == "node" and not term.coefficient.is_zero()}
            dependencies[index] = references
            pending.extend(references)
        users = dict.fromkeys(needed, 0)
        for references in dependencies.values():
            for index in references:
                users[index] += 1
        expanded: dict = {}
        retained = 0
        for index in self.order:
            if index not in needed:
                continue
            output: dict = {}
            for term in graph.derivations[index]:
                coefficient = arithmetic.from_element(term.coefficient)
                if coefficient == 0:
                    continue
                children = {((), term.index, ()): arithmetic.one} if term.kind == "input" else expanded[term.index]
                for (left, relation, right), value in children.items():
                    self.budget.consume()
                    context = term.left + left, relation, right + term.right
                    _add(output, context, coefficient*value, arithmetic.zero)
                    self.peak = max(self.peak, retained + len(output))
                    if retained + len(output) > self.maximum:
                        raise VerificationLimit("cofactor expansion exceeds its context-storage bound")
            expanded[index] = output
            retained += len(output)
            self.expanded_nodes += 1
            for reference in dependencies[index]:
                users[reference] -= 1
                if users[reference] == 0:
                    retained -= len(expanded.pop(reference))
        return expanded[root]

    def identity(self, root, maximum):
        contexts = self.expand(root)
        algebra = self.presentation.algebra
        terms = tuple((algebra.element({left: self.arithmetic.to_element(value)}), relation,
                       algebra.element({right: 1}))
                      for (left, relation, right), value in sorted(contexts.items(),
                          key=lambda item: (len(item[0][0]) + len(item[0][2]), item[0])))
        identity = CofactorIdentity(self.presentation.content_hash(), self.graph.nodes[root], terms)
        if not verify_cofactor(self.presentation, identity, max_operations=maximum):
            raise ValueError("expanded cofactor fails its exact original-input identity")
        return identity


def expand_derivation(presentation: AlgebraPresentation, graph: DerivationGraph, node: int, *,
                      backend="literal", max_operations=10000000, max_contexts=100000,
                      max_nodes=100000, max_matrix_entries=100000, max_primes=5,
                      max_degree=24, threads=1, start_prime=0, shift=None) -> CofactorIdentity:
    """Expand one graph node into a checked explicit original-input cofactor.

    All raw graph references and cycles are checked, including zero-coefficient
    and unrelated references. Only dependencies needed by the selected node are
    expanded. The final explicit identity is checked directly against the original
    presentation; no truth claim is made about unrelated node polynomials.
    ``backend="literal"`` propagates the existing path coefficients.
    ``backend="finiteflow"`` reconstructs one witness from structural input
    contexts over Q or a rational-function field. Both check the complete exact
    original-input identity. Native reconstruction options apply to FiniteFlow
    only; its matrix-entry bound does not bound native fill, time or RAM.

    ``max_contexts`` bounds peak retained intermediate context entries. Production
    and final identity verification each use their own operation budget.
    """
    _validate_backend(backend)
    identity, _ = _produce_identity(presentation, graph, node, backend=backend,
        max_operations=max_operations, max_contexts=max_contexts, max_nodes=max_nodes,
        native_options={"max_matrix_entries": max_matrix_entries, "max_primes": max_primes,
                        "max_degree": max_degree, "threads": threads,
                        "start_prime": start_prime, "shift": shift})
    return identity


def _validate_backend(backend):
    if backend not in ("literal", "finiteflow"):
        raise ValueError("cofactor backend must be 'literal' or 'finiteflow'")


def _produce_identity(presentation, graph, node, *, backend, max_operations,
                      max_contexts, max_nodes, native_options):
    if backend == "literal":
        expansion = _Expansion(presentation, graph, max_operations, max_contexts, max_nodes)
        identity = expansion.identity(node, max_operations)
        return identity, {"expanded_nodes": expansion.expanded_nodes,
                          "contexts": len(identity.terms), "peak_contexts": expansion.peak,
                          "expansion_operations": expansion.budget.used}
    from .reconstruction import reconstruct_derivation

    result = reconstruct_derivation(presentation, graph, node, max_operations=max_operations,
        max_contexts=max_contexts, max_nodes=max_nodes, **native_options)
    return result.identity, result.statistics


@dataclass(frozen=True, slots=True)
class CofactorReduction:
    """A representative and checked cofactor for target minus representative."""

    remainder: AlgebraElement
    identity: CofactorIdentity
    _statistics: tuple[tuple[str, int], ...]

    @property
    def statistics(self):
        return dict(self._statistics)


def cofactor_representation(presentation: AlgebraPresentation, certificate: WordBasisCertificate,
                            target: AlgebraElement, *, backend="literal", max_operations=10000000,
                            max_contexts=100000, max_nodes=100000, max_dimension=10000,
                            max_matrix_entries=100000, max_primes=5, max_degree=24,
                            threads=1, start_prime=0, shift=None) -> CofactorReduction:
    """Reduce an element and explicitly prove its difference from its representative.

    The finite-basis certificate is independently verified first. A zero remainder
    gives an explicit ideal-membership cofactor; a nonzero remainder represents a
    nonzero element of the certified quotient. Left generator closure identities
    recursively reduce suffixes, retaining word order on both sides of each input.
    ``backend="literal"`` expands coefficients along the graph. The optional
    ``"finiteflow"`` backend reconstructs a particular solution on its structural
    original-input contexts. Both produce the same kind of exactly checked
    identity; they may produce different witnesses. FiniteFlow supports Q and
    rational-function fields and requires the reconstruction dependencies.
    Its additional options bound matrix entries, degree/primes and worker count.

    Verification, reduction construction, cofactor production and final explicit checking
    have separate operation budgets. Limits do not bound process RAM or time.
    """
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(target, AlgebraElement):
        raise TypeError("cofactor reduction requires a presentation and free-polynomial target")
    if target.algebra != presentation.algebra:
        raise TypeError("target and presentation have different free algebras")
    _validate_backend(backend)
    _positive_limit(max_contexts, "max_contexts")
    algebra = verify_word_basis(presentation, certificate, max_operations=max_operations,
                                max_nodes=max_nodes, max_dimension=max_dimension)
    budget = _Budget(max_operations)
    arithmetic = WorkingField(presentation.algebra.field)
    nodes, rows = list(certificate.proof.nodes), list(certificate.proof.derivations)
    parent, field = presentation.algebra, presentation.algebra.field

    def polynomial(vector):
        budget.consume(len(vector))
        return parent.element({word: arithmetic.to_element(value) for word, value in vector.items()})

    def append(polynomial, terms):
        if len(nodes) >= max_nodes:
            raise VerificationLimit("cofactor reduction exceeds its proof-node bound")
        budget.consume(len(polynomial.terms) + len(terms) + 1)
        index = len(nodes)
        nodes.append(polynomial)
        rows.append(tuple(terms))
        return index

    remainder: dict = {}
    terms = []
    if not algebra.dimension:
        for word, coefficient in target.terms:
            budget.consume()
            terms.append(DerivationTerm(coefficient, word, "node", cast(int, certificate.unit_node), ()))
    else:
        actions = tuple(tuple(tuple((row, arithmetic.from_element(coefficient))
                                    for row, coefficient in column) for column in action.columns)
                        for action in certificate.actions)
        memo: dict[tuple[int, ...], tuple[dict[int, Any], int | None]] = {(): ({0: arithmetic.one}, None)}
        for word, coefficient in target.terms:
            for offset in range(len(word)-1, -1, -1):
                current = word[offset:]
                if current in memo:
                    continue
                vector, previous = memo[current[1:]]
                generator = current[0]-1
                result: dict = {}
                proof = [] if previous is None else [DerivationTerm(field.one, (generator+1,), "node", previous, ())]
                for column, value in vector.items():
                    budget.consume()
                    reference = certificate.action_nodes[generator][column]
                    if not nodes[reference].is_zero():
                        proof.append(DerivationTerm(arithmetic.to_element(value), (), "node", reference, ()))
                    for row, scalar in actions[generator][column]:
                        budget.consume()
                        _add(result, row, scalar*value, arithmetic.zero)
                difference = {current: arithmetic.one}
                for row, value in result.items():
                    _add(difference, certificate.basis_words[row], -value, arithmetic.zero)
                memo[current] = result, append(polynomial(difference), proof) if difference else None
            vector, word_reference = memo[word]
            scalar = arithmetic.from_element(coefficient)
            for row, value in vector.items():
                budget.consume()
                _add(remainder, certificate.basis_words[row], scalar*value, arithmetic.zero)
            if word_reference is not None:
                terms.append(DerivationTerm(coefficient, (), "node", word_reference, ()))
    representative = polynomial(remainder)
    difference = target - representative
    root = append(difference, terms)
    graph = DerivationGraph(presentation.content_hash(), tuple(nodes), tuple(rows))
    identity, counts = _produce_identity(presentation, graph, root, backend=backend,
        max_operations=max_operations, max_contexts=max_contexts, max_nodes=max_nodes,
        native_options={"max_matrix_entries": max_matrix_entries, "max_primes": max_primes,
                        "max_degree": max_degree, "threads": threads,
                        "start_prime": start_prime, "shift": shift})
    statistics = {"basis_dimension": algebra.dimension, "graph_nodes": len(nodes),
                  "construction_operations": budget.used, **counts}
    return CofactorReduction(representative, identity, tuple(sorted(statistics.items())))
