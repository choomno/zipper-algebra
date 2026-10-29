"""Verify complete noncommutative Groebner bases over explicit exact fields.

The completion criterion is Bergman's Diamond Lemma, Theorem 1.2:
G. M. Bergman, Advances in Mathematics 29 (1978), 178--218,
https://doi.org/10.1016/0001-8708(78)90010-5 . We check both ideal inclusions
and every inclusion and proper overlap ambiguity under degree-lex order.
Derivations are checked as an acyclic graph of local identities, without
expanding shared subproofs into large cofactor lists.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from graphlib import CycleError, TopologicalSorter
import re

from ..algebra import AlgebraElement, AlgebraPresentation, FreeAlgebra, Word
from ..fields import FieldElement

__all__ = ["ProofTerm", "GroebnerCertificate", "VerifiedReduction",
           "VerificationLimit", "verify_groebner"]


class VerificationLimit(RuntimeError):
    """The verification budget ended before the requested computation finished."""


def _word(value: object) -> None:
    if (not isinstance(value, tuple)
            or any(type(index) is not int or index < 1 for index in value)):
        raise ValueError("proof contexts must be tuples of positive integer indices")


@dataclass(frozen=True, slots=True)
class ProofTerm:
    """One coefficient * left * referenced polynomial * right derivation term."""

    coefficient: FieldElement
    left: Word
    kind: str
    index: int
    right: Word

    def __post_init__(self) -> None:
        if not isinstance(self.coefficient, FieldElement):
            raise TypeError("proof coefficients must be exact field elements")
        _word(self.left)
        _word(self.right)
        if self.kind not in ("input", "basis"):
            raise ValueError("proof references must name input or basis polynomials")
        if type(self.index) is not int or self.index < 0:
            raise ValueError("proof reference indices must be nonnegative integers")

    def to_json(self) -> dict:
        return {"coefficient": self.coefficient.to_json(), "left": list(self.left),
                "kind": self.kind, "index": self.index, "right": list(self.right)}

    @classmethod
    def from_json(cls, value: object) -> ProofTerm:
        if (not isinstance(value, dict)
                or set(value) != {"coefficient", "left", "kind", "index", "right"}
                or not isinstance(value["left"], list)
                or not isinstance(value["right"], list)):
            raise ValueError("invalid Groebner derivation term")
        return cls(FieldElement.from_json(value["coefficient"]), tuple(value["left"]),
                   value["kind"], value["index"], tuple(value["right"]))


@dataclass(frozen=True, slots=True)
class GroebnerCertificate:
    """An ordered basis and local input/basis derivations bound to a presentation.

    A basis reference may point forward, but its graph must be acyclic. Raw
    references are validated even when their coefficient is zero.
    """

    presentation_hash: str
    basis: tuple[AlgebraElement, ...]
    derivations: tuple[tuple[ProofTerm, ...], ...]

    def __post_init__(self) -> None:
        if (not isinstance(self.presentation_hash, str)
                or re.fullmatch(r"[0-9a-f]{64}", self.presentation_hash) is None):
            raise ValueError("a Groebner certificate needs the full presentation SHA-256")
        if (not isinstance(self.basis, tuple)
                or any(not isinstance(polynomial, AlgebraElement) for polynomial in self.basis)):
            raise TypeError("the candidate basis must be a tuple of free polynomials")
        if any(polynomial.is_zero() for polynomial in self.basis):
            raise ValueError("candidate basis polynomials must be nonzero")
        if (not isinstance(self.derivations, tuple)
                or any(not isinstance(row, tuple)
                       or any(not isinstance(term, ProofTerm) for term in row)
                       for row in self.derivations)):
            raise TypeError("derivations must be a tuple of tuples of ProofTerm values")
        if len(self.basis) != len(self.derivations):
            raise ValueError("each candidate basis polynomial needs one derivation")
        if self.basis:
            algebra = self.basis[0].algebra
            if any(polynomial.algebra != algebra for polynomial in self.basis):
                raise TypeError("candidate basis polynomials have different free algebras")
            for row in self.derivations:
                for term in row:
                    if term.coefficient.field != algebra.field:
                        raise TypeError("a proof coefficient has the wrong field")
                    if any(index > len(algebra.generators) for index in term.left + term.right):
                        raise ValueError("a proof context uses an undeclared generator")

    def to_json(self) -> dict:
        return {"schema": "zipper-groebner-certificate-v1",
                "presentation_hash": self.presentation_hash,
                "basis": [polynomial.to_json() for polynomial in self.basis],
                "derivations": [[term.to_json() for term in row] for row in self.derivations]}

    @classmethod
    def from_json(cls, value: object) -> GroebnerCertificate:
        if (not isinstance(value, dict)
                or set(value) != {"schema", "presentation_hash", "basis", "derivations"}
                or value["schema"] != "zipper-groebner-certificate-v1"
                or not isinstance(value["basis"], list)
                or not isinstance(value["derivations"], list)
                or any(not isinstance(row, list) for row in value["derivations"])):
            raise ValueError("invalid Groebner certificate data")
        return cls(value["presentation_hash"],
                   tuple(AlgebraElement.from_json(item) for item in value["basis"]),
                   tuple(tuple(ProofTerm.from_json(term) for term in row)
                         for row in value["derivations"]))


def _order(word: Word) -> tuple:
    return len(word), word


def _occurrences(word: Word, divisor: Word):
    for position in range(len(word) - len(divisor) + 1):
        if word[position:position + len(divisor)] == divisor:
            yield position


def _accumulate(polynomial: dict, word: Word, coefficient) -> None:
    value = polynomial.get(word, 0) + coefficient
    if value == 0:
        polynomial.pop(word, None)
    else:
        polynomial[word] = value


class _Reducer:
    """Private sparse arithmetic; Q without parameters uses Fraction directly."""

    def __init__(self, algebra: FreeAlgebra, basis: tuple[AlgebraElement, ...],
                 maximum: int, prune: bool) -> None:
        self.algebra = algebra
        self.rational = not algebra.field.parameters and not algebra.field.characteristic
        self.one = Fraction(1) if self.rational else algebra.field.one
        self.maximum = maximum
        self.work = 0
        self.reductions = 0
        self.memo: dict = {}
        original = [self.native(polynomial) for polynomial in basis]
        leading = [max(polynomial, key=_order) for polynomial in original]
        selected = sorted(range(len(basis)), key=lambda i: (*_order(leading[i]), i))
        if prune:
            kept: list[int] = []
            for index in selected:
                if not any(next(_occurrences(leading[index], leading[other]), None) is not None
                           for other in kept):
                    kept.append(index)
            selected = kept
        # Leading words are computed, never accepted as caller assertions.
        # Consequently every tail word is strictly smaller in degree-lex.
        self.rules = []
        self.monic = []
        for index in selected:
            polynomial = original[index]
            word = leading[index]
            inverse = self.one / polynomial[word]
            normalized = {monomial: coefficient * inverse
                          for monomial, coefficient in polynomial.items()}
            self.monic.append(normalized)
            self.rules.append((word, tuple((monomial, -coefficient)
                                           for monomial, coefficient in normalized.items()
                                           if monomial != word)))

    def scalar(self, coefficient: FieldElement):
        if self.rational:
            return Fraction(coefficient.numerator[0][1]) if coefficient.numerator else Fraction(0)
        return coefficient

    def native(self, polynomial: AlgebraElement) -> dict:
        return {word: self.scalar(coefficient) for word, coefficient in polynomial.terms}

    def tick(self, *, reduction=False) -> None:
        if self.work >= self.maximum:
            raise VerificationLimit("Groebner verification exhausted its reduction/ambiguity budget")
        self.work += 1
        self.reductions += int(reduction)

    def match(self, word: Word):
        for leading, tail in self.rules:
            position = next(_occurrences(word, leading), None)
            if position is not None:
                return word[:position], tail, word[position + len(leading):]
        return None

    def word_normal_form(self, word: Word) -> dict:
        if word in self.memo:
            return self.memo[word]
        pending: dict = {}
        stack = [(word, False)]
        while stack:
            current, ready = stack.pop()
            if current in self.memo:
                continue
            if ready:
                result: dict = {}
                for coefficient, child in pending.pop(current):
                    for monomial, value in self.memo[child].items():
                        _accumulate(result, monomial, coefficient * value)
                self.memo[current] = result
                continue
            split = self.match(current)
            if split is None:
                self.memo[current] = {current: self.one}
                continue
            self.tick(reduction=True)
            left, tail, right = split
            children = [(coefficient, left + monomial + right) for monomial, coefficient in tail]
            pending[current] = children
            stack.append((current, True))
            stack.extend((child, False) for _, child in children if child not in self.memo)
        return self.memo[word]

    def normal_form(self, polynomial: dict) -> dict:
        result: dict = {}
        for word, coefficient in polynomial.items():
            for monomial, value in self.word_normal_form(word).items():
                _accumulate(result, monomial, coefficient * value)
        return result

    def contextual(self, polynomial: dict, left: Word, right: Word) -> dict:
        return {left + word + right: coefficient for word, coefficient in polynomial.items()}


def _check_derivations(presentation: AlgebraPresentation, certificate: GroebnerCertificate,
                       reducer: _Reducer) -> None:
    inputs = [reducer.native(polynomial) for polynomial in presentation.relations]
    basis = [reducer.native(polynomial) for polynomial in certificate.basis]
    graph: dict[int, set[int]] = {}
    for index, row in enumerate(certificate.derivations):
        graph[index] = set()
        for term in row:
            bound = len(inputs) if term.kind == "input" else len(basis)
            if term.index >= bound:
                raise ValueError("a raw proof reference names no input or basis polynomial")
            if term.kind == "basis":
                graph[index].add(term.index)
    try:
        order = tuple(TopologicalSorter(graph).static_order())
    except CycleError as error:
        raise ValueError("cyclic Groebner derivation graph") from error
    for index in order:
        result: dict = {}
        for term in certificate.derivations[index]:
            referenced = inputs[term.index] if term.kind == "input" else basis[term.index]
            coefficient = reducer.scalar(term.coefficient)
            for word, value in referenced.items():
                _accumulate(result, term.left + word + term.right, coefficient * value)
        if result != basis[index]:
            raise ValueError(f"basis entry {index} fails its local input-ideal derivation")


def _ambiguity_checks(reducer: _Reducer) -> tuple[int, int]:
    inclusions = overlaps = 0
    for i, (leading_i, _) in enumerate(reducer.rules):
        for j, (leading_j, _) in enumerate(reducer.rules):
            for position in _occurrences(leading_i, leading_j):
                if i == j and position == 0:
                    continue
                inclusions += 1
                reducer.tick()
                difference = dict(reducer.monic[i])
                contextual = reducer.contextual(reducer.monic[j], leading_i[:position],
                                                leading_i[position + len(leading_j):])
                for word, coefficient in contextual.items():
                    _accumulate(difference, word, -coefficient)
                if reducer.normal_form(difference):
                    raise ValueError(f"unresolved inclusion ambiguity between rules {i} and {j}")
            for length in range(1, min(len(leading_i), len(leading_j))):
                if leading_i[-length:] != leading_j[:length]:
                    continue
                overlaps += 1
                reducer.tick()
                difference = reducer.contextual(reducer.monic[i], (), leading_j[length:])
                contextual = reducer.contextual(reducer.monic[j], leading_i[:-length], ())
                for word, coefficient in contextual.items():
                    _accumulate(difference, word, -coefficient)
                if reducer.normal_form(difference):
                    raise ValueError(f"unresolved overlap ambiguity between rules {i} and {j}")
    return inclusions, overlaps


def _standard_words(reducer: _Reducer, maximum: int) -> tuple[Word, ...] | None:
    if any(not leading for leading, _ in reducer.rules):
        return ()
    words: list[Word] = [()]
    frontier: deque[Word] = deque([()])
    while frontier:
        word = frontier.popleft()
        for index in range(1, len(reducer.algebra.generators) + 1):
            candidate = word + (index,)
            if reducer.match(candidate) is None:
                words.append(candidate)
                if len(words) > maximum:
                    return None
                frontier.append(candidate)
    return tuple(words)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class VerifiedReduction:
    """Immutable defining data with a verified complete degree-lex reduction.

    Later normal-form calls share the verifier's memo and remaining work budget.
    Exhaustion raises VerificationLimit; it does not revoke the proved completion.
    """

    algebra: FreeAlgebra
    basis: tuple[AlgebraElement, ...]
    standard_words: tuple[Word, ...] | None
    _reducer: _Reducer
    _counts: tuple
    __hash__ = None  # type: ignore[assignment]

    def __init__(self) -> None:
        raise TypeError("VerifiedReduction values are returned by verify_groebner")

    @property
    def dimension(self) -> int | None:
        return None if self.standard_words is None else len(self.standard_words)

    @property
    def report(self) -> dict:
        """A copy of verification counts; later reductions do not change it."""
        return dict(self._counts)

    def normal_form(self, element: AlgebraElement) -> AlgebraElement:
        if not isinstance(element, AlgebraElement) or element.algebra != self.algebra:
            raise TypeError("normal forms require an element of the verified free algebra")
        result = self._reducer.normal_form(self._reducer.native(element))
        return self.algebra.element(result)


def verify_groebner(presentation: AlgebraPresentation, certificate: GroebnerCertificate, *,
                    max_reductions: int = 2_000_000, max_words: int = 100_000,
                    prune: bool = False) -> VerifiedReduction:
    """Check ideal equality, all ambiguities, then enumerate standard words.

    The work budget counts one unit per uncached word rewrite and per tested
    ambiguity. The word-count bound applies only after completion has been
    proved: exceeding it returns dimension=None, never an infinity verdict.
    Malformed or false certificates raise TypeError/ValueError. A work-budget
    limit raises VerificationLimit rather than claiming failure of membership.
    """
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(certificate, GroebnerCertificate):
        raise TypeError("verification requires a presentation and GroebnerCertificate")
    if certificate.presentation_hash != presentation.content_hash():
        raise ValueError("Groebner certificate presentation hash mismatch")
    if any(polynomial.algebra != presentation.algebra for polynomial in certificate.basis):
        raise TypeError("the candidate basis and presentation have different free algebras")
    if (type(max_reductions) is not int or max_reductions < 1
            or type(max_words) is not int or max_words < 1 or type(prune) is not bool):
        raise ValueError("verification limits must be positive integers and prune must be Boolean")
    reducer = _Reducer(presentation.algebra, certificate.basis, max_reductions, prune)
    _check_derivations(presentation, certificate, reducer)
    for index, polynomial in enumerate(presentation.relations):
        if reducer.normal_form(reducer.native(polynomial)):
            raise ValueError(f"input relation {index} has a nonzero remainder against the candidate basis")
    # A verified nonzero scalar proves the unit ideal; no ambiguities remain.
    inclusions, overlaps = ((0, 0) if any(not word for word, _ in reducer.rules)
                            else _ambiguity_checks(reducer))
    words = _standard_words(reducer, max_words)
    counts = {"basis_complete": True, "basis_entries": len(certificate.basis),
              "rules": len(reducer.rules), "input_relations": len(presentation.relations),
              "inclusion_ambiguities": inclusions, "overlap_ambiguities": overlaps,
              "reductions": reducer.reductions, "work_units": reducer.work,
              "standard_words_complete": words is not None,
              "dimension": None if words is None else len(words)}
    verified = object.__new__(VerifiedReduction)
    for name, value in (("algebra", presentation.algebra), ("basis", certificate.basis),
                        ("standard_words", words), ("_reducer", reducer),
                        ("_counts", tuple(counts.items()))):
        object.__setattr__(verified, name, value)
    return verified
