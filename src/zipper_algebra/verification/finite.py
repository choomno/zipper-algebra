"""Certified finite word bases, sparse generator actions and exact products.

Spanning is proved by ideal identities for every generator/basis-word product.
Independence is proved separately by a representation with a cyclic coordinate
basis. Neither a representation alone nor closure identities alone suffice.
"""
from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice
import re
from typing import Any, cast

from ..algebra import AlgebraElement, AlgebraPresentation, Word
from ..fields import FieldElement, RationalFunctionField
from .derivations import (DerivationGraph, VerificationLimit, _Budget,
                          _equal_sparse, _positive_limit, _verify_graph, _word)
from .working_field import WorkingField

__all__ = ["GeneratorAction", "WordBasisCertificate", "FiniteElement",
           "VerifiedFiniteAlgebra", "verify_word_basis", "VerificationLimit"]

_DEFAULT_OPERATIONS = 10_000_000


@dataclass(frozen=True, slots=True)
class GeneratorAction:
    """A square matrix in sparse columns; row and column indices start at zero."""

    field: RationalFunctionField
    columns: tuple[tuple[tuple[int, FieldElement], ...], ...]

    def __post_init__(self):
        if not isinstance(self.field, RationalFunctionField):
            raise TypeError("a generator action requires an exact coefficient field")
        if not isinstance(self.columns, tuple) or any(not isinstance(column, tuple) for column in self.columns):
            raise TypeError("matrix columns must be immutable tuples")
        size = len(self.columns)
        normalized = []
        for column in self.columns:
            seen = set()
            entries = []
            for entry in column:
                if (not isinstance(entry, tuple) or len(entry) != 2
                        or type(entry[0]) is not int or not 0 <= entry[0] < size):
                    raise ValueError("matrix entries need an in-range integer row and exact coefficient")
                row, coefficient = entry
                if row in seen:
                    raise ValueError("duplicate sparse matrix row")
                seen.add(row)
                if not isinstance(coefficient, FieldElement) or coefficient.field != self.field:
                    raise TypeError("matrix coefficients must belong to the declared field")
                if not coefficient.is_zero():
                    entries.append((row, coefficient))
            normalized.append(tuple(sorted(entries)))
        object.__setattr__(self, "columns", tuple(normalized))

    @property
    def dimension(self):
        return len(self.columns)

    def to_json(self):
        return {"schema": "zipper-generator-action-v1", "field": self.field.to_json(),
                "columns": [[[row, coefficient.to_json()] for row, coefficient in column]
                            for column in self.columns]}

    @classmethod
    def from_json(cls, value):
        if (not isinstance(value, dict) or set(value) != {"schema", "field", "columns"}
                or value["schema"] != "zipper-generator-action-v1"
                or not isinstance(value["columns"], list)
                or any(not isinstance(column, list) for column in value["columns"])):
            raise ValueError("invalid generator-action data")
        columns = []
        for column in value["columns"]:
            entries = []
            for entry in column:
                if not isinstance(entry, list) or len(entry) != 2:
                    raise ValueError("invalid sparse matrix entry")
                entries.append((entry[0], FieldElement.from_json(entry[1])))
            columns.append(tuple(entries))
        return cls(RationalFunctionField.from_json(value["field"]), tuple(columns))


@dataclass(frozen=True, slots=True)
class WordBasisCertificate:
    """A finite word basis proposed over a literal presentation.

    Actions follow the presentation's exact generator order. action_nodes[s][j]
    names a graph node equal to g_s*b_j minus the supplied column expression.
    The zero algebra instead has no words and a unit_node proving 1 in the ideal.
    """

    presentation_hash: str
    basis_words: tuple[Word, ...]
    actions: tuple[GeneratorAction, ...]
    proof: DerivationGraph
    action_nodes: tuple[tuple[int, ...], ...]
    unit_node: int | None = None

    def __post_init__(self):
        if (not isinstance(self.presentation_hash, str)
                or re.fullmatch(r"[0-9a-f]{64}", self.presentation_hash) is None):
            raise ValueError("a word-basis certificate needs the full presentation SHA-256")
        if not isinstance(self.basis_words, tuple):
            raise TypeError("basis words must be an immutable tuple")
        for word in self.basis_words:
            _word(word)
        if len(set(self.basis_words)) != len(self.basis_words):
            raise ValueError("basis words must be distinct")
        if self.basis_words:
            if self.basis_words[0] != () or self.unit_node is not None:
                raise ValueError("a nonempty word basis must start with the unit and have no unit_node")
        elif type(self.unit_node) is not int or self.unit_node < 0:
            raise ValueError("the zero algebra needs a nonnegative unit_node")
        if (not isinstance(self.actions, tuple)
                or any(not isinstance(action, GeneratorAction) for action in self.actions)):
            raise TypeError("generator actions must be an immutable tuple")
        if not isinstance(self.proof, DerivationGraph):
            raise TypeError("a word-basis certificate requires a DerivationGraph")
        if self.proof.presentation_hash != self.presentation_hash:
            raise ValueError("word-basis and derivation presentation hashes differ")
        if (not isinstance(self.action_nodes, tuple)
                or any(not isinstance(row, tuple)
                       or any(type(index) is not int or index < 0 for index in row)
                       for row in self.action_nodes)):
            raise ValueError("action node references must be immutable nonnegative integer rows")
        if (len(self.action_nodes) != len(self.actions)
                or any(len(row) != len(self.basis_words) for row in self.action_nodes)
                or any(action.dimension != len(self.basis_words) for action in self.actions)):
            raise ValueError("action matrices and proof-reference rows must have the basis dimension")

    def to_json(self):
        return {"schema": "zipper-word-basis-certificate-v1",
                "presentation_hash": self.presentation_hash,
                "basis_words": [list(word) for word in self.basis_words],
                "actions": [action.to_json() for action in self.actions],
                "proof": self.proof.to_json(), "action_nodes": [list(row) for row in self.action_nodes],
                "unit_node": self.unit_node}

    @classmethod
    def from_json(cls, value):
        if (not isinstance(value, dict)
                or set(value) != {"schema", "presentation_hash", "basis_words", "actions",
                                  "proof", "action_nodes", "unit_node"}
                or value["schema"] != "zipper-word-basis-certificate-v1"
                or not isinstance(value["basis_words"], list)
                or any(not isinstance(word, list) for word in value["basis_words"])
                or not isinstance(value["actions"], list)
                or not isinstance(value["action_nodes"], list)
                or any(not isinstance(row, list) for row in value["action_nodes"])):
            raise ValueError("invalid word-basis certificate data")
        return cls(value["presentation_hash"], tuple(tuple(word) for word in value["basis_words"]),
                   tuple(GeneratorAction.from_json(action) for action in value["actions"]),
                   DerivationGraph.from_json(value["proof"]),
                   tuple(tuple(row) for row in value["action_nodes"]), value["unit_node"])


def _add(vector, row, value):
    total = vector.get(row, 0) + value
    if total == 0:
        vector.pop(row, None)
    else:
        vector[row] = total


def _add_scaled(destination, source, scale, budget):
    if scale == 0:
        return
    for row, value in source.items():
        budget.consume()
        _add(destination, row, scale * value)


def _apply(action, vector, budget):
    result: dict = {}
    for column, coefficient in vector.items():
        for row, value in action[column]:
            budget.consume()
            _add(result, row, value * coefficient)
    return result


def _coefficient_key(value):
    if isinstance(value, FieldElement):
        return "field", value.numerator, value.denominator
    return value


class _ClearedActions:
    """Polynomial matrices A_g and nonzero denominators d_g with M_g=A_g/d_g.

    Denominators may be arbitrary parameter polynomials. Their least common
    multiples and exact quotients are computed once, outside matrix-vector
    arithmetic; nonmonic denominators over finite fields are allowed.
    """

    def __init__(self, arithmetic, actions, budget):
        self.field = arithmetic._native
        self.ring = self.field.ring
        self.denominators = []
        matrices = []
        for action in actions:
            denominator = self.ring.one
            for divisor in {value.denom for column in action for _, value in column}:
                budget.consume()
                denominator = denominator.lcm(divisor)
            self.denominators.append(denominator)
            columns = []
            for column in action:
                budget.consume(len(column))
                columns.append(tuple((row, value.numer * denominator.exquo(value.denom))
                                     for row, value in column))
            matrices.append(tuple(columns))
        self.matrices = tuple(matrices)

    def seed(self, vector, budget):
        denominator = self.ring.one
        for value in vector.values():
            budget.consume()
            denominator = denominator.lcm(value.denom)
        budget.consume(len(vector))
        return denominator, {row: value.numer * denominator.exquo(value.denom)
                             for row, value in vector.items()}


class _PolynomialDAG:
    """Normalized noncommutative Horner polynomials, sharing right factors.

    p = c + sum(g * p_g). Evaluating each p_g before multiplying permits its
    cancellations to happen before potentially dense generator actions.
    """

    def __init__(self, arithmetic, budget, cleared_actions=None):
        self.arithmetic = arithmetic
        self.budget = budget
        self.nodes: list[tuple[Any, tuple]] = [(arithmetic.zero, ())]
        self.indices = {(): 0}
        self.cleared_actions = cleared_actions
        self.cleared_nodes: list[tuple[Any, tuple]] = []
        self.denominators: list[Any] = []

    def _prepare(self, polynomial):
        if not polynomial:
            return self.arithmetic.zero, (), {}
        scale = polynomial[min(polynomial)]
        self.budget.consume(len(polynomial))
        normalized = {word: coefficient / scale for word, coefficient in polynomial.items()}
        key = tuple((word, _coefficient_key(coefficient)) for word, coefficient in sorted(normalized.items()))
        return scale, key, normalized

    def _frame(self, key, polynomial):
        groups: dict = {}
        for word, coefficient in polynomial.items():
            if word:
                groups.setdefault(word[0] - 1, {})[word[1:]] = coefficient
        children = []
        for generator, child in sorted(groups.items()):
            scale, child_key, normalized = self._prepare(child)
            children.append((generator, scale, child_key, normalized))
        return key, polynomial.get((), self.arithmetic.zero), children

    def add(self, polynomial):
        scale, key, normalized = self._prepare(polynomial)
        if key not in self.indices:
            stack = [self._frame(key, normalized)]
            while stack:
                current, constant, children = stack[-1]
                missing = next((child for child in children if child[2] not in self.indices), None)
                if missing is not None:
                    stack.append(self._frame(missing[2], missing[3]))
                    continue
                self.indices[current] = len(self.nodes)
                self.nodes.append((constant, tuple((generator, factor, self.indices[child_key])
                                                   for generator, factor, child_key, _ in children)))
                stack.pop()
        return scale, self.indices[key]

    def _compile_cleared(self):
        cleared = self.cleared_actions
        ring, field = cleared.ring, cleared.field
        for constant, children in self.nodes[len(self.cleared_nodes):]:
            denominator = constant.denom
            ratios = []
            for generator, factor, child in children:
                self.budget.consume()
                ratio = factor / field.new(cleared.denominators[generator] * self.denominators[child], ring.one)
                ratios.append((generator, ratio, child))
                denominator = denominator.lcm(ratio.denom)
            self.budget.consume(len(ratios) + 1)
            self.denominators.append(denominator)
            self.cleared_nodes.append((constant.numer * denominator.exquo(constant.denom),
                                      tuple((generator, ratio.numer * denominator.exquo(ratio.denom), child)
                                            for generator, ratio, child in ratios)))

    def _evaluate_node(self, index, seed, actions, budget, cache, nodes):
        stack = [index]
        while stack:
            current = stack[-1]
            if current in cache:
                stack.pop()
                continue
            constant, children = nodes[current]
            missing = next((child for _, _, child in children if child not in cache), None)
            if missing is not None:
                stack.append(missing)
                continue
            result: dict = {}
            _add_scaled(result, seed, constant, budget)
            for generator, factor, child in children:
                if cache[child]:
                    product = _apply(actions[generator], cache[child], budget)
                    _add_scaled(result, product, factor, budget)
            cache[current] = result
            stack.pop()
        return cache[index]

    def evaluate_node(self, index, seed, actions, budget, cache):
        if self.cleared_actions is None:
            return self._evaluate_node(index, seed, actions, budget, cache, self.nodes)
        self._compile_cleared()
        cleared = self.cleared_actions
        seed_denominator, polynomial_seed = cleared.seed(seed, budget)
        # If P_n=c+sum(f M_g P_child), choose D_n clearing c and
        # f/(d_g D_child). The cached polynomial vector is exactly
        # V_n=D_n P_n(M) U, where U=S*seed. Thus all cancellations in the
        # matrix loop use polynomial +/*, with no fraction GCD operations.
        # Division by D_n*S occurs only at the requested output boundary.
        vector = self._evaluate_node(index, polynomial_seed, cleared.matrices, budget, cache,
                                     self.cleared_nodes)
        denominator = self.denominators[index] * seed_denominator
        budget.consume(len(vector))
        return {row: cleared.field.new(value, denominator) for row, value in vector.items()}

    def evaluate(self, root, seed, actions, budget, cache):
        factor, index = root
        value = self.evaluate_node(index, seed, actions, budget, cache)
        result: dict = {}
        _add_scaled(result, value, factor, budget)
        return result


def _native_actions(actions, arithmetic, budget):
    result = []
    for action in actions:
        columns = []
        for column in action.columns:
            budget.consume(len(column))
            columns.append(tuple((row, arithmetic.from_element(coefficient)) for row, coefficient in column))
        result.append(tuple(columns))
    return tuple(result)


def _native_polynomial(element, arithmetic, budget):
    budget.consume(len(element.terms))
    return {word: arithmetic.from_element(coefficient) for word, coefficient in element.terms}


@dataclass(frozen=True, slots=True, init=False, eq=False)
class VerifiedFiniteAlgebra:
    """A verified word basis with sparse regular generator actions."""

    presentation: AlgebraPresentation
    basis_words: tuple[Word, ...]
    _actions: tuple[GeneratorAction, ...]
    _native_actions: tuple
    _cleared_actions: Any
    _arithmetic: WorkingField
    _counts: tuple
    __hash__ = None  # type: ignore[assignment]

    def __init__(self):
        raise TypeError("VerifiedFiniteAlgebra values are returned by verify_word_basis")

    @property
    def dimension(self):
        return len(self.basis_words)

    @property
    def field(self):
        return self.presentation.algebra.field

    @property
    def report(self):
        return dict(self._counts)

    def element(self, coordinates: Iterable[int | Fraction | FieldElement]):
        return FiniteElement(self, coordinates)

    def scalar(self, coefficient):
        value = self.field.scalar(coefficient)
        zero = self.field.zero
        return self.element((value,) + (zero,) * (self.dimension - 1) if self.dimension else ())

    @property
    def zero(self):
        return self.element((self.field.zero,) * self.dimension)

    @property
    def one(self):
        return self.scalar(1)

    def generator_action(self, name):
        return self._actions[self.presentation.algebra.generators.index(name)]

    def _native_vector(self, element):
        if not isinstance(element, FiniteElement) or element.algebra is not self:
            raise TypeError("finite elements belong to different verified algebra parents")
        return {row: self._arithmetic.from_element(value)
                for row, value in enumerate(element.coordinates) if not value.is_zero()}

    def _element_from_native(self, vector):
        zero = self.field.zero
        return self.element(tuple(self._arithmetic.to_element(vector[row]) if row in vector else zero
                                  for row in range(self.dimension)))

    def coordinates(self, polynomial: AlgebraElement, *, max_operations=_DEFAULT_OPERATIONS):
        if not isinstance(polynomial, AlgebraElement) or polynomial.algebra != self.presentation.algebra:
            raise TypeError("coordinates require a polynomial of the declared free algebra")
        budget = _Budget(max_operations)
        dag = _PolynomialDAG(self._arithmetic, budget, self._cleared_actions)
        root = dag.add(_native_polynomial(polynomial, self._arithmetic, budget))
        seed = {0: self._arithmetic.one} if self.dimension else {}
        vector = dag.evaluate(root, seed, self._native_actions, budget, {0: {}})
        return self._element_from_native(vector).coordinates

    def multiply(self, left, right, *, max_operations=_DEFAULT_OPERATIONS):
        left_vector, right_vector = self._native_vector(left), self._native_vector(right)
        budget = _Budget(max_operations)
        budget.consume(len(left_vector) + len(right_vector))
        dag = _PolynomialDAG(self._arithmetic, budget, self._cleared_actions)
        root = dag.add({self.basis_words[index]: value for index, value in left_vector.items()})
        result = dag.evaluate(root, right_vector, self._native_actions, budget, {0: {}})
        return self._element_from_native(result)

    def _linear(self, left, right, sign=1):
        a, b = self._native_vector(left), self._native_vector(right)
        budget = _Budget(_DEFAULT_OPERATIONS)
        _add_scaled(a, b, self._arithmetic.one * sign, budget)
        return self._element_from_native(a)


@dataclass(frozen=True, slots=True, init=False, eq=False)
class FiniteElement:
    """Exact coordinates in one immutable verified word-basis parent."""

    algebra: VerifiedFiniteAlgebra
    coordinates: tuple[FieldElement, ...]
    __hash__ = None  # type: ignore[assignment]

    def __init__(self, algebra, coordinates):
        if not isinstance(algebra, VerifiedFiniteAlgebra):
            raise TypeError("finite elements require a verified algebra parent")
        values = tuple(algebra.field.scalar(value) for value in islice(coordinates, algebra.dimension + 1))
        if len(values) != algebra.dimension:
            raise ValueError("coordinate count differs from the verified dimension")
        object.__setattr__(self, "algebra", algebra)
        object.__setattr__(self, "coordinates", values)

    def _lift(self, other):
        if isinstance(other, FiniteElement):
            if other.algebra is not self.algebra:
                raise TypeError("finite elements belong to different verified algebra parents")
            return other
        if isinstance(other, (int, Fraction, FieldElement)):
            return self.algebra.scalar(other)
        return NotImplemented

    def __add__(self, other):
        rhs = self._lift(other)
        return NotImplemented if rhs is NotImplemented else self.algebra._linear(self, rhs)

    __radd__ = __add__

    def __neg__(self):
        return self.algebra._linear(self.algebra.zero, self, -1)

    def __sub__(self, other):
        rhs = self._lift(other)
        return NotImplemented if rhs is NotImplemented else self.algebra._linear(self, rhs, -1)

    def __rsub__(self, other):
        lhs = self._lift(other)
        return NotImplemented if lhs is NotImplemented else self.algebra._linear(lhs, self, -1)

    def __mul__(self, other):
        rhs = self._lift(other)
        return NotImplemented if rhs is NotImplemented else self.algebra.multiply(self, rhs)

    def __rmul__(self, other):
        lhs = self._lift(other)
        return NotImplemented if lhs is NotImplemented else self.algebra.multiply(lhs, self)

    def __eq__(self, other):
        if not isinstance(other, FiniteElement):
            return NotImplemented
        return self.algebra is other.algebra and self.coordinates == other.coordinates


def verify_word_basis(presentation: AlgebraPresentation, certificate: WordBasisCertificate, *,
                      max_operations: int = _DEFAULT_OPERATIONS, max_dimension: int = 10_000,
                      max_nodes: int = 100_000) -> VerifiedFiniteAlgebra:
    """Verify spanning identities, all relation matrices, and a cyclic basis.

    Every nonzero sparse coefficient multiply/accumulate consumes work budget.
    Matrix relations use normalized polynomial Horner graphs, not expanded
    full matrix products. An exhausted bound raises VerificationLimit.
    """
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(certificate, WordBasisCertificate):
        raise TypeError("verification requires an AlgebraPresentation and WordBasisCertificate")
    if certificate.presentation_hash != presentation.content_hash():
        raise ValueError("word-basis certificate presentation hash mismatch")
    _positive_limit(max_dimension, "max_dimension")
    _positive_limit(max_nodes, "max_nodes")
    dimension = len(certificate.basis_words)
    if dimension > max_dimension or len(certificate.proof.nodes) > max_nodes:
        raise VerificationLimit("word-basis certificate exceeds its dimension or node bound")
    algebra = presentation.algebra
    if len(certificate.actions) != len(algebra.generators):
        raise ValueError("one action matrix is required for every declared generator")
    if any(action.field != algebra.field for action in certificate.actions):
        raise TypeError("action matrices and presentation have different coefficient fields")
    if any(index > len(algebra.generators) for word in certificate.basis_words for index in word):
        raise ValueError("a basis word uses an undeclared generator")
    budget = _Budget(max_operations)
    arithmetic = WorkingField(algebra.field)
    proof_counts = _verify_graph(presentation, certificate.proof, arithmetic.from_element, budget)
    if not dimension:
        unit_index = cast(int, certificate.unit_node)
        if (unit_index >= len(certificate.proof.nodes)
                or certificate.proof.nodes[unit_index] != algebra.one):
            raise ValueError("the zero algebra requires a node proving the unit belongs to the ideal")
    else:
        for generator, references in enumerate(certificate.action_nodes, 1):
            for column, node_index in enumerate(references):
                if node_index >= len(certificate.proof.nodes):
                    raise ValueError("an action identity references a missing derivation node")
                expected: dict = {}
                expected[(generator,) + certificate.basis_words[column]] = arithmetic.one
                for row, coefficient in certificate.actions[generator - 1].columns[column]:
                    budget.consume()
                    word = certificate.basis_words[row]
                    value = expected.get(word, 0) - arithmetic.from_element(coefficient)
                    if value == 0:
                        expected.pop(word, None)
                    else:
                        expected[word] = value
                actual = _native_polynomial(certificate.proof.nodes[node_index], arithmetic, budget)
                if not _equal_sparse(expected, actual):
                    raise ValueError("an action derivation does not prove its required generator-column identity")
    actions = _native_actions(certificate.actions, arithmetic, budget)
    cleared_actions = _ClearedActions(arithmetic, actions, budget) if algebra.field.parameters else None
    dag = _PolynomialDAG(arithmetic, budget, cleared_actions)
    relation_roots = [dag.add(_native_polynomial(relation, arithmetic, budget))
                      for relation in presentation.relations]
    cyclic_roots = [dag.add({word: arithmetic.one}) for word in certificate.basis_words]
    for column in range(dimension):
        seed = {column: arithmetic.one}
        cache: dict[int, dict] = {0: {}}
        for _, root in relation_roots:
            if dag.evaluate_node(root, seed, actions, budget, cache):
                raise ValueError("an original defining relation does not vanish on the supplied matrices")
        if column == 0:
            for index, root in enumerate(cyclic_roots):
                if not _equal_sparse(dag.evaluate(root, seed, actions, budget, cache), {index: arithmetic.one}):
                    raise ValueError("the proposed words do not give a cyclic coordinate basis")
    counts = {"verified": True, "dimension": dimension, "generators": len(actions),
              "action_identities": len(actions) * dimension, "derivation_nodes": proof_counts["nodes"],
              "derivation_terms": proof_counts["terms"],
              "relation_columns": len(presentation.relations) * dimension,
              "cyclic_words": dimension, "horner_nodes": len(dag.nodes), "operations": budget.used}
    result = object.__new__(VerifiedFiniteAlgebra)
    for name, value in (("presentation", presentation), ("basis_words", certificate.basis_words),
                        ("_actions", certificate.actions), ("_native_actions", actions),
                        ("_cleared_actions", cleared_actions),
                        ("_arithmetic", arithmetic), ("_counts", tuple(counts.items()))):
        object.__setattr__(result, name, value)
    return result
