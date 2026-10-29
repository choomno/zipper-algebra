"""Domain-bound free algebras, literal presentations and exact cofactor identities.

This interface supports arbitrary parameter fields without widening the legacy
q,t certificate format. Equality of coefficients is exact; presentation identity
records the literal coefficient expressions, including their stored denominators.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
from types import NotImplementedType

from .fields import FieldElement, RationalFunctionField

__all__ = ["FreeAlgebra", "AlgebraElement", "AlgebraPresentation", "CofactorIdentity"]

Scalar = int | Fraction | FieldElement
Word = tuple[int, ...]


@dataclass(frozen=True, slots=True)
class FreeAlgebra:
    """A free associative algebra over an explicit field and ordered alphabet."""

    field: RationalFunctionField
    generators: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.field, RationalFunctionField):
            raise TypeError("a free algebra requires a RationalFunctionField")
        if (not isinstance(self.generators, tuple)
                or any(not isinstance(name, str)
                       or re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", name) is None
                       for name in self.generators)
                or len(set(self.generators)) != len(self.generators)):
            raise ValueError("generators must be a tuple of distinct ASCII identifiers")

    def element(self, terms: Mapping[Word, Scalar]) -> AlgebraElement:
        return AlgebraElement(self, terms)

    def scalar(self, value: Scalar) -> AlgebraElement:
        return self.element({(): value})

    @property
    def zero(self) -> AlgebraElement:
        return self.element({})

    @property
    def one(self) -> AlgebraElement:
        return self.scalar(1)

    def generator(self, name: str) -> AlgebraElement:
        """Return a named generator; word indices in stored data are one-based."""
        return self.element({(self.generators.index(name) + 1,): 1})

    def to_json(self) -> dict:
        return {"schema": "zipper-free-algebra-v1", "field": self.field.to_json(),
                "generators": list(self.generators)}

    @classmethod
    def from_json(cls, value: object) -> FreeAlgebra:
        if (not isinstance(value, dict) or set(value) != {"schema", "field", "generators"}
                or value["schema"] != "zipper-free-algebra-v1"
                or not isinstance(value["generators"], list)):
            raise ValueError("invalid free-algebra data")
        return cls(RationalFunctionField.from_json(value["field"]), tuple(value["generators"]))


@dataclass(frozen=True, slots=True, init=False, eq=False)
class AlgebraElement:
    """Immutable sparse free polynomial; the parent survives cancellation to zero."""

    algebra: FreeAlgebra
    terms: tuple[tuple[Word, FieldElement], ...]
    __hash__ = None  # type: ignore[assignment]

    def __init__(self, algebra: FreeAlgebra, terms: Mapping[Word, Scalar]) -> None:
        if not isinstance(algebra, FreeAlgebra):
            raise TypeError("an element requires a FreeAlgebra")
        if not isinstance(terms, Mapping):
            raise TypeError("terms must map words to exact coefficients")
        data = []
        for word, value in terms.items():
            if (not isinstance(word, tuple)
                    or any(type(index) is not int or not 1 <= index <= len(algebra.generators)
                           for index in word)):
                raise ValueError("words must contain valid one-based generator indices")
            coefficient = algebra.field.scalar(value)
            if not coefficient.is_zero():
                data.append((word, coefficient))
        object.__setattr__(self, "algebra", algebra)
        object.__setattr__(self, "terms", tuple(sorted(data)))

    def is_zero(self) -> bool:
        return not self.terms

    def _lift(self, other: object) -> AlgebraElement | NotImplementedType:
        if isinstance(other, AlgebraElement):
            if other.algebra != self.algebra:
                raise TypeError("free polynomials have different coefficient fields or alphabets")
            return other
        if isinstance(other, (int, Fraction, FieldElement)):
            return self.algebra.scalar(other)
        return NotImplemented

    def __add__(self, other: object) -> AlgebraElement:
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        terms = dict(self.terms)
        for word, coefficient in rhs.terms:
            terms[word] = terms.get(word, self.algebra.field.zero) + coefficient
        return self.algebra.element(terms)

    __radd__ = __add__

    def __neg__(self) -> AlgebraElement:
        return self.algebra.element({word: -coefficient for word, coefficient in self.terms})

    def __sub__(self, other: object) -> AlgebraElement:
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        return self + (-rhs)

    def __rsub__(self, other: object) -> AlgebraElement:
        lhs = self._lift(other)
        if lhs is NotImplemented:
            return NotImplemented
        return lhs + (-self)

    def __mul__(self, other: object) -> AlgebraElement:
        rhs = self._lift(other)
        if rhs is NotImplemented:
            return NotImplemented
        terms: dict[Word, FieldElement] = {}
        for left_word, left_coefficient in self.terms:
            for right_word, right_coefficient in rhs.terms:
                word = left_word + right_word
                terms[word] = (terms.get(word, self.algebra.field.zero)
                               + left_coefficient * right_coefficient)
        return self.algebra.element(terms)

    __rmul__ = __mul__

    def __pow__(self, exponent: int) -> AlgebraElement:
        if type(exponent) is not int or exponent < 0:
            raise ValueError("free-polynomial powers require a nonnegative integer")
        result = self.algebra.one
        power = self
        while exponent:
            if exponent % 2:
                result = result * power
            exponent //= 2
            if exponent:
                power = power * power
        return result

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AlgebraElement):
            return NotImplemented
        return self.algebra == other.algebra and self.terms == other.terms

    def to_json(self) -> dict:
        return {"schema": "zipper-algebra-element-v1", "algebra": self.algebra.to_json(),
                "terms": [[list(word), coefficient.to_json()]
                          for word, coefficient in self.terms]}

    @classmethod
    def from_json(cls, value: object) -> AlgebraElement:
        if (not isinstance(value, dict) or set(value) != {"schema", "algebra", "terms"}
                or value["schema"] != "zipper-algebra-element-v1"
                or not isinstance(value["terms"], list)):
            raise ValueError("invalid free-polynomial data")
        algebra = FreeAlgebra.from_json(value["algebra"])
        terms = {}
        for item in value["terms"]:
            if (not isinstance(item, list) or len(item) != 2
                    or not isinstance(item[0], list)):
                raise ValueError("a free-polynomial term needs a word and a coefficient")
            word = tuple(item[0])
            if any(type(index) is not int for index in word):
                raise ValueError("word indices must be integers")
            if word in terms:
                raise ValueError("duplicate free-polynomial word")
            terms[word] = FieldElement.from_json(item[1])
        return algebra.element(terms)


@dataclass(frozen=True, slots=True)
class AlgebraPresentation:
    """Literal ordered defining relations over a declared parameter field.

    No completion, dimension, monicity or unit-constant-term claim is inferred
    from this data. A content hash binds a witness to these precise inputs.
    """

    algebra: FreeAlgebra
    relations: tuple[AlgebraElement, ...]
    labels: tuple[str, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.algebra, FreeAlgebra):
            raise TypeError("a presentation requires a FreeAlgebra")
        if (not isinstance(self.relations, tuple)
                or any(not isinstance(relation, AlgebraElement)
                       or relation.algebra != self.algebra for relation in self.relations)):
            raise TypeError("relations must be a tuple of elements of the declared free algebra")
        if (not isinstance(self.labels, tuple) or len(self.labels) != len(self.relations)
                or any(not isinstance(label, str) or not label for label in self.labels)
                or len(set(self.labels)) != len(self.labels)):
            raise ValueError("relation labels must be distinct nonempty strings, one per relation")

    @classmethod
    def zipper(cls, n: int, *, field: RationalFunctionField, q: Scalar,
               untwistings: Mapping[int, Iterable[Scalar]],
               max_zipper_level: int | None = None) -> AlgebraPresentation:
        """Construct raw zipper relations with independently chosen field polynomials.

        The existing bare symbolic emitter supplies the zipper polynomials; its
        q is substituted exactly. No higher untwisting is inserted implicitly.
        At level i>2 the twist is the thesis's Theta_{i-1,j}.
        """
        from .coefficients import SymbolicRing
        from .presentation import Presentation

        bare = Presentation(n, untwist="none", max_zipper_level=max_zipper_level)
        algebra = FreeAlgebra(field, tuple(bare.var_names()))
        q_value = field.scalar(q)
        if not isinstance(untwistings, Mapping):
            raise TypeError("untwistings must map levels to coefficient sequences")
        polynomials = {}
        for level, values in untwistings.items():
            if type(level) is not int or not 2 <= level <= n:
                raise ValueError("untwisting levels must be integers from two to n")
            coefficients = tuple(field.scalar(value) for value in values)
            last = len(coefficients)
            while last and coefficients[last - 1].is_zero():
                last -= 1
            if not last:
                raise ValueError("an untwisting polynomial must be nonzero; omit absent levels")
            polynomials[level] = coefficients[:last]
        _, legacy_relations, legacy_sliders, tags = bare.build(SymbolicRing())

        def convert(polynomial):
            terms = {}
            for word, coefficient in polynomial.terms.items():
                # The bare emitter has no t terms. Substitution t=0 also makes
                # the complete bivariate substitution explicit here.
                terms[word] = sum((field.scalar(c) * q_value**i * field.zero**j
                                   for (i, j), c in coefficient.v.items()), field.zero)
            return algebra.element(terms)

        relations = [convert(relation) for relation in legacy_relations]
        labels = [tag.label(n) for tag in tags]
        for level, coefficients in sorted(polynomials.items()):
            for anchor in range(1, n - level + 2):
                if level == 2:
                    operator = algebra.generator(f"x{anchor}")
                else:
                    ascending = tuple(range(anchor, anchor + level - 1))
                    operator = algebra.element({ascending + ascending[::-1]: 1})
                polynomial = algebra.scalar(coefficients[-1])
                for coefficient in reversed(coefficients[:-1]):
                    polynomial = polynomial * operator + coefficient
                relations.append(polynomial * convert(legacy_sliders[level, anchor]))
                labels.append(f"T_{level}" + (f"_j{anchor}" if anchor != 1 else ""))
        return cls(algebra, tuple(relations), tuple(labels))

    def to_json(self) -> dict:
        return {"schema": "zipper-algebra-presentation-v1", "algebra": self.algebra.to_json(),
                "relations": [relation.to_json() for relation in self.relations],
                "labels": list(self.labels)}

    @classmethod
    def from_json(cls, value: object) -> AlgebraPresentation:
        if (not isinstance(value, dict)
                or set(value) != {"schema", "algebra", "relations", "labels"}
                or value["schema"] != "zipper-algebra-presentation-v1"
                or not isinstance(value["relations"], list)
                or not isinstance(value["labels"], list)):
            raise ValueError("invalid algebra-presentation data")
        return cls(FreeAlgebra.from_json(value["algebra"]),
                   tuple(AlgebraElement.from_json(item) for item in value["relations"]),
                   tuple(value["labels"]))

    def content_hash(self) -> str:
        """Full SHA-256 of the versioned literal data, including domain and alphabet."""
        data = json.dumps(self.to_json(), sort_keys=True, separators=(",", ":"), allow_nan=False)
        return hashlib.sha256(data.encode("utf-8")).hexdigest()

    def specialize(self, values: Mapping[str, int | Fraction], *,
                   characteristic: int | None = None) -> AlgebraPresentation:
        """Specialize stored expressions where their denominators survive.

        This operation does not define a map on the whole rational-function
        field. A zero defining relation is retained, with its original label.
        """
        field = self.algebra.field.one.specialize(values, characteristic=characteristic).field
        algebra = FreeAlgebra(field, self.algebra.generators)
        relations = tuple(algebra.element({
            word: coefficient.specialize(values, characteristic=characteristic)
            for word, coefficient in relation.terms}) for relation in self.relations)
        return AlgebraPresentation(algebra, relations, self.labels)


@dataclass(frozen=True, slots=True)
class CofactorIdentity:
    """An exact identity target = sum(left * defining_relation * right).

    Verification proves ideal membership in the explicitly declared field. It
    makes no assertion about a smaller coefficient ring or an excluded fibre.
    Relation references are zero-based indices in the ordered presentation.
    """

    presentation_hash: str
    target: AlgebraElement
    terms: tuple[tuple[AlgebraElement, int, AlgebraElement], ...]

    def __post_init__(self) -> None:
        if (not isinstance(self.presentation_hash, str)
                or re.fullmatch(r"[0-9a-f]{64}", self.presentation_hash) is None):
            raise ValueError("a cofactor identity requires a full presentation SHA-256")
        if not isinstance(self.target, AlgebraElement):
            raise TypeError("a cofactor target must be a free-polynomial element")
        if not isinstance(self.terms, tuple):
            raise TypeError("cofactor terms must be an immutable tuple")
        for term in self.terms:
            if (not isinstance(term, tuple) or len(term) != 3
                    or type(term[1]) is not int or term[1] < 0):
                raise ValueError("a cofactor term needs left, nonnegative relation index, right")
            if any(not isinstance(factor, AlgebraElement)
                   or factor.algebra != self.target.algebra for factor in (term[0], term[2])):
                raise TypeError("all cofactor factors must belong to the target's free algebra")

    def verify(self, presentation: AlgebraPresentation) -> bool:
        if not isinstance(presentation, AlgebraPresentation):
            raise TypeError("verification requires an AlgebraPresentation")
        if self.presentation_hash != presentation.content_hash():
            raise ValueError("cofactor presentation hash mismatch")
        if self.target.algebra != presentation.algebra:
            raise TypeError("cofactor target and presentation have different free algebras")
        result = presentation.algebra.zero
        for left, index, right in self.terms:
            if index >= len(presentation.relations):
                raise ValueError("cofactor relation index is out of range")
            result = result + left * presentation.relations[index] * right
        return result == self.target

    def to_json(self) -> dict:
        return {"schema": "zipper-cofactor-identity-v1", "presentation_hash": self.presentation_hash,
                "target": self.target.to_json(),
                "terms": [[left.to_json(), index, right.to_json()]
                          for left, index, right in self.terms]}

    @classmethod
    def from_json(cls, value: object) -> CofactorIdentity:
        if (not isinstance(value, dict)
                or set(value) != {"schema", "presentation_hash", "target", "terms"}
                or value["schema"] != "zipper-cofactor-identity-v1"
                or not isinstance(value["terms"], list)):
            raise ValueError("invalid cofactor-identity data")
        terms = []
        for item in value["terms"]:
            if not isinstance(item, list) or len(item) != 3:
                raise ValueError("cofactor data need left, relation index, right")
            terms.append((AlgebraElement.from_json(item[0]), item[1],
                          AlgebraElement.from_json(item[2])))
        return cls(value["presentation_hash"], AlgebraElement.from_json(value["target"]), tuple(terms))

    def save(self, path: str | Path) -> None:
        """Write a JSON identity with automatic source/output provenance."""
        from .provenance import capture_provenance, write_output_provenance

        producer = capture_provenance()
        Path(path).write_text(json.dumps(self.to_json(), sort_keys=True, indent=2) + "\n",
                              encoding="utf-8")
        write_output_provenance(path, producer=producer)

    @classmethod
    def load(cls, path: str | Path) -> CofactorIdentity:
        return cls.from_json(json.loads(Path(path).read_text(encoding="utf-8")))
