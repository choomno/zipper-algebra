"""Bounded independent checks in Z[q,q^-1][B_n] of local zipper identities.

The faithful Artin action checks braid-group equality exactly. An accepted
certificate proves lhs-rhs = sum c L T(Z_i,j) R, where Z is built here from
its original defining polynomial and T is I, M, phi or Mphi. Symmetry images
use the symmetry of the defining ideal. This is a braid-group-algebra
cofactor identity; braid/inverse equalities are absorbed by the group keys.
It is not an explicit cofactor over the maintained free associative input
presentation. This checker is optional and does not normalize an execution.

Derived BM rules additionally invert q and crossing signs. After canceling
their common slider unit, the cofactor has coefficient -q^(1-i)*a(q^-1),
bar-mirrored contexts, and the original relator's M image. The frozen catalog
is unchanged; the derived rule hash binds its source record and proof version.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

from .catalog import Catalog, CatalogError, expand_polynomial, expand_word, integer, instantiate_rule
from .semilinear import semilinear_zip_rule

Polynomial = dict[tuple[int, tuple[int, ...]], int]
SYMMETRIES = {"I": (False, False), "M": (False, True),
              "phi": (True, True), "Mphi": (True, False)}
SCOPE = "exact_braid_group_algebra_original_zipper_relator_identity"


class LocalVerificationError(ValueError):
    """A false, malformed or unsupported local identity certificate."""


class LocalVerificationLimit(LocalVerificationError):
    """The explicit optional verification budget was exhausted."""


@dataclass(frozen=True)
class VerificationLimits:
    max_terms: int = 100000
    max_word_size: int = 10000
    max_image_letters: int = 100000
    max_operations: int = 10000000
    max_integer_bits: int = 100000

    def __post_init__(self) -> None:
        if any(type(x) is not int or x < 1 for x in vars(self).values()):
            raise LocalVerificationError("verification limits must be positive integers")


class _Budget:
    def __init__(self, limits: VerificationLimits):
        self.limits = limits
        self.used = 0

    def consume(self, count: int = 1) -> None:
        if count > self.limits.max_operations - self.used:
            raise LocalVerificationLimit("Artin verification operation budget exceeded")
        self.used += count

    def integer(self, value: Any) -> int:
        if type(value) is not int:
            raise LocalVerificationError("exact integer required")
        if abs(value).bit_length() > self.limits.max_integer_bits:
            raise LocalVerificationLimit("integer encoding budget exceeded")
        return value

    def word(self, value: Any, n: int) -> tuple[int, ...]:
        if not isinstance(value, (tuple, list)):
            raise LocalVerificationError("signed braid word required")
        if len(value) > self.limits.max_word_size:
            raise LocalVerificationLimit("braid word budget exceeded")
        if any(type(x) is not int or not 1 <= abs(x) < n for x in value):
            raise LocalVerificationError("signed generator outside B_n")
        return tuple(value)


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _hash(value: Any) -> str:
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _source_hash() -> str:
    return _hash({name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                  for name in ("verification.py", "catalog.py", "semilinear.py")})


def _fields(value: Any, keys: set[str]) -> None:
    if not isinstance(value, dict) or set(value) != keys:
        raise LocalVerificationError("unexpected or missing certificate fields")


def _catalog_call(function: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    try:
        return function(*args, **kwargs)
    except CatalogError as exc:
        if "budget exceeded" in str(exc):
            raise LocalVerificationLimit(str(exc)) from exc
        raise LocalVerificationError("invalid catalog instantiation: " + str(exc)) from exc


def _add(out: Polynomial, q: int, word: tuple[int, ...], c: int, budget: _Budget) -> None:
    budget.consume()
    budget.integer(q)
    budget.integer(c)
    value = budget.integer(out.get((q, word), 0) + c)
    if value:
        out[q, word] = value
    else:
        out.pop((q, word), None)
    if len(out) > budget.limits.max_terms:
        raise LocalVerificationLimit("polynomial term budget exceeded")


def _multiply(a: Polynomial, b: Polynomial, budget: _Budget) -> Polynomial:
    out: Polynomial = {}
    for (q, w), c in a.items():
        for (p, v), d in b.items():
            if len(w) + len(v) > budget.limits.max_word_size:
                raise LocalVerificationLimit("braid word budget exceeded")
            _add(out, q + p, w + v, c * d, budget)
    return out


def _slider(size: int, start: int, budget: _Budget) -> Polynomial:
    # Independent defining expansion; no catalog slider implementation is used.
    out = {(1, (-start,)): 1, (0, ()): 1, (1, ()): -1, (0, (start,)): -1}
    for i in range(3, size + 1):
        positive = tuple(range(start, start + i - 1))
        out = _multiply({(i - 1, tuple(-x for x in positive)): 1, (0, positive): -1}, out, budget)
    return out


def _relator(n: int, size: int, start: int, budget: _Budget) -> Polynomial:
    if type(size) is not int or type(start) is not int or not 2 <= size < n or not 1 <= start <= n - size:
        raise LocalVerificationError("original zipper relator requires 2 <= i < n and a valid anchor")
    full = tuple(range(start, start + size))
    if size == 2:
        factor = _slider(2, start + 1, budget)
        _add(factor, 1, tuple(-x for x in full), -1, budget)
        _add(factor, 0, full, 1, budget)
    else:
        shorter = tuple(range(start + 1, start + size))
        factor = {(size - 1, tuple(-x for x in shorter)): 1, (0, shorter): -1,
                  (size - 1, tuple(-x for x in full)): -1, (0, full): 1}
    return _multiply(factor, _slider(size, start, budget), budget)


def _transform(word: tuple[int, ...], n: int, mode: str) -> tuple[int, ...]:
    if mode not in SYMMETRIES:
        raise LocalVerificationError("unknown symmetry")
    reverse, reflect = SYMMETRIES[mode]
    if reverse:
        word = tuple(reversed(word))
    if reflect:
        word = tuple((1 if x > 0 else -1) * (n - abs(x)) for x in word)
    return word


def _serialize(poly: Polynomial) -> list[dict[str, Any]]:
    return [{"integer": c, "q_power": q, "word": list(w)} for (q, w), c in sorted(poly.items())]


def _parse(value: Any, n: int, budget: _Budget) -> Polynomial:
    if not isinstance(value, list):
        raise LocalVerificationError("expanded polynomial must be a list")
    if len(value) > budget.limits.max_terms:
        raise LocalVerificationLimit("polynomial term budget exceeded")
    out: Polynomial = {}
    for term in value:
        _fields(term, {"integer", "q_power", "word"})
        _add(out, budget.integer(term["q_power"]), budget.word(term["word"], n),
             budget.integer(term["integer"]), budget)
    if _serialize(out) != value:
        raise LocalVerificationError("expanded polynomial is not canonically collected")
    return out


def _artin(word: tuple[int, ...], n: int, budget: _Budget) -> tuple[tuple[int, ...], ...]:
    # The unchanged ambient n has already passed _binding's image budget.
    images: list[tuple[int, ...]] = [(i,) for i in range(1, n + 1)]
    total = n
    for letter in word:
        budget.consume()
        j = abs(letter) - 1
        a, b = images[j], images[j + 1]
        pieces = (a, b, tuple(-x for x in reversed(a))) if letter > 0 else (tuple(-x for x in reversed(b)), a, b)
        stack: list[int] = []
        for piece in pieces:
            for x in piece:
                budget.consume()
                if stack and stack[-1] == -x:
                    stack.pop()
                else:
                    stack.append(x)
                if len(stack) > budget.limits.max_image_letters:
                    raise LocalVerificationLimit("Artin image budget exceeded")
        changed = tuple(stack)
        total += len(changed) - (len(b) if letter > 0 else len(a))
        if total > budget.limits.max_image_letters:
            raise LocalVerificationLimit("Artin image budget exceeded")
        if letter > 0:
            images[j], images[j + 1] = changed, a
        else:
            images[j], images[j + 1] = b, changed
    return tuple(images)


def _zero(poly: Polynomial, n: int, budget: _Budget) -> bool:
    collected: dict[Any, int] = defaultdict(int)
    for (q, word), coefficient in poly.items():
        key = q, _artin(word, n, budget)
        collected[key] = budget.integer(collected[key] + coefficient)
    return all(c == 0 for c in collected.values())


def _binding(parameters: Any, words: Any, budget: _Budget) -> tuple[dict[str, int], dict[str, tuple[int, ...]]]:
    if (not isinstance(parameters, dict) or any(not isinstance(k, str) for k in parameters)
            or "n" not in parameters or not isinstance(words, dict)
            or any(not isinstance(k, str) for k in words)):
        raise LocalVerificationError("explicit parameter and named-word bindings required")
    env = {k: budget.integer(v) for k, v in parameters.items()}
    n = env["n"]
    if n < 1:
        raise LocalVerificationError("positive ambient rank required")
    if n > budget.limits.max_image_letters:
        raise LocalVerificationLimit("Artin image budget exceeded")
    return env, {k: budget.word(w, n) for k, w in words.items()}


def _expanded_rule(rule: dict[str, Any], env: dict[str, int], words: dict[str, tuple[int, ...]], budget: _Budget) -> dict[str, Any]:
    return {side: _serialize(_catalog_call(expand_polynomial, rule[side], env, words=words, expand_sliders=True,
                                          max_terms=min(budget.limits.max_terms, budget.limits.max_word_size)))
            for side in ("lhs", "rhs")}


def _get_rule(catalog: Catalog, rule_id: str, env: dict[str, int]) -> dict[str, Any]:
    if not isinstance(rule_id, str):
        raise LocalVerificationError("string rule identifier required")
    if rule_id.endswith(".BM"):
        return _catalog_call(semilinear_zip_rule, catalog.get(rule_id[:-3] + ".M"), env)
    return catalog.get(rule_id)


def _check_accepted(rule: dict[str, Any], env: dict[str, int], words: dict[str, tuple[int, ...]], budget: _Budget) -> None:
    _catalog_call(instantiate_rule, rule, env, words=words,
                  max_terms=budget.limits.max_terms, require_verified=True)


def certify_rule(rule_id: str, parameters: dict[str, int], *, words: dict[str, list[int] | tuple[int, ...]] | None = None,
                 catalog: Catalog | None = None, limits: VerificationLimits | None = None) -> dict[str, Any]:
    """Produce and exactly check a finite local cofactor witness, or refuse it."""
    catalog = Catalog.bundled() if catalog is None else catalog
    budget = _Budget(VerificationLimits() if limits is None else limits)
    env, named = _binding(parameters, {} if words is None else words, budget)
    rule = _get_rule(catalog, rule_id, env)
    # This invokes the catalog guard and acceptance check before certificate production.
    _check_accepted(rule, env, named, budget)
    expanded = _expanded_rule(rule, env, named, budget)
    mode = rule["symmetry"]["name"]
    base = catalog.get(rule["base_id"] + ".I")
    cofactors = []
    metadata = base.get("original_relator_cofactor")
    if metadata is not None:
        rel = metadata["relator"]
        if rel["family"] != "Z":
            raise LocalVerificationError("unsupported original relator family")
        size, start = integer(rel["level"], env), integer(rel["start"], env)
        original = _relator(env["n"], size, start, budget)
        contexts = []
        for side in ("left", "right"):
            poly = _catalog_call(expand_word, metadata[side], env, words=named, max_terms=budget.limits.max_word_size)
            # Compact word AST expansion is one monomial of coefficient one.
            # The signed-word check below rejects any preserved slider token.
            contexts.append(budget.word(next(iter(poly))[1], env["n"]))
        if mode != "BM" and SYMMETRIES[mode][0]:
            contexts.reverse()
        relator_mode = "M" if mode == "BM" else mode
        left, right = (_transform(w, env["n"], relator_mode) for w in contexts)
        if mode == "BM":
            left, right = tuple(-x for x in left), tuple(-x for x in right)
        for coefficient in metadata["coefficient"]:
            scalar, power = coefficient["integer"], integer(coefficient["q_power"], env)
            if mode == "BM":
                # B(X_i)'s common unit is canceled in the derived local rule.
                # B(Z_i)/B(X_i)'s unit contributes -q^(1-i) times Z_i.
                scalar, power = -scalar, 1 - env["i"] - power
            cofactors.append({"integer": scalar, "q_power": power,
                              "left": list(left), "right": list(right),
                              "relator": {"family": "Z", "level": size, "start": start, "symmetry": relator_mode,
                                          "definition_sha256": _hash({"n": env["n"], "original_Z": _serialize(original)})}})
    certificate = {"schema": "zipper-local-original-relator-v1", "scope": SCOPE,
                   "coefficient_ring": "Z[q,q^-1]", "source_sha256": _source_hash(),
                   "catalog_sha256": catalog.sha256, "rule_sha256": _hash(rule), "rule_id": rule_id,
                   "parameters": env, "words": {k: list(v) for k, v in named.items()},
                   **expanded, "cofactors": cofactors}
    check_local_certificate(certificate, catalog=catalog, limits=limits)
    return certificate


def check_local_certificate(certificate: Any, *, catalog: Catalog | None = None,
                            limits: VerificationLimits | None = None) -> bool:
    """Check reloaded explicit cofactors directly against original definitions.

    A nonzero residual is refused, including identities lacking cofactor
    metadata. No membership search or matrix representation is substituted.
    """
    _fields(certificate, {"schema", "scope", "coefficient_ring", "source_sha256", "catalog_sha256", "rule_sha256",
                          "rule_id", "parameters", "words", "lhs", "rhs", "cofactors"})
    if (certificate["schema"] != "zipper-local-original-relator-v1" or certificate["scope"] != SCOPE
            or certificate["coefficient_ring"] != "Z[q,q^-1]" or certificate["source_sha256"] != _source_hash()):
        raise LocalVerificationError("unsupported certificate scope or source identity")
    catalog = Catalog.bundled() if catalog is None else catalog
    budget = _Budget(VerificationLimits() if limits is None else limits)
    env, words = _binding(certificate["parameters"], certificate["words"], budget)
    rule = _get_rule(catalog, certificate["rule_id"], env)
    if certificate["catalog_sha256"] != catalog.sha256 or certificate["rule_sha256"] != _hash(rule):
        raise LocalVerificationError("catalog or rule hash mismatch")
    left_side = _parse(certificate["lhs"], env["n"], budget)
    right_side = _parse(certificate["rhs"], env["n"], budget)
    _check_accepted(rule, env, words, budget)
    expected = _expanded_rule(rule, env, words, budget)
    if any(certificate[side] != expected[side] for side in ("lhs", "rhs")):
        raise LocalVerificationError("expanded sides differ from bound catalog rule")
    residual = left_side
    for (q, w), c in right_side.items():
        _add(residual, q, w, -c, budget)
    if not isinstance(certificate["cofactors"], list):
        raise LocalVerificationError("explicit cofactor list required")
    if len(certificate["cofactors"]) > budget.limits.max_terms:
        raise LocalVerificationLimit("cofactor term budget exceeded")
    for term in certificate["cofactors"]:
        _fields(term, {"integer", "q_power", "left", "right", "relator"})
        tag = term["relator"]
        _fields(tag, {"family", "level", "start", "symmetry", "definition_sha256"})
        if tag["family"] != "Z":
            raise LocalVerificationError("unsupported original relator family")
        original = _relator(env["n"], tag["level"], tag["start"], budget)
        if tag["definition_sha256"] != _hash({"n": env["n"], "original_Z": _serialize(original)}):
            raise LocalVerificationError("original relator definition hash mismatch")
        left, right = budget.word(term["left"], env["n"]), budget.word(term["right"], env["n"])
        coefficient, power = budget.integer(term["integer"]), budget.integer(term["q_power"])
        for (q, w), c in original.items():
            middle = _transform(w, env["n"], tag["symmetry"])
            whole = budget.word(left + middle + right, env["n"])
            _add(residual, q + power, whole, -coefficient * c, budget)
    if not _zero(residual, env["n"], budget):
        raise LocalVerificationError("nonzero braid-group-algebra cofactor residual; local identity unverified")
    return True
