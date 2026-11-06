"""Safe symbolic syntax and bounded expansion for the zipper rewrite catalog.

Words use algebraic left-to-right order, NOT chronological braid drawing order.
This module checks syntax and computes formal Laurent expansions. It does not
certify equality in a braid group or any quotient algebra.
"""

from __future__ import annotations
import json
import hashlib
import sys
from threading import RLock
from copy import deepcopy
from importlib.resources import files
from typing import Any
from collections import defaultdict, OrderedDict

PARAMETERS = frozenset(("n", "i", "j", "k", "l", "eps"))


class CatalogError(ValueError):
    pass


def integer(expr: Any, env: dict[str, int]) -> int:
    if type(expr) is int:
        return expr
    if isinstance(expr, str) and expr in PARAMETERS:
        value = env.get(expr)
        if type(value) is not int:
            raise CatalogError("Missing integer parameter: " + expr)
        return value
    if isinstance(expr, dict) and set(expr) == {"op", "args"}:
        op, args = expr["op"], expr["args"]
        if (
            op not in ("add", "sub", "mul")
            or not isinstance(args, list)
            or len(args) < 2
        ):
            raise CatalogError("Invalid integer operation")
        values = [integer(x, env) for x in args]
        if op == "add":
            return sum(values)
        if op == "sub":
            if len(values) != 2:
                raise CatalogError("sub has arity two")
            return values[0] - values[1]
        result = 1
        for x in values:
            result *= x
        return result
    raise CatalogError("Invalid integer expression: " + repr(expr))


def condition(expr: Any, env: dict[str, int]) -> bool:
    if type(expr) is bool:
        return expr
    if not isinstance(expr, dict):
        raise CatalogError("Invalid condition")
    op = expr.get("op")
    if (
        op in ("and", "or")
        and set(expr) == {"op", "args"}
        and isinstance(expr["args"], list)
    ):
        values = [condition(x, env) for x in expr["args"]]
        return all(values) if op == "and" else any(values)
    if op == "not" and set(expr) == {"op", "arg"}:
        return not condition(expr["arg"], env)
    if op in ("eq", "ne", "lt", "le", "gt", "ge") and set(expr) == {
        "op",
        "left",
        "right",
    }:
        a, b = integer(expr["left"], env), integer(expr["right"], env)
        return {
            "eq": a == b,
            "ne": a != b,
            "lt": a < b,
            "le": a <= b,
            "gt": a > b,
            "ge": a >= b,
        }[op]
    raise CatalogError("Invalid condition operation")


def _keys(node: Any, required: tuple[str, ...], optional: tuple[str, ...] = ()) -> None:
    if (
        not isinstance(node, dict)
        or set(node) - set(required) - set(optional)
        or not set(required) <= set(node)
    ):
        raise CatalogError("Unexpected/missing fields in " + repr(node))


def _clean(poly: dict[Any, int]) -> dict[Any, int]:
    return {key: value for key, value in poly.items() if value}


def _multiply(a: dict[Any, int], b: dict[Any, int], limit: int) -> dict[Any, int]:
    if len(a) * len(b) > limit:
        raise CatalogError("Expansion term budget exceeded")
    if a and b and max(len(w) for q, w in a) + max(len(w) for q, w in b) > limit:
        raise CatalogError("Expansion word budget exceeded")
    out: dict[Any, int] = defaultdict(int)
    for (qa, wa), ca in a.items():
        for (qb, wb), cb in b.items():
            out[(qa + qb, wa + wb)] += ca * cb
    # Collection has at most len(a)*len(b) keys, already bounded above.
    return _clean(out)


def _slider(size: int, start: int, env: dict[str, int], limit: int) -> dict[Any, int]:
    if size == 1:
        return {(0, ()): 1}
    base = {(1, (-start,)): 1, (0, ()): 1, (1, ()): -1, (0, (start,)): -1}
    for s in range(3, size + 1):
        positive = tuple(range(start, start + s - 1))
        factor = {(s - 1, tuple(-x for x in positive)): 1, (0, positive): -1}
        base = _multiply(factor, base, limit)
    return base


def expand_word(
    node: Any,
    env: dict[str, int],
    *,
    words: dict[str, Any] | None = None,
    expand_sliders: bool = False,
    max_terms: int = 100000,
) -> dict[Any, int]:
    """Return {(q_power, tuple(tokens)): integer}; braid tokens are signed ints.

    Preserved slider tokens are tuples ('slider',size,start,mirrored,reversed,
    barred). The latter flags describe actual polynomial transformations; they
    are never silently identified with the original slider.
    """
    words = {} if words is None else words
    if type(max_terms) is not int or max_terms < 1:
        raise CatalogError("Expansion budget must be positive")
    if not isinstance(node, dict):
        raise CatalogError("Word must be an AST object")
    kind = node.get("type")
    unit = {(0, ()): 1}

    def expand(child):
        return expand_word(
            child, env, words=words, expand_sliders=expand_sliders, max_terms=max_terms
        )

    def iv(key):
        return integer(node[key], env)

    n = env.get("n")
    if type(n) is not int or n < 1:
        raise CatalogError("n must be a positive integer")

    def braid(seq):
        seq = tuple(seq)
        if len(seq) > max_terms:
            raise CatalogError("Expansion word budget exceeded")
        if any(type(a) is not int or a == 0 or abs(a) >= n for a in seq):
            raise CatalogError("Braid generator outside B_n")
        return {(0, seq): 1}

    if kind == "identity":
        _keys(node, ("type",))
        return unit
    if kind == "sigma":
        _keys(node, ("type", "index"), ("sign",))
        idx = iv("index")
        sign = integer(node.get("sign", 1), env)
        if idx < 1 or sign not in (-1, 1):
            raise CatalogError("Invalid signed generator")
        return braid([sign * idx])
    if kind in ("fat", "block"):
        _keys(
            node,
            ("type", "length", "start")
            + (("count", "order") if kind == "block" else ()),
            ("sign",),
        )
        length, start = iv("length"), iv("start")
        count = iv("count") if kind == "block" else 1
        sign = integer(node.get("sign", 1), env)
        if length < 0 or count < 0 or start < 1 or sign not in (-1, 1):
            raise CatalogError("Invalid block parameters")
        order = node.get("order", "ascending")
        if order not in ("ascending", "descending"):
            raise CatalogError("Invalid block order")
        if length == 0 or count == 0:
            return unit
        if length * count > max_terms:
            raise CatalogError("Expansion word budget exceeded")
        offsets = range(count) if order == "ascending" else reversed(range(count))
        return braid(sign * (start + t + x) for t in offsets for x in range(length))
    if kind == "half_twist":
        _keys(node, ("type", "size", "start"))
        size, start = iv("size"), iv("start")
        if size < 1 or start < 1 or start + size - 1 > n:
            raise CatalogError("Invalid half-twist interval")
        if size * (size - 1) // 2 > max_terms:
            raise CatalogError("Expansion word budget exceeded")
        return braid(
            x for length in range(size - 1, 0, -1) for x in range(start, start + length)
        )
    if kind == "slider":
        _keys(node, ("type", "size", "start"), ("orientation",))
        size, start = iv("size"), iv("start")
        orientation = node.get("orientation", "normal")
        if (
            size < 1
            or start < 1
            or start + size - 1 > n
            or orientation not in ("normal", "mirrored")
        ):
            raise CatalogError("Invalid slider")
        mirrored = orientation == "mirrored"
        if expand_sliders and size > max_terms:
            raise CatalogError("Expansion word budget exceeded")
        if not expand_sliders:
            return {(0, (("slider", size, start, mirrored, False, False),)): 1}
        result = _slider(size, start, env, max_terms)
        if mirrored:
            result = {
                (
                    q,
                    tuple(
                        (1 if x > 0 else -1) * (2 * start + size - 2 - abs(x))
                        for x in w
                    ),
                ): c
                for (q, w), c in result.items()
            }
        return result
    if kind == "concat":
        _keys(node, ("type", "items"))
        if not isinstance(node["items"], list):
            raise CatalogError("concat.items must be a list")
        out = unit
        for child in node["items"]:
            out = _multiply(out, expand(child), max_terms)
        return out
    if kind == "word_ref":
        _keys(node, ("type", "name"))
        name = node["name"]
        if not isinstance(name, str) or name not in words:
            raise CatalogError("Unbound word reference")
        return braid(words[name])
    if kind in ("inverse", "bar", "reverse", "reflect"):
        _keys(node, ("type", "word"), ("n",) if kind == "reflect" else ())
        original = expand(node["word"])
        if kind == "inverse" and (
            len(original) != 1
            or next(iter(original.values())) != 1
            or next(iter(original))[0] != 0
            or any(not isinstance(x, int) for x in next(iter(original))[1])
        ):
            raise CatalogError(
                "inverse is only defined for braid words, not slider polynomials"
            )
        width = integer(node.get("n", n), env)
        if kind == "reflect" and width != n:
            raise CatalogError("reflect.n must be the ambient n")
        out = defaultdict(int)
        for (q, w), c in original.items():
            if kind in ("inverse", "reverse"):
                w = tuple(reversed(w))
            transformed = []
            for token in w:
                if isinstance(token, int):
                    if kind in ("inverse", "bar"):
                        token = -token
                    elif kind == "reflect":
                        token = (1 if token > 0 else -1) * (n - abs(token))
                else:
                    _, size, start, mirrored, reversed_, barred = token
                    if kind == "bar":
                        barred = not barred
                    elif kind == "reverse":
                        reversed_ = not reversed_
                    else:  # Reflection; inverse rejected slider tokens above.
                        start, mirrored = n - start - size + 2, not mirrored
                    token = ("slider", size, start, mirrored, reversed_, barred)
                transformed.append(token)
            out[(q, tuple(transformed))] += c
        return _clean(out)
    raise CatalogError("Unknown word node: " + repr(kind))


def expand_polynomial(
    terms: Any,
    env: dict[str, int],
    *,
    words: dict[str, Any] | None = None,
    expand_sliders: bool = False,
    max_terms: int = 100000,
) -> dict[Any, int]:
    if not isinstance(terms, list):
        raise CatalogError("Polynomial must be a list")
    out: dict[Any, int] = defaultdict(int)
    for term in terms:
        if not isinstance(term, dict):
            raise CatalogError("Invalid polynomial term")
        _keys(term, ("coefficient", "word"))
        if not isinstance(term["coefficient"], list):
            raise CatalogError("Coefficient must be Laurent term list")
        word = expand_word(
            term["word"],
            env,
            words=words,
            expand_sliders=expand_sliders,
            max_terms=max_terms,
        )
        for coeff in term["coefficient"]:
            _keys(coeff, ("integer", "q_power"))
            if type(coeff["integer"]) is not int:
                raise CatalogError("Coefficient numerator must be integer")
            power = integer(coeff["q_power"], env)
            for (q, w), c in word.items():
                out[(q + power, w)] += coeff["integer"] * c
        if len(out) > max_terms:
            raise CatalogError("Expansion term budget exceeded")
    return _clean(out)


def serialize(poly: dict[Any, int]) -> list[dict[str, Any]]:
    """JSON-safe exact output, sorted without assuming token comparability."""
    return [
        {
            "integer": c,
            "q_power": q,
            "word": [list(x) if isinstance(x, tuple) else x for x in w],
        }
        for (q, w), c in sorted(poly.items(), key=lambda item: repr(item[0]))
    ]


def instantiate_rule(
    rule: dict[str, Any],
    env: dict[str, int],
    *,
    words: dict[str, Any] | None = None,
    expand_sliders: bool = False,
    max_terms: int = 100000,
    require_verified: bool = False,
) -> dict[str, Any]:
    verified = rule.get("execution", {}).get("algebra_verified", False) is True
    if require_verified and not verified:
        raise CatalogError("Rule lacks reviewed algebraic acceptance")
    if not condition(rule["condition"], env):
        raise CatalogError("Rule condition is false")
    result = {
        "id": rule["id"],
        "parameters": dict(env),
        "algebra_verified": verified,
        "scope": "Formal instantiation only; control preconditions are not checked.",
    }
    for side in ("lhs", "rhs"):
        result[side] = serialize(
            expand_polynomial(
                rule[side],
                env,
                words=words,
                expand_sliders=expand_sliders,
                max_terms=max_terms,
            )
        )
    return result


def validate_catalog(catalog: Any) -> bool:
    """Validate against bundled Schema without requiring external packages.

    The small interpreter implements exactly the keywords used by schema.json.
    Full Draft 2020-12 validators can also consume that schema directly.
    """
    schema = json.loads(_data("schema.json").decode("utf-8"))

    def check(value, s, path):
        if "$ref" in s:
            return check(value, schema["$defs"][s["$ref"].split("/")[-1]], path)
        if "oneOf" in s:
            passed = 0
            for child in s["oneOf"]:
                try:
                    check(value, child, path)
                    passed += 1
                except CatalogError:
                    pass
            if passed != 1:
                raise CatalogError(path + ": no unique schema alternative")
            return
        if "const" in s and (value != s["const"] or type(value) != type(s["const"])):
            raise CatalogError(path + ": wrong constant")
        if "enum" in s and not any(
            value == x and type(value) == type(x) for x in s["enum"]
        ):
            raise CatalogError(path + ": invalid enum value")
        expected = s.get("type")
        types = {
            "object": dict,
            "array": list,
            "integer": int,
            "boolean": bool,
            "string": str,
        }
        if expected and type(value) is not types[expected]:
            raise CatalogError(path + ": expected " + expected)
        if expected == "object":
            for key in s.get("required", []):
                if key not in value:
                    raise CatalogError(path + ": missing " + key)
            props = s.get("properties", {})
            if s.get("additionalProperties") is False and set(value) - set(props):
                raise CatalogError(path + ": unexpected properties")
            for key, child in props.items():
                if key in value:
                    check(value[key], child, path + "." + key)
        if expected == "array":
            if len(value) < s.get("minItems", 0) or len(value) > s.get(
                "maxItems", len(value)
            ):
                raise CatalogError(path + ": invalid length")
            if "items" in s:
                for i, item in enumerate(value):
                    check(item, s["items"], path + f"[{i}]")
        if expected == "string" and len(value) < s.get("minLength", 0):
            raise CatalogError(path + ": empty string")

    check(catalog, schema, "catalog")
    ids = [r["id"] for r in catalog["rules"]]
    if len(ids) != len(set(ids)):
        raise CatalogError("Duplicate rule id")
    return True


CATALOG_VERSION = "2026-11-05.1"
CATALOG_SHA256 = "27c5dd71d09718a5561d18c802f714e07abac87890dc87169215369541462470"
SCHEMA_SHA256 = "ec586a1b9c4ff3ebc21d17688f46bd683a90e88b082e0f1ac11d945cd8334906"


def _data(name: str) -> bytes:
    return files("zipper_algebra.zipping_algorithm").joinpath("data", name).read_bytes()


class Catalog:
    """Versioned exact local identities, with separate author-review status.

    Parameter guards and formal expansions are checked here. Geometric labels,
    admissible processed bands and scheduler coverage require control evidence.
    Records returned to callers are independent copies of the frozen data.
    """

    def __init__(self, catalog: dict[str, Any], *, source_sha256: str | None = None):
        validate_catalog(catalog)
        self._catalog = deepcopy(catalog)
        self._rules = {rule["id"]: rule for rule in self._catalog["rules"]}
        self.source_sha256 = source_sha256
        self.sha256 = source_sha256
        self.author_review = "partial"
        self.version = CATALOG_VERSION
        self._rule_cache: OrderedDict[bytes, bytes] = OrderedDict()
        self._rule_cache_lock = RLock()
        self._rule_cache_max_bytes = 4 * 1024 * 1024
        self._rule_cache_max_entries = 256
        self._rule_cache_bytes = 0
        self._rule_cache_hits = 0
        self._rule_cache_misses = 0
        self._rule_cache_evictions = 0

    def configure_rule_cache(self, *, max_bytes: int, max_entries: int) -> None:
        """Bound retained derived-rule JSON; eviction never discards a calculation."""
        if any(type(x) is not int or x < 0 for x in (max_bytes, max_entries)):
            raise CatalogError("nonnegative integer rule-cache budgets required")
        with self._rule_cache_lock:
            self._rule_cache_max_bytes, self._rule_cache_max_entries = max_bytes, max_entries
            self._trim_rule_cache()

    def _trim_rule_cache(self) -> None:
        while self._rule_cache and (
            len(self._rule_cache) > self._rule_cache_max_entries
            or self._rule_cache_bytes > self._rule_cache_max_bytes
        ):
            key, value = self._rule_cache.popitem(last=False)
            self._rule_cache_bytes -= len(key) + len(value)
            self._rule_cache_evictions += 1

    def rule_cache_clear(self) -> None:
        """Release derived rules and reset their counters."""
        with self._rule_cache_lock:
            self._rule_cache.clear()
            self._rule_cache_bytes = 0
            self._rule_cache_hits = self._rule_cache_misses = self._rule_cache_evictions = 0

    def rule_cache_info(self) -> dict[str, int]:
        """Payload accounting and an object-size estimate, not a process RSS limit."""
        with self._rule_cache_lock:
            return self._rule_cache_info()

    def _rule_cache_info(self) -> dict[str, int]:
        return {
            "hits": self._rule_cache_hits, "misses": self._rule_cache_misses,
            "evictions": self._rule_cache_evictions, "entries": len(self._rule_cache),
            "retained_bytes": self._rule_cache_bytes,
            "max_bytes": self._rule_cache_max_bytes, "max_entries": self._rule_cache_max_entries,
            "resident_estimate_bytes": sys.getsizeof(self._rule_cache) + sys.getsizeof(self._rule_cache_lock) + sum(
                sys.getsizeof(key) + sys.getsizeof(value) for key, value in self._rule_cache.items()
            ),
        }

    def semilinear_rule(self, source_id: str, env: dict[str, int]) -> dict[str, Any]:
        """Reuse a parameterized BM derivation with fresh, validated bindings.

        The cache belongs to this frozen catalog. Widths, exact q powers and
        all bindings remain part of the derivation; no slider sizes are
        identified. Immutable JSON backing prevents mutation through callers.
        Expansion and resource checks still run when the rule is instantiated.
        """
        bindings = dict(env)
        with self._rule_cache_lock:
            return self._semilinear_rule(source_id, bindings)

    def _semilinear_rule(self, source_id: str, env: dict[str, int]) -> dict[str, Any]:
        from .semilinear import (
            _bind_semilinear_template, _semilinear_template, validate_semilinear_bindings,
        )

        if source_id not in self._rules:
            raise CatalogError("Unknown rule id: " + source_id)
        rule = self._rules[source_id]
        validate_semilinear_bindings(rule, env)
        if not self._rule_cache_max_entries or not self._rule_cache_max_bytes:
            self._rule_cache_misses += 1
            return _bind_semilinear_template(_semilinear_template(rule), env)
        key = source_id.encode()
        if key in self._rule_cache:
            self._rule_cache_hits += 1
            self._rule_cache.move_to_end(key)
            return _bind_semilinear_template(json.loads(self._rule_cache[key]), env)
        self._rule_cache_misses += 1
        result = _semilinear_template(rule)
        payload = json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
        if self._rule_cache_max_entries and len(key) + len(payload) <= self._rule_cache_max_bytes:
            self._rule_cache[key] = payload
            self._rule_cache_bytes += len(key) + len(payload)
            self._trim_rule_cache()
        return _bind_semilinear_template(result, env)

    @classmethod
    def bundled(cls) -> Catalog:
        """Load the immutable catalog and check both execution-data hashes."""
        payload = _data("full-catalog.json")
        schema = _data("schema.json")
        if hashlib.sha256(payload).hexdigest() != CATALOG_SHA256:
            raise CatalogError("Bundled catalog hash mismatch")
        if hashlib.sha256(schema).hexdigest() != SCHEMA_SHA256:
            raise CatalogError("Bundled schema hash mismatch")
        return cls(json.loads(payload), source_sha256=CATALOG_SHA256)

    @property
    def rules(self) -> tuple[dict[str, Any], ...]:
        return tuple(deepcopy(list(self._rules.values())))

    def select_rules(
        self,
        env: dict[str, int],
        *,
        prefix: str = "",
        suffix: str = "",
        exclude_prefixes: tuple[str, ...] = (),
    ) -> tuple[dict[str, Any], ...]:
        """Return independent copies of applicable rules in catalog order.

        Test the identifiers and conditions on the validated internal records
        before copying. Controllers need only a few applicable rules, while
        ``rules`` continues to provide a defensive copy of the full catalog.
        """
        return tuple(deepcopy(self._rules[rule_id]) for rule_id in self.select_rule_ids(
            env, prefix=prefix, suffix=suffix, exclude_prefixes=exclude_prefixes
        ))

    def select_rule_ids(
        self, env: dict[str, int], *, prefix: str = "", suffix: str = "",
        exclude_prefixes: tuple[str, ...] = (),
    ) -> tuple[str, ...]:
        """Select immutable family IDs before materializing any rule record."""
        return tuple(
            rule_id
            for rule_id, rule in self._rules.items()
            if rule_id.startswith(prefix)
            and rule_id.endswith(suffix)
            and not rule_id.startswith(exclude_prefixes)
            and condition(rule["condition"], env)
        )

    def get(self, rule_id: str) -> dict[str, Any]:
        try:
            return deepcopy(self._rules[rule_id])
        except KeyError as exc:
            raise CatalogError("Unknown rule id: " + rule_id) from exc

    def instantiate(
        self,
        rule_id: str,
        env: dict[str, int],
        *,
        words: dict[str, list[int] | tuple[int, ...]] | None = None,
        expand_sliders: bool = False,
        max_terms: int = 100000,
        require_verified: bool = True,
    ) -> dict[str, Any]:
        return instantiate_rule(
            self.get(rule_id),
            env,
            words=words,
            expand_sliders=expand_sliders,
            max_terms=max_terms,
            require_verified=require_verified,
        )
