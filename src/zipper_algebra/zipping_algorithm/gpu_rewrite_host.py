"""Admission and independent CPU oracle for packed catalog expansion.

Admission binds coordinates and computes lengths without materializing RHS
words. Full scheduler geometry remains outside this local expansion boundary.
The CPU oracle calls the unchanged interpreter and does not use the compiler.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Iterable

from .catalog import Catalog, CatalogError, condition, expand_polynomial
from .gpu_rewrite_program import (
    BLOCK, HALF_TWIST, SIGMA, SLIDER, WORD_REF,
    CompiledRule, Primitive, checked_i64, evaluate,
)
from .interpreter import Branch, Limits, LimitExceeded, ReplayError, Request, RunState, apply_rule


@dataclass(frozen=True)
class AdmittedCase:
    raw: dict[str, Any]
    branch: Branch
    request: Request
    limits: Limits
    program: CompiledRule
    values: tuple[int, ...]
    end: int
    lengths: tuple[int, ...]


def inputs(case: dict[str, Any]) -> tuple[dict[str, Any], Branch, Request, Limits]:
    """Normalize common input contracts without invoking packed compilation."""
    required = {"root", "id", "batch", "n", "rule", "parameters", "offset", "predecessor"}
    optional = {"words", "successor_controls", "limits"}
    if not isinstance(case, dict) or not required <= case.keys() or case.keys() - required - optional:
        raise ReplayError("Unexpected/missing packed rewrite case fields")
    raw = deepcopy(case)
    if any(type(raw[k]) is not int or raw[k] < 0 for k in ("root", "id", "batch")):
        raise ReplayError("Packed rewrite root/id/batch must be nonnegative integers")
    rule = raw["rule"]
    if not isinstance(rule, dict) or not isinstance(rule.get("id"), str):
        raise ReplayError("Packed rewrite requires an identified rule AST")
    try:
        branch = Branch.from_json(raw["predecessor"])
        request = Request(0, raw["offset"], rule["id"], raw["parameters"],
                          words=raw.get("words"),
                          successor_controls=raw.get("successor_controls"))
        limits = Limits(**raw["limits"]) if "limits" in raw else Limits(
            max_steps=100000, max_terms=100000, expansion_budget=1000000)
        limits.check(RunState(raw["n"], (branch,)))
    except (TypeError, KeyError) as error:
        raise ReplayError("Malformed packed rewrite case") from error
    return raw, branch, request, limits


def _primitive_length(
    primitive: Primitive, values: tuple[int, ...], n: int,
    program: CompiledRule, words: dict[str, tuple[int, ...]], budget: int,
) -> int:
    kind, a, b, c, d, _ = primitive
    if kind == WORD_REF:
        name = program.refs[a]
        if name not in words:
            raise CatalogError("Unbound packed word reference: " + name)
        word = words[name]
        if any(type(x) is not int or not 0 < abs(x) < n for x in word):
            raise CatalogError("Packed named words must be pure braid words in B_n")
        count = len(word)
    elif kind == SIGMA:
        index, sign = values[a], values[b]
        if not 0 < index < n or sign not in (-1, 1):
            raise CatalogError("Invalid packed signed generator")
        count = 1
    elif kind == BLOCK:
        length, start, repetitions, sign = values[a], values[b], values[c], values[d]
        if min(length, repetitions) < 0 or start < 1 or sign not in (-1, 1):
            raise CatalogError("Invalid packed block parameters")
        count = checked_i64(length * repetitions)
        if count and checked_i64(start + repetitions + length - 2) >= n:
            raise CatalogError("Packed block generator outside B_n")
    elif kind in (HALF_TWIST, SLIDER):
        size, start = values[a], values[b]
        if size < 1 or start < 1 or checked_i64(start + size - 1) > n:
            raise CatalogError("Invalid packed slider/half-twist interval")
        count = 1 if kind == SLIDER else checked_i64(size * (size - 1)) // 2
    else:
        raise CatalogError("Unknown packed primitive opcode")
    if count > budget:
        raise LimitExceeded("Packed expansion word budget exceeded")
    return count


def admit(case: dict[str, Any], program: CompiledRule) -> AdmittedCase:
    """Check an actual contextual occurrence and allocate only length metadata."""
    raw, branch, request, limits = inputs(case)
    env = dict(request.parameters)
    n = checked_i64(raw["n"])
    if env.get("n") != n:
        raise ReplayError("Packed rule ambient rank differs from predecessor")
    if raw["rule"].get("execution", {}).get("algebra_verified") is not True:
        raise ReplayError("Packed rule lacks accepted algebraic verification")
    if not condition(raw["rule"]["condition"], env):
        raise ReplayError("Packed rule guard is false")
    lhs = expand_polynomial(raw["rule"]["lhs"], env,
                            words=dict(request.words or {}), max_terms=limits.expansion_budget)
    if len(lhs) != 1:
        raise ReplayError("Packed contextual replay requires a monomial left side")
    (power, occurrence), scalar = next(iter(lhs.items()))
    if power != 0 or scalar != 1 or not occurrence:
        raise ReplayError("Packed left side must be a nonempty word with coefficient one")
    end = request.offset + len(occurrence)
    if branch.word[request.offset:end] != occurrence:
        raise ReplayError("Packed literal left side does not match occurrence")
    for token in branch.word:
        if isinstance(token, int):
            checked_i64(token)
        else:
            checked_i64(token[1])
            checked_i64(token[2])
    values = evaluate(program, env)
    # Only coordinates used by templates are required; unused bindings may
    # remain absent. Device parameter slots for unused names can be zero.
    words = {key: tuple(value) for key, value in (request.words or {}).items()}
    for word in words.values():
        for token in word:
            checked_i64(token)
    lengths = []
    contributions = 0
    exterior = len(branch.word) - len(occurrence)
    for term in program.terms:
        local_length = 0
        for primitive in term.primitives:
            local_length = checked_i64(local_length + _primitive_length(
                primitive, values, n, program, words, limits.expansion_budget))
        if local_length > limits.expansion_budget:
            raise LimitExceeded("Packed concatenation word budget exceeded")
        following = checked_i64(exterior + local_length)
        if any(factor for _, factor in term.coefficients) and following > limits.max_word_size:
            raise LimitExceeded("Packed successor word size limit exceeded")
        lengths.append(following)
        for q, _ in branch.coefficient:
            checked_i64(q)
            for power_register, factor in term.coefficients:
                if factor:
                    checked_i64(q + values[power_register])
                contributions += 1
    if contributions > limits.expansion_budget:
        raise LimitExceeded("Packed coefficient expansion budget exceeded")
    return AdmittedCase(raw, branch, request, limits, program, values, end, tuple(lengths))


class _RuleView(Catalog):
    """Independent interpreter view of the exact captured source rule."""

    def __init__(self, rule: dict[str, Any]):
        # inputs() has already isolated this rule; get() retains the ordinary
        # catalog's fresh-return contract without a redundant second snapshot.
        self.rule = rule

    def get(self, rule_id: str) -> dict[str, Any]:
        if rule_id != self.rule["id"]:
            raise ReplayError("Captured oracle rule identity mismatch")
        return deepcopy(self.rule)


def expand_cpu(cases: Iterable[dict[str, Any]]) -> tuple[dict[str, Any], ...]:
    """Run unchanged contextual expansion independently of the packed compiler."""
    result = []
    for case in cases:
        raw, branch, request, limits = inputs(case)
        following, _ = apply_rule(
            RunState(raw["n"], (branch,)), request,
            catalog=_RuleView(raw["rule"]), limits=limits, trace=False)
        result.append({"root": raw["root"], "id": raw["id"],
                       "branches": [item.to_json() for item in following.branches]})
    return tuple(result)
