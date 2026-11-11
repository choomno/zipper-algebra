"""Exact contextual catalog replay with bounded, resumable execution.

This layer executes an explicitly supplied schedule of algebraic identities.
It checks literal contexts and catalog guards. Control evidence is recorded data;
its presence does not prove that a state is reachable by the zipper scheduler.
"""
from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from .catalog import Catalog, CatalogError, expand_polynomial

Token = int | tuple[Any, ...]
Word = tuple[Token, ...]
Laurent = tuple[tuple[int, int], ...]


class ReplayError(ValueError):
    """An invalid state, contextual application, or replay certificate."""


class LimitExceeded(ReplayError):
    """A refusal to complete a transition; the predecessor is retained."""


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _digest(value: Any) -> str:
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _source_sha256() -> str:
    """Bind both transition semantics and the formal expansion implementation."""
    return _digest({name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                    for name in ("interpreter.py", "catalog.py")})


def _fields(value: Any, keys: set[str]) -> None:
    if not isinstance(value, dict) or set(value) != keys:
        raise ReplayError("unexpected or missing record fields")


def _word(value: Iterable[Any]) -> Word:
    return tuple(tuple(x) if isinstance(x, list) else x for x in value)


def _coefficient(value: Iterable[Iterable[int]]) -> Laurent:
    pairs = tuple(tuple(x) for x in value)
    if not pairs or any(len(x) != 2 or any(type(v) is not int for v in x) or x[1] == 0 for x in pairs):
        raise ReplayError("a nonzero Laurent coefficient is required")
    if tuple(sorted(pairs)) != pairs or len({x[0] for x in pairs}) != len(pairs):
        raise ReplayError("Laurent exponents must be distinct and increasing")
    return tuple((x[0], x[1]) for x in pairs)


def _multiply(a: Laurent, b: Laurent) -> Laurent:
    out: Counter[int] = Counter()
    for p, c in a:
        for q, d in b:
            out[p + q] += c * d
    return tuple((p, c) for p, c in sorted(out.items()) if c)


@dataclass(frozen=True)
class Branch:
    """One uncollected word and its exact Laurent coefficient.

    Control data are copied into canonical JSON, so sibling branches cannot
    share mutable algorithm state. No braid cancellation is implicit.
    """

    word: Word
    coefficient: Laurent = ((0, 1),)
    control_json: str = "{}"

    def __post_init__(self) -> None:
        object.__setattr__(self, "word", _word(self.word))
        object.__setattr__(self, "coefficient", _coefficient(self.coefficient))
        control = json.loads(self.control_json)
        if not isinstance(control, dict):
            raise ReplayError("branch control must be a JSON object")
        object.__setattr__(self, "control_json", _json(control))

    @property
    def control(self) -> dict[str, Any]:
        return json.loads(self.control_json)

    def to_json(self) -> dict[str, Any]:
        return {"word": [list(x) if isinstance(x, tuple) else x for x in self.word],
                "coefficient": [list(x) for x in self.coefficient], "control": self.control}

    @classmethod
    def from_json(cls, value: Any) -> Branch:
        _fields(value, {"word", "coefficient", "control"})
        return cls(_word(value["word"]), _coefficient(value["coefficient"]), _json(value["control"]))


@dataclass(frozen=True)
class RunState:
    n: int
    branches: tuple[Branch, ...]

    def __post_init__(self) -> None:
        if type(self.n) is not int or self.n < 1:
            raise ReplayError("rank must be a positive integer")
        object.__setattr__(self, "branches", tuple(self.branches))
        for branch in self.branches:
            if not isinstance(branch, Branch):
                raise ReplayError("state entries must be Branch records")
            for token in branch.word:
                if type(token) is int:
                    valid = 0 < abs(token) < self.n
                else:
                    valid = (isinstance(token, tuple) and len(token) == 6 and token[0] == "slider"
                             and type(token[1]) is int and type(token[2]) is int
                             and 1 <= token[1] <= self.n and 1 <= token[2] <= self.n - token[1] + 1
                             and all(type(x) is bool for x in token[3:]))
                if not valid:
                    raise ReplayError("word contains an invalid braid or slider token")

    def to_json(self) -> dict[str, Any]:
        return {"n": self.n, "branches": [x.to_json() for x in self.branches]}

    @classmethod
    def from_json(cls, value: Any) -> RunState:
        _fields(value, {"n", "branches"})
        return cls(value["n"], tuple(Branch.from_json(x) for x in value["branches"]))

    @property
    def sha256(self) -> str:
        return _digest(self.to_json())


@dataclass(frozen=True)
class Request:
    """A literal occurrence and explicit branch-control updates.

    ``successor_controls`` are ordered by first occurrence of each RHS word in
    the expanded catalog polynomial. Omission copies predecessor control data.
    These updates require separate scheduler verification.
    """

    branch: int
    offset: int
    rule_id: str
    parameters: Mapping[str, int]
    words: Mapping[str, Iterable[int]] | None = None
    evidence: Mapping[str, Any] | None = None
    successor_controls: Iterable[Mapping[str, Any]] | None = None
    _evidence_json: str = field(init=False, repr=False)
    _controls_json: tuple[str, ...] | None = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if any(type(x) is not int or x < 0 for x in (self.branch, self.offset)):
            raise ReplayError("branch and offset must be nonnegative integers")
        if not isinstance(self.rule_id, str) or not self.rule_id:
            raise ReplayError("rule id must be nonempty")
        parameters = dict(self.parameters)
        if any(not isinstance(k, str) or type(v) is not int for k, v in parameters.items()):
            raise ReplayError("parameters must map names to integers")
        words = {k: tuple(v) for k, v in (self.words or {}).items()}
        if any(not isinstance(k, str) or any(type(x) is not int for x in v) for k, v in words.items()):
            raise ReplayError("named words must contain signed braid integers")
        # Mapping proxies preserve a convenient mapping API while preventing edits.
        from types import MappingProxyType
        object.__setattr__(self, "parameters", MappingProxyType(parameters))
        object.__setattr__(self, "words", MappingProxyType(words))
        object.__setattr__(self, "_evidence_json", _json(dict(self.evidence or {})))
        object.__setattr__(self, "_controls_json", None)
        if self.successor_controls is not None:
            object.__setattr__(self, "_controls_json", tuple(_json(dict(x)) for x in self.successor_controls))

    def to_json(self) -> dict[str, Any]:
        return {"branch": self.branch, "offset": self.offset, "rule_id": self.rule_id,
                "parameters": dict(self.parameters), "words": {k: list(v) for k, v in (self.words or {}).items()},
                "evidence": json.loads(self._evidence_json),
                "successor_controls": None if self._controls_json is None else
                [json.loads(x) for x in self._controls_json]}

    @classmethod
    def from_json(cls, value: Any) -> Request:
        _fields(value, {"branch", "offset", "rule_id", "parameters", "words", "evidence", "successor_controls"})
        return cls(**value)


@dataclass(frozen=True)
class Limits:
    max_steps: int = 10000
    max_terms: int = 10000
    max_word_size: int = 100000
    max_coefficient_support: int = 10000
    max_integer_bits: int = 100000
    expansion_budget: int = 100000

    def __post_init__(self) -> None:
        if any(type(x) is not int or x < 1 for x in vars(self).values()):
            raise ReplayError("execution limits must be positive integers")

    def check(self, state: RunState) -> None:
        if len(state.branches) > self.max_terms:
            raise LimitExceeded("term limit exceeded")
        for branch in state.branches:
            if len(branch.word) > self.max_word_size:
                raise LimitExceeded("word size limit exceeded")
            if len(branch.coefficient) > self.max_coefficient_support:
                raise LimitExceeded("coefficient support limit exceeded")
            if any(abs(c).bit_length() > self.max_integer_bits for _, c in branch.coefficient):
                raise LimitExceeded("coefficient integer bit limit exceeded")


def apply_rule(state: RunState, request: Request, *, catalog: Catalog | None = None,
               limits: Limits | None = None, trace: bool = True) -> tuple[RunState, dict[str, Any]]:
    """Apply one accepted identity, retaining both complete exterior words."""
    if type(trace) is not bool:
        raise ReplayError("trace must be a boolean")
    catalog = Catalog.bundled() if catalog is None else catalog
    limits = Limits() if limits is None else limits
    limits.check(state)
    if request.branch >= len(state.branches):
        raise ReplayError("selected branch does not exist")
    if request.parameters.get("n") != state.n:
        raise ReplayError("rule ambient rank differs from state rank")
    rule = catalog.get(request.rule_id)
    if rule.get("execution", {}).get("algebra_verified") is not True:
        raise ReplayError("rule lacks accepted algebraic verification")
    from .catalog import condition
    if not condition(rule["condition"], dict(request.parameters)):
        raise ReplayError("rule guard is false")
    try:
        sides = [expand_polynomial(rule[side], dict(request.parameters), words=dict(request.words or {}),
                                   max_terms=limits.expansion_budget) for side in ("lhs", "rhs")]
    except CatalogError as exc:
        if "budget exceeded" in str(exc):
            raise LimitExceeded(str(exc)) from exc
        raise ReplayError("invalid catalog instantiation: " + str(exc)) from exc
    lhs, rhs = sides
    if len(lhs) != 1:
        raise ReplayError("contextual replay requires a monomial left side")
    (power, occurrence), scalar = next(iter(lhs.items()))
    if power != 0 or scalar != 1 or not occurrence:
        raise ReplayError("left side must be a nonempty word with coefficient one")
    predecessor = state.branches[request.branch]
    end = request.offset + len(occurrence)
    if predecessor.word[request.offset:end] != occurrence:
        raise ReplayError("literal left side does not match selected occurrence")
    grouped: dict[Word, Counter[int]] = {}
    for (q, word), c in rhs.items():
        grouped.setdefault(word, Counter())[q] += c
    rhs_terms = [(word, tuple((q, c) for q, c in sorted(coef.items()) if c))
                 for word, coef in grouped.items() if any(coef.values())]
    controls = request._controls_json
    if controls is not None and len(controls) != len(rhs_terms):
        raise ReplayError("one control update is required per distinct RHS word")
    successors = tuple(Branch(predecessor.word[:request.offset] + word + predecessor.word[end:],
                              _multiply(predecessor.coefficient, coefficient),
                              predecessor.control_json if controls is None else controls[idx])
                       for idx, (word, coefficient) in enumerate(rhs_terms))
    result = RunState(state.n, state.branches[:request.branch] + successors + state.branches[request.branch + 1:])
    limits.check(result)
    if not trace:
        return result, {}
    transition = {"request": request.to_json(), "rule_sha256": _digest(rule),
                  "predecessor_state_sha256": state.sha256, "successor_state_sha256": result.sha256,
                  "predecessor": predecessor.to_json(), "successors": [x.to_json() for x in successors],
                  "equality_reference": {"catalog_rule": request.rule_id, "kind": "accepted_catalog_identity"}}
    return result, transition


def metrics(state: RunState) -> dict[str, Any]:
    """Current uncollected size and exact coefficient encoding metrics."""
    collected: dict[Word, Counter[int]] = {}
    for branch in state.branches:
        coefficient = collected.setdefault(branch.word, Counter())
        for q, c in branch.coefficient:
            coefficient[q] += c
    exponents = [q for branch in state.branches for q, _ in branch.coefficient]
    return {"uncollected_terms": len(state.branches),
            "collected_terms": sum(any(c.values()) for c in collected.values()),
            "word_tokens": sum(len(b.word) for b in state.branches),
            "max_word_tokens": max((len(b.word) for b in state.branches), default=0),
            "coefficient_support": sum(len(b.coefficient) for b in state.branches),
            "q_exponent_min": min(exponents, default=0), "q_exponent_max": max(exponents, default=0),
            "max_integer_bits": max((abs(c).bit_length() for b in state.branches for _, c in b.coefficient), default=0)}


def execute(state: RunState, requests: Sequence[Request], *, catalog: Catalog | None = None,
            limits: Limits | None = None) -> dict[str, Any]:
    """Replay a finite supplied schedule, saving an atomic predecessor on limits.

    The schedule must be a finite Sequence, because resumable output records
    every pending request. Arbitrary iterators are rejected without consuming
    them. ``schedule_exhausted`` means that the supplied list is exhausted; it makes
    no assertion about zipping completion. ``incomplete`` retains pending work.
    """
    catalog = Catalog.bundled() if catalog is None else catalog
    limits = Limits() if limits is None else limits
    if not isinstance(requests, Sequence):
        raise ReplayError("a finite request sequence is required for resumable execution")
    schedule = requests
    initial = state
    transitions = []
    history = [metrics(state)]
    status, diagnostic = "schedule_exhausted", None
    try:
        limits.check(state)
    except LimitExceeded as exc:
        status, diagnostic = "incomplete", str(exc)
    for idx, request in enumerate(schedule):
        if status == "incomplete":
            break
        try:
            if idx >= limits.max_steps:
                raise LimitExceeded("step limit exceeded")
            following, transition = apply_rule(state, request, catalog=catalog, limits=limits)
        except LimitExceeded as exc:
            status, diagnostic = "incomplete", str(exc)
            break
        state = following
        transitions.append(transition)
        history.append(metrics(state))
    return {"schema": "zipper-contextual-replay-v1", "scope": "algebraic_catalog_replay",
            "source_sha256": _source_sha256(),
            "catalog_sha256": catalog.sha256,
            "initial": initial.to_json(), "final": state.to_json(), "transitions": transitions,
            "pending_requests": [x.to_json() for x in schedule[len(transitions):]],
            "status": status, "diagnostic": diagnostic, "limits": vars(limits), "metrics": history}


def resume(certificate: Mapping[str, Any], *, catalog: Catalog | None = None,
           limits: Limits | None = None) -> dict[str, Any]:
    """Verify saved work and continue its pending schedule under new limits."""
    check_certificate(certificate, catalog=catalog)
    return execute(RunState.from_json(certificate["final"]),
                   [Request.from_json(x) for x in certificate["pending_requests"]],
                   catalog=catalog, limits=limits)


def check_certificate(certificate: Mapping[str, Any], *, catalog: Catalog | None = None) -> bool:
    """Independently reinstantiate every saved local equality after JSON reload.

    Acceptance certifies contextual applications of the accepted catalog
    identities and exact coefficient propagation. It does not certify recorded
    control evidence, the identities' derivations, or the scheduler's coverage.
    """
    _fields(certificate, {"schema", "scope", "source_sha256", "catalog_sha256", "initial", "final", "transitions",
                          "pending_requests", "status", "diagnostic", "limits", "metrics"})
    if certificate["schema"] != "zipper-contextual-replay-v1" or certificate["scope"] != "algebraic_catalog_replay":
        raise ReplayError("unsupported certificate schema or scope")
    if certificate["source_sha256"] != _source_sha256():
        raise ReplayError("interpreter source hash mismatch")
    catalog = Catalog.bundled() if catalog is None else catalog
    if certificate["catalog_sha256"] != catalog.sha256:
        raise ReplayError("catalog source hash mismatch")
    limits = Limits(**certificate["limits"])
    initial = RunState.from_json(certificate["initial"])
    requests = [Request.from_json(x["request"]) for x in certificate["transitions"]]
    requests.extend(Request.from_json(x) for x in certificate["pending_requests"])
    expected = execute(initial, requests, catalog=catalog, limits=limits)
    if _json(expected) != _json(certificate):
        raise ReplayError("certificate differs from exact contextual replay")
    return True
