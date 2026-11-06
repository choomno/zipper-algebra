"""Bounded reuse of intrinsic unit transitions, including complement steps.

The active diagram and all of its control data are retained in cache keys.
Already extracted outer components, endpoint label values, coefficients and
diagnostic depth are supplied afresh by each caller.  Cached results describe
one local transition, not a completed complement or a global formal sum.
Only successful transitions are retained. Limits and trace policy are keyed.

The LRU stores private immutable pickle bytes, bounded by both entry count and the
combined encoded key/payload bytes. These byte counts are not process RSS.
Sources and catalog must remain fixed for a session's lifetime.
The independent verifier lazily retains a second parsed catalog, but keeps
both of its result caches disabled; process RSS includes that catalog.
"""

from __future__ import annotations

from collections import Counter, OrderedDict
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from dataclasses import asdict
import hashlib
import pickle
import sys
from pathlib import Path
from multiprocessing import get_context
from typing import Any, Iterable

from . import collected
from . import algorithm as reference
from .interpreter import Limits


def _source() -> dict[str, str]:
    return {
        **collected._source(),
        "memoized.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def _initialize_worker(rule_bytes: int, rule_entries: int) -> None:
    collected._initialize_worker()
    assert collected._CATALOG is not None
    collected._CATALOG.configure_rule_cache(
        max_bytes=rule_bytes, max_entries=rule_entries
    )


def _intrinsic(work: dict[str, Any]) -> dict[str, Any]:
    unit = dict(work)
    if "branch" in unit:
        unit["branch"] = dict(unit["branch"])
    unit["contexts"] = []
    unit["bottom_labels"] = list(range(1, work["rank"] + 1))
    unit["top_labels"] = list(range(1, work["rank"] + 1))
    unit["depth"] = 0
    collected._set_coefficient(unit, [[0, 1]])
    return unit


def _relabel_record(
    record: dict[str, Any], work: dict[str, Any], seen: set[int]
) -> None:
    if id(record) in seen:
        return
    seen.add(id(record))
    for side in ("bottom", "top"):
        name = "original_" + side + "_ports"
        record[name] = [work[side + "_labels"][p - 1] for p in record[name]]


def _contexts(
    contexts: list[dict[str, Any]], work: dict[str, Any], seen: set[int]
) -> list[dict[str, Any]]:
    for context in contexts:
        _relabel_record(context["record"], work, seen)
    return deepcopy(work["contexts"]) + contexts


def _rebase(template: dict[str, Any], work: dict[str, Any]) -> dict[str, Any]:
    result = deepcopy(template)
    if "error" in result:
        return result
    seen: set[int] = set()
    for child in result["children"]:
        child["contexts"] = _contexts(child["contexts"], work, seen)
        for side in ("bottom", "top"):
            child[side + "_labels"] = [
                work[side + "_labels"][p - 1] for p in child[side + "_labels"]
            ]
        child["depth"] = work["depth"] + 1
        collected._set_coefficient(
            child, reference._multiply(work["coefficient"], child["coefficient"])
        )
    for output in result["output"]:
        output["components"] = _contexts(output["components"], work, seen)
        output["word"] = reference._encoded(
            reference._lift(reference._word(output["word"]), work["contexts"])
        )
        collected._set_coefficient(
            output, reference._multiply(work["coefficient"], output["coefficient"])
        )
    event = result["unit_evidence"]
    if "component" in event:
        _relabel_record(event["component"], work, seen)
    return result


class MemoizedSession(collected.CollectedSession):
    """A collected executor with a bounded, per-session transition LRU."""

    def __init__(
        self,
        *,
        cache_bytes: int = 64 * 1024**2,
        cache_entries: int = 4096,
        rule_cache_bytes: int = 4 * 1024**2,
        rule_cache_entries: int = 256,
    ) -> None:
        if any(
            type(x) is not int or x < 0
            for x in (cache_bytes, cache_entries, rule_cache_bytes, rule_cache_entries)
        ):
            raise collected.CollectedError(
                "cache byte and entry budgets must be nonnegative integers"
            )
        super().__init__()
        self._catalog.configure_rule_cache(
            max_bytes=rule_cache_bytes, max_entries=rule_cache_entries
        )
        self._max_bytes, self._max_entries = cache_bytes, cache_entries
        self._rule_bytes, self._rule_entries = rule_cache_bytes, rule_cache_entries
        self._cache: OrderedDict[str, bytes] = OrderedDict()
        self._checker: _UnfactoredReplay | None = None
        self._input = 0
        self._identity = hashlib.sha256(
            reference._json(
                {"sources": _source(), "catalog": self._catalog.sha256}
            ).encode()
        ).hexdigest()
        self.clear_cache()

    def _source_identity(self) -> dict[str, str]:
        return _source()

    def _new_pool(self, workers: int, batch_size: int) -> ProcessPoolExecutor:
        return ProcessPoolExecutor(
            max_workers=min(workers, batch_size),
            initializer=_initialize_worker,
            initargs=(self._rule_bytes, self._rule_entries),
            mp_context=get_context("spawn"),
        )

    def clear_cache(self) -> None:
        """Release all retained templates and reset cache telemetry."""
        self._cache.clear()
        self._bytes = 0
        self._resident_payload_bytes = 0
        self._stats: Counter[str] = Counter()
        self._phase_hits: Counter[str] = Counter()
        self._phase_misses: Counter[str] = Counter()
        self._rule_evaluations: Counter[str] = Counter()
        self._rule_hits: Counter[str] = Counter()
        self._rule_batch_reuses: Counter[str] = Counter()
        self._catalog.rule_cache_clear()

    def cache_info(self) -> dict[str, Any]:
        """Runtime-only counters; encoded storage is distinct from resident RAM."""
        fields = (
            "hits",
            "misses",
            "evictions",
            "oversized",
            "batch_reuses",
            "complement_hits",
            "cross_input_hits",
            "evaluations",
            "peak_bytes",
            "peak_entries",
        )
        rule_cache = self._catalog.rule_cache_info()
        resident = sys.getsizeof(self._cache) + self._resident_payload_bytes
        return {
            **{name: self._stats[name] for name in fields},
            "entries": len(self._cache),
            "bytes": self._bytes,
            "max_bytes": self._max_bytes,
            "max_entries": self._max_entries,
            "byte_measure": "encoded-key-and-payload",
            "phase_hits": dict(self._phase_hits),
            "phase_misses": dict(self._phase_misses),
            "actual_unit_evaluations": dict(self._rule_evaluations),
            "unit_lru_hits": dict(self._rule_hits),
            "unit_batch_reuses": dict(self._rule_batch_reuses),
            "evaluation_scope": "attempted-local-transitions-including-rolled-back-batches",
            "resident_estimate_bytes": resident,
            "rule_cache": rule_cache,
            "rule_cache_scope": "coordinator-counters; identical-budget-separate-worker-caches",
            "total_encoded_bytes": self._bytes + rule_cache["retained_bytes"],
            "total_resident_estimate_bytes": resident
            + rule_cache["resident_estimate_bytes"],
        }

    @staticmethod
    def _rule(result: dict[str, Any]) -> str:
        event = result.get("unit_evidence", {})
        return str(
            event.get(
                "rule_id",
                event.get("catalog_transition", {})
                .get("request", {})
                .get("rule_id", "refused-local-transition"),
            )
        )

    def _key(self, work: dict[str, Any], limits: Limits, trace: bool) -> str:
        return reference._json(
            {
                "identity": self._identity,
                "limits": asdict(limits),
                "trace": trace,
                "active": collected.continuation_key(work),
            }
        )

    def _lookup(self, key: str, work: dict[str, Any]) -> dict[str, Any] | None:
        payload = self._cache.get(key)
        if payload is None:
            self._stats["misses"] += 1
            self._phase_misses[work["phase"]] += 1
            return None
        self._cache.move_to_end(key)
        saved = pickle.loads(payload)
        self._stats["hits"] += 1
        self._phase_hits[work["phase"]] += 1
        self._stats["cross_input_hits"] += int(saved["input"] != self._input)
        self._stats["complement_hits"] += bool(work["contexts"])
        self._rule_hits[self._rule(saved["result"])] += 1
        return saved["result"]

    def _remember(self, key: str, result: dict[str, Any]) -> None:
        if not self._max_entries or not self._max_bytes or "error" in result:
            return
        payload = pickle.dumps({"input": self._input, "result": result}, protocol=5)
        cost = len(key.encode()) + len(payload)
        if cost > self._max_bytes:
            self._stats["oversized"] += 1
            return
        while self._cache and (
            len(self._cache) >= self._max_entries
            or self._bytes + cost > self._max_bytes
        ):
            old_key, old_payload = self._cache.popitem(last=False)
            self._bytes -= len(old_key.encode()) + len(old_payload)
            self._resident_payload_bytes -= sys.getsizeof(old_key) + sys.getsizeof(
                old_payload
            )
            self._stats["evictions"] += 1
        self._cache[key] = payload
        self._bytes += cost
        self._resident_payload_bytes += sys.getsizeof(key) + sys.getsizeof(payload)
        self._stats["peak_bytes"] = max(self._stats["peak_bytes"], self._bytes)
        self._stats["peak_entries"] = max(self._stats["peak_entries"], len(self._cache))

    def _execute_batch(
        self,
        selected: list[tuple[str, dict[str, Any]]],
        limits: Limits,
        trace: bool,
        workers: int,
        batch_size: int,
        pool: ProcessPoolExecutor | None,
    ) -> tuple[list[dict[str, Any]], ProcessPoolExecutor | None]:
        templates: dict[str, dict[str, Any] | None] = {}
        keys = []
        missing = []
        aliases = []
        for _, work in selected:
            active = _intrinsic(work)
            key = self._key(active, limits, trace)
            keys.append(key)
            if key in templates:
                self._stats["batch_reuses"] += 1
                self._stats["complement_hits"] += bool(work["contexts"])
                aliases.append(key)
                continue
            result = self._lookup(key, work)
            templates[key] = result
            if result is None:
                missing.append((key, active))
        if missing:
            evaluated, pool = super()._execute_batch(
                missing, limits, trace, workers, batch_size, pool
            )
            self._stats["evaluations"] += len(missing)
            for (key, _), result in zip(missing, evaluated, strict=True):
                templates[key] = result
                self._rule_evaluations[self._rule(result)] += 1
                self._remember(key, result)
        for key in aliases:
            template = templates[key]
            assert template is not None
            self._rule_batch_reuses[self._rule(template)] += 1
        results = []
        for (_, work), key in zip(selected, keys, strict=True):
            template = templates[key]
            assert template is not None
            results.append(_rebase(template, work))
        return results, pool

    def reduce(
        self,
        word: Iterable[int],
        *,
        n: int,
        limits: Limits | None = None,
        trace: bool = False,
        workers: int = 1,
        batch_size: int = 64,
    ) -> dict[str, Any]:
        self._input += 1
        result = super().reduce(
            word,
            n=n,
            limits=limits,
            trace=trace,
            workers=workers,
            batch_size=batch_size,
        )
        result["schema"] = "zipper-memoized-run-v1"
        result["scope"] = "literal_word_memoized_collected_spanning_procedure"
        result["execution"]["memoization"] = "intrinsic-unit-transition-v1"
        return result

    def check(self, certificate: dict[str, Any]) -> bool:
        """Replay unfactored caller states with both cache layers disabled."""
        if self._checker is None:
            self._checker = _UnfactoredReplay()
        return collected.CollectedSession.check(
            self._checker,
            certificate,
        )


class _UnfactoredReplay(MemoizedSession):
    def __init__(self) -> None:
        super().__init__(
            cache_bytes=0, cache_entries=0, rule_cache_bytes=0, rule_cache_entries=0
        )

    def _execute_batch(
        self,
        selected: list[tuple[str, dict[str, Any]]],
        limits: Limits,
        trace: bool,
        workers: int,
        batch_size: int,
        pool: ProcessPoolExecutor | None,
    ) -> tuple[list[dict[str, Any]], ProcessPoolExecutor | None]:
        return collected.CollectedSession._execute_batch(
            self, selected, limits, trace, workers, batch_size, pool
        )


def reduce_word(
    word: Iterable[int],
    *,
    n: int,
    limits: Limits | None = None,
    trace: bool = False,
    workers: int = 1,
    batch_size: int = 64,
    cache_bytes: int = 64 * 1024**2,
    cache_entries: int = 4096,
    rule_cache_bytes: int = 4 * 1024**2,
    rule_cache_entries: int = 256,
) -> dict[str, Any]:
    return MemoizedSession(
        cache_bytes=cache_bytes,
        cache_entries=cache_entries,
        rule_cache_bytes=rule_cache_bytes,
        rule_cache_entries=rule_cache_entries,
    ).reduce(
        word, n=n, limits=limits, trace=trace, workers=workers, batch_size=batch_size
    )


def check_certificate(certificate: dict[str, Any]) -> bool:
    return collected.CollectedSession.check(_UnfactoredReplay(), certificate)


def resume(
    certificate: dict[str, Any],
    *,
    limits: Limits | None = None,
    workers: int = 1,
    cache_bytes: int = 64 * 1024**2,
    cache_entries: int = 4096,
    rule_cache_bytes: int = 4 * 1024**2,
    rule_cache_entries: int = 256,
) -> dict[str, Any]:
    session = MemoizedSession(
        cache_bytes=cache_bytes,
        cache_entries=cache_entries,
        rule_cache_bytes=rule_cache_bytes,
        rule_cache_entries=rule_cache_entries,
    )
    session.check(certificate)
    old = Limits(**certificate["limits"])
    following = old if limits is None else limits
    if any(getattr(following, name) < value for name, value in asdict(old).items()):
        raise collected.CollectedError("continuation limits must not decrease")
    return session.reduce(
        certificate["input_word"],
        n=certificate["n"],
        limits=following,
        trace=certificate["trace_enabled"],
        workers=workers,
        batch_size=certificate["execution"]["batch_size"],
    )
