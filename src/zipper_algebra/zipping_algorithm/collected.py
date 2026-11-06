"""Collected, deterministic batches of independent zipper transitions.

Only identical continuation states are combined.  The scalar coefficient
and the diagnostic maximum depth are excluded from their identity; every
other control, boundary, label and extraction field remains in the key.
Workers execute unit-coefficient transitions and scale their resulting edges.
The coordinator collects the formal sum and commits a whole batch.

``max_steps`` counts committed term transitions.  A batch which would cross
that bound is refused in its entirety.  Frontier term and accumulated scalar
bounds apply after collection; primitive transitions retain word, expansion
and unit-scalar bit bounds, with local term/support floors of four/two for
the verified local rules.
Worker count affects execution only, and is absent from the certificate.
Execution requires the bound package sources to remain fixed during a run.
"""

from __future__ import annotations

from collections import Counter, OrderedDict
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from copy import deepcopy
from dataclasses import asdict, replace
import hashlib
from itertools import islice
from multiprocessing import get_context
from pathlib import Path
from typing import Any, Iterable

from . import algorithm as reference
from .catalog import Catalog
from .interpreter import Limits, ReplayError


class CollectedError(reference.AlgorithmError):
    """Invalid collected-executor input or replay certificate."""


_CATALOG: Catalog | None = None
_EXPECTED = (
    reference.LimitExceeded,
    reference.ControlObligation,
    reference.upward.ControlError,
    reference.RibbonError,
    reference.FrameError,
    reference.TurnaroundError,
    reference.TransportError,
    reference.CatalogError,
)


def continuation_key(item: dict[str, Any], *, pending: bool = True) -> str:
    """Canonical identity of a complete continuation, or a completed term."""
    payload = dict(item)
    del payload["coefficient"]
    if pending:
        payload.pop("depth", None)
        if "branch" in payload:
            payload["branch"] = dict(payload["branch"])
            del payload["branch"]["coefficient"]
    return reference._json(payload)


def _set_coefficient(item: dict[str, Any], coefficient: list[list[int]]) -> None:
    item["coefficient"] = deepcopy(coefficient)
    if "branch" in item:
        item["branch"]["coefficient"] = deepcopy(coefficient)


def _add(
    left: Iterable[Iterable[int]], right: Iterable[Iterable[int]]
) -> list[list[int]]:
    values: Counter[int] = Counter()
    for exponent, scalar in (*left, *right):
        values[exponent] += scalar
    return [[q, c] for q, c in sorted(values.items()) if c]


def _stage(
    current: dict[str, dict[str, Any]],
    incoming: Iterable[dict[str, Any]],
    *,
    pending: bool,
    removed: set[str],
) -> tuple[dict[str, dict[str, Any] | None], int, int]:
    """Collect affected keys without changing the predecessor maps.

    Merge counts count additions to nonzero accumulators. Zero counts count
    keys absent at batch end; cancellation followed by re-addition is retained.
    """
    updates: dict[str, dict[str, Any] | None] = {}
    merged = 0
    for item in incoming:
        key = continuation_key(item, pending=pending)
        previous = updates.get(key)
        if key not in updates and key not in removed:
            previous = current.get(key)
        combined = deepcopy(item)
        if previous is not None:
            merged += 1
            coefficient = _add(previous["coefficient"], item["coefficient"])
            if pending:
                combined["depth"] = max(previous["depth"], item["depth"])
        else:
            coefficient = item["coefficient"]
        _set_coefficient(combined, coefficient)
        updates[key] = combined if coefficient else None
    zeroes = sum(item is None for item in updates.values())
    return updates, merged, zeroes


def collect_terms(
    items: Iterable[dict[str, Any]], *, pending: bool = True
) -> list[dict[str, Any]]:
    """Add exact Laurent coefficients of identical saved terms; remove zeros."""
    updates, _, _ = _stage({}, items, pending=pending, removed=set())
    return [item for item in updates.values() if item is not None]


def _initialize_worker() -> None:
    global _CATALOG
    _CATALOG = Catalog.bundled()


def _transition(
    work: dict[str, Any], limits: Limits, trace: bool, catalog: Catalog
) -> dict[str, Any]:
    unit = deepcopy(work)
    _set_coefficient(unit, [[0, 1]])
    local = replace(
        limits,
        max_terms=max(4, limits.max_terms),
        max_coefficient_support=max(2, limits.max_coefficient_support),
    )
    try:
        added: list[dict[str, Any]] = []
        phase = unit["phase"]
        if phase == "single":
            children, event = reference._single_step(unit, catalog, local)
        elif phase in ("zip_up", "zip_down"):
            children, event = reference._zip_step(unit, catalog, local, trace)
        elif phase == "layer":
            children, event = reference._layer_step(unit, local)
        else:
            children, added, event = reference._extract_step(
                unit, limits=local, trace=trace
            )
    except _EXPECTED as exc:
        bounded = isinstance(exc, (reference.LimitExceeded, reference.RibbonLimit)) or (
            isinstance(exc, reference.CatalogError) and "budget" in str(exc)
        )
        return {
            "error": "incomplete" if bounded else "control_obligation",
            "message": str(exc),
        }
    for item in children + added:
        _set_coefficient(
            item, reference._multiply(work["coefficient"], item["coefficient"])
        )
    for child in children:
        child["depth"] = work["depth"] + 1
    return {"children": children, "output": added, "unit_evidence": event}


def _worker(job: tuple[dict[str, Any], Limits, bool]) -> dict[str, Any]:
    work, limits, trace = job
    # Initialized once per process, reused for all its bounded batch jobs.
    assert _CATALOG is not None
    return _transition(work, limits, trace, _CATALOG)


def _check_items(
    items: Iterable[dict[str, Any]], limits: Limits, *, pending: bool
) -> None:
    for item in items:
        size = reference._size(item) if pending else len(item["word"])
        if size > limits.max_word_size:
            raise reference.LimitExceeded("word size limit exceeded")
        if len(item["coefficient"]) > limits.max_coefficient_support:
            raise reference.LimitExceeded("coefficient support limit exceeded")
        if any(
            abs(c).bit_length() > limits.max_integer_bits
            for _, c in item["coefficient"]
        ):
            raise reference.LimitExceeded("coefficient integer bit limit exceeded")


def _commit(
    current: dict[str, Any], updates: dict[str, Any], removed: Iterable[str] = ()
) -> None:
    for key in removed:
        del current[key]
    for key, value in updates.items():
        if value is None:
            current.pop(key, None)
        else:
            current[key] = value


def _count_after(
    current: dict[str, Any], updates: dict[str, Any], removed: set[str]
) -> int:
    return (
        len(current)
        - len(removed)
        + sum(
            int(value is not None) - int(key in current and key not in removed)
            for key, value in updates.items()
        )
    )


def _source() -> dict[str, str]:
    return {
        **reference._source(),
        "collected.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


class CollectedSession:
    """Reuse local identities; isolate each collected input calculation."""

    def __init__(self) -> None:
        self._catalog = Catalog.bundled()

    def _source_identity(self) -> dict[str, str]:
        return _source()

    def _new_pool(self, workers: int, batch_size: int) -> ProcessPoolExecutor:
        return ProcessPoolExecutor(
            max_workers=min(workers, batch_size),
            initializer=_initialize_worker,
            mp_context=get_context("spawn"),
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
        if workers == 1 or len(selected) == 1:
            return [
                _transition(item, limits, trace, self._catalog) for _, item in selected
            ], pool
        if pool is None:
            pool = self._new_pool(workers, batch_size)
        try:
            results = list(
                pool.map(
                    _worker,
                    ((item, limits, trace) for _, item in selected),
                    chunksize=max(1, len(selected) // (4 * min(workers, batch_size))),
                )
            )
        except BaseException:
            pool.shutdown(wait=True, cancel_futures=True)
            raise
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
        if (
            type(n) is not int
            or n < 1
            or type(trace) is not bool
            or type(workers) is not int
            or workers < 1
            or type(batch_size) is not int
            or batch_size < 1
        ):
            raise CollectedError(
                "positive integer rank, workers and batch size, and boolean trace required"
            )
        if not isinstance(word, (list, tuple)) or any(
            type(x) is not int or not 1 <= abs(x) < n for x in word
        ):
            raise CollectedError(
                "a literal signed-generator list or tuple in the stated rank is required"
            )
        limits = Limits() if limits is None else limits
        initial = reference._single(
            {
                "rank": n,
                "coefficient": [[0, 1]],
                "contexts": [],
                "world_phi": False,
                "boundary_left": [],
                "boundary_right": [],
                "boundary_left_factors": [],
                "boundary_right_factors": [],
                "depth": 0,
                "bottom_labels": list(range(1, n + 1)),
                "top_labels": list(range(1, n + 1)),
            },
            word,
        )
        frontier = OrderedDict(((continuation_key(initial), initial),))
        output: dict[str, dict[str, Any]] = {}
        transitions: list[dict[str, Any]] = []
        counts: Counter[str] = Counter()
        metrics = {
            "steps": 0,
            "attempted_steps": 0,
            "batches": 0,
            "max_terms": 1,
            "max_local_word_tokens": len(word) + 1,
            "max_coefficient_support": 1,
            "max_integer_bits": 1,
            "q_exponent_min": 0,
            "q_exponent_max": 0,
            "maximum_depth": 0,
            "max_branch_count": 0,
            "merged_terms": 0,
            "zero_terms": 0,
        }
        status, diagnostic = "complete", None
        pool: ProcessPoolExecutor | None = None
        try:
            while frontier:
                selected = list(islice(frontier.items(), batch_size))
                before = selected[0][1]
                try:
                    if metrics["steps"] + len(selected) > limits.max_steps:
                        raise reference.LimitExceeded(
                            "step limit would be exceeded by the next atomic batch"
                        )
                    _check_items((item for _, item in selected), limits, pending=True)
                    metrics["attempted_steps"] += len(selected)
                    results, pool = self._execute_batch(
                        selected, limits, trace, workers, batch_size, pool
                    )
                    failure = next(
                        (
                            (item, result)
                            for (_, item), result in zip(selected, results)
                            if "error" in result
                        ),
                        None,
                    )
                    if failure is not None:
                        before, result = failure
                        status = result["error"]
                        diagnostic = {
                            "message": result["message"],
                            "phase": before["phase"],
                            "state": deepcopy(before),
                            "batch_rolled_back": True,
                        }
                        break
                    removed = {key for key, _ in selected}
                    children = [
                        item for result in results for item in result["children"]
                    ]
                    added = [item for result in results for item in result["output"]]
                    pending_updates, pm, pz = _stage(
                        frontier, children, pending=True, removed=removed
                    )
                    output_updates, om, oz = _stage(
                        output, added, pending=False, removed=set()
                    )
                    term_count = _count_after(
                        frontier, pending_updates, removed
                    ) + _count_after(output, output_updates, set())
                    if term_count > limits.max_terms:
                        raise reference.LimitExceeded("collected term limit exceeded")
                    _check_items(
                        (x for x in pending_updates.values() if x is not None),
                        limits,
                        pending=True,
                    )
                    _check_items(
                        (x for x in output_updates.values() if x is not None),
                        limits,
                        pending=False,
                    )
                except reference.LimitExceeded as exc:
                    status = "incomplete"
                    diagnostic = {
                        "message": str(exc),
                        "phase": "batch",
                        "batch_representative": deepcopy(before),
                        "batch_rolled_back": True,
                    }
                    break
                _commit(frontier, pending_updates, removed)
                _commit(output, output_updates)
                metrics["steps"] += len(selected)
                metrics["batches"] += 1
                metrics["merged_terms"] += pm + om
                metrics["zero_terms"] += pz + oz
                metrics["max_terms"] = max(metrics["max_terms"], term_count)
                for (_, parent), result in zip(selected, results):
                    event = result["unit_evidence"]
                    name = event.get(
                        "rule_id",
                        event.get("catalog_transition", {})
                        .get("request", {})
                        .get("rule_id"),
                    )
                    counts[str(name)] += 1
                    metrics["maximum_depth"] = max(
                        metrics["maximum_depth"], parent["depth"] + 1
                    )
                    metrics["max_branch_count"] = max(
                        metrics["max_branch_count"],
                        len(result["children"]) + len(result["output"]),
                    )
                for pending, updates in (
                    (True, pending_updates),
                    (False, output_updates),
                ):
                    for item in updates.values():
                        if item is None:
                            continue
                        size = reference._size(item) if pending else len(item["word"])
                        metrics["max_local_word_tokens"] = max(
                            metrics["max_local_word_tokens"], size
                        )
                        metrics["max_coefficient_support"] = max(
                            metrics["max_coefficient_support"], len(item["coefficient"])
                        )
                        for q, c in item["coefficient"]:
                            metrics["q_exponent_min"] = min(
                                metrics["q_exponent_min"], q
                            )
                            metrics["q_exponent_max"] = max(
                                metrics["q_exponent_max"], q
                            )
                            metrics["max_integer_bits"] = max(
                                metrics["max_integer_bits"], abs(c).bit_length()
                            )
                if trace:
                    transitions.append(
                        {
                            "batch": metrics["batches"] - 1,
                            "predecessors": [deepcopy(item) for _, item in selected],
                            "rewrites": results,
                            "pending_updates": deepcopy(pending_updates),
                            "output_updates": deepcopy(output_updates),
                        }
                    )
        except BrokenProcessPool as exc:
            raise CollectedError(
                f"worker process failed after {metrics['batches']} committed batches "
                f"and {metrics['steps']} committed term transitions"
            ) from exc
        finally:
            if pool is not None:
                pool.shutdown(wait=True, cancel_futures=True)
        return {
            "schema": "zipper-collected-run-v1",
            "scope": "literal_word_collected_spanning_procedure",
            "n": n,
            "input_word": list(word),
            "coefficient_ring": "Z[q,q^-1]",
            "source_sha256": self._source_identity(),
            "catalog_sha256": self._catalog.sha256,
            "trace_enabled": trace,
            "limits": asdict(limits),
            "execution": {
                "batch_size": batch_size,
                "schedule": "stable-prefix-atomic-v1",
                "collection": "full-continuation-except-scalar-and-depth-v1",
                "local_term_floor": 4,
                "local_support_floor": 2,
            },
            "status": status,
            "diagnostic": diagnostic,
            "pending": list(frontier.values()),
            "output": list(output.values()),
            "transitions": transitions,
            "metrics": {
                **metrics,
                "rule_applications": dict(counts),
                "output_terms": len(output),
                "pending_terms": len(frontier),
            },
        }

    def check(self, certificate: dict[str, Any]) -> bool:
        try:
            expected = self.reduce(
                certificate["input_word"],
                n=certificate["n"],
                limits=Limits(**certificate["limits"]),
                trace=certificate["trace_enabled"],
                batch_size=certificate["execution"]["batch_size"],
            )
        except (KeyError, TypeError, ReplayError) as exc:
            raise CollectedError("malformed collected certificate") from exc
        if reference._json(expected) != reference._json(certificate):
            raise CollectedError(
                "collected certificate differs from deterministic replay"
            )
        return True


def reduce_word(
    word: Iterable[int],
    *,
    n: int,
    limits: Limits | None = None,
    trace: bool = False,
    workers: int = 1,
    batch_size: int = 64,
) -> dict[str, Any]:
    return CollectedSession().reduce(
        word, n=n, limits=limits, trace=trace, workers=workers, batch_size=batch_size
    )


def check_certificate(certificate: dict[str, Any]) -> bool:
    return CollectedSession().check(certificate)


def resume(
    certificate: dict[str, Any], *, limits: Limits | None = None, workers: int = 1
) -> dict[str, Any]:
    """Verify a checkpoint, then replay and extend its deterministic batches."""
    session = CollectedSession()
    session.check(certificate)
    old = Limits(**certificate["limits"])
    following = old if limits is None else limits
    if any(getattr(following, name) < value for name, value in asdict(old).items()):
        raise CollectedError("continuation limits must not decrease")
    return session.reduce(
        certificate["input_word"],
        n=certificate["n"],
        limits=following,
        trace=certificate["trace_enabled"],
        workers=workers,
        batch_size=certificate["execution"]["batch_size"],
    )
