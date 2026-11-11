"""Independent Python replay with a bounded cache of complete local states.

Every cached value is computed by the Python reference transition function.
Keys retain full caller geometry, contexts, labels, resource limits and
trace setting. Coefficients are normalized to one and depth to zero, then restored
exactly. Unlike the producer cache, no outer context or endpoint label is removed. No native transition, cache key, result, or witness populates this
cache. Pickle is used only for bytes created within this process; external
certificates are never unpickled.
"""
from __future__ import annotations

from collections import Counter, OrderedDict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import hashlib
from pathlib import Path
import pickle
import sys
from threading import Lock
from typing import Any, TYPE_CHECKING, cast

from . import collected
from .interpreter import Limits
from .catalog import Catalog, CatalogError

if TYPE_CHECKING:
    from .native import NativeSession


def source_identity() -> dict[str, str]:
    """Bind the verifier and all Python reference implementation sources."""
    return {**collected._source(),
            "independent_replay.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "native.py": hashlib.sha256(Path(__file__).with_name("native.py").read_bytes()).hexdigest()}


class _ReplayCatalog(Catalog):
    """Clone frozen rule trees using trusted local bytes, preserving isolation."""

    def __init__(self, catalog: dict[str, Any], *, source_sha256: str | None = None) -> None:
        super().__init__(catalog, source_sha256=source_sha256)
        self._frozen_rules = {key: pickle.dumps(value, protocol=5)
                              for key, value in self._rules.items()}

    def get(self, rule_id: str) -> dict[str, Any]:
        try:
            return pickle.loads(self._frozen_rules[rule_id])
        except KeyError as error:
            raise CatalogError("Unknown rule id: " + rule_id) from error

    def select_rules(self, env: dict[str, int], *, prefix: str = "", suffix: str = "",
                     exclude_prefixes: tuple[str, ...] = ()) -> tuple[dict[str, Any], ...]:
        return tuple(self.get(key) for key in self.select_rule_ids(
            env, prefix=prefix, suffix=suffix, exclude_prefixes=exclude_prefixes))


class IndependentVerifier(collected.CollectedSession):
    """Retain independent full-state transition evaluations across inputs."""

    _catalog: _ReplayCatalog

    def __init__(self, *, cache_bytes: int = 256 * 1024**2,
                 cache_entries: int = 262144, rule_cache_bytes: int = 4 * 1024**2,
                 rule_cache_entries: int = 256) -> None:
        if any(type(x) is not int or x < 0 for x in
               (cache_bytes, cache_entries, rule_cache_bytes, rule_cache_entries)):
            raise collected.CollectedError("Verifier cache budgets must be nonnegative integers")
        self._catalog = cast(_ReplayCatalog, _ReplayCatalog.bundled())
        self._catalog.configure_rule_cache(max_bytes=rule_cache_bytes,
                                          max_entries=rule_cache_entries)
        self._max_bytes, self._max_entries = cache_bytes, cache_entries
        self._identity = source_identity()
        self._lock = Lock()
        self._cache: OrderedDict[bytes, tuple[bytes, int]] = OrderedDict()
        self.clear_cache()

    def clear_cache(self) -> None:
        """Drop all locally computed cached values and reset counters."""
        with self._lock:
            self._cache.clear()
            self._bytes = 0
            self._input = 0
            self._counts: Counter[str] = Counter()
            self._catalog.rule_cache_clear()

    def cache_info(self) -> dict[str, Any]:
        """Return cache costs separately from native producer telemetry."""
        with self._lock:
            return self._cache_info()

    def _cache_info(self) -> dict[str, Any]:
        fields = ("hits", "misses", "evaluations", "evictions", "oversized",
                  "cross_input_hits", "complement_hits", "peak_bytes", "peak_entries")
        return {**{key: self._counts[key] for key in fields},
                "entries": len(self._cache), "bytes": self._bytes,
                "max_entries": self._max_entries, "max_bytes": self._max_bytes,
                "byte_measure": "retained-key-and-payload-object-sizes",
                "resident_estimate_bytes": self._bytes + sys.getsizeof(self._cache)
                + sum(sys.getsizeof(entry) + sys.getsizeof(entry[1]) for entry in self._cache.values()),
                "frozen_rule_payload_bytes": sum(map(len, self._catalog._frozen_rules.values())),
                "rule_cache": self._catalog.rule_cache_info(),
                "mode": "independent-python-unfactored-unit-replay-v1"}

    def _execute_batch(
        self, selected: list[tuple[str, dict[str, Any]]], limits: Limits,
        trace: bool, workers: int, batch_size: int, pool: ProcessPoolExecutor | None,
    ) -> tuple[list[dict[str, Any]], ProcessPoolExecutor | None]:
        results = []
        for _, work in selected:
            # Dict insertion order may cause a harmless cache miss. It cannot
            # identify unequal states; bytes retain every input field and type.
            unit = dict(work, coefficient=[[0, 1]], depth=0)
            if "branch" in work:
                unit["branch"] = dict(work["branch"], coefficient=[[0, 1]])
            key = pickle.dumps((asdict(limits), trace, unit), protocol=5)
            stored = self._cache.get(key)
            if stored is not None:
                payload, origin = stored
                self._cache.move_to_end(key)
                self._counts["hits"] += 1
                self._counts["cross_input_hits"] += origin != self._input
                self._counts["complement_hits"] += work["phase"] == "extract"
                result = pickle.loads(payload)
            else:
                self._counts["misses"] += 1
                self._counts["evaluations"] += 1
                result = collected._transition(unit, limits, trace, self._catalog)
                if self._max_entries and self._max_bytes:
                    payload = pickle.dumps(result, protocol=5)
                    size = sys.getsizeof(key) + sys.getsizeof(payload)
                    if size > self._max_bytes:
                        self._counts["oversized"] += 1
                    else:
                        while self._cache and (len(self._cache) >= self._max_entries
                                               or self._bytes + size > self._max_bytes):
                            old_key, (old_payload, _) = self._cache.popitem(last=False)
                            self._bytes -= sys.getsizeof(old_key) + sys.getsizeof(old_payload)
                            self._counts["evictions"] += 1
                        self._cache[key] = (payload, self._input)
                        self._bytes += size
                        self._counts["peak_bytes"] = max(self._counts["peak_bytes"], self._bytes)
                        self._counts["peak_entries"] = max(self._counts["peak_entries"], len(self._cache))
            # Cache payloads precede scaling. All mutable returned values are
            # isolated from the retained bytes, including the cold miss path.
            if "error" not in result:
                for item in result["children"] + result["output"]:
                    collected._set_coefficient(item, collected.reference._multiply(
                        work["coefficient"], item["coefficient"]))
                for child in result["children"]:
                    child["depth"] = work["depth"] + 1
            results.append(result)
        return results, pool

    def check_native(self, session: NativeSession, certificate: dict[str, Any]) -> bool:
        """Recompute the entire certificate; producer outputs are comparison data only."""
        from .native import NativeError, _json_arrays, _strict_equal, reference_certificate

        with self._lock:
            session._validate_certificate_identity(certificate)
            if source_identity() != self._identity:
                raise NativeError("Independent verifier sources changed during the session")
            self._input += 1
            try:
                expected = self.reduce(
                    certificate["input_word"], n=certificate["n"],
                    limits=Limits(**certificate["limits"]), trace=certificate["trace_enabled"],
                    batch_size=certificate["execution"]["batch_size"], workers=1)
                actual = reference_certificate(certificate)
            except (KeyError, TypeError, ValueError) as error:
                raise NativeError("Malformed native certificate") from error
            if not _strict_equal(_json_arrays(expected), actual):
                raise NativeError("Native certificate differs from independent Python replay")
            return True
