"""Whole-engine timings and independent exact replay, consolidated in SQLite.

Run from a checkout with ``PYTHONPATH=src:.``. Native call time includes request
encoding, IPC, execution, final certificate encoding, and response decoding.
Each timed producer is followed by a separate uncached Python replay. Both
producer and verification times are retained, as is their sum. Certificates
are compressed SQLite blobs; no per-word files are emitted.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import resource
import sqlite3
import time
from typing import Any
import zlib

from zipper_algebra.zipping_algorithm import collected, memoized
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import _strict_equal


CASES = {
    "b3_alternating_l10": (3, [1, -2] * 5),
    "b5_positive_l6": (5, [1, 2, 3, 4, 1, 2]),
    "b8_positive_l7": (8, [1, 2, 3, 4, 5, 6, 7]),
    "b5_positive_l10": (5, [1, 2, 3, 4, 1, 2, 3, 4, 1, 2]),
    "b8_positive_l10": (8, [1, 2, 3, 4, 5, 6, 7, 1, 2, 3]),
    "b3_negative_l10": (3, [-1, -2] * 5),
    "b2_unreduced_l10": (2, [1, -1] * 5),
    "smoke": (3, [1, 2]),
}
BACKENDS = ("native", "python_memoized", "python_collected")


def encoded(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def as_arrays(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return [as_arrays(item) for item in value]
    if isinstance(value, dict):
        return {key: as_arrays(item) for key, item in value.items()}
    return value


def digest(value: Any) -> str:
    return hashlib.sha256(encoded(value)).hexdigest()


def reference_form(certificate: dict[str, Any], backend: str) -> dict[str, Any]:
    """Remove implementation provenance only; retain all mathematical fields."""
    if backend == "native":
        from zipper_algebra.zipping_algorithm.native import reference_certificate
        return reference_certificate(certificate)
    result = dict(certificate)
    if backend == "python_memoized":
        result["schema"] = "zipper-collected-run-v1"
        result["scope"] = "literal_word_collected_spanning_procedure"
        result["source_sha256"] = collected._source()
        result["execution"] = dict(result["execution"])
        result["execution"].pop("memoization")
    return result


def source_identity() -> dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    native = root / "native" / "zipper_engine"
    return {
        "python": memoized._source(),
        "native": {
            str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(native.rglob("*"))
            if path.is_file() and path.suffix in (".cpp", ".hpp", ".h")
        },
        "harness": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "native_wrapper": hashlib.sha256(
            (root / "src/zipper_algebra/zipping_algorithm/native.py").read_bytes()).hexdigest(),
    }


def rss_snapshot(pid: int | str = "self") -> dict[str, int]:
    result = {}
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            key = line.split(":", 1)[0]
            if key in ("VmRSS", "VmHWM", "VmSize"):
                result[key + "_kib"] = int(line.split()[1])
    except OSError:
        pass
    return result


def initialize_database(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=WAL")
    connection.executescript("""
        CREATE TABLE IF NOT EXISTS campaigns (
            id INTEGER PRIMARY KEY, started_utc TEXT NOT NULL,
            configuration_json TEXT NOT NULL, sources_json TEXT NOT NULL,
            setup_json TEXT NOT NULL, status TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS measurements (
            id INTEGER PRIMARY KEY, campaign_id INTEGER NOT NULL,
            round INTEGER NOT NULL, position INTEGER NOT NULL,
            backend TEXT NOT NULL, case_name TEXT NOT NULL,
            n INTEGER NOT NULL, word_json TEXT NOT NULL,
            producer_seconds REAL NOT NULL, verification_seconds REAL NOT NULL,
            verified_total_seconds REAL NOT NULL, status TEXT NOT NULL,
            certificate_sha256 TEXT NOT NULL, mathematical_sha256 TEXT NOT NULL,
            certificate_zlib BLOB NOT NULL, metrics_json TEXT NOT NULL,
            telemetry_json TEXT NOT NULL, verified INTEGER NOT NULL CHECK(verified=1),
            FOREIGN KEY(campaign_id) REFERENCES campaigns(id),
            UNIQUE(campaign_id, round, backend, case_name)
        );
    """)
    connection.commit()
    return connection


def run(
    database: Path,
    *,
    cases: list[str],
    rounds: int,
    backends: list[str],
    limits: Limits,
    batch_size: int,
    trace: bool = False,
    cache_policy: str = "cold_each_input",
) -> int:
    if rounds < 1 or not cases or not backends or len(set(backends)) != len(backends):
        raise ValueError("positive rounds and nonempty distinct backends required")
    if any(name not in CASES for name in cases) or any(name not in BACKENDS for name in backends):
        raise ValueError("unknown workload or backend")
    if cache_policy not in ("cold_each_input", "retained"):
        raise ValueError("unknown cache policy")
    sources = source_identity()
    configuration = {
        "cases": cases, "rounds": rounds, "backends": backends,
        "limits": asdict(limits), "batch_size": batch_size, "trace": trace,
        "ordering": "cyclic rotations followed by reversed rotations; each backend occupies each position per block",
        "native_timing": "entire reduce call including boundary serialization and IPC",
        "verification": "fresh full-state Python replay with both result caches disabled",
        "cache_policy": cache_policy,
        "cache_budgets": {"max_entries": 262144, "max_bytes": 256 * 1024**2,
                          "byte_measure": "native estimated typed storage; Python encoded key/payload"},
        "cache_clear_timing": "recorded separately before producer timing for cold_each_input",
        "memory": "Linux proc status KiB; parent peak includes preceding runs and verifier",
    }
    setup: dict[str, Any] = {}
    sessions: dict[str, Any] = {}
    connection = initialize_database(database)
    campaign = None
    try:
        for backend in backends:
            start = time.perf_counter()
            if backend == "native":
                from zipper_algebra.zipping_algorithm.native import NativeSession
                session = NativeSession(cache_entries=262144, cache_bytes=256 * 1024**2)
            elif backend == "python_memoized":
                session = memoized.MemoizedSession(cache_bytes=256 * 1024**2, cache_entries=262144)
            else:
                session = collected.CollectedSession()
            sessions[backend] = session
            setup[backend] = {"seconds": time.perf_counter() - start}
            if backend == "native":
                setup[backend]["build"] = session.metadata
        cursor = connection.execute(
            "INSERT INTO campaigns(started_utc,configuration_json,sources_json,setup_json,status) VALUES(?,?,?,?,?)",
            (datetime.now(timezone.utc).isoformat(), encoded(configuration).decode(),
             encoded(sources).decode(), encoded(setup).decode(), "running"),
        )
        campaign = int(cursor.lastrowid)
        connection.commit()
        orderings = [tuple(backends[i:] + backends[:i]) for i in range(len(backends))]
        if len(backends) > 2:
            reverse = list(reversed(backends))
            orderings += [tuple(reverse[i:] + reverse[:i]) for i in range(len(reverse))]
        for round_index in range(rounds):
            order = orderings[round_index % len(orderings)]
            for case_name in cases:
                n, word = CASES[case_name]
                for position, backend in enumerate(order):
                    if source_identity() != sources:
                        raise RuntimeError("source changed during benchmark; refusing mixed-source measurement")
                    session = sessions[backend]
                    cache_clear_start = time.perf_counter()
                    if cache_policy == "cold_each_input":
                        if hasattr(session, "clear_cache"):
                            session.clear_cache()
                        elif hasattr(session, "_catalog"):
                            session._catalog.rule_cache_clear()
                    cache_clear_seconds = time.perf_counter() - cache_clear_start
                    before = rss_snapshot()
                    start = time.perf_counter()
                    certificate = session.reduce(word, n=n, limits=limits, trace=trace, batch_size=batch_size)
                    producer_seconds = time.perf_counter() - start
                    verification_start = time.perf_counter()
                    oracle = collected.CollectedSession()
                    oracle._catalog.configure_rule_cache(max_bytes=0, max_entries=0)
                    expected = oracle.reduce(word, n=n, limits=limits, trace=trace, batch_size=batch_size)
                    normalized = reference_form(certificate, backend)
                    if not _strict_equal(as_arrays(normalized), as_arrays(expected)):
                        raise AssertionError(f"exact replay mismatch: {backend}/{case_name}/{round_index}")
                    mathematical_sha = digest(normalized)
                    verification_seconds = time.perf_counter() - verification_start
                    blob = encoded(certificate)
                    telemetry = {
                        "parent_before": before, "parent_after": rss_snapshot(),
                        "parent_process_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        "certificate_bytes": len(blob), "cache_clear_seconds": cache_clear_seconds,
                    }
                    if backend == "native":
                        telemetry["native_process"] = rss_snapshot(session.pid)
                        telemetry["native_runtime"] = session.last_runtime_info
                    if hasattr(session, "cache_info"):
                        telemetry["cache_info"] = session.cache_info()
                    if source_identity() != sources:
                        raise RuntimeError("source changed during measurement; result not stored")
                    connection.execute("""
                        INSERT INTO measurements(campaign_id,round,position,backend,case_name,n,
                          word_json,producer_seconds,verification_seconds,verified_total_seconds,
                          status,certificate_sha256,mathematical_sha256,certificate_zlib,
                          metrics_json,telemetry_json,verified)
                        VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,1)
                    """, (campaign, round_index, position, backend, case_name, n, encoded(word).decode(),
                          producer_seconds, verification_seconds, producer_seconds + verification_seconds,
                          certificate["status"], hashlib.sha256(blob).hexdigest(), mathematical_sha,
                          zlib.compress(blob), encoded(certificate["metrics"]).decode(), encoded(telemetry).decode()))
                    connection.commit()
                    print(json.dumps({"campaign": campaign, "round": round_index, "backend": backend,
                                      "case": case_name, "producer_seconds": producer_seconds,
                                      "verification_seconds": verification_seconds, "status": certificate["status"],
                                      "steps": certificate["metrics"]["steps"]}), flush=True)
        connection.execute("UPDATE campaigns SET status='complete' WHERE id=?", (campaign,))
        connection.commit()
        return campaign
    except BaseException as error:
        if campaign is not None:
            status = "interrupted" if isinstance(error, KeyboardInterrupt) else "failed"
            connection.execute("UPDATE campaigns SET status=? WHERE id=?", (status, campaign))
            connection.commit()
        raise
    finally:
        for session in sessions.values():
            if hasattr(session, "close"):
                session.close()
        connection.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--cases", nargs="+", choices=CASES, default=["b3_alternating_l10", "b5_positive_l6", "b8_positive_l7"])
    parser.add_argument("--backends", nargs="+", choices=BACKENDS, default=list(BACKENDS))
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--max-steps", type=int, default=20000)
    parser.add_argument("--max-terms", type=int, default=10000)
    parser.add_argument("--expansion-budget", type=int, default=1000000)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--trace", action="store_true")
    parser.add_argument("--cache-policy", choices=["cold_each_input", "retained"], default="cold_each_input")
    args = parser.parse_args()
    run(args.database, cases=args.cases, rounds=args.rounds, backends=args.backends,
        limits=Limits(max_steps=args.max_steps, max_terms=args.max_terms, expansion_budget=args.expansion_budget),
        batch_size=args.batch_size, trace=args.trace, cache_policy=args.cache_policy)


if __name__ == "__main__":
    main()
