"""Retained-session benchmark on every literal word in a finite corpus.

Default input: all 341 signed words in B3 through length 4 and all 259 in B4
through length 3. Inputs are never freely reduced or identified by isotopy.
Each backend retains its own transition cache across distinct inputs. Producer
order rotates per input. Every returned certificate is separately checked by
an uncached Python replay, and all records live in one SQLite database.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from itertools import product
import json
from pathlib import Path
import resource
import time
from typing import Any
import zlib

from experiments import zipper_native_engine_benchmark as shared
from zipper_algebra.zipping_algorithm import collected, memoized
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession, _strict_equal


DEFAULT_BOUNDS = ((3, 4), (4, 3))


def literal_corpus(bounds: tuple[tuple[int, int], ...] = DEFAULT_BOUNDS) -> list[dict[str, Any]]:
    if not bounds or len({n for n, _ in bounds}) != len(bounds):
        raise ValueError("nonempty distinct ranks required")
    result = []
    for n, cutoff in bounds:
        if type(n) is not int or n < 2 or type(cutoff) is not int or cutoff < 0:
            raise ValueError("integer rank at least two and nonnegative cutoff required")
        alphabet = tuple(range(-(n - 1), 0)) + tuple(range(1, n))
        for length in range(cutoff + 1):
            for ordinal, word in enumerate(product(alphabet, repeat=length)):
                result.append({"case_name": f"b{n}_l{length}_ordinal{ordinal}",
                               "n": n, "length": length, "ordinal": ordinal,
                               "word": list(word)})
    return result


def corpus_source_identity(corpus: list[dict[str, Any]]) -> dict[str, Any]:
    return {**shared.source_identity(),
            "corpus_harness": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "literal_input_manifest": corpus,
            "literal_input_manifest_sha256": shared.digest(corpus)}


def run(
    database: Path,
    *,
    bounds: tuple[tuple[int, int], ...] = DEFAULT_BOUNDS,
    backends: tuple[str, ...] = ("native", "python_memoized"),
    limits: Limits | None = None,
    batch_size: int = 64,
) -> int:
    if not backends or len(set(backends)) != len(backends) or any(b not in shared.BACKENDS for b in backends):
        raise ValueError("distinct supported backends required")
    corpus = literal_corpus(bounds)
    limits = limits or Limits(max_steps=100000, max_terms=100000, expansion_budget=1000000)
    sources = corpus_source_identity(corpus)
    config = {
        "type": "retained_literal_corpus", "bounds": bounds, "inputs": len(corpus),
        "literal_input_manifest_sha256": shared.digest(corpus), "backends": backends,
        "limits": asdict(limits), "batch_size": batch_size, "trace": False,
        "cache_policy": "retained across every literal input; fresh session per backend",
        "cache_budgets": {"max_entries": 262144, "max_bytes": 256 * 1024**2,
                          "byte_measure": "native estimated typed storage; Python encoded key/payload"},
        "ordering": "rotate producer order per input, retain identical corpus order per backend",
        "verification": "separate fresh uncached full-state Python replay for every backend and input",
        "timing": "entire reduce call including serialization and IPC; SQL compression/storage excluded",
        "scope": "finite performance corpus, not an authoritative census extension",
    }
    db = shared.initialize_database(database)
    sessions: dict[str, Any] = {}
    setup = {}
    campaign = None
    try:
        for backend in backends:
            start = time.perf_counter()
            if backend == "native":
                session = NativeSession(cache_entries=262144, cache_bytes=256 * 1024**2)
            elif backend == "python_memoized":
                session = memoized.MemoizedSession(cache_entries=262144, cache_bytes=256 * 1024**2)
            else:
                session = collected.CollectedSession()
            sessions[backend] = session
            setup[backend] = {"seconds": time.perf_counter() - start}
            if backend == "native":
                setup[backend]["build"] = session.metadata
        cursor = db.execute(
            "INSERT INTO campaigns(started_utc,configuration_json,sources_json,setup_json,status) VALUES(?,?,?,?,?)",
            (datetime.now(timezone.utc).isoformat(), shared.encoded(config).decode(),
             shared.encoded(sources).decode(), shared.encoded(setup).decode(), "running"),
        )
        campaign = int(cursor.lastrowid)
        db.commit()
        for index, case in enumerate(corpus):
            n, word, name = case["n"], case["word"], case["case_name"]
            at = index % len(backends)
            order = backends[at:] + backends[:at]
            for position, backend in enumerate(order):
                if corpus_source_identity(corpus) != sources:
                    raise RuntimeError("source changed during corpus benchmark")
                session = sessions[backend]
                before = shared.rss_snapshot()
                started = time.perf_counter()
                cert = session.reduce(word, n=n, limits=limits, trace=False, batch_size=batch_size)
                producer = time.perf_counter() - started
                started = time.perf_counter()
                oracle = collected.CollectedSession()
                oracle._catalog.configure_rule_cache(max_bytes=0, max_entries=0)
                expected = oracle.reduce(word, n=n, limits=limits, trace=False, batch_size=batch_size)
                normalized = shared.reference_form(cert, backend)
                if not _strict_equal(shared.as_arrays(normalized), shared.as_arrays(expected)):
                    raise AssertionError(f"exact replay mismatch: {backend}/{name}")
                mathematical_sha = shared.digest(normalized)
                verification = time.perf_counter() - started
                blob = shared.encoded(cert)
                telemetry: dict[str, Any] = {
                    "input_index": index, "parent_before": before,
                    "parent_after": shared.rss_snapshot(), "certificate_bytes": len(blob),
                    "parent_process_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                }
                if backend == "native":
                    telemetry["native_process"] = shared.rss_snapshot(session.pid)
                    telemetry["native_runtime"] = session.last_runtime_info
                if hasattr(session, "cache_info"):
                    telemetry["cache_info"] = session.cache_info()
                if corpus_source_identity(corpus) != sources:
                    raise RuntimeError("source changed during corpus measurement; result not stored")
                db.execute("""
                    INSERT INTO measurements(campaign_id,round,position,backend,case_name,n,
                      word_json,producer_seconds,verification_seconds,verified_total_seconds,
                      status,certificate_sha256,mathematical_sha256,certificate_zlib,
                      metrics_json,telemetry_json,verified)
                    VALUES(?,0,?,?,?,?,?,?,?,?,?,?,?,?,?,?,1)
                """, (campaign, position, backend, name, n, shared.encoded(word).decode(),
                      producer, verification, producer + verification, cert["status"],
                      hashlib.sha256(blob).hexdigest(), mathematical_sha, zlib.compress(blob),
                      shared.encoded(cert["metrics"]).decode(), shared.encoded(telemetry).decode()))
                db.commit()
            if (index + 1) % 25 == 0 or index + 1 == len(corpus):
                print(json.dumps({"campaign": campaign, "verified_inputs": index + 1,
                                  "total_inputs": len(corpus), "backends": backends}), flush=True)
        db.execute("UPDATE campaigns SET status='complete' WHERE id=?", (campaign,))
        db.commit()
        return campaign
    except BaseException as error:
        if campaign is not None:
            db.execute("UPDATE campaigns SET status=? WHERE id=?",
                       ("interrupted" if isinstance(error, KeyboardInterrupt) else "failed", campaign))
            db.commit()
        raise
    finally:
        for session in sessions.values():
            if hasattr(session, "close"):
                session.close()
        db.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--n3-length", type=int, default=4)
    parser.add_argument("--n4-length", type=int, default=3)
    parser.add_argument("--backends", nargs="+", choices=shared.BACKENDS, default=["native", "python_memoized"])
    args = parser.parse_args()
    run(args.database, bounds=((3, args.n3_length), (4, args.n4_length)), backends=tuple(args.backends))


if __name__ == "__main__":
    main()
