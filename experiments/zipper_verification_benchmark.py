"""Compare cold and retained independent Python checks on actual native results.

Each explicitly supplied literal input is produced once in a retained C++
session. Both verifiers check the complete certificate, with alternating order.
SQLite retains exact compressed certificates, source/build identities, cache
counters and separate producer/checker timings. Compilation, setup and SQLite
writes are excluded from the producer-plus-checker totals.
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

from zipper_algebra.zipping_algorithm import native
from zipper_algebra.zipping_algorithm.independent_replay import IndependentVerifier, source_identity
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession, _encode_coefficients


def encoded(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def validate_cases(cases: Any) -> list[dict[str, Any]]:
    if not isinstance(cases, list) or not cases:
        raise ValueError("A nonempty list of literal rank/word cases is required")
    for case in cases:
        if (not isinstance(case, dict) or type(case.get("n")) is not int or case["n"] < 2
                or not isinstance(case.get("word"), list)
                or any(type(x) is not int or not 1 <= abs(x) < case["n"] for x in case["word"])):
            raise ValueError("Each case requires integer n >= 2 and a literal signed word")
    if len({(case["n"], tuple(case["word"])) for case in cases}) != len(cases):
        raise ValueError("Duplicate cases are not a distinct-input benchmark")
    return cases


def benchmark_sources() -> dict[str, Any]:
    return {"reference": source_identity(), "native": native.source_identity(),
            "wrapper": hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
            "benchmark": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def run(database: Path, cases: list[dict[str, Any]], *, limits: Limits,
        verifier_cache_bytes: int = 128 * 1024**2, cache_entries: int = 262144) -> int:
    cases = validate_cases(cases)
    sources = benchmark_sources()
    config = {"cases": cases, "limits": asdict(limits), "trace": False, "batch_size": 64,
              "producer_cache_bytes": 256 * 1024**2, "verifier_cache_bytes": verifier_cache_bytes,
              "cache_entries": cache_entries, "checker_order": "alternate by case index",
              "cold": "Python collected replay with both cache layers disabled",
              "retained": "independent Python unit replay retaining full contexts and labels",
              "certificate_codec": "zlib-json-decimal-coefficients-v1"}
    database.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(database) as db, NativeSession(cache_entries=cache_entries) as producer:
        db.execute("PRAGMA foreign_keys=ON")
        db.executescript("""
            CREATE TABLE IF NOT EXISTS verification_campaigns(
              id INTEGER PRIMARY KEY, started_utc TEXT NOT NULL, configuration_json TEXT NOT NULL,
              sources_json TEXT NOT NULL, native_build_json TEXT NOT NULL, status TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS verification_measurements(
              campaign INTEGER NOT NULL REFERENCES verification_campaigns(id), ordinal INTEGER NOT NULL,
              input_json TEXT NOT NULL, certificate BLOB NOT NULL, certificate_sha256 TEXT NOT NULL,
              producer_seconds REAL NOT NULL, cold_seconds REAL NOT NULL, retained_seconds REAL NOT NULL,
              cache_json TEXT NOT NULL, python_peak_rss_kib INTEGER NOT NULL,
              PRIMARY KEY(campaign,ordinal));
        """)
        cursor = db.execute("INSERT INTO verification_campaigns VALUES(NULL,?,?,?,?,?)", (
            datetime.now(timezone.utc).isoformat(), encoded(config).decode(), encoded(sources).decode(),
            encoded(producer.metadata).decode(), "running"))
        assert cursor.lastrowid is not None
        campaign = cursor.lastrowid
        db.commit()
        verifier = IndependentVerifier(cache_bytes=verifier_cache_bytes, cache_entries=cache_entries)
        try:
            for index, case in enumerate(cases):
                if benchmark_sources() != sources:
                    raise RuntimeError("Reference sources changed during benchmark")
                started = time.perf_counter()
                certificate = producer.reduce(case["word"], n=case["n"], limits=limits)
                generation = time.perf_counter() - started
                timings = {}
                for mode in (("retained", "cold") if index % 2 else ("cold", "retained")):
                    started = time.perf_counter()
                    if mode == "cold":
                        producer.check(certificate, cold=True)
                    else:
                        verifier.check_native(producer, certificate)
                    timings[mode] = time.perf_counter() - started
                payload = encoded(_encode_coefficients(certificate))
                cache = verifier.cache_info()
                db.execute("INSERT INTO verification_measurements VALUES(?,?,?,?,?,?,?,?,?,?)", (
                    campaign, index, encoded(case).decode(), zlib.compress(payload),
                    hashlib.sha256(payload).hexdigest(), generation, timings["cold"], timings["retained"],
                    encoded(cache).decode(), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
                db.commit()
                print(json.dumps({"campaign": campaign, "ordinal": index, "n": case["n"],
                                  "status": certificate["status"], "producer": generation, **timings}), flush=True)
        except BaseException:
            db.execute("UPDATE verification_campaigns SET status='interrupted' WHERE id=?", (campaign,))
            db.commit()
            raise
        db.execute("UPDATE verification_campaigns SET status='complete' WHERE id=?", (campaign,))
        db.commit()
        return campaign


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("cases", type=Path)
    parser.add_argument("--max-steps", type=int, default=4096)
    args = parser.parse_args()
    run(args.database, json.loads(args.cases.read_text()),
        limits=Limits(max_steps=args.max_steps, max_terms=100000, expansion_budget=1000000))


if __name__ == "__main__":
    main()
