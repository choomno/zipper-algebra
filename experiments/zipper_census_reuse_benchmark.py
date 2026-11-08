"""Paired exact standalone, completed-subproblem, and SQLite census timings.

The source census is read-only. Captured inputs are frozen before measurement.
Each backend executes in a separate process. Timed work includes authentication,
process transfers, production, independent checking, exact reconstruction,
canonical output encoding, compression and committed result storage. Session
startup and diagnostic telemetry writes are recorded separately. Only complete,
independently checked outputs are compared; a limit or timeout is never counted
as an equal completed computation. Retained sessions are shared across distinct
inputs and passes; ``fresh`` creates a new native session pair for every input
but retains the SQLite sidecar within that run. It measures process restart
with persisted results, not an empty persistent cache for every input.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
import os
from pathlib import Path
import random
import signal
import sqlite3
import subprocess
import sys
from threading import Event, Thread
import time
from typing import Any, Iterator
import zlib

from zipper_algebra.zipping_algorithm.catalog import Catalog
from zipper_algebra.zipping_algorithm import output_normal_form
from zipper_algebra.zipping_algorithm.census_reuse import ReuseCensus, verifier_binding
from zipper_algebra.zipping_algorithm.input_normal_form import input_normal_form, leaf_requests
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession, _encode_coefficients
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier, VerifierIncomplete
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import execution_identity
from zipper_algebra.zipping_algorithm.output_normal_form import canonical_bytes, canonical_output

BACKENDS = ("standalone", "dag", "census")


class BenchmarkDeadline(BaseException):
    """An interrupt that transport wrappers must not relabel as an IO failure."""


def interface_sources() -> dict[str, str]:
    directory = Path(output_normal_form.__file__).parent
    return {path.name: digest(path.read_bytes()) for path in sorted(directory.glob("*.py"))}


def encoded(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


_IMPORTED_INTERFACE_SOURCES = interface_sources()


def assert_interface_sources() -> None:
    if interface_sources() != _IMPORTED_INTERFACE_SOURCES:
        raise RuntimeError("benchmark Python interface source changed after import")


def schema(db: sqlite3.Connection) -> None:
    db.executescript("""
    CREATE TABLE IF NOT EXISTS reuse_bench_metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS reuse_bench_cases(
      id INTEGER PRIMARY KEY,label TEXT UNIQUE NOT NULL,n INTEGER NOT NULL,
      word TEXT NOT NULL,origin TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS reuse_bench_runs(
      id INTEGER PRIMARY KEY,backend TEXT NOT NULL,cache_mode TEXT NOT NULL,
      started_utc TEXT NOT NULL,status TEXT NOT NULL,config TEXT NOT NULL,
      identity TEXT,setup_seconds REAL,error TEXT);
    CREATE TABLE IF NOT EXISTS reuse_bench_payloads(
      run INTEGER,pass INTEGER,case_id INTEGER,result BLOB,output BLOB,output_sha TEXT,
      PRIMARY KEY(run,pass,case_id));
    CREATE TABLE IF NOT EXISTS reuse_bench_measurements(
      run INTEGER,pass INTEGER,case_id INTEGER,status TEXT NOT NULL,
      seconds REAL NOT NULL,producer_seconds REAL,checker_seconds REAL,
      reconstruction_seconds REAL,store_seconds REAL,peak_tree_rss_bytes INTEGER,
      telemetry TEXT NOT NULL,error TEXT,PRIMARY KEY(run,pass,case_id));
    """)


def frozen_cases(db: sqlite3.Connection, *, establish: bool = False) -> list[dict[str, Any]]:
    cases = [{"label": label, "n": n, "word": json.loads(word), "origin": json.loads(origin)}
             for label, n, word, origin in db.execute(
                 "SELECT label,n,word,origin FROM reuse_bench_cases ORDER BY id")]
    actual = digest(encoded(cases))
    row = db.execute("SELECT value FROM reuse_bench_metadata WHERE key='cases_sha256'").fetchone()
    if row is None and establish and not db.execute("SELECT 1 FROM reuse_bench_runs LIMIT 1").fetchone():
        db.execute("INSERT INTO reuse_bench_metadata VALUES('cases_sha256',?)", (encoded(actual).decode(),))
    elif row is None or json.loads(row[0]) != actual:
        raise ValueError("frozen benchmark input hash differs or is absent")
    return cases


def selected_cases(db: sqlite3.Connection, labels: list[str] | None) -> list[tuple[Any, ...]]:
    rows = list(db.execute("SELECT id,label,n,word FROM reuse_bench_cases ORDER BY id"))
    if labels and set(labels) - {row[1] for row in rows}:
        raise ValueError("unknown benchmark labels: " + repr(sorted(set(labels) - {row[1] for row in rows})))
    selected = [row for row in rows if not labels or row[1] in labels]
    if not selected:
        raise ValueError("no selected benchmark inputs")
    return selected


def synthetic_cases() -> list[dict[str, Any]]:
    """Deterministic full-support inputs, plus explicitly tagged reuse families."""
    result = []
    for n in range(3, 9):
        for length in (8, 9, 10):
            rng = random.Random(1107000 + 100 * n + length)
            word = list(range(1, n))
            word.extend(rng.randrange(1, n) for _ in range(length - len(word)))
            rng.shuffle(word)
            if length == 9:
                word = [letter * rng.choice((-1, 1)) for letter in word]
            elif length == 10:
                word = [-letter for letter in word]
            result.append({"label": f"full-support-b{n}-l{length}", "n": n,
                           "word": word, "origin": {"kind": "constructed-full-support",
                                                     "seed": 1107000 + 100 * n + length}})
    base = [1, -2, 1, -2, 1, -2, 1, -2]
    for name, n, word in (
        ("base", 3, base),
        ("ambient-inclusion", 8, base),
        ("index-translation", 8, [(1 if a > 0 else -1) * (abs(a) + 4) for a in base]),
        ("disjoint-right-twist", 6, base + [5, 5]),
        ("disjoint-interleaved", 6, base[:4] + [5] + base[4:] + [5]),
    ):
        result.append({"label": "reuse-family-" + name, "n": n, "word": word,
                       "origin": {"kind": "constructed-related-context", "family": "b3-alternating"}})
    return result


def input_classification(word: list[int], n: int) -> dict[str, Any]:
    requests = leaf_requests(input_normal_form(word, n=n))
    return {"full_support": {abs(a) for a in word} == set(range(1, n)),
            "unchanged_primitive": len(requests) == 1 and requests[0]["n"] == n
            and requests[0]["word"] == word,
            "unique_primitive_requests": requests}


def capture(source: Path, target: Path, *, recent: int = 256) -> int:
    """Freeze one hard recent input per rank, a known hard word, and longer inputs."""
    if source.resolve() == target.resolve():
        raise ValueError("benchmark and source census must differ")
    target.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(target) as db:
        schema(db)
        previous = dict(db.execute("SELECT key,value FROM reuse_bench_metadata"))
        if previous:
            if json.loads(previous.get("source", "null")) != str(source.resolve()):
                raise ValueError("frozen cases belong to another census")
            frozen_cases(db)
            return db.execute("SELECT COUNT(*) FROM reuse_bench_cases").fetchone()[0]
        with sqlite3.connect(source.resolve().as_uri() + "?mode=ro", uri=True) as src:
            src.row_factory = sqlite3.Row
            src.execute("BEGIN")
            metadata = dict(src.execute("SELECT key,value FROM metadata"))
            config = json.loads(metadata["run_state"])["config_sha"]
            cases, missing = [], []
            for n in range(3, 9):
                rows = src.execute("""SELECT n,length,ordinal,word,metrics,config_sha,certificate_sha
                  FROM words WHERE n=? AND config_sha=? AND status='complete'
                  ORDER BY length DESC,ordinal DESC LIMIT ?""", (n, config, recent)).fetchall()
                if rows:
                    row = max(rows, key=lambda r: json.loads(r["metrics"])["steps"])
                    cases.append({"label": f"census-recent-hard-b{n}", "n": n,
                                  "word": json.loads(row["word"]),
                                  "origin": {"kind": "stored-census", **dict(row)}})
                else:
                    missing.append(f"recent-complete-b{n}")
                stopped = src.execute("""SELECT n,length,ordinal,word,metrics,config_sha,certificate_sha,status
                  FROM words WHERE n=? AND config_sha=? AND status!='complete'
                  ORDER BY length DESC,ordinal DESC LIMIT 1""", (n, config)).fetchone()
                if stopped:
                    cases.append({"label": f"census-stopped-b{n}", "n": n,
                                  "word": json.loads(stopped["word"]),
                                  "origin": {"kind": "stored-census-noncomplete", **dict(stopped)}})
                else:
                    missing.append(f"noncomplete-b{n}")
            row = src.execute("""SELECT n,length,ordinal,word,metrics,config_sha,certificate_sha
                FROM words WHERE n=5 AND length=6 AND ordinal=6519 AND status='complete'
                ORDER BY config_sha,certificate_sha LIMIT 1""").fetchone()
            if row:
                cases.append({"label": "census-known-hard-b5-l6-6519", "n": 5,
                              "word": json.loads(row["word"]),
                              "origin": {"kind": "stored-census", **dict(row)}})
            else:
                missing.append("known-hard-b5-l6-6519")
            cases.extend(synthetic_cases())
            for case in cases:
                db.execute("INSERT INTO reuse_bench_cases VALUES(NULL,?,?,?,?)",
                           (case["label"], case["n"], encoded(case["word"]).decode(),
                            encoded(case["origin"]).decode()))
            for key, value in {"source": str(source.resolve()), "source_metadata": metadata,
                               "captured_utc": datetime.now(timezone.utc).isoformat(),
                               "selection_recent": recent, "method": __doc__,
                               "missing_selections": missing,
                               "cases_sha256": digest(encoded(cases))}.items():
                db.execute("INSERT INTO reuse_bench_metadata VALUES(?,?)", (key, encoded(value).decode()))
        return len(cases)


def tree_rss(pid: int) -> int:
    total, queue = 0, [pid]
    while queue:
        current = queue.pop()
        try:
            total += int(Path(f"/proc/{current}/statm").read_text().split()[1]) * os.sysconf("SC_PAGE_SIZE")
            queue.extend(map(int, Path(f"/proc/{current}/task/{current}/children").read_text().split()))
        except (FileNotFoundError, ProcessLookupError):
            pass
    return total


class PeakMemory:
    """Sample summed Linux process-tree RSS; shared pages can be counted twice."""
    def __init__(self) -> None:
        self.peak = 0
        self.stop = Event()
        self.thread = Thread(target=self._sample, daemon=True)

    def _sample(self) -> None:
        while not self.stop.is_set():
            self.peak = max(self.peak, tree_rss(os.getpid()))
            self.stop.wait(0.05)

    def __enter__(self) -> PeakMemory:
        self.thread.start()
        return self

    def __exit__(self, *args: Any) -> None:
        self.stop.set()
        self.thread.join()
        self.peak = max(self.peak, tree_rss(os.getpid()))


@contextmanager
def deadline(seconds: float) -> Iterator[None]:
    def expired(signum: int, frame: Any) -> None:
        raise BenchmarkDeadline("per-input wall deadline exceeded")
    previous = signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


class Backend:
    def __init__(self, name: str, sidecar: Path, config: dict[str, Any]) -> None:
        if name not in BACKENDS:
            raise ValueError("unknown benchmark backend")
        assert_interface_sources()
        self.name, self.census = name, None
        self.limits = Limits(**config["limits"])
        self.producer = NativeSession(cache_bytes=config["cache_bytes"], cache_entries=config["cache_entries"])
        try:
            self.verifier = CppVerifier(producer_metadata=self.producer.metadata,
                                       cache_bytes=config["cache_bytes"], cache_entries=config["cache_entries"],
                                       allocator_library=Path(config["allocator"]) if config["allocator"] else None)
            if name == "census":
                identity_result = self.producer.reduce_reused([], n=1, limits=self.limits)
                self.census = ReuseCensus(
                    sidecar, produce=lambda n, word: self.producer.reduce_reused(word, n=n, limits=self.limits),
                    verify=self.verifier.check_reused, producer_identity=execution_identity(identity_result),
                    verifier_identity=verifier_binding(self.verifier.metadata), limits=asdict(self.limits),
                    cache_entries=config["cache_entries"])
            self.identity = {"producer": self.producer.metadata, "verifier": self.verifier.metadata,
                             "catalog_sha256": Catalog.bundled().sha256,
                             "comparison_sources": dict(_IMPORTED_INTERFACE_SOURCES),
                             "harness_sha256": digest(Path(__file__).read_bytes())}
        except BaseException:
            self.close()
            raise

    def close(self) -> None:
        if self.census is not None:
            self.census.close()
        if hasattr(self, "verifier"):
            self.verifier.close()
        self.producer.close()

    def evaluate(self, word: list[int], n: int) -> tuple[dict[str, Any], dict[str, Any] | None, dict[str, Any]]:
        if self.census is not None:
            result = self.census.evaluate(word, n=n)
            complete = result["status"] == "complete"
            normal = result.get("output") if complete else None
            if complete and (result.get("provenance") not in (
                    "checked-reconstruction-from-independently-verified-leaves",
                    "independently-verified-leaf-with-checked-identity-embedding")
                             or not isinstance(normal, dict) or normal.get("completion") != "complete"
                             or normal.get("n") != n):
                raise ValueError("census completion lacks independently checked output provenance")
            stats = dict(result.get("stats", {}))
            stats["independently_checked"] = complete
            stats["verification_policy"] = "independent leaf replay once per process plus checked reconstruction"
            return result, normal, stats
        started = time.perf_counter()
        result = (self.producer.reduce(word, n=n, limits=self.limits) if self.name == "standalone"
                  else self.producer.reduce_reused(word, n=n, limits=self.limits))
        producer_seconds = time.perf_counter() - started
        checked, checker_seconds = False, 0.0
        if result["status"] == "complete" or self.name == "standalone":
            started = time.perf_counter()
            try:
                checked = (self.verifier.check(result) if self.name == "standalone"
                           else self.verifier.check_reused(result))
            except VerifierIncomplete as failure:
                checker_seconds = time.perf_counter() - started
                return ({"status": "verification_incomplete", "producer_result": result,
                         "replay_status": failure.replay_status,
                         "replay_diagnostic": failure.replay_diagnostic}, None,
                        {"producer_seconds": producer_seconds, "checker_seconds": checker_seconds,
                         "reconstruction_seconds": 0.0, "store_seconds": 0.0,
                         "independently_checked": False})
            checker_seconds = time.perf_counter() - started
            if checked is not True:
                raise ValueError("independent verifier did not accept result")
        normal = None
        if result["status"] == "complete":
            normal = canonical_output(result["output"], n=n, complete=True)
        return result, normal, {"producer_seconds": producer_seconds,
                                "checker_seconds": checker_seconds,
                                "reconstruction_seconds": 0.0, "store_seconds": 0.0,
                                "independently_checked": checked}

    def telemetry(self) -> dict[str, Any]:
        result = {"producer": self.producer.cache_info(), "verifier": self.verifier.cache_info()}
        if self.census is not None:
            result["census_counts"] = self.census.counts()
        return result


def create_run(database: Path, backend: str, cache_mode: str, config: dict[str, Any]) -> int:
    if backend not in BACKENDS or cache_mode not in ("fresh", "retained"):
        raise ValueError("invalid backend or cache mode")
    with sqlite3.connect(database) as db:
        schema(db)
        frozen_cases(db, establish=True)
        cases = selected_cases(db, config.get("labels"))
        config = dict(config, harness_sha256=digest(Path(__file__).read_bytes()),
                      selected_cases_sha256=digest(encoded(cases)))
        cursor = db.execute("""INSERT INTO reuse_bench_runs
          (backend,cache_mode,started_utc,status,config) VALUES(?,?,?,'pending',?)""",
          (backend, cache_mode, datetime.now(timezone.utc).isoformat(), encoded(config).decode()))
        assert cursor.lastrowid is not None
        return cursor.lastrowid


def worker(database: Path, run: int) -> None:
    backend: Backend | None = None
    initial_hash = digest(Path(__file__).read_bytes())
    with sqlite3.connect(database, timeout=60) as db:
        row = db.execute("SELECT backend,cache_mode,config FROM reuse_bench_runs WHERE id=?", (run,)).fetchone()
        if row is None:
            raise ValueError("unknown benchmark run")
        name, cache_mode, config_json = row
        config = json.loads(config_json)
        frozen_cases(db)
        cases = selected_cases(db, config.get("labels"))
        sidecar = database.with_name(database.stem + f"-run{run}-census.sqlite3")
        setup_seconds = 0.0
        session_generation = 0
        try:
            if config.get("harness_sha256") != initial_hash:
                raise RuntimeError("benchmark source changed before measurement")
            if config.get("selected_cases_sha256") != digest(encoded(cases)):
                raise ValueError("selected benchmark inputs changed")
            if name == "census" and any(Path(str(sidecar) + suffix).exists() for suffix in ("", "-wal", "-shm")):
                raise ValueError("benchmark requires a new sidecar for each run")
            db.execute("UPDATE reuse_bench_runs SET status='running' WHERE id=?", (run,))
            db.commit()
            for repeat in range(config["passes"]):
                for case_id, label, n, word_json in cases:
                    assert_interface_sources()
                    classification = input_classification(json.loads(word_json), n)
                    if backend is None:
                        start = time.perf_counter()
                        backend = Backend(name, sidecar, config)
                        session_generation += 1
                        setup_seconds += time.perf_counter() - start
                        db.execute("UPDATE reuse_bench_runs SET identity=COALESCE(identity,?),setup_seconds=? WHERE id=?",
                                   (encoded(backend.identity).decode(), setup_seconds, run))
                        db.commit()
                    result: dict[str, Any] = {}
                    stats: dict[str, Any] = {}
                    output = None
                    error = None
                    timed_out = False
                    fatal: BaseException | None = None
                    telemetry = {"harness_sha256": initial_hash, "identity": backend.identity,
                                 "session_generation": session_generation,
                                 "input_classification": classification,
                                 "persistent_sidecar_retained": name == "census"}
                    db.execute("INSERT INTO reuse_bench_measurements VALUES(?,?,?,'running',0,NULL,NULL,NULL,NULL,NULL,?,NULL)",
                               (run, repeat, case_id, encoded(telemetry).decode()))
                    db.commit()
                    with PeakMemory() as memory:
                        started = time.perf_counter()
                        try:
                            with deadline(config["case_seconds"]):
                                result, normal, stats = backend.evaluate(json.loads(word_json), n)
                                stats = dict(stats)
                                if result["status"] == "complete" and normal is None:
                                    raise ValueError("completed result omitted its canonical output")
                                if normal is not None:
                                    output = canonical_bytes(normal)
                                status = result["status"]
                                storage_started = time.perf_counter()
                                raw = encoded(_encode_coefficients(result))
                                db.execute("INSERT INTO reuse_bench_payloads VALUES(?,?,?,?,?,?)",
                                           (run, repeat, case_id, zlib.compress(raw),
                                            zlib.compress(output) if output is not None else None,
                                            digest(output) if output is not None else None))
                                db.commit()
                                stats["store_seconds"] = stats.get("store_seconds", 0.0) + time.perf_counter() - storage_started
                        except BenchmarkDeadline as failure:
                            db.rollback()
                            status, error, timed_out = "timeout", repr(failure), True
                        except BaseException as failure:
                            db.rollback()
                            status, error, fatal = "failed", repr(failure), failure
                        elapsed = time.perf_counter() - started
                    telemetry["stats"] = stats
                    if not timed_out and fatal is None:
                        telemetry["cache"] = backend.telemetry()
                    if digest(Path(__file__).read_bytes()) != initial_hash:
                        status = "failed"
                        fatal = RuntimeError("benchmark source changed during measurement")
                        error = repr(fatal)
                    try:
                        assert_interface_sources()
                    except RuntimeError as failure:
                        status, error, fatal = "failed", repr(failure), failure
                    telemetry["sidecar_bytes"] = sum(Path(str(sidecar) + suffix).stat().st_size
                        for suffix in ("", "-wal", "-shm") if Path(str(sidecar) + suffix).exists())
                    db.execute("""UPDATE reuse_bench_measurements SET status=?,seconds=?,producer_seconds=?,
                               checker_seconds=?,reconstruction_seconds=?,store_seconds=?,peak_tree_rss_bytes=?,
                               telemetry=?,error=? WHERE run=? AND pass=? AND case_id=?""",
                               (status, elapsed, stats.get("producer_seconds"),
                                stats.get("checker_seconds"), stats.get("reconstruction_seconds"),
                                stats.get("store_seconds"), memory.peak, encoded(telemetry).decode(), error,
                                run, repeat, case_id))
                    db.commit()
                    if fatal is not None:
                        raise fatal
                    print(encoded({"run": run, "pass": repeat, "case": label, "status": status,
                                   "seconds": elapsed}).decode(), flush=True)
                    if cache_mode == "fresh" or timed_out:
                        backend.close()
                        backend = None
            db.execute("UPDATE reuse_bench_runs SET status='complete' WHERE id=?", (run,))
            db.commit()
        except BaseException as failure:
            db.execute("UPDATE reuse_bench_runs SET status='failed',error=? WHERE id=?", (repr(failure), run))
            db.commit()
            raise
        finally:
            if backend is not None:
                backend.close()


def compare(database: Path) -> dict[str, Any]:
    """Compare backend pairs and time only their commonly completed inputs.

    Same-backend repetitions are reported as cache/reproducibility rechecks,
    never added to cross-backend comparisons. Historical rows without a session
    identity use their retained run identity and are explicitly counted.
    """
    comparisons = rechecks = legacy_identities = legacy_unbound = 0
    complete: dict[tuple[int, int], dict[int, dict[str, Any]]] = {}
    all_rows: dict[tuple[int, int], dict[int, str]] = {}
    summaries, pairs, paired_inputs, standalone_inputs = [], [], set(), set()
    with sqlite3.connect(database) as db:
        frozen_cases(db)
        runs = {row[0]: {"backend": row[1], "cache_mode": row[2], "status": row[3],
                         "setup_seconds": row[4], "config": json.loads(row[5])}
                for row in db.execute("SELECT id,backend,cache_mode,status,setup_seconds,config FROM reuse_bench_runs")}
        for run, repeat, case, status in db.execute("SELECT run,pass,case_id,status FROM reuse_bench_measurements"):
            all_rows.setdefault((run, repeat), {})[case] = status
        rows = db.execute("""SELECT p.run,p.pass,p.case_id,p.output,p.output_sha,r.config,r.identity,
                             m.seconds,m.telemetry
            FROM reuse_bench_payloads p JOIN reuse_bench_measurements m
            ON (p.run,p.pass,p.case_id)=(m.run,m.pass,m.case_id)
            JOIN reuse_bench_runs r ON r.id=p.run
            WHERE m.status='complete' ORDER BY p.run,p.pass,p.case_id""").fetchall()
        same_backend: dict[tuple[Any, ...], bytes] = {}
        for run, repeat, case, blob, sha, config_json, identity_json, seconds, telemetry_json in rows:
            if blob is None:
                raise ValueError("completed measurement has no output payload")
            raw = zlib.decompress(blob)
            if digest(raw) != sha:
                raise ValueError("stored output hash differs")
            if canonical_bytes(json.loads(raw)) != raw or json.loads(raw)["completion"] != "complete":
                raise ValueError("stored completed output is not canonical and complete")
            config, telemetry = json.loads(config_json), json.loads(telemetry_json)
            identity = telemetry.get("identity")
            if identity is None:
                legacy_identities += 1
                identity = json.loads(identity_json)
            if not identity.get("comparison_sources") or not identity.get("harness_sha256"):
                legacy_unbound += 1
                continue
            binding = (encoded(config["limits"]),
                       encoded(identity["producer"]["source_sha256"]),
                       identity["producer"]["wrapper_source_sha256"],
                       encoded(identity["verifier"]["source_sha256"]),
                       encoded(identity["verifier"].get("reference_source_sha256")),
                       encoded(identity["verifier"].get("allocator_library")),
                       identity.get("catalog_sha256"), encoded(identity.get("comparison_sources")),
                       identity.get("harness_sha256"))
            same_key = (case, runs[run]["backend"], *binding)
            if same_key in same_backend:
                if same_backend[same_key] != raw:
                    raise AssertionError(f"exact completed outputs differ at run {run}, pass {repeat}, case {case}")
                rechecks += 1
            else:
                same_backend[same_key] = raw
            complete.setdefault((run, repeat), {})[case] = {"bytes": raw, "binding": binding, "seconds": seconds}
        for left, right in combinations(sorted(all_rows), 2):
            lr, rr = runs[left[0]], runs[right[0]]
            if (lr["backend"] == rr["backend"] or lr["cache_mode"] != rr["cache_mode"]
                    or left[1] != right[1]
                    or any(lr["config"].get(k) != rr["config"].get(k)
                           for k in ("cache_bytes", "cache_entries", "case_seconds", "allocator"))):
                continue
            common = sorted(set(complete.get(left, {})) & set(complete.get(right, {})))
            common = [c for c in common if complete[left][c]["binding"] == complete[right][c]["binding"]]
            if not common:
                continue
            for case in common:
                if complete[left][case]["bytes"] != complete[right][case]["bytes"]:
                    raise AssertionError(f"exact completed outputs differ between {left} and {right}, case {case}")
                comparisons += 1
                paired_inputs.add(case)
                if "standalone" in (lr["backend"], rr["backend"]):
                    standalone_inputs.add(case)
            lt = sum(complete[left][c]["seconds"] for c in common)
            rt = sum(complete[right][c]["seconds"] for c in common)
            union = set(all_rows[left]) | set(all_rows[right])
            pairs.append({"left_run": left[0], "right_run": right[0], "pass": left[1],
                          "left_backend": lr["backend"], "right_backend": rr["backend"],
                          "includes_standalone": "standalone" in (lr["backend"], rr["backend"]),
                          "cache_mode": lr["cache_mode"],
                          "census_sidecar_policy": ("new per run, retained across native session restarts"
                                                     if "census" in (lr["backend"], rr["backend"]) else None),
                          "session_generations": {
                              side: sorted({json.loads(row[0]).get("session_generation", 1)
                                  for row in db.execute("SELECT telemetry FROM reuse_bench_measurements WHERE run=? AND pass=?", key)})
                              for side, key in (("left", left), ("right", right))},
                          "common_complete_case_ids": common, "left_seconds": lt, "right_seconds": rt,
                          "left_over_right_seconds": lt / rt if rt else None,
                          "excluded_cases": [{"case_id": c, "left": all_rows[left].get(c, "not_selected"),
                                              "right": all_rows[right].get(c, "not_selected")}
                                             for c in sorted(union - set(common))]})
        for run, record in sorted(runs.items()):
            measurements = db.execute("""SELECT pass,status,COUNT(*),SUM(seconds),MAX(peak_tree_rss_bytes)
                FROM reuse_bench_measurements WHERE run=? GROUP BY pass,status ORDER BY pass,status""", (run,)).fetchall()
            sessions = sorted({json.loads(row[0]).get("session_generation", 1)
                               for row in db.execute("SELECT telemetry FROM reuse_bench_measurements WHERE run=?", (run,))})
            summaries.append({"run": run, **record, "session_generations": sessions,
                              "persistent_sidecar_retained": record["backend"] == "census",
                              "totals_scope": "unpaired diagnostic totals; use backend_pairs for timing comparisons",
                              "measurements": [{"pass": p, "status": s, "cases": c, "seconds": t,
                                                "peak_tree_rss_bytes": m} for p, s, c, t, m in measurements]})
    completed_inputs = {c for values in complete.values() for c in values}
    return {"exact_comparisons": comparisons, "same_backend_rechecks": rechecks,
            "distinct_inputs_matching_standalone": len(standalone_inputs),
            "standalone_matched_case_ids": sorted(standalone_inputs),
            "distinct_complete_inputs": len(completed_inputs), "paired_complete_inputs": len(paired_inputs),
            "unpaired_complete_case_ids": sorted(completed_inputs - paired_inputs),
            "legacy_measurements_with_run_identity": legacy_identities,
            "legacy_complete_outputs_excluded_from_source_bound_pairs": legacy_unbound,
            "backend_pairs": pairs, "runs": summaries}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("--capture", type=Path)
    parser.add_argument("--worker", type=int)
    parser.add_argument("--run", choices=BACKENDS)
    parser.add_argument("--cache-mode", choices=("fresh", "retained"), default="retained",
                        help="Native session lifetime; both modes retain this run's SQLite sidecar")
    parser.add_argument("--passes", type=int, default=2)
    parser.add_argument("--case-seconds", type=float, default=60)
    parser.add_argument("--cache-mib", type=int, default=256)
    parser.add_argument("--cache-entries", type=int, default=262144)
    parser.add_argument("--max-steps", type=int, default=250000)
    parser.add_argument("--allocator", type=Path)
    parser.add_argument("--labels", nargs="+")
    args = parser.parse_args()
    if args.capture:
        print(capture(args.capture, args.database))
    elif args.worker is not None:
        worker(args.database, args.worker)
    elif args.run:
        if args.passes < 1 or args.case_seconds <= 0 or args.cache_mib < 0 or args.cache_entries < 0:
            parser.error("passes/deadline must be positive and cache budgets nonnegative")
        config = {"passes": args.passes, "case_seconds": args.case_seconds,
                  "cache_bytes": args.cache_mib * 1024**2, "cache_entries": args.cache_entries,
                  "limits": asdict(Limits(max_steps=args.max_steps, max_terms=100000,
                                          expansion_budget=1000000)),
                  "allocator": str(args.allocator.resolve()) if args.allocator else None,
                  "labels": args.labels}
        run = create_run(args.database, args.run, args.cache_mode, config)
        subprocess.run([sys.executable, str(Path(__file__).resolve()), str(args.database),
                        "--worker", str(run)], check=True)
        print(json.dumps(compare(args.database), indent=2))
    else:
        print(json.dumps(compare(args.database), indent=2))


if __name__ == "__main__":
    main()
