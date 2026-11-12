"""Compare independent Python and C++ replay of saved census certificates.

The census is opened read-only. A separate SQLite database freezes the exact
compressed certificates, their original identifiers and hashes, implementation
identities, complete-check timings, cache counters, and sampled process memory.
Each backend runs in a fresh process; cases stream from SQLite. The first pass
starts with empty caches and retains them across distinct inputs. The second
pass repeats the identical sequence with those caches retained. Timings include
source authentication, encoding, transfers and exact comparison, but exclude
setup, decompression, memory sampling, cache telemetry and SQLite writes.
Recapture reuses frozen cases only for the same census path and preserved source
identity/version; later census progress does not change the selected inputs.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import signal
from pathlib import Path
import sqlite3
import subprocess
import sys
import time
from typing import Any
import zlib


def encoded(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def load_certificate(blob: bytes, expected: str, size: int) -> dict[str, Any]:
    from zipper_algebra.zipping_algorithm.native import _decimal_integer
    raw = zlib.decompress(blob)
    if len(raw) != size or sha(raw) != expected:
        raise ValueError("Captured certificate size or SHA256 differs")
    value = json.loads(raw, parse_int=_decimal_integer)
    if not isinstance(value, dict):
        raise ValueError("Certificate must be an object")
    return value


def capture(source: Path, target: Path, recent: int = 512) -> int:
    """Freeze three workload quantiles per rank plus one known hard case."""
    if source.resolve() == target.resolve():
        raise ValueError("The benchmark database must differ from the census")
    target.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(f"file:{source.resolve()}?mode=ro", uri=True) as src, sqlite3.connect(target) as db:
        src.row_factory = sqlite3.Row
        src.execute("BEGIN")  # one coherent read snapshot, no census writes
        db.executescript("""
        CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS cases(
          id INTEGER PRIMARY KEY,n INTEGER,length INTEGER,ordinal INTEGER,
          word TEXT,input_sha TEXT,config_sha TEXT,certificate BLOB,
          certificate_sha TEXT,certificate_bytes INTEGER,metrics TEXT,
          source_checker_ns INTEGER,selection TEXT,UNIQUE(n,length,ordinal));
        CREATE TABLE IF NOT EXISTS runs(
          id INTEGER PRIMARY KEY,backend TEXT,started_utc TEXT,status TEXT,
          setup_seconds REAL,identity TEXT,cache_bytes INTEGER,cache_entries INTEGER,
          benchmark_sha TEXT,peak_tree_rss_bytes INTEGER,peak_process_rss TEXT,error TEXT);
        CREATE TABLE IF NOT EXISTS measurements(
          run INTEGER,pass INTEGER,case_id INTEGER,seconds REAL,cache TEXT,
          PRIMARY KEY(run,pass,case_id));
        """)
        metadata = dict(src.execute("SELECT key,value FROM metadata"))
        if db.execute("SELECT COUNT(*) FROM cases").fetchone()[0]:
            saved = dict(db.execute("SELECT key,value FROM metadata"))
            if json.loads(saved.get("source", "null")) != str(source.resolve()):
                raise ValueError("Frozen benchmark belongs to a different census source")
            previous = json.loads(saved["source_metadata"])
            if previous.get("identity") != metadata.get("identity"):
                raise ValueError("Census identity differs from the frozen benchmark")
            version = previous.get("active_version")
            if version is not None:
                has_versions = src.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='versions'").fetchone()
                if has_versions:
                    row = src.execute("SELECT identity FROM versions WHERE hash=?", (version,)).fetchone()
                    if row is None or sha(encoded(json.loads(row[0])).encode()) != version:
                        raise ValueError("Frozen benchmark source version is absent or differs")
                elif version != metadata.get("active_version"):
                    raise ValueError("Frozen benchmark source version differs")
            return db.execute("SELECT COUNT(*) FROM cases").fetchone()[0]
        config = json.loads(metadata["last_config"])
        config_sha = json.loads(metadata["run_state"])["config_sha"]
        selected = []
        for n in range(3, 9):
            rows = list(src.execute("""SELECT * FROM words WHERE n=? AND config_sha=?
                AND status='complete' ORDER BY length DESC,ordinal DESC LIMIT ?""", (n, config_sha, recent)))
            rows.sort(key=lambda row: (json.loads(row["metrics"])["steps"], row["length"], row["ordinal"]))
            if len(rows) < 3:
                raise ValueError(f"Insufficient current census cases for rank {n}")
            for fraction, label in ((0.5, "recent-median"), (0.9, "recent-p90"), (1.0, "recent-maximum")):
                selected.append((rows[round((len(rows) - 1) * fraction)], label))
        hard = src.execute("SELECT * FROM words WHERE n=5 AND length=6 AND ordinal=6519").fetchone()
        if hard is not None and hard["status"] == "complete":
            selected.append((hard, "known-hard"))
        for row, label in selected:
            certificate = load_certificate(row["certificate"], row["certificate_sha"], row["certificate_bytes"])
            if certificate["schema"] != "zipper-native-run-v1":
                if label == "known-hard":
                    continue
                raise ValueError("Captured case does not use the native producer")
            if certificate["n"] != row["n"] or certificate["input_word"] != json.loads(row["word"]):
                raise ValueError("Census input differs from captured certificate")
            db.execute("INSERT OR IGNORE INTO cases VALUES(NULL,?,?,?,?,?,?,?,?,?,?,?,?)", (
                row["n"], row["length"], row["ordinal"], row["word"], row["input_sha"],
                row["config_sha"], row["certificate"], row["certificate_sha"], row["certificate_bytes"],
                row["metrics"], row["checker_ns"], label))
        for key, value in {"source": str(source.resolve()), "captured_utc": datetime.now(timezone.utc).isoformat(),
                           "source_config": config, "source_metadata": metadata, "selection_recent": recent,
                           "method": __doc__}.items():
            db.execute("INSERT INTO metadata VALUES(?,?)", (key, encoded(value)))
        db.commit()
        return db.execute("SELECT COUNT(*) FROM cases").fetchone()[0]


def worker(database: Path, backend: str, passes: int, budget: int, entries: int,
           allocator_library: Path | None = None) -> int:
    from zipper_algebra.zipping_algorithm.independent_replay import IndependentVerifier, source_identity
    from zipper_algebra.zipping_algorithm.native import NativeSession
    from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier
    with sqlite3.connect(database, timeout=30) as db:
        db.execute("PRAGMA journal_mode=WAL")
        cursor = db.execute("""INSERT INTO runs(backend,started_utc,status,cache_bytes,cache_entries,benchmark_sha)
                               VALUES(?,?,'initializing',?,?,?)""", (
            backend, datetime.now(timezone.utc).isoformat(), budget, entries, sha(Path(__file__).read_bytes())))
        run = cursor.lastrowid
        db.commit()
        print(encoded({"run": run, "backend": backend, "event": "initializing"}), flush=True)
        verifier: CppVerifier | IndependentVerifier | None = None
        try:
            # Producer initialization authenticates metadata only. Its process is
            # closed before either verifier starts, and never computes a result.
            with NativeSession(cache_bytes=0, cache_entries=0) as producer:
                producer_metadata = producer.metadata
            started = time.perf_counter()
            if backend == "cpp":
                verifier = CppVerifier(producer_metadata=producer_metadata, cache_bytes=budget, cache_entries=entries,
                                       allocator_library=allocator_library)
                check = verifier.check
                identity = verifier.metadata
            else:
                python_verifier = IndependentVerifier(cache_bytes=budget, cache_entries=entries)
                verifier = python_verifier
                check = lambda certificate: python_verifier.check_native(producer, certificate)
                identity = {"reference": source_identity(), "producer": producer_metadata}
            setup = time.perf_counter() - started
            db.execute("UPDATE runs SET setup_seconds=?,identity=?,status='running' WHERE id=?",
                       (setup, encoded(identity), run))
            db.commit()
            print(encoded({"run": run, "event": "ready", "setup_seconds": setup}), flush=True)
            case_ids = [row[0] for row in db.execute("SELECT id FROM cases ORDER BY id").fetchall()]
            for repeat in range(passes):
                for case_id in case_ids:
                    case = db.execute("SELECT id,n,length,ordinal,certificate,certificate_sha,certificate_bytes FROM cases WHERE id=?",
                                      (case_id,)).fetchone()
                    case_id, n, length, ordinal, blob, digest, size = case
                    certificate = load_certificate(blob, digest, size)
                    started = time.perf_counter()
                    if check(certificate) is not True:
                        raise ValueError("Verifier failed to accept an actual census certificate")
                    elapsed = time.perf_counter() - started
                    db.execute("INSERT INTO measurements VALUES(?,?,?,?,?)", (
                        run, repeat, case_id, elapsed, encoded(verifier.cache_info())))
                    db.commit()
                    print(encoded({"run": run, "pass": repeat, "case": case_id, "n": n,
                                   "length": length, "ordinal": ordinal, "seconds": elapsed}), flush=True)
            db.execute("UPDATE runs SET status='complete' WHERE id=?", (run,))
            db.commit()
        except BaseException as error:
            db.execute("UPDATE runs SET status='failed',error=? WHERE id=?", (repr(error), run))
            db.commit()
            raise
        finally:
            if isinstance(verifier, CppVerifier):
                verifier.close()
    return 0


def tree_rss(pid: int) -> dict[int, int]:
    """Linux RSS by process, retaining separate components to expose shared pages."""
    result = {}
    queue = [pid]
    while queue:
        current = queue.pop()
        try:
            fields = Path(f"/proc/{current}/statm").read_text().split()
            result[current] = int(fields[1]) * os.sysconf("SC_PAGE_SIZE")
            children = Path(f"/proc/{current}/task/{current}/children").read_text().split()
            queue.extend(map(int, children))
        except (FileNotFoundError, ProcessLookupError):
            pass
    return result


def run_backend(database: Path, backend: str, passes: int, budget: int, entries: int,
                allocator_library: Path | None = None) -> None:
    command = [sys.executable, str(Path(__file__).resolve()), str(database), "--worker", backend,
               "--passes", str(passes), "--cache-bytes", str(budget), "--cache-entries", str(entries)]
    if allocator_library is not None:
        command.extend(["--allocator-library", str(allocator_library.resolve())])
    with sqlite3.connect(database, timeout=30) as db:
        before = db.execute("SELECT COALESCE(MAX(id),0) FROM runs").fetchone()[0]
    peak = 0
    process_peaks: dict[int, int] = {}
    with subprocess.Popen(command, start_new_session=True) as process:
        try:
            while process.poll() is None:
                rss = tree_rss(process.pid)
                peak = max(peak, sum(rss.values()))
                for pid, value in rss.items():
                    process_peaks[pid] = max(process_peaks.get(pid, 0), value)
                time.sleep(0.05)
        except BaseException as error:
            # This session belongs solely to this benchmark worker and its
            # verifier children; census workers are never in this process group.
            try:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass
            finally:
                # A second interrupt during the graceful wait still reaches
                # forced cleanup and cannot strand a verifier descendant.
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                with sqlite3.connect(database, timeout=30) as db:
                    db.execute("UPDATE runs SET status='interrupted',error=COALESCE(error,?),peak_tree_rss_bytes=?,peak_process_rss=? WHERE id>? AND backend=? AND status IN ('initializing','running')",
                               (repr(error), peak, encoded(process_peaks), before, backend))
                    db.commit()
            raise
        finally:
            # Also handle an unexpected worker death and stubborn descendants.
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    with sqlite3.connect(database, timeout=30) as db:
        db.execute("UPDATE runs SET peak_tree_rss_bytes=?,peak_process_rss=? WHERE id>? AND backend=?",
                   (peak, encoded(process_peaks), before, backend))
        if process.returncode:
            db.execute("UPDATE runs SET status='failed',error=COALESCE(error,?) WHERE id>? AND backend=?",
                       (f"Worker exit code {process.returncode}", before, backend))
        db.commit()
    if process.returncode:
        raise RuntimeError(f"{backend} benchmark failed with exit code {process.returncode}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--capture-only", action="store_true")
    parser.add_argument("--recent", type=int, default=512)
    parser.add_argument("--passes", type=int, default=2)
    parser.add_argument("--cache-bytes", type=int, default=128 * 1024**2)
    parser.add_argument("--cache-entries", type=int, default=262144)
    parser.add_argument("--allocator-library", type=Path, help="Optional allocator shared library for the C++ verifier child")
    parser.add_argument("--worker", choices=("python", "cpp"))
    parser.add_argument("--order", choices=("python-first", "cpp-first"), default="python-first")
    args = parser.parse_args()
    if min(args.recent, args.passes) < 1 or min(args.cache_bytes, args.cache_entries) < 0:
        parser.error("Passes and recent count must be positive; cache budgets must be nonnegative")
    if args.worker:
        return worker(args.database, args.worker, args.passes, args.cache_bytes, args.cache_entries, args.allocator_library)
    if args.source:
        print(encoded({"captured_cases": capture(args.source, args.database, args.recent)}), flush=True)
    if args.capture_only:
        return 0
    backends = ("python", "cpp") if args.order == "python-first" else ("cpp", "python")
    for backend in backends:
        run_backend(args.database, backend, args.passes, args.cache_bytes, args.cache_entries, args.allocator_library)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
