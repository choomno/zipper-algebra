"""Capture real zipper collection batches and benchmark exact optional GPU backends.

Capture uses a separate collected session; it never modifies a census database.
The GPU prototype accelerates coefficient collection, not the full zipper run.
"""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from dataclasses import asdict
import hashlib
from itertools import groupby
import json
from pathlib import Path
import platform
import sqlite3
import time
from typing import Any, Iterable
import zlib

from zipper_algebra.zipping_algorithm import collected
from zipper_algebra.zipping_algorithm.interpreter import Limits

Contribution = tuple[int, int, int, int]
DEFAULT_WORKLOADS = [
    {"name": "negative-cycle-B3-L10", "n": 3, "word": [-1, -2] * 5},
    {"name": "positive-cycle-B3-L8", "n": 3, "word": [1, 2] * 4},
    {"name": "alternating-cycle-B3-L8", "n": 3, "word": [1, -2] * 4},
    {"name": "positive-chain-B5-L4", "n": 5, "word": [1, 2, 3, 4]},
    {"name": "negative-chain-B8-L7", "n": 8, "word": [-1, -2, -3, -4, -5, -6, -7]},
]


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def digest(value: Any) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


def encode(value: Any) -> tuple[bytes, str]:
    raw = canonical(value)
    return zlib.compress(raw), hashlib.sha256(raw).hexdigest()


def decode(payload: bytes, expected: str) -> Any:
    raw = zlib.decompress(payload)
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("compressed artifact hash mismatch")
    return json.loads(raw)


def reference_collection(records: Iterable[Contribution]) -> tuple[Contribution, ...]:
    """Independent sorted-group Python-integer sum, preserving each input root."""
    result = []
    for key, items in groupby(sorted(records), key=lambda row: row[:3]):
        coefficient = sum(row[3] for row in items)
        if coefficient:
            result.append((*key, coefficient))
    return tuple(result)


class CaptureSession(collected.CollectedSession):
    """Observe real incoming accumulators through the protected dispatch hook."""

    def __init__(self, root_id: int):
        super().__init__()
        self.root_id = root_id
        self.batches: list[dict[str, Any]] = []
        self.frontier: dict[str, Any] = {}
        self.output: dict[str, Any] = {}
        self.next_state = 0
        self.before_latest: tuple[dict[str, Any], dict[str, Any]] | None = None

    def _execute_batch(self, selected, limits, trace, workers, batch_size, pool):
        if not self.batches:
            self.frontier = {key: deepcopy(item) for key, item in selected}
        results, pool = super()._execute_batch(
            selected, limits, trace, workers, batch_size, pool
        )
        if any("error" in result for result in results):
            return results, pool
        removed = {key for key, _ in selected}
        for key, item in selected:
            if key not in self.frontier or self.frontier[key] != item:
                raise ValueError("shadow frontier disagrees with scheduler")
        self.before_latest = (deepcopy(self.frontier), deepcopy(self.output))
        children = [item for result in results for item in result["children"]]
        added = [item for result in results for item in result["output"]]
        records: list[Contribution] = []
        states = []
        expected: list[Contribution] = []
        batch = len(self.batches)
        for pending, current, incoming, excluded in (
            (True, self.frontier, children, removed),
            (False, self.output, added, set()),
        ):
            affected = {collected.continuation_key(x, pending=pending) for x in incoming}
            sources = [*incoming, *(current[key] for key in sorted(affected)
                                   if key in current and key not in excluded)]
            identities: dict[str, int] = {}
            for item in sources:
                key = collected.continuation_key(item, pending=pending)
                if key not in identities:
                    identities[key] = self.next_state
                    self.next_state += 1
                    states.append({"state_id": identities[key], "batch": batch,
                                   "pending": pending, "continuation": key})
                records.extend((self.root_id, identities[key], exponent, coefficient)
                               for exponent, coefficient in item["coefficient"])
            updates, _, _ = collected._stage(
                current, incoming, pending=pending, removed=excluded
            )
            for key, item in updates.items():
                if item is not None:
                    expected.extend((self.root_id, identities[key], exponent, coefficient)
                                    for exponent, coefficient in item["coefficient"])
            collected._commit(current, updates, excluded)
        expected = sorted(expected)
        if reference_collection(records) != tuple(expected):
            raise ValueError("captured contributions disagree with scheduler collection")
        rules = Counter()
        for result in results:
            event = result["unit_evidence"]
            rule = event.get("rule_id", event.get("catalog_transition", {}).get(
                "request", {}).get("rule_id"))
            rules[str(rule)] += 1
        self.batches.append({"batch": batch, "records": records, "states": states,
                             "expected": expected, "selected": len(selected),
                             "rule_applications": dict(rules)})
        return results, pool


    def reconcile(self, certificate: dict[str, Any]) -> None:
        """Match the committed shadow against the complete saved scheduler state."""
        committed = certificate["metrics"]["batches"]
        if len(self.batches) == committed:
            frontier, output = self.frontier, self.output
        elif len(self.batches) == committed + 1 and self.before_latest is not None:
            # The dispatcher ran, but the enclosing atomic batch hit a later limit.
            frontier, output = self.before_latest
        else:
            raise ValueError("capture/scheduler committed batch counts disagree")
        if canonical(list(frontier.values())) != canonical(certificate["pending"]):
            raise ValueError("final shadow frontier differs from scheduler certificate")
        if canonical(list(output.values())) != canonical(certificate["output"]):
            raise ValueError("final shadow output differs from scheduler certificate")


def initialize_database(db: sqlite3.Connection) -> None:
    db.executescript("""
        CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS workloads(
            root INTEGER PRIMARY KEY,name TEXT NOT NULL,n INTEGER NOT NULL,
            word TEXT NOT NULL,status TEXT NOT NULL,certificate BLOB NOT NULL,
            certificate_sha TEXT NOT NULL,capture_wall REAL NOT NULL,
            replay_wall REAL NOT NULL,metadata TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS frontiers(
            root INTEGER NOT NULL,batch INTEGER NOT NULL,payload BLOB NOT NULL,
            payload_sha TEXT NOT NULL,contributions INTEGER NOT NULL,
            PRIMARY KEY(root,batch),FOREIGN KEY(root) REFERENCES workloads(root));
        CREATE TABLE IF NOT EXISTS benchmarks(
            id INTEGER PRIMARY KEY,created_utc TEXT DEFAULT CURRENT_TIMESTAMP,
            backend TEXT NOT NULL,result TEXT NOT NULL);
    """)
    db.execute("INSERT OR IGNORE INTO metadata VALUES('schema',?)",
               (json.dumps("zipper-gpu-frontiers-v1"),))
    schema = json.loads(db.execute("SELECT value FROM metadata WHERE key='schema'").fetchone()[0])
    if schema != "zipper-gpu-frontiers-v1":
        raise ValueError("unexpected capture database schema")


def capture(database: Path, workloads: list[dict[str, Any]], *,
            max_steps: int = 4096, batch_size: int = 64) -> dict[str, Any]:
    database = database.resolve()
    if database.exists() and database.stat().st_size:
        raise ValueError("capture requires a new database path")
    if type(max_steps) is not int or max_steps < 1:
        raise ValueError("capture requires a positive step budget")
    if not workloads:
        raise ValueError("at least one workload is required")
    database.parent.mkdir(parents=True, exist_ok=True)
    limits = Limits(max_steps=max_steps, max_terms=10000, expansion_budget=1000000)
    with sqlite3.connect(database) as db:
        initialize_database(db)
        if db.execute("SELECT count(*) FROM workloads").fetchone()[0]:
            raise ValueError("capture database already contains workloads; choose a new path")
        intended = [{"root": root, "name": work.get("name", str(root)),
                     "n": work["n"], "word": work["word"]}
                    for root, work in enumerate(workloads)]
        db.execute("INSERT INTO metadata VALUES('intended_workloads',?)", (json.dumps(intended),))
        db.execute("INSERT INTO metadata VALUES('capture_complete','false')")
        db.commit()
        manifest = []
        for root, workload in enumerate(workloads):
            session = CaptureSession(root)
            started = time.perf_counter()
            cert = session.reduce(workload["word"], n=workload["n"], limits=limits,
                                  trace=False, workers=1, batch_size=batch_size)
            capture_wall = time.perf_counter() - started
            session.reconcile(cert)
            cert = json.loads(json.dumps(cert))
            started = time.perf_counter()
            collected.CollectedSession().check(cert)
            replay_wall = time.perf_counter() - started
            batches = session.batches[:cert["metrics"]["batches"]]
            payload, sha = encode(cert)
            metadata = {"limits": asdict(limits), "batch_size": batch_size,
                        "source_sha256": cert["source_sha256"],
                        "catalog_sha256": cert["catalog_sha256"],
                        "capture_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        "certificate_reloaded_and_replayed": True,
                        "final_shadow_reconciled": True,
                        "capture_scope": "actual committed batch incoming terms plus affected prior accumulators",
                        "metrics": cert["metrics"]}
            db.execute("INSERT INTO workloads VALUES(?,?,?,?,?,?,?,?,?,?)", (
                root, workload.get("name", str(root)), workload["n"],
                json.dumps(workload["word"]), cert["status"], payload, sha,
                capture_wall, replay_wall, json.dumps(metadata)))
            contribution_count = 0
            for batch in batches:
                compressed, hashed = encode(batch)
                count = len(batch["records"])
                contribution_count += count
                db.execute("INSERT INTO frontiers VALUES(?,?,?,?,?)", (
                    root, batch["batch"], compressed, hashed, count))
            row = {"root": root, "name": workload.get("name", str(root)),
                   "n": workload["n"], "word": workload["word"], "status": cert["status"],
                   "batches": len(batches), "contributions": contribution_count,
                   "capture_wall": capture_wall, "replay_wall": replay_wall,
                   "certificate_sha": sha}
            manifest.append(row)
            db.commit()
            print(json.dumps(row), flush=True)
        db.execute("UPDATE metadata SET value='true' WHERE key='capture_complete'")
        db.commit()
        return {"database": str(database), "workloads": manifest}


def load_frontiers(database: Path, roots: set[int] | None = None) -> tuple[list[Contribution], dict[str, Any]]:
    records = []
    expected = []
    workload_rows = []
    states_seen = set()
    expected_batches: dict[int, int] = {}
    seen_batches: dict[int, list[int]] = {}
    with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as db:
        db.execute("BEGIN")
        capture_metadata = {key: json.loads(value) for key, value in db.execute("SELECT key,value FROM metadata")}
        if "capture_complete" in capture_metadata and capture_metadata["capture_complete"] is not True:
            raise ValueError("capture is not marked complete")
        available_roots = {row[0] for row in db.execute("SELECT root FROM workloads")}
        intended = capture_metadata.get("intended_workloads")
        if intended is not None:
            intended_rows = {row["root"]: (row["n"], row["word"]) for row in intended}
            if len(intended_rows) != len(intended) or set(intended_rows) != available_roots:
                raise ValueError("workload roots differ from intended capture manifest")
            saved_rows = {root: (n, json.loads(word)) for root, n, word in db.execute(
                "SELECT root,n,word FROM workloads")}
            if intended_rows != saved_rows:
                raise ValueError("workload inputs differ from intended capture manifest")
        if roots is not None and not roots <= available_roots:
            raise ValueError("requested workload root is absent")
        selected_roots = available_roots if roots is None else roots
        for root, n, word, status, cert_payload, cert_sha in db.execute(
            "SELECT root,n,word,status,certificate,certificate_sha FROM workloads ORDER BY root"
        ):
            if roots is not None and root not in roots:
                continue
            cert = decode(cert_payload, cert_sha)
            if cert["input_word"] != json.loads(word) or cert["n"] != n or cert["status"] != status:
                raise ValueError("workload/certificate input binding mismatch")
            count = cert["metrics"]["batches"]
            if type(count) is not int or count < 0:
                raise ValueError("invalid certificate batch count")
            expected_batches[root] = count
            seen_batches[root] = []
            workload_rows.append({"root": root, "n": n, "word": json.loads(word),
                                  "status": status, "certificate_sha": cert_sha})
        for root, batch, payload, sha, count in db.execute(
            "SELECT root,batch,payload,payload_sha,contributions FROM frontiers ORDER BY root,batch"
        ):
            if root not in available_roots:
                raise ValueError("orphan frontier workload root")
            if root not in selected_roots:
                continue
            if type(batch) is not int or batch != len(seen_batches[root]):
                raise ValueError("frontier batches must be contiguous from zero")
            seen_batches[root].append(batch)
            frontier = decode(payload, sha)
            incoming = [tuple(row) for row in frontier["records"]]
            if frontier["batch"] != batch or len(incoming) != count:
                raise ValueError("frontier batch/count binding mismatch")
            local_ids = {row["state_id"] for row in frontier["states"]}
            if any((root, state) in states_seen for state in local_ids):
                raise ValueError("a state ID was reused across scheduler times")
            states_seen.update((root, state) for state in local_ids)
            if any(len(row) != 4 or any(type(x) is not int for x in row)
                   or row[0] != root or row[1] not in local_ids for row in incoming):
                raise ValueError("malformed or unbound contribution")
            if reference_collection(incoming) != tuple(tuple(row) for row in frontier["expected"]):
                raise ValueError("saved frontier collection mismatch")
            records.extend(incoming)
            expected.extend(tuple(row) for row in frontier["expected"])
        for root in selected_roots:
            if len(seen_batches[root]) != expected_batches[root]:
                raise ValueError("frontier batch coverage differs from certificate")
    if not records:
        raise ValueError("capture database contains no contributions")
    if reference_collection(records) != tuple(sorted(expected)):
        raise ValueError("aggregate collection crossed capture boundaries")
    metadata = {"workloads": workload_rows, "contributions": len(records),
                "groups_before_zero_removal": len({row[:3] for row in records}),
                "input_sha256": digest(records), "expected_output_sha256": digest(sorted(expected)),
                "max_coefficient_bits": max(abs(row[3]).bit_length() for row in records),
                "roots": len(workload_rows), "frontier_states": len(states_seen),
                "capture_completeness_binding": "intended manifest and completed capture marker" if intended is not None
                else "legacy: saved workload rows and certificate batch counts only"}
    return records, metadata


def benchmark(database: Path, *, backend: str, repeats: int = 3,
              device: int = 0, roots: set[int] | None = None,
              mode: str = "sequential", collector: Any = None) -> dict[str, Any]:
    if backend not in {"cpu", "cuda", "opencl"}:
        raise ValueError("unknown backend")
    if mode not in {"sequential", "pooled"}:
        raise ValueError("unknown collection mode")
    if repeats < 1:
        raise ValueError("repeats must be positive")
    from zipper_algebra.zipping_algorithm.gpu_collection import collect_cpu

    started = time.perf_counter()
    records, workload = load_frontiers(database, roots)
    if mode == "sequential":
        with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as db:
            chunks = [
                [tuple(row) for row in decode(payload, sha)["records"]]
                for root, payload, sha in db.execute(
                    "SELECT root,payload,payload_sha FROM frontiers ORDER BY root,batch"
                ) if roots is None or root in roots
            ]
    else:
        chunks = [records]
    loading = time.perf_counter() - started
    if backend == "cpu" and collector is not None:
        raise ValueError("CPU baseline cannot use a GPU collector")
    reused_collector = collector is not None
    if reused_collector and getattr(collector, "metadata", {}).get("backend") != backend:
        raise ValueError("reused collector backend differs from requested backend")
    setup = time.perf_counter()
    if backend != "cpu" and collector is None:
        from zipper_algebra.zipping_algorithm.gpu_collection import CudaCollector, OpenClCollector
        collector = (CudaCollector if backend == "cuda" else OpenClCollector)(device=device)
    setup_wall = 0.0 if reused_collector else time.perf_counter() - setup
    rows = []
    for iteration in range(repeats):
        started = time.perf_counter()
        output_rows = []
        total_timings: Counter[str] = Counter()
        chunk_rows = []
        for index, chunk in enumerate(chunks):
            chunk_started = time.perf_counter()
            if collector is None:
                output = collect_cpu(chunk)
                computed = time.perf_counter()
                reference = reference_collection(chunk)
                if output != reference:
                    raise ValueError("CPU result differs from independent sorted-group verification")
                ended = time.perf_counter()
                timings = {"collection": computed - chunk_started,
                           "verification": ended - computed, "total": ended - chunk_started}
                metadata = {"device": platform.processor(),
                            "arithmetic": "Python arbitrary-precision integers"}
            else:
                collected_result = collector.collect(chunk, verify=False)
                output = collected_result.output
                timings, metadata = dict(collected_result.timings), dict(collected_result.metadata)
                verify_started = time.perf_counter()
                reference = reference_collection(chunk)
                if output != reference:
                    raise ValueError("GPU result differs from independent sorted-group verification")
                timings["verification"] = time.perf_counter() - verify_started
                timings["total"] = time.perf_counter() - chunk_started
                metadata["external_verification"] = "same sorted-group Python-integer oracle as CPU baseline"
                metadata["verified"] = True
            output_rows.extend(output)
            total_timings.update(timings)
            chunk_rows.append({"chunk": index, "contributions": len(chunk),
                               "nonzero_groups": len(output), "timings": timings,
                               "metadata": metadata})
        verify_started = time.perf_counter()
        output_rows.sort()
        output_sha = digest(output_rows)
        if output_sha != workload["expected_output_sha256"]:
            raise ValueError("collected output differs from captured scheduler output")
        total_timings["captured_output_verification"] = time.perf_counter() - verify_started
        observed = time.perf_counter() - started
        total_timings["total"] = observed
        rows.append({"iteration": iteration, "timings": dict(total_timings),
                     "observed_call_wall": observed,
                     "output_sha256": output_sha, "nonzero_groups": len(output_rows),
                     "chunks": chunk_rows})
    result = {"schema": "zipper-gpu-benchmark-v1", "backend": backend, "mode": mode,
              "chunk_count": len(chunks), "largest_chunk_contributions": max(map(len, chunks)),
              "workload": workload, "loading_and_capture_integrity_wall": loading,
              "backend_setup_wall": setup_wall, "collector_reused": reused_collector,
              "timing_scope": "actual scheduler batches processed sequentially" if mode == "sequential" else "offline pooled captured batches; successive scheduler times are not simultaneously available online",
              "timing_includes": "packing, all transfers, coefficient collection, reconstruction, independent exact verification; capture/rewrite, full-continuation canonicalization/state-ID assignment, and database loading separate",
              "first_use_with_setup_wall": None if reused_collector else rows[0]["observed_call_wall"] + setup_wall,
              "retained_collector_initialization_seconds": getattr(collector, "initialization_seconds", None),
              "hardware": {"platform": platform.platform(), "python": platform.python_version()},
              "benchmark_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "collector_source_sha256": hashlib.sha256((Path(collected.__file__).parent / "gpu_collection.py").read_bytes()).hexdigest(),
              "rows": rows}
    with sqlite3.connect(database) as db:
        db.execute("INSERT INTO benchmarks(backend,result) VALUES(?,?)", (backend, json.dumps(result)))
    return result


def status(database: Path) -> dict[str, Any]:
    records, manifest = load_frontiers(database)
    with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as db:
        runs = [{"id": row[0], "backend": row[1], "created_utc": row[2]}
                for row in db.execute("SELECT id,backend,created_utc FROM benchmarks ORDER BY id")]
    return {"database": str(database.resolve()), **manifest,
            "nonzero_groups": len(reference_collection(records)), "benchmarks": runs}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    cap = sub.add_parser("capture")
    cap.add_argument("database", type=Path)
    cap.add_argument("--workloads", type=Path, help="JSON list of {name,n,word} objects")
    cap.add_argument("--max-steps", type=int, default=4096)
    cap.add_argument("--batch-size", type=int, default=64)
    bench = sub.add_parser("benchmark")
    bench.add_argument("database", type=Path)
    bench.add_argument("--backend", choices=["cpu", "cuda", "opencl"], required=True)
    bench.add_argument("--repeats", type=int, default=3)
    bench.add_argument("--device", type=int, default=0)
    bench.add_argument("--mode", choices=["sequential", "pooled"], default="sequential")
    bench.add_argument("--roots", nargs="+", type=int, help="Optional workload root IDs")
    stat = sub.add_parser("status")
    stat.add_argument("database", type=Path)
    args = parser.parse_args()
    if args.command == "capture":
        workloads = json.loads(args.workloads.read_text()) if args.workloads else DEFAULT_WORKLOADS
        result = capture(args.database, workloads, max_steps=args.max_steps, batch_size=args.batch_size)
    elif args.command == "benchmark":
        result = benchmark(args.database, backend=args.backend, repeats=args.repeats, device=args.device,
                           roots=set(args.roots) if args.roots else None, mode=args.mode)
    else:
        result = status(args.database)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
