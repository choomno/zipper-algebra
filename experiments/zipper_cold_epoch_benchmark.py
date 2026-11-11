"""Run the reviewed reuse benchmark against an immutable source epoch.

A run imports only its selected epoch. Cross-epoch reports retain both source
manifests and compare exact independently checked output bytes; differences in
source are intentional and never treated as within-epoch equivalence evidence.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import sqlite3
import sys
import time
from typing import Any, Callable
import zlib


_IMPORTED_DRIVER_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def assert_driver_source() -> None:
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != _IMPORTED_DRIVER_SHA256:
        raise RuntimeError("epoch driver source changed after import")


def encode(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def source_manifest(epoch: Path) -> dict[str, str]:
    paths = [epoch / "experiments/zipper_census_reuse_benchmark.py"]
    for directory in (epoch / "src", epoch / "native"):
        paths.extend(p for p in directory.rglob("*")
                     if p.is_file() and p.suffix in (".py", ".cpp", ".hpp", ".h", ".json"))
    return {str(p.relative_to(epoch)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(paths)}


def conditions() -> dict[str, Any]:
    result: dict[str, Any] = {"affinity": sorted(os.sched_getaffinity(0)), "loadavg": os.getloadavg()}
    for name in ("cpu", "memory", "io"):
        path = Path("/proc/pressure") / name
        result[f"pressure_{name}"] = path.read_text() if path.exists() else None
    return result


def cpu_seconds() -> float:
    return sum(getattr(resource.getrusage(kind), field)
               for kind in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)
               for field in ("ru_utime", "ru_stime"))


def validate_cold_cases(cases: list[dict[str, Any]], config: dict[str, Any], cache_mode: str,
                        classify: Callable[[list[int], int], dict[str, Any]]) -> None:
    if cache_mode != "fresh" or config["passes"] != 1:
        raise ValueError("cold cohort requires fresh sessions and one pass")
    seen_requests: set[str] = set()
    for case in cases:
        if config.get("labels") and case["label"] not in config["labels"]:
            continue
        classification = classify(case["word"], case["n"])
        requests = classification["unique_primitive_requests"]
        if not classification["unchanged_primitive"] or len(requests) != 1:
            raise ValueError("cold cohort requires unchanged primitive inputs")
        request = requests[0]["key"]
        if request in seen_requests:
            raise ValueError("cold cohort repeats a primitive request")
        seen_requests.add(request)


def run_epoch(epoch: Path, database: Path, cases: list[dict[str, Any]],
              backend: str, cache_mode: str, config: dict[str, Any], *,
              epoch_name: str, cohort: str, round_index: int, cpu: int) -> int:
    """Run one epoch in a fresh Python process (one call per invocation)."""
    assert_driver_source()
    epoch = epoch.resolve()
    if any(name == "zipper_algebra" or name.startswith("zipper_algebra.")
           or name == "experiments.zipper_census_reuse_benchmark" for name in sys.modules):
        raise RuntimeError("an epoch was already imported; use a fresh process")
    package = sys.modules.get("experiments")
    if package is not None and any(Path(path).resolve() != epoch / "experiments"
                                   for path in getattr(package, "__path__", ())):
        raise RuntimeError("foreign experiments package was imported; invoke the driver by absolute file path")
    if cpu not in os.sched_getaffinity(0):
        raise ValueError("requested CPU is outside allowed affinity")
    os.sched_setaffinity(0, {cpu})
    resource.setrlimit(resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3))
    sys.path[:0] = [str(epoch), str(epoch / "src")]
    frozen = source_manifest(epoch)
    bench = importlib.import_module("experiments.zipper_census_reuse_benchmark")
    if bench.__file__ is None or Path(bench.__file__).resolve() != epoch / "experiments/zipper_census_reuse_benchmark.py":
        raise RuntimeError("benchmark was imported from a different epoch")
    for name, module in list(sys.modules.items()):
        if name == "zipper_algebra" or name.startswith("zipper_algebra."):
            filename = getattr(module, "__file__", None)
            if filename is not None and not Path(filename).resolve().is_relative_to(epoch / "src"):
                raise RuntimeError("zipper module was imported from a different epoch: " + name)
    if cohort == "cold":
        validate_cold_cases(cases, config, cache_mode, bench.input_classification)
    database.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(database) as db:
        bench.schema(db)
        db.execute("""CREATE TABLE IF NOT EXISTS cold_epoch_runs(
            run INTEGER PRIMARY KEY, epoch_name TEXT NOT NULL, epoch_path TEXT NOT NULL,
            cohort TEXT NOT NULL, round_index INTEGER NOT NULL, source_manifest TEXT NOT NULL,
            driver_sha256 TEXT NOT NULL, conditions_before TEXT NOT NULL,
            conditions_after TEXT, process_seconds REAL, process_cpu_seconds REAL)""")
        expected = [(i, case["label"], case["n"], encode(case["word"]), encode(case.get("origin", {})))
                    for i, case in enumerate(cases, 1)]
        actual = db.execute("SELECT * FROM reuse_bench_cases ORDER BY id").fetchall()
        if not actual:
            db.executemany("INSERT INTO reuse_bench_cases VALUES(?,?,?,?,?)", expected)
        elif actual != expected:
            raise ValueError("epoch input set differs from frozen SQLite inputs")
    config = dict(config, epoch_name=epoch_name, cohort=cohort, round_index=round_index,
                  process_address_space_bytes=4 * 1024**3, cpu=cpu)
    run = bench.create_run(database, backend, cache_mode, config)
    before, start, cpu_start = conditions(), time.perf_counter(), cpu_seconds()
    with sqlite3.connect(database) as db:
        db.execute("INSERT INTO cold_epoch_runs VALUES(?,?,?,?,?,?,?,?,NULL,NULL,NULL)",
                   (run, epoch_name, str(epoch), cohort, round_index, encode(frozen),
                    _IMPORTED_DRIVER_SHA256, encode(before)))
    try:
        bench.worker(database, run)
        if cohort == "cold" and backend in ("census", "census-standalone"):
            with sqlite3.connect(database) as db:
                for (telemetry,) in db.execute("SELECT telemetry FROM reuse_bench_measurements WHERE run=?", (run,)):
                    if json.loads(telemetry).get("stats", {}).get("leaf_hits", 0):
                        db.execute("UPDATE reuse_bench_runs SET status='failed',error='cold sidecar hit' WHERE id=?", (run,))
                        db.commit()
                        raise ValueError("cold cohort unexpectedly reused a persistent leaf")
    finally:
        with sqlite3.connect(database) as db:
            db.execute("UPDATE cold_epoch_runs SET conditions_after=?,process_seconds=?,process_cpu_seconds=? WHERE run=?",
                       (encode(conditions()), time.perf_counter() - start, cpu_seconds() - cpu_start, run))
            if (source_manifest(epoch) != frozen
                    or hashlib.sha256(Path(__file__).read_bytes()).hexdigest() != _IMPORTED_DRIVER_SHA256):
                db.execute("UPDATE reuse_bench_runs SET status='failed',error='epoch source drift' WHERE id=?", (run,))
                db.commit()
                raise RuntimeError("epoch source changed while measuring")
    return run


def pair_report(database: Path, left_run: int, right_run: int) -> dict[str, Any]:
    """Compare completed intersections of two explicitly selected epoch runs."""
    with sqlite3.connect(f"file:{database}?mode=ro", uri=True) as db:
        cases = [{"label": label, "n": n, "word": json.loads(word), "origin": json.loads(origin)}
                 for label, n, word, origin in db.execute(
                     "SELECT label,n,word,origin FROM reuse_bench_cases ORDER BY id")]
        frozen = db.execute("SELECT value FROM reuse_bench_metadata WHERE key='cases_sha256'").fetchone()
        if frozen is None or json.loads(frozen[0]) != hashlib.sha256(encode(cases).encode()).hexdigest():
            raise ValueError("frozen epoch benchmark inputs differ")
        case_ids = {row[0] for row in db.execute("SELECT id FROM reuse_bench_cases")}
        runs = {}
        measurements = {}
        for run in (left_run, right_run):
            row = db.execute("""SELECT r.backend,r.cache_mode,r.status,r.config,r.identity,e.epoch_name,
                e.source_manifest,e.driver_sha256,e.cohort,e.round_index,e.conditions_before,e.conditions_after,
                r.setup_seconds,e.process_seconds,e.process_cpu_seconds
                FROM reuse_bench_runs r JOIN cold_epoch_runs e ON e.run=r.id WHERE r.id=?""", (run,)).fetchone()
            if row is None or row[2] != "complete":
                raise ValueError("epoch run is absent or did not finish successfully")
            runs[run] = dict(zip(("backend", "cache_mode", "status", "config", "identity", "epoch_name",
                "source_manifest", "driver_sha256", "cohort", "round_index", "conditions_before",
                "conditions_after", "setup_seconds", "process_seconds", "process_cpu_seconds"), row))
            for key in ("config", "identity", "source_manifest", "conditions_before", "conditions_after"):
                runs[run][key] = json.loads(runs[run][key])
            measurements[run] = {(repeat, case): {"status": status, "seconds": seconds,
                "producer_seconds": producer, "checker_seconds": checker, "store_seconds": storage,
                "reconstruction_seconds": reconstruction, "peak_tree_rss_bytes": memory,
                "telemetry": json.loads(telemetry), "output": output, "sha": sha}
                for repeat, case, status, seconds, producer, checker, storage, reconstruction, memory,
                    telemetry, output, sha in db.execute("""SELECT m.pass,m.case_id,m.status,m.seconds,
                    m.producer_seconds,m.checker_seconds,m.store_seconds,m.reconstruction_seconds,
                    m.peak_tree_rss_bytes,m.telemetry,p.output,p.output_sha
                    FROM reuse_bench_measurements m LEFT JOIN reuse_bench_payloads p
                    ON (m.run,m.pass,m.case_id)=(p.run,p.pass,p.case_id) WHERE m.run=?""", (run,))}
        left, right = runs[left_run], runs[right_run]
        for measurement_set in measurements.values():
            if any(case_id not in case_ids for _, case_id in measurement_set):
                raise ValueError("measurement references an unknown frozen case")
        worker_source = "experiments/zipper_census_reuse_benchmark.py"
        if (not left["source_manifest"].get(worker_source)
                or left["source_manifest"].get(worker_source) != right["source_manifest"].get(worker_source)):
            raise ValueError("epoch benchmark worker sources differ or are unbound")
        if left["cache_mode"] != right["cache_mode"]:
            raise ValueError("cache modes differ")
        for key in ("limits", "cache_bytes", "cache_entries", "case_seconds", "allocator", "cpu",
                    "process_address_space_bytes", "labels", "passes"):
            if left["config"].get(key) != right["config"].get(key):
                raise ValueError(f"comparison configuration differs: {key}")
        if left["driver_sha256"] != right["driver_sha256"]:
            raise ValueError("epoch driver sources differ")
        matches, excluded = [], []
        for case_key in sorted(set(measurements[left_run]) | set(measurements[right_run])):
            a, b = measurements[left_run].get(case_key), measurements[right_run].get(case_key)
            if a is None or b is None or a["status"] != "complete" or b["status"] != "complete":
                excluded.append({"pass": case_key[0], "case_id": case_key[1],
                    "left": None if a is None else a["status"], "right": None if b is None else b["status"]})
                continue
            for record in (a, b):
                if record["output"] is None or not record["telemetry"]["stats"].get("independently_checked"):
                    raise ValueError("complete output lacks independently checked evidence")
                raw = zlib.decompress(record["output"])
                if hashlib.sha256(raw).hexdigest() != record["sha"]:
                    raise ValueError("stored output hash differs")
                if json.loads(raw).get("completion") != "complete":
                    raise ValueError("complete output has incomplete canonical payload")
                record["raw"] = raw
            if a["raw"] != b["raw"]:
                raise AssertionError(f"exact cross-epoch output mismatch at {case_key}")
            if a["telemetry"]["input_classification"] != b["telemetry"]["input_classification"]:
                raise AssertionError(f"input normal-form classification changed at {case_key}")
            metrics = ("seconds", "producer_seconds", "checker_seconds", "store_seconds",
                       "reconstruction_seconds", "peak_tree_rss_bytes")
            matches.append({"pass": case_key[0], "case_id": case_key[1], "output_sha256": a["sha"],
                            "left": {k: a[k] for k in metrics}, "right": {k: b[k] for k in metrics}})
    lt = sum(item["left"]["seconds"] for item in matches)
    rt = sum(item["right"]["seconds"] for item in matches)
    comparable = all(left[key] == right[key] for key in ("backend", "cohort", "round_index"))
    return {"scope": ("explicit cross-epoch optimization comparison" if comparable else
                       "exact reference comparison across different backend/cohort/round axes"),
            "optimization_timing_comparable": comparable,
            "left_run": left_run, "right_run": right_run, "epochs": runs,
            "matches": matches, "excluded": excluded, "left_seconds": lt, "right_seconds": rt,
            "left_over_right": lt / rt if comparable and rt else None}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("--epoch", type=Path)
    parser.add_argument("--epoch-name")
    parser.add_argument("--cases", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--cohort", default="cold")
    parser.add_argument("--round", type=int, default=0)
    parser.add_argument("--cpu", type=int, default=6)
    parser.add_argument("--backend", choices=("standalone", "dag", "census", "census-standalone"), default="census")
    parser.add_argument("--cache-mode", choices=("fresh", "retained"), default="fresh")
    parser.add_argument("--pair", type=int, nargs=2)
    args = parser.parse_args()
    if args.pair:
        print(encode(pair_report(args.database, *args.pair)))
    else:
        if not all((args.epoch, args.epoch_name, args.cases, args.config)):
            parser.error("epoch execution requires --epoch, --epoch-name, --cases, and --config")
        run_epoch(args.epoch, args.database, json.loads(args.cases.read_text()), args.backend,
                  args.cache_mode, json.loads(args.config.read_text()), epoch_name=args.epoch_name,
                  cohort=args.cohort, round_index=args.round, cpu=args.cpu)


if __name__ == "__main__":
    main()
