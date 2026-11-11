"""Measure exact zipper execution in one isolated process.

Example::

    python experiments/zipper_memoization_benchmark.py --execution memoized \
        --case negative_l10 --output results/negative --check-after-first
    python experiments/zipper_benchmark_artifacts.py results \
        --remove-verified-originals

Select the package version with PYTHONPATH before starting this process.
RSS is measured; optional tracemalloc runs should be separate from timing runs.
The archive command verifies bytes before removing loose certificates.
"""

import argparse, gc, hashlib, itertools, json, os, resource, time, tracemalloc
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--output", required=True)
p.add_argument("--case", required=True)
p.add_argument("--session-kwargs", default="{}")
p.add_argument("--passes", type=int, default=2)
p.add_argument("--trace-memory", action="store_true")
p.add_argument("--check-after-first", action="store_true")
p.add_argument("--rule-cache-bytes", type=int)
p.add_argument("--rule-cache-entries", type=int, default=256)
p.add_argument("--execution", choices=["collected", "memoized"], default="collected")
args = p.parse_args()
if args.execution == "memoized":
    from zipper_algebra.zipping_algorithm.memoized import (
        MemoizedSession as CollectedSession,
        _source,
    )
else:
    from zipper_algebra.zipping_algorithm.collected import CollectedSession, _source
from zipper_algebra.zipping_algorithm.interpreter import Limits


def resident_memory():
    values = {}
    for path in ["/proc/self/status", "/proc/self/smaps_rollup"]:
        try:
            for line in Path(path).read_text().splitlines():
                key = line.split(":")[0]
                if key in [
                    "VmRSS",
                    "VmHWM",
                    "Rss",
                    "Pss",
                    "Private_Clean",
                    "Private_Dirty",
                ]:
                    values[key + "_kib"] = int(line.split()[1])
        except OSError:
            pass
    return values


cases = {
    "negative_l10": [(3, [-1, -2] * 5)],
    "positive_l10": [(3, [1, 2] * 5)],
    "alternating_l10": [(3, [1, -2] * 5)],
    "positive_n5": [(5, [1, 2, 3, 4])],
    "corpus_n3_l3": [
        (3, list(w))
        for l in range(4)
        for w in itertools.product([-2, -1, 1, 2], repeat=l)
    ],
}
words = cases[args.case]
out = Path(args.output)
out.mkdir(parents=True, exist_ok=True)
if args.trace_memory:
    tracemalloc.start()
t0 = time.perf_counter()
session = CollectedSession(**json.loads(args.session_kwargs))
setup = time.perf_counter() - t0
if args.rule_cache_bytes is not None:
    session._catalog.configure_rule_cache(
        max_bytes=args.rule_cache_bytes, max_entries=args.rule_cache_entries
    )
result = {
    "case": args.case,
    "source": _source(),
    "session_kwargs": json.loads(args.session_kwargs),
    "setup_wall": setup,
    "pid": os.getpid(),
    "rss_unit": "KiB on Linux",
    "instrumented_tracemalloc": args.trace_memory,
    "passes": [],
}
for run in range(args.passes):
    records = []
    start = time.perf_counter()
    cpu = time.process_time()
    for index, (n, word) in enumerate(words):
        w = time.perf_counter()
        t = time.process_time()
        cert = session.reduce(
            word,
            n=n,
            limits=Limits(max_steps=100000, max_terms=100000, expansion_budget=1000000),
            trace=False,
            workers=1,
            batch_size=64,
        )
        rec = {
            "index": index,
            "n": n,
            "word": word,
            "status": cert["status"],
            "wall": time.perf_counter() - w,
            "cpu": time.process_time() - t,
            "metrics": cert["metrics"],
        }
        encoded = json.dumps(
            cert["output"], sort_keys=True, separators=(",", ":")
        ).encode()
        rec["output_sha256"] = hashlib.sha256(encoded).hexdigest()
        (out / f"pass{run}-case{index}.json").write_text(
            json.dumps(cert, separators=(",", ":"))
        )
        records.append(rec)
    row = {
        "pass": run,
        "wall_including_certificate_writes": time.perf_counter() - start,
        "cpu_including_certificate_writes": time.process_time() - cpu,
        "producer_wall": sum(r["wall"] for r in records),
        "producer_cpu": sum(r["cpu"] for r in records),
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "cases": records,
        "resident_memory": resident_memory(),
    }
    if args.trace_memory:
        row["tracemalloc_current_peak"] = tracemalloc.get_traced_memory()
    for key in ["cache_info", "cache_stats", "memo_stats", "statistics"]:
        val = getattr(session, key, None)
        if val is not None:
            try:
                row[key] = val() if callable(val) else val
            except TypeError:
                pass
    if hasattr(session._catalog, "rule_cache_info"):
        row["rule_cache_info"] = session._catalog.rule_cache_info()
    if args.check_after_first and run == 0:
        ct = time.perf_counter()
        cc = time.process_time()
        for cert_path in sorted(out.glob("pass0-case*.json")):
            assert session.check(json.loads(cert_path.read_text()))
        row["independent_check"] = {
            "wall": time.perf_counter() - ct,
            "cpu": time.process_time() - cc,
            "resident_memory": resident_memory(),
            "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        }
    result["passes"].append(row)
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
if hasattr(session, "clear_cache"):
    session.clear_cache()
    if hasattr(session._catalog, "clear_rule_cache"):
        session._catalog.clear_rule_cache()
    gc.collect()
    result["after_clear_resident_memory"] = resident_memory()
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
print(
    json.dumps(
        {
            "case": args.case,
            "passes": [
                {k: v for k, v in row.items() if k != "cases"}
                for row in result["passes"]
            ],
        }
    )
)
