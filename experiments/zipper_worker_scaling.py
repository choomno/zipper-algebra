"""Fixed distinct-word worker scaling; run while the production pool is stopped."""

import argparse
import concurrent.futures
import hashlib
import json
import multiprocessing
import os
import random
import resource
import sqlite3
import threading
import time
import zlib
from pathlib import Path

BASE = Path(".")
DB = Path("words.sqlite3")


def canonical(output):
    return hashlib.sha256(
        json.dumps(
            sorted(
                json.dumps(x, sort_keys=True, separators=(",", ":")) for x in output
            ),
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def prepare(
    database=None,
    output=None,
    *,
    strata=((3, 5), (3, 6), (4, 4), (5, 4)),
    measured_per_stratum=64,
    warmup_per_stratum=8,
    seed=20261106,
):
    global BASE, DB
    if database is not None:
        DB = Path(database).resolve()
    if output is not None:
        BASE = Path(output).resolve()
    BASE.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    groups = []
    with sqlite3.connect(DB.as_uri() + "?mode=ro", uri=True) as db:
        db.execute("BEGIN")
        for n, L in strata:
            rows = db.execute(
                'SELECT ordinal,word,certificate,certificate_sha FROM words WHERE n=? AND length=? AND status="complete" ORDER BY ordinal',
                (n, L),
            ).fetchall()
            if len(rows) != (2 * (n - 1)) ** L:
                raise ValueError(f"incomplete exact layer n={n}, L={L}")
            cases = []
            for ordinal, word, payload, digest in rng.sample(
                rows, measured_per_stratum + warmup_per_stratum
            ):
                raw = zlib.decompress(payload)
                if hashlib.sha256(raw).hexdigest() != digest:
                    raise ValueError("saved certificate hash mismatch")
                cert = json.loads(raw)
                cases.append(
                    {
                        "n": n,
                        "L": L,
                        "ordinal": ordinal,
                        "word": json.loads(word),
                        "output_sha": canonical(cert["output"]),
                        "steps": cert["metrics"]["steps"],
                        "source_certificate_sha": digest,
                    }
                )
            groups.append(cases)
    result = {
        "seed": seed,
        "database": str(DB),
        "warmup": [g[i] for i in range(warmup_per_stratum) for g in groups],
        "measured": [
            g[i]
            for i in range(
                warmup_per_stratum, warmup_per_stratum + measured_per_stratum
            )
            for g in groups
        ],
    }
    (BASE / "worker-scaling-inputs.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    return result


SESSION = None


def init(cache_bytes, cache_entries):
    global SESSION
    from zipper_algebra.zipping_algorithm.memoized import MemoizedSession

    SESSION = MemoizedSession(cache_bytes=cache_bytes, cache_entries=cache_entries)


def job(case):
    from zipper_algebra.zipping_algorithm.interpreter import Limits

    before = SESSION.cache_info()
    w = time.perf_counter()
    cpu = time.process_time()
    cert = SESSION.reduce(
        case["word"],
        n=case["n"],
        limits=Limits(max_steps=100000, max_terms=100000, expansion_budget=1000000),
        trace=False,
        workers=1,
        batch_size=64,
    )
    generated = time.perf_counter()
    gen_cpu = time.process_time() - cpu
    if not (
        cert["status"] == "complete"
        and cert["metrics"]["steps"] == case["steps"]
        and canonical(cert["output"]) == case["output_sha"]
    ):
        raise ValueError("benchmark result differs from saved exact output")
    SESSION.check(json.loads(json.dumps(cert)))
    after = SESSION.cache_info()
    return {
        "n": case["n"],
        "L": case["L"],
        "ordinal": case["ordinal"],
        "output_sha": case["output_sha"],
        "steps": case["steps"],
        "generate_wall": generated - w,
        "generate_cpu": gen_cpu,
        "checked_wall": time.perf_counter() - w,
        "checked_cpu": time.process_time() - cpu,
        "worker_pid": os.getpid(),
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "cache_delta": {
            k: after[k] - before[k]
            for k in [
                "evaluations",
                "hits",
                "batch_reuses",
                "evictions",
                "cross_input_hits",
            ]
        },
        "cache_bytes": after["total_encoded_bytes"],
        "cache_entries": after["entries"],
    }


def memory(pids):
    totals = {"rss_kib": 0, "pss_kib": 0, "processes": 0}
    for pid in pids:
        try:
            text = Path(f"/proc/{pid}/smaps_rollup").read_text()
            row = {
                line.split(":")[0]: int(line.split()[1])
                for line in text.splitlines()
                if line.startswith(("Rss:", "Pss:"))
            }
            totals["rss_kib"] += row["Rss"]
            totals["pss_kib"] += row["Pss"]
            totals["processes"] += 1
        except (OSError, KeyError):
            pass
    return totals


def run(count, inputs, *, cache_bytes=512 * 1024**2, cache_entries=262144):
    from zipper_algebra.zipping_algorithm.memoized import _source

    setup = time.perf_counter()
    ctx = multiprocessing.get_context("spawn")
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=count,
        mp_context=ctx,
        initializer=init,
        initargs=(cache_bytes, cache_entries),
    ) as pool:
        warm = list(pool.map(job, inputs["warmup"], chunksize=1))
        warmtime = time.perf_counter() - setup
        pids = [os.getpid(), *pool._processes]
        tracker = multiprocessing.resource_tracker._resource_tracker._pid
        if tracker:
            pids.append(tracker)
        stop = threading.Event()
        samples = []

        def sample():
            while not stop.is_set():
                samples.append(memory(pids))
                stop.wait(0.25)

        monitor = threading.Thread(target=sample)
        monitor.start()
        t = time.perf_counter()
        cpu = time.process_time()
        try:
            rows = list(pool.map(job, inputs["measured"], chunksize=1))
        finally:
            elapsed = time.perf_counter() - t
            parentcpu = time.process_time() - cpu
            stop.set()
            monitor.join()
        samples.append(memory(pids))
    result = {
        "workers": count,
        "source": _source(),
        "cache_bytes_per_worker": cache_bytes,
        "cache_entries_per_worker": cache_entries,
        "warmup_words": len(warm),
        "warmup_and_startup_wall": warmtime,
        "measured_words": len(rows),
        "measured_wall": elapsed,
        "words_per_second": len(rows) / elapsed,
        "logical_steps_per_second": sum(x["steps"] for x in rows) / elapsed,
        "checked_worker_cpu": sum(x["checked_cpu"] for x in rows),
        "parent_cpu": parentcpu,
        "max_sampled_rss_mib": max(x["rss_kib"] for x in samples) / 1024,
        "max_sampled_pss_mib": max(x["pss_kib"] for x in samples) / 1024,
        "memory_sample_count": len(samples),
        "cache_delta": {
            k: sum(x["cache_delta"][k] for x in rows) for k in rows[0]["cache_delta"]
        },
        "rows": rows,
        "warmup": warm,
        "memory_scope": "Parent, resource tracker and worker pool; sampled every250ms; not an exact peak.",
    }
    path = BASE / f"worker-scaling-w{count}.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    with sqlite3.connect(BASE / "benchmark-results.sqlite3") as db:
        db.execute(
            "CREATE TABLE IF NOT EXISTS worker_scaling(workers INTEGER PRIMARY KEY, result TEXT NOT NULL)"
        )
        db.execute(
            "INSERT OR REPLACE INTO worker_scaling VALUES (?,?)",
            (count, json.dumps(result)),
        )
    print(
        json.dumps(
            {k: v for k, v in result.items() if k not in ["rows", "warmup", "source"]}
        ),
        flush=True,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--database", type=Path)
    p.add_argument(
        "--strata", type=json.loads, default=[[3, 5], [3, 6], [4, 4], [5, 4]]
    )
    p.add_argument("--measured-per-stratum", type=int, default=64)
    p.add_argument("--warmup-per-stratum", type=int, default=8)
    p.add_argument("--seed", type=int, default=20261106)
    p.add_argument("--cache-mib", type=int, default=512)
    p.add_argument("--cache-entries", type=int, default=262144)
    p.add_argument("--prepare-only", action="store_true")
    p.add_argument("--workers", nargs="+", type=int, default=[4, 8, 12, 16, 24])
    a = p.parse_args()
    BASE = a.output.resolve()
    BASE.mkdir(parents=True, exist_ok=True)
    if a.prepare_only:
        if a.database is None:
            p.error("--database is required with --prepare-only")
        prepare(
            a.database,
            BASE,
            strata=a.strata,
            measured_per_stratum=a.measured_per_stratum,
            warmup_per_stratum=a.warmup_per_stratum,
            seed=a.seed,
        )
    else:
        inputs = json.loads((BASE / "worker-scaling-inputs.json").read_text())
        for count in a.workers:
            run(
                count,
                inputs,
                cache_bytes=a.cache_mib * 1024**2,
                cache_entries=a.cache_entries,
            )
