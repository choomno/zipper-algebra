"""Compare synchronous and split certification using the production queue.

Workers initialize before a barrier. Timed work includes exact native calls,
transfers, candidate compression/storage, certification and queue coordination.
Setup, total CPU and sampled process-tree RSS are reported separately. The split
arm uses two concurrent workers; it is not evidence of reduced mathematical work.
Use paired2 as a two-worker synchronous comparison, or split12 for one producer
and two certifiers under the same inherited CPU affinity. No result is accepted merely
because its output digest agrees: comparison checks full canonical outputs.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import queue as queue_module
import resource
import sqlite3
import time
from typing import Any
import zlib

from experiments.zipper_async_census import initialize, read_policy, verify_result
from experiments.zipper_census_reuse_benchmark import PeakMemory
from zipper_algebra.zipping_algorithm.census_queue import CensusQueue
from zipper_algebra.zipping_algorithm.census_reuse import _canonical_digest, verifier_binding
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession, _decode_coefficient_object
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import completed_output


def cpu_seconds() -> float:
    return sum(getattr(resource.getrusage(kind), field)
               for kind in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN)
               for field in ("ru_utime", "ru_stime"))


def manifest(root: Path) -> dict[str, str]:
    paths = [p for folder in ("src", "native") for p in (root / folder).rglob("*")
             if p.is_file() and p.suffix in (".py", ".cpp", ".hpp", ".h", ".json")]
    paths += [root / "experiments" / name for name in
              ("zipper_async_benchmark.py", "zipper_async_census.py", "zipper_census_reuse_benchmark.py")]
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def validate_requests(requests: Any) -> list[list[Any]]:
    if not isinstance(requests, list) or not requests:
        raise ValueError("a nonempty request list is required")
    seen: set[str] = set()
    for item in requests:
        if (not isinstance(item, list) or len(item) != 2 or type(item[0]) is not int or item[0] < 3
                or not isinstance(item[1], list)
                or any(type(a) is not int or not 0 < abs(a) < item[0] for a in item[1])):
            raise ValueError("requests must be [rank>=3, signed word]")
        key = json.dumps(item)
        if key in seen:
            raise ValueError("benchmark requests must be distinct")
        seen.add(key)
    return requests


def _worker(path: Path, policy: dict[str, Any], role: str, config: dict[str, Any],
            ready: Any, start: Any, reports: Any) -> None:
    before = time.perf_counter()
    try:
        with ExitStack() as stack:
            queue = CensusQueue(path, policy)
            stack.callback(queue.close)
            queue.assert_sources()
            produce = verify = None
            verifier_identity: dict[str, Any] = {}
            if role in ("pair", "compute"):
                producer = stack.enter_context(NativeSession(cache_bytes=config["cache_bytes"]))
                reduce = producer.reduce_reused if policy["mode"] == "reuse" else producer.reduce
                limits = Limits(**policy["limits"])
                produce = lambda n, word: reduce(word, n=n, limits=limits)
            if role in ("pair", "certify"):
                verifier = stack.enter_context(CppVerifier(producer_metadata=policy["producer_metadata"],
                    cache_bytes=config["cache_bytes"], allocator_library=config["allocator"]))
                verify = lambda result: verify_result(verifier, result, policy["mode"])
                verifier_identity = verifier_binding(verifier.metadata)
            setup = time.perf_counter() - before
            ready.put({"role": role, "setup_seconds": setup})
            if not start.wait(config["timeout"]):
                raise TimeoutError("benchmark start barrier timed out")
            begin = time.perf_counter()
            deadline = begin + config["timeout"]
            while time.perf_counter() < deadline:
                did_work = False
                if produce is not None:
                    did_work |= queue.compute_one(produce, max_pending=config["max_pending"],
                                                  max_pending_bytes=config["max_pending_bytes"])
                if verify is not None:
                    did_work |= queue.certify_one(verify, verifier_identity)
                counts = queue.counts()
                states = ("queued", "computing") if role == "compute" else (
                    "queued", "computing", "candidate", "certifying")
                if not any(counts.get(s, 0) for s in states):
                    break
                if not did_work:
                    time.sleep(0.01)
            elapsed = time.perf_counter() - begin
            counts = queue.counts()
        reports.put({"role": role, "setup_seconds": setup, "work_seconds": elapsed,
                     "total_worker_cpu_seconds": cpu_seconds(), "counts_at_exit": counts})
    except BaseException as error:
        reports.put({"role": role, "error": repr(error)})
        ready.put({"role": role, "error": repr(error)})
        raise


def summarize_store(path: Path) -> dict[str, Any]:
    with sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True) as db:
        counts = dict(db.execute("SELECT state,count(*) FROM async_jobs GROUP BY state"))
        stages: dict[str, dict[str, float]] = {}
        for stage, status, record in db.execute("SELECT stage,status,record FROM async_attempts"):
            data = json.loads(record)
            stats = stages.setdefault(stage, {"attempts": 0, "seconds": 0})
            stats["attempts"] += 1
            stats["seconds"] += data.get("seconds", 0)
        size = db.execute("SELECT coalesce(sum(bytes),0) FROM async_candidates").fetchone()[0]
    return {"counts": counts, "stages": stages, "candidate_bytes": size,
            "database_bytes": sum(p.stat().st_size for p in path.parent.glob(path.name + "*") if p.is_file())}


def read_outputs(path: Path) -> dict[tuple[int, str], tuple[str, Any]]:
    outputs = {}
    with sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True) as db:
        for policy, key, n, word, state, payload, digest, output_digest in db.execute(
                "SELECT j.policy,j.key,j.n,j.word,j.state,c.payload,c.digest,c.output_digest FROM async_jobs j "
                "LEFT JOIN async_candidates c ON j.policy=c.policy AND j.key=c.key"):
            normal = None
            if state == "certified":
                if payload is None:
                    raise ValueError("certified candidate payload is missing")
                records = db.execute("SELECT record FROM async_attempts WHERE policy=? AND key=? "
                                     "AND stage='certify' AND status='certified' ORDER BY id DESC LIMIT 1",
                                     (policy, key)).fetchall()
                if not records:
                    raise ValueError("certified candidate has no certification attempt")
                record = json.loads(records[0][0])
                if record.get("digest") != digest:
                    raise ValueError("certification attempt digest mismatch")
                verifier_binding(record.get("verifier", {}))
                raw = zlib.decompress(payload)
                if hashlib.sha256(raw).hexdigest() != digest:
                    raise ValueError("candidate payload digest mismatch")
                certificate = json.loads(raw, object_hook=_decode_coefficient_object)
                if certificate.get("n") != n or certificate.get("input_word") != json.loads(word):
                    raise ValueError("candidate input binding mismatch")
                normal = completed_output(certificate)
                if _canonical_digest(normal) != output_digest:
                    raise ValueError("candidate output digest mismatch")
            key = (n, word)
            if key in outputs:
                raise ValueError("benchmark database contains multiple policies for one input")
            outputs[key] = (state, normal)
    return outputs


def compare(left: Path, right: Path) -> dict[str, Any]:
    a, b = read_outputs(left), read_outputs(right)
    if a.keys() != b.keys():
        raise ValueError("benchmark input sets differ")
    matches, unresolved = [], []
    for key in sorted(a):
        if a[key][0] == b[key][0] == "certified":
            if a[key][1] != b[key][1]:
                raise ValueError("independently certified full outputs differ")
            matches.append([key[0], json.loads(key[1])])
        else:
            unresolved.append({"n": key[0], "word": json.loads(key[1]), "left": a[key][0], "right": b[key][0]})
    return {"exact_output_matches": matches, "unresolved": unresolved,
            "all_inputs_certified_and_equal": bool(a) and not unresolved}


def receive(channel: Any, workers: list[Any], timeout: float) -> Any:
    """Notice a crashed worker while waiting for a different worker's message."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            return channel.get(timeout=min(0.2, max(0.001, deadline - time.monotonic())))
        except queue_module.Empty:
            if any(worker.exitcode not in (None, 0) for worker in workers):
                raise RuntimeError("benchmark worker exited unsuccessfully")
    raise queue_module.Empty


def run_arm(path: Path, requests: list[list[Any]], arm: str, config: dict[str, Any]) -> dict[str, Any]:
    if path.exists():
        raise ValueError("benchmark database must be new")
    roles = {"sequential": ["pair"], "paired2": ["pair", "pair"],
             "split": ["compute", "certify"], "split12": ["compute", "certify", "certify"]}[arm]
    setup_start, cpu_start = time.perf_counter(), cpu_seconds()
    initialize(path, Limits(**config["limits"]), config["mode"], config["cache_bytes"])
    policy = read_policy(path, None)
    queue = CensusQueue(path, policy)
    try:
        for n, word in requests:
            queue.enqueue(n, word)
    finally:
        queue.close()
    context = multiprocessing.get_context("spawn")
    ready, reports, start = context.Queue(), context.Queue(), context.Event()
    workers = [context.Process(target=_worker, args=(path, policy, role, config, ready, start, reports))
               for role in roles]
    try:
        with PeakMemory() as memory:
            for worker in workers:
                worker.start()
            readiness = [receive(ready, workers, config["timeout"]) for _ in workers]
            if any("error" in item for item in readiness):
                raise RuntimeError(str(readiness))
            setup = time.perf_counter() - setup_start
            begin = time.perf_counter()
            start.set()
            results = [receive(reports, workers, config["timeout"] + 60) for _ in workers]
            for worker in workers:
                worker.join(timeout=10)
            elapsed = time.perf_counter() - begin
            if any("error" in item for item in results) or any(w.exitcode != 0 for w in workers):
                raise RuntimeError(str(results))
        return {"arm": arm, "database": str(path), "setup_seconds": setup, "wall_seconds": elapsed,
                "wall_includes_worker_teardown": True,
                "max_worker_work_seconds": max(r["work_seconds"] for r in results),
                "total_cpu_seconds_including_setup": cpu_seconds() - cpu_start,
                "peak_tree_rss_bytes": memory.peak, "workers": results, "policy": policy,
                **summarize_store(path)}
    except queue_module.Empty as error:
        raise TimeoutError("benchmark worker did not finish within guard") from error
    finally:
        for worker in workers:
            if worker.is_alive():
                worker.terminate()
            if worker.pid is not None:
                worker.join(timeout=10)
        ready.close()
        reports.close()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--requests", type=Path, required=True)
    parser.add_argument("--order", default="sequential,split,split,sequential")
    parser.add_argument("--limits", default='{"max_steps":250000,"max_terms":20000,"expansion_budget":1000000}')
    parser.add_argument("--cache-mib", type=int, default=256)
    parser.add_argument("--mode", choices=("reuse", "standalone"), default="reuse")
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--allocator", type=Path)
    parser.add_argument("--cpus", help="optional comma-separated inherited CPU affinity")
    args = parser.parse_args(argv)
    order = args.order.split(",")
    if any(arm not in ("sequential", "paired2", "split", "split12") for arm in order):
        parser.error("order must contain sequential, paired2, split or split12")
    if args.cache_mib < 0 or args.timeout <= 0:
        parser.error("cache must be nonnegative and timeout positive")
    requests = validate_requests(json.loads(args.requests.read_text()))
    if args.cpus:
        os.sched_setaffinity(0, set(map(int, args.cpus.split(","))))
    config = {"mode": args.mode, "limits": json.loads(args.limits), "cache_bytes": args.cache_mib * 2**20,
              "allocator": args.allocator, "timeout": args.timeout,
              "max_pending": len(requests) + 1, "max_pending_bytes": 2**30}
    root = Path(__file__).resolve().parents[1]
    frozen = manifest(root)
    args.directory.mkdir(parents=True, exist_ok=True)
    destination = args.directory / "report.json"
    if destination.exists():
        raise ValueError("benchmark report already exists")
    report: dict[str, Any] = {"source_manifest": frozen, "requests": requests,
        "config": {**config, "allocator": str(args.allocator) if args.allocator else None},
        "affinity": sorted(os.sched_getaffinity(0)), "runs": [], "comparisons": []}
    for i, arm in enumerate(order):
        if manifest(root) != frozen:
            raise RuntimeError("source changed during benchmark")
        run = run_arm(args.directory / f"{i}-{arm}.sqlite3", requests, arm, config)
        report["runs"].append(run)
        if manifest(root) != frozen:
            raise RuntimeError("source changed during benchmark")
        destination.write_text(json.dumps(report, indent=2) + "\n")
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            report["comparisons"].append({"left": i, "right": j,
                **compare(Path(report["runs"][i]["database"]), Path(report["runs"][j]["database"]))})
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"report": str(destination), "runs": [
        {key: run[key] for key in ("arm", "wall_seconds", "counts")} for run in report["runs"]]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
