"""Durable asynchronous native computation and independent certification.

Run ``python -m experiments.zipper_async_census --help``. Multiple worker
processes may share one sidecar. Historical stores are never imported implicitly.
"""
from __future__ import annotations

import argparse
from contextlib import closing
from dataclasses import asdict, replace
import hashlib
import json
import multiprocessing
from collections import deque
from pathlib import Path
import sqlite3
import time
import uuid
from typing import Any

from zipper_algebra.zipping_algorithm.census_queue import CensusQueue, SCHEMA, LostLease
from zipper_algebra.zipping_algorithm.census_reuse import reconstruction_sources, verifier_binding, _sha
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier, VerifierError, VerifierIncomplete
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import execution_identity
from zipper_algebra.zipping_algorithm.peeling_enumeration import NonPeelableWords


def read_policy(path: Path, key: str | None) -> dict[str, Any]:
    with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)) as db:
        rows = db.execute("SELECT key,payload FROM async_policies" + (" WHERE key=?" if key else ""),
                          (key,) if key else ()).fetchall()
    if len(rows) != 1:
        raise ValueError("select exactly one existing policy with --policy")
    policy = json.loads(rows[0][1])
    if policy.get("schema") != SCHEMA or _sha(policy) != rows[0][0]:
        raise ValueError("stored policy schema or digest mismatch")
    return policy


def read_counts(path: Path, policy: dict[str, Any]) -> dict[str, int]:
    """Read historical coverage without authenticating a current runtime."""
    with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)) as db:
        counts = dict(db.execute("SELECT state,count(*) FROM async_jobs WHERE policy=? GROUP BY state", (_sha(policy),)))
        counts["stored_candidate_bytes"] = db.execute(
            "SELECT coalesce(sum(bytes),0) FROM async_candidates WHERE policy=?", (_sha(policy),)).fetchone()[0]
    return counts


def initialize(path: Path, limits: Limits, mode: str, cache_bytes: int) -> str:
    with NativeSession(cache_bytes=cache_bytes) as producer:
        reduce = producer.reduce_reused if mode == "reuse" else producer.reduce
        seed = reduce([], n=1, limits=replace(limits, max_steps=max(2, limits.max_steps),
                                            max_terms=max(2, limits.max_terms)))
        metadata = {k: v for k, v in producer.metadata.items()
                    if k not in ("compiled_during_initialization", "setup_seconds")}
        policy = {"schema": SCHEMA, "producer": execution_identity(seed),
                  "producer_metadata": metadata, "limits": asdict(limits), "mode": mode,
                  "package_sources": reconstruction_sources(),
                  "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        queue = CensusQueue(path, policy)
        try:
            return queue.policy_key
        finally:
            queue.close()


def verify_result(verifier: Any, result: dict[str, Any], mode: str) -> bool:
    """Preserve resource/operational exceptions; classify completed mismatches."""
    check = verifier.check_reused if mode == "reuse" else verifier.check
    try:
        return check(result)
    except VerifierIncomplete:
        raise
    except VerifierError:
        report = verifier.last_result
        if (report.get("ok") is False and (report.get("replay_status") == "complete"
                or report.get("error") == "native certificate differs from independent C++ replay")):
            return False
        raise


def worker(path: Path, policy: dict[str, Any], role: str, settings: dict[str, Any]) -> None:
    queue = CensusQueue(path, policy)
    session: Any = None
    attempt_metadata: dict[str, Any] = {"worker_session": str(uuid.uuid4()), "cache_bytes": settings["cache_bytes"],
                        "cache_entries": settings.get("cache_entries", 262144) if role == "certify" else 65536,
                        "cache_lifetime": "worker_session", "attempt_index": 0}
    try:
        queue.assert_sources()
        if policy["runner_sha256"] != hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
            raise ValueError("runner differs from policy")
        if role == "compute":
            session = NativeSession(cache_bytes=settings["cache_bytes"],
                                    cache_entries=attempt_metadata["cache_entries"])
            if any(session.metadata.get(key) != policy["producer_metadata"].get(key)
                   for key in ("source_sha256", "binary_sha256", "wrapper_source_sha256", "abi",
                               "compiler", "compiler_version", "flags", "link_libraries")):
                raise ValueError("producer build differs from queue policy")
            reduce = session.reduce_reused if policy["mode"] == "reuse" else session.reduce
            limits = Limits(**policy["limits"])
            def step() -> bool:
                return queue.compute_one(lambda n, word: reduce(word, n=n, limits=limits),
                                         attempt_metadata=attempt_metadata,
                                         max_pending=settings["max_pending"],
                                         max_pending_bytes=settings["max_pending_bytes"])
        else:
            session = CppVerifier(producer_metadata=policy["producer_metadata"],
                                  cache_bytes=settings["cache_bytes"],
                                  cache_entries=settings.get("cache_entries", 262144),
                                  allocator_library=settings["allocator"])
            def verify(result: dict[str, Any]) -> bool:
                return verify_result(session, result, policy["mode"])
            identity = verifier_binding(session.metadata)
            def step() -> bool:
                return queue.certify_one(verify, identity, attempt_metadata=attempt_metadata,
                                         verification_telemetry=lambda: {
                                             "report": session.last_result,
                                             "native_call_seconds": session.last_call_seconds})
        deadline = time.monotonic() + settings["wall_seconds"]
        while time.monotonic() < deadline and not (settings.get("stop") and settings["stop"].is_set()):
            try:
                did_work = step()
                if did_work:
                    attempt_metadata["attempt_index"] += 1
            except LostLease:
                attempt_metadata["attempt_index"] += 1
                did_work = True
            if role == "certify" and session.pid is None:
                raise RuntimeError("verifier process exited; remaining candidates preserved")
            if settings["once"]:
                break
            if not did_work:
                counts = queue.work_presence()
                if settings["drain"]:
                    remaining = (sum(counts.get(s, False) for s in ("queued", "computing")) if role == "compute"
                                 else sum(counts.get(s, False) for s in ("queued", "computing", "candidate", "certifying")))
                    if remaining == 0:
                        break
                    if (role == "compute" and counts.get("queued", False)
                            and not sum(counts.get(s, False) for s in ("computing", "candidate", "certifying"))
                            and sum(counts.get(s, False) for s in ("verification_incomplete", "verification_error"))):
                        print(json.dumps({"blocked": "uncertified backlog requires retry", "present": counts}), flush=True)
                        if settings.get("stop"):
                            settings["stop"].set()
                        break
                time.sleep(0.2)
    finally:
        if session is not None:
            session.close()
        queue.close()


def admit(queue: CensusQueue, targets: dict[str, int], *, maximum: int, minimum_length: int = 0) -> int:
    """Round-robin bounded admission of primitive rank >=3 layers."""
    if type(minimum_length) is not int or minimum_length < 0:
        raise ValueError("minimum length must be a nonnegative integer")
    layers: list[NonPeelableWords] = []
    for rank, cutoff in sorted((int(n), value) for n, value in targets.items()):
        if rank < 3 or type(cutoff) is not int or cutoff < 0:
            raise ValueError("census targets require rank >=3 and nonnegative integer cutoffs")
        layers.extend(NonPeelableWords(rank, length) for length in range(max(rank - 1, minimum_length), cutoff + 1))
    rotation_key = "admission-rotation:" + _sha({"targets": targets, "minimum_length": minimum_length})
    row = queue.db.execute("SELECT ordinal FROM async_cursors WHERE policy=? AND layer=?",
                           (queue.policy_key, rotation_key)).fetchone()
    order = deque(range(len(layers)))
    if row and order:
        order.rotate(-((int(row[0]) + 1) % len(layers)))
    admitted = 0
    while order and admitted < maximum:
        index = order.popleft()
        layer = layers[index]
        def request_at(i: int, layer: NonPeelableWords = layer) -> tuple[int, list[int]]:
            return layer.n, list(layer.word_at(i))
        count = queue.enqueue_interval(f"nonpeelable-v1:{layer.n}:{layer.length}", layer.count,
                  request_at, batch=min(64, maximum - admitted),
                  rotation=(rotation_key, str(index)))
        admitted += count
        if count:
            order.append(index)
    return admitted


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("init", "enqueue", "compute", "certify", "run", "status", "retry"))
    parser.add_argument("database", type=Path)
    parser.add_argument("--policy")
    parser.add_argument("--limits", default='{}')
    parser.add_argument("--mode", choices=("reuse", "standalone"), default="reuse")
    parser.add_argument("--words", type=Path, help="JSON array of [rank, literal word] requests")
    parser.add_argument("--targets", help='primitive rank>=3 cutoffs, e.g. {"3":6,"4":4}')
    parser.add_argument("--min-length", type=int, default=0, help="lowest length to admit from --targets")
    parser.add_argument("--admit-max", type=int, default=10000)
    parser.add_argument("--cache-mib", type=int, default=256)
    parser.add_argument("--cache-entries", type=int, default=262144, help="verifier cache entry budget")
    parser.add_argument("--max-pending", type=int, default=256)
    parser.add_argument("--max-pending-mib", type=int, default=1024)
    parser.add_argument("--compute-workers", type=int, default=1)
    parser.add_argument("--certify-workers", type=int, default=1)
    parser.add_argument("--wall-seconds", type=float, default=3600)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--drain", action="store_true", help="exit once admitted work drains; time limit still applies")
    parser.add_argument("--allocator", type=Path)
    parser.add_argument("--state", choices=("compute_incomplete", "compute_error", "verification_incomplete", "verification_error"))
    args = parser.parse_args(argv)
    if (args.min_length < 0 or args.cache_mib < 0 or args.cache_entries < 0 or min(args.max_pending, args.max_pending_mib, args.compute_workers,
                                  args.certify_workers, args.admit_max) < 1 or args.wall_seconds <= 0):
        parser.error("budgets and worker counts must be positive; cache and minimum length may be zero")
    if args.command == "init":
        key = initialize(args.database, Limits(**json.loads(args.limits)), args.mode, args.cache_mib * 2**20)
        print(json.dumps({"policy": key}))
        return 0
    policy = read_policy(args.database, args.policy)
    if args.command == "status":
        print(json.dumps({"policy": _sha(policy), "counts": read_counts(args.database, policy)}, sort_keys=True))
        return 0
    queue = CensusQueue(args.database, policy)
    try:
        if args.command == "enqueue":
            queue.assert_sources()
            if policy["runner_sha256"] != hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
                raise ValueError("runner differs from policy")
            if bool(args.words) == bool(args.targets):
                parser.error("enqueue requires exactly one of --words or --targets")
            if args.words:
                requests = json.loads(args.words.read_text())
                for n, word in requests:
                    queue.enqueue(n, word)
            else:
                admit(queue, json.loads(args.targets), maximum=args.admit_max, minimum_length=args.min_length)
        elif args.command == "retry":
            if args.state is None:
                parser.error("retry requires --state")
            queue.retry(args.state)
        elif args.command in ("compute", "certify", "run"):
            settings = {"cache_bytes": args.cache_mib * 2**20, "cache_entries": args.cache_entries, "max_pending": args.max_pending,
                        "max_pending_bytes": args.max_pending_mib * 2**20, "wall_seconds": args.wall_seconds,
                        "once": args.once, "drain": args.drain, "allocator": args.allocator}
            if args.command != "run":
                worker(args.database, policy, args.command, settings)
            else:
                context = multiprocessing.get_context("spawn")
                settings["stop"] = context.Event()
                workers = [context.Process(target=worker, args=(args.database, policy, role, settings))
                           for role, count in (("compute", args.compute_workers), ("certify", args.certify_workers))
                           for _ in range(count)]
                for process in workers:
                    process.start()
                while any(process.is_alive() for process in workers):
                    for process in workers:
                        process.join(timeout=0.1)
                    if any(process.exitcode not in (None, 0) for process in workers):
                        settings["stop"].set()
                if any(process.exitcode != 0 for process in workers):
                    print(json.dumps({"policy": queue.policy_key, "counts": queue.counts(), "worker_failure": True}))
                    return 1
        print(json.dumps({"policy": queue.policy_key, "counts": queue.counts()}, sort_keys=True))
    finally:
        queue.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
