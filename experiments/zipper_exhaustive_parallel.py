"""Parallel literal-word enumeration with one ordered SQLite writer.

Workers reduce and JSON-replay-check independent literal inputs. The writer
commits results in the original length/rank/ordinal order, so interruption,
worker failure, or an incomplete word leaves an exact consecutive prefix.
At most twice the worker count is in flight; later speculative results are
never committed past the first unresolved word. Existing sequential-driver
rows are preserved and fully replay-audited before each invocation.

The database retains the algorithm/catalog identity used by exhaust_general.
This runner's SHA-256 is additionally bound in metadata and every new config;
changing it requires a separate campaign rather than silently reusing results.
Only the parent opens SQLite; workers use spawn. WAL allows progress queries
during execution. A Linux flock rejects another parallel writer for this
campaign. The sequential driver does not use this lock and must not run
concurrently on the same database.

Example (from the repository with its package installed):
  python experiments/zipper_exhaustive_parallel.py results --ranks 2 3 4 5 6 7 8 \
      --workers 4 --wall-seconds 900 --max-words 8463
"""
from __future__ import annotations

import argparse
import fcntl
import multiprocessing
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time

from zipper_algebra.zipping_algorithm.algorithm import AlgorithmSession
from zipper_algebra.zipping_algorithm.enumeration import _advance, _integer, _json, _sha, _word_at
from zipper_algebra.zipping_algorithm.general_enumeration import _audit, _identity, _layers
from zipper_algebra.zipping_algorithm.interpreter import Limits

_session = None
_worker_limits = None
_worker_trace = False
_worker_identity = None


def _runner_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _assert_source(identity, runner_hash):
    if _identity(identity["ranks"], identity["trace"]) != identity:
        raise ValueError("enumeration source or catalog changed during the calculation")
    if _runner_hash() != runner_hash:
        raise ValueError("parallel runner source changed during the calculation")


def _initialize_worker(limits, trace, identity):
    global _session, _worker_limits, _worker_trace, _worker_identity
    _session = AlgorithmSession()
    _worker_limits = Limits(**limits)
    _worker_trace = trace
    _worker_identity = identity


def _compute(n, word):
    produced = time.perf_counter_ns()
    certificate = _session.reduce(word, n=n, limits=_worker_limits, trace=_worker_trace)
    producer_ns = time.perf_counter_ns() - produced
    encoded = _json(certificate)
    checked_at = time.perf_counter_ns()
    if _session.check(json.loads(encoded)) is not True:
        raise ValueError("produced word certificate failed exact checking")
    checker_ns = time.perf_counter_ns() - checked_at
    if (certificate["source_sha256"] != _worker_identity["algorithm_source_sha256"]
            or certificate["catalog_sha256"] != _worker_identity["catalog_sha256"]):
        raise ValueError("algorithm source or catalog changed during the calculation")
    return encoded, producer_ns, checker_ns


def _process_pool(**kwargs):
    return ProcessPoolExecutor(mp_context=multiprocessing.get_context("spawn"), **kwargs)


def exhaust_parallel(output, *, ranks=tuple(range(2, 9)), workers=4,
                     wall_seconds=900, max_length=None, max_words=None,
                     limits=None, trace=False, _executor_factory=_process_pool):
    """Extend one checked prefix; max_words caps total rows, not new rows.

    Wall time starts after the complete saved-prefix audit and is soft: a word
    already being awaited can finish after the deadline. Shutdown also waits
    for the bounded number of running speculative words, which are discarded
    if they lie beyond the committed prefix. No input is removed by isotopy.
    """
    if (not isinstance(ranks, (tuple, list)) or not ranks
            or len(set(ranks)) != len(ranks)):
        raise ValueError("distinct positive ranks required")
    for rank in ranks:
        _integer(rank, "rank", 1)
    _integer(workers, "workers", 1)
    if type(wall_seconds) not in (int, float) or not math.isfinite(wall_seconds) or wall_seconds <= 0:
        raise ValueError("wall_seconds must be finite and positive")
    if max_length is not None:
        _integer(max_length, "max_length", 0)
    if max_words is not None:
        _integer(max_words, "max_words", 0)
    if type(trace) is not bool:
        raise ValueError("trace must be boolean")
    limits = Limits(max_steps=200000, max_terms=100000, expansion_budget=1000000) if limits is None else limits
    if not isinstance(limits, Limits):
        raise ValueError("limits must be a Limits record")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    database = output / "words.sqlite3"
    ranks = list(ranks)
    identity = _identity(ranks, trace)
    runner_hash = _runner_hash()
    config = {"wall_seconds": wall_seconds, "max_length": max_length,
              "max_words": max_words, "workers": workers, "max_in_flight": 2 * workers,
              "limits": asdict(limits), "identity_sha256": _sha(_json(identity)),
              "runner_sha256": runner_hash}
    config_sha = _sha(_json(config))
    lock = (output / ".parallel.lock").open("a+")
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        lock.close()
        raise ValueError("another parallel runner owns this campaign") from error
    try:
        connection = sqlite3.connect(database)
    except BaseException:
        lock.close()
        raise
    executor = None
    pending = deque()
    try:
        connection.execute("PRAGMA journal_mode=WAL")
        connection.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY,value TEXT NOT NULL)")
        connection.execute("CREATE TABLE IF NOT EXISTS words (n INTEGER,length INTEGER,ordinal INTEGER,rank_order INTEGER,word TEXT,input_sha TEXT,status,diagnostic TEXT,certificate TEXT,certificate_sha TEXT,metrics TEXT,producer_ns INTEGER,checker_ns INTEGER,config_sha TEXT, PRIMARY KEY(n,length,ordinal))")
        connection.execute("CREATE TABLE IF NOT EXISTS configurations (hash TEXT PRIMARY KEY,config TEXT NOT NULL)")
        metadata = dict(connection.execute("SELECT key,value FROM metadata"))
        if metadata:
            if metadata.get("identity") != _json(identity):
                raise ValueError("enumeration source, catalog, rank order, or trace configuration changed")
            previous_runner = metadata.get("parallel_runner_sha256")
            if previous_runner is not None and previous_runner != runner_hash:
                raise ValueError("parallel runner source changed")
            cursor = json.loads(metadata["cursor"])
        else:
            cursor = [0, 0, 0]
            with connection:
                connection.executemany("INSERT INTO metadata VALUES (?,?)",
                                       [("identity", _json(identity)), ("cursor", _json(cursor))])
        audit_started = time.perf_counter()
        audited = _audit(connection, ranks, cursor, identity, AlgorithmSession())
        audit_seconds = time.perf_counter() - audit_started
        with connection:
            connection.execute("INSERT OR IGNORE INTO configurations VALUES (?,?)", (config_sha, _json(config)))
            connection.execute("INSERT OR REPLACE INTO metadata VALUES ('parallel_runner_sha256',?)", (runner_hash,))
        previous = connection.execute("SELECT certificate FROM words WHERE n=? AND length=? AND ordinal=?",
                                      (ranks[cursor[1]], cursor[0], cursor[2])).fetchone()
        reason = None
        if previous:
            checkpoint = json.loads(previous[0])
            old, new = checkpoint["limits"], asdict(limits)
            if checkpoint["status"] == "control_obligation":
                reason = "control_obligation"
            elif any(new[key] < old[key] for key in old) or new == old:
                reason = "word_limit_needs_larger_limits"
        start = time.perf_counter()
        new_words = 0
        total_rows = audited
        peak_in_flight = 0
        next_to_submit = cursor[:]
        if reason is None:
            executor = _executor_factory(max_workers=workers, initializer=_initialize_worker,
                                         initargs=(asdict(limits), trace, identity))
        while reason is None:
            if ranks == [1] and cursor[0] > 0:
                reason = "all_nonempty_layers_exhausted"
                break
            if max_length is not None and cursor[0] > max_length:
                reason = "max_length"
                break
            if max_words is not None and total_rows >= max_words and not previous:
                reason = "max_words"
                break
            if time.perf_counter() - start >= wall_seconds:
                reason = "wall_budget_between_words"
                break
            while len(pending) < 2 * workers:
                length, rank_index, ordinal = next_to_submit
                if ((ranks == [1] and length > 0)
                        or (max_length is not None and length > max_length)
                        or (max_words is not None and total_rows + len(pending) >= max_words
                            and not previous)):
                    break
                n = ranks[rank_index]
                word = list(_word_at(n, length, ordinal))
                pending.append((next_to_submit[:], n, word, executor.submit(_compute, n, word)))
                next_to_submit = _advance(next_to_submit, ranks)
                peak_in_flight = max(peak_in_flight, len(pending))
            if not pending:
                raise RuntimeError("no job available before the requested stopping boundary")
            committed_cursor, n, word, future = pending.popleft()
            if committed_cursor != cursor:
                raise ValueError("parallel result is not the next literal input")
            encoded, producer_ns, checker_ns = future.result()
            certificate = json.loads(encoded)
            if (certificate["n"] != n or certificate["input_word"] != word
                    or certificate["source_sha256"] != identity["algorithm_source_sha256"]
                    or certificate["catalog_sha256"] != identity["catalog_sha256"]
                    or certificate["limits"] != asdict(limits)
                    or certificate["trace_enabled"] != trace):
                raise ValueError("worker result identity does not match its literal input")
            _assert_source(identity, runner_hash)
            status = certificate["status"]
            if status not in ("complete", "incomplete", "control_obligation"):
                raise ValueError("unknown worker status")
            next_cursor = _advance(cursor, ranks) if status == "complete" else cursor
            length, rank_index, ordinal = cursor
            with connection:
                connection.execute("INSERT OR REPLACE INTO words VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (n,length,ordinal,rank_index,_json(word),_sha(_json({"n":n,"word":word})),status,
                     _json(certificate["diagnostic"]),encoded,_sha(encoded),_json(certificate["metrics"]),
                     producer_ns,checker_ns,config_sha))
                connection.execute("UPDATE metadata SET value=? WHERE key='cursor'", (_json(next_cursor),))
                connection.execute("INSERT OR REPLACE INTO metadata VALUES ('last_config',?)", (_json(config),))
            total_rows += not bool(previous)
            previous = None
            cursor = next_cursor
            new_words += 1
            if status != "complete":
                reason = "control_obligation" if status == "control_obligation" else "word_limit"
        for _, _, _, future in pending:
            future.cancel()
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
            executor = None
        _assert_source(identity, runner_hash)
        summary = {"schema": identity["schema"], "identity": identity, "config": config,
                   "config_sha256": config_sha, "database": str(database.resolve()),
                   "runner_sha256": runner_hash, "audited_words": audited,
                   "new_or_resumed_words": new_words, "audit_seconds": audit_seconds,
                   "new_work_seconds": time.perf_counter() - start, "stop_reason": reason,
                   "peak_in_flight": peak_in_flight,
                   "next_cursor": {"length":cursor[0],"rank_index":cursor[1],"ordinal":cursor[2]},
                   "layers": _layers(connection)}
        temporary = output / "summary.json.tmp"
        temporary.write_text(_json(summary) + "\n")
        temporary.replace(output / "summary.json")
        return summary
    finally:
        for _, _, _, future in pending:
            future.cancel()
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
        connection.close()
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--ranks", nargs="+", type=int, default=list(range(2, 9)))
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--wall-seconds", type=float, default=900)
    parser.add_argument("--max-length", type=int)
    parser.add_argument("--max-words", type=int)
    parser.add_argument("--limits", type=json.loads)
    parser.add_argument("--trace", action="store_true")
    args = parser.parse_args()
    options = vars(args)
    options["limits"] = Limits(**args.limits) if args.limits is not None else None
    result = exhaust_parallel(**options)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
