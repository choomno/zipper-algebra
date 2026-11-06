"""Paired historical C++ epochs with wall time and native process CPU time.

This diagnostic validates each executable against an explicit archived source
snapshot and its build manifest. One independent uncached Python reference per
case is constructed before timing; every native result is compared in full with
that reference. Reference construction, per-trial comparison, startup and cold
cache clearing are separate measurements. Per-trial checked totals include only
the comparison, not the separately recorded reusable Python replay.
"""
from __future__ import annotations

import argparse
from contextlib import suppress
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import tarfile
import tempfile
import time
from typing import Any
import zlib

from zipper_algebra.zipping_algorithm import collected, memoized
from zipper_algebra.zipping_algorithm.catalog import Catalog
from zipper_algebra.zipping_algorithm.interpreter import Limits

ABI = "zipper-native-engine-v1"
PREFIX = "native/zipper_engine/"
CASES = {
    "b3_alternating_l10": (3, [1, -2] * 5),
    "b5_positive_l6": (5, [1, 2, 3, 4, 1, 2]),
    "b8_positive_l10": (8, [1, 2, 3, 4, 5, 6, 7, 1, 2, 3]),
    "smoke": (3, [1, 2]),
}


def encoded(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def arrays(value: Any) -> Any:
    if isinstance(value, (tuple, list)):
        return [arrays(item) for item in value]
    if isinstance(value, dict):
        return {key: arrays(item) for key, item in value.items()}
    return value


def equal(left: Any, right: Any) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(equal(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(equal(a, b) for a, b in zip(left, right))
    return bool(left == right)


def decimal_integer(value: Any) -> int:
    if type(value) is int:
        return value
    if not isinstance(value, str) or not value or value == "-":
        raise ValueError("invalid decimal coefficient")
    negative = value.startswith("-")
    digits = value[1:] if negative else value
    if any(c < "0" or c > "9" for c in digits):
        raise ValueError("invalid decimal coefficient")
    result = 0
    for start in range(0, len(digits), 9):
        piece = digits[start:start + 9]
        result = result * 10**len(piece) + int(piece)
    return -result if negative else result


def decode(value: Any) -> Any:
    if isinstance(value, list):
        return [decode(item) for item in value]
    if isinstance(value, dict):
        return {key: ([[q, decimal_integer(c)] for q, c in item] if key == "coefficient"
                      else decode(item)) for key, item in value.items()}
    return value


def validate_epoch(specification: dict[str, str]) -> dict[str, Any]:
    """Bind every historical translation unit/header to the retained archive."""
    if set(specification) != {"label", "binary", "source_identity"} or not specification["label"]:
        raise ValueError("epoch needs label, binary and source_identity")
    identity_path = Path(specification["source_identity"]).resolve()
    binary = Path(specification["binary"]).resolve()
    identity = json.loads(identity_path.read_text())
    archive = Path(identity["archive"])
    if not archive.is_absolute():
        archive = identity_path.parent / archive
    if sha(archive) != identity["archive_sha256"]:
        raise ValueError("source archive hash differs")
    record = json.loads(binary.with_suffix(".json").read_text())
    if record["abi"] != ABI or sha(binary) != record["binary_sha256"]:
        raise ValueError("executable ABI or binary hash differs")
    native = identity["sources"]["native"]
    expected = {name[len(PREFIX):]: value for name, value in native.items()
                if name.startswith(PREFIX)}
    if len(expected) != len(native) or not expected or expected != record["source_sha256"]:
        raise ValueError("archive source identity differs from executable build")
    found = {}
    with tarfile.open(archive, "r:gz") as saved:
        for member in saved.getmembers():
            if (member.name.startswith(PREFIX)
                    and Path(member.name).suffix in (".cpp", ".hpp", ".h")):
                name = member.name[len(PREFIX):]
                if not member.isfile() or name in found:
                    raise ValueError("nonregular or duplicate archived source")
                payload = saved.extractfile(member)
                if payload is None:
                    raise ValueError("archived source is unreadable")
                found[name] = hashlib.sha256(payload.read()).hexdigest()
    if found != expected:
        raise ValueError("archived source bytes differ")
    if identity["sources"]["python"] != memoized._source():
        raise ValueError("historical epoch requires another Python replay source")
    return {"specification": specification, "source_identity_sha256": sha(identity_path),
            "archive_sha256": identity["archive_sha256"], "archive": str(archive),
            "build": record, "build_manifest_sha256": sha(binary.with_suffix(".json"))}


def cpu_snapshot(pid: int) -> tuple[int, int]:
    # comm may contain spaces and parentheses. The fields after its final ')'
    # begin at Linux stat field 3 (state), so utime/stime are offsets 11/12.
    fields = Path(f"/proc/{pid}/stat").read_text().rpartition(")")[2].split()
    return int(fields[11]) + int(fields[12]), int(fields[19])


def memory_snapshot(pid: int) -> dict[str, int]:
    result = {}
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        key = line.partition(":")[0]
        if key in ("VmRSS", "VmHWM", "VmSize"):
            result[key + "_kib"] = int(line.split()[1])
    return result


class EpochSession:
    """A historical executable with archive validation, without source monkeypatches."""

    def __init__(self, specification: dict[str, str], *, cpu: int | None = None) -> None:
        self.identity = validate_epoch(specification)
        self._error = tempfile.TemporaryFile(mode="w+t")
        self.process: subprocess.Popen[str] | None = None
        start = time.perf_counter()
        try:
            self.process = subprocess.Popen([specification["binary"]], stdin=subprocess.PIPE,
                                            stdout=subprocess.PIPE, stderr=self._error,
                                            text=True, encoding="utf-8", bufsize=1)
            if cpu is not None:
                os.sched_setaffinity(self.process.pid, {cpu})
            self.affinity = sorted(os.sched_getaffinity(self.process.pid))
            self.catalog_sha = Catalog.bundled().sha256
            payload = files("zipper_algebra.zipping_algorithm").joinpath("data/full-catalog.json").read_bytes()
            if hashlib.sha256(payload).hexdigest() != self.catalog_sha:
                raise ValueError("catalog changed during initialization")
            ready = self.request({"op": "initialize", "abi": ABI, "catalog": json.loads(payload),
                                  "catalog_sha256": self.catalog_sha, "reference_sources": collected._source(),
                                  "cache_entries": 262144, "cache_bytes": 256 * 1024**2})
            if ready != {"abi": ABI, "ready": True}:
                raise ValueError("unexpected historical initialization response")
        except BaseException:
            self.close()
            raise
        self.startup_seconds = time.perf_counter() - start

    @property
    def pid(self) -> int:
        if self.process is None:
            raise RuntimeError("session closed")
        return self.process.pid

    def request(self, request: dict[str, Any]) -> dict[str, Any]:
        process = self.process
        if process is None or process.stdin is None or process.stdout is None:
            raise RuntimeError("session closed")
        try:
            process.stdin.write(encoded(request).decode() + "\n")
            process.stdin.flush()
            line = process.stdout.readline()
            if not line:
                raise RuntimeError("historical engine stopped without response")
            value = json.loads(line)
            if not isinstance(value, dict) or "protocol_error" in value:
                raise RuntimeError("historical engine rejected request: " + repr(value))
            return value
        except BaseException:
            self.close()
            raise

    def reduce(self, n: int, word: list[int], limits: Limits, batch_size: int,
               trace: bool) -> tuple[dict[str, Any], dict[str, Any]]:
        body = decode(self.request({"op": "reduce", "n": n, "word": word,
                                    "limits": asdict(limits), "batch_size": batch_size, "trace": trace}))
        runtime = body.pop("native_runtime")
        body.update(n=n, input_word=list(word), coefficient_ring="Z[q,q^-1]",
                    trace_enabled=trace, limits=asdict(limits),
                    schema="zipper-collected-run-v1", scope="literal_word_collected_spanning_procedure",
                    source_sha256=collected._source(), catalog_sha256=self.catalog_sha)
        return body, runtime

    def close(self) -> None:
        process, self.process = self.process, None
        if process is not None:
            with suppress(OSError):
                process.kill()
            with suppress(OSError, subprocess.TimeoutExpired):
                process.wait(timeout=5)
            for stream in (process.stdin, process.stdout):
                if stream is not None:
                    with suppress(OSError, ValueError):
                        stream.close()
        self._error.close()


def initialize_database(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute("PRAGMA journal_mode=WAL")
    connection.executescript('''
        CREATE TABLE IF NOT EXISTS native_epoch_campaigns(
            id INTEGER PRIMARY KEY, started_utc TEXT NOT NULL, config_json TEXT NOT NULL,
            sources_json TEXT NOT NULL, setup_json TEXT NOT NULL, status TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS native_epoch_references(
            campaign_id INTEGER NOT NULL REFERENCES native_epoch_campaigns(id), case_name TEXT NOT NULL,
            n INTEGER NOT NULL, word_json TEXT NOT NULL, replay_seconds REAL NOT NULL,
            certificate_sha256 TEXT NOT NULL, certificate_zlib BLOB NOT NULL,
            PRIMARY KEY(campaign_id,case_name));
        CREATE TABLE IF NOT EXISTS native_epoch_measurements(
            campaign_id INTEGER NOT NULL REFERENCES native_epoch_campaigns(id),
            round INTEGER NOT NULL, position INTEGER NOT NULL, epoch TEXT NOT NULL,
            case_name TEXT NOT NULL, producer_wall_seconds REAL NOT NULL,
            native_cpu_ticks INTEGER NOT NULL, native_cpu_seconds REAL NOT NULL,
            comparison_seconds REAL NOT NULL, producer_plus_comparison_seconds REAL NOT NULL,
            certificate_sha256 TEXT NOT NULL, certificate_zlib BLOB NOT NULL,
            telemetry_json TEXT NOT NULL, verified INTEGER NOT NULL CHECK(verified=1),
            PRIMARY KEY(campaign_id,round,epoch,case_name));
    ''')
    connection.commit()
    return connection


def run(database: Path, specifications: list[dict[str, str]], *, cases: list[str], rounds: int,
        limits: Limits, batch_size: int = 64, trace: bool = False, cpu: int | None = None) -> int:
    labels = [spec["label"] for spec in specifications]
    if (not labels or len(set(labels)) != len(labels) or type(rounds) is not int or rounds < 1
            or not cases or len(set(cases)) != len(cases) or any(c not in CASES for c in cases)):
        raise ValueError("distinct epochs, distinct known cases, and positive rounds required")
    if type(batch_size) is not int or batch_size < 1 or type(trace) is not bool:
        raise ValueError("positive batch size and boolean trace required")
    if cpu is not None and (type(cpu) is not int or cpu not in os.sched_getaffinity(0)):
        raise ValueError("CPU must be an allowed integer CPU id")
    source = {"harness": sha(Path(__file__)), "python": memoized._source()}
    ticks = os.sysconf("SC_CLK_TCK")
    config = {"epochs": labels, "cases": {name: CASES[name] for name in cases},
              "limits": asdict(limits), "batch_size": batch_size, "trace": trace, "rounds": rounds,
              "cache_policy": "cold_each_input", "cache_entries": 262144, "cache_bytes": 256 * 1024**2,
              "ordering": "cyclic rotations then reversed rotations; adjacent epochs per input",
              "native_cpu_clock_ticks_per_second": ticks,
              "native_child_cpu": cpu, "parent_allowed_cpus": sorted(os.sched_getaffinity(0)),
              "affinity_scope": "only benchmark-owned native children; census unchanged",
              "verification": "one uncached Python replay per case; strict full comparison and hash each trial",
              "checked_total_scope": "producer plus comparison only; reusable reference construction separate",
              "producer_scope": "request encoding, IPC, native computation, response decode and reference projection"}
    connection = initialize_database(database)
    sessions: dict[str, EpochSession] = {}
    campaign = None
    try:
        for spec in specifications:
            sessions[spec["label"]] = EpochSession(spec, cpu=cpu)
        setup = {name: {"identity": session.identity, "startup_seconds": session.startup_seconds,
                        "native_affinity": session.affinity}
                 for name, session in sessions.items()}
        cur = connection.execute("INSERT INTO native_epoch_campaigns(started_utc,config_json,sources_json,setup_json,status) VALUES(?,?,?,?,?)",
                                 (datetime.now(timezone.utc).isoformat(), encoded(config).decode(), encoded(source).decode(),
                                  encoded(setup).decode(), "running"))
        campaign = int(cur.lastrowid)
        connection.commit()
        references = {}
        for name in cases:
            n, word = CASES[name]
            started = time.perf_counter()
            oracle = collected.CollectedSession()
            oracle._catalog.configure_rule_cache(max_bytes=0, max_entries=0)
            expected = arrays(oracle.reduce(word, n=n, limits=limits, trace=trace, batch_size=batch_size))
            duration = time.perf_counter() - started
            payload = encoded(expected)
            references[name] = expected, hashlib.sha256(payload).hexdigest()
            connection.execute("INSERT INTO native_epoch_references VALUES(?,?,?,?,?,?,?)",
                               (campaign, name, n, encoded(word).decode(), duration,
                                references[name][1], zlib.compress(payload)))
            connection.commit()
        orders = [labels[i:] + labels[:i] for i in range(len(labels))]
        if len(labels) > 2:
            reverse = list(reversed(labels))
            orders += [reverse[i:] + reverse[:i] for i in range(len(labels))]
        for round_index in range(rounds):
            for name in cases:
                n, word = CASES[name]
                for position, label in enumerate(orders[round_index % len(orders)]):
                    session = sessions[label]
                    if validate_epoch(session.identity["specification"]) != session.identity:
                        raise RuntimeError("historical executable identity changed")
                    if {"harness": sha(Path(__file__)), "python": memoized._source()} != source:
                        raise RuntimeError("benchmark or Python replay source changed")
                    started = time.perf_counter()
                    session.request({"op": "clear_cache"})
                    clear_seconds = time.perf_counter() - started
                    before_ticks, process_start = cpu_snapshot(session.pid)
                    started = time.perf_counter()
                    actual, runtime = session.reduce(n, word, limits, batch_size, trace)
                    producer = time.perf_counter() - started
                    after_ticks, following_start = cpu_snapshot(session.pid)
                    if process_start != following_start or after_ticks < before_ticks:
                        raise RuntimeError("native process CPU clock identity changed")
                    started = time.perf_counter()
                    expected, reference_sha = references[name]
                    if not equal(actual, expected):
                        raise AssertionError(f"exact reference differs: {label}/{name}/{round_index}")
                    payload = encoded(actual)
                    actual_sha = hashlib.sha256(payload).hexdigest()
                    if actual_sha != reference_sha:
                        raise AssertionError("equal certificate serialization hash differs")
                    comparison = time.perf_counter() - started
                    telemetry = {"native_runtime": runtime, "native_process": memory_snapshot(session.pid),
                                 "cache_clear_seconds": clear_seconds, "certificate_bytes": len(payload),
                                 "native_process_start_ticks": process_start}
                    if validate_epoch(session.identity["specification"]) != session.identity:
                        raise RuntimeError("historical executable identity changed during measurement")
                    if {"harness": sha(Path(__file__)), "python": memoized._source()} != source:
                        raise RuntimeError("benchmark or Python replay source changed during measurement")
                    connection.execute("INSERT INTO native_epoch_measurements VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,1)",
                                       (campaign, round_index, position, label, name, producer, after_ticks-before_ticks,
                                        (after_ticks-before_ticks)/ticks, comparison, producer+comparison,
                                        actual_sha, zlib.compress(payload), encoded(telemetry).decode()))
                    connection.commit()
                    print(json.dumps({"epoch_campaign": campaign, "round": round_index, "epoch": label, "case": name,
                                      "producer_seconds": producer, "native_cpu_seconds": (after_ticks-before_ticks)/ticks,
                                      "comparison_seconds": comparison}), flush=True)
        connection.execute("UPDATE native_epoch_campaigns SET status='complete' WHERE id=?", (campaign,))
        connection.commit()
        return campaign
    except BaseException as error:
        if campaign is not None:
            connection.execute("UPDATE native_epoch_campaigns SET status=? WHERE id=?",
                               ("interrupted" if isinstance(error, KeyboardInterrupt) else "failed", campaign))
            connection.commit()
        raise
    finally:
        for session in sessions.values():
            session.close()
        connection.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--epochs", type=Path, required=True, help="JSON list of label,binary,source_identity objects")
    parser.add_argument("--cases", nargs="+", choices=CASES, default=list(CASES)[:3])
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--max-steps", type=int, default=4096)
    parser.add_argument("--max-terms", type=int, default=10000)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--trace", action="store_true")
    parser.add_argument("--cpu", type=int, help="Pin only the benchmark-owned native processes to this allowed CPU")
    args = parser.parse_args()
    run(args.database, json.loads(args.epochs.read_text()), cases=args.cases, rounds=args.rounds,
        limits=Limits(max_steps=args.max_steps, max_terms=args.max_terms, expansion_budget=1000000),
        batch_size=args.batch_size, trace=args.trace, cpu=args.cpu)


if __name__ == "__main__":
    main()
