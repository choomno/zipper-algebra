"""Durable provisional computations and independent certification in a sidecar.

No candidate is reusable as certified evidence. Native process-local caches may
share computations, but certification replays from the literal input. Existing
synchronous and historical census databases are never opened for writing.
"""
from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import contextmanager, closing
from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import threading
import time
from typing import Any
import uuid
import zlib

from .census_reuse import _json, _sha, _canonical_digest, reconstruction_sources, verifier_binding, _IMPORT_SOURCES
from .native import _decode_coefficient_object
from .native_verifier import VerifierIncomplete
from .normal_form_reconstruction import completed_output, _identity

SCHEMA = "zipper-async-census-v1"
_BACKLOG = "('computing','candidate','certifying','verification_incomplete','verification_error')"


@dataclass(frozen=True)
class Claim:
    key: str
    n: int
    word: list[int]
    token: str
    stage: str


class SourceChanged(ValueError):
    """Source drift is operational, not evidence of a mathematical mismatch."""


class LostLease(RuntimeError):
    """A superseded worker cannot commit its result."""


class CensusQueue:
    """Concurrent SQLite queue; instantiate one connection per worker/thread.

    Policy namespaces are immutable. Every retry and certification records an
    append-only attempt. Lease tokens prevent a reclaimed worker from committing.
    """

    def __init__(self, path: str | Path, policy: dict[str, Any]) -> None:
        self.path = Path(path)
        self.policy = deepcopy(policy)
        if (policy.get("schema") != SCHEMA or not policy.get("producer")
                or not policy.get("producer_metadata") or not policy.get("limits")
                or not policy.get("package_sources")):
            raise ValueError("explicit async policy and source identities required")
        if policy["package_sources"] != _IMPORT_SOURCES:
            raise ValueError("queue package sources differ from imported code")
        self.policy_key = _sha(policy)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.path.exists():
            with closing(sqlite3.connect(self.path.resolve().as_uri() + "?mode=ro", uri=True)) as probe:
                names = {r[0] for r in probe.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                if names and "async_policies" not in names:
                    raise ValueError("async census requires a separate sidecar database")
        self.db = sqlite3.connect(self.path, timeout=60, isolation_level=None)
        self.db.row_factory = sqlite3.Row
        self.db.execute("PRAGMA foreign_keys=ON")
        self.db.execute("PRAGMA journal_mode=WAL")
        self.db.executescript("""
            CREATE TABLE IF NOT EXISTS async_policies(key TEXT PRIMARY KEY,payload TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS async_jobs(
              policy TEXT NOT NULL REFERENCES async_policies(key),key TEXT NOT NULL,
              n INTEGER NOT NULL,word TEXT NOT NULL,state TEXT NOT NULL,
              token TEXT,lease REAL,created REAL NOT NULL,updated REAL NOT NULL,
              PRIMARY KEY(policy,key));
            CREATE INDEX IF NOT EXISTS async_work ON async_jobs(policy,state,created);
            CREATE TABLE IF NOT EXISTS async_candidates(
              policy TEXT NOT NULL,key TEXT NOT NULL,payload BLOB NOT NULL,digest TEXT NOT NULL,
              output_digest TEXT NOT NULL,bytes INTEGER NOT NULL,
              PRIMARY KEY(policy,key),FOREIGN KEY(policy,key) REFERENCES async_jobs(policy,key));
            CREATE TABLE IF NOT EXISTS async_attempts(
              id INTEGER PRIMARY KEY,policy TEXT NOT NULL,key TEXT NOT NULL,stage TEXT NOT NULL,
              status TEXT NOT NULL,record TEXT NOT NULL,created REAL NOT NULL,
              FOREIGN KEY(policy,key) REFERENCES async_jobs(policy,key));
            CREATE TABLE IF NOT EXISTS async_cursors(
              policy TEXT NOT NULL,layer TEXT NOT NULL,ordinal TEXT NOT NULL,
              PRIMARY KEY(policy,layer));
        """)
        self.db.execute("INSERT OR IGNORE INTO async_policies VALUES (?,?)", (self.policy_key, _json(policy)))
        stored = self.db.execute("SELECT payload FROM async_policies WHERE key=?", (self.policy_key,)).fetchone()[0]
        if stored != _json(policy):
            self.close()
            raise ValueError("policy binding mismatch")

    @contextmanager
    def transaction(self) -> Iterator[None]:
        self.db.execute("BEGIN IMMEDIATE")
        try:
            yield
            self.db.execute("COMMIT")
        except BaseException:
            self.db.execute("ROLLBACK")
            raise

    def close(self) -> None:
        self.db.close()

    def assert_sources(self) -> None:
        if (_IMPORT_SOURCES != self.policy["package_sources"]
                or reconstruction_sources() != self.policy["package_sources"]):
            raise SourceChanged("queue package source changed; use the recorded source checkout")

    def _enqueue(self, n: int, word: list[int]) -> str:
        if (type(n) is not int or n < 1 or type(word) is not list
                or any(type(a) is not int or not 1 <= abs(a) < n for a in word)):
            raise ValueError("invalid literal braid request")
        key = _sha({"n": n, "word": word})
        now = time.time()
        self.db.execute("INSERT OR IGNORE INTO async_jobs VALUES (?,?,?,?,?,?,?,?,?)",
                        (self.policy_key, key, n, _json(word), "queued", None, None, now, now))
        return key

    def enqueue(self, n: int, word: list[int]) -> str:
        self.assert_sources()
        with self.transaction():
            return self._enqueue(n, word)

    def enqueue_interval(self, layer: str, count: int, word_at: Callable[[int], tuple[int, list[int]]],
                         *, batch: int = 256, rotation: tuple[str, str] | None = None) -> int:
        """Atomically admit a bounded ordinal interval and advance its cursor."""
        if type(count) is not int or count < 0 or type(batch) is not int or batch < 1:
            raise ValueError("nonnegative count and positive batch required")
        self.assert_sources()
        with self.transaction():
            row = self.db.execute("SELECT ordinal FROM async_cursors WHERE policy=? AND layer=?",
                                  (self.policy_key, layer)).fetchone()
            start = int(row[0]) if row else 0
            stop = min(count, start + batch)
            for ordinal in range(start, stop):
                self._enqueue(*word_at(ordinal))
            self.db.execute("INSERT OR REPLACE INTO async_cursors VALUES (?,?,?)",
                            (self.policy_key, layer, str(stop)))
            if rotation is not None:
                self.db.execute("INSERT OR REPLACE INTO async_cursors VALUES (?,?,?)",
                                (self.policy_key, rotation[0], rotation[1]))
        return stop - start

    def claim(self, stage: str, *, lease_seconds: float = 300,
              max_pending: int = 256, max_pending_bytes: int = 2**30) -> Claim | None:
        if stage not in ("compute", "certify") or not math.isfinite(lease_seconds) or lease_seconds <= 0:
            raise ValueError("valid stage and positive lease required")
        if (type(max_pending) is not int or type(max_pending_bytes) is not int
                or max_pending < 1 or max_pending_bytes < 1):
            raise ValueError("positive backlog limits required")
        active, ready = ("computing", "queued") if stage == "compute" else ("certifying", "candidate")
        with self.transaction():
            now = time.time()
            # Recovery is conditional on expiry; live workers heartbeat separately.
            expired = self.db.execute("SELECT key FROM async_jobs WHERE policy=? AND state=? AND lease<?",
                                      (self.policy_key, active, now)).fetchall()
            for row in expired:
                self._attempt(row[0], stage, "lease_expired", {})
            self.db.execute("UPDATE async_jobs SET state=?,token=NULL,lease=NULL,updated=? "
                            "WHERE policy=? AND state=? AND lease<?", (ready, now, self.policy_key, active, now))
            if stage == "compute":
                pending = self.db.execute("SELECT count(*) FROM async_jobs WHERE policy=? AND state IN " + _BACKLOG,
                                          (self.policy_key,)).fetchone()[0]
                size = self.db.execute("SELECT coalesce(sum(c.bytes),0) FROM async_candidates c JOIN async_jobs j "
                                       "ON c.policy=j.policy AND c.key=j.key WHERE j.policy=? AND j.state IN " + _BACKLOG,
                                       (self.policy_key,)).fetchone()[0]
                if pending >= max_pending or size >= max_pending_bytes:
                    return None
            row = self.db.execute("SELECT key,n,word FROM async_jobs WHERE policy=? AND state=? "
                                  "ORDER BY created,rowid LIMIT 1", (self.policy_key, ready)).fetchone()
            if row is None:
                return None
            token = uuid.uuid4().hex
            self.db.execute("UPDATE async_jobs SET state=?,token=?,lease=?,updated=? WHERE policy=? AND key=?",
                            (active, token, now + lease_seconds, now, self.policy_key, row[0]))
            return Claim(row[0], row[1], json.loads(row[2]), token, stage)

    def heartbeat(self, claim: Claim, lease_seconds: float) -> bool:
        return self.db.execute("UPDATE async_jobs SET lease=? WHERE policy=? AND key=? AND token=?",
                               (time.time() + lease_seconds, self.policy_key, claim.key, claim.token)).rowcount == 1

    @contextmanager
    def keepalive(self, claim: Claim, lease_seconds: float = 300) -> Iterator[None]:
        stop = threading.Event()
        def beat() -> None:
            try:
                with closing(sqlite3.connect(self.path, timeout=60)) as db:
                    while not stop.wait(lease_seconds / 3):
                        changed = db.execute("UPDATE async_jobs SET lease=? WHERE policy=? AND key=? AND token=?",
                                             (time.time() + lease_seconds, self.policy_key, claim.key, claim.token)).rowcount
                        db.commit()
                        if changed != 1:
                            return
            except Exception:
                # Commit-time token comparison remains authoritative. A dead
                # heartbeat permits recovery; it never permits a stale commit.
                return
        thread = threading.Thread(target=beat, daemon=True)
        thread.start()
        try:
            yield
        finally:
            stop.set()
            thread.join()

    def _owned(self, claim: Claim) -> None:
        row = self.db.execute("SELECT token,state,n,word FROM async_jobs WHERE policy=? AND key=?",
                              (self.policy_key, claim.key)).fetchone()
        state = "computing" if claim.stage == "compute" else "certifying"
        if row is None or tuple(row) != (claim.token, state, claim.n, _json(claim.word)):
            raise LostLease("claim was superseded or altered")

    def _attempt(self, key: str, stage: str, status: str, record: dict[str, Any]) -> None:
        self.db.execute("INSERT INTO async_attempts(policy,key,stage,status,record,created) VALUES (?,?,?,?,?,?)",
                        (self.policy_key, key, stage, status, _json(record), time.time()))

    def finish(self, claim: Claim, status: str, record: dict[str, Any]) -> None:
        allowed = ({"compute_incomplete", "compute_error"} if claim.stage == "compute" else
                   {"certified", "rejected", "verification_incomplete", "verification_error"})
        if status not in allowed:
            raise ValueError("invalid stage outcome")
        with self.transaction():
            self._owned(claim)
            self._attempt(claim.key, claim.stage, status, record)
            self.db.execute("UPDATE async_jobs SET state=?,token=NULL,lease=NULL,updated=? WHERE policy=? AND key=?",
                            (status, time.time(), self.policy_key, claim.key))

    def _validate(self, result: dict[str, Any], claim: Claim) -> str:
        if (result.get("n") != claim.n or type(result.get("n")) is not int
                or _json(result.get("input_word")) != _json(claim.word)
                or _json(_identity({"schema": result.get("schema"),
                    "source_sha256": result.get("source_sha256"),
                    "catalog_sha256": result.get("catalog_sha256"),
                    "execution": result.get("execution")})) != _json(self.policy["producer"])
                or _json(result.get("limits")) != _json(self.policy["limits"])):
            raise ValueError("candidate request or source binding mismatch")
        if result.get("schema") in ("zipper-native-run-v1", "zipper-native-reuse-v1"):
            fields = ("abi", "binary_sha256", "compiler", "compiler_version", "flags",
                      "wrapper_source_sha256", "link_libraries")
            build = {key: self.policy["producer_metadata"].get(key) for key in fields}
            if any(value is None for value in build.values()) or _json(result.get("native_build")) != _json(build):
                raise ValueError("candidate native build binding mismatch")
        return _canonical_digest(completed_output(result))

    def publish(self, claim: Claim, result: dict[str, Any], *, seconds: float) -> None:
        if claim.stage != "compute":
            raise ValueError("only a compute claim can publish")
        self.assert_sources()
        output_digest = self._validate(result, claim)
        encoded = _json(result).encode()
        payload, digest = zlib.compress(encoded), hashlib.sha256(encoded).hexdigest()
        with self.transaction():
            self._owned(claim)
            self.db.execute("INSERT INTO async_candidates VALUES (?,?,?,?,?,?)",
                            (self.policy_key, claim.key, payload, digest, output_digest, len(payload)))
            self._attempt(claim.key, "compute", "candidate", {"seconds": seconds, "digest": digest,
                          "compressed_bytes": len(payload), "metrics": result.get("metrics", {})})
            self.db.execute("UPDATE async_jobs SET state='candidate',token=NULL,lease=NULL,updated=? "
                            "WHERE policy=? AND key=?", (time.time(), self.policy_key, claim.key))

    def candidate(self, claim: Claim) -> tuple[dict[str, Any], str]:
        self.assert_sources()
        self._owned(claim)
        row = self.db.execute("SELECT payload,digest,output_digest FROM async_candidates WHERE policy=? AND key=?",
                              (self.policy_key, claim.key)).fetchone()
        if row is None:
            raise ValueError("candidate missing")
        encoded = zlib.decompress(row[0])
        if hashlib.sha256(encoded).hexdigest() != row[1]:
            raise ValueError("candidate payload digest mismatch")
        result = json.loads(encoded, object_hook=_decode_coefficient_object)
        if self._validate(result, claim) != row[2]:
            raise ValueError("candidate output digest mismatch")
        return result, row[1]

    def compute_one(self, produce: Callable[[int, list[int]], dict[str, Any]], **claim_options: Any) -> bool:
        self.assert_sources()
        claim = self.claim("compute", **claim_options)
        if claim is None:
            return False
        started = time.perf_counter()
        with self.keepalive(claim, claim_options.get("lease_seconds", 300)):
            try:
                result = produce(claim.n, claim.word)
                if result.get("status") != "complete":
                    self.finish(claim, "compute_incomplete", {"seconds": time.perf_counter() - started,
                                "status": result.get("status"), "diagnostic": result.get("diagnostic")})
                else:
                    self.publish(claim, result, seconds=time.perf_counter() - started)
            except LostLease:
                raise
            except Exception as exc:
                self.finish(claim, "compute_error", {"seconds": time.perf_counter() - started, "error": str(exc)})
                raise
        return True

    def certify_one(self, verify: Callable[[dict[str, Any]], bool], verifier_identity: dict[str, Any],
                    **claim_options: Any) -> bool:
        verifier_identity = verifier_binding(verifier_identity)
        self.assert_sources()
        claim = self.claim("certify", **claim_options)
        if claim is None:
            return False
        record: dict[str, Any] = {"verifier": deepcopy(verifier_identity)}
        started = time.perf_counter()
        status = "certified"
        with self.keepalive(claim, claim_options.get("lease_seconds", 300)):
            try:
                result, record["digest"] = self.candidate(claim)
            except (LostLease, SourceChanged):
                raise
            except (ValueError, zlib.error) as exc:
                status, record["error"] = "rejected", str(exc)
                record["reason"] = "candidate_integrity"
            else:
                try:
                    if verify(result) is not True:
                        status = "rejected"
                        record["reason"] = "replay_mismatch"
                        record["error"] = "independent verifier rejected output"
                except VerifierIncomplete as exc:
                    status, record["error"] = "verification_incomplete", str(exc)
                    record["reason"] = "resource_limit"
                except Exception as exc:
                    status, record["error"] = "verification_error", str(exc)
                    record["reason"] = "operational_error"
            self.assert_sources()
            record["seconds"] = time.perf_counter() - started
            self.finish(claim, status, record)
        return True

    def retry(self, state: str) -> int:
        destinations = {"compute_incomplete": "queued", "compute_error": "queued",
                        "verification_incomplete": "candidate", "verification_error": "candidate"}
        if state not in destinations:
            raise ValueError("only incomplete or operationally failed jobs are retryable")
        with self.transaction():
            for row in self.db.execute("SELECT key FROM async_jobs WHERE policy=? AND state=?",
                                       (self.policy_key, state)).fetchall():
                self._attempt(row[0], "retry", destinations[state], {"previous_state": state})
            return self.db.execute("UPDATE async_jobs SET state=?,updated=? WHERE policy=? AND state=?",
                                   (destinations[state], time.time(), self.policy_key, state)).rowcount

    def counts(self) -> dict[str, Any]:
        rows = dict(self.db.execute("SELECT state,count(*) FROM async_jobs WHERE policy=? GROUP BY state",
                                    (self.policy_key,)).fetchall())
        rows["stored_candidate_bytes"] = self.db.execute("SELECT coalesce(sum(bytes),0) FROM async_candidates WHERE policy=?",
                                                        (self.policy_key,)).fetchone()[0]
        return rows
