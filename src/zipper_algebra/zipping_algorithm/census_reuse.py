"""Persistent exact census reuse with independently replayed primitive results.

The sidecar never imports or alters the legacy literal census. Only completed
leaves enter the reusable table. On each new process a stored leaf is replayed
once before being trusted, with a bounded LRU of checked content digests.
Reconstruction records identify outputs, not standalone execution traces.
"""
from __future__ import annotations

from collections import OrderedDict
from collections.abc import Callable
from copy import deepcopy
import hashlib
import fcntl
import json
from pathlib import Path
import sqlite3
from time import perf_counter
from typing import Any, NamedTuple
import zlib

from .input_normal_form import input_normal_form, leaf_requests, reconstruction_plan
from .native import _decode_coefficients, _encode_coefficients
from .native_verifier import VerifierIncomplete
from .normal_form_reconstruction import (
    completed_output, execution_cache_key, _identity, reconstruct_output,
)
from .output_normal_form import output_fingerprint, engine_output, _exact
from .reconstruction_check import check_reconstruction, identity_leaf_request

SCHEMA = "zipper-reuse-census-v1"


def _json(value: Any) -> str:
    return json.dumps(_encode_coefficients(value), sort_keys=True, separators=(",", ":"),
                      allow_nan=False)


def _sha(value: Any) -> str:
    return hashlib.sha256(_json(value).encode()).hexdigest()


class _ValidatedLeaf(NamedTuple):
    """Per-call validated representation; never retained across census calls."""
    certificate: dict[str, Any]
    normal: dict[str, Any]
    output_digest: str


def _canonical_digest(normal: dict[str, Any]) -> str:
    # Only called on the direct return value of completed_output: its output
    # already has canonical ordering and decimal coefficient strings.
    encoded = json.dumps(normal, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=True, allow_nan=False).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def reconstruction_sources() -> dict[str, str]:
    """Bind the complete Python package, including all witness dependencies."""
    directory = Path(__file__).parent
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.glob("*.py"))}


_IMPORT_SOURCES = reconstruction_sources()


def verifier_binding(metadata: dict[str, Any]) -> dict[str, Any]:
    """Keep authenticated verifier identity, excluding session timing and paths."""
    required = {"abi", "source_sha256", "binary_sha256"}
    if not required <= metadata.keys() or any(not metadata[key] for key in required):
        raise ValueError("authenticated verifier ABI, source, and binary identities required")
    fields = ("abi", "source_sha256", "reference_source_sha256", "wrapper_source_sha256",
              "binary_sha256", "compiler", "compiler_version", "flags", "link_libraries",
              "libraries", "allocator_library")
    return {key: deepcopy(metadata[key]) for key in fields if key in metadata}


class ReuseCensus:
    """Single-writer sidecar; one instance owns one producer and verifier pair.

    Callbacks must produce authentic results and independently replay them.
    The official CLI supplies NativeSession and CppVerifier directly. Different
    source, verifier, limit, and execution policies occupy distinct namespaces.
    """

    def __init__(self, path: str | Path, *, produce: Callable[[int, list[int]], dict[str, Any]],
                 verify: Callable[[dict[str, Any]], bool], producer_identity: dict[str, Any],
                 verifier_identity: dict[str, Any], limits: dict[str, int],
                 cache_entries: int = 4096) -> None:
        if type(cache_entries) is not int or cache_entries < 0:
            raise ValueError("cache_entries must be a nonnegative integer")
        if not verifier_identity or not limits or any(type(v) is not int for v in limits.values()):
            raise ValueError("explicit verifier identity and integer resource limits required")
        self.produce, self.verify = produce, verify
        self.producer_identity = deepcopy(producer_identity)
        self.limits = deepcopy(limits)
        self._sources = reconstruction_sources()
        if self._sources != _IMPORT_SOURCES:
            raise ValueError("census package sources changed since import")
        self.policy = {"schema": SCHEMA, "producer": self.producer_identity,
                       "verifier": deepcopy(verifier_identity), "limits": self.limits,
                       "reconstruction_sources": self._sources}
        self.policy_key = _sha(self.policy)
        self.cache_entries = cache_entries
        self._checked: OrderedDict[str, None] = OrderedDict()
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        # Probe existing files before a writable open: SQLite may otherwise
        # recover a legacy WAL/journal even if we subsequently reject it.
        if self.path.exists():
            probe = sqlite3.connect(self.path.resolve().as_uri() + "?mode=ro", uri=True)
            try:
                tables = {r[0] for r in probe.execute("SELECT name FROM sqlite_master WHERE type='table'")}
                if tables and "reuse_policies" not in tables:
                    raise ValueError("reuse census requires its own sidecar SQLite database")
            finally:
                probe.close()
        self._lock_file = self.path.with_name(self.path.name + ".lock").open("a")
        try:
            fcntl.flock(self._lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BaseException:
            self._lock_file.close()
            raise
        try:
            self.db = sqlite3.connect(self.path, timeout=60)
            # Refuse historical census files even though our table names differ.
            names = {r[0] for r in self.db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if names and "reuse_policies" not in names:
                self.db.close()
                self._lock_file.close()
                raise ValueError("reuse census requires its own sidecar SQLite database")
            self.db.execute("PRAGMA foreign_keys=ON")
            self.db.execute("PRAGMA journal_mode=WAL")
            self.db.executescript("""
                CREATE TABLE IF NOT EXISTS reuse_policies(
                    key TEXT PRIMARY KEY, payload TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS reuse_leaves(
                    policy TEXT NOT NULL REFERENCES reuse_policies(key),
                    key TEXT NOT NULL, request TEXT NOT NULL,
                    payload BLOB NOT NULL, digest TEXT NOT NULL,
                    output_digest TEXT NOT NULL,
                    PRIMARY KEY(policy,key));
                CREATE TABLE IF NOT EXISTS reuse_roots(
                    policy TEXT NOT NULL REFERENCES reuse_policies(key),
                    key TEXT NOT NULL, witness TEXT NOT NULL,
                    output_digest TEXT NOT NULL, PRIMARY KEY(policy,key));
                CREATE TABLE IF NOT EXISTS reuse_dependencies(
                    policy TEXT NOT NULL, root TEXT NOT NULL, leaf TEXT NOT NULL,
                    PRIMARY KEY(policy,root,leaf),
                    FOREIGN KEY(policy,root) REFERENCES reuse_roots(policy,key),
                    FOREIGN KEY(policy,leaf) REFERENCES reuse_leaves(policy,key));
                CREATE TABLE IF NOT EXISTS reuse_attempts(
                    id INTEGER PRIMARY KEY, policy TEXT NOT NULL REFERENCES reuse_policies(key),
                    root TEXT NOT NULL, status TEXT NOT NULL, payload TEXT NOT NULL);
            """)
            with self.db:
                self.db.execute("INSERT OR IGNORE INTO reuse_policies VALUES (?,?)",
                                (self.policy_key, _json(self.policy)))
            existing = self.db.execute("SELECT payload FROM reuse_policies WHERE key=?",
                                       (self.policy_key,)).fetchone()[0]
            if existing != _json(self.policy):
                self.db.close()
                self._lock_file.close()
                raise ValueError("stored reuse policy identity mismatch")
        except BaseException:
            connection = getattr(self, "db", None)
            if connection is not None:
                connection.close()
            self._lock_file.close()
            raise

    def _assert_sources(self) -> None:
        if reconstruction_sources() != self._sources:
            raise ValueError("census reconstruction sources changed during the run")

    def _remember(self, digest: str) -> None:
        self._checked[digest] = None
        self._checked.move_to_end(digest)
        while len(self._checked) > self.cache_entries:
            self._checked.popitem(last=False)

    def _validate_leaf(self, result: dict[str, Any], request: dict[str, Any]) -> _ValidatedLeaf:
        # This is metadata only. The following completed_output call performs
        # exactly one whole-output structural validation and canonicalization.
        metadata = _identity({"schema": result.get("schema"),
                              "source_sha256": result.get("source_sha256"),
                              "catalog_sha256": result.get("catalog_sha256"),
                              "execution": result.get("execution", {})})
        if (type(result.get("n")) is not int or type(result.get("input_word")) is not list
                or any(type(a) is not int for a in result["input_word"])
                or result.get("n") != request["n"] or result.get("input_word") != request["word"]
                or _json(metadata) != _json(self.producer_identity)
                or _json(result.get("limits")) != _json(self.limits)
                or result.get("status") != "complete"):
            raise ValueError("cached primitive result binding mismatch or incomplete result")
        normal = completed_output(result)
        return _ValidatedLeaf(result, normal, _canonical_digest(normal))

    def _check(self, result: dict[str, Any], stats: dict[str, Any]) -> bool:
        started = perf_counter()
        try:
            return self.verify(result)
        finally:
            stats["checker_seconds"] += perf_counter() - started

    def _leaf(self, request: dict[str, Any], stats: dict[str, Any]) -> _ValidatedLeaf | dict[str, Any]:
        key = execution_cache_key(request, self.producer_identity)
        started = perf_counter()
        row = self.db.execute("SELECT request,payload,digest,output_digest FROM reuse_leaves "
                              "WHERE policy=? AND key=?", (self.policy_key, key)).fetchone()
        if row is not None:
            encoded = zlib.decompress(row[1])
            digest = hashlib.sha256(encoded).hexdigest()
            if row[0] != _json(request) or digest != row[2]:
                raise ValueError("cached primitive payload or request digest mismatch")
            result = _decode_coefficients(json.loads(encoded))
            validated = self._validate_leaf(result, request)
            if validated.output_digest != row[3]:
                raise ValueError("cached primitive output digest mismatch")
            stats["store_seconds"] += perf_counter() - started
            stats["leaf_hits"] += 1
            if digest not in self._checked:
                if self._check(result, stats) is not True:
                    raise ValueError("independent verifier rejected stored primitive result")
                self._remember(digest)
            else:
                self._checked.move_to_end(digest)
            return validated
        stats["store_seconds"] += perf_counter() - started
        stats["leaf_misses"] += 1
        started = perf_counter()
        result = self.produce(request["n"], request["word"])
        stats["producer_seconds"] += perf_counter() - started
        if result.get("status") != "complete":
            # No cached result and no coverage claim for a resource-limited leaf.
            return result
        validated = self._validate_leaf(result, request)
        if self._check(result, stats) is not True:
            raise ValueError("independent verifier rejected primitive result")
        self._assert_sources()
        started = perf_counter()
        encoded = _json(result).encode()
        digest = hashlib.sha256(encoded).hexdigest()
        values = (self.policy_key, key, _json(request), zlib.compress(encoded), digest,
                  validated.output_digest)
        with self.db:
            self.db.execute("BEGIN IMMEDIATE")
            old = self.db.execute("SELECT request,digest,output_digest FROM reuse_leaves "
                                   "WHERE policy=? AND key=?", values[:2]).fetchone()
            if old is None:
                self.db.execute("INSERT INTO reuse_leaves VALUES (?,?,?,?,?,?)", values)
            elif old != (values[2], values[4], values[5]):
                raise ValueError("conflicting immutable primitive result")
        self._remember(digest)
        stats["store_seconds"] += perf_counter() - started
        return validated

    def _reconstruction_budget(self, form: Any, leaves: dict[str, Any]) -> str | None:
        """Conservative allocation bounds, before forming any Cartesian product."""
        plan = reconstruction_plan(form)
        terms, support, bits, word_size = 1, 1, 0, 0
        decoded = {key: engine_output(completed_output(value)) for key, value in leaves.items()}
        for atom in plan["atoms"]:
            if atom["kind"] == "singleton":
                continue
            choices = decoded[atom["key"]]
            terms *= len(choices)
            support *= max((len(t["coefficient"]) for t in choices), default=0)
            bits += max((abs(c).bit_length() for t in choices for _, c in t["coefficient"]), default=0)
            word_size += max((sum(len(c["component_word"]) for c in t["components"])
                              for t in choices), default=0)
            if terms > min(self.limits["max_terms"], self.limits["expansion_budget"]):
                return "reconstruction Cartesian product exceeds term budget"
            if support > min(self.limits["max_coefficient_support"], self.limits["expansion_budget"]):
                return "reconstruction polynomial product exceeds support budget"
        # Two stable-partition routes per extraction cross each pair at most twice.
        word_size += plan["n"] * (plan["n"] - 1)
        if terms and word_size > self.limits["max_word_size"]:
            return "reconstruction route bound exceeds word budget"
        if bits + max(0, (terms * support - 1).bit_length()) > self.limits["max_integer_bits"]:
            return "reconstruction coefficient bound exceeds integer budget"
        return None

    def _incomplete(self, root: str, payload: dict[str, Any], stats: dict[str, Any],
                    started: float) -> dict[str, Any]:
        with self.db:
            self.db.execute("INSERT INTO reuse_attempts(policy,root,status,payload) VALUES (?,?,?,?)",
                            (self.policy_key, root, "incomplete", _json(payload)))
        stats["wall_seconds"] = perf_counter() - started
        return {"status": "incomplete", "root": root, "policy": self.policy_key,
                "diagnostic": payload, "stats": stats}

    def evaluate(self, word: list[int] | tuple[int, ...], *, n: int) -> dict[str, Any]:
        """Return a checked canonical output plus separately identified provenance.

        Partial leaves return incomplete without storing a reusable leaf/root.
        Calls may be retried with the same or a larger policy-bound budget.
        """
        self._assert_sources()
        started = perf_counter()
        stats: dict[str, Any] = {"producer_seconds": 0.0, "checker_seconds": 0.0,
                                 "reconstruction_seconds": 0.0, "store_seconds": 0.0,
                                 "leaf_hits": 0, "leaf_misses": 0}
        form = input_normal_form(word, n=n)
        if (type(form) is not dict or type(form.get("root")) is not dict
                or type(form["root"].get("rank")) is not int or form["root"]["rank"] != n
                or _json(form["root"].get("word")) != _json(list(word))):
            raise ValueError("input reconstruction witness differs from the requested literal input")
        root = _sha({"n": n, "word": list(word)})
        leaves = {}
        validated_leaves = {}
        requests = leaf_requests(form)
        for request in requests:
            try:
                result = self._leaf(request, stats)
            except VerifierIncomplete as error:
                return self._incomplete(root, {"stage": "verification", "request": request,
                                              "replay_status": error.replay_status,
                                              "diagnostic": error.replay_diagnostic}, stats, started)
            if not isinstance(result, _ValidatedLeaf):
                payload = {"request": request, "leaf_status": result.get("status"),
                           "diagnostic": result.get("diagnostic"), "stats": stats}
                return self._incomplete(root, payload, stats, started)
            validated_leaves[request["key"]] = result
            leaves[request["key"]] = result.certificate
        before = perf_counter()
        identity_key = identity_leaf_request(form, requests)
        stats["identity_reconstruction"] = identity_key is not None
        if identity_key is not None:
            output = validated_leaves[identity_key].normal
        else:
            diagnostic = self._reconstruction_budget(form, leaves)
            if diagnostic:
                stats["reconstruction_seconds"] += perf_counter() - before
                return self._incomplete(root, {"stage": "reconstruction", "message": diagnostic}, stats, started)
            output = reconstruct_output(form, leaves, source_identity=self.producer_identity)
        decoded = output["terms"]
        if (len(decoded) > self.limits["max_terms"]
                or any(len(term["word"]) > self.limits["max_word_size"]
                       or len(term["coefficient"]) > self.limits["max_coefficient_support"]
                       or any(abs(_exact(c)).bit_length() > self.limits["max_integer_bits"]
                              for _, c in term["coefficient"]) for term in decoded)):
            return self._incomplete(root, {"stage": "reconstruction", "message": "ambient output limit"}, stats, started)
        if identity_key is None and not check_reconstruction(form, leaves, output, source_identity=self.producer_identity):
            raise ValueError("independent reconstruction check failed")
        stats["reconstruction_seconds"] += perf_counter() - before
        self._assert_sources()
        before = perf_counter()
        fingerprint = (validated_leaves[identity_key].output_digest if identity_key is not None
                       else output_fingerprint(output))
        with self.db:
            self.db.execute("BEGIN IMMEDIATE")
            old = self.db.execute("SELECT witness,output_digest FROM reuse_roots WHERE policy=? AND key=?",
                                  (self.policy_key, root)).fetchone()
            if old is None:
                self.db.execute("INSERT INTO reuse_roots VALUES (?,?,?,?)",
                                (self.policy_key, root, _json(form), fingerprint))
                for request in requests:
                    self.db.execute("INSERT INTO reuse_dependencies VALUES (?,?,?)",
                                    (self.policy_key, root,
                                     execution_cache_key(request, self.producer_identity)))
            elif old != (_json(form), fingerprint):
                raise ValueError("conflicting immutable reconstruction result")
            expected = {execution_cache_key(r, self.producer_identity) for r in requests}
            actual = {r[0] for r in self.db.execute("SELECT leaf FROM reuse_dependencies "
                                                   "WHERE policy=? AND root=?", (self.policy_key, root))}
            if expected != actual:
                raise ValueError("stored reconstruction dependency mismatch")
        stats["store_seconds"] += perf_counter() - before
        stats["wall_seconds"] = perf_counter() - started
        return {"schema": SCHEMA, "status": "complete", "root": root, "policy": self.policy_key,
                "output": output, "output_digest": fingerprint,
                "provenance": ("independently-verified-leaf-with-checked-identity-embedding"
                               if identity_key is not None else
                               "checked-reconstruction-from-independently-verified-leaves"),
                "stats": stats}

    def counts(self) -> dict[str, int]:
        return {name: self.db.execute(f"SELECT count(*) FROM reuse_{name} WHERE policy=?",
                                     (self.policy_key,)).fetchone()[0]
                for name in ("leaves", "roots", "attempts")}

    def close(self) -> None:
        self.db.close()
        self._lock_file.close()

    def __enter__(self) -> ReuseCensus:
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()
