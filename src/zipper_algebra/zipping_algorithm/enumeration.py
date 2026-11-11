"""Exhaust literal braid words in length layers with resumable exact checking.

The signed alphabet is enumerated without cancellation or isotopy pruning.
The maintained complete reducer currently supports ranks one and two. Every
produced compact checkpoint is reloaded and checked before its SQLite commit.
The wall budget is soft and checked between bounded single-word computations.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path
import sqlite3
import time
from typing import Any, Iterator, Sequence

from .catalog import Catalog
from .two_strand import TwoStrandSession


def _integer(value: int, name: str, minimum: int) -> None:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer at least {minimum}")


def alphabet(n: int) -> tuple[int, ...]:
    _integer(n, "n", 1)
    return tuple(range(1-n, 0)) + tuple(range(1, n))


def word_count(n: int, length: int) -> int:
    _integer(length, "length", 0)
    return len(alphabet(n)) ** length


def iter_words(n: int, length: int) -> Iterator[tuple[int, ...]]:
    """Yield every literal signed-generator word; adjacent inverses are retained."""
    _integer(length, "length", 0)
    yield from itertools.product(alphabet(n), repeat=length)


def _word_at(n: int, length: int, ordinal: int) -> tuple[int, ...]:
    letters = alphabet(n)
    digits = []
    for _ in range(length):
        ordinal, digit = divmod(ordinal, len(letters))
        digits.append(letters[digit])
    return tuple(reversed(digits))


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def _identity(ranks: Sequence[int], trace: bool) -> dict[str, Any]:
    return {"schema": "zipper-literal-enumeration-v1", "ranks": list(ranks), "trace": trace,
            "coefficient_ring": "Z[q,q^-1]", "catalog_sha256": Catalog.bundled().sha256,
            "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                              for name in ("enumeration.py", "two_strand.py", "catalog.py")}}


def _normalize(cursor: list[int], supported: list[int]) -> list[int]:
    length, rank_index, ordinal = cursor
    while supported and (rank_index >= len(supported) or word_count(supported[rank_index], length) == 0):
        if rank_index >= len(supported):
            length, rank_index, ordinal = length+1, 0, 0
        else:
            rank_index += 1
        # B_1 has only its empty literal word; no further layer is nonempty.
        if supported == [1] and length > 0:
            return [length, 0, 0]
    return [length, rank_index, ordinal]


def _advance(cursor: list[int], supported: list[int]) -> list[int]:
    length, rank_index, ordinal = cursor
    ordinal += 1
    if ordinal == word_count(supported[rank_index], length):
        rank_index, ordinal = rank_index+1, 0
    return _normalize([length, rank_index, ordinal], supported)


def _audit(connection: sqlite3.Connection, supported: list[int], stored_cursor: list[int], identity: dict[str, Any], session: TwoStrandSession) -> int:
    """Check that the database is an exact checked prefix, including its last checkpoint."""
    cursor, count, incomplete = [0, 0, 0], 0, False
    for row in connection.execute("SELECT n,length,ordinal,rank_order,word,input_sha,status,certificate,certificate_sha,config_sha,metrics,diagnostic FROM words ORDER BY length,rank_order,ordinal"):
        n, length, ordinal, rank_index, word_json, ident, status, encoded, digest, configuration_sha, metrics_json, diagnostic = row
        if incomplete or not supported or [length, rank_index, ordinal] != cursor or n != supported[rank_index]:
            raise ValueError("saved words are not a consecutive literal enumeration")
        word = list(_word_at(n, length, ordinal))
        if word_json != _json(word) or ident != _sha(_json({"n": n, "word": word})) or digest != _sha(encoded):
            raise ValueError("saved literal word or certificate hash mismatch")
        certificate = json.loads(encoded)
        configuration_row = connection.execute("SELECT config FROM configurations WHERE hash=?", (configuration_sha,)).fetchone()
        if not configuration_row:
            raise ValueError("saved word configuration is missing")
        configuration = json.loads(configuration_row[0])
        if (_sha(configuration_row[0]) != configuration_sha or configuration["identity_sha256"] != _sha(_json(identity))
                or certificate["trace_enabled"] != identity["trace"]
                or certificate["limits"] != {key: configuration[key] for key in ("max_steps", "max_terms")}):
            raise ValueError("saved word configuration hash or limits mismatch")
        checked = session.check(certificate)
        if (checked["pass"] is not True or certificate["input_word"] != word or certificate["n"] != n
                or certificate["status"] != status or certificate["diagnostic"] != diagnostic
                or _json(certificate["metrics"]) != metrics_json):
            raise ValueError("saved literal certificate does not match its input or status")
        count += 1
        incomplete = status != "complete"
        if not incomplete:
            cursor = _advance(cursor, supported)
    if cursor != stored_cursor:
        raise ValueError("saved cursor does not match the checked word prefix")
    return count


def _layers(connection: sqlite3.Connection) -> list[dict[str, Any]]:
    return [{"n": n, "length": length, "expected_words": word_count(n, length),
             "visited_words": count, "complete_words": complete,
             "layer_complete": complete == count == word_count(n, length)}
            for n, length, count, complete in connection.execute(
                "SELECT n,length,count(*),sum(status='complete') FROM words GROUP BY n,length ORDER BY length,rank_order")]


def exhaust(output: str | Path, *, ranks: Sequence[int] = (1, 2), wall_seconds: float = 60,
            max_length: int | None = None, max_steps: int = 5000, max_terms: int = 1024,
            trace: bool = False) -> dict[str, Any]:
    """Extend an exact literal prefix until a wall budget, word limit, or length cap.

    ``words.sqlite3`` stores exact inputs, ordinals, compact certificates, hashes,
    metrics, and timings. A limit-hit word remains the next cursor and requires
    larger per-word limits on resume. Source, rank order, and trace mode are bound
    to the database. Wall time and limits may change on later invocations.
    """
    if not isinstance(ranks, (list, tuple)) or not ranks or len(set(ranks)) != len(ranks):
        raise ValueError("ranks must be a nonempty sequence of distinct positive integers")
    for n in ranks:
        _integer(n, "rank", 1)
    if type(wall_seconds) not in (int, float) or not math.isfinite(wall_seconds) or wall_seconds <= 0:
        raise ValueError("wall_seconds must be finite and positive")
    if max_length is not None:
        _integer(max_length, "max_length", 0)
    _integer(max_steps, "max_steps", 0)
    _integer(max_terms, "max_terms", 1)
    if type(trace) is not bool:
        raise ValueError("trace must be a boolean")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    supported = [n for n in ranks if n <= 2]
    identity = _identity(ranks, trace)
    session = TwoStrandSession()
    config = {"wall_seconds": wall_seconds, "max_length": max_length,
              "max_steps": max_steps, "max_terms": max_terms, "identity_sha256": _sha(_json(identity))}
    config_sha = _sha(_json(config))
    database = output / "words.sqlite3"
    connection = sqlite3.connect(database)
    try:
        connection.execute("CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY,value TEXT NOT NULL)")
        connection.execute("CREATE TABLE IF NOT EXISTS words (n INTEGER,length INTEGER,ordinal INTEGER,rank_order INTEGER,word TEXT,input_sha TEXT,status TEXT,diagnostic TEXT,certificate TEXT,certificate_sha TEXT,metrics TEXT,producer_ns INTEGER,checker_ns INTEGER,config_sha TEXT, PRIMARY KEY(n,length,ordinal))")
        connection.execute("CREATE TABLE IF NOT EXISTS configurations (hash TEXT PRIMARY KEY,config TEXT NOT NULL)")
        metadata = dict(connection.execute("SELECT key,value FROM metadata"))
        if metadata:
            if metadata.get("identity") != _json(identity):
                raise ValueError("enumeration source, catalog, rank order, or trace configuration changed")
            cursor = json.loads(metadata["cursor"])
        else:
            cursor = [0, 0, 0]
            with connection:
                connection.executemany("INSERT INTO metadata VALUES (?,?)", [("identity", _json(identity)), ("cursor", _json(cursor))])
        audit_start = time.perf_counter()
        audited = _audit(connection, supported, cursor, identity, session)
        audit_seconds = time.perf_counter()-audit_start
        with connection:
            connection.execute("INSERT OR IGNORE INTO configurations VALUES (?,?)", (config_sha, _json(config)))
        start = time.perf_counter()
        new_words = 0
        while True:
            length, rank_index, ordinal = cursor
            if not supported:
                reason = "unsupported_ranks"
                break
            if supported == [1] and length > 0:
                reason = "all_nonempty_layers_exhausted"
                break
            if max_length is not None and length > max_length:
                reason = "max_length"
                break
            if time.perf_counter()-start >= wall_seconds:
                reason = "wall_budget_between_words"
                break
            n = supported[rank_index]
            word = list(_word_at(n, length, ordinal))
            previous = connection.execute("SELECT certificate FROM words WHERE n=? AND length=? AND ordinal=?", (n,length,ordinal)).fetchone()
            produced = time.perf_counter_ns()
            if previous:
                checkpoint = json.loads(previous[0])
                old = checkpoint["limits"]
                if max_steps < old["max_steps"] or max_terms < old["max_terms"] or (max_steps == old["max_steps"] and max_terms == old["max_terms"]):
                    reason = "word_limit_needs_larger_limits"
                    break
                certificate = session.resume(checkpoint, max_steps=max_steps, max_terms=max_terms)
            else:
                certificate = session.reduce(word, n=n, max_steps=max_steps, max_terms=max_terms, trace=trace)
            producer_ns = time.perf_counter_ns()-produced
            encoded = _json(certificate)
            checked_at = time.perf_counter_ns()
            checked = session.check(json.loads(encoded))
            checker_ns = time.perf_counter_ns()-checked_at
            if checked["pass"] is not True:
                raise ValueError("produced word certificate failed exact checking")
            status = certificate["status"]
            next_cursor = _advance(cursor, supported) if status == "complete" else cursor
            with connection:
                connection.execute("INSERT OR REPLACE INTO words VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (n,length,ordinal,rank_index,_json(word),_sha(_json({"n":n,"word":word})),status,
                     certificate["diagnostic"],encoded,_sha(encoded),_json(certificate["metrics"]),
                     producer_ns,checker_ns,config_sha))
                connection.execute("UPDATE metadata SET value=? WHERE key='cursor'", (_json(next_cursor),))
                connection.execute("INSERT OR REPLACE INTO metadata VALUES ('last_config',?)", (_json(config),))
            cursor = next_cursor
            new_words += 1
            if status != "complete":
                reason = "word_limit"
                break
        summary = {"schema": identity["schema"], "identity": identity, "config": config,
                   "config_sha256": config_sha, "database": str(database.resolve()),
                   "audited_words": audited, "new_or_resumed_words": new_words,
                   "audit_seconds": audit_seconds, "new_work_seconds": time.perf_counter()-start,
                   "stop_reason": reason, "next_cursor": {"length":cursor[0],"rank_index":cursor[1],"ordinal":cursor[2]},
                   "unsupported_ranks": [n for n in ranks if n > 2], "layers": _layers(connection)}
        temporary = output / "summary.json.tmp"
        temporary.write_text(_json(summary)+"\n")
        temporary.replace(output / "summary.json")
        return summary
    finally:
        connection.close()
