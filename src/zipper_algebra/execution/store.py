"""Transactional storage of exact solver artifacts in one SQLite database.

The database records proof-validation results, not basis completeness. Readers
check byte identities; the execution layer independently rechecks proof identities.
"""
from __future__ import annotations

from contextlib import closing, contextmanager
import json
import os
from pathlib import Path
import sqlite3
import tempfile

from .artifacts import (ARTIFACT_NAMES, SCHEMA as RECORD_SCHEMA, ArtifactError,
                        canonical_json, digest_bytes, prepare_output_root,
                        validate_source_snapshot)
from ..provenance import validate_provenance


SCHEMA = "zipper-sample-store-v1"
APPLICATION_ID = 1514758226
VERSION = 1
_REQUEST_V3_FIELDS = {"schema", "presentation", "point", "specialized_domain",
                      "monomial_order", "input_sha256", "producer", "solver"}
_TABLES = {
    "metadata": "CREATE TABLE metadata (schema TEXT NOT NULL)",
    "generations": ("CREATE TABLE generations (key TEXT NOT NULL, generation TEXT NOT NULL, "
                    "record BLOB NOT NULL, PRIMARY KEY (key, generation))"),
    "artifacts": ("CREATE TABLE artifacts (key TEXT NOT NULL, generation TEXT NOT NULL, "
                  "name TEXT NOT NULL, content BLOB NOT NULL, PRIMARY KEY (key, generation, name), "
                  "FOREIGN KEY (key, generation) REFERENCES generations (key, generation) ON DELETE CASCADE)"),
    "selections": ("CREATE TABLE selections (key TEXT PRIMARY KEY NOT NULL, generation TEXT NOT NULL, "
                   "FOREIGN KEY (key, generation) REFERENCES generations (key, generation))"),
    "f4ncgb_io": ("CREATE TABLE f4ncgb_io (job_id TEXT PRIMARY KEY,input BLOB NOT NULL,basis BLOB,"
                  "proof BLOB,state TEXT NOT NULL CHECK(state IN('ready','complete')))"),
}


def _hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _require_hash(value, name):
    if not _hash(value):
        raise ArtifactError(f"invalid sample {name}")


def _request_key(request):
    if (not isinstance(request, dict) or request.get("schema") not in
            ("zipper-solver-request-v2", "zipper-solver-request-v3")
            or not _hash(request.get("input_sha256"))):
        raise ArtifactError("invalid sample request")
    if request["schema"] == "zipper-solver-request-v3" and set(request) != _REQUEST_V3_FIELDS:
        raise ArtifactError("sample request-v3 must contain exactly the presentation computation fields")
    try:
        validate_provenance(request.get("producer"))
        # Source metadata and configuration must be finite, serializable JSON.
        json.dumps(request, allow_nan=False)
    except (ValueError, TypeError) as error:
        raise ArtifactError(f"invalid sample request: {error}") from error
    return digest_bytes(canonical_json(request))


def _file_identities(files):
    if not isinstance(files, dict) or set(files) != set(ARTIFACT_NAMES):
        raise ArtifactError("sample artifact names differ")
    if any(not isinstance(content, bytes) for content in files.values()):
        raise ArtifactError("sample artifacts must be exact bytes")
    return {name: {"bytes": len(content), "sha256": digest_bytes(content)}
            for name, content in files.items()}


def _validate_schema(connection):
    if (connection.execute("PRAGMA application_id").fetchone()[0] != APPLICATION_ID
            or connection.execute("PRAGMA user_version").fetchone()[0] != VERSION):
        raise ArtifactError("not a supported zipper-algebra sample database")
    objects = connection.execute(
        "SELECT name, type, sql FROM sqlite_master WHERE name NOT GLOB 'sqlite_*'").fetchall()
    if {name for name, _, _ in objects} != set(_TABLES):
        raise ArtifactError("sample database objects differ from its schema")
    if any(kind != "table" or " ".join(sql.split()) != " ".join(_TABLES[name].split())
           for name, kind, sql in objects):
        raise ArtifactError("sample database table definitions differ")
    if connection.execute("SELECT schema FROM metadata").fetchall() != [(SCHEMA,)]:
        raise ArtifactError("sample database schema marker differs")


def _recover_hot_journal(path):
    """Let SQLite restore committed pages before retrying read-only access."""
    try:
        with closing(sqlite3.connect(Path(path).resolve().as_uri() + "?mode=rw", uri=True,
                                     timeout=30, isolation_level=None)) as connection:
            with connection:
                connection.execute("BEGIN IMMEDIATE")
                connection.execute("SELECT name FROM sqlite_master").fetchall()
    except sqlite3.Error as error:
        raise ArtifactError(f"sample database has a hot journal and needs writable recovery: {error}") from error


@contextmanager
def _connection(path, *, create=False):
    path = Path(path)
    if path.is_symlink():
        raise ArtifactError("sample database must not be a symbolic link")
    if not create and not path.exists():
        yield None
        return
    if create:
        prepare_output_root(path.parent)
    connection = None
    try:
        connection = (sqlite3.connect(path, timeout=30, isolation_level=None) if create else
                      sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True,
                                      timeout=30, isolation_level=None))
        connection.execute("PRAGMA foreign_keys=ON")
        if not create:
            connection.execute("PRAGMA query_only=ON")
        connection.execute("BEGIN IMMEDIATE" if create else "BEGIN")
        try:
            objects = connection.execute("SELECT name FROM sqlite_master").fetchall()
        except sqlite3.Error as error:
            if create or getattr(error, "sqlite_errorcode", None) != sqlite3.SQLITE_READONLY_ROLLBACK:
                raise
            connection.close()
            connection = None
            _recover_hot_journal(path)
            connection = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True,
                                         timeout=30, isolation_level=None)
            connection.execute("PRAGMA foreign_keys=ON")
            connection.execute("PRAGMA query_only=ON")
            connection.execute("BEGIN")
            objects = connection.execute("SELECT name FROM sqlite_master").fetchall()
        if not objects:
            uninitialized = (connection.execute("PRAGMA application_id").fetchone()[0] == 0
                             and connection.execute("PRAGMA user_version").fetchone()[0] == 0)
            if not uninitialized:
                raise ArtifactError("empty database has an unrelated schema marker")
            if not create:
                # A concurrent initial publisher, or an interrupted schema
                # initialization, may leave an empty SQLite file temporarily.
                yield None
                connection.commit()
                return
            for statement in _TABLES.values():
                connection.execute(statement)
            connection.execute("INSERT INTO metadata VALUES (?)", (SCHEMA,))
            connection.execute(f"PRAGMA application_id={APPLICATION_ID}")
            connection.execute(f"PRAGMA user_version={VERSION}")
            # Schema creation has its own atomic commit. A later failed sample
            # publication therefore leaves a valid, reusable empty database.
            connection.commit()
            connection.execute("BEGIN IMMEDIATE")
        _validate_schema(connection)
        yield connection
        connection.commit()
    except sqlite3.Error as error:
        raise ArtifactError(f"invalid or unavailable sample database: {error}") from error
    finally:
        if connection is not None:
            if connection.in_transaction:
                connection.rollback()
            connection.close()


def _read(connection, key, generation=None):
    if generation is None:
        selected = connection.execute("SELECT generation FROM selections WHERE key=?", (key,)).fetchone()
        if selected is None:
            return None
        generation = selected[0]
    _require_hash(generation, "generation")
    row = connection.execute("SELECT record FROM generations WHERE key=? AND generation=?",
                             (key, generation)).fetchone()
    if row is None:
        raise ArtifactError("selected sample generation is missing")
    try:
        if not isinstance(row[0], bytes):
            raise ValueError("record is not a JSON byte string")
        record = json.loads(row[0])
        fields = {"schema", "request", "status", "basis_complete", "basis_entries", "generation", "files"}
        if (not isinstance(record, dict) or set(record) != fields
                or record["schema"] != RECORD_SCHEMA or record["status"] != "proof-validated"
                or record["basis_complete"] is not False or type(record["basis_entries"]) is not int
                or record["basis_entries"] < 1 or record["generation"] != generation
                or _request_key(record["request"]) != key):
            raise ValueError("record identity or proof status differs")
        files = dict(connection.execute(
            "SELECT name, content FROM artifacts WHERE key=? AND generation=?", (key, generation)))
        identities = _file_identities(files)
        if (canonical_json(identities) != canonical_json(record["files"])
                or digest_bytes(canonical_json(identities)) != generation
                or identities["input.ms"]["sha256"] != record["request"]["input_sha256"]):
            raise ValueError("artifact bytes do not match their recorded identities")
    except (ValueError, TypeError, KeyError) as error:
        raise ArtifactError(f"invalid stored sample: {error}") from error
    return files, record


@contextmanager
def transaction(path, *, write=False):
    """Open a short, schema-checked transport transaction.

    A write transaction initializes an absent store. A read transaction is
    read-only and yields None if the store is absent or not yet initialized.
    Native solver jobs must release this context before starting a subprocess.
    """
    with _connection(path, create=write) as connection:
        yield connection


def load(path, key, request):
    """Read a selected sample for this exact request; absent samples return None."""
    _require_hash(key, "key")
    if _request_key(request) != key:
        raise ArtifactError("sample key does not identify the requested computation")
    with _connection(path) as connection:
        result = None if connection is None else _read(connection, key)
    return result


def read(path, key, generation=None):
    """Read a selected or historical generation, checking metadata and all bytes."""
    _require_hash(key, "key")
    if generation is not None:
        _require_hash(generation, "generation")
    with _connection(path) as connection:
        result = None if connection is None else _read(connection, key, generation)
    if result is None:
        raise ArtifactError("sample is not present in the database")
    return result


def publish(path, key, files, request, basis_entries):
    """Atomically retain a validated generation and select it for its request key."""
    _require_hash(key, "key")
    if _request_key(request) != key:
        raise ArtifactError("sample key does not identify the requested computation")
    if type(basis_entries) is not int or basis_entries < 1:
        raise ArtifactError("sample basis count must be a positive integer")
    identities = _file_identities(files)
    files = dict(files)
    if identities["input.ms"]["sha256"] != request["input_sha256"]:
        raise ArtifactError("sample input differs from the requested input")
    generation = digest_bytes(canonical_json(identities))
    record = {"schema": RECORD_SCHEMA, "request": request, "status": "proof-validated",
              "basis_complete": False, "basis_entries": basis_entries,
              "generation": generation, "files": identities}
    with _connection(path, create=True) as connection:
        connection.execute("INSERT OR REPLACE INTO generations VALUES (?, ?, ?)",
                           (key, generation, canonical_json(record)))
        for name, content in files.items():
            connection.execute("INSERT OR REPLACE INTO artifacts VALUES (?, ?, ?, ?)",
                               (key, generation, name, content))
        validate_source_snapshot(request["producer"])
        connection.execute("INSERT OR REPLACE INTO selections VALUES (?, ?)", (key, generation))
    return record


def export(path, key, output, generation=None):
    """Export exact sample bytes and their record, preserving differing files."""
    files, record = read(path, key, generation)
    files = {**files, "record.json": canonical_json(record) + b"\n"}
    output = Path(output)
    prepare_output_root(output)
    for name, content in files.items():
        destination = output / name
        if destination.is_symlink() or destination.exists() and (
                not destination.is_file() or destination.read_bytes() != content):
            raise ArtifactError(f"export would overwrite a different file: {name}")
    for name, content in files.items():
        destination = output / name
        if destination.exists():
            continue
        stream = tempfile.NamedTemporaryFile(dir=output, prefix=".export-", delete=False)
        temporary = Path(stream.name)
        try:
            with stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
            try:
                os.link(temporary, destination)
            except FileExistsError:
                if destination.is_symlink() or destination.read_bytes() != content:
                    raise ArtifactError(f"export encountered a differing file: {name}")
        finally:
            temporary.unlink(missing_ok=True)
    return record
