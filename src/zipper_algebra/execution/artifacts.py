"""Versioned proof-validated solver artifacts and atomic publication records."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile

from ..cofactor.proof import Basis, Proof


SCHEMA = "zipper-solver-run-v1"
ARTIFACT_NAMES = ("input.ms", "basis.txt", "proof.txt", "stdout.txt", "stderr.txt")


class ArtifactError(ValueError):
    """A recorded computation does not match its declared inputs or proof."""


def _contained(path: Path, parent: Path):
    if not path.resolve().is_relative_to(parent.resolve()):
        raise ArtifactError("artifact location leaves its configured directory")


def canonical_json(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def digest_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def file_identity(path: Path) -> dict:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        while chunk := stream.read(131072):
            size += len(chunk)
            digest.update(chunk)
    return {"bytes": size, "sha256": digest.hexdigest()}


def _check_syntax(path: Path, names, *, proof: bool):
    """Consume complete solver terms before using the historical tolerant parser."""
    variable = "(?:" + "|".join(re.escape(name) for name in names) + ")"
    scalar = r"\d+(?:/\d+)?"
    word = variable + r"(?:\*" + variable + ")*"
    if proof:
        term = r"[+-]?(?:" + scalar + r"\*)?(?:" + variable + r"\*)*\[i?\d+\](?:\*" + variable + ")*"
    else:
        term = r"[+-]?(?:" + scalar + r"(?:\*" + word + ")?|" + word + ")"
    pattern = re.compile(term)
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.endswith(","):
            line = line[:-1]
        tokens = line.split()
        pending = False
        for index, token in enumerate(tokens):
            if token in ("+", "-"):
                if pending:
                    raise ValueError("repeated standalone sign in solver output")
                pending = True
            else:
                if pending and token.startswith(("+", "-")):
                    raise ValueError("repeated sign in solver output")
                if index and not pending and not token.startswith(("+", "-")):
                    raise ValueError("missing sign between solver output terms")
                if pattern.fullmatch(token) is None:
                    raise ValueError("malformed solver output term")
                pending = False
        if pending:
            raise ValueError("incomplete solver output term")


def validate_proofs(directory: Path, inputs, names) -> int:
    """Check every output identity against the actual ordered input relations.

    This establishes that the output polynomials belong to the input ideal.
    It does not establish that they form a complete Gröbner basis.
    """
    try:
        _check_syntax(directory / "basis.txt", names, proof=False)
        _check_syntax(directory / "proof.txt", names, proof=True)
        basis = Basis.load(directory / "basis.txt", 0, names)
        proof = Proof.load(directory / "proof.txt", names)
        if not len(basis) or len(basis) != len(proof):
            raise ValueError("basis and proof must contain the same positive number of entries")
        for terms in proof.terms:
            for _coefficient, _left, kind, index, _right in terms:
                bound = len(inputs) if kind == "input" else len(basis)
                if not 0 <= index < bound:
                    raise ValueError("a raw proof reference does not name an existing relation")
        for index in range(len(basis)):
            if not proof.reconstructs(index, basis, inputs):
                raise ValueError(f"basis entry {index} fails its input-ideal identity")
    except (OSError, ValueError, IndexError, ZeroDivisionError) as error:
        raise ArtifactError(f"invalid solver proof: {error}") from error
    return len(basis)


def load_record(root: Path, key: str, request: dict, inputs, names):
    """Load a completed generation; absent records return None, corrupt ones fail."""
    slot = root / key
    _contained(slot, root)
    path = slot / "record.json"
    if not path.exists():
        return None
    try:
        record = json.loads(path.read_text())
        if (record["schema"] != SCHEMA
                or canonical_json(record["request"]) != canonical_json(request)
                or record["status"] != "proof-validated" or record["basis_complete"] is not False):
            raise ValueError("schema, request identity or proof status differs")
        generation = record["generation"]
        if (not isinstance(generation, str) or len(generation) != 64
                or any(c not in "0123456789abcdef" for c in generation)):
            raise ValueError("invalid generation identifier")
        directory = root / key / generation
        _contained(directory, slot)
        files = record["files"]
        if set(files) != set(ARTIFACT_NAMES):
            raise ValueError("artifact list differs")
        if digest_bytes(canonical_json(files)) != generation:
            raise ValueError("generation does not identify the artifact hashes")
        for name in ARTIFACT_NAMES:
            if file_identity(directory / name) != files[name]:
                raise ValueError(f"artifact hash differs: {name}")
        if files["input.ms"]["sha256"] != request["input_sha256"]:
            raise ValueError("saved input differs from the requested input")
        count = validate_proofs(directory, inputs, names)
        if type(record["basis_entries"]) is not int or record["basis_entries"] != count:
            raise ValueError("basis entry count differs")
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise ArtifactError(f"invalid run record at {path}: {error}") from error
    return directory, record


def _same_files(directory, files):
    try:
        return all(file_identity(directory / name) == files[name] for name in ARTIFACT_NAMES)
    except OSError:
        return False


def publish(root: Path, key: str, source: Path, request: dict, basis_entries: int):
    """Publish immutable bytes, then atomically select their completed generation.

    Interruption before replacing record.json leaves an older selected generation
    intact. Unselected generations cannot be discovered as completed samples.
    """
    root.mkdir(parents=True, exist_ok=True)
    slot = root / key
    _contained(slot, root)
    slot.mkdir(exist_ok=True)
    files = {name: file_identity(source / name) for name in ARTIFACT_NAMES}
    generation = digest_bytes(canonical_json(files))
    final = slot / generation
    _contained(final, slot)
    with tempfile.TemporaryDirectory(prefix=".publish-", dir=slot) as temporary:
        stage = Path(temporary) / "generation"
        stage.mkdir()
        for name in ARTIFACT_NAMES:
            shutil.copy2(source / name, stage / name)
        if final.exists() and not _same_files(final, files):
            shutil.rmtree(final)
        if not final.exists():
            try:
                os.replace(stage, final)
            except OSError:
                # A concurrent writer can publish the same immutable generation
                # after our existence check. Accept only its exact verified bytes.
                _contained(final, slot)
                if not final.exists() or not _same_files(final, files):
                    raise
        record = {"schema": SCHEMA, "request": request, "status": "proof-validated",
                  "basis_complete": False, "basis_entries": basis_entries,
                  "generation": generation, "files": files}
        marker = Path(temporary) / "record.json"
        marker.write_bytes(canonical_json(record) + b"\n")
        os.replace(marker, slot / "record.json")
    return final, record
