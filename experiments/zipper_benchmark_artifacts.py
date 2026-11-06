"""Lossless SQLite storage for per-case zipper benchmark certificates."""

from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import zlib


def _alias(value: str) -> str:
    path = Path(value)
    if path.is_absolute() or not path.parts or ".." in path.parts:
        raise ValueError("artifact alias must be a relative path without traversal")
    return path.as_posix()


def load(database: Path, alias: str) -> bytes:
    """Load and hash-check a certificate by its original relative path."""
    with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as db:
        row = db.execute(
            "SELECT sha256,raw_bytes,payload FROM certificates WHERE alias=?",
            (_alias(alias),),
        ).fetchone()
    if row is None:
        raise KeyError(alias)
    digest, size, payload = row
    raw = zlib.decompress(payload)
    if len(raw) != size or hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError("certificate hash or size mismatch")
    return raw


def aliases(database: Path, run_alias: str, pass_index: int | None = None) -> list[str]:
    """List certificate aliases in stable pass/case order."""
    query = "SELECT alias FROM certificates WHERE run_alias=?"
    args: list[object] = [_alias(run_alias)]
    if pass_index is not None:
        query += " AND pass=?"
        args.append(pass_index)
    with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as db:
        return [row[0] for row in db.execute(query + " ORDER BY pass,case_index", args)]


def archive(directory: Path, *, remove: bool = False) -> dict[str, object]:
    """Archive generated files; optionally remove only byte-verified originals."""
    directory = directory.resolve()
    database = directory / "benchmark-results.sqlite3"
    files = sorted(directory.glob("*/pass*-case*.json"))
    with sqlite3.connect(database) as db:
        db.execute(
            "CREATE TABLE IF NOT EXISTS runs(alias TEXT PRIMARY KEY, summary TEXT NOT NULL)"
        )
        db.execute(
            "CREATE TABLE IF NOT EXISTS certificates(alias TEXT PRIMARY KEY,run_alias TEXT,pass INTEGER,case_index INTEGER,sha256 TEXT,raw_bytes INTEGER,payload BLOB)"
        )
        for path in sorted(directory.glob("*/summary.json")):
            json.loads(path.read_text())
            db.execute(
                "INSERT OR REPLACE INTO runs VALUES (?,?)",
                (str(path.parent.relative_to(directory)), path.read_text()),
            )
        for path in files:
            match = re.fullmatch(r"pass(\d+)-case(\d+)\.json", path.name)
            if match is None:
                raise ValueError(f"invalid certificate filename: {path.name}")
            raw = path.read_bytes()
            json.loads(raw)
            alias = _alias(str(path.relative_to(directory)))
            digest = hashlib.sha256(raw).hexdigest()
            previous = db.execute(
                "SELECT sha256 FROM certificates WHERE alias=?", (alias,)
            ).fetchone()
            if previous is not None and previous[0] != digest:
                raise ValueError(f"conflicting artifact alias: {alias}")
            db.execute(
                "INSERT OR REPLACE INTO certificates VALUES (?,?,?,?,?,?,?)",
                (
                    alias,
                    str(path.parent.relative_to(directory)),
                    int(match[1]),
                    int(match[2]),
                    digest,
                    len(raw),
                    zlib.compress(raw, 9),
                ),
            )
    for path in files:
        if load(database, str(path.relative_to(directory))) != path.read_bytes():
            raise ValueError("archive differs from source certificate")
    if remove:
        for path in files:
            path.unlink()
    return {
        "database": str(database),
        "verified_files": len(files),
        "removed_files": len(files) if remove else 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--remove-verified-originals", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            archive(args.directory, remove=args.remove_verified_originals), indent=2
        )
    )


if __name__ == "__main__":
    main()
