"""Restore optional historical research files from a checked archive.

Run from a source checkout with ``python -m experiments.restore_data``.
Files are copied as bytes; archived Python objects are never deserialized.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
from urllib.parse import urlparse
from urllib.request import Request, urlopen

from experiments.output import prepare_output_dir


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "research-archive.json"
GROUPS = ("all", "inputs", "x3-history", "variants", "gap", "documents")


class RestoreError(ValueError):
    """The manifest, archive or destination failed validation."""


@dataclass(frozen=True)
class Entry:
    path: str
    bytes: int
    sha256: str
    group: str


def _checked_path(name: str) -> str:
    if (not isinstance(name, str) or not name or "\\" in name
            or any(ord(character) < 32 or 127 <= ord(character) < 160 for character in name)
            or name.startswith("/") or re.match(r"^[A-Za-z]:", name)
            or any(part in ("", ".", "..") for part in name.split("/"))):
        raise RestoreError(f"Unsafe archive path: {name!r}")
    return name


def _load_manifest(path: Path, group: str) -> tuple[dict, list[Entry]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(document, dict) or document.get("schema") != "zipper-research-archive-v1":
        raise RestoreError("Unsupported research archive manifest schema")
    if not isinstance(document.get("groups"), dict):
        raise RestoreError("Manifest groups must be an object")
    if group != "all" and group not in document["groups"]:
        raise RestoreError(f"Unknown archive group: {group}")
    revision = document.get("source_revision")
    if not isinstance(revision, str) or not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise RestoreError("Manifest must pin a full Git source revision")
    if not isinstance(document.get("files"), list):
        raise RestoreError("Manifest files must be an array")
    entries = []
    seen = set()
    for record in document["files"]:
        if not isinstance(record, dict):
            raise RestoreError("Manifest file entries must be objects")
        name = _checked_path(record.get("path"))
        if name in seen:
            raise RestoreError(f"Duplicate manifest path: {name}")
        seen.add(name)
        size, digest, category = (record.get(key) for key in ("bytes", "sha256", "group"))
        if (type(size) is not int or size < 0 or not isinstance(digest, str)
                or not re.fullmatch(r"[0-9a-f]{64}", digest)
                or not isinstance(category, str) or category not in document["groups"]):
            raise RestoreError(f"Invalid manifest entry: {name}")
        if group == "all" or group == category:
            entries.append(Entry(name, size, digest, category))
    return document, entries


def _matches(path: Path, entry: Entry) -> bool:
    if not path.is_file() or path.stat().st_size != entry.bytes:
        return False
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest() == entry.sha256


def _target(root: Path, name: str) -> Path:
    target = root.joinpath(*name.split("/"))
    current = root
    for part in name.split("/"):
        current = current / part
        if current.is_symlink():
            raise RestoreError(f"Destination contains a symbolic link: {current}")
        if current != target and current.exists() and not current.is_dir():
            raise RestoreError(f"Destination parent is not a directory: {current}")
    return target


def _archive_git_blobs(source_root: Path, revision: str, entries: list[Entry], stream) -> None:
    # cat-file returns the stored bytes without checkout line-ending filters or
    # archive attributes such as export-subst and export-ignore.
    queries = "".join(f"{revision}:{entry.path}\n" for entry in entries).encode("utf-8")
    result = subprocess.run(
        ["git", "-C", str(source_root), "cat-file", "--batch"], input=queries,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if result.returncode:
        raise RestoreError("Cannot read pinned Git blobs: "
                           + result.stderr.decode(errors="replace").strip())
    blobs = io.BytesIO(result.stdout)
    with tarfile.open(fileobj=stream, mode="w") as bundle:
        for entry in entries:
            header = blobs.readline().split()
            if len(header) != 3 or header[1] != b"blob" or not header[2].isdigit():
                raise RestoreError(f"Cannot read pinned Git blob: {entry.path}")
            size = int(header[2])
            if size != entry.bytes:
                raise RestoreError(f"Size mismatch: {entry.path}")
            content = blobs.read(size)
            if len(content) != size or blobs.read(1) != b"\n":
                raise RestoreError(f"Incomplete pinned Git blob: {entry.path}")
            member = tarfile.TarInfo(entry.path)
            member.size = size
            bundle.addfile(member, io.BytesIO(content))
    if blobs.read(1):
        raise RestoreError("Unexpected data after pinned Git blobs")


@contextmanager
def _source(document: dict, entries: list[Entry], archive: Path | None, source_root: Path):
    if archive is not None:
        yield archive
        return
    revision = document["source_revision"]
    try:
        available = subprocess.run(
            ["git", "-C", str(source_root), "cat-file", "-e", revision + "^{commit}"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
        ).returncode == 0
    except FileNotFoundError:
        available = False
    with tempfile.TemporaryDirectory(prefix="zipper-archive-") as directory:
        bundle = Path(directory) / "source.tar"
        with bundle.open("wb") as stream:
            if available:
                _archive_git_blobs(source_root, revision, entries, stream)
            else:
                url = document.get("source_archive_url", "")
                parsed = urlparse(url)
                if parsed.scheme != "https" or not parsed.netloc or revision not in parsed.path:
                    raise RestoreError("Manifest download URL must use HTTPS and the pinned revision")
                request = Request(url, headers={"User-Agent": "zipper-algebra-data-restore"})
                with urlopen(request, timeout=60) as response:
                    shutil.copyfileobj(response, stream)
        yield bundle


def _verify_archive(archive: Path, entries: list[Entry], staging: Path) -> dict[str, Path]:
    wanted = {entry.path: entry for entry in entries}
    staged = {}
    with tarfile.open(archive, mode="r:*") as bundle:
        for member in bundle:
            # A GitHub source tarball adds one directory above repository paths.
            name = member.name.rstrip("/") if member.isdir() else member.name
            _checked_path(name)
            if not (member.isfile() or member.isdir()):
                raise RestoreError(f"Archive member is not a regular file or directory: {name}")
            candidates = [name]
            if "/" in name:
                candidates.append(name.split("/", 1)[1])
            matches = [candidate for candidate in candidates if candidate in wanted]
            if len(matches) > 1:
                raise RestoreError(f"Ambiguous archive member: {name}")
            if not matches:
                continue
            relative = matches[0]
            if relative in staged:
                raise RestoreError(f"Duplicate archive member for: {relative}")
            entry = wanted[relative]
            if not member.isfile():
                raise RestoreError(f"Expected a regular file: {relative}")
            if member.size != entry.bytes:
                raise RestoreError(f"Size mismatch: {relative}")
            temporary = staging / str(len(staged))
            digest = hashlib.sha256()
            with bundle.extractfile(member) as source, temporary.open("wb") as target:
                while chunk := source.read(1024 * 1024):
                    digest.update(chunk)
                    target.write(chunk)
            if digest.hexdigest() != entry.sha256:
                raise RestoreError(f"SHA256 mismatch: {relative}")
            staged[relative] = temporary
    missing = wanted.keys() - staged.keys()
    if missing:
        raise RestoreError("Archive is missing selected files: " + ", ".join(sorted(missing)))
    return staged


def restore_data(output: Path | str, group: str = "all", archive: Path | str | None = None,
                 *, manifest: Path | str = MANIFEST, source_root: Path | str = ROOT) -> list[Path]:
    """Restore selected files, preserving existing files and their original paths.

    Existing identical files are reused. Any differing destination is an error.
    All selected archive contents are verified before files are installed.
    """
    document, entries = _load_manifest(Path(manifest), group)
    root = Path(output).absolute()
    if root.is_symlink():
        raise RestoreError(f"Destination is a symbolic link: {root}")
    root = root.resolve()
    if root.exists() and not root.is_dir():
        raise RestoreError(f"Destination is not a directory: {root}")
    pending = []
    for entry in entries:
        target = _target(root, entry.path)
        if target.exists():
            if not _matches(target, entry):
                raise RestoreError(f"Existing destination differs; preserved: {target}")
        else:
            pending.append(entry)
    prepare_output_dir(root)
    if not pending:
        return [_target(root, entry.path) for entry in entries]
    with _source(document, pending, Path(archive) if archive is not None else None,
                 Path(source_root)) as bundle, \
            tempfile.TemporaryDirectory(prefix=".zipper-restore-", dir=root) as directory:
        staged = _verify_archive(bundle, pending, Path(directory))
        for entry in pending:
            target = _target(root, entry.path)
            target.parent.mkdir(parents=True, exist_ok=True)
            try:
                # Linking a verified temporary file installs it atomically and
                # fails if another writer has created the destination meanwhile.
                os.link(staged[entry.path], target)
            except FileExistsError:
                if not _matches(target, entry):
                    raise RestoreError(f"Existing destination differs; preserved: {target}") from None
    return [_target(root, entry.path) for entry in entries]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(".research-data"),
                        help="destination directory (default: .research-data)")
    parser.add_argument("--group", choices=GROUPS, default="all")
    parser.add_argument("--archive", type=Path, help="use a local tar or tar.gz bundle")
    args = parser.parse_args(argv)
    try:
        paths = restore_data(args.output, args.group, args.archive)
    except (RestoreError, OSError, tarfile.TarError, json.JSONDecodeError) as error:
        parser.exit(1, f"Cannot restore research data: {error}\n")
    print(f"Verified {len(paths)} research files in {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
