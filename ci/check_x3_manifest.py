"""Check the paths and SHA-256 hashes in the X3 reproduction inventory."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "docs/x3-reproduction-files.json"


def check(root=ROOT, manifest=MANIFEST):
    root = Path(root).resolve()
    data = json.loads(Path(manifest).read_text(encoding="utf-8"))
    if data.get("schema") != "zipper-x3-reproduction-files-v1":
        raise ValueError("unknown reproduction inventory format")
    files = data["files"]
    if not isinstance(files, list) or not files:
        raise ValueError("the reproduction inventory must list files")
    seen = set()
    for entry in files:
        name = entry["path"]
        relative = PurePosixPath(name)
        if relative.is_absolute() or ".." in relative.parts or "\\" in name or name in seen:
            raise ValueError(f"invalid or repeated inventory path: {name}")
        seen.add(name)
        path = root / name
        if not path.resolve().is_relative_to(root):
            raise ValueError(f"inventory path leaves the checkout: {name}")
        content = path.read_bytes()
        if len(content) != entry["bytes"]:
            raise ValueError(f"file size differs: {name}")
        if hashlib.sha256(content).hexdigest() != entry["sha256"]:
            raise ValueError(f"SHA-256 differs: {name}")
        if not isinstance(entry.get("role"), str) or not entry["role"].strip():
            raise ValueError(f"missing reproduction role: {name}")
    print(f"X3 reproduction inventory passed: {len(files)} files")
    return len(files)


if __name__ == "__main__":
    try:
        check()
    except (ValueError, KeyError, OSError, TypeError) as error:
        raise SystemExit(f"X3 reproduction inventory failed: {error}")
