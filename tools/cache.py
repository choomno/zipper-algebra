#!/usr/bin/env python3
"""Content-addressed cache for pipeline artifacts.

An artifact is cached under a key derived from **everything that determines it**:
the input file's content, the tool flags, and the tool's own identity. If the key
matches, the artifact is reused; if any input changes, the key changes and the
artifact is recomputed. There is no staleness to reason about.

Why this exists: the expensive phases are 1b (f4ncgb at an anchor) and 3 (the
multi-point lift, ~226 specialisations, multi-hour). Phases 4-10 are seconds.
Re-running the tail of the pipeline should not re-run the head. Before this,
`reproduce_X3_zero.sh` hard-coded a single cached Phase 3 file and recomputed
Phase 1 every time.

Location comes from `[paths].cache` in settings.toml, so it is a machine fact
like every other path.

Shell use (the driver does this):
    KEY=$(python3 cache.py key --file input.ms --tag f4ncgb -- -e 0 -v 2 -r 0)
    if python3 cache.py get --key "$KEY" --into basis.txt proof.txt; then
        echo "[cache] hit"
    else
        ... compute ...
        python3 cache.py put --key "$KEY" -- basis.txt proof.txt
    fi

Python use:
    from cache import Cache
    c = Cache()
    k = c.key(files=[ms_path], tag="f4ncgb", args=["-e","0"])
    if not c.get(k, [basis, proof]):
        run_it(); c.put(k, [basis, proof])
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]          # tools/ -> zipper-algebra/ -> repo root
sys.path.insert(0, str(REPO_ROOT))


def _cache_dir() -> Path:
    try:
        from repo_settings import settings
        return settings().path("cache")
    except Exception:
        d = REPO_ROOT / "zipper-algebra" / ".cache"
        d.mkdir(parents=True, exist_ok=True)
        return d


class Cache:
    def __init__(self, root: Path | None = None):
        self.root = Path(root) if root else _cache_dir()
        self.root.mkdir(parents=True, exist_ok=True)

    def key(self, files=(), tag: str = "", args=(), tools=()) -> str:
        """Hash of every input that determines the artifact — INCLUDING the code
        that produces it.

        `tools` are the producing binaries/scripts. Without them a cache hit only
        means "same inputs", and a changed tool would silently return an artifact
        the current revision would never produce. That is the failure class that
        gave this repo a Jun-2 symbolic lift feeding today's phases, and a
        37-sample cofactor that disagreed with the exact answer on all 695 slots.
        With them, a hit means "this IS what the current revision produces".
        """
        h = hashlib.sha256()
        h.update(tag.encode())
        for a in args:
            h.update(b"\0"); h.update(str(a).encode())
        for f in files:
            p = Path(f)
            h.update(b"\0file\0"); h.update(p.name.encode())
            h.update(p.read_bytes() if p.exists() else b"<missing>")
        for t in tools:
            p = Path(t)
            h.update(b"\0tool\0"); h.update(p.name.encode())
            if p.exists():
                # hash content; for a large binary this is still well under a second
                hh = hashlib.sha256()
                with open(p, "rb") as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        hh.update(chunk)
                h.update(hh.digest())
            else:
                h.update(b"<missing>")
        return h.hexdigest()[:16]

    def _slot(self, key: str) -> Path:
        return self.root / key

    def get(self, key: str, dests) -> bool:
        """Copy cached artifacts to dests. True on a complete hit."""
        slot = self._slot(key)
        man = slot / "manifest.json"
        if not man.exists():
            return False
        names = json.loads(man.read_text())["names"]
        if len(names) != len(dests):
            return False
        if not all((slot / n).exists() for n in names):
            return False
        for n, d in zip(names, dests):
            Path(d).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(slot / n, d)
        return True

    # NOTE: artifacts are stored under an index-prefixed name, so two inputs with
    # the same basename cannot overwrite each other inside a slot -- that would
    # have returned one artifact's bytes in another's place.

    def put(self, key: str, srcs) -> None:
        slot = self._slot(key)
        slot.mkdir(parents=True, exist_ok=True)
        names = []
        for i, s in enumerate(srcs):
            p = Path(s)
            if not p.exists():
                raise SystemExit(f"cache put: missing artifact {p}")
            stored = f"{i:02d}-{p.name}"      # index prefix: basenames may collide
            shutil.copy2(p, slot / stored)
            names.append(stored)
        (slot / "manifest.json").write_text(json.dumps({"names": names}, indent=1))

    def entries(self):
        for slot in sorted(self.root.iterdir()):
            man = slot / "manifest.json"
            if man.is_file():
                names = json.loads(man.read_text())["names"]
                size = sum((slot / n).stat().st_size for n in names if (slot / n).exists())
                yield slot.name, names, size


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    k = sub.add_parser("key"); k.add_argument("--file", action="append", default=[])
    k.add_argument("--tool", action="append", default=[],
                   help="producing binary/script; its content enters the key")
    k.add_argument("--tag", default=""); k.add_argument("rest", nargs="*")

    g = sub.add_parser("get"); g.add_argument("--key", required=True)
    g.add_argument("--into", nargs="+", required=True)

    p = sub.add_parser("put"); p.add_argument("--key", required=True)
    p.add_argument("rest", nargs="+")

    sub.add_parser("list")
    c = sub.add_parser("clear"); c.add_argument("--key", default=None)

    a = ap.parse_args()
    cache = Cache()

    if a.cmd == "key":
        print(cache.key(files=a.file, tag=a.tag, args=a.rest, tools=a.tool)); return 0
    if a.cmd == "get":
        return 0 if cache.get(a.key, a.into) else 1
    if a.cmd == "put":
        cache.put(a.key, a.rest); return 0
    if a.cmd == "list":
        total = 0
        for key, names, size in cache.entries():
            print(f"  {key}  {size:>10,} B  {', '.join(names)}"); total += size
        print(f"  {'':16}  {total:>10,} B total  in {cache.root}")
        return 0
    if a.cmd == "clear":
        if a.key:
            shutil.rmtree(cache.root / a.key, ignore_errors=True)
        else:
            shutil.rmtree(cache.root, ignore_errors=True)
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
