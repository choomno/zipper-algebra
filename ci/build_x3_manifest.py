"""Rebuild the repository's X3 provenance inventory from explicit file families.

This is an audit of source, retained scientific data and test/build support. It
is not a list of files a user must manipulate to run a fresh computation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

# Keep selection rules visible here rather than relying on the current Git
# index, caches, an enclosing workspace or an unversioned generator.
FAMILIES = (
    ("src/zipper_algebra/**/*.py", "Maintained implementation used by presentation, production or verification"),
    ("tests/**/*.py", "Mathematical regression, installation or independent implementation comparison"),
    ("ci/*.py", "Installation, test, coverage, documentation or inventory support"),
    ("examples/x3/*", "Sparse current X3 output example and its reproduction instructions"),
)

EXACT = {
    ".coveragerc": "Complete line and branch measurement, including Python subprocesses",
    ".gitattributes": "Canonical line endings and solver-input classification",
    ".gitignore": "Generated computation files and caches stay untracked",
    ".github/workflows/ci.yml": "Hosted software tests, small solver integration and documentation checks",
    "MANIFEST.in": "Source-distribution inclusion rules",
    "pyproject.toml": "Package dependencies and installed commands",
    "setup.py": "Build-time source provenance embedded in wheels and source archives",
    "experiments/output.py": "Ignore markers for generated research output directories",
    "experiments/__init__.py": "Research recipe package marker",
    "patches/f4ncgb/msolve-parser-fix.patch": "Pinned external-solver parser correction",
    "patches/f4ncgb/LICENSE": "License accompanying the solver correction",
    "singular_run/x3-membership-Q.sing": "Independent rational-specialization membership cross-check in Singular",
}

REQUIRED = (
    "src/zipper_algebra/presentation.py",
    "src/zipper_algebra/cofactor/verify.py",
    "examples/x3/certificate.txt",
    "singular_run/x3-membership-Q.sing",
)

def build(root=ROOT):
    root = Path(root).resolve()
    roles = {}
    for pattern, role in FAMILIES:
        for path in root.glob(pattern):
            if path.is_file() and "__pycache__" not in path.parts:
                roles[path.relative_to(root).as_posix()] = role
    for name, role in EXACT.items():
        if (root / name).is_file():
            roles[name] = role
    for name in REQUIRED:
        if name not in roles:
            raise ValueError(f"Required reproduction input or source is absent: {name}")
    files = []
    for name, role in sorted(roles.items()):
        path = root / name
        if not path.resolve().is_relative_to(root):
            raise ValueError(f"Inventory path leaves the checkout: {name}")
        content = path.read_bytes()
        files.append({"path": name, "bytes": len(content),
                      "sha256": hashlib.sha256(content).hexdigest(), "role": role})
    return {
        "schema": "zipper-x3-reproduction-files-v1",
        "scope": "Curated X3 source and fixture inventory. Fresh production uses the installed discovery command.",
        "excluded": "This inventory, reader documentation, private working records, workspace configuration, caches and retired duplicate files.",
        "files": files,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true", help="require an already current inventory")
    args = parser.parse_args(argv)
    output = args.output or args.root / "docs/x3-reproduction-files.json"
    data = build(args.root)
    if args.check:
        if json.loads(output.read_text(encoding="utf-8")) != data:
            raise SystemExit("X3 inventory differs from its source families; rebuild with ci/build_x3_manifest.py")
        print(f"X3 inventory regeneration agrees: {len(data['files'])} entries")
        return
    output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(data['files'])} audit entries to {output}")


if __name__ == "__main__":
    main()
