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
HISTORICAL_REVISION = "0db1dd9e2a33fb0e5d07e216c8a5354911446e97"

# Keep selection rules visible here rather than relying on the current Git
# index, caches, an enclosing workspace or an unversioned generator.
FAMILIES = (
    ("src/zipper_algebra/**/*.py", "Maintained implementation used by presentation, production or verification"),
    ("tests/**/*.py", "Mathematical regression, installation or independent implementation comparison"),
    ("tests/fixtures/legacy_x3/*", "Curated historical X3 regression fixture or its provenance manifest"),
    ("ci/*.py", "Installation, test, coverage, documentation or inventory support"),
    ("experiments/x3/*.py", "Compact X3 transformations and checks"),
    ("examples/x3/*", "Sparse current X3 output example and its reproduction instructions"),
)

EXACT = {
    ".coveragerc": "Complete line and branch measurement, including Python subprocesses",
    ".gitattributes": "Stable source line endings and byte-preserving fixture rules",
    ".gitignore": "Generated runs, caches and historical datasets stay untracked",
    ".github/workflows/ci.yml": "Hosted software tests, small solver integration and documentation checks",
    "MANIFEST.in": "Source-distribution inclusion rules",
    "pyproject.toml": "Package dependencies and installed commands",
    "experiments/output.py": "Ignore markers for generated research output directories",
    "experiments/__init__.py": "Research recipe package marker",
    "patches/f4ncgb/msolve-parser-fix.patch": "Pinned external-solver parser correction",
    "patches/f4ncgb/LICENSE": "License accompanying the solver correction",
    "singular_run/x3-membership-Q.sing": "Independent specialized Singular producer for the canonical presentation",
}

REQUIRED = (
    "src/zipper_algebra/presentation.py",
    "src/zipper_algebra/cofactor/verify.py",
    "tests/fixtures/legacy_x3/proof-n4-Q-noex.txt",
    "tests/fixtures/legacy_x3/basis-symbolic-Qqt.json",
    "tests/fixtures/legacy_x3/skeleton-n4-Q.json",
    "tests/fixtures/legacy_x3/basis-n4-Q.txt",
    "tests/fixtures/legacy_x3/cofactor-X3-n4-Qqt-symbolic.txt",
    "examples/x3/certificate.txt",
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
    retired = (
        "extract_cofactor.py", "extract_cofactor_qt.py", "symbolic_lift_x3.py",
        "compose_x3_cofactor_qt.py", "lift_until_verified.py", "run_lift_grid.py",
        "run_lift_grid_canonical.py", "reproduce_X3_zero.sh",
    )
    prefix = "archive/cofactor-pre-consolidation/"
    return {
        "schema": "zipper-x3-reproduction-files-v1",
        "scope": "Curated X3 source and fixture inventory. Fresh production uses the installed discovery command.",
        "excluded": "This inventory, reader documentation, private working records, workspace configuration, caches and retired duplicate files.",
        "historical_sources": [
            {"path": prefix + name, "revision": HISTORICAL_REVISION,
             "url": f"https://github.com/choomno/zipper-algebra/blob/{HISTORICAL_REVISION}/{prefix}{name}",
             "role": "Retired historical producer; source retained in the phase 1 Git baseline"}
            for name in retired
        ],
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
