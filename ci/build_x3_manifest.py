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
    ("src/zipper_algebra/**/*.py", "Maintained presentation, arithmetic, production or verification source"),
    ("tests/**/*.py", "Regression, integration, reproduction or frozen implementation comparison"),
    ("tests/fixtures/legacy_x3/*", "Immutable symbolic lift, reduction skeleton or canonical anchor proof"),
    ("ci/*.py", "Installation, test, coverage, documentation or inventory validation support"),
    ("examples/*.py", "Executable mathematical example used by regression checks"),
    ("experiments/**/*.py", "Retained research recipe or independent computational comparison"),
    ("ms_inputs/qt_sweep/*.ms", "Retained specialized solver input or reference presentation"),
    ("logs_and_traces/f4ncgb/basis-n4-Q*.txt", "Retained specialized solver basis"),
    ("logs_and_traces/f4ncgb/cofactor-X3-*.txt", "Retained numeric cofactor or symbolic certificate stage"),
    ("archive/gap/traces/grobner_zipper_bigelow-untraced-[34]-*.gap", "Saved finite-field basis for the GBNP comparison"),
)

EXACT = {
    ".coveragerc": "Complete line and branch measurement, including Python subprocesses",
    ".github/workflows/ci.yml": "Hosted test, solver integration, coverage and documentation execution",
    "MANIFEST.in": "Source-distribution inclusion rules",
    "pyproject.toml": "Package dependencies and installed command entry point",
    "data/aliases.json": "Hashes, canonical destinations and roles of deduplicated historical data",
    "patches/f4ncgb/msolve-parser-fix.patch": "Pinned external-solver parser correction",
    "patches/f4ncgb/LICENSE": "License accompanying the solver correction",
    "singular_run/x3-membership-Q.sing": "Preserved canonical 31-relation Singular reduction producer at q=5, t=7",
    "logs_and_traces/reproduce-X3-zero-transcript.txt": "Historical corrected 31-relation production transcript",
    "logs_and_traces/f4ncgb/zipper-n4-Q-q5-t7.ms": "Original bare 31-relation canonical anchor input",
    "logs_and_traces/f4ncgb/atoms-n4-Q.json": "Historical denominator-family reconstruction prior",
    "logs_and_traces/f4ncgb/skeleton-n4-Q.json": "Recorded canonical anchor reduction skeleton",
    "logs_and_traces/f4ncgb/basis-symbolic-Qqt.json": "Retained newer-schema symbolic basis-proof reconstruction",
    "logs_and_traces/f4ncgb/basis-symbolic-Qqt.txt": "Historical symbolic lift rendering from a different revision",
    "logs_and_traces/f4ncgb/lift_samples/coefficients.json": "212 numeric cofactor samples with 695 coefficient slots",
    "logs_and_traces/f4ncgb/reproduce-zipper-n4-Q-q5-t7.ms": "Older augmented input: 31 canonical relations plus X3 itself",
    "logs_and_traces/f4ncgb/reproduce-basis-n4-Q.txt": "131-entry basis of the older target-augmented presentation",
    "logs_and_traces/f4ncgb/reproduce-proof-n4-Q-noex.txt": "131-entry augmented-presentation proof using appended target input",
    "logs_and_traces/f4ncgb/X3-membership-singular-independent.txt": "Independent specialized Singular membership summary",
}

for name in ("Q_d12", "Q_d14", "GFp_d10", "GFp_d14"):
    EXACT[f"logs_and_traces/f4ncgb/singular_X3_{name}.log"] = "Independent specialized Singular execution output"
for name in ("zipper", "traced_X3_witness_n3", "reduce_X3_against_B", "trace_replay", "test_n4_presentation"):
    EXACT[f"archive/gap/{name}.gap"] = "Historical independent GBNP presentation, witness or replay route"

REQUIRED = (
    "src/zipper_algebra/presentation.py",
    "src/zipper_algebra/cofactor/verify.py",
    "tests/fixtures/legacy_x3/proof-n4-Q-noex.txt",
    "tests/fixtures/legacy_x3/basis-symbolic-Qqt.json",
    "tests/fixtures/legacy_x3/skeleton-n4-Q.json",
    "logs_and_traces/f4ncgb/zipper-n4-Q-q5-t7.ms",
    "logs_and_traces/f4ncgb/basis-n4-Q.txt",
    "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt",
    "data/aliases.json",
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
    aliases = json.loads((root / "data/aliases.json").read_text())
    if aliases.get("schema") != "zipper-data-aliases-v1":
        raise ValueError("Unknown historical-data alias schema")
    for alias in aliases["aliases"]:
        name = alias["canonical_path"]
        path = root / name
        if not path.resolve().is_relative_to(root):
            raise ValueError(f"Canonical data alias leaves the checkout: {name}")
        content = path.read_bytes()
        if (len(content) != alias["bytes"]
                or hashlib.sha256(content).hexdigest() != alias["sha256"]):
            raise ValueError(f"Canonical data differs from its recorded historical alias: {name}")
        roles.setdefault(name, alias["role"])
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
        "scope": "Repository audit inventory: maintained source, retained historical evidence, and test/build support. Fresh production is managed by the installed discovery command.",
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
