"""Inspect release archives before installation; needs only the standard library."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile
import zipfile


LEGACY_X3_FILES = {
    "basis-symbolic-Qqt.json", "skeleton-n4-Q.json", "proof-n4-Q-noex.txt",
    "basis-n4-Q.txt", "basis-n4-Q-q3-t11.txt", "basis-n4-Q-q7-t13.txt",
    "cofactor-X3-n4-Qqt-symbolic.txt",
}


def single(directory, pattern):
    found = list(Path(directory).glob(pattern))
    if len(found) != 1:
        raise SystemExit(f"expected one {pattern} in {directory}, found {len(found)}")
    return found[0]


def build_provenance(content, files, read, prefix):
    """Check the packaged build origin against the packaged implementation bytes."""
    metadata = json.loads(content)
    assert set(metadata) == {"schema", "origin"}
    assert metadata["schema"] == "zipper-build-provenance-v1"
    origin = metadata["origin"]
    assert origin["schema"] == "zipper-source-provenance-v1"
    canonical_sources = json.dumps(origin["source_files"], sort_keys=True,
                                   separators=(",", ":"), allow_nan=False).encode()
    assert hashlib.sha256(canonical_sources).hexdigest() == origin["source_sha256"]
    package_prefix = "src/zipper_algebra/"
    sources = {name[len(prefix):] for name in files
               if name.startswith(prefix) and name.endswith(".py")}
    recorded = {name[len(package_prefix):] for name in origin["source_files"]
                if name.startswith(package_prefix)}
    assert sources == recorded, sources ^ recorded
    for name in sources:
        actual = read(prefix + name)
        identity = origin["source_files"][package_prefix + name]
        assert len(actual) == identity["bytes"], name
        assert hashlib.sha256(actual).hexdigest() == identity["sha256"], name
    return metadata


def inspect(directory):
    wheel = single(directory, "*.whl")
    source = single(directory, "*.tar.gz")
    distribution, version, *_ = wheel.name.split("-")
    metadata = f"{distribution}-{version}.dist-info/"
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        wheel_names = list(names)
        if not all(name.startswith("zipper_algebra/") or name.startswith(metadata)
                   for name in names):
            raise AssertionError("the wheel contains material outside the public package")
        for needed in ("__init__.py", "cli.py", "presentation.py", "provenance.py", "core.py",
                       "coefficients.py", "polynomial.py", "relations.py", "msolve.py",
                       "linear.py", "cofactor/verify.py", "cofactor/lift.py",
                       "cofactor/inference.py", "cofactor/univariate.py",
                       "cofactor/artifacts.py", "execution/config.py",
                       "execution/runner.py", "execution/artifacts.py",
                       "execution/sampling.py", "execution/pipeline.py"):
            assert f"zipper_algebra/{needed}" in names, needed
        assert metadata + "licenses/LICENSE" in names, "the wheel must include the license"
        assert b"License-Expression: MIT" in archive.read(metadata + "METADATA")
        wheel_origin = build_provenance(
            archive.read("zipper_algebra/_build_provenance.json"), names,
            archive.read, "zipper_algebra/")
    with tarfile.open(source) as archive:
        names = archive.getnames()
        roots = {PurePosixPath(name).parts[0] for name in names}
        assert len(roots) == 1, "source distribution must have one enclosing directory"
        relative = {str(PurePosixPath(*PurePosixPath(name).parts[1:])) for name in names}
        for needed in ("pyproject.toml", "setup.py", "README.md", "INDEX.md", "WORKFLOW.md",
                       ".github/workflows/ci.yml", "singular_run/x3-membership-Q.sing",
                       "singular_run/README.md", "reconstruction_theory.tex", "docs/TESTING.md",
                       "docs/X3_REPRODUCTION.md", "docs/PROVENANCE.md", "docs/x3-reproduction-files.json",
                       "patches/f4ncgb/msolve-parser-fix.patch", "patches/f4ncgb/LICENSE",
                       "docs/mathematical-specification/main.tex",
                       "docs/mathematical-specification/mathematical-specification.pdf",
                       "tests/reproduction/test_x3.py",
                       "tests/reproduction/test_legacy_proofs.py",
                       "tests/reference/legacy_x3.py", "examples/infer_denominators.py",
                       "examples/x3/certificate.txt", "examples/x3/README.md",
                       "tests/fixtures/legacy_x3/manifest.json",
                       "experiments/survey/pointalg.py",
                       "experiments/research/theta_min_poly.py",
                       "experiments/variants/square_zero_regular_rep_qt.py",
                       "ci/run_tests.py", "ci/build_x3_manifest.py"):
            assert needed in relative, needed
        fixture_prefix = "tests/fixtures/legacy_x3/"
        expected_fixtures = {fixture_prefix + name for name in LEGACY_X3_FILES | {"manifest.json"}}
        actual_fixtures = {name for name in relative if name.startswith(fixture_prefix)}
        assert actual_fixtures == expected_fixtures, actual_fixtures ^ expected_fixtures
        root = next(iter(roots))
        def read(name):
            member = archive.extractfile(name)
            assert member is not None, name
            return member.read()
        source_origin = build_provenance(
            read(f"{root}/src/zipper_algebra/_build_provenance.json"), names,
            read, f"{root}/src/zipper_algebra/")
        assert source_origin == wheel_origin, "wheel and source archive build origins differ"
        recorded = source_origin["origin"]["source_files"]
        recipes = {name for name in relative if name.startswith("experiments/")
                   and name.endswith(".py")}
        recorded_recipes = {name for name in recorded if name.startswith("experiments/")}
        assert recipes == recorded_recipes, recipes ^ recorded_recipes
        for name in recipes | {"setup.py", "pyproject.toml", "MANIFEST.in"}:
            payload = read(f"{root}/{name}")
            assert {"bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest()} \
                == recorded[name], name
        manifest_file = archive.extractfile(f"{root}/{fixture_prefix}manifest.json")
        assert manifest_file is not None
        manifest = json.load(manifest_file)
        assert manifest["source_commit"] == "d30b36aa2da0819dfa28954bd9223ef72163edf5"
        assert {entry["path"] for entry in manifest["files"]} == LEGACY_X3_FILES
        assert len(manifest["files"]) == len(LEGACY_X3_FILES)
        for entry in manifest["files"]:
            payload = archive.extractfile(f"{root}/{fixture_prefix}{entry['path']}")
            assert payload is not None
            assert hashlib.sha256(payload.read()).hexdigest() == entry["sha256"], entry["path"]
            source_path = PurePosixPath(entry["source_path"])
            assert not source_path.is_absolute() and ".." not in source_path.parts
            assert "source_record" not in entry
    forbidden = {".git", ".cache", "__pycache__", "private-baselines", "research-records",
                 "recovery", "manuscripts", "historical-reconstruction", "logs_and_traces",
                 "logs", "ms_inputs",
                 ".claude", ".codex", ".agents", ".agent-notes", "agent-notes",
                 "CONTRIBUTING.md", "CONTRIBUTING.MD", "AGENTS.md", "CLAUDE.md"}
    forbidden.update({"PRIVATE_BASELINE.md", "PUBLIC_SCOPE.md", "VALIDATION.md", "HISTORY.md",
                      "baseline.json", "public-scope.json", "PLAN_FOR_REVIEW.md", "REPORT.md",
                      "AUDIT_spec_vs_generic.md", "GENERALITY_CENSUS.md", "LEDGER_claims.md",
                      "NIGHT_LOG.md", "REPTHEORY_STRATEGY.md", "RESEARCH_AGENDA.md",
                      "SOURCE_IMPORT.md", "SURVEY_WORKFLOW.md"})
    for name in wheel_names + names:
        path = PurePosixPath(name)
        assert not path.is_absolute() and ".." not in path.parts and "\\" not in name, name
        assert not forbidden.intersection(path.parts), name
        assert path.suffix not in {".pyc", ".pyo", ".pkl", ".pickle"}, name
    print(f"Distribution contents passed: {wheel.name}, {source.name}")


if __name__ == "__main__":
    inspect(sys.argv[1] if len(sys.argv) > 1 else "dist")
