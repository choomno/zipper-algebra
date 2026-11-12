"""Inspect release archives before installation; needs only the standard library."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import tarfile
import zipfile


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
    for native_prefix, package_bundle in (
            ("native/zipper_engine/", "zipping_algorithm/_native_sources/"),
            ("native/zipper_verifier/", "zipping_algorithm/_verifier_sources/")):
        payload_prefix = (prefix + package_bundle if prefix == "zipper_algebra/"
                          else prefix.removesuffix("src/zipper_algebra/") + native_prefix)
        native_files = {name[len(payload_prefix):] for name in files
                        if name.startswith(payload_prefix)
                        and (name.endswith((".cpp", ".hpp", ".h"))
                             or name.endswith("vendor/nlohmann/LICENSE.MIT"))}
        recorded_native = {name[len(native_prefix):] for name in origin["source_files"]
                           if name.startswith(native_prefix)}
        assert native_files == recorded_native, native_files ^ recorded_native
        for name in native_files:
            actual = read(payload_prefix + name)
            identity = origin["source_files"][native_prefix + name]
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
                       "execution/sampling.py", "execution/pipeline.py",
                       "verification/gap.py",
                       "zipping_algorithm/__init__.py", "zipping_algorithm/__main__.py",
                       "zipping_algorithm/catalog.py", "zipping_algorithm/interpreter.py",
                       "zipping_algorithm/transport.py", "zipping_algorithm/upward.py",
                       "zipping_algorithm/two_strand.py", "zipping_algorithm/native.py",
                       "zipping_algorithm/census_reuse.py",
                       "zipping_algorithm/reconstruction_check.py",
                       "zipping_algorithm/normal_form_reconstruction.py",
                       "zipping_algorithm/_native_sources/main.cpp",
                       "zipping_algorithm/_native_sources/engine.hpp",
                       "zipping_algorithm/_native_sources/vendor/nlohmann/json.hpp",
                       "zipping_algorithm/_native_sources/vendor/nlohmann/LICENSE.MIT",
                       "zipping_algorithm/native_verifier.py",
                       "zipping_algorithm/_verifier_sources/main.cpp",
                       "zipping_algorithm/_verifier_sources/core.hpp",
                       "zipping_algorithm/_verifier_sources/typed_state.hpp",
                       "zipping_algorithm/_verifier_sources/typed_algorithm.cpp",
                       "zipping_algorithm/_verifier_sources/vendor/nlohmann/json.hpp",
                       "zipping_algorithm/_verifier_sources/vendor/nlohmann/LICENSE.MIT",
                       "zipping_algorithm/data/full-catalog.json",
                       "zipping_algorithm/data/schema.json", "zipping_algorithm/data/manifest.json"):
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
                       "singular_run/README.md",
                       "docs/TESTING.md",
                       "docs/X3_REPRODUCTION.md", "docs/PROVENANCE.md", "docs/x3-reproduction-files.json",
                       "patches/f4ncgb/msolve-parser-fix.patch", "patches/f4ncgb/LICENSE",
                       "docs/mathematical-specification/main.tex",
                       "docs/mathematical-specification/mathematical-specification.pdf",
                       "tests/test_proof_graphs.py", "examples/infer_denominators.py",
                       "examples/x3/certificate.txt", "examples/x3/README.md",
                       "experiments/survey/pointalg.py",
                       "experiments/research/theta_min_poly.py",
                       "experiments/variants/square_zero_regular_rep_qt.py",
                       "ci/build_x3_manifest.py", "docs/ZIPPING_REFERENCE.md", "docs/ZIPPING_NATIVE.md",
                       "docs/zipper-census-reuse.md", "docs/zipper-computation-graph.md",
                       "docs/zipper-reuse-performance.md",
                       "docs/zipper-async-census.md", "docs/zipper-census-plan.md",
                       "docs/zipper-async-performance.md",
                       "docs/zipper-scheduling-performance.md",
                       "docs/benchmarks/schedule-2026-11-09.json",
                       "docs/zipper-derivation-certificates.md",
                       "experiments/zipper_exhaustive_reuse.py",
                       "experiments/zipper_census_reuse_benchmark.py",
                       "experiments/zipper_cold_epoch_benchmark.py",
                       "experiments/zipper_schedule_report.py",
                       "experiments/zipper_async_census.py", "experiments/zipper_async_benchmark.py",
                       "experiments/zipper_derivation_dag.py",
                       "native/zipper_engine/main.cpp", "native/zipper_engine/engine.hpp",
                       "native/zipper_engine/vendor/nlohmann/json.hpp",
                       "native/zipper_engine/vendor/nlohmann/LICENSE.MIT",
                       "docs/ZIPPING_NATIVE_VERIFIER.md",
                       "docs/ZIPPING_GPU_PROTOTYPE.md", "docs/ZIPPING_GPU_REWRITE.md",
                       "experiments/native/zipper_rewrite_cpu.cpp",
                       "native/zipper_verifier/main.cpp", "native/zipper_verifier/core.hpp",
                       "native/zipper_verifier/work_size.hpp",
                       "native/zipper_verifier/typed_state.hpp",
                       "native/zipper_verifier/typed_algorithm.cpp",
                       "tests/native/frozen_verifier_9c477dd.tar.gz",
                       "tests/native/frozen_verifier_9c477dd.json",
                       "docs/zipper-typed-verifier.md",
                       "native/zipper_verifier/vendor/nlohmann/json.hpp",
                       "native/zipper_verifier/vendor/nlohmann/LICENSE.MIT",
                       "src/zipper_algebra/zipping_algorithm/data/full-catalog.json",
                       "src/zipper_algebra/zipping_algorithm/data/schema.json",
                       "src/zipper_algebra/zipping_algorithm/data/manifest.json"):
            assert needed in relative, needed
        root = next(iter(roots))
        def read(name):
            member = archive.extractfile(name)
            assert member is not None, name
            return member.read()
        # Bind the separately bundled oracle source to this checkout's curated
        # inventory, as well as the archive's self-contained fixture manifest.
        inventory = json.loads((Path(__file__).resolve().parents[1]
                                / "docs/x3-reproduction-files.json").read_text())
        indexed = {item["path"]: item for item in inventory["files"]}
        for name in ("tests/native/frozen_verifier_9c477dd.json",
                     "tests/native/frozen_verifier_9c477dd.tar.gz"):
            assert name in indexed, name
            content = read(f"{root}/{name}")
            assert len(content) == indexed[name]["bytes"], name
            assert hashlib.sha256(content).hexdigest() == indexed[name]["sha256"], name
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
    forbidden = {".git", ".cache", "__pycache__", "private-baselines", "research-records",
                 "recovery", "manuscripts", "historical-reconstruction", "logs_and_traces",
                 "logs", "ms_inputs", "archive", "gap_scripts",
                 "legacy_x3", "zipper_in_ms.py",
                 "research-archive.json", "restore_data.py", "reconstruction_theory.tex",
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
        # The maintained subpackage shares the historical experiment's name.
        # Keep the retired top-level tree outside both release formats.
        if "zipping_algorithm" in path.parts:
            assert (path.parts[:2] == ("zipper_algebra", "zipping_algorithm")
                    or path.parts[:4] == (root, "src", "zipper_algebra", "zipping_algorithm")), name
        assert path.suffix not in {".pyc", ".pyo", ".pkl", ".pickle"}, name
    print(f"Distribution contents passed: {wheel.name}, {source.name}")


if __name__ == "__main__":
    inspect(sys.argv[1] if len(sys.argv) > 1 else "dist")
