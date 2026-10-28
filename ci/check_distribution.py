"""Inspect release archives before installation; needs only the standard library."""
from __future__ import annotations

from pathlib import Path, PurePosixPath
import sys
import tarfile
import zipfile


def single(directory, pattern):
    found = list(Path(directory).glob(pattern))
    if len(found) != 1:
        raise SystemExit(f"expected one {pattern} in {directory}, found {len(found)}")
    return found[0]


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
        for needed in ("__init__.py", "cli.py", "presentation.py", "core.py",
                       "coefficients.py", "polynomial.py", "relations.py", "msolve.py",
                       "linear.py", "cofactor/verify.py", "cofactor/lift.py",
                       "cofactor/inference.py", "cofactor/univariate.py",
                       "cofactor/artifacts.py", "execution/config.py",
                       "execution/runner.py", "execution/artifacts.py",
                       "execution/sampling.py", "execution/pipeline.py"):
            assert f"zipper_algebra/{needed}" in names, needed
        assert metadata + "licenses/LICENSE" in names, "the wheel must include the license"
        assert b"License-Expression: MIT" in archive.read(metadata + "METADATA")
    with tarfile.open(source) as archive:
        names = archive.getnames()
        roots = {PurePosixPath(name).parts[0] for name in names}
        assert len(roots) == 1, "source distribution must have one enclosing directory"
        relative = {str(PurePosixPath(*PurePosixPath(name).parts[1:])) for name in names}
        for needed in ("pyproject.toml", "README.md", "docs/TESTING.md",
                       "docs/X3_REPRODUCTION.md", "docs/x3-reproduction-files.json",
                       "patches/f4ncgb/msolve-parser-fix.patch", "patches/f4ncgb/LICENSE",
                       "docs/mathematical-specification/main.tex",
                       "docs/mathematical-specification/mathematical-specification.pdf",
                       "tests/reproduction/test_x3.py",
                       "tests/reproduction/test_legacy_proofs.py",
                       "tests/reference/legacy_x3.py", "examples/infer_denominators.py",
                       "experiments/survey/pointalg.py",
                       "experiments/research/theta_min_poly.py",
                       "experiments/variants/square_zero_regular_rep_qt.py",
                       "ci/run_tests.py", "ci/build_x3_manifest.py"):
            assert needed in relative, needed
    forbidden = {".git", ".cache", "__pycache__", "private-baselines", "research-records",
                 "recovery", "manuscripts", "historical-reconstruction",
                 ".claude", ".codex", ".agents", ".agent-notes", "agent-notes",
                 "CONTRIBUTING.md", "CONTRIBUTING.MD", "AGENTS.md", "CLAUDE.md"}
    for name in wheel_names + names:
        path = PurePosixPath(name)
        assert not path.is_absolute() and ".." not in path.parts and "\\" not in name, name
        assert not forbidden.intersection(path.parts), name
        assert path.suffix not in {".pyc", ".pyo", ".pkl", ".pickle"}, name
    print(f"Distribution contents passed: {wheel.name}, {source.name}")


if __name__ == "__main__":
    inspect(sys.argv[1] if len(sys.argv) > 1 else "dist")
