"""Rebuild and test the source distribution in an independent directory."""
from __future__ import annotations

from pathlib import Path
import argparse
import os
import subprocess
import sys
import tarfile
import tempfile

from check_distribution import single


def check(directory, no_isolation=False):
    source = single(directory, "*.tar.gz")
    with tempfile.TemporaryDirectory(prefix="zipper-sdist-") as directory:
        with tarfile.open(source) as archive:
            for entry in archive:
                target = Path(directory, entry.name).resolve()
                assert target.is_relative_to(Path(directory).resolve()), entry.name
                assert not entry.issym() and not entry.islnk(), entry.name
            archive.extractall(directory)
        trees = [p for p in Path(directory).iterdir() if p.is_dir()]
        assert len(trees) == 1
        tree = trees[0]
        for script, arguments in (("ci/check_docs.py", ()),
                                  ("ci/check_x3_manifest.py", ()),
                                  ("ci/build_x3_manifest.py", ("--check",))):
            subprocess.run([sys.executable, script, *arguments], cwd=tree, check=True)
        # Local ignored working files must never leak into a rebuilt release.
        # These are synthetic privacy-regression fixtures, not instructions.
        for name in ("docs/PRIVATE_BASELINE.md", "docs/baseline.json", "docs/HISTORY.md",
                     "AGENTS.md", "CONTRIBUTING.md", ".agents/notes.md"):
            private = tree / name
            private.parent.mkdir(parents=True, exist_ok=True)
            private.write_text("Private build regression fixture.\n", encoding="utf-8")
        command = [sys.executable, "-m", "build", "--wheel", "--sdist"]
        if no_isolation:
            command.append("--no-isolation")
        subprocess.run(command, cwd=tree, check=True)
        subprocess.run([sys.executable, "ci/check_distribution.py", "dist"], cwd=tree, check=True)
        wheel = single(tree / "dist", "*.whl")
        subprocess.run([sys.executable, "-m", "pip", "install", "--force-reinstall",
                        "--no-deps", str(wheel)], check=True)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1")
        env.pop("PYTHONPATH", None)
        subprocess.run([sys.executable, "ci/run_tests.py", "--suite", "fast"],
                       cwd=tree, env=env, check=True)
        # The full historical comparison suite runs once in the coverage job.
        # Here one complete certificate recomposition proves that its immutable
        # data and reference adapters survive extraction and wheel rebuilding.
        subprocess.run([sys.executable, "-m", "unittest",
                        "tests.reproduction.test_x3.TestCompose."
                        "test_composed_body_is_byte_identical_to_the_shipped_certificate",
                        "-v"], cwd=tree, env=env, check=True)
    print("Source-distribution build and regressions passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="dist")
    parser.add_argument("--no-isolation", action="store_true",
                        help="use already-installed build tools for an offline local check")
    args = parser.parse_args()
    check(args.directory, no_isolation=args.no_isolation)
