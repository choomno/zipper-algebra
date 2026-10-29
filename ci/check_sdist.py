"""Rebuild and test the source distribution in an independent directory."""
from __future__ import annotations

from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import zipfile

from check_distribution import single


def check(directory, no_isolation=False):
    source = single(directory, "*.tar.gz")
    with tempfile.TemporaryDirectory(prefix="zipper-sdist-") as directory:
        # A source archive may be rebuilt beneath another checkout. Its origin
        # must stay tied to the archived source, not to that unrelated parent.
        workspace = Path(directory)
        subprocess.run(["git", "init", "-q", str(workspace)], check=True)
        (workspace / "unrelated.txt").write_text("Unrelated source archive build.\n")
        subprocess.run(["git", "-C", str(workspace), "add", "unrelated.txt"], check=True)
        subprocess.run(["git", "-C", str(workspace), "-c", "user.name=Build fixture",
                        "-c", "user.email=build-fixture@example.invalid", "commit", "-q",
                        "-m", "Unrelated build environment"], check=True)
        extracted = workspace / "source"
        extracted.mkdir()
        with tarfile.open(source) as archive:
            for entry in archive:
                target = Path(extracted, entry.name).resolve()
                assert target.is_relative_to(extracted.resolve()), entry.name
                assert not entry.issym() and not entry.islnk(), entry.name
            archive.extractall(extracted)
        trees = [p for p in extracted.iterdir() if p.is_dir()]
        assert len(trees) == 1
        tree = trees[0]
        origin = json.loads((tree / "src/zipper_algebra/_build_provenance.json").read_text())
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
        with zipfile.ZipFile(wheel) as archive:
            assert json.loads(archive.read("zipper_algebra/_build_provenance.json")) == origin, \
                "rebuilding a source archive must preserve its recorded source origin"
        subprocess.run([sys.executable, "-m", "pip", "install", "--force-reinstall",
                        "--no-deps", str(wheel)], check=True)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1")
        env.pop("PYTHONPATH", None)
        subprocess.run([sys.executable, "ci/run_tests.py"],
                       cwd=tree, env=env, check=True)
        subprocess.run([sys.executable, "ci/check_installed.py"],
                       cwd=tree, env=env, check=True)
    print("Source-distribution build and regressions passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="dist")
    parser.add_argument("--no-isolation", action="store_true",
                        help="use already-installed build tools for an offline local check")
    args = parser.parse_args()
    check(args.directory, no_isolation=args.no_isolation)
