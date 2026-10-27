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
        command = [sys.executable, "-m", "build", "--wheel"]
        if no_isolation:
            command.append("--no-isolation")
        subprocess.run(command, cwd=tree, check=True)
        wheel = single(tree / "dist", "*.whl")
        subprocess.run([sys.executable, "-m", "pip", "install", "--force-reinstall",
                        "--no-deps", str(wheel)], check=True)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1")
        env.pop("PYTHONPATH", None)
        for args in (["discover", "-s", "tests", "-v"],
                     ["tools.cofactor.tests_proof_skeleton", "-v"]):
            subprocess.run([sys.executable, "-m", "unittest", *args], cwd=tree, env=env,
                           check=True)
    print("Source-distribution build and regressions passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="dist")
    parser.add_argument("--no-isolation", action="store_true",
                        help="use already-installed build tools for an offline local check")
    args = parser.parse_args()
    check(args.directory, no_isolation=args.no_isolation)
