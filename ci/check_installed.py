"""Check wheel imports and command behaviour from a fresh working directory."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def run(*args, cwd, expected=0):
    completed = subprocess.run([sys.executable, "-I", "-m", "zipper_algebra", *args],
                               cwd=cwd, text=True, capture_output=True)
    if completed.returncode != expected:
        raise AssertionError(f"{args}: exit {completed.returncode}\n"
                             f"{completed.stdout}\n{completed.stderr}")
    return json.loads(completed.stdout)


def check():
    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="zipper-installed-") as directory:
        probe = subprocess.run(
            [sys.executable, "-I", "-c",
             "import pathlib, zipper_algebra; print(pathlib.Path(zipper_algebra.__file__).resolve())"],
            cwd=directory, text=True, capture_output=True, check=True,
        )
        package = Path(probe.stdout.strip())
        assert not package.is_relative_to(root), "this check requires an installed wheel"
        report = run("describe", "--n", "3", "--json", cwd=directory)
        assert report["n"] == 3 and "Z_2" in report["relations"]
        certificate = Path(directory) / "historical.txt"
        shutil.copyfile(root / "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt",
                        certificate)
        report = run("verify", str(certificate), "--n", "4", "--target", "X_3", "--json",
                     cwd=directory)
        assert report["pass"] is True
        # An empty, well-formed certificate cannot express the nonzero free polynomial X3.
        false = Path(directory) / "false.txt"
        false.write_text("# target: X_3\n", encoding="utf-8")
        report = run("verify", str(false), "--n", "4", "--target", "X_3", "--json",
                     cwd=directory, expected=1)
        assert report["pass"] is False
        false.write_text("this is not a certificate\n", encoding="utf-8")
        report = run("verify", str(false), "--n", "4", "--json", cwd=directory, expected=2)
        assert report["status"] == "invalid-input"
    print("Installed-wheel imports and CLI checks passed outside the checkout")


if __name__ == "__main__":
    check()
