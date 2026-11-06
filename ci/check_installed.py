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
             "import pathlib, zipper_algebra; "
             "from zipper_algebra.invariants import minimal_polynomial; "
             "from zipper_algebra.production import reconstruct_basis, reconstruct_groebner; "
             "from zipper_algebra.reconstruction import reconstruct_kernel_candidate; "
             "from zipper_algebra.cofactor.general import "
             "cofactor_representation, expand_derivation, verify_cofactor; "
             "from zipper_algebra.action_production import word_basis_from_rewrites; "
             "from zipper_algebra.execution.production import produce_algebra; "
             "from zipper_algebra.verification.krylov import MinimalPolynomialCertificate; "
             "from zipper_algebra.zipping_algorithm import Catalog; "
             "from zipper_algebra.zipping_algorithm.two_strand import reduce_two_strand, check_two_strand; "
             "assert len(Catalog.bundled().rules) == 128; "
             "assert check_two_strand(reduce_two_strand([1, -1, 1]))['pass']; "
              "from zipper_algebra.cofactor.reconstruction import reconstruct_derivation; "
             "import sys; "
             "assert not {'sympy', 'fflow', 'sage'}.intersection(sys.modules); "
             "print(pathlib.Path(zipper_algebra.__file__).resolve())"],
            cwd=directory, text=True, capture_output=True, check=True,
        )
        package = Path(probe.stdout.strip())
        assert not package.is_relative_to(root), "this check requires an installed wheel"
        subprocess.run(
            [sys.executable, "-I", "-c",
             "from zipper_algebra.zipping_algorithm.native import NativeSession; "
             "session = NativeSession(); "
             "result = session.reduce([1, -2, 1], n=3, trace=True, batch_size=4); "
             "assert result['status'] == 'complete'; "
             "assert session.check(result); session.close()"],
            cwd=directory, check=True,
        )
        report = run("describe", "--n", "3", "--json", cwd=directory)
        assert report["n"] == 3 and "Z_2" in report["relations"]
        certificate = Path(directory) / "certificate.txt"
        shutil.copyfile(root / "examples/x3/certificate.txt",
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
