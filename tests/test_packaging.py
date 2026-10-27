"""Exercise the installed public interface, CLI and legacy import aliases."""
from __future__ import annotations

import importlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import zipper_algebra
from zipper_algebra import AtomSet, Cofactor, Presentation


class TestInstalledInterface(unittest.TestCase):
    def test_legacy_imports_alias_the_maintained_modules(self):
        for suffix in ("core", "zipper_in_ms", "presentation", "cofactor.atoms",
                       "cofactor.cofactor_io", "cofactor.compose", "cofactor.lift",
                       "cofactor.ncpoly", "cofactor.proof", "cofactor.skeleton",
                       "cofactor.verify"):
            with self.subTest(module=suffix):
                self.assertIs(importlib.import_module(f"tools.{suffix}"),
                              importlib.import_module(f"zipper_algebra.{suffix}"))

    def test_version_and_runtime_dependencies_match_distribution_metadata(self):
        self.assertEqual(zipper_algebra.__version__, importlib.metadata.version("zipper-algebra"))
        requirements = importlib.metadata.requires("zipper-algebra") or []
        self.assertEqual([requirement for requirement in requirements if ";" not in requirement], [])

    def test_public_imports_need_no_workspace_configuration_or_optional_inference(self):
        code = """
import json, sys
before = list(sys.path)
from zipper_algebra import Presentation, AtomSet, Cofactor
from zipper_algebra.cofactor import compose, lift, verify
assert sys.path == before
assert 'repo_settings' not in sys.modules
assert 'sympy' not in sys.modules
print(json.dumps({'relations': len(Presentation(4).labels())}))
"""
        result = subprocess.run([sys.executable, "-I", "-c", code], cwd=tempfile.gettempdir(),
                                capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(result.stdout), {"relations": 31})


class TestCommandLine(unittest.TestCase):
    def run_cli(self, *arguments):
        return subprocess.run([sys.executable, "-I", "-m", "zipper_algebra", *arguments],
                              cwd=tempfile.gettempdir(), capture_output=True, text=True)

    def test_describe_outputs_the_public_presentation_identity(self):
        result = self.run_cli("describe", "--n", "4", "--json")
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads(result.stdout)
        self.assertEqual(report["presentation_hash"], Presentation(4).content_hash())
        self.assertEqual(report["relations"], Presentation(4).labels())

    def test_verify_distinguishes_true_false_and_invalid_input(self):
        presentation = Presentation(4, extra=("X_4",))
        with tempfile.TemporaryDirectory(prefix="zipper-cli-") as directory:
            path = Path(directory) / "certificate.txt"
            for coefficient, expected_status in ((1, 0), (2, 1)):
                certificate = Cofactor.from_presentation(
                    presentation, "EXTRA[X_4]", AtomSet([]),
                    [((), "EXTRA[X_4]", (), {(0, 0): coefficient}, ())], "test_packaging",
                )
                certificate.save(str(path))
                result = self.run_cli("verify", str(path), "--n", "4", "--extra", "X_4",
                                      "--pit", "--json")
                self.assertEqual(result.returncode, expected_status, result.stderr)
                self.assertEqual(json.loads(result.stdout)["pass"], expected_status == 0)
            path.write_text("unparseable certificate\n")
            result = self.run_cli("verify", str(path), "--n", "4", "--json")
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertEqual(json.loads(result.stdout)["status"], "invalid-input")

    def test_verify_requires_an_explicit_algebra_and_matching_presentation(self):
        result = self.run_cli("verify", "missing-certificate.txt", "--json")
        self.assertEqual(result.returncode, 2)
        with tempfile.TemporaryDirectory(prefix="zipper-cli-") as directory:
            path = Path(directory) / "certificate.txt"
            Cofactor.from_presentation(Presentation(4), "X_3", AtomSet([]), [], "test").save(str(path))
            result = self.run_cli("verify", str(path), "--n", "4", "--untwist", "none", "--json")
            self.assertEqual(result.returncode, 2)
            self.assertIn("hash", json.loads(result.stdout)["error"])

    def test_invalid_presentation_options_fail(self):
        for arguments in (("--n", "1"), ("--untwist", "unrecognised"),
                          ("--max-untwist-level", "1"), ("--zero-from", "8")):
            with self.subTest(arguments=arguments):
                self.assertEqual(self.run_cli("describe", *arguments).returncode, 2)
