"""Exercise the installed public interface and CLI."""
from __future__ import annotations

import importlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import pytest

import zipper_algebra
from zipper_algebra import AtomSet, Cofactor, Presentation


class TestInstalledInterface:
    def test_version_and_runtime_dependencies_match_distribution_metadata(self):
        assert (zipper_algebra.__version__) == (importlib.metadata.version('zipper-algebra'))
        requirements = importlib.metadata.requires("zipper-algebra") or []
        assert ([requirement for requirement in requirements if ';' not in requirement]) == ([])

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
        assert (json.loads(result.stdout)) == ({'relations': 31})


class TestCommandLine:
    def run_cli(self, *arguments):
        return subprocess.run([sys.executable, "-I", "-m", "zipper_algebra", *arguments],
                              cwd=tempfile.gettempdir(), capture_output=True, text=True)

    def test_describe_outputs_the_public_presentation_identity(self):
        result = self.run_cli("describe", "--n", "4", "--json")
        assert (result.returncode) == (0), result.stderr
        report = json.loads(result.stdout)
        assert (report['presentation_hash']) == (Presentation(4).content_hash())
        assert (report['relations']) == (Presentation(4).labels())

    def test_verify_distinguishes_true_false_and_invalid_input(self, tmp_path):
        presentation = Presentation(4, extra=("X_4",))
        directory = tmp_path
        path = Path(directory) / "certificate.txt"
        for coefficient, expected_status in ((1, 0), (2, 1)):
            certificate = Cofactor.from_presentation(
                presentation, "EXTRA[X_4]", AtomSet([]),
                [((), "EXTRA[X_4]", (), {(0, 0): coefficient}, ())], "test_packaging",
            )
            certificate.save(str(path))
            result = self.run_cli("verify", str(path), "--n", "4", "--extra", "X_4",
                                  "--pit", "--json")
            assert (result.returncode) == (expected_status), result.stderr
            assert (json.loads(result.stdout)['pass']) == (expected_status == 0)
        path.write_text("unparseable certificate\n")
        result = self.run_cli("verify", str(path), "--n", "4", "--json")
        assert (result.returncode) == (2), result.stderr
        assert (json.loads(result.stdout)['status']) == ('invalid-input')

    def test_verify_requires_an_explicit_algebra_and_matching_presentation(self, tmp_path):
        result = self.run_cli("verify", "missing-certificate.txt", "--json")
        assert (result.returncode) == (2)
        directory = tmp_path
        path = Path(directory) / "certificate.txt"
        Cofactor.from_presentation(Presentation(4), "X_3", AtomSet([]), [], "test").save(str(path))
        result = self.run_cli("verify", str(path), "--n", "4", "--untwist", "none", "--json")
        assert (result.returncode) == (2)
        assert ('hash') in (json.loads(result.stdout)['error'])

    @pytest.mark.parametrize('arguments', [('--n', '1'), ('--untwist', 'unrecognised'), ('--max-untwist-level', '1'), ('--zero-from', '8')])
    def test_invalid_presentation_options_fail(self, arguments):
        assert (self.run_cli('describe', *arguments).returncode) == (2)
