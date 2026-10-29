"""CLI reports and exit statuses for small independently checkable identities."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import runpy
import sys
import pytest
from unittest.mock import patch

import zipper_algebra
from zipper_algebra import AtomSet, Cofactor, Presentation
from zipper_algebra import cli


class TestCommandLineRequirements:
    def invoke(self, arguments):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            code = cli.main(arguments)
        return code, stdout.getvalue(), stderr.getvalue()

    def certificate(self, directory, *, coefficient=1, target="X_2"):
        # ZERO_2_1 is precisely X_2, so this one-term identity is transparent.
        presentation = Presentation(2, zero_from=2)
        certificate = Cofactor.from_presentation(
            presentation, target, AtomSet([]),
            [((), "ZERO_2_1", (), {(0, 0): coefficient}, ())], "test_cli",
        )
        path = Path(directory) / "certificate.txt"
        certificate.save(str(path))
        return path

    @staticmethod
    def verify_arguments(path, *extra):
        return ["verify", str(path), "--n", "2", "--zero-from", "2", *extra]

    def test_describe_text_lists_relations_in_presentation_order(self):
        code, output, errors = self.invoke(["describe", "--n", "2"])
        assert (code) == (0)
        assert (errors) == ('')
        assert ('n2-bigelow') in (output)
        assert ('3 relations') in (output)
        assert ([line.split()[-1] for line in output.splitlines()[1:]]) == (['INV_1_p', 'INV_1_m', 'T_2'])

    @pytest.mark.parametrize('value, message', [('two', 'expected an integer'), ('1', 'expected an integer at least 2')])
    def test_integer_arguments_explain_invalid_and_out_of_range_values(self, value, message):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            with pytest.raises(SystemExit) as error:
                cli.main(["describe", "--n", value])
        assert (error.value.code) == (2)
        assert (message) in (stderr.getvalue())
        assert (stdout.getvalue()) == ('')

    @pytest.mark.parametrize('option', ['--max-zipper-level', '--max-untwist-level', '--zero-from'])
    def test_level_limits_cannot_exceed_the_supplied_strand_count(self, option):
        code, output, errors = self.invoke(
            ["describe", "--n", "2", option, "3", "--json"])
        assert (code) == (2)
        assert (errors) == ('')
        report = json.loads(output)
        assert (report['status']) == ('invalid-input')
        assert (option + ' cannot exceed --n') in (report['error'])

    @pytest.mark.parametrize('extra', [pytest.param([], id='saved-target'),
                                      pytest.param(['--target', 'X_2'], id='explicit-target')])
    def test_saved_target_and_explicit_target_verify_the_same_small_identity(self, extra, tmp_path):
        directory = tmp_path
        path = self.certificate(directory)
        code, output, errors = self.invoke(
            self.verify_arguments(path, "--json", "--pit", *extra))
        assert (code) == (0)
        assert (errors) == ('')
        report = json.loads(output)
        assert report['pass']
        assert report['exact']['pass']
        assert report['pit']['pass']
        assert (report['exact']['tuples']) == (1)

    def test_text_verification_displays_exact_identity_result_without_grid(self, tmp_path):
        directory = tmp_path
        path = self.certificate(directory)
        code, output, errors = self.invoke(self.verify_arguments(path))
        assert (code) == (0)
        assert (errors) == ('')
        assert ('PASS:') in (output)
        assert ('FREE algebra') in (output)

    def test_wrong_coefficient_is_a_false_identity_with_both_checks(self, tmp_path):
        directory = tmp_path
        path = self.certificate(directory, coefficient=2)
        code, output, errors = self.invoke(
            self.verify_arguments(path, "--json", "--pit"))
        assert (code) == (1)
        assert (errors) == ('')
        report = json.loads(output)
        assert not (report['pass'])
        assert not (report['exact']['pass'])
        assert not (report['pit']['pass'])

    def test_missing_target_requires_explicit_target_and_then_succeeds(self, tmp_path):
        directory = tmp_path
        path = self.certificate(directory, target="")
        code, output, errors = self.invoke(self.verify_arguments(path, "--json"))
        assert (code) == (2)
        assert (errors) == ('')
        report = json.loads(output)
        assert (report['status']) == ('invalid-input')
        assert ('supply --target') in (report['error'])
        code, output, errors = self.invoke(
            self.verify_arguments(path, "--target", "X_2", "--json"))
        assert (code) == (0)
        assert json.loads(output)['pass']

    def test_file_errors_go_to_stderr_for_text_mode(self, tmp_path):
        directory = tmp_path
        path = Path(directory) / "absent.txt"
        code, output, errors = self.invoke(self.verify_arguments(path))
        assert (code) == (2)
        assert (output) == ('')
        assert ('zipper-algebra:') in (errors)
        assert ('absent.txt') in (errors)

    def test_json_description_defaults_to_four_strands(self):
        code, output, errors = self.invoke(["describe", "--json"])
        assert (code) == (0)
        assert (errors) == ('')
        report = json.loads(output)
        assert (report['n']) == (4)
        assert (report['variables']) == (['x1', 'x2', 'x3', 'y1', 'y2', 'y3'])
        assert (len(report['relations'])) == (31)

    def test_module_entrypoint_uses_the_same_cli(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(sys, "argv", ["zipper-algebra", "describe", "--n", "2", "--json"]):
            with redirect_stdout(stdout), redirect_stderr(stderr):
                with pytest.raises(SystemExit) as error:
                    runpy.run_module("zipper_algebra", run_name="__main__")
        assert (error.value.code) == (0)
        assert (json.loads(stdout.getvalue())['n']) == (2)
        assert (stderr.getvalue()) == ('')


class TestPackageAttributeRequirements:
    def test_unknown_public_attribute_has_a_normal_attribute_error(self):
        with pytest.raises(AttributeError, match='zipper_algebra.*not_a_public_name'):
            getattr(zipper_algebra, "not_a_public_name")

    def test_lazy_public_class_is_cached_with_its_ordinary_identity(self):
        cached = zipper_algebra.__dict__.pop("Presentation", None)
        try:
            resolved = getattr(zipper_algebra, "Presentation")
            assert (resolved) is (Presentation)
            assert (zipper_algebra.__dict__['Presentation']) is (Presentation)
        finally:
            if cached is not None:
                zipper_algebra.__dict__["Presentation"] = cached
