"""CLI reports and exit statuses for small independently checkable identities."""
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import runpy
import sys
import tempfile
import unittest
from unittest.mock import patch

import zipper_algebra
from zipper_algebra import AtomSet, Cofactor, Presentation
from zipper_algebra import cli


class TestCommandLineRequirements(unittest.TestCase):
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
        self.assertEqual(code, 0)
        self.assertEqual(errors, "")
        self.assertIn("n2-bigelow", output)
        self.assertIn("3 relations", output)
        self.assertEqual([line.split()[-1] for line in output.splitlines()[1:]],
                         ["INV_1_p", "INV_1_m", "T_2"])

    def test_integer_arguments_explain_invalid_and_out_of_range_values(self):
        for value, message in (("two", "expected an integer"),
                               ("1", "expected an integer at least 2")):
            with self.subTest(value=value):
                stdout, stderr = io.StringIO(), io.StringIO()
                with redirect_stdout(stdout), redirect_stderr(stderr):
                    with self.assertRaises(SystemExit) as error:
                        cli.main(["describe", "--n", value])
                self.assertEqual(error.exception.code, 2)
                self.assertIn(message, stderr.getvalue())
                self.assertEqual(stdout.getvalue(), "")

    def test_level_limits_cannot_exceed_the_supplied_strand_count(self):
        for option in ("--max-zipper-level", "--max-untwist-level", "--zero-from"):
            with self.subTest(option=option):
                code, output, errors = self.invoke(
                    ["describe", "--n", "2", option, "3", "--json"])
                self.assertEqual(code, 2)
                self.assertEqual(errors, "")
                report = json.loads(output)
                self.assertEqual(report["status"], "invalid-input")
                self.assertIn(option + " cannot exceed --n", report["error"])

    def test_saved_target_and_explicit_target_verify_the_same_small_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.certificate(directory)
            for extra in ([], ["--target", "X_2"]):
                with self.subTest(extra=extra):
                    code, output, errors = self.invoke(
                        self.verify_arguments(path, "--json", "--pit", *extra))
                    self.assertEqual(code, 0)
                    self.assertEqual(errors, "")
                    report = json.loads(output)
                    self.assertTrue(report["pass"])
                    self.assertTrue(report["exact"]["pass"])
                    self.assertTrue(report["pit"]["pass"])
                    self.assertEqual(report["exact"]["tuples"], 1)

    def test_text_verification_displays_exact_identity_result_without_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.certificate(directory)
            code, output, errors = self.invoke(self.verify_arguments(path))
            self.assertEqual(code, 0)
            self.assertEqual(errors, "")
            self.assertIn("PASS:", output)
            self.assertIn("FREE algebra", output)

    def test_wrong_coefficient_is_a_false_identity_with_both_checks(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.certificate(directory, coefficient=2)
            code, output, errors = self.invoke(
                self.verify_arguments(path, "--json", "--pit"))
            self.assertEqual(code, 1)
            self.assertEqual(errors, "")
            report = json.loads(output)
            self.assertFalse(report["pass"])
            self.assertFalse(report["exact"]["pass"])
            self.assertFalse(report["pit"]["pass"])

    def test_missing_target_requires_explicit_target_and_then_succeeds(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.certificate(directory, target="")
            code, output, errors = self.invoke(self.verify_arguments(path, "--json"))
            self.assertEqual(code, 2)
            self.assertEqual(errors, "")
            report = json.loads(output)
            self.assertEqual(report["status"], "invalid-input")
            self.assertIn("supply --target", report["error"])
            code, output, errors = self.invoke(
                self.verify_arguments(path, "--target", "X_2", "--json"))
            self.assertEqual(code, 0)
            self.assertTrue(json.loads(output)["pass"])

    def test_file_errors_go_to_stderr_for_text_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "absent.txt"
            code, output, errors = self.invoke(self.verify_arguments(path))
            self.assertEqual(code, 2)
            self.assertEqual(output, "")
            self.assertIn("zipper-algebra:", errors)
            self.assertIn("absent.txt", errors)

    def test_json_description_defaults_to_four_strands(self):
        code, output, errors = self.invoke(["describe", "--json"])
        self.assertEqual(code, 0)
        self.assertEqual(errors, "")
        report = json.loads(output)
        self.assertEqual(report["n"], 4)
        self.assertEqual(report["variables"], ["x1", "x2", "x3", "y1", "y2", "y3"])
        self.assertEqual(len(report["relations"]), 31)

    def test_module_entrypoint_uses_the_same_cli(self):
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(sys, "argv", ["zipper-algebra", "describe", "--n", "2", "--json"]):
            with redirect_stdout(stdout), redirect_stderr(stderr):
                with self.assertRaises(SystemExit) as error:
                    runpy.run_module("zipper_algebra", run_name="__main__")
        self.assertEqual(error.exception.code, 0)
        self.assertEqual(json.loads(stdout.getvalue())["n"], 2)
        self.assertEqual(stderr.getvalue(), "")


class TestPackageAttributeRequirements(unittest.TestCase):
    def test_unknown_public_attribute_has_a_normal_attribute_error(self):
        with self.assertRaisesRegex(AttributeError, "zipper_algebra.*not_a_public_name"):
            getattr(zipper_algebra, "not_a_public_name")

    def test_lazy_public_class_is_cached_with_its_ordinary_identity(self):
        cached = zipper_algebra.__dict__.pop("Presentation", None)
        try:
            resolved = getattr(zipper_algebra, "Presentation")
            self.assertIs(resolved, Presentation)
            self.assertIs(zipper_algebra.__dict__["Presentation"], Presentation)
        finally:
            if cached is not None:
                zipper_algebra.__dict__["Presentation"] = cached
