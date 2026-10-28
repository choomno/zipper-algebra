"""Small exact identities and observable diagnostics for both certificate checks."""
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction
import io
from pathlib import Path
import runpy
import sys
import tempfile
import unittest
from unittest.mock import patch
import warnings

from zipper_algebra import Presentation
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.cofactor_io import Cofactor
from zipper_algebra.cofactor import verify


ONE = {(0, 0): Fraction(1)}


class TestVerifierExactExamples(unittest.TestCase):
    def certificate(self, pres=None, terms=(), atoms=None, target="X_2"):
        pres = pres or Presentation(2)
        return Cofactor.from_presentation(pres, target, atoms or AtomSet([]), terms,
                                          "exact-verifier-examples")

    def test_untagged_certificates_require_an_explicit_presentation(self):
        cof = self.certificate()
        cof.presentation_hash = ""
        for checker in (verify.check, verify.pit_check):
            with self.subTest(checker=checker.__name__):
                with self.assertRaisesRegex(ValueError, "explicit pres"):
                    checker(cof, target={}, verbose=False)
                self.assertTrue(checker(cof, Presentation(2), {}, verbose=False)["pass"])

    def test_default_presentation_and_stored_target(self):
        pres = Presentation(3)
        # The target is the canonical lasso untwisting relation. Its reverse-word
        # representative differs in the free algebra by explicit inverse terms.
        cof = self.lasso_certificate(pres)
        self.assertTrue(verify.check(cof, verbose=False)["pass"])
        self.assertTrue(verify.pit_check(cof, verbose=False)["pass"])

    @staticmethod
    def lasso_certificate(pres):
        sliders = pres.symbolic_ncpolys()[1]
        slider = sliders[(3, 1)]
        terms = [((1, 2), "T_3", (), ONE, ())]
        # x1*x2*y2*y1 - 1 = x1*(x2*y2 - 1)*y1 + (x1*y1 - 1).
        for word, coefficient in slider.items():
            q_squared = {(i + 2, j): c for (i, j), c in coefficient.items()}
            terms.append(((1,), "INV_2_p", (pres.n,) + word, q_squared, ()))
            terms.append(((), "INV_1_p", word, q_squared, ()))
        return Cofactor.from_presentation(pres, "(Theta_3-q^2)X_3", AtomSet([]),
                                          terms, "lasso-unit-transport")

    def test_many_small_terms_produce_progress_and_correct_identity(self):
        # 100 copies of one hundredth of x1*y1-1 reproduce that relator exactly.
        pres = Presentation(2)
        target = {(1, 2): ONE, (): {(0, 0): Fraction(-1)}}
        terms = [((), "INV_1_p", (), {(0, 0): Fraction(1, 100)}, ())] * 100
        terms.append(((), "INV_1_p", (), {}, ()))
        cof = self.certificate(pres, terms)
        output = io.StringIO()
        with redirect_stdout(output):
            result = verify.check(cof, pres, target, verbose=True)
        self.assertTrue(result["pass"])
        self.assertEqual(result["tuples"], 101)
        self.assertEqual(result["labels_used"], ["INV_1_p"])
        self.assertEqual(result["residual_monomials"], 0)
        self.assertIn("[acc]   100/101", output.getvalue())
        self.assertIn("identity verified in the FREE algebra", output.getvalue())

    def test_failure_diagnostics_include_scalar_and_nonempty_word(self):
        cof = self.certificate()
        target = {(): ONE, (1,): {(0, 0): Fraction(3)}}
        output = io.StringIO()
        with redirect_stdout(output):
            result = verify.check(cof, Presentation(2), target, verbose=True)
        self.assertFalse(result["pass"])
        self.assertEqual(result["residual_monomials"], 2)
        self.assertIn("FAIL: residual has 2", output.getvalue())
        self.assertIn("1: 1 q,t-terms", output.getvalue())
        self.assertIn("x1: 1 q,t-terms", output.getvalue())

    def test_fractional_identity_remains_defined_on_clearing_factor_zero_set(self):
        pres = Presentation(2, "none")
        atoms = AtomSet([Atom("q-t", {(1, 0): 1, (0, 1): -1})])
        terms = [((), "INV_1_p", (), {(1, 0): 1, (0, 1): -1}, (1,))]
        cof = self.certificate(pres, terms, atoms)
        target = {(1, 2): ONE, (): {(0, 0): Fraction(-1)}}
        output = io.StringIO()
        with redirect_stdout(output):
            exact = verify.check(cof, pres, target, verbose=True)
            grid = verify.pit_check(cof, pres, target, verbose=True)
        self.assertTrue(exact["pass"])
        self.assertTrue(grid["pass"])
        self.assertEqual(grid["failures"], [])
        self.assertIn("exact by interpolation", output.getvalue())

    def test_grid_failure_report_is_bounded_but_covers_entire_grid(self):
        pres = Presentation(2, "none")
        cof = self.certificate(pres)
        # q**11+1 never vanishes on the positive grid, yielding 12 failures.
        target = {(): {(11, 0): 1, (0, 0): 1}}
        output = io.StringIO()
        with redirect_stdout(output):
            result = verify.pit_check(cof, pres, target, verbose=True)
        self.assertFalse(result["pass"])
        self.assertEqual(result["bidegree"], [11, 0])
        self.assertEqual(result["grid"], [12, 1])
        self.assertEqual(len(result["failures"]), 10)
        self.assertEqual(result["failures"][0], (2, 2, 1))
        self.assertEqual(result["failures"][-1], (11, 2, 1))
        self.assertIn("FAIL at 12 point(s)", output.getvalue())

    def test_compatibility_target_resolver_retains_word_order(self):
        result = verify._symbolic_target(Presentation(2), {}, "X_2")
        self.assertEqual(result, {(2,): {(1, 0): Fraction(1)},
                                  (): {(0, 0): Fraction(1), (1, 0): Fraction(-1)},
                                  (1,): {(0, 0): Fraction(-1)}})

    def test_atom_powers_and_grid_have_independent_expected_values(self):
        atoms = AtomSet([Atom("q-2", {(1, 0): 1, (0, 0): -2})])
        powers = verify._atom_powers(atoms, (3,))[0]
        self.assertEqual(powers[0], ONE)
        self.assertEqual(powers[3], {(3, 0): 1, (2, 0): -6,
                                   (1, 0): 12, (0, 0): -8})
        self.assertEqual(verify._grid_points(0), [])
        self.assertEqual(verify._grid_points(4), [2, 3, 4, 5])

    def invoke_legacy_cli(self, path, arguments):
        output = io.StringIO()
        with patch.object(sys, "argv", ["verify", str(path), *arguments]), \
             warnings.catch_warnings(), redirect_stdout(output), redirect_stderr(output):
            warnings.simplefilter("ignore", RuntimeWarning)
            with self.assertRaises(SystemExit) as caught:
                runpy.run_module("zipper_algebra.cofactor.verify", run_name="__main__")
        return caught.exception.code, output.getvalue()

    def test_module_entrypoint_verifies_free_algebra_lasso_certificate(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cofactor.txt"
            cof = self.lasso_certificate(Presentation(4))
            cof.save(path)
            status, output = self.invoke_legacy_cli(path, ["--pit"])
            self.assertEqual(status, 0)
            self.assertIn("PASS: identity holds on the full grid", output)
            status, output = self.invoke_legacy_cli(path, ["--n", "4"])
            self.assertEqual(status, 0)
            self.assertIn("identity verified in the FREE algebra", output)

    def test_module_entrypoint_false_identity_and_explicit_target(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cofactor.txt"
            self.certificate(Presentation(2)).save(path)
            status, output = self.invoke_legacy_cli(path, ["--n", "2", "--pit", "--target", "X_2"])
            self.assertEqual(status, 1)
            self.assertIn("FAIL: residual", output)

    def test_module_entrypoint_requests_missing_target(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cofactor.txt"
            self.certificate(Presentation(2), target="").save(path)
            status, _output = self.invoke_legacy_cli(path, ["--n", "2"])
            self.assertIn("carries no target", status)


if __name__ == "__main__":
    unittest.main()
