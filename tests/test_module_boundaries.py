"""Public arithmetic imports and the separated solver-input entry point."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from zipper_algebra import coefficients, msolve, polynomial, relations, zipper_in_ms


class TestEmitterModuleBoundaries(unittest.TestCase):
    def test_legacy_imports_refer_to_the_same_classes_and_functions(self):
        for module, names in (
            (coefficients, ("Coef", "QTCoef", "NumericRing", "SymbolicRing")),
            (polynomial, ("Poly",)),
            (relations, ("RelTag", "build_relations", "build_relations_ring",
                         "append_zero_sliders", "append_target")),
            (msolve, ("main",)),
        ):
            for name in names:
                with self.subTest(name=name):
                    self.assertIs(getattr(zipper_in_ms, name), getattr(module, name))

    def test_separated_arithmetic_and_relations_share_their_coefficient_parent(self):
        names, defining_relations, sliders, tags = relations.build_relations_ring(
            2, coefficients.NumericRing(5, 7, 0), "bigelow")
        self.assertEqual(names, ["x1", "y1"])
        self.assertEqual([tag.label(2) for tag in tags], ["INV_1_p", "INV_1_m", "T_2"])
        relation = defining_relations[-1]
        self.assertIsInstance(relation, polynomial.Poly)
        # (x-7)(5y-4-x) = -x*x + 5*x*y + 3*x - 35*y + 28.
        self.assertEqual({word: coefficient.v for word, coefficient in relation.terms.items()},
                         {(1, 1): -1, (1, 2): 5, (1,): 3, (2,): -35, (): 28})
        self.assertEqual({word: coefficient.v for word, coefficient in sliders[(2, 1)].terms.items()},
                         {(1,): -1, (2,): 5, (): -4})

    def test_new_msolve_module_emits_the_exact_two_generator_input(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "input.ms"
            result = subprocess.run(
                [sys.executable, "-I", "-m", "zipper_algebra.msolve", "--n", "2",
                 "--untwist", "none", "--out", str(output)],
                cwd=directory, capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            expected = "x1,y1\n0\nx1*y1-1,\ny1*x1-1\n"
            self.assertEqual(output.read_bytes(), expected.encode("ascii"))
            self.assertIn("Wrote 2 relations in 2 variables", result.stdout)

    def test_new_msolve_module_rejects_the_removed_membership_shortcut(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "input.ms"
            result = subprocess.run(
                [sys.executable, "-I", "-m", "zipper_algebra.msolve", "--n", "3",
                 "--target", "3", "--out", str(output)],
                cwd=directory, capture_output=True, text=True, check=False,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("does not prove membership in the original ideal", result.stderr)
            self.assertFalse(output.exists())
