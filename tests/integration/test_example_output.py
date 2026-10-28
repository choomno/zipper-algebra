"""The one public output example remains a checkable mathematical witness."""
import unittest

from tests.support import ROOT
from zipper_algebra import Cofactor, Presentation
from zipper_algebra.cofactor.verify import check, pit_check


class TestExampleOutput(unittest.TestCase):
    def test_curated_x3_witness_passes_both_exact_checks(self):
        presentation = Presentation(4)
        certificate = Cofactor.load(ROOT / "examples/x3/certificate.txt", pres=presentation)
        symbolic = check(certificate, presentation, "X_3", verbose=False)
        grid = pit_check(certificate, presentation, "X_3", verbose=False)
        self.assertEqual(len(certificate.terms), 129)
        self.assertTrue(symbolic["pass"])
        self.assertEqual(symbolic["residual_monomials"], 0)
        self.assertTrue(grid["pass"])
        self.assertEqual(grid["grid"], [25, 7])
