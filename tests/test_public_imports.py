"""Use the certificate interface with ordinary package imports and exact arithmetic."""
from fractions import Fraction
import unittest

from tools import zipper_in_ms
from tools.presentation import Presentation
from tools.cofactor.atoms import AtomSet
from tools.cofactor.cofactor_io import Cofactor
from tools.cofactor.verify import check, pit_check


class TestPackageImports(unittest.TestCase):
    def test_symbolic_and_numeric_objects_share_the_public_coefficient_classes(self):
        pres = Presentation(3)
        target = pres.relations(zipper_in_ms.SymbolicRing())[0]
        self.assertIsInstance(target, zipper_in_ms.Poly)
        cof = Cofactor.from_presentation(
            pres, "inverse relation", AtomSet([]),
            [((), pres.labels()[0], (), {(0, 0): Fraction(1)}, ())], "test",
        )
        for verifier in (check, pit_check):
            self.assertTrue(verifier(cof, pres, target, verbose=False)["pass"])
        self.assertEqual(cof.evaluate(3, 5, pres).terms,
                         pres.relations(zipper_in_ms.NumericRing(3, 5, 0))[0].terms)
