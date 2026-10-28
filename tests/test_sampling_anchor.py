"""Discovery anchors obey the same sampling restrictions as grown points."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from zipper_algebra import Atom, AtomSet, Presentation
from zipper_algebra.execution.config import RunConfig
from zipper_algebra.execution.sampling import SampleSet


class TestSamplingAnchor(unittest.TestCase):
    def test_excluded_or_singular_anchor_is_rejected_before_solver_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "runs"
            config = RunConfig(Path(sys.executable), output)
            atom = Atom("q-5", {(1, 0): 1, (0, 0): -5})
            for point, atoms in (((2, 7, 0), AtomSet([])),
                                 ((4, 7, 0), AtomSet([])),
                                 ((5, 7, 0), AtomSet([atom]))):
                with self.subTest(point=point), patch(
                    "zipper_algebra.execution.sampling.run_point"
                ) as execute:
                    samples = SampleSet(Presentation(2), "X_2", point, config, atoms)
                    with self.assertRaisesRegex(ValueError, "anchor"):
                        samples.run_anchor()
                    execute.assert_not_called()
            self.assertFalse(output.exists())
