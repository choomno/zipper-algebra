"""Discovery anchors obey the same sampling restrictions as grown points."""
from pathlib import Path
import sys
import tempfile
import pytest
from unittest.mock import patch

from zipper_algebra import Atom, AtomSet, Presentation
from zipper_algebra.execution.config import RunConfig
from zipper_algebra.execution.sampling import SampleSet


class TestSamplingAnchor:
    def test_excluded_or_singular_anchor_is_rejected_before_solver_execution(self, subtests):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "runs"
            config = RunConfig(Path(sys.executable), output)
            atom = Atom("q-5", {(1, 0): 1, (0, 0): -5})
            for point, atoms in (((2, 7, 0), AtomSet([])),
                                 ((4, 7, 0), AtomSet([])),
                                 ((5, 7, 0), AtomSet([atom]))):
                with subtests.test(point=point), patch(
                    "zipper_algebra.execution.sampling.run_point"
                ) as execute:
                    samples = SampleSet(Presentation(2), "X_2", point, config, atoms)
                    with pytest.raises(ValueError, match="anchor"):
                        samples.run_anchor()
                    execute.assert_not_called()
            assert not (output.exists())


@pytest.mark.parametrize("point", [(0, 5, 0), (-1, 5, 0), (3, 3, 0), (3, -9, 0), (5, 0, 0)])
def test_unseeded_anchor_does_not_assume_historical_denominator_factors(tmp_path, point):
    # Relations are polynomials over Q[q,t]. Only a supplied factor family
    # imposes its pole exclusions; generic certificate localizations are separate.
    config = RunConfig(Path(sys.executable), tmp_path / "runs")
    samples = SampleSet(Presentation(2), "X_2", point, config)
    with patch.object(samples, "run_point", return_value="accepted") as execute:
        assert samples.run_anchor() == "accepted"
    execute.assert_called_once_with(*point[:2], refresh=False)
