"""Bounded discovery outcomes with real validated process data and exact checks."""
from __future__ import annotations
import pytest

from contextlib import redirect_stdout
from fractions import Fraction
import io
import sys
from unittest.mock import patch

from tests.test_discovery_pipeline import DiscoveryFixture
from zipper_algebra import cli
from zipper_algebra.coefficients import NumericRing
from zipper_algebra.cofactor.lift import Lift, NeedMoreSamples
from zipper_algebra.execution import pipeline
from zipper_algebra.execution.sampling import SampleSet


class TestPipelineReconstructionOutcome(DiscoveryFixture):
    def fitted(self, coefficient=1):
        # The synthetic process emits g=X2=[i3]. This is an exact symbolic lift;
        # multiplying its input coefficient by two creates an actually false candidate.
        return Lift(self.atoms,
                    {'y1': [((), {(0,0): Fraction(coefficient)}, (), 'input', 3, ())]},
                    7, 4, (5,7,0), {'y1': 0})

    def controlled_run(self, fitted=None, *, multiplier_failure=False, unstable=False):
        effects = fitted if fitted is not None else [self.fitted(), self.fitted()]
        growth = []
        original_grow = SampleSet.grow
        def grow(samples, amount, **kwargs):
            growth.append(amount)
            return original_grow(samples, amount, **kwargs)
        actual_multipliers = pipeline.multiplier_samples
        attempts = 0
        def multipliers(*arguments):
            nonlocal attempts
            attempts += 1
            if multiplier_failure and attempts == 1:
                raise NeedMoreSamples('multiplier fit unresolved')
            return actual_multipliers(*arguments)
        with patch.object(pipeline, 'lift', side_effect=effects), \
             patch.object(pipeline, 'multiplier_samples', side_effect=multipliers), \
             patch.object(SampleSet, 'grow', grow):
            if unstable:
                # Compatibility is a sampling policy. Test its unavailable-modal
                # outcome without changing proof arithmetic or verification results.
                with patch.object(SampleSet, 'modal_shapes', return_value=({}, ['y1'])), \
                     patch.object(SampleSet, 'stable', return_value=[]):
                    result = self.discover(start=7, maxpts=9)
            else:
                result = self.discover(start=7, maxpts=9)
        return result, growth

    def test_failed_reconstruction_grows_to_cap_and_returns_inconclusive(self):
        result, growth = self.controlled_run([NeedMoreSamples('failed held-out fit')] * 2)
        assert (growth) == ([6, 2])
        assert (result['status']) == ('inconclusive')
        assert not (result['verified'])
        assert result['in_ideal'] is None

    def test_failed_exact_candidate_can_retry_with_more_data(self):
        result, growth = self.controlled_run([self.fitted(2), self.fitted()])
        assert (growth) == ([6, 2])
        assert (result['status']) == ('verified')
        assert result['verified']

    def test_failed_multiplier_fit_can_retry_with_more_data(self):
        result, growth = self.controlled_run(multiplier_failure=True)
        assert (growth) == ([6, 2])
        assert result['verified']

    def test_no_verified_candidate_at_cap_is_inconclusive(self):
        result, _ = self.controlled_run([self.fitted(2)] * 2)
        assert (result['status']) == ('inconclusive')
        assert not (result['verified'])
        assert ('failed exact verification') in (result['reason'])

    def test_cli_returns_nonzero_for_inconclusive_reconstruction(self):
        with patch.object(pipeline, 'run', return_value={'status':'inconclusive',
                                                        'verified':False, 'reason':'sample cap'}), \
             redirect_stdout(io.StringIO()):
            assert (cli.main(['discover', '--n', '2', '--target', 'X_2', '--solver', sys.executable, '--output', str(self.root)])) == (3)

    def test_unavailable_modal_shape_grows_to_cap_and_is_inconclusive(self):
        result, growth = self.controlled_run(unstable=True)
        assert (growth) == ([6, 2])
        assert (result['status']) == ('inconclusive')
        assert ('modal proof shapes') in (result['reason'])

    def test_defining_relation_nonzero_against_truncated_basis_is_inconclusive(self):
        # X2 is itself an imposed relation. The process returns only an inverse
        # relation, so its nonzero remainder cannot establish nonmembership.
        relation = self.pres.target(2,1,NumericRing(5,7,0))
        result = self.discover(config=self.configuration('omit'))
        assert (result['status']) == ('inconclusive')
        assert result['in_ideal'] is None
        assert (result['remainder_terms']) == (len(relation.terms))
        assert not (result['basis_complete'])

    @pytest.mark.parametrize('parameters', [
        pytest.param({'start': 0}, id='nonpositive-start'),
        pytest.param({'step': 0}, id='nonpositive-step'),
        pytest.param({'maxpts': 0}, id='nonpositive-cap'),
        pytest.param({'start': 6, 'maxpts': 5}, id='cap-below-start'),
    ])
    def test_invalid_sampling_limits_fail_before_execution(self, parameters):
        with pytest.raises(ValueError):
            self.discover(**parameters)
