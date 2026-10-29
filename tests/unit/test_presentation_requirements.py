"""Algebra choices must be explicit and must not silently select another presentation."""
import pytest

from zipper_algebra import Presentation
from zipper_algebra.coefficients import SymbolicRing


class TestPresentationRequirements:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    @pytest.mark.parametrize("value", ["typo", "Bigelow", "indeterminate", None, []])
    def test_unknown_untwisting_cannot_fall_back_to_canonical(self, value):
        with pytest.raises(ValueError):
            Presentation(4, untwist=value)

    def test_integer_domains_and_level_limits(self):
        for n in (True, 1, 4.0, "4"):
            with self.subtests.test(n=n):
                with pytest.raises(ValueError):
                    Presentation(n)
        for field in ("framing_power", "theta_power"):
            for value in (True, 0, -1, 1.0):
                with self.subtests.test(field=field, value=value):
                    with pytest.raises(ValueError):
                        Presentation(4, **{field: value})
        for field in ("max_zipper_level", "max_untwist_level", "zero_from"):
            for value in (True, 1, 5, 2.0):
                with self.subtests.test(field=field, value=value):
                    with pytest.raises(ValueError):
                        Presentation(4, **{field: value})

    def test_extra_relations_exist_and_have_unique_labels(self):
        for extra in (("unknown",), ("X_4", "X_4"), ["X_4"], (1,), ("X_5",)):
            with self.subtests.test(extra=extra):
                with pytest.raises(ValueError):
                    Presentation(4, extra=extra)
        pres = Presentation(4, extra=("X_4", "(Theta_3-q^2)X_3"))
        assert len(pres.labels()) == len(set(pres.labels()))

    def test_indices_are_not_coerced(self):
        pres = Presentation(4)
        for indices in ((2, True), (2.0, 1), (1, 1), (4, 2), (5, 1)):
            with self.subtests.test(indices=indices):
                with pytest.raises(ValueError):
                    pres.target(*indices, SymbolicRing())

    def test_default_relations_preserve_historical_hash(self):
        # Computed independently from the preserved pre-packaging checkout.
        assert Presentation(4).content_hash() == '091bbb92271ab870'
        assert Presentation(4).n_relations() == 31
        assert 'T_2' in Presentation(2, max_untwist_level=2).labels()
