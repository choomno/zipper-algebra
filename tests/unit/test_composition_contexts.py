"""Synthetic regressions for contextual composition and step multiplier sampling.

Run with ``python -m pytest tests/unit/test_composition_contexts.py``.
No solver, precomputed basis, sibling research directory or cache is used.
"""
from __future__ import annotations

import sys
import pytest
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tests.support import ROOT
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.compose import Rat, compose, lift_multipliers
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.presentation import Presentation
from zipper_algebra.cofactor.proof import Basis
from zipper_algebra.cofactor.skeleton import Skeleton
from zipper_algebra.coefficients import Coef, NumericRing
from zipper_algebra.polynomial import Poly
from zipper_algebra.execution import pipeline


def polynomial(terms):
    out = Poly(0)
    for word, coefficient in terms:
        out.add_term(tuple(word), Coef(coefficient, 0))
    return out


class TestContextualComposition:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    @pytest.fixture(autouse=True)
    def setup_algebra(self):
        self.pres = Presentation(3)
        self.atoms = AtomSet([Atom("q", {(1, 0): 1}), Atom("t", {(0, 1): 1})])
        self.lift = Lift(self.atoms, {
            "inner": [((1, 0), {(0, 0): Fraction(1)}, (2,), "input", 0, (3,))],
            "outer": [((0, 1), {(0, 0): Fraction(1)}, (4,), "basis", "inner", (1,))],
        }, samples_used=0, holdout=0, anchor=(5, 7, 0), anchor_index={})

    def skeleton(self, steps):
        return Skeleton(self.pres.slug(), self.pres.content_hash(), "test target", (5, 7, 0),
                        [{"basis_key": key, "left": left, "right": right,
                          "mult": Fraction(1)} for key, left, right in steps])

    def constant(self, value, denominator=(0, 0)):
        return ({(0, 0): Fraction(value)} if value else {}, denominator)

    def test_noncommuting_contexts_and_multiple_denominator_atoms(self):
        skel = self.skeleton([("inner", (1,), (4,)), ("inner", (4,), (1,))])
        cof = compose(self.lift, skel, self.pres,
                      multipliers=[self.constant(1, (0, 1)), self.constant(2, (0, 1))],
                      verbose=False)
        label = self.pres.labels()[0]
        by_context = {(left, lab, right): (num, den)
                      for left, lab, right, num, den in cof.terms}
        assert by_context == {((1, 2), label, (3, 4)): self.constant(1, (1, 1)), ((4, 2), label, (3, 1)): self.constant(2, (1, 1))}
        for q, t in ((5, 7), (11, 3)):
            relation = self.pres.relations(NumericRing(q, t, 0))[0]
            expected = Poly(0)
            for left, right, coefficient in (((1, 2), (3, 4), 1), ((4, 2), (3, 1), 2)):
                expected = expected + (
                    polynomial([(left, 1)]) * relation * polynomial([(right, 1)])
                ).scalar_mul(Coef(Fraction(coefficient, q * t), 0))
            assert cof.evaluate(q, t, self.pres).terms == expected.terms

    def test_nested_basis_proof_preserves_context_order(self):
        skel = self.skeleton([("outer", (1, 2), (3, 4))])
        cof = compose(self.lift, skel, self.pres, multipliers=[self.constant(1)], verbose=False)
        assert cof.terms == [((1, 2, 4, 2), self.pres.labels()[0], (3, 1, 3, 4), {(0, 0): Fraction(1)}, (1, 1))]

    def test_repeated_identical_contexts_add_distinct_multipliers(self):
        skel = self.skeleton([("inner", (1,), (4,)), ("inner", (1,), (4,))])
        cof = compose(self.lift, skel, self.pres,
                      multipliers=[self.constant(1, (0, 1)), self.constant(2, (0, 1))],
                      verbose=False)
        assert cof.terms == [((1, 2), self.pres.labels()[0], (3, 4), {(0, 0): Fraction(3)}, (1, 1))]

    def test_repeated_contexts_can_cancel(self):
        skel = self.skeleton([("inner", (1,), (4,)), ("inner", (1,), (4,))])
        cof = compose(self.lift, skel, self.pres,
                      multipliers=[self.constant(1, (0, 1)), self.constant(-1, (0, 1))],
                      verbose=False)
        assert cof.terms == []

    @pytest.mark.parametrize("count", [0, 1, 3])
    def test_wrong_multiplier_count_is_rejected(self, count):
        skel = self.skeleton([("inner", (1,), (4,)), ("inner", (4,), (1,))])
        with pytest.raises(ValueError, match="multipliers"):
            compose(self.lift, skel, self.pres,
                    multipliers=[self.constant(1)] * count, verbose=False)

    def test_nonzero_skeleton_remainder_is_rejected(self):
        skel = self.skeleton([("inner", (), ())])
        skel.remainder_terms = 1
        with pytest.raises(ValueError, match='nonzero remainder'):
            compose(self.lift, skel, self.pres, multipliers=[self.constant(1)], verbose=False)

    def test_mismatched_presentation_metadata_is_rejected(self):
        for owner in ("skeleton", "lift"):
            for field, bad_value in (("presentation_hash", "not-the-presentation-hash"),
                                     ("presentation", "different-presentation")):
                with self.subtests.test(owner=owner, field=field):
                    skel = self.skeleton([("inner", (), ())])
                    data = skel if owner == "skeleton" else self.lift
                    previous = getattr(data, field, None)
                    setattr(data, field, bad_value)
                    try:
                        with pytest.raises(ValueError):
                            compose(self.lift, skel, self.pres,
                                    multipliers=[self.constant(1)], verbose=False)
                    finally:
                        if previous is None:
                            delattr(data, field)
                        else:
                            setattr(data, field, previous)

    def test_lift_repeated_step_multipliers_independently(self):
        skel = self.skeleton([("inner", (1,), (4,)), ("inner", (4,), (1,))])
        samples = {(q, t): [Fraction(1, q), Fraction(2, t)]
                   for q, t in ((3, 5), (5, 7), (7, 11), (11, 13))}
        fitted = lift_multipliers(skel, samples, self.atoms, max_pq=0, max_pt=0)
        rat = Rat(self.atoms)
        for q, t in ((13, 17), (19, 23)):
            assert [rat.eval(value, q, t) for value in fitted] == [Fraction(1, q), Fraction(2, t)]

    def test_legacy_mapping_is_supported_only_for_unique_keys(self):
        skel = self.skeleton([("inner", (), ()), ("outer", (), ())])
        samples = {(q, t): {"inner": Fraction(1), "outer": Fraction(2)}
                   for q, t in ((3, 5), (5, 7), (7, 11))}
        fitted = lift_multipliers(skel, samples, self.atoms, max_pq=0, max_pt=0)
        assert [Rat(self.atoms).eval(value, 13, 17) for value in fitted] == [1, 2]
        skel.steps[1]["basis_key"] = "inner"
        with pytest.raises(ValueError, match='ambiguous'):
            lift_multipliers(skel, samples, self.atoms)

    def test_malformed_multiplier_samples_are_rejected(self):
        skel = self.skeleton([("inner", (), ())])
        for values in ([], [1, 2], {"outer": 1}):
            with self.subtests.test(values=values):
                with pytest.raises(ValueError):
                    lift_multipliers(skel, {(3, 5): values}, self.atoms)

    def test_invalid_input_indices_cannot_select_python_negative_indices(self):
        skel = self.skeleton([("inner", (), ())])
        for index in (-1, len(self.pres.labels()), True, 0.5):
            with self.subtests.test(index=index):
                self.lift.entries["inner"] = [((0, 0), {}, (), "input", index, ())]
                with pytest.raises(ValueError, match='input relation index'):
                    compose(self.lift, skel, self.pres,
                            multipliers=[self.constant(0)], verbose=False)

    def test_unknown_kinds_and_missing_dependencies_are_rejected_before_cancellation(self):
        skel = self.skeleton([("outer", (), ())])
        self.lift.entries["outer"] = [((0, 0), {}, (), "typo", "inner", ())]
        with pytest.raises(ValueError, match='unknown proof kind'):
            compose(self.lift, skel, self.pres, multipliers=[self.constant(0)], verbose=False)
        self.lift.entries["outer"] = [((0, 0), {}, (), "basis", "missing", ())]
        with pytest.raises(KeyError, match='not in the lift'):
            compose(self.lift, skel, self.pres, multipliers=[self.constant(0)], verbose=False)

    def test_recursive_dependency_cycle_is_rejected(self):
        skel = self.skeleton([("inner", (), ())])
        self.lift.entries["inner"] = [((0, 0), {}, (), "basis", "outer", ())]
        with pytest.raises(ValueError, match='cycle'):
            compose(self.lift, skel, self.pres, multipliers=[self.constant(0)], verbose=False)

    def test_every_lift_entry_is_validated_before_zero_multipliers_hide_it(self):
        skel = self.skeleton([("inner", (), ())])
        malformed = [((0,), {}, (), "input", 0, ()),
                     ((-1, 0), {}, (), "input", 0, ()),
                     ((0, 0), {(-1, 0): 1}, (), "input", 0, ()),
                     ((0, 0), {(0, 0): 0.5}, (), "input", 0, ()),
                     ((0, 0), {}, (0,), "input", 0, ()),
                     ((0, 0), {}, (), "input", 0, (5,))]
        for entry in malformed:
            with self.subtests.test(entry=entry):
                # 'outer' is not requested by this skeleton, but malformed stored entries
                # are rejected consistently rather than disappearing with coefficient zero.
                self.lift.entries["outer"] = [entry]
                with pytest.raises(ValueError):
                    compose(self.lift, skel, self.pres,
                            multipliers=[self.constant(0)], verbose=False)

    def test_invalid_multipliers_are_rejected_even_when_the_result_is_zero(self):
        skel = self.skeleton([("inner", (), ())])
        malformed = [({}, (0,)), ({}, (-1, 0)), ({}, (0.5, 0)),
                     ({(-1, 0): 1}, (0, 0)), ({(0, 0): 0.5}, (0, 0)), ({},)]
        for multiplier in malformed:
            with self.subtests.test(multiplier=multiplier):
                with pytest.raises(ValueError):
                    compose(self.lift, skel, self.pres, multipliers=[multiplier], verbose=False)

    def test_invalid_skeleton_contexts_are_rejected_with_zero_multipliers(self):
        for left, right in (((0,), ()), ((), (5,)), ((True,), ())):
            with self.subtests.test(left=left, right=right):
                skel = self.skeleton([("inner", left, right)])
                with pytest.raises(ValueError):
                    compose(self.lift, skel, self.pres,
                            multipliers=[self.constant(0)], verbose=False)


class TestMultiplierSampling:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    @pytest.fixture(autouse=True)
    def setup_algebra(self):
        self.pres = Presentation(3)
        self.steps = [((4,), 0, (2,), Coef(7, 0)), ((2,), 0, (4,), Coef(5, 0))]
        self.skel = Skeleton(self.pres.slug(), self.pres.content_hash(), "test target", (5, 7, 0),
                             [{"left": left, "basis_key": "x1", "right": right,
                               "mult": Fraction(mult.v)} for left, _, right, mult in self.steps])

    def test_actual_reductions_match_by_context_across_basis_renumbering(self):
        points = ((5, 7), (11, 13))
        g = polynomial([((1,), 1)])
        h = polynomial([((3,), 1)])
        bases = {(5, 7): Basis([g, h], 0, self.pres.var_names()),
                 (11, 13): Basis([h, g], 0, self.pres.var_names())}
        samples = SimpleNamespace(char=0, stable=lambda sk: points,
                                  basis=lambda q, t: bases[(q, t)])
        # Only target construction is substituted: sampling uses the real reducer and Basis.
        def target(_pres, _name, ring):
            return polynomial([((4, 1, 2), ring.t.v), ((2, 1, 4), ring.q.v)])
        with patch.object(pipeline, "resolve_target", side_effect=target):
            result = pipeline.multiplier_samples(self.skel, samples, self.pres)
        assert result == {(5, 7): [Fraction(7), Fraction(5)], (11, 13): [Fraction(13), Fraction(11)]}

    def test_repeated_identical_context_occurrences_are_retained(self):
        self.skel.steps[1] = dict(self.skel.steps[0])
        repeated = [(self.steps[0][0], 0, self.steps[0][2], Coef(mult, 0)) for mult in (2, 3)]
        assert self.sample_mocked(repeated) == {(5, 7): [Fraction(2), Fraction(3)]}

    def sample_mocked(self, steps, remainder=None):
        samples = SimpleNamespace(char=0, stable=lambda sk: [(5, 7)],
                                  basis=lambda q, t: SimpleNamespace(key=["x1", "other"]))
        if remainder is None:
            remainder = Poly(0)
        with patch.object(pipeline, "resolve_target", return_value=Poly(0)), \
             patch.object(pipeline, "reduce", return_value=(steps, remainder)):
            return pipeline.multiplier_samples(self.skel, samples, self.pres)

    def test_incompatible_reduction_schedules_are_skipped(self):
        changed_context = [((3,), 0, (2,), Coef(7, 0)), self.steps[1]]
        changed_key = [(self.steps[0][0], 1, self.steps[0][2], Coef(7, 0)), self.steps[1]]
        for steps in (list(reversed(self.steps)), self.steps[:1], self.steps + self.steps[:1],
                      changed_context, changed_key):
            with self.subtests.test(steps=steps):
                assert self.sample_mocked(steps) == {}

    def test_nonzero_remainder_is_skipped(self):
        assert self.sample_mocked(self.steps, polynomial([((), 1)])) == {}
