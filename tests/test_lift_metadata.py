"""Exact numeric metadata at the explicit reader for historical unbound lifts."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from zipper_algebra import Presentation
from zipper_algebra.cofactor.atoms import AtomSet
from zipper_algebra.cofactor.compose import _composition_inputs, expand_over_inputs
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.cofactor.skeleton import Skeleton


class LiftMetadata(unittest.TestCase):
    def setUp(self):
        self.lift = Lift(AtomSet([]), {
            "g": [((), {(0, 0): Fraction(-2, 3)}, (1,), "input", 0, (2,))],
        }, 7, 2, (5, 7, 0), {"g": 0})
        self.data = self.lift.to_json()

    def test_exact_roundtrip_preserves_the_expanded_identity(self):
        loaded = Lift.from_json(json.loads(json.dumps(self.data)))
        self.assertEqual(loaded.to_json(), self.data)
        self.assertEqual(expand_over_inputs(loaded), expand_over_inputs(self.lift))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "lift.json"
            loaded.save(path)
            self.assertEqual(Lift.load(path).to_json(), self.data)

    def test_lossy_values_are_rejected_through_the_legacy_compose_reader(self):
        pres = Presentation(2)
        skeleton = Skeleton(pres.slug(), pres.content_hash(), "S[2][1]", (5, 7, 0))
        with tempfile.TemporaryDirectory() as directory:
            skeleton_path = Path(directory) / "skeleton.json"
            lift_path = Path(directory) / "lift.json"
            skeleton.save(skeleton_path)
            args = SimpleNamespace(skeleton=skeleton_path, lift=lift_path,
                                   target=None, solver_order="degree-lex")
            for field, value in (("anchor", 5.9), ("samples_used", 7.9),
                                 ("holdout", True), ("num", 0.1)):
                data = deepcopy(self.data)
                if field == "anchor":
                    data[field][0] = value
                elif field == "num":
                    data["entries"]["g"]["terms"][0]["num"][0][1] = value
                else:
                    data[field] = value
                lift_path.write_text(json.dumps(data))
                with self.subTest(field=field), self.assertRaisesRegex(ValueError, "lift"):
                    _composition_inputs(args, pres)

    def test_every_anchor_coordinate_rejects_floats_booleans_and_nonintegers(self):
        for index in range(3):
            for value in (5.0, 5.9, True, False, None, [], {}, "5.9", ""):
                data = deepcopy(self.data)
                data["anchor"][index] = value
                with self.subTest(index=index, value=value), self.assertRaisesRegex(ValueError, "integer"):
                    Lift.from_json(data)
        self.data["anchor"][2] = -1
        with self.assertRaisesRegex(ValueError, "characteristic"):
            Lift.from_json(self.data)

    def test_counts_and_indices_require_nonnegative_exact_integers(self):
        for field in ("samples_used", "holdout", "anchor_index"):
            for value in (0.0, 0.9, True, False, [], {}, "0.9", "", -1, "-1"):
                data = deepcopy(self.data)
                if field == "anchor_index":
                    data["entries"]["g"][field] = value
                else:
                    data[field] = value
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    Lift.from_json(data)
        for field in ("samples_used", "holdout"):
            self.data[field] = None
            with self.subTest(field=field), self.assertRaises(ValueError):
                Lift.from_json(self.data)
            self.data[field] = 0

    def test_integer_strings_and_missing_traceability_index_remain_supported(self):
        self.data.update(anchor=["-5", "+7", "0"], samples_used="7", holdout="2")
        self.data["entries"]["g"]["anchor_index"] = "0"
        loaded = Lift.from_json(self.data)
        self.assertEqual(loaded.anchor, (-5, 7, 0))
        self.assertEqual((loaded.samples_used, loaded.holdout), (7, 2))
        self.assertEqual(loaded.anchor_index, {"g": 0})
        self.data["entries"]["g"]["anchor_index"] = None
        self.data["anchor"][2] = "101"
        self.assertEqual(Lift.from_json(self.data).anchor_index, {"g": None})
        self.assertEqual(Lift.from_json(self.data).anchor, (-5, 7, 101))

    def test_numerator_coefficients_keep_exact_rational_strings(self):
        for value, expected in ((0, 0), (-3, -3), ("-2/3", Fraction(-2, 3)),
                                ("0.1", Fraction(1, 10)), ("1e-2", Fraction(1, 100)),
                                (Fraction(7, 11), Fraction(7, 11))):
            self.data["entries"]["g"]["terms"][0]["num"][0][1] = value
            with self.subTest(value=value):
                loaded = Lift.from_json(self.data)
                self.assertEqual(loaded.entries["g"][0][1], {(0, 0): expected})

    def test_numerator_coefficients_reject_floats_booleans_and_invalid_rationals(self):
        for value in (0.1, 1.0, True, False, None, [], {}, "1/0", "invalid", "NaN"):
            self.data["entries"]["g"]["terms"][0]["num"][0][1] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "exact rational"):
                Lift.from_json(self.data)

    def test_object_and_anchor_shapes_are_explicit(self):
        for value in (None, [], "lift", 3):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "JSON object"):
                Lift.from_json(value)
        for value in (None, {}, "570", 5, [], [5, 7], [5, 7, 0, 1]):
            self.data["anchor"] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "anchor.*three"):
                Lift.from_json(self.data)
        self.data["anchor"] = (5, 7, 0)
        self.assertEqual(Lift.from_json(self.data), self.lift)

    def test_empty_lift_with_zero_counts_roundtrips(self):
        empty = Lift(AtomSet([]), {}, 0, 0, (0, -1, 0), {})
        self.assertEqual(Lift.from_json(empty.to_json()).to_json(), empty.to_json())


if __name__ == "__main__":
    unittest.main()
