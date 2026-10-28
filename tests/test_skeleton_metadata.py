"""Exact metadata at the explicit reader for historical unbound skeletons."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from zipper_algebra import Presentation
from zipper_algebra.coefficients import Coef, NumericRing
from zipper_algebra.polynomial import Poly
from zipper_algebra.cofactor.compose import _composition_inputs
from zipper_algebra.cofactor.proof import Basis, Proof
from zipper_algebra.cofactor.skeleton import Skeleton, compose_at, evaluate, reduce


class SkeletonMetadata(unittest.TestCase):
    def setUp(self):
        self.pres = Presentation(2)
        self.inputs = self.pres.relations(NumericRing(5, 7, 0))
        self.basis = Basis([self.inputs[0]], 0, self.pres.var_names())
        self.proof = Proof([[(Fraction(1), (), "input", 0, ())]])
        left = Poly(0)
        left.add_term((1,), Coef.one(0))
        right = Poly(0)
        right.add_term((2,), Coef.one(0))
        self.target = (left * self.inputs[0] * right).scalar_mul(Coef(Fraction(-2, 3), 0))
        steps, remainder = reduce(self.target, self.basis)
        self.skeleton = Skeleton.from_reduction(
            self.pres, "test polynomial", (5, 7, 0), steps, remainder, self.basis, self.proof)
        self.data = self.skeleton.to_json()

    def test_fractional_reduction_roundtrip_reconstructs_its_target(self):
        self.assertEqual(self.skeleton.steps[0]["mult"], Fraction(-2, 3))
        loaded = Skeleton.from_json(json.loads(json.dumps(self.data)))
        self.assertEqual(loaded, self.skeleton)
        self.assertEqual(loaded.steps[0]["left"], (1,))
        self.assertEqual(loaded.steps[0]["right"], (2,))
        result = evaluate(compose_at(loaded, self.basis, self.proof), self.inputs, 0)
        self.assertEqual(result.terms, self.target.terms)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "skeleton.json"
            loaded.save(path)
            self.assertEqual(Skeleton.load(path), loaded)

    def test_three_reported_lossy_values_are_rejected_by_legacy_compose_reader(self):
        with tempfile.TemporaryDirectory() as directory:
            skeleton_path = Path(directory) / "skeleton.json"
            lift_path = Path(directory) / "lift.json"
            lift_path.write_text("{}")
            args = SimpleNamespace(skeleton=skeleton_path, lift=lift_path,
                                   target=None, solver_order="degree-lex")
            for field, value in (("anchor", 5.9), ("remainder_terms", 0.9), ("mult", 0.1)):
                data = deepcopy(self.data)
                if field == "anchor":
                    data["anchor"]["q"] = value
                elif field == "mult":
                    data["steps"][0]["mult"] = value
                else:
                    data[field] = value
                skeleton_path.write_text(json.dumps(data))
                with self.subTest(field=field), self.assertRaisesRegex(ValueError, "skeleton"):
                    _composition_inputs(args, self.pres)

    def test_every_anchor_coordinate_rejects_float_boolean_and_other_nonintegers(self):
        for field in ("q", "t", "char"):
            for value in (5.0, 5.9, True, False, None, [], {}, "5.9", ""):
                data = deepcopy(self.data)
                data["anchor"][field] = value
                with self.subTest(field=field, value=value), self.assertRaisesRegex(ValueError, "integer"):
                    Skeleton.from_json(data)

    def test_counts_reject_lossy_values_negative_values_and_wrong_step_totals(self):
        for field in ("remainder_terms", "n_steps"):
            for value in (0.0, 0.9, True, False, None, [], {}, "0.9", "", -1, "-1"):
                data = deepcopy(self.data)
                data[field] = value
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    Skeleton.from_json(data)
        self.data["n_steps"] += 1
        with self.assertRaisesRegex(ValueError, "step count"):
            Skeleton.from_json(self.data)

    def test_legacy_integer_strings_and_exact_rational_inputs(self):
        self.data["anchor"] = {"q": "5", "t": "+7", "char": "0"}
        self.data["remainder_terms"] = "0"
        self.data["n_steps"] = str(len(self.skeleton.steps))
        self.assertEqual(Skeleton.from_json(self.data), self.skeleton)
        del self.data["n_steps"]  # Older readers did not require this redundant count.
        self.assertEqual(Skeleton.from_json(self.data), self.skeleton)
        for value, expected in ((0, 0), (-3, -3), ("-2/3", Fraction(-2, 3)),
                                ("0.1", Fraction(1, 10)), ("1e-2", Fraction(1, 100)),
                                (Fraction(7, 11), Fraction(7, 11))):
            self.data["steps"][0]["mult"] = value
            with self.subTest(value=value):
                self.assertEqual(Skeleton.from_json(self.data).steps[0]["mult"], expected)

    def test_multiplier_rejects_float_bool_nonrational_and_zero_denominator(self):
        for value in (0.1, 1.0, True, False, None, [], {}, "1/0", "invalid", "NaN"):
            self.data["steps"][0]["mult"] = value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "exact rational"):
                Skeleton.from_json(self.data)

    def test_contexts_require_arrays_of_positive_integer_indices(self):
        for field in ("left", "right"):
            for value in ("12", 1, None, {}, [True], [1.0], [0], [-1], ["1"], [None]):
                data = deepcopy(self.data)
                data["steps"][0][field] = value
                with self.subTest(field=field, value=value), self.assertRaisesRegex(ValueError, "context"):
                    Skeleton.from_json(data)
        self.data["steps"][0].update(left=[], right=())
        loaded = Skeleton.from_json(self.data)
        self.assertEqual(loaded.steps[0]["left"], ())
        self.assertEqual(loaded.steps[0]["right"], ())

    def test_object_and_array_fields_reject_scalar_or_ambiguous_containers(self):
        for value in (None, [], "skeleton"):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "JSON object"):
                Skeleton.from_json(value)
        for field in ("anchor", "steps", "basis_keys"):
            for value in (None, 1, "data"):
                data = deepcopy(self.data)
                data[field] = value
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    Skeleton.from_json(data)
        for value in (None, [], "step", 1):
            self.data["steps"] = [value]
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "step.*JSON object"):
                Skeleton.from_json(self.data)

    def test_text_fields_and_basis_keys_require_strings(self):
        for field in ("presentation", "presentation_hash", "target", "basis_keys", "basis_key"):
            for value in (None, True, 1, [], {}):
                data = deepcopy(self.data)
                if field == "basis_key":
                    data["steps"][0][field] = value
                elif field == "basis_keys":
                    data[field] = [value]
                else:
                    data[field] = value
                with self.subTest(field=field, value=value), self.assertRaisesRegex(ValueError, "string"):
                    Skeleton.from_json(data)
        for field in ("basis_keys", "basis_key"):
            data = deepcopy(self.data)
            if field == "basis_keys":
                data[field] = [""]
            else:
                data["steps"][0][field] = ""
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "nonempty"):
                Skeleton.from_json(data)

    def test_empty_zero_skeleton_and_scalar_remainder_remain_valid(self):
        empty = Skeleton("", "", "", (-5, 0, 0))
        self.assertEqual(Skeleton.from_json(empty.to_json()), empty)
        scalar = Poly(0)
        scalar.add_term((), Coef(7, 0))
        steps, remainder = reduce(scalar, self.basis)
        skeleton = Skeleton.from_reduction(self.pres, "7", (5, 7, 0), steps, remainder,
                                           self.basis, self.proof)
        self.assertEqual(skeleton.steps, [])
        self.assertEqual(skeleton.remainder_terms, 1)
        self.assertFalse(skeleton.in_ideal)
        self.assertEqual(Skeleton.from_json(skeleton.to_json()), skeleton)

    def test_negative_characteristic_rejected_and_positive_characteristic_preserved(self):
        self.data["anchor"]["char"] = -1
        with self.assertRaisesRegex(ValueError, "characteristic"):
            Skeleton.from_json(self.data)
        self.data["anchor"]["char"] = "101"
        self.assertEqual(Skeleton.from_json(self.data).anchor, (5, 7, 101))


if __name__ == "__main__":
    unittest.main()
