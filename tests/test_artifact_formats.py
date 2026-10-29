"""Exact data, context binding and atomic publication of mathematical artifacts."""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from zipper_algebra import Presentation
from zipper_algebra.cofactor import artifacts as A
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.cofactor_io import Cofactor
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.cofactor.skeleton import Skeleton
from zipper_algebra.cofactor.verify import check


class MathematicalArtifacts(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.path = self.root / "artifact.json"
        self.pres = Presentation(2, zero_from=2)
        self.atoms = AtomSet([Atom("q", {(1, 0): 1}), Atom("t", {(0, 1): 1})], [2, 3])
        label = "ZERO_2_1"
        index = self.pres.labels().index(label)
        self.cofactor = Cofactor.from_presentation(
            self.pres, "X_2", self.atoms, [((), label, (), {(0, 0): F(1)}, (0, 0))],
            producer="test", anchor=(5, 7, 0), notes=["exact contextual data"])
        self.skeleton = Skeleton(
            self.pres.slug(), self.pres.content_hash(), "X_2", (5, 7, 0),
            [{"basis_key": "g", "left": (1,), "right": (), "mult": F(1, 2)},
             {"basis_key": "g", "left": (), "right": (2,), "mult": F(-3)}],
            ["g", "h"], 0)
        self.lift = Lift(self.atoms,
                         {"g": [((0, 0), {(0, 0): F(1)}, (), "input", index, ())],
                          "h": [((1, 0), {(1, 0): F(1, 2)}, (1,), "basis", "g", ()),
                                ((0, 0), {(0, 1): F(-3)}, (), "basis", "g", (2,))]},
                         12, 4, (5, 7, 0), {"g": 0, "h": None}, timing={"fit": 0.5})
        self.order = {"solver": "degree-lex", "reduction": "length-lex-index-v1"}

    def save(self, value=None, **kwargs):
        A.save_artifact(self.path, self.cofactor if value is None else value,
                        self.pres, "X_2", **kwargs)

    def load(self, **kwargs):
        return A.load_artifact(self.path, self.pres, "X_2", **kwargs)

    def mutate(self, change, *, redigest=False):
        data = json.loads(self.path.read_text())
        change(data)
        if redigest:
            payload = json.dumps(data["payload"], sort_keys=True, separators=(",", ":"),
                                 ensure_ascii=True, allow_nan=False).encode()
            data["payload_sha256"] = hashlib.sha256(payload).hexdigest()
        self.path.write_text(json.dumps(data))

    def test_presentation_identity_records_actual_ordered_definition(self):
        identity = A.presentation_identity(self.pres)
        self.assertEqual(identity["coefficient_domain"], "Q[q,t]")
        self.assertEqual(identity["generators"], ["x1", "y1"])
        self.assertEqual(identity["relation_labels"], self.pres.labels())
        self.assertEqual(identity["canonical_relations"], self.pres.canonical_text())
        self.assertEqual(identity["sha256"], hashlib.sha256(self.pres.canonical_text().encode()).hexdigest())
        self.assertEqual(identity["legacy_hash"], self.pres.content_hash())
        self.assertNotIn("order", identity)

    def test_target_identity_agrees_for_aliases_and_direct_polynomials(self):
        expected = {"coefficient_domain": "Q[q,t]", "terms": [
            [[], [[[0, 0], "1"], [[1, 0], "-1"]]],
            [[1], [[[0, 0], "-1"]]], [[2], [[[1, 0], "1"]]]],
                    "sha256": ""}
        actual = A.target_identity(self.pres, "X_2")
        self.assertEqual(actual["terms"], expected["terms"])
        self.assertEqual(actual, A.target_identity(self.pres, "S[2][1]"))
        self.assertEqual(actual, A.target_identity(self.pres, {(1,): {(0, 0): -1},
                        (2,): {(1, 0): 1}, (): {(0, 0): 1, (1, 0): -1}}))
        self.assertEqual(A.target_identity(self.pres, {})["terms"], [])

    def test_roundtrip_cofactor_preserves_exact_terms_and_legacy_reader(self):
        self.save()
        loaded = self.load(kind="cofactor", atoms=self.atoms)
        self.assertEqual(loaded.terms, self.cofactor.terms)
        self.assertEqual(loaded.header, self.cofactor.header)
        self.assertTrue(check(loaded, self.pres, verbose=False)["pass"])
        self.assertIsNone(json.loads(self.path.read_text())["order"])
        legacy = self.root / "legacy.txt"
        self.cofactor.save(legacy)
        self.assertEqual(Cofactor.load(legacy, pres=self.pres).terms, loaded.terms)
        with self.assertRaises(json.JSONDecodeError):
            A.load_artifact(legacy, self.pres, "X_2")

    def test_roundtrip_cofactor_optional_header_fields_and_no_verification_claim(self):
        false = Cofactor([], 2, self.atoms, target="X_2")
        self.save(false)
        data = json.loads(self.path.read_text())
        self.assertNotIn("verified", data)
        self.assertIsNone(data["payload"]["anchor"])
        self.assertEqual(data["payload"]["notes"], [])
        self.assertFalse(check(self.load(), self.pres, verbose=False)["pass"])

    def test_skeleton_roundtrip_retains_repeated_contextual_occurrences(self):
        self.save(self.skeleton, order=self.order)
        loaded = self.load(kind="skeleton", order=deepcopy(self.order))
        self.assertEqual(loaded, self.skeleton)
        self.assertEqual(loaded.step_keys(), ["g", "g"])
        self.assertEqual(loaded.steps[0]["left"], (1,))
        self.assertEqual(loaded.steps[1]["right"], (2,))
        with self.assertRaisesRegex(ValueError, "atoms"):
            self.load(order=self.order, atoms=self.atoms)

    def test_lift_roundtrip_retains_rationals_references_and_atom_order(self):
        self.save(self.lift, order="length-lex-index-v1")
        loaded = self.load(kind="lift", order="length-lex-index-v1", atoms=self.atoms)
        self.assertEqual(loaded, self.lift)
        self.assertEqual(loaded.to_json(), self.lift.to_json())
        self.assertEqual([term[4] for term in loaded.entries["h"]], ["g", "g"])
        reversed_atoms = AtomSet(list(reversed(self.atoms.atoms)), [3, 2])
        with self.assertRaisesRegex(ValueError, "atoms"):
            self.load(order="length-lex-index-v1", atoms=reversed_atoms)
        mutated_atoms = deepcopy(self.atoms)
        mutated_atoms.max_exp[0] = True
        with self.assertRaisesRegex(ValueError, "exponents"):
            self.load(order="length-lex-index-v1", atoms=mutated_atoms)

    def test_unknown_objects_cannot_be_saved(self):
        with self.assertRaises(TypeError):
            self.save({})

    def test_order_is_required_for_order_dependent_objects(self):
        for value in (self.skeleton, self.lift):
            for order in (None, "", {}, 12):
                with self.subTest(kind=type(value), order=order), self.assertRaises(ValueError):
                    self.save(value, order=order)
        with self.assertRaises(ValueError):
            self.save(order=self.order)

    def test_order_and_kind_mismatch_are_rejected(self):
        self.save(self.lift, order=self.order)
        with self.assertRaisesRegex(ValueError, "kind"):
            self.load(kind="skeleton", order=self.order)
        with self.assertRaisesRegex(ValueError, "order"):
            self.load(order={"solver": "other", "reduction": "length-lex-index-v1"})

    def test_cross_presentation_and_target_are_rejected(self):
        self.save()
        with self.assertRaisesRegex(ValueError, "presentation or target"):
            A.load_artifact(self.path, Presentation(2), "X_2")
        with self.assertRaisesRegex(ValueError, "presentation or target"):
            A.load_artifact(self.path, self.pres, {})
        self.assertEqual(A.load_artifact(self.path, self.pres, "S[2][1]").terms, self.cofactor.terms)

    def test_schema_rejects_unknown_fields_unbound_and_nonobject_data(self):
        for data in ([], {}, {"schema": "other"}, self.lift.to_json()):
            self.path.write_text(json.dumps(data))
            with self.subTest(data=data), self.assertRaisesRegex(ValueError, "schema"):
                self.load()
        self.save()
        self.mutate(lambda data: data.update(verified=True))
        with self.assertRaisesRegex(ValueError, "schema"):
            self.load()

    def test_domain_alphabet_relation_and_target_metadata_are_checked(self):
        for field, replacement in (("coefficient_domain", "GF(5)"), ("generators", ["y1", "x1"]),
                                   ("relation_labels", []), ("canonical_relations", "different"),
                                   ("sha256", "0" * 64), ("legacy_hash", "0" * 16)):
            self.save()
            self.mutate(lambda data: data["presentation"].update({field: replacement}))
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.load()
        self.save()
        self.mutate(lambda data: data["target"].update(coefficient_domain="Q"))
        with self.assertRaises(ValueError):
            self.load()

    def test_payload_digest_detects_coefficient_and_atom_mutation(self):
        for mutate in (lambda data: data["payload"]["terms"][0].update(num=[[[0, 0], "2"]]),
                       lambda data: data["payload"]["atoms"]["atoms"].reverse()):
            self.save()
            self.mutate(mutate)
            with self.assertRaisesRegex(ValueError, "digest"):
                self.load()

    def test_boolean_and_float_metadata_do_not_equal_integer_indices(self):
        for replacement in (True, 1.0):
            self.save()
            self.mutate(lambda data: data["target"]["terms"][1][0].__setitem__(0, replacement))
            with self.assertRaises(ValueError):
                self.load()
        order = {"alphabet_origin": 1}
        self.save(self.lift, order=order)
        self.mutate(lambda data: data["order"].update(alphabet_origin=True))
        with self.assertRaisesRegex(ValueError, "order"):
            self.load(order=order)

    def test_recomputed_digest_does_not_hide_invalid_or_noncanonical_data(self):
        changes = [lambda data: data["payload"].update(unexpected="ignored"),
                   lambda data: data["payload"]["terms"][0].update(num=[[[0, 0], 0.5]]),
                   lambda data: data["payload"]["terms"][0].update(num=[[[0, 0], "1"], [[0, 0], "1"]])]
        for mutate in changes:
            self.save()
            self.mutate(mutate, redigest=True)
            with self.assertRaisesRegex(ValueError, "canonical"):
                self.load()
        self.save()
        self.mutate(lambda data: data["payload"]["terms"][0].update(denom=[-1, 0]), redigest=True)
        with self.assertRaises(ValueError):
            self.load()

    def test_unknown_kind_cannot_select_an_unchecked_decoder(self):
        self.save(self.lift, order=self.order)
        self.mutate(lambda data: data.update(kind="unknown"))
        with self.assertRaisesRegex(ValueError, "unknown.*kind"):
            self.load(order=self.order)

    def test_new_skeleton_rejects_inconsistent_counts_or_contexts(self):
        cases = [lambda value: setattr(value, "anchor", (5, 7)),
                 lambda value: setattr(value, "anchor", (True, 7, 0)),
                 lambda value: setattr(value, "anchor", (5, 7, 5)),
                 lambda value: setattr(value, "remainder_terms", -1),
                 lambda value: setattr(value, "remainder_terms", True),
                 lambda value: value.basis_keys.append("g"),
                 lambda value: value.basis_keys.append([]),
                 lambda value: value.basis_keys.append(""),
                 lambda value: value.steps[0].update(basis_key="missing"),
                 lambda value: value.steps[0].update(left=(3,)),
                 lambda value: value.steps[0].update(mult=0.5),
                 lambda value: value.steps[0].update(mult=True)]
        for change in cases:
            value = deepcopy(self.skeleton)
            change(value)
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.save(value, order=self.order)

    def test_recomputed_skeleton_digest_cannot_conceal_implicit_coercions(self):
        self.save(self.skeleton, order=self.order)
        self.mutate(lambda data: data["payload"]["anchor"].update(q=5.0), redigest=True)
        with self.assertRaisesRegex(ValueError, "canonical"):
            self.load(order=self.order)

    def test_lift_reference_domain_and_metadata_are_validated(self):
        changes = [lambda value: value.entries.update({"": []}),
                   lambda value: value.anchor_index.update(g=-1),
                   lambda value: setattr(value, "samples_used", True),
                   lambda value: setattr(value, "holdout", -1),
                   lambda value: value.entries.update(h=[((0, 0), {(0, 0): F(1)}, (), "basis", "missing", ())]),
                   lambda value: value.entries.update(g=[((0, 0), {(0, 0): F(1)}, (), "basis", "h", ())]),
                   lambda value: value.entries.update(h=[((0, 0), {(0, 0): F(1)}, (), "input", -1, ())])]
        for change in changes:
            value = deepcopy(self.lift)
            change(value)
            with self.subTest(change=change), self.assertRaises((ValueError, KeyError)):
                self.save(value, order=self.order)

    def test_recorded_object_presentation_or_target_cannot_disagree(self):
        skeleton = deepcopy(self.skeleton)
        skeleton.presentation_hash = "other"
        with self.assertRaisesRegex(ValueError, "presentation"):
            self.save(skeleton, order=self.order)
        cofactor = deepcopy(self.cofactor)
        cofactor.target = "X_3"
        with self.assertRaises(ValueError):
            self.save(cofactor)
        cofactor.target = "X_2"
        with self.assertRaisesRegex(ValueError, "different target"):
            A.save_artifact(self.path, cofactor, self.pres, {})

    def test_empty_skeleton_and_lift_are_valid_data_for_a_zero_target(self):
        skel = Skeleton(self.pres.slug(), self.pres.content_hash(), "", (5, 7, 0))
        A.save_artifact(self.path, skel, self.pres, {}, order="length-lex-index-v1")
        self.assertEqual(A.load_artifact(self.path, self.pres, {}, order="length-lex-index-v1"), skel)
        lift = Lift(AtomSet([], []), {}, 0, 0, (5, 7, 0), {})
        A.save_artifact(self.path, lift, self.pres, {}, order="length-lex-index-v1")
        self.assertEqual(A.load_artifact(self.path, self.pres, {}, order="length-lex-index-v1"), lift)

    def test_serialization_errors_do_not_create_or_replace_files(self):
        self.path.write_text("existing")
        with self.assertRaises(ValueError):
            A.write_json_atomic(self.path, {"bad": float("nan")})
        self.assertEqual(self.path.read_text(), "existing")
        self.assertEqual(list(self.root.iterdir()), [self.path])

    def test_failed_promotion_preserves_old_file_and_removes_temporary_output(self):
        self.path.write_text("existing")
        with patch.object(A.os, "replace", side_effect=OSError("promotion interrupted")):
            with self.assertRaises(OSError):
                self.save()
        self.assertEqual(self.path.read_text(), "existing")
        self.assertEqual(list(self.root.iterdir()), [self.path])

    def test_failed_temporary_write_is_cleaned(self):
        opened = []
        make_temporary = A.tempfile.NamedTemporaryFile

        def remember_stream(*args, **kwargs):
            stream = make_temporary(*args, **kwargs)
            opened.append(stream)
            return stream

        with patch.object(A.os, "fsync", side_effect=OSError("disk error")), \
                patch.object(A.tempfile, "NamedTemporaryFile", side_effect=remember_stream):
            with self.assertRaises(OSError):
                self.save()
        self.assertTrue(opened[0].closed)
        self.assertFalse(self.path.exists())
        self.assertEqual(list(self.root.iterdir()), [])

    def test_nested_destination_is_created_and_output_is_stable_json(self):
        nested = self.root / "new" / "data.json"
        A.write_json_atomic(nested, {"b": 2, "a": "exact"})
        self.assertEqual(nested.read_text(), '{\n  "a": "exact",\n  "b": 2\n}\n')


if __name__ == "__main__":
    unittest.main()
