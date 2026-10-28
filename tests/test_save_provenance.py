"""Public output writers preserve data bytes and identify their producing sources."""
from contextlib import redirect_stdout
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

from zipper_algebra import Presentation
from zipper_algebra import msolve
from zipper_algebra.cofactor.atoms import AtomSet
from zipper_algebra.cofactor.cofactor_io import Cofactor
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.cofactor.skeleton import Skeleton


class TestPublicSaveProvenance(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.pres = Presentation(2, zero_from=2)
        self.atoms = AtomSet([])
        self.certificate = Cofactor.from_presentation(
            self.pres, "X_2", self.atoms,
            [((), "ZERO_2_1", (), {(0, 0): Fraction(1)}, (0,) * len(self.atoms))],
            producer="test")
        self.skeleton = Skeleton(self.pres.slug(), self.pres.content_hash(),
                                 "X_2", (5, 7, 0), [], [], 0)
        self.lift = Lift(self.atoms, {}, 0, 0, (5, 7, 0), {})

    def assert_tagged(self, path):
        record = json.loads(Path(str(path) + ".provenance.json").read_text())
        self.assertEqual(record["schema"], "zipper-output-provenance-v1")
        content = path.read_bytes()
        self.assertEqual(record["outputs"], {
            path.name: {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}})
        self.assertEqual(record["producer"]["schema"], "zipper-source-provenance-v1")
        self.assertIn("src/zipper_algebra/presentation.py", record["producer"]["source_files"])

    def writers(self):
        def standalone(path, *, provenance=True):
            with patch.object(sys, "argv", ["msolve", "--n", "2", "--out", str(path)]), \
                    redirect_stdout(io.StringIO()):
                return msolve.main(provenance=provenance)

        return [("atoms", self.atoms.save), ("lift", self.lift.save),
                ("skeleton", self.skeleton.save), ("certificate", self.certificate.save),
                ("presentation", lambda path, **kwargs: self.pres.emit_ms(path, 5, 7, 0, **kwargs)),
                ("standalone", standalone)]

    def test_producer_snapshot_precedes_serialization_for_every_public_writer(self):
        for name, writer in self.writers():
            with self.subTest(name=name):
                output = self.root / name
                snapshot = {"state": "before serialization"}

                def capture():
                    self.assertFalse(output.exists())
                    return snapshot

                def tag(path, *, producer):
                    self.assertEqual(Path(path), output)
                    self.assertTrue(output.is_file())
                    self.assertIs(producer, snapshot)

                with patch("zipper_algebra.provenance.capture_provenance", side_effect=capture) as collect, \
                        patch("zipper_algebra.provenance.write_output_provenance", side_effect=tag) as write:
                    writer(output)
                collect.assert_called_once_with()
                write.assert_called_once()

    def test_suppression_skips_source_capture_as_well_as_sidecar_creation(self):
        with patch("zipper_algebra.provenance.capture_provenance") as capture, \
                patch("zipper_algebra.provenance.write_output_provenance") as tag:
            for name, writer in self.writers():
                with self.subTest(name=name):
                    writer(self.root / name, provenance=False)
            capture.assert_not_called()
            tag.assert_not_called()

    def test_legacy_objects_keep_the_same_bytes_and_readers(self):
        for name, value in (("atoms.json", self.atoms), ("lift.json", self.lift),
                            ("skeleton.json", self.skeleton), ("certificate.txt", self.certificate)):
            with self.subTest(name=name):
                tagged, untagged = self.root / name, self.root / ("raw-" + name)
                value.save(tagged)
                value.save(untagged, provenance=False)
                self.assertEqual(tagged.read_bytes(), untagged.read_bytes())
                self.assert_tagged(tagged)
                self.assertFalse(Path(str(untagged) + ".provenance.json").exists())
        self.assertEqual(AtomSet.load(self.root / "atoms.json").to_json(), self.atoms.to_json())
        self.assertEqual(Lift.load(self.root / "lift.json"), self.lift)
        self.assertEqual(Skeleton.load(self.root / "skeleton.json"), self.skeleton)
        self.assertEqual(Cofactor.load(self.root / "certificate.txt", pres=self.pres).terms,
                         self.certificate.terms)

    def test_presentation_emission_preserves_solver_grammar_and_count(self):
        tagged, raw = self.root / "input.ms", self.root / "raw.ms"
        count = self.pres.emit_ms(tagged, 5, 7, 0)
        self.assertEqual(self.pres.emit_ms(raw, 5, 7, 0, provenance=False), count)
        self.assertEqual(count, len(self.pres.labels()))
        self.assertEqual(tagged.read_bytes(), raw.read_bytes())
        self.assertEqual(tagged.read_text().splitlines()[:2], ["x1,y1", "0"])
        self.assert_tagged(tagged)
        self.assertFalse(Path(str(raw) + ".provenance.json").exists())

    def test_standalone_emitter_tags_by_default_and_can_suppress_managed_writes(self):
        tagged, raw = self.root / "standalone.ms", self.root / "standalone-raw.ms"
        for path, enabled in ((tagged, True), (raw, False)):
            with patch.object(sys, "argv", ["msolve", "--n", "2", "--out", str(path)]), \
                    redirect_stdout(io.StringIO()):
                self.assertEqual(msolve.main(provenance=enabled), 0)
        self.assertEqual(tagged.read_bytes(), raw.read_bytes())
        self.assert_tagged(tagged)
        self.assertFalse(Path(str(raw) + ".provenance.json").exists())

    def test_failed_data_writes_do_not_create_provenance_records(self):
        invalid = self.root / "missing" / "output"
        with patch("zipper_algebra.provenance.write_output_provenance") as tag:
            for value in (self.atoms, self.lift, self.skeleton, self.certificate):
                with self.subTest(value=type(value).__name__), self.assertRaises(OSError):
                    value.save(invalid)
            with self.assertRaises(OSError):
                self.pres.emit_ms(invalid, 5, 7, 0)
            with patch.object(sys, "argv", ["msolve", "--n", "2", "--out", str(invalid)]), \
                    self.assertRaises(OSError):
                msolve.main()
            tag.assert_not_called()


if __name__ == "__main__":
    unittest.main()
