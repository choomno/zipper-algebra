"""Public output writers preserve data bytes and identify their producing sources."""
from contextlib import redirect_stdout
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import pytest
from unittest.mock import patch

from zipper_algebra import Presentation
from zipper_algebra import msolve
from zipper_algebra.cofactor.atoms import AtomSet
from zipper_algebra.cofactor.cofactor_io import Cofactor
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.cofactor.skeleton import Skeleton


class TestPublicSaveProvenance:
    @pytest.fixture(autouse=True)
    def resources(self, request):
        temporary = tempfile.TemporaryDirectory()
        request.addfinalizer(temporary.cleanup)
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

    def check_tagged_output(self, path):
        record = json.loads(Path(str(path) + ".provenance.json").read_text())
        assert record['schema'] == 'zipper-output-provenance-v1'
        content = path.read_bytes()
        assert record['outputs'] == {path.name: {'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()}}
        assert record['producer']['schema'] == 'zipper-source-provenance-v1'
        assert 'src/zipper_algebra/presentation.py' in record['producer']['source_files']

    def writers(self):
        def standalone(path, *, provenance=True):
            with patch.object(sys, "argv", ["msolve", "--n", "2", "--out", str(path)]), \
                    redirect_stdout(io.StringIO()):
                return msolve.main(provenance=provenance)

        return [("atoms", self.atoms.save), ("lift", self.lift.save),
                ("skeleton", self.skeleton.save), ("certificate", self.certificate.save),
                ("presentation", lambda path, **kwargs: self.pres.emit_ms(path, 5, 7, 0, **kwargs)),
                ("standalone", standalone)]

    def test_producer_snapshot_precedes_serialization_for_every_public_writer(self, subtests):
        for name, writer in self.writers():
            with subtests.test(name=name):
                output = self.root / name
                snapshot = {"state": "before serialization"}

                def capture():
                    assert not output.exists()
                    return snapshot

                def tag(path, *, producer):
                    assert Path(path) == output
                    assert output.is_file()
                    assert producer is snapshot

                with patch("zipper_algebra.provenance.capture_provenance", side_effect=capture) as collect, \
                        patch("zipper_algebra.provenance.write_output_provenance", side_effect=tag) as write:
                    writer(output)
                collect.assert_called_once_with()
                write.assert_called_once()

    def test_suppression_skips_source_capture_as_well_as_sidecar_creation(self, subtests):
        with patch("zipper_algebra.provenance.capture_provenance") as capture, \
                patch("zipper_algebra.provenance.write_output_provenance") as tag:
            for name, writer in self.writers():
                with subtests.test(name=name):
                    writer(self.root / name, provenance=False)
            capture.assert_not_called()
            tag.assert_not_called()

    def test_legacy_objects_keep_the_same_bytes_and_readers(self, subtests):
        for name, value in (("atoms.json", self.atoms), ("lift.json", self.lift),
                            ("skeleton.json", self.skeleton), ("certificate.txt", self.certificate)):
            with subtests.test(name=name):
                tagged, untagged = self.root / name, self.root / ("raw-" + name)
                value.save(tagged)
                value.save(untagged, provenance=False)
                assert tagged.read_bytes() == untagged.read_bytes()
                self.check_tagged_output(tagged)
                assert not Path(str(untagged) + '.provenance.json').exists()
        assert AtomSet.load(self.root / 'atoms.json').to_json() == self.atoms.to_json()
        assert Lift.load(self.root / 'lift.json') == self.lift
        assert Skeleton.load(self.root / 'skeleton.json') == self.skeleton
        assert Cofactor.load(self.root / 'certificate.txt', pres=self.pres).terms == self.certificate.terms

    def test_presentation_emission_preserves_solver_grammar_and_count(self):
        tagged, raw = self.root / "input.ms", self.root / "raw.ms"
        count = self.pres.emit_ms(tagged, 5, 7, 0)
        assert self.pres.emit_ms(raw, 5, 7, 0, provenance=False) == count
        assert count == len(self.pres.labels())
        assert tagged.read_bytes() == raw.read_bytes()
        assert tagged.read_text().splitlines()[:2] == ['x1,y1', '0']
        self.check_tagged_output(tagged)
        assert not Path(str(raw) + '.provenance.json').exists()

    def test_standalone_emitter_tags_by_default_and_can_suppress_managed_writes(self):
        tagged, raw = self.root / "standalone.ms", self.root / "standalone-raw.ms"
        for path, enabled in ((tagged, True), (raw, False)):
            with patch.object(sys, "argv", ["msolve", "--n", "2", "--out", str(path)]), \
                    redirect_stdout(io.StringIO()):
                assert msolve.main(provenance=enabled) == 0
        assert tagged.read_bytes() == raw.read_bytes()
        self.check_tagged_output(tagged)
        assert not Path(str(raw) + '.provenance.json').exists()

    def test_failed_data_writes_do_not_create_provenance_records(self, subtests):
        invalid = self.root / "missing" / "output"
        with patch("zipper_algebra.provenance.write_output_provenance") as tag:
            for value in (self.atoms, self.lift, self.skeleton, self.certificate):
                with subtests.test(value=type(value).__name__), pytest.raises(OSError):
                    value.save(invalid)
            with pytest.raises(OSError):
                self.pres.emit_ms(invalid, 5, 7, 0)
            with patch.object(sys, "argv", ["msolve", "--n", "2", "--out", str(invalid)]), \
                    pytest.raises(OSError):
                msolve.main()
            tag.assert_not_called()
