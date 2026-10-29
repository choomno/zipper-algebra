"""Provenance hashes identify the exact byte snapshots consumed by readers."""
from contextlib import redirect_stdout
from copy import deepcopy
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from unittest.mock import patch

from zipper_algebra import Presentation
from zipper_algebra.cofactor import artifacts, compose
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.cofactor_io import Cofactor
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.cofactor.skeleton import Skeleton


class TestInputSnapshots:
    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        self.root = tmp_path
        self.pres = Presentation(2, zero_from=2)
        self.atoms = AtomSet([Atom("q", {(1, 0): 1})], [1])
        self.cofactor = Cofactor.from_presentation(
            self.pres, "X_2", self.atoms,
            [((), "ZERO_2_1", (), {(0, 0): Fraction(1)}, (0,))],
            producer="snapshot-reader-test", notes=["line\vseparator"])

    def test_text_decoder_matches_file_readers_across_current_and_legacy_formats(self, subtests):
        path = self.root / "certificate.txt"
        self.cofactor.save(path, provenance=False)
        texts = [path.read_text(),
                 "# (q, t) = (-2, 3)\n# note: exact numeric input\n+1/2 * x1 * [0] * y1\n",
                 "# historical symbolic witness\n() * INV_1_p * () :: (1*q^1*t^0) / (q)\n"]
        for text in texts:
            for newline in ("\n", "\r\n", "\r"):
                raw = text.replace("\n", newline)
                path.write_bytes(raw.encode("utf-8"))
                from_file = Cofactor.load(path, pres=self.pres, atoms=self.atoms, target="X_2")
                with patch("builtins.open", side_effect=AssertionError("text decoder reopened a file")):
                    from_text = Cofactor.from_text(raw, pres=self.pres, atoms=self.atoms,
                                                  target="X_2", name=str(path))
                with subtests.test(format=text[:35], newline=repr(newline)):
                    assert (from_text.terms) == (from_file.terms)
                    assert (from_text.header) == (from_file.header)
        with pytest.raises(ValueError, match="<memory> has no presentation"):
            Cofactor.from_text(texts[1])

    def test_bound_json_decoder_checks_current_schema_without_reopening(self):
        path = self.root / "artifact.json"
        artifacts.save_artifact(path, self.cofactor, self.pres, "X_2")
        document = json.loads(path.read_bytes())
        loaded = artifacts.load_artifact(path, self.pres, "X_2", kind="cofactor", atoms=self.atoms)
        path.write_text("not the original artifact")
        with patch.object(Path, "read_bytes", side_effect=AssertionError("decoder reopened its source")):
            decoded = artifacts.artifact_from_json(document, self.pres, "X_2",
                                                   kind="cofactor", atoms=self.atoms)
        assert (decoded.terms) == (loaded.terms)
        assert (decoded.header) == (loaded.header)
        retired = deepcopy(document)
        retired["schema"] = "zipper-mathematical-artifact-v1"
        del retired["provenance"]
        path.write_text(json.dumps(retired))
        with pytest.raises(ValueError, match="schema"):
            artifacts.load_artifact(path, self.pres, "X_2")
        with pytest.raises(ValueError, match="schema"):
            artifacts.artifact_from_json(retired, self.pres, "X_2")
        document["payload"]["terms"][0]["num"] = [[[0, 0], "2"]]
        with pytest.raises(ValueError, match="digest"):
            artifacts.artifact_from_json(document, self.pres, "X_2")

    def composition_inputs(self, bound):
        index = self.pres.labels().index("ZERO_2_1")
        lift = Lift(self.atoms, {"g": [((1,), {(0, 0): Fraction(1)}, (), "input", index, ())]},
                    3, 1, (5, 7, 0), {"g": 0})
        skeleton = Skeleton(self.pres.slug(), self.pres.content_hash(), "X_2", (5, 7, 0),
                            [{"basis_key": "g", "left": (), "right": (), "mult": Fraction(5)}],
                            basis_keys=["g"])
        order = {"solver": "degree-lex", "reduction": "length-lex-index-v1"}
        paths = [self.root / "skeleton.json", self.root / "lift.json", self.root / "multipliers.json"]
        if bound:
            artifacts.save_artifact(paths[0], skeleton, self.pres, "X_2", order=order)
            artifacts.save_artifact(paths[1], lift, self.pres, "X_2", order=order)
        else:
            paths[0].write_text(json.dumps(skeleton.to_json()))
            paths[1].write_text(json.dumps(lift.to_json()))
        multipliers = [({(1, 0): Fraction(1)}, (0,))]
        paths[2].write_text(json.dumps(compose.multipliers_to_json(
            multipliers, skeleton, self.atoms, self.pres)))
        return paths, skeleton, multipliers

    def test_multiplier_decoder_consumes_its_snapshot(self):
        paths, skeleton, expected = self.composition_inputs(bound=False)
        document = json.loads(paths[2].read_bytes())
        assert (compose.load_multipliers(paths[2], skeleton, self.atoms, self.pres)) == (expected)
        paths[2].write_text("changed after the snapshot")
        assert (compose.multipliers_from_json(document, skeleton, self.atoms, self.pres)) == (expected)
        document["target"] = "different"
        with pytest.raises(ValueError, match="different presentation or target"):
            compose.multipliers_from_json(document, skeleton, self.atoms, self.pres)

    def test_direct_composition_reader_reads_each_artifact_once(self):
        paths, skeleton, _ = self.composition_inputs(bound=True)
        reads = dict.fromkeys(paths[:2], 0)
        actual_read = Path.read_bytes

        def count_read(path):
            if path in reads:
                reads[path] += 1
            return actual_read(path)

        arguments = SimpleNamespace(skeleton=str(paths[0]), lift=str(paths[1]),
                                    target="X_2", solver_order="degree-lex")
        with patch.object(Path, "read_bytes", count_read):
            fitted, decoded = compose._composition_inputs(arguments, self.pres)
        assert (reads) == (dict.fromkeys(paths[:2], 1))
        assert (decoded) == (skeleton)
        assert (set(fitted.entries)) == ({"g"})

    def test_composition_hashes_and_decodes_the_same_bytes_after_source_files_change(self, subtests):
        for bound in (False, True):
            paths, _, _ = self.composition_inputs(bound)
            snapshots = {path: path.read_bytes() for path in paths}
            reads = {path: 0 for path in paths}
            actual_read = Path.read_bytes

            def change_after_read(path):
                data = actual_read(path)
                if path in snapshots:
                    reads[path] += 1
                    path.write_bytes(b'{"changed":"after capture"}')
                return data

            output = self.root / "composed.txt"
            arguments = ["--n", "2", "--zero-from", "2", "--target", "X_2", "--json",
                         "--skeleton", str(paths[0]), "--lift", str(paths[1]),
                         "--multipliers", str(paths[2]), "--out", str(output)]
            stream = io.StringIO()
            with patch.object(Path, "read_bytes", change_after_read), redirect_stdout(stream):
                assert (compose.main(arguments)) == (0), (stream.getvalue())
            result = json.loads(stream.getvalue())
            sidecar = json.loads(output.with_name(output.name + ".provenance.json").read_text())
            expected = [{"name": path.name, "bytes": len(snapshots[path]),
                         "sha256": hashlib.sha256(snapshots[path]).hexdigest()} for path in paths]
            with subtests.test(bound=bound):
                assert (reads) == (dict.fromkeys(paths, 1))
                assert (result["inputs"]) == (expected)
                assert (sidecar["inputs"]) == (expected)
                assert (result["pass"])
                assert (Cofactor.load(output, pres=self.pres).terms) == ([((), "ZERO_2_1", (), {(1, 0): Fraction(1)}, (1,))])
                assert not (list(self.root.glob(".*.provenance.json")))
                assert not (list(self.root.glob("*.tmp*")))

    def test_snapshot_composition_rejects_a_wellformed_different_target(self):
        pres = Presentation(3)
        skeleton = Skeleton(pres.slug(), pres.content_hash(), "X_2", (5, 7, 0))
        fitted = Lift(AtomSet([]), {}, 0, 0, (5, 7, 0), {})
        paths = [self.root / "skeleton.json", self.root / "lift.json", self.root / "multipliers.json"]
        for path, data in zip(paths, [skeleton.to_json(), fitted.to_json(),
                                    compose.multipliers_to_json([], skeleton, fitted.atoms, pres)]):
            path.write_text(json.dumps(data))
        output = self.root / "wrong-target.txt"
        arguments = ["--n", "3", "--target", "X_3", "--json", "--skeleton", str(paths[0]),
                     "--lift", str(paths[1]), "--multipliers", str(paths[2]), "--out", str(output)]
        stream = io.StringIO()
        with redirect_stdout(stream):
            assert (compose.main(arguments)) == (2)
        assert ("requested target differs") in (json.loads(stream.getvalue())["error"])
        assert not (output.exists())
