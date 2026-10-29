"""Source and input identities accompany experimental files without changing them."""
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from experiments.output import (capture_experiment_provenance, input_identity,
                                prepare_output_dir, read_input, save_pickle_output,
                                tag_solver_input, tag_solver_basis)
from experiments.survey import run as survey_run
from zipper_algebra.provenance import validate_provenance


ROOT = Path(__file__).resolve().parents[1]


def identity(content):
    return {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def sidecar(path):
    return json.loads(path.with_name(path.name + ".provenance.json").read_text())


class TestExperimentProvenance(unittest.TestCase):
    def test_recipe_snapshot_includes_family_and_shared_sources(self):
        producer = capture_experiment_provenance(survey_run.__file__)
        self.assertIs(validate_provenance(producer), producer)
        for name in ("experiments/survey/run.py", "experiments/survey/pointalg.py",
                     "experiments/survey/spec.py", "experiments/output.py",
                     "src/zipper_algebra/presentation.py"):
            self.assertEqual(producer["source_files"][name], identity((ROOT / name).read_bytes()))

    def test_pickle_format_and_prehashed_input_survive_in_place_changes(self):
        producer = capture_experiment_provenance(survey_run.__file__)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.pkl"
            path.write_bytes(b"original input")
            original = input_identity(path)
            payload = {"basis_complete": False, "words": [("x1",), ()]}
            save_pickle_output(path, payload, producer=producer, inputs=[original])
            self.assertEqual(path.read_bytes(), pickle.dumps(payload))
            record = sidecar(path)
            self.assertEqual(record["producer"], producer)
            self.assertEqual(record["inputs"], [{"name": "candidate.pkl", **identity(b"original input")}])
            self.assertEqual(record["outputs"], {path.name: identity(pickle.dumps(payload))})

    def test_failed_pickle_serialization_preserves_previous_output_and_sidecar(self):
        producer = capture_experiment_provenance(survey_run.__file__)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.pkl"
            save_pickle_output(path, {"old": True}, producer=producer)
            old, old_record = path.read_bytes(), sidecar(path)
            def fail_dump(payload, stream):
                stream.write(b"partial pickle")
                raise RuntimeError("serialization failed")
            with patch.object(pickle, "dump", side_effect=fail_dump):
                with self.assertRaisesRegex(RuntimeError, "serialization failed"):
                    save_pickle_output(path, {"new": True}, producer=producer)
            self.assertEqual(path.read_bytes(), old)
            self.assertEqual(sidecar(path), old_record)
            self.assertEqual(sorted(p.name for p in Path(directory).iterdir()),
                             ["candidate.pkl", "candidate.pkl.provenance.json"])

    def test_failed_pickle_replacement_preserves_previous_output_and_sidecar(self):
        producer = capture_experiment_provenance(survey_run.__file__)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidate.pkl"
            save_pickle_output(path, {"old": True}, producer=producer)
            old, old_record = path.read_bytes(), sidecar(path)
            with patch('experiments.output.os.replace', side_effect=OSError("replacement failed")):
                with self.assertRaisesRegex(OSError, "replacement failed"):
                    save_pickle_output(path, {"new": True}, producer=producer)
            self.assertEqual(path.read_bytes(), old)
            self.assertEqual(sidecar(path), old_record)
            self.assertEqual(sorted(p.name for p in Path(directory).iterdir()),
                             ["candidate.pkl", "candidate.pkl.provenance.json"])

    def test_existing_checkout_examples_keep_unrelated_files_visible(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "-q", str(root)], check=True)
            examples = root / "examples/x3"
            examples.mkdir(parents=True)
            source = examples / "human-example.txt"
            source.write_text("retain visibility\n")
            prepare_output_dir(examples)
            self.assertFalse((examples / ".gitignore").exists())
            visible = subprocess.check_output(["git", "-C", str(root), "status", "--porcelain",
                                               "--untracked-files=all"], text=True)
            self.assertIn("examples/x3/human-example.txt", visible)
            generated = root / "runs/existing"
            generated.mkdir(parents=True)
            prepare_output_dir(generated)
            self.assertEqual((generated / ".gitignore").read_text(), "*\n!.gitignore\n")
            new_output = examples / "new-run"
            prepare_output_dir(new_output)
            self.assertTrue((new_output / ".gitignore").is_file())

    @unittest.skipUnless(importlib.util.find_spec('sympy'), 'trace cache needs Sympy')
    def test_trace_cache_is_bound_to_representation_and_exact_consumed_bytes(self):
        import sympy as sp
        from experiments.variants.square_zero_traceform_qt import load_trace_cache
        producer = capture_experiment_provenance(survey_run.__file__)
        with tempfile.TemporaryDirectory() as directory:
            representation, cache = Path(directory) / "representation.pkl", Path(directory) / "trace.pkl"
            representation.write_bytes(b"current representation bytes")
            _, parent = read_input(representation)
            payload = {"standard": [()], "T": sp.Matrix([[3]])}
            save_pickle_output(cache, payload, producer=producer, inputs=[parent])
            content = cache.read_bytes()
            real_read, reads = Path.read_bytes, []
            def replace_after_read(path):
                result = real_read(path)
                if path == cache and not reads:
                    reads.append(path)
                    path.write_bytes(b"changed cache after read")
                return result
            with patch.object(Path, "read_bytes", replace_after_read):
                matrix, used = load_trace_cache(cache, [()], parent)
            self.assertEqual(matrix, sp.Matrix([[3]]))
            self.assertEqual(used, {"name": cache.name, **identity(content)})
            self.assertEqual(load_trace_cache(cache, [()], parent), (None, None))
            save_pickle_output(cache, payload, producer=producer, inputs=[parent])
            representation.write_bytes(b"different representation")
            _, changed_parent = read_input(representation)
            self.assertEqual(load_trace_cache(cache, [()], changed_parent), (None, None))

    def test_solver_sidecars_preserve_raw_grammar_and_bind_preexecution_input(self):
        producer = capture_experiment_provenance(survey_run.__file__)
        with tempfile.TemporaryDirectory() as directory:
            source, basis = Path(directory) / "input.ms", Path(directory) / "basis.txt"
            content = b"x1,y1\n0\nx1*y1-1\n"
            source.write_bytes(content)
            inputs = tag_solver_input(source, producer=producer)
            self.assertEqual(source.read_bytes(), content)
            source.write_bytes(b"subsequent edit")
            tag_solver_basis(basis, producer=producer, inputs=inputs)
            self.assertFalse(basis.with_name(basis.name + ".provenance.json").exists())
            basis.write_bytes(b"x1*y1-1\n")
            tag_solver_basis(basis, producer=producer, inputs=inputs)
            self.assertEqual(basis.read_bytes(), b"x1*y1-1\n")
            self.assertEqual(sidecar(basis)["inputs"], [{"name": source.name, **identity(content)}])

    def test_survey_json_has_provenance_without_changing_its_observation_payload(self):
        from experiments.survey import run
        from experiments.survey.spec import Spec
        observation = {"points": [], "coefficient_domain": "Q"}
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(run.ss, "semisimplicity", return_value=observation):
                result = run.run_spec(Spec("linear", 3), directory, directory,
                                      do_basis=False, char=0, log=lambda text: None)
            path = Path(directory) / "linear-n3-mul2.json"
            self.assertEqual(json.loads(path.read_text()), result)
            record = sidecar(path)
            self.assertEqual(record["outputs"], {path.name: identity(path.read_bytes())})
            self.assertIn("experiments/research/theta_min_poly.py", record["producer"]["source_files"])
            self.assertEqual(record["inputs"], [])

    @unittest.skipUnless(importlib.util.find_spec('sympy'), 'survey reconstruction needs Sympy')
    def test_survey_sibling_outputs_are_not_reported_as_input_dependencies(self):
        import sympy as sp
        from experiments.survey import reconstruct, run
        from experiments.survey.spec import Spec
        candidate = {"LMS": [("x1",)], "atoms": [], "D": sp.Integer(1),
                     "GEN": {("x1",): {(): sp.Integer(1)}}}
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(run.ss, "semisimplicity", return_value={"points": [], "coefficient_domain": "Q"}), \
                 patch.object(reconstruct, "reconstruct_gb", return_value=candidate), \
                 patch.object(reconstruct, "validate", return_value=(False, "candidate rejected")):
                result = run.run_spec(Spec("linear", 3), directory, directory, log=lambda text: None)
            self.assertFalse(result["basis_candidate"]["sample_check_passed"])
            summary = Path(directory) / "linear-n3-mul2.json"
            candidate_path = Path(result["basis_candidate"]["artifact"])
            self.assertEqual(sidecar(summary)["inputs"], [])
            self.assertEqual(sidecar(candidate_path)["inputs"], [])


if __name__ == "__main__":
    unittest.main()
