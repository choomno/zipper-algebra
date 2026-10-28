"""Process and artifact boundary tests using an explicitly synthetic producer."""
from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
from threading import Barrier
import unittest
from unittest.mock import patch

from zipper_algebra import Presentation
from zipper_algebra.execution import ArtifactError, ExecutionError, RunConfig, run_point
from zipper_algebra.execution.artifacts import (
    ARTIFACT_NAMES, canonical_json, digest_bytes, file_identity, load_record,
    publish, validate_proofs,
)
from zipper_algebra.execution.config import DEFAULT_FLAGS
from zipper_algebra.zipper_in_ms import NumericRing


# This writes a known subset of the input relations with tautological input
# proofs. It exercises process control and validation; it is not a GB solver.
PRODUCER = r'''
import argparse
import json
from pathlib import Path
import re
import sys
import time

p = argparse.ArgumentParser()
p.add_argument("-i", required=True)
p.add_argument("-o", required=True)
p.add_argument("-p", required=True)
p.add_argument("-m", type=int)
p.add_argument("-d", type=int)
p.add_argument("--mode", default="success")
p.add_argument("--counter")
a, _ = p.parse_known_args()
if a.counter:
    counter = Path(a.counter)
    n = int(counter.read_text()) if counter.exists() else 0
    counter.write_text(str(n + 1))
if a.mode == "timeout":
    time.sleep(10)
if a.mode == "exit":
    print("test producer rejected its input", file=sys.stderr)
    sys.exit(7)
relations = [line.rstrip(",") for line in Path(a.i).read_text().splitlines()[2:]]
spaced = [re.sub(r"(?<!^)([+-])", r" \1 ", line) for line in relations]
basis = spaced[0] + "\n"
proof = "[i0]\n"
if a.mode == "all":
    basis = "\n".join(spaced) + "\n"
    proof = "\n".join("[i%d]" % i for i in range(len(spaced))) + "\n"
elif a.mode == "forward":
    basis = (spaced[0] + "\n") * 2
    proof = "[1]\n[i0]\n"
elif a.mode == "contexts":
    basis = "1/2*x1*x1*y1*y1 - 1/2*x1*y1\n"
    proof = "1/2*x1*[i0]*y1\n"
elif a.mode == "missing":
    Path(a.o).write_text(basis)
    sys.exit(0)
elif a.mode == "empty":
    basis = proof = ""
elif a.mode == "empty-proof":
    proof = "\n"
elif a.mode == "count":
    proof += "[i0]\n"
elif a.mode == "identity":
    basis = spaced[0] + " + 1\n"
elif a.mode == "cycle":
    proof = "[0]\n"
elif a.mode == "zero-reference":
    proof = "[i0] +0*[i999]\n"
elif a.mode == "cancelled-reference":
    proof = "[i0] +[i999] -[i999]\n"
elif a.mode == "basis-reference":
    proof = "[999]\n"
elif a.mode == "dangling":
    proof = "[i0] +\n"
elif a.mode == "double-standalone":
    proof = "[i0] + - [i0]\n"
elif a.mode == "double-sign":
    proof = "[i0] + +[i0]\n"
elif a.mode == "missing-sign":
    proof = "[i0] [i0]\n"
elif a.mode == "malformed":
    proof = "[i0]junk\n"
elif a.mode == "zero-basis":
    basis = "0\n"
    proof = "0*[i0]\n"
elif a.mode == "zero-denominator":
    basis = "1/0*x1*y1\n"
elif a.mode == "comma":
    basis = basis.rstrip() + ",\n"
    proof = proof.rstrip() + ",\n"
elif a.mode == "standalone":
    basis = "+ " + basis
    proof = "+ [i0]\n"
elif a.mode == "blank-basis":
    basis = "\n" + basis + "\n"
elif a.mode == "mutate-input":
    Path(a.i).write_text("not the canonical input\n")
elif a.mode == "extra":
    Path("unrequested-solver-output.tmp").write_text("not a declared artifact")
Path(a.o).write_text(basis)
if a.mode == "bad-bytes":
    Path(a.p).write_bytes(b"\xff")
else:
    Path(a.p).write_text(proof)
print(json.dumps({"iterations": a.m, "degree": a.d}))
'''


class RunnerFixture(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="zipper-runner-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.producer = self.root / "producer.py"
        self.producer.write_text(PRODUCER)
        self.counter = self.root / "counter"
        self.pres = Presentation(2)
        self.point = (5, 7, 0)

    def config(self, mode="success", **kwargs):
        values = dict(
            binary=Path(sys.executable), output_root=self.root / "output",
            timeout_seconds=3, max_iterations=13, command_prefix=(self.producer,),
            flags=(*DEFAULT_FLAGS, "--mode", mode, "--counter", str(self.counter)),
            solver_revision="deterministic test producer; not f4ncgb",
        )
        values.update(kwargs)
        return RunConfig(**values)

    def run_solver(self, mode="success", **kwargs):
        return run_point(self.pres, "X_2", self.point, self.config(mode, **kwargs))

    def marker(self, run, root=None):
        return (root or self.root / "output") / run.key / "record.json"

    def restore_record(self, run, record):
        self.marker(run).write_bytes(canonical_json(record))


class TestRunConfiguration(RunnerFixture):
    def test_configuration_is_immutable_and_has_explicit_limits(self):
        config = self.config(max_degree=9)
        with self.assertRaises(FrozenInstanceError):
            config.max_degree = 10
        self.assertEqual(config.arguments()[-4:], ("-m", "13", "-d", "9"))
        self.assertEqual(config.binary, Path(sys.executable).resolve())
        self.assertIsNone(config.cache_root)
        self.assertEqual(RunConfig(Path(sys.executable), self.root).flags, DEFAULT_FLAGS)

    def test_invalid_configuration_is_refused(self):
        invalid = [
            {"binary": None}, {"output_root": None}, {"binary": self.root / "absent"},
            *({"timeout_seconds": value} for value in (True, "3", 0, -1, float("inf"), float("nan"))),
            *({"max_iterations": value} for value in (None, True, 0, -1, 1.5)),
            *({"max_degree": value} for value in (True, 0, -1, 1.5)),
            {"flags": []}, {"flags": ("",)}, {"flags": (3,)},
            *({"flags": (value,)} for value in ("-i", "--input=file", "-ooutput",
                                               "--proof=file", "-m10", "--maxdeg=3")),
            {"solver_revision": ""}, {"solver_revision": None},
            {"monomial_order": "unknown"}, {"command_prefix": []},
            {"command_prefix": (self.root / "absent",)},
        ]
        for change in invalid:
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.config(**change)

    def test_named_target_and_characteristic_are_explicit(self):
        for point in ([5, 7, 0], (5, 7), (True, 7, 0), (5, 7, 5)):
            with self.subTest(point=point), self.assertRaises(ValueError):
                run_point(self.pres, "X_2", point, self.config())
        for target in ((2, 1), "X_3"):
            with self.subTest(target=target), self.assertRaises(ValueError):
                run_point(self.pres, target, self.point, self.config())
        self.assertFalse(self.counter.exists())


class TestSolverExecution(RunnerFixture):
    def test_exact_input_and_all_proof_identities_are_retained(self):
        run = self.run_solver("all", max_degree=11)
        self.assertFalse(run.hit)
        self.assertGreaterEqual(run.seconds, 0)
        expected = self.root / "expected.ms"
        self.pres.emit_ms(str(expected), *self.point)
        self.assertEqual((run.directory / "input.ms").read_bytes(), expected.read_bytes())
        self.assertEqual(run.record["status"], "proof-validated")
        self.assertIs(run.record["basis_complete"], False)
        self.assertEqual(run.record["basis_entries"], self.pres.n_relations())
        self.assertEqual(run.record["request"]["presentation"]["legacy_hash"],
                         self.pres.content_hash())
        self.assertEqual(run.record["request"]["presentation"]["coefficient_domain"], "Q[q,t]")
        self.assertEqual(run.record["request"]["point"], [5, 7, 0])
        self.assertEqual(run.record["request"]["target"], "X_2")
        self.assertEqual(run.record["request"]["monomial_order"], "degree-lex")
        self.assertEqual(run.basis_path, run.directory / "basis.txt")
        self.assertEqual(run.proof_path, run.directory / "proof.txt")
        stdout = json.loads((run.directory / "stdout.txt").read_text())
        self.assertEqual(stdout, {"iterations": 13, "degree": 11})

    def test_forward_references_and_supported_output_spacing(self):
        for mode in ("forward", "contexts", "comma", "standalone", "blank-basis", "extra"):
            with self.subTest(mode=mode):
                run = self.run_solver(mode)
                self.assertEqual(set(p.name for p in run.directory.iterdir()), set(ARTIFACT_NAMES))
                self.assertEqual(run.record["basis_entries"], 2 if mode == "forward" else 1)

    def test_invalid_and_partial_outputs_never_have_completion_records(self):
        for mode in ("missing", "empty", "empty-proof", "count", "identity", "cycle",
                     "zero-reference", "cancelled-reference", "basis-reference",
                     "dangling", "double-standalone", "double-sign", "missing-sign",
                     "malformed", "zero-basis", "zero-denominator", "bad-bytes", "mutate-input"):
            output = self.root / mode
            with self.subTest(mode=mode), self.assertRaises(ArtifactError):
                self.run_solver(mode, output_root=output)
            self.assertFalse(list(output.rglob("record.json")))
            self.assertFalse(list(output.glob(".solver-*")))

    def test_nonzero_exit_and_timeout_are_distinct_process_failures(self):
        with self.assertRaisesRegex(ExecutionError, "code 7.*test producer"):
            self.run_solver("exit")
        with self.assertRaisesRegex(ExecutionError, "exceeded"):
            self.run_solver("timeout", timeout_seconds=0.03)
        self.assertFalse(list((self.root / "output").rglob("record.json")))

    def test_missing_binary_after_configuration_is_a_start_failure(self):
        executable = self.root / "copied-executable"
        shutil.copy2(sys.executable, executable)
        config = self.config(binary=executable)
        executable.unlink()
        # Identity hashing catches this before a process is attempted.
        with self.assertRaises(FileNotFoundError):
            run_point(self.pres, "X_2", self.point, config)
        nonexecutable = self.root / "not-an-executable"
        nonexecutable.write_bytes(b"invalid executable contents")
        with self.assertRaisesRegex(ExecutionError, "could not start"):
            self.run_solver(binary=nonexecutable)

    def test_output_cache_and_refresh_do_not_conflate_computation_with_completeness(self):
        config = self.config(cache_root=self.root / "cache")
        first = run_point(self.pres, "X_2", self.point, config)
        hit = run_point(self.pres, "S[2][1]", self.point, config)
        self.assertTrue(hit.hit)
        self.assertEqual(first.key, hit.key)
        self.assertEqual(self.counter.read_text(), "1")
        other = replace(config, output_root=self.root / "other-output")
        copied = run_point(self.pres, "X_2", self.point, other)
        self.assertTrue(copied.hit)
        self.assertEqual(copied.key, first.key)
        self.assertEqual(copied.basis_path.read_bytes(), first.basis_path.read_bytes())
        self.assertEqual(self.counter.read_text(), "1")
        refreshed = run_point(self.pres, "X_2", self.point, config, refresh=True)
        self.assertFalse(refreshed.hit)
        self.assertEqual(self.counter.read_text(), "2")
        self.assertEqual(refreshed.directory, first.directory)
        self.assertIs(refreshed.record["basis_complete"], False)

    def test_changed_inputs_binary_source_flags_and_limits_change_keys(self):
        config = self.config()
        first = run_point(self.pres, "X_2", self.point, config)
        points = run_point(self.pres, "X_2", (13, 17, 0), config)
        different = run_point(Presentation(2, untwist="none"), "X_2", self.point, config)
        limit = run_point(self.pres, "X_2", self.point, replace(config, max_iterations=14))
        flag = self.run_solver("comma")
        self.producer.write_text(PRODUCER + "\n# changed executable source\n")
        changed = run_point(self.pres, "X_2", self.point, config)
        self.assertEqual(len({r.key for r in (first, points, different, limit, flag, changed)}), 6)


class TestCompletedRunRecords(RunnerFixture):
    def test_absent_marker_is_not_a_completed_run(self):
        self.assertIsNone(load_record(self.root, "a"*64, {}, [], []))

    def test_record_schema_request_and_metadata_tampering_are_rejected(self):
        run = self.run_solver()
        original = json.loads(self.marker(run).read_text())
        changes = [
            lambda d: d.update(schema="unknown"),
            lambda d: d.update(status="verified"),
            lambda d: d.update(basis_complete=True),
            lambda d: d.update(basis_entries=True),
            lambda d: d.update(basis_entries=99),
            lambda d: d.update(generation=None),
            lambda d: d.update(generation="short"),
            lambda d: d.update(generation="z"*64),
            lambda d: d.update(files={}),
            lambda d: d.update(generation="0"*64),
            lambda d: d["request"].update(target="X_3"),
            lambda d: d["request"].update(point=[5, 7, False]),
            lambda d: d.pop("schema"),
        ]
        for mutate in changes:
            record = json.loads(json.dumps(original))
            mutate(record)
            self.restore_record(run, record)
            with self.subTest(record=record), self.assertRaises(ArtifactError):
                self.run_solver()
        self.marker(run).write_text("{not json")
        with self.assertRaises(ArtifactError):
            self.run_solver()
        self.restore_record(run, original)
        self.assertTrue(self.run_solver().hit)

    def test_output_byte_tampering_is_rejected_and_refresh_repairs_it(self):
        run = self.run_solver()
        run.proof_path.write_text("[i1]\n")
        with self.assertRaisesRegex(ArtifactError, "hash differs"):
            self.run_solver()
        fixed = run_point(self.pres, "X_2", self.point, self.config(), refresh=True)
        self.assertEqual(fixed.proof_path.read_text(), "[i0]\n")
        fixed.proof_path.unlink()
        with self.assertRaises(ArtifactError):
            self.run_solver()
        fixed = run_point(self.pres, "X_2", self.point, self.config(), refresh=True)
        self.assertEqual(fixed.proof_path.read_text(), "[i0]\n")

    def _replace_generation(self, run, filename, content):
        """Model an internally consistent but semantically false external record."""
        record = json.loads(self.marker(run).read_text())
        new = self.root / "changed"
        shutil.copytree(run.directory, new)
        (new / filename).write_bytes(content)
        record["files"] = {name: file_identity(new / name) for name in ARTIFACT_NAMES}
        generation = digest_bytes(canonical_json(record["files"]))
        destination = self.marker(run).parent / generation
        shutil.move(new, destination)
        record["generation"] = generation
        self.restore_record(run, record)

    def test_a_rehashed_false_proof_still_fails_numeric_validation(self):
        run = self.run_solver()
        self._replace_generation(run, "proof.txt", b"[i1]\n")
        with self.assertRaisesRegex(ArtifactError, "fails its input-ideal identity"):
            self.run_solver()

    def test_a_rehashed_different_input_is_not_the_requested_input(self):
        run = self.run_solver()
        self._replace_generation(run, "input.ms", b"x1,y1\n0\nx1*y1-2\n")
        with self.assertRaisesRegex(ArtifactError, "saved input differs"):
            self.run_solver()

    def test_corrupt_shared_cache_is_refused_for_a_new_output(self):
        config = self.config(cache_root=self.root / "cache")
        run = run_point(self.pres, "X_2", self.point, config)
        self.marker(run, self.root / "cache").write_text("{}")
        with self.assertRaises(ArtifactError):
            run_point(self.pres, "X_2", self.point,
                      replace(config, output_root=self.root / "another"))

    def test_interrupted_publication_preserves_a_prior_record(self):
        run = self.run_solver()
        before = self.marker(run).read_bytes()
        stage = self.root / "new-generation"
        shutil.copytree(run.directory, stage)
        (stage / "stdout.txt").write_text("a new valid run")
        real_replace = os.replace
        def interrupt_marker(source, destination):
            if Path(destination).name == "record.json":
                raise OSError("interrupted before selecting the new generation")
            return real_replace(source, destination)
        with patch("zipper_algebra.execution.artifacts.os.replace",
                   side_effect=interrupt_marker):
            with self.assertRaisesRegex(OSError, "interrupted"):
                publish(self.root / "output", run.key, stage,
                        run.record["request"], run.record["basis_entries"])
        self.assertEqual(self.marker(run).read_bytes(), before)
        self.assertTrue(self.run_solver().hit)

    def test_concurrent_identical_generations_both_publish_successfully(self):
        run = self.run_solver()
        destination = self.root / "concurrent-output"
        collision = Barrier(2)
        failures = []
        real_replace = os.replace

        def simultaneous_replace(source, final):
            if Path(source).name == "generation":
                collision.wait(timeout=10)
                try:
                    return real_replace(source, final)
                except OSError as error:
                    failures.append(error)
                    raise
            return real_replace(source, final)

        def publish_one():
            return publish(destination, run.key, run.directory,
                           run.record["request"], run.record["basis_entries"])

        with patch("zipper_algebra.execution.artifacts.os.replace",
                   side_effect=simultaneous_replace):
            with ThreadPoolExecutor(max_workers=2) as workers:
                futures = [workers.submit(publish_one) for _ in range(2)]
                results = [future.result(timeout=15) for future in futures]
        self.assertEqual(len(failures), 1)
        self.assertEqual(results[0], results[1])
        names, inputs, _, _ = self.pres.build(NumericRing(*self.point))
        loaded = load_record(destination, run.key, run.record["request"], inputs, names)
        self.assertEqual(loaded, results[0])

    def test_publication_does_not_hide_an_unrelated_rename_failure(self):
        run = self.run_solver()
        destination = self.root / "failed-publication"
        failure = OSError("storage unavailable")
        with patch("zipper_algebra.execution.artifacts.os.replace", side_effect=failure):
            with self.assertRaises(OSError) as caught:
                publish(destination, run.key, run.directory,
                        run.record["request"], run.record["basis_entries"])
        self.assertIs(caught.exception, failure)
        self.assertFalse((destination / run.key / "record.json").exists())

    def test_publication_does_not_accept_a_different_generation_after_failure(self):
        run = self.run_solver()
        destination = self.root / "different-publication"
        failure = OSError("rename failed after another writer appeared")

        def different_writer(source, final):
            Path(final).mkdir()
            (Path(final) / "input.ms").write_bytes(b"different bytes")
            raise failure

        with patch("zipper_algebra.execution.artifacts.os.replace", side_effect=different_writer):
            with self.assertRaises(OSError) as caught:
                publish(destination, run.key, run.directory,
                        run.record["request"], run.record["basis_entries"])
        self.assertIs(caught.exception, failure)
        self.assertFalse((destination / run.key / "record.json").exists())

    def test_unselected_generation_is_not_resumed(self):
        run = self.run_solver()
        self.marker(run).unlink()
        rerun = self.run_solver()
        self.assertFalse(rerun.hit)
        self.assertEqual(self.counter.read_text(), "2")
        self.assertEqual(rerun.directory, run.directory)

    def test_external_slot_and_generation_symlinks_are_rejected(self):
        run = self.run_solver()
        external = self.root / "external-generation"
        shutil.copytree(run.directory, external)
        shutil.rmtree(run.directory)
        try:
            run.directory.symlink_to(external, target_is_directory=True)
        except OSError:
            self.skipTest("directory symlinks are unavailable on this platform")
        with self.assertRaisesRegex(ArtifactError, "leaves"):
            self.run_solver()
        with self.assertRaisesRegex(ArtifactError, "leaves"):
            run_point(self.pres, "X_2", self.point, self.config(), refresh=True)
        run.directory.unlink()
        shutil.copytree(external, run.directory)
        slot = self.marker(run).parent
        moved = self.root / "external-slot"
        shutil.move(slot, moved)
        slot.symlink_to(moved, target_is_directory=True)
        with self.assertRaisesRegex(ArtifactError, "leaves"):
            self.run_solver()
        with self.assertRaisesRegex(ArtifactError, "leaves"):
            run_point(self.pres, "X_2", self.point, self.config(), refresh=True)

    def test_file_hashing_streams_more_than_one_chunk(self):
        path = self.root / "large-artifact"
        content = b"proof output" * 20000
        path.write_bytes(content)
        self.assertEqual(file_identity(path),
                         {"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()})


if __name__ == "__main__":
    unittest.main()
