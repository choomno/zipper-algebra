"""Process and artifact boundary tests using an explicitly synthetic producer."""
from __future__ import annotations

from contextlib import closing
from dataclasses import FrozenInstanceError, replace
import hashlib
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import pytest
from unittest.mock import patch

from zipper_algebra import Presentation
from zipper_algebra.execution import ArtifactError, ExecutionError, RunConfig, run_point
from zipper_algebra.execution import pipeline
from zipper_algebra.execution.artifacts import (
    ARTIFACT_NAMES, canonical_json, digest_bytes, file_identity,
    prepare_output_root,
)
from zipper_algebra.execution.config import DEFAULT_FLAGS
from zipper_algebra.execution import store
from zipper_algebra.cofactor.proof import Basis, Proof
from zipper_algebra.cofactor.artifacts import target_identity
from zipper_algebra.execution.runner import run_basis_point
from zipper_algebra.provenance import capture_provenance


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
    Path(a.o).write_text(basis, newline="\n")
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
    Path(a.i).write_text("not the canonical input\n", newline="\n")
elif a.mode == "extra":
    Path("unrequested-solver-output.tmp").write_text("not a declared artifact")
Path(a.o).write_text(basis, newline="\n")
if a.mode == "bad-bytes":
    Path(a.p).write_bytes(b"\xff")
else:
    Path(a.p).write_text(proof, newline="\n")
print(json.dumps({"iterations": a.m, "degree": a.d}))
'''


class RunnerFixture:
    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        self.root = tmp_path
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

    def replace_record(self, run, record, database=None):
        with closing(sqlite3.connect(database or run.store)) as connection, connection:
            connection.execute("UPDATE generations SET record=? WHERE key=? AND generation=?",
                               (canonical_json(record), run.key, run.record["generation"]))

    def alter_file(self, run, name, content):
        with closing(sqlite3.connect(run.store)) as connection, connection:
            connection.execute("UPDATE artifacts SET content=? WHERE key=? AND generation=? AND name=?",
                               (content, run.key, run.record["generation"], name))


class TestRunConfiguration(RunnerFixture):
    def test_configuration_is_immutable_and_has_explicit_limits(self):
        config = self.config(max_degree=9)
        with pytest.raises(FrozenInstanceError):
            config.max_degree = 10
        assert (config.arguments()[-4:]) == (("-m", "13", "-d", "9"))
        assert (config.binary) == (Path(sys.executable).resolve())
        assert (config.cache_root) is None
        assert (RunConfig(Path(sys.executable), self.root).flags) == (DEFAULT_FLAGS)

    def test_invalid_configuration_is_refused(self, subtests):
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
            with subtests.test(change=change), pytest.raises(ValueError):
                self.config(**change)

    def test_named_target_and_characteristic_are_explicit(self, subtests):
        for point in ([5, 7, 0], (5, 7), (True, 7, 0), (5, 7, 5)):
            with subtests.test(point=point), pytest.raises(ValueError):
                run_point(self.pres, "X_2", point, self.config())
        for target in ((2, 1), "X_3"):
            with subtests.test(target=target), pytest.raises(ValueError):
                run_point(self.pres, target, self.point, self.config())
        assert not (self.counter.exists())
        assert not ((self.root / "output").exists())


class TestOutputRootProtection(RunnerFixture):
    def test_custom_output_and_cache_roots_ignore_generated_artifacts_in_git(self, subtests):
        subprocess.run(["git", "init", "--quiet", str(self.root)], check=True)
        config = self.config(output_root=self.root / "custom-computation",
                             cache_root=self.root / "custom-shared-cache")
        run = run_point(self.pres, "X_2", self.point, config)
        for root in (config.output_root, config.cache_root):
            with subtests.test(root=root):
                assert ((root / ".gitignore").read_text()) == ("*\n!.gitignore\n")
                assert (subprocess.run(
                    ["git", "check-ignore", "--quiet", str(root / "samples.sqlite3")],
                    cwd=self.root).returncode) == (0)
                assert (subprocess.run(
                    ["git", "check-ignore", "--quiet", str(root / ".gitignore")],
                    cwd=self.root).returncode) == (1)
        assert (set(run.files)) == (set(ARTIFACT_NAMES))
        assert (set(path.name for path in config.output_root.iterdir())) == ({".gitignore", "samples.sqlite3"})

    def test_reused_runs_prepare_missing_output_and_cache_markers(self):
        config = self.config(cache_root=self.root / "cache")
        first = run_point(self.pres, "X_2", self.point, config)
        for root in (config.output_root, config.cache_root):
            (root / ".gitignore").unlink()
        reused = run_point(self.pres, "X_2", self.point, config)
        assert (reused.hit)
        assert (reused.store) == (first.store)
        for root in (config.output_root, config.cache_root):
            assert ((root / ".gitignore").read_text()) == ("*\n!.gitignore\n")
        other = replace(config, output_root=self.root / "from-cache")
        copied = run_point(self.pres, "X_2", self.point, other)
        assert (copied.hit)
        assert ((other.output_root / ".gitignore").read_text()) == ("*\n!.gitignore\n")
        assert (self.counter.read_text()) == ("1")

    def test_existing_custom_ignore_policy_survives_repeated_preparation(self):
        config = self.config(cache_root=self.root / "cache")
        custom = b"# project-specific output rules\n*.txt\n"
        for root in (config.output_root, config.cache_root):
            root.mkdir()
            (root / ".gitignore").write_bytes(custom)
        run_point(self.pres, "X_2", self.point, config)
        run_point(self.pres, "X_2", self.point, config)
        for root in (config.output_root, config.cache_root):
            assert ((root / ".gitignore").read_bytes()) == (custom)

    def test_invalid_root_or_marker_propagates_filesystem_errors(self):
        occupied = self.root / "occupied"
        occupied.write_text("existing file")
        with pytest.raises(FileExistsError):
            prepare_output_root(occupied)
        root = self.root / "blocked-marker"
        (root / ".gitignore").mkdir(parents=True)
        with pytest.raises(OSError):
            prepare_output_root(root)
        assert ((root / ".gitignore").is_dir())

    def test_pipeline_prepares_root_and_sample_markers_before_execution(self):
        config = self.config("exit", output_root=self.root / "custom-pipeline")
        result = pipeline.run(self.pres, "X_2", self.point, config,
                              start=1, step=1, maxpts=1, verbose=False)
        assert (result["execution_status"]) == ("failed")
        for root in (config.output_root,):
            assert ((root / ".gitignore").read_text()) == ("*\n!.gitignore\n")


class TestSolverExecution(RunnerFixture):
    def test_in_memory_parsers_preserve_file_line_and_whitespace_semantics(self):
        basis_text = "x1*y1 -1\v +x1\r\nx1\r"
        proof_text = "[i0]\v +[i1]\r\n\r"
        basis_path, proof_path = self.root / "basis.txt", self.root / "proof.txt"
        basis_path.write_bytes(basis_text.encode())
        proof_path.write_bytes(proof_text.encode())
        names = self.pres.var_names()
        actual = Basis.from_text(basis_text, 0, names)
        expected = Basis.load(basis_path, 0, names)
        assert ([p.terms for p in actual.polys]) == ([p.terms for p in expected.polys])
        assert (Proof.from_text(proof_text, names).terms) == (Proof.load(proof_path, names).terms)

    def test_exact_input_and_all_proof_identities_are_retained(self):
        run = self.run_solver("all", max_degree=11)
        assert not (run.hit)
        assert (run.seconds) >= (0)
        expected = self.root / "expected.ms"
        self.pres.emit_ms(str(expected), *self.point)
        assert (run.files["input.ms"]) == (expected.read_bytes())
        assert (run.record["status"]) == ("proof-validated")
        assert (run.record["basis_complete"]) is (False)
        assert (run.record["basis_entries"]) == (self.pres.n_relations())
        assert (run.record["request"]["presentation"]["legacy_hash"]) == (self.pres.content_hash())
        assert (run.record["request"]["presentation"]["coefficient_domain"]) == ("Q[q,t]")
        assert (run.record["request"]["point"]) == ([5, 7, 0])
        assert run.record["request"]["schema"] == "zipper-solver-request-v3"
        assert {"target", "target_polynomial"}.isdisjoint(run.record["request"])
        assert (run.record["request"]["monomial_order"]) == ("degree-lex")
        names = self.pres.var_names()
        assert (len(Basis.from_text(run.files["basis.txt"].decode(), 0, names))) == (run.record["basis_entries"])
        assert (len(Proof.from_text(run.files["proof.txt"].decode(), names))) == (run.record["basis_entries"])
        stdout = json.loads(run.files["stdout.txt"])
        assert (stdout) == ({"iterations": 13, "degree": 11})

    def test_forward_references_and_supported_output_spacing(self, subtests):
        for mode in ("forward", "contexts", "comma", "standalone", "blank-basis", "extra"):
            with subtests.test(mode=mode):
                run = self.run_solver(mode)
                assert (set(run.files)) == (set(ARTIFACT_NAMES))
                assert (run.record["basis_entries"]) == (2 if mode == "forward" else 1)

    def test_invalid_and_partial_outputs_never_have_completion_records(self, subtests):
        for mode in ("missing", "empty", "empty-proof", "count", "identity", "cycle",
                     "zero-reference", "cancelled-reference", "basis-reference",
                     "dangling", "double-standalone", "double-sign", "missing-sign",
                     "malformed", "zero-basis", "zero-denominator", "bad-bytes", "mutate-input"):
            output = self.root / mode
            with subtests.test(mode=mode), pytest.raises(ArtifactError):
                self.run_solver(mode, output_root=output)
            assert not ((output / "samples.sqlite3").exists())
            assert not (list(output.glob(".solver-*")))

    def test_nonzero_exit_and_timeout_are_distinct_process_failures(self):
        with pytest.raises(ExecutionError, match="code 7.*test producer"):
            self.run_solver("exit")
        with pytest.raises(ExecutionError, match="exceeded"):
            self.run_solver("timeout", timeout_seconds=0.03)
        assert not ((self.root / "output" / "samples.sqlite3").exists())

    def test_missing_binary_after_configuration_is_a_start_failure(self):
        executable = self.root / "copied-executable"
        shutil.copy2(sys.executable, executable)
        config = self.config(binary=executable)
        executable.unlink()
        # Identity hashing catches this before a process is attempted.
        with pytest.raises(FileNotFoundError):
            run_point(self.pres, "X_2", self.point, config)
        nonexecutable = self.root / "not-an-executable"
        nonexecutable.write_bytes(b"invalid executable contents")
        with pytest.raises(ExecutionError, match="could not start"):
            self.run_solver(binary=nonexecutable)

    def test_output_cache_and_refresh_do_not_conflate_computation_with_completeness(self):
        config = self.config(cache_root=self.root / "cache")
        first = run_point(self.pres, "X_2", self.point, config)
        hit = run_point(self.pres, "S[2][1]", self.point, config)
        assert (hit.hit)
        assert (first.key) == (hit.key)
        assert (self.counter.read_text()) == ("1")
        other = replace(config, output_root=self.root / "other-output")
        copied = run_point(self.pres, "X_2", self.point, other)
        assert (copied.hit)
        assert (copied.key) == (first.key)
        assert (copied.files) == (first.files)
        assert (self.counter.read_text()) == ("1")
        refreshed = run_point(self.pres, "X_2", self.point, config, refresh=True)
        assert not (refreshed.hit)
        assert (self.counter.read_text()) == ("2")
        assert (refreshed.store) == (first.store)
        assert (refreshed.record["generation"]) == (first.record["generation"])
        assert (refreshed.record["basis_complete"]) is (False)

    def test_changed_inputs_binary_source_flags_and_limits_change_keys(self):
        config = self.config()
        first = run_point(self.pres, "X_2", self.point, config)
        points = run_point(self.pres, "X_2", (13, 17, 0), config)
        different = run_point(Presentation(2, untwist="none"), "X_2", self.point, config)
        limit = run_point(self.pres, "X_2", self.point, replace(config, max_iterations=14))
        flag = self.run_solver("comma")
        self.producer.write_text(PRODUCER + "\n# changed executable source\n")
        changed = run_point(self.pres, "X_2", self.point, config)
        assert (len({r.key for r in (first, points, different, limit, flag, changed)})) == (6)


class TestCompletedRunRecords(RunnerFixture):
    def test_corruption_helpers_commit_and_close_their_connections(self):
        run = self.run_solver()
        connect = sqlite3.connect
        connections = []
        def track_connection(*args, **kwargs):
            connection = connect(*args, **kwargs)
            connections.append(connection)
            return connection
        changed = {**run.record, 'basis_entries': 99}
        with patch('tests.test_execution_runner.sqlite3.connect', side_effect=track_connection):
            self.replace_record(run, changed)
            self.alter_file(run, 'proof.txt', b'changed\n')
        assert (len(connections)) == (2)
        for connection in connections:
            with pytest.raises(sqlite3.ProgrammingError):
                connection.execute('SELECT 1')
        with closing(connect(run.store)) as connection:
            record = connection.execute('SELECT record FROM generations').fetchone()[0]
            proof = connection.execute("SELECT content FROM artifacts WHERE name='proof.txt'").fetchone()[0]
        assert (json.loads(record)['basis_entries']) == (99)
        assert (proof) == (b'changed\n')

    def test_absent_database_is_not_a_completed_run(self):
        run = self.run_solver()
        absent = self.root / "absent.sqlite3"
        assert (store.load(absent, run.key, run.record["request"])) is None
        assert not (absent.exists())

    def test_record_schema_request_and_metadata_tampering_are_rejected(self, subtests):
        run = self.run_solver()
        original = run.record
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
            self.replace_record(run, record)
            with subtests.test(record=record), pytest.raises(ArtifactError):
                self.run_solver()
        with closing(sqlite3.connect(run.store)) as connection, connection:
            connection.execute("UPDATE generations SET record=?", (b"{not json",))
        with pytest.raises(ArtifactError):
            self.run_solver()
        self.replace_record(run, original)
        assert (self.run_solver().hit)

    def test_output_byte_tampering_is_rejected_and_refresh_repairs_it(self):
        run = self.run_solver()
        self.alter_file(run, "proof.txt", b"[i1]\n")
        with pytest.raises(ArtifactError, match="recorded identities"):
            self.run_solver()
        fixed = run_point(self.pres, "X_2", self.point, self.config(), refresh=True)
        assert (fixed.files["proof.txt"]) == (b"[i0]\n")
        with closing(sqlite3.connect(run.store)) as connection, connection:
            connection.execute("DELETE FROM artifacts WHERE name='proof.txt'")
        with pytest.raises(ArtifactError):
            self.run_solver()
        fixed = run_point(self.pres, "X_2", self.point, self.config(), refresh=True)
        assert (fixed.files["proof.txt"]) == (b"[i0]\n")

    def _replace_generation(self, run, filename, content):
        """Model a coherently rehashed but mathematically false external record."""
        files = dict(run.files)
        files[filename] = content
        record = json.loads(json.dumps(run.record))
        record["files"] = {name: {"bytes": len(data), "sha256": digest_bytes(data)}
                           for name, data in files.items()}
        generation = digest_bytes(canonical_json(record["files"]))
        record["generation"] = generation
        with closing(sqlite3.connect(run.store)) as connection, connection:
            connection.execute("INSERT INTO generations VALUES(?,?,?)",
                               (run.key, generation, canonical_json(record)))
            connection.executemany("INSERT INTO artifacts VALUES(?,?,?,?)",
                                   [(run.key, generation, name, data)
                                    for name, data in files.items()])
            connection.execute("UPDATE selections SET generation=? WHERE key=?",
                               (generation, run.key))

    def test_a_rehashed_false_proof_still_fails_numeric_validation(self):
        run = self.run_solver()
        self._replace_generation(run, "proof.txt", b"[i1]\n")
        with pytest.raises(ArtifactError, match="fails its input-ideal identity"):
            self.run_solver()

    def test_another_target_cannot_reuse_a_rehashed_false_proof(self):
        presentation = Presentation(3)
        config = self.config()
        first = run_point(presentation, "X_2", self.point, config)
        self._replace_generation(first, "proof.txt", b"[i1]\n")
        with pytest.raises(ArtifactError, match="fails its input-ideal identity"):
            run_point(presentation, "X_3", self.point, config)
        assert self.counter.read_text() == "1"

    def test_a_rehashed_different_input_is_not_the_requested_input(self):
        run = self.run_solver()
        self._replace_generation(run, "input.ms", b"x1,y1\n0\nx1*y1-2\n")
        with pytest.raises(ArtifactError, match="recorded identities"):
            self.run_solver()

    def test_corrupt_shared_cache_is_refused_for_a_new_output(self):
        config = self.config(cache_root=self.root / "cache")
        run = run_point(self.pres, "X_2", self.point, config)
        self.replace_record(run, {}, self.root / "cache" / "samples.sqlite3")
        with pytest.raises(ArtifactError):
            run_point(self.pres, "X_2", self.point,
                      replace(config, output_root=self.root / "another"))

    def test_unselected_generation_is_not_resumed(self):
        run = self.run_solver()
        with closing(sqlite3.connect(run.store)) as connection, connection:
            connection.execute("DELETE FROM selections WHERE key=?", (run.key,))
        rerun = self.run_solver()
        assert not (rerun.hit)
        assert (self.counter.read_text()) == ("2")
        assert (rerun.record["generation"]) == (run.record["generation"])

    def test_multiple_points_use_one_database_without_loose_sample_files(self):
        config = self.config()
        runs = [run_point(self.pres, "X_2", (q, 7, 0), config) for q in (3, 5, 7, 11)]
        assert (len({run.store for run in runs})) == (1)
        assert (set(path.name for path in config.output_root.iterdir())) == ({".gitignore", "samples.sqlite3"})
        with closing(sqlite3.connect(runs[0].store)) as connection, connection:
            assert (connection.execute("SELECT COUNT(*) FROM selections").fetchone()[0]) == (4)
        assert not (list(config.output_root.glob(".solver-*")))

    def test_file_hashing_streams_more_than_one_chunk(self):
        path = self.root / "large-artifact"
        content = b"proof output" * 20000
        path.write_bytes(content)
        assert (file_identity(path)) == ({"bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()})


class TestTargetIndependentSolverReuse(RunnerFixture):
    def test_v3_storage_rejects_target_fields_or_missing_computation_data(self, subtests):
        sample = run_basis_point(self.pres, self.point, self.config())
        request = sample.record['request']
        without_point = {key: value for key, value in request.items() if key != 'point'}
        for malformed in (without_point, {**request, 'target': 'X_2'},
                          {**request, 'target_polynomial': target_identity(self.pres, 'X_2')},
                          {**request, 'unrecognized': 'data'}):
            key = digest_bytes(canonical_json(malformed))
            with subtests.test(fields=set(malformed)), pytest.raises(ArtifactError, match='request-v3'):
                store.publish(sample.store, key, sample.files, malformed, sample.record['basis_entries'])

    def test_different_targets_and_target_free_call_share_one_calculation(self):
        presentation = Presentation(3)
        config = self.config("all", cache_root=self.root / "cache")
        first = run_point(presentation, "X_2", self.point, config)
        another = run_point(presentation, "X_3", self.point, config)
        direct = run_basis_point(presentation, self.point, config)
        assert target_identity(presentation, "X_2") != target_identity(presentation, "X_3")
        assert not first.hit and another.hit and direct.hit
        assert first.key == another.key == direct.key
        assert first.files == another.files == direct.files
        assert first.record == another.record == direct.record
        assert first.record["request"]["schema"] == "zipper-solver-request-v3"
        assert {"target", "target_polynomial"}.isdisjoint(first.record["request"])
        assert not first.record["basis_complete"]
        assert self.counter.read_text() == "1"
        with closing(sqlite3.connect(first.store)) as connection:
            assert connection.execute("SELECT count(*) FROM selections").fetchone() == (1,)

    def test_different_target_in_a_new_output_reuses_the_shared_cache(self):
        presentation = Presentation(3)
        config = self.config(cache_root=self.root / "cache")
        first = run_point(presentation, "X_2", self.point, config)
        other_config = replace(config, output_root=self.root / "another-target")
        other = run_point(presentation, "X_3", self.point, other_config)
        assert not first.hit and other.hit
        assert first.key == other.key
        assert first.record["generation"] == other.record["generation"]
        assert first.store != other.store
        assert store.read(first.store, first.key) == store.read(other.store, other.key)
        assert self.counter.read_text() == "1"

    def test_invalid_targets_fail_before_any_cache_lookup(self, subtests):
        presentation = Presentation(3)
        config = self.config()
        run_point(presentation, "X_2", self.point, config)
        with patch("zipper_algebra.execution.runner.load") as lookup:
            for target in (None, (2, 1), "unknown", "X_4"):
                with subtests.test(target=target), pytest.raises(ValueError):
                    run_point(presentation, target, self.point, config)
        lookup.assert_not_called()
        assert self.counter.read_text() == "1"

    def test_new_requests_preserve_executable_bytes_and_captured_source_identity(self):
        executable = self.root / "solver-launcher"
        executable.write_text(
            f"#!{sys.executable}\nimport os, sys\n"
            "os.execv(sys.executable, [sys.executable, *sys.argv[1:]])\n")
        executable.chmod(0o755)
        config = self.config(binary=executable)
        producer = capture_provenance()
        first = run_basis_point(self.pres, self.point, config, producer=producer)
        # The executable's actual identity must change even though its path
        # and behaviour do not.
        with executable.open("ab") as stream:
            stream.write(b"\n# changed executable bytes\n")
        changed_binary = run_basis_point(self.pres, self.point, config, producer=producer)
        recipe = Path(__file__).resolve().parents[1] / "experiments" / "output.py"
        with_recipe = capture_provenance(extra_sources=(recipe,))
        changed_sources = run_basis_point(self.pres, self.point, config, producer=with_recipe)
        assert len({run.key for run in (first, changed_binary, changed_sources)}) == 3
        assert all(not run.hit for run in (first, changed_binary, changed_sources))
        assert first.record["request"]["solver"]["identity"] != changed_binary.record["request"]["solver"]["identity"]
        assert changed_binary.record["request"]["producer"] != changed_sources.record["request"]["producer"]
        assert self.counter.read_text() == "3"

    def test_old_target_bound_requests_remain_readable_without_becoming_v3_hits(self):
        config = self.config()
        seed = run_basis_point(self.pres, self.point, config)
        legacy_config = replace(config, output_root=self.root / "legacy-output")
        legacy_store = legacy_config.output_root / "samples.sqlite3"
        legacy_request = {**seed.record["request"], "schema": "zipper-solver-request-v2",
                          "target": "X_2", "target_polynomial": target_identity(self.pres, "X_2")}
        legacy_key = digest_bytes(canonical_json(legacy_request))
        legacy_record = store.publish(legacy_store, legacy_key, seed.files, legacy_request,
                                      seed.record["basis_entries"])
        before = store.read(legacy_store, legacy_key)
        fresh = run_point(self.pres, "X_2", self.point, legacy_config)
        assert not fresh.hit
        assert fresh.key != legacy_key and fresh.key == seed.key
        assert self.counter.read_text() == "2"
        assert store.read(legacy_store, legacy_key) == before == (seed.files, legacy_record)
        exported = self.root / "old-generation-export"
        assert store.export(legacy_store, legacy_key, exported) == legacy_record
        assert json.loads((exported / "record.json").read_text()) == legacy_record
        assert (exported / "proof.txt").read_bytes() == seed.files["proof.txt"]
        with closing(sqlite3.connect(legacy_store)) as connection:
            assert connection.execute("SELECT count(*) FROM selections").fetchone() == (2,)
