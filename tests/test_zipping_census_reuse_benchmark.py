"""Reproducible capture and exact acceptance rules for the reuse benchmark."""
from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
from threading import Lock
from types import SimpleNamespace

import pytest

from experiments import zipper_census_reuse_benchmark as bench
from zipper_algebra.zipping_algorithm.interpreter import Limits


def config(**extra):
    return dict(passes=1, case_seconds=1, cache_bytes=0, cache_entries=0,
                limits=asdict(Limits()), allocator=None, labels=None, **extra)


def store(tmp_path):
    path = tmp_path / "benchmark.sqlite3"
    with sqlite3.connect(path) as db:
        bench.schema(db)
        db.execute("INSERT INTO reuse_bench_cases VALUES(1,'case',2,'[1]','{}')")
    return path


def fake_backend(monkeypatch, *, status="complete", wrong=False, timeout=False):
    opened, closed = [], []

    class Backend:
        def __init__(self, name, sidecar, config):
            opened.append(name)
            self.identity = {"producer": {"source_sha256": {"a": "hash"},
                                          "wrapper_source_sha256": "wrapper"},
                             "verifier": {"source_sha256": {"b": "hash"}},
                             "comparison_sources": {"test-interface": "hash"},
                             "harness_sha256": bench.digest(Path(bench.__file__).read_bytes())}

        def evaluate(self, word, n):
            if timeout:
                raise bench.BenchmarkDeadline("bounded test")
            normal = bench.canonical_output([], n=n + int(wrong), complete=True)
            result = {"status": status, "word": word, "n": n}
            return result, normal if status == "complete" else None, {
                "producer_seconds": 0.1, "checker_seconds": 0.2,
                "reconstruction_seconds": 0.3, "store_seconds": 0.0}

        def telemetry(self):
            return {"test": True}

        def close(self):
            closed.append(True)

    monkeypatch.setattr(bench, "Backend", Backend)
    return opened, closed


def test_synthetic_long_cases_have_full_support_and_valid_letters():
    cases = bench.synthetic_cases()
    assert cases == bench.synthetic_cases()
    assert len(cases) == 23
    for case in cases[:18]:
        assert len(case["word"]) in (8, 9, 10)
        assert {abs(x) for x in case["word"]} == set(range(1, case["n"]))
    for case in cases:
        assert all(0 < abs(x) < case["n"] for x in case["word"])


def test_classification_distinguishes_full_support_from_primitive():
    cases = bench.synthetic_cases()[:18]
    classifications = [bench.input_classification(c["word"], c["n"]) for c in cases]
    assert all(c["full_support"] for c in classifications)
    assert sum(c["unchanged_primitive"] for c in classifications) == 16
    reduced = [c for c in classifications if not c["unchanged_primitive"]]
    assert all(c["unique_primitive_requests"] for c in reduced)


def test_capture_preserves_source_and_freezes_cases(tmp_path):
    source = tmp_path / "source.sqlite3"
    with sqlite3.connect(source) as db:
        db.executescript("""CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
          CREATE TABLE words(n,length,ordinal,word,metrics,config_sha,certificate_sha,status);""")
        db.execute("INSERT INTO metadata VALUES('run_state',?)", (json.dumps({"config_sha": "epoch"}),))
        db.execute("INSERT INTO words VALUES(3,8,4,'[1,2]',?, 'epoch','proof','complete')",
                   (json.dumps({"steps": 20}),))
    before = source.read_bytes()
    target = tmp_path / "captured.sqlite3"
    assert bench.capture(source, target) == 24
    assert source.read_bytes() == before
    assert bench.capture(source, target, recent=1) == 24
    with pytest.raises(ValueError, match="another census"):
        bench.capture(tmp_path / "different.sqlite3", target)
    with pytest.raises(ValueError, match="must differ"):
        bench.capture(source, source)
    with sqlite3.connect(target) as db:
        row = db.execute("SELECT label,origin FROM reuse_bench_cases WHERE id=1").fetchone()
        assert row[0] == "census-recent-hard-b3"
        assert json.loads(row[1])["certificate_sha"] == "proof"


def test_worker_commits_payload_and_compares_complete_bytes(tmp_path, monkeypatch):
    path = store(tmp_path)
    opened, closed = fake_backend(monkeypatch)
    for name in ("standalone", "dag", "census"):
        run = bench.create_run(path, name, "retained", config())
        bench.worker(path, run)
    result = bench.compare(path)
    assert result["exact_comparisons"] == 3
    assert sum(p["includes_standalone"] for p in result["backend_pairs"]) == 2
    assert result["paired_complete_inputs"] == 1
    assert result["distinct_inputs_matching_standalone"] == 1
    assert result["distinct_complete_inputs"] == 1
    assert len(opened) == len(closed) == 3
    assert all(run["status"] == "complete" for run in result["runs"])
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT COUNT(*) FROM reuse_bench_payloads").fetchone()[0] == 3
        assert db.execute("SELECT MIN(store_seconds) FROM reuse_bench_measurements").fetchone()[0] > 0


def test_fresh_restarts_between_cases_and_passes(tmp_path, monkeypatch):
    path = store(tmp_path)
    opened, closed = fake_backend(monkeypatch)
    settings = config()
    settings["passes"] = 2
    run = bench.create_run(path, "dag", "fresh", settings)
    bench.worker(path, run)
    assert len(opened) == len(closed) == 2
    assert bench.compare(path)["exact_comparisons"] == 0
    assert bench.compare(path)["same_backend_rechecks"] == 1


@pytest.mark.parametrize("status", ["incomplete", "control_obligation"])
def test_partial_results_are_not_claimed_equal(tmp_path, monkeypatch, status):
    path = store(tmp_path)
    fake_backend(monkeypatch, status=status)
    run = bench.create_run(path, "dag", "retained", config())
    bench.worker(path, run)
    result = bench.compare(path)
    assert result["exact_comparisons"] == result["distinct_complete_inputs"] == 0
    assert result["runs"][0]["measurements"][0]["status"] == status


def test_timeout_records_no_completed_payload_and_closes_session(tmp_path, monkeypatch):
    path = store(tmp_path)
    opened, closed = fake_backend(monkeypatch, timeout=True)
    run = bench.create_run(path, "dag", "retained", config())
    bench.worker(path, run)
    assert len(opened) == len(closed) == 1
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT status FROM reuse_bench_measurements").fetchone() == ("timeout",)
        assert db.execute("SELECT COUNT(*) FROM reuse_bench_payloads").fetchone() == (0,)


def test_disagreeing_complete_outputs_raise(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)
    bench.worker(path, bench.create_run(path, "standalone", "retained", config()))
    fake_backend(monkeypatch, wrong=True)
    bench.worker(path, bench.create_run(path, "dag", "retained", config()))
    with pytest.raises(AssertionError, match="exact completed outputs differ"):
        bench.compare(path)


def test_corrupt_stored_output_fails_integrity(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)
    bench.worker(path, bench.create_run(path, "standalone", "retained", config()))
    with sqlite3.connect(path) as db:
        db.execute("UPDATE reuse_bench_payloads SET output_sha='wrong'")
    with pytest.raises(ValueError, match="hash differs"):
        bench.compare(path)


def test_different_semantic_sources_do_not_form_comparison(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)
    bench.worker(path, bench.create_run(path, "standalone", "retained", config()))
    bench.worker(path, bench.create_run(path, "dag", "retained", config()))
    with sqlite3.connect(path) as db:
        telemetry = json.loads(db.execute("SELECT telemetry FROM reuse_bench_measurements WHERE run=2").fetchone()[0])
        telemetry["identity"]["producer"]["source_sha256"]["a"] = "new"
        db.execute("UPDATE reuse_bench_measurements SET telemetry=? WHERE run=2", (json.dumps(telemetry),))
    assert bench.compare(path)["exact_comparisons"] == 0


def test_measurement_rejects_harness_source_drift(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)
    original = Path.read_bytes
    calls = 0

    def changed(self):
        nonlocal calls
        if self == Path(bench.__file__):
            calls += 1
            return b"old" if calls == 1 else b"new"
        return original(self)

    monkeypatch.setattr(Path, "read_bytes", changed)
    run = bench.create_run(path, "dag", "retained", config())
    with pytest.raises(RuntimeError, match="source changed"):
        bench.worker(path, run)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT status FROM reuse_bench_runs").fetchone() == ("failed",)


def test_backend_bad_mode_rejected(tmp_path):
    with pytest.raises(ValueError, match="unknown"):
        bench.Backend("bad", tmp_path / "sidecar", {})
    with pytest.raises(ValueError, match="invalid backend"):
        bench.create_run(tmp_path / "results", "bad", "fresh", {})


def test_deadline_restores_handler():
    before = bench.signal.getsignal(bench.signal.SIGALRM)
    with bench.deadline(1):
        pass
    assert bench.signal.getsignal(bench.signal.SIGALRM) == before


def test_native_backend_requires_independent_acceptance(monkeypatch, tmp_path):
    result = {"status": "complete", "n": 2, "output": []}

    class Producer:
        metadata = {"source_sha256": {}}
        def __init__(self, **kwargs):
            pass
        def reduce(self, word, **kwargs):
            return deepcopy(result)
        reduce_reused = reduce
        def close(self):
            pass

    class Verifier:
        metadata = {}
        def __init__(self, **kwargs):
            pass
        def check(self, certificate):
            return False
        check_reused = check
        def close(self):
            pass

    monkeypatch.setattr(bench, "NativeSession", Producer)
    monkeypatch.setattr(bench, "CppVerifier", Verifier)
    for name in ("standalone", "dag"):
        backend = bench.Backend(name, tmp_path / name, config())
        with pytest.raises(ValueError, match="independent verifier"):
            backend.evaluate([1], 2)
        backend.close()


def test_verifier_resource_refusal_is_not_completed_or_a_mismatch():
    def refuse(result):
        raise bench.VerifierIncomplete("budget", replay_diagnostic={"kind": "steps"})

    backend = object.__new__(bench.Backend)
    backend.census, backend.name, backend.limits = None, "dag", Limits()
    backend.producer = SimpleNamespace(reduce_reused=lambda *a, **kw: {
        "status": "complete", "output": []})
    backend.verifier = SimpleNamespace(check_reused=refuse)
    result, normal, stats = backend.evaluate([1], 2)
    assert result["status"] == "verification_incomplete"
    assert result["producer_result"]["status"] == "complete"
    assert result["replay_diagnostic"] == {"kind": "steps"}
    assert normal is None and stats["independently_checked"] is False


def test_frozen_inputs_and_unknown_labels_are_rejected(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)
    settings = config()
    settings["labels"] = ["case", "typo"]
    with pytest.raises(ValueError, match="unknown benchmark labels"):
        bench.create_run(path, "dag", "retained", settings)
    run = bench.create_run(path, "dag", "retained", config())
    with sqlite3.connect(path) as db:
        db.execute("UPDATE reuse_bench_cases SET word='[-1]'")
    with pytest.raises(ValueError, match="frozen benchmark input hash"):
        bench.worker(path, run)
    with pytest.raises(ValueError, match="frozen benchmark input hash"):
        bench.compare(path)


def test_existing_sidecar_is_never_silently_reused(tmp_path, monkeypatch):
    path = store(tmp_path)
    opened, _ = fake_backend(monkeypatch)
    run = bench.create_run(path, "census", "retained", config())
    sidecar = path.with_name(path.stem + f"-run{run}-census.sqlite3")
    sidecar.write_bytes(b"prior data")
    with pytest.raises(ValueError, match="new sidecar"):
        bench.worker(path, run)
    assert sidecar.read_bytes() == b"prior data" and opened == []


def test_hard_failure_retains_provisional_case_row(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)

    def fail(self, word, n):
        raise ValueError("independent mismatch")

    monkeypatch.setattr(bench.Backend, "evaluate", fail)
    run = bench.create_run(path, "dag", "retained", config())
    with pytest.raises(ValueError, match="independent mismatch"):
        bench.worker(path, run)
    with sqlite3.connect(path) as db:
        status, error = db.execute("SELECT status,error FROM reuse_bench_measurements").fetchone()
        assert status == "failed" and "independent mismatch" in error
        assert db.execute("SELECT status FROM reuse_bench_runs").fetchone() == ("failed",)


def test_pair_timings_only_use_common_completed_inputs(tmp_path, monkeypatch):
    path = store(tmp_path)
    with sqlite3.connect(path) as db:
        db.execute("INSERT INTO reuse_bench_cases VALUES(2,'second',2,'[-1]','{}')")
    fake_backend(monkeypatch)
    bench.worker(path, bench.create_run(path, "standalone", "retained", config()))
    bench.worker(path, bench.create_run(path, "dag", "retained", config()))
    with sqlite3.connect(path) as db:
        db.execute("UPDATE reuse_bench_measurements SET status='incomplete' WHERE run=2 AND case_id=2")
    summary = bench.compare(path)
    pair = summary["backend_pairs"][0]
    assert pair["common_complete_case_ids"] == [1]
    assert pair["excluded_cases"] == [{"case_id": 2, "left": "complete", "right": "incomplete"}]
    assert summary["paired_complete_inputs"] == 1
    assert summary["unpaired_complete_case_ids"] == [2]


def test_session_identity_is_retained_per_measurement(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)
    original = bench.Backend.__init__
    generation = 0

    def init(self, *args, **kwargs):
        nonlocal generation
        original(self, *args, **kwargs)
        generation += 1
        self.identity["producer"]["source_sha256"]["a"] = str(generation)

    monkeypatch.setattr(bench.Backend, "__init__", init)
    settings = config()
    settings["passes"] = 2
    run = bench.create_run(path, "dag", "fresh", settings)
    bench.worker(path, run)
    with sqlite3.connect(path) as db:
        first = json.loads(db.execute("SELECT identity FROM reuse_bench_runs").fetchone()[0])
        identities = [json.loads(row[0])["identity"] for row in db.execute(
            "SELECT telemetry FROM reuse_bench_measurements ORDER BY pass")]
    assert first == identities[0] and identities[0] != identities[1]
    assert bench.compare(path)["same_backend_rechecks"] == 0


def test_signal_deadline_records_timeout(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)

    def block(self, word, n):
        bench.time.sleep(1)

    monkeypatch.setattr(bench.Backend, "evaluate", block)
    settings = config()
    settings["case_seconds"] = 0.01
    bench.worker(path, bench.create_run(path, "dag", "retained", settings))
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT status FROM reuse_bench_measurements").fetchone() == ("timeout",)


def test_census_complete_requires_verified_output_provenance():
    backend = object.__new__(bench.Backend)
    backend.census = SimpleNamespace(evaluate=lambda *a, **kw: {"status": "complete", "stats": {}})
    with pytest.raises(ValueError, match="independently checked output provenance"):
        backend.evaluate([1], 2)


def test_actual_census_backend_identity_embedding_and_reuse(tmp_path):
    backend = bench.Backend("census", tmp_path / "sidecar.sqlite3", config())
    try:
        result, normal, stats = backend.evaluate([1], 2)
        assert result["status"] == "complete" and normal["completion"] == "complete"
        assert stats["independently_checked"] is True
        assert result["provenance"] == "independently-verified-leaf-with-checked-identity-embedding"
        again, second, repeated = backend.evaluate([1], 2)
        assert second == normal and again["output_digest"] == result["output_digest"]
        assert repeated["leaf_hits"] == 1 and repeated["leaf_misses"] == 0
    finally:
        backend.close()


def test_signal_interrupts_actual_native_transport_pipe_and_reaps_child():
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(5)"],
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    session = object.__new__(bench.NativeSession)
    session._process, session._stderr, session._lock = child, None, Lock()
    try:
        with pytest.raises(bench.BenchmarkDeadline), bench.deadline(0.02):
            session._request({"op": "unanswered"})
        assert session.pid is None and child.poll() is not None
    finally:
        session.close()
        if child.poll() is None:
            child.kill()
            child.wait()


def test_timing_pair_rejects_incompatible_resource_configuration(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)
    bench.worker(path, bench.create_run(path, "standalone", "retained", config()))
    changed = config()
    changed["case_seconds"] = 2
    bench.worker(path, bench.create_run(path, "dag", "retained", changed))
    assert bench.compare(path)["backend_pairs"] == []


def test_legacy_incomplete_identity_is_explicitly_excluded_from_pairing(tmp_path, monkeypatch):
    path = store(tmp_path)
    fake_backend(monkeypatch)
    for name in ("standalone", "dag"):
        bench.worker(path, bench.create_run(path, name, "retained", config()))
    with sqlite3.connect(path) as db:
        for run, payload in db.execute("SELECT run,telemetry FROM reuse_bench_measurements").fetchall():
            value = json.loads(payload)
            value["identity"].pop("comparison_sources")
            db.execute("UPDATE reuse_bench_measurements SET telemetry=? WHERE run=?", (json.dumps(value), run))
    result = bench.compare(path)
    assert result["exact_comparisons"] == 0
    assert result["legacy_complete_outputs_excluded_from_source_bound_pairs"] == 2


def test_fresh_session_refuses_python_drift_since_import(tmp_path, monkeypatch):
    path = store(tmp_path)
    opened, _ = fake_backend(monkeypatch)
    settings = config()
    settings["passes"] = 2
    original = bench.Backend.close
    changed = False

    def close(self):
        nonlocal changed
        original(self)
        changed = True

    monkeypatch.setattr(bench.Backend, "close", close)
    monkeypatch.setattr(bench, "interface_sources", lambda:
                        {"changed": "source"} if changed else bench._IMPORTED_INTERFACE_SOURCES)
    run = bench.create_run(path, "dag", "fresh", settings)
    with pytest.raises(RuntimeError, match="changed after import"):
        bench.worker(path, run)
    assert len(opened) == 1
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT status FROM reuse_bench_runs").fetchone() == ("failed",)
        assert db.execute("SELECT COUNT(*) FROM reuse_bench_measurements").fetchone() == (1,)


@pytest.mark.parametrize("name,expected", [("census", "dag"), ("census-standalone", "standalone")])
def test_census_leaf_backend_uses_matching_independent_checker(monkeypatch, tmp_path, name, expected):
    calls = []

    class Producer:
        metadata = {}

        def __init__(self, **kwargs):
            pass

        def reduce(self, word, **kwargs):
            calls.append(("standalone", word))
            return {"kind": "standalone"}

        def reduce_reused(self, word, **kwargs):
            calls.append(("dag", word))
            return {"kind": "dag"}

        def close(self):
            pass

    class Verifier:
        metadata = {}

        def __init__(self, **kwargs):
            pass

        def check(self, result):
            return result == {"kind": "standalone"}

        def check_reused(self, result):
            return result == {"kind": "dag"}

        def close(self):
            pass

    class Census:
        def __init__(self, path, *, produce, verify, **kwargs):
            self.produce, self.verify = produce, verify

        def close(self):
            pass

    monkeypatch.setattr(bench, "NativeSession", Producer)
    monkeypatch.setattr(bench, "CppVerifier", Verifier)
    monkeypatch.setattr(bench, "ReuseCensus", Census)
    monkeypatch.setattr(bench, "execution_identity", lambda result: result)
    monkeypatch.setattr(bench, "verifier_binding", lambda metadata: metadata)
    backend = bench.Backend(name, tmp_path / "sidecar.sqlite3", config())
    try:
        leaf = backend.census.produce(3, [1, 2])
        assert leaf == {"kind": expected}
        assert backend.census.verify(leaf)
        assert calls == [(expected, []), (expected, [1, 2])]
    finally:
        backend.close()
