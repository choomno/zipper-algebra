"""Real zipper captures bind collection inputs to replay-checked scheduler output."""

import json
import sqlite3

import pytest

from experiments import zipper_gpu_benchmark as bench


def test_reference_keeps_input_roots_states_and_exact_large_coefficients():
    huge = 1 << 240
    rows = [(1, 2, -3, huge), (1, 2, -3, 7), (1, 2, -3, -huge),
            (2, 2, -3, -7), (1, 3, -3, 2), (1, 3, -3, -2), (1, 2, 4, -9)]
    assert bench.reference_collection(rows) == ((1, 2, -3, 7), (1, 2, 4, -9), (2, 2, -3, -7))


def test_real_capture_replays_and_binds_every_batch_to_saved_scheduler_output(tmp_path):
    path = tmp_path / "frontiers.sqlite3"
    manifest = bench.capture(path, [
        {"name": "actual-positive", "n": 3, "word": [1, 2]},
        {"name": "actual-negative", "n": 3, "word": [-1, -2]},
    ], max_steps=1000, batch_size=4)
    assert all(row["status"] == "complete" for row in manifest["workloads"])
    records, metadata = bench.load_frontiers(path)
    assert records and {row[0] for row in records} == {0, 1}
    assert metadata["roots"] == 2
    assert bench.digest(bench.reference_collection(records)) == metadata["expected_output_sha256"]
    filtered, summary = bench.load_frontiers(path, {1})
    assert {row[0] for row in filtered} == {1} and summary["roots"] == 1
    assert len(filtered) < len(records)
    with sqlite3.connect(path) as db:
        rows = db.execute("SELECT metadata FROM workloads").fetchall()
    assert all(json.loads(row[0])["certificate_reloaded_and_replayed"] for row in rows)
    original = path.read_bytes()
    with pytest.raises(ValueError, match="new database"):
        bench.capture(path, [{"n": 2, "word": [1]}])
    assert path.read_bytes() == original


def test_capture_discards_uncommitted_limit_batch_and_rejects_tampering(tmp_path):
    path = tmp_path / "frontiers.sqlite3"
    manifest = bench.capture(path, [{"n": 3, "word": [1, -2] * 3}],
                             max_steps=10, batch_size=4)
    assert manifest["workloads"][0]["status"] == "incomplete"
    with sqlite3.connect(path) as db:
        certificate = db.execute("SELECT certificate,certificate_sha FROM workloads").fetchone()
        cert = bench.decode(*certificate)
        assert db.execute("SELECT count(*) FROM frontiers").fetchone()[0] == cert["metrics"]["batches"]
        db.execute("UPDATE frontiers SET payload_sha='bad' WHERE batch=0")
    with pytest.raises(ValueError, match="hash mismatch"):
        bench.load_frontiers(path)


def test_cpu_and_gpu_benchmark_use_same_external_oracle_and_full_timing(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from zipper_algebra.zipping_algorithm import gpu_collection

    path = tmp_path / "frontiers.sqlite3"
    bench.capture(path, [{"n": 3, "word": [1, 2]}], max_steps=1000, batch_size=4)
    calls = []

    class FakeCollector:
        def __init__(self, *, device):
            assert device == 0

        def collect(self, records, *, verify):
            calls.append(verify)
            return SimpleNamespace(output=gpu_collection.collect_cpu(records),
                                   timings={"packing": 0.0, "h2d": 0.0, "reduce": 0.0,
                                            "d2h": 0.0, "total": 0.0},
                                   metadata={"verified": verify, "device": "test double"})

    monkeypatch.setattr(gpu_collection, "OpenClCollector", FakeCollector)
    cpu = bench.benchmark(path, backend="cpu", repeats=1)
    gpu = bench.benchmark(path, backend="opencl", repeats=1)
    assert calls and not any(calls)
    assert cpu["rows"][0]["output_sha256"] == gpu["rows"][0]["output_sha256"]
    assert cpu["mode"] == gpu["mode"] == "sequential"
    for result in (cpu, gpu):
        row = result["rows"][0]
        assert row["timings"]["total"] == row["observed_call_wall"]
        assert row["timings"]["total"] >= row["timings"]["verification"] > 0
        assert row["timings"]["captured_output_verification"] > 0
    assert all(chunk["metadata"]["verified"] for chunk in gpu["rows"][0]["chunks"])
    pooled = bench.benchmark(path, backend="cpu", repeats=1, mode="pooled")
    assert pooled["chunk_count"] == 1
    assert "offline" in pooled["timing_scope"]


@pytest.mark.parametrize("mutation", ["drop_root_batches", "drop_middle_batch", "orphan_root"])
def test_loader_rejects_missing_batches_and_orphan_roots(tmp_path, mutation):
    path = tmp_path / "frontiers.sqlite3"
    bench.capture(path, [{"n": 3, "word": [1, 2]}], max_steps=1000, batch_size=4)
    with sqlite3.connect(path) as db:
        if mutation == "drop_root_batches":
            db.execute("DELETE FROM frontiers WHERE root=0")
        elif mutation == "drop_middle_batch":
            db.execute("DELETE FROM frontiers WHERE root=0 AND batch=1")
        else:
            db.execute("UPDATE frontiers SET root=99")
    with pytest.raises(ValueError, match="batch|orphan"):
        bench.load_frontiers(path)


def test_loader_requires_every_requested_root(tmp_path):
    path = tmp_path / "frontiers.sqlite3"
    bench.capture(path, [{"n": 3, "word": [1, 2]}], max_steps=1000, batch_size=4)
    with pytest.raises(ValueError, match="requested workload root"):
        bench.load_frontiers(path, {0, 99})


def test_reused_collector_backend_and_initialization_labels(tmp_path):
    from types import SimpleNamespace
    from zipper_algebra.zipping_algorithm import gpu_collection

    path = tmp_path / "frontiers.sqlite3"
    bench.capture(path, [{"n": 3, "word": [1, 2]}], max_steps=1000, batch_size=4)
    wrong = SimpleNamespace(metadata={"backend": "cuda"})
    with pytest.raises(ValueError, match="backend differs"):
        bench.benchmark(path, backend="opencl", repeats=1, collector=wrong)

    class ReusedCollector:
        metadata = {"backend": "opencl"}
        initialization_seconds = 0.125

        def collect(self, records, *, verify):
            assert verify is False
            return SimpleNamespace(output=gpu_collection.collect_cpu(records),
                                   timings={"total": 0.0}, metadata=self.metadata)

    result = bench.benchmark(path, backend="opencl", repeats=1, collector=ReusedCollector())
    assert result["collector_reused"] and result["backend_setup_wall"] == 0
    assert result["first_use_with_setup_wall"] is None
    assert result["retained_collector_initialization_seconds"] == 0.125


def test_shadow_reconciliation_includes_rolled_back_attempt_and_detects_wrong_final_output():
    from copy import deepcopy
    from zipper_algebra.zipping_algorithm.interpreter import Limits

    session = bench.CaptureSession(0)
    cert = session.reduce([1, 2], n=3, limits=Limits(max_terms=1), batch_size=4)
    assert cert["status"] == "incomplete"
    session.reconcile(cert)
    assert len(session.batches) == cert["metrics"]["batches"] + 1
    session = bench.CaptureSession(0)
    cert = session.reduce([1, 2], n=3, batch_size=4)
    session.reconcile(cert)
    bad = deepcopy(cert)
    bad["output"][0]["coefficient"][0][1] += 1
    with pytest.raises(ValueError, match="final shadow output"):
        session.reconcile(bad)


def test_incomplete_capture_and_deleted_workload_are_rejected(tmp_path):
    path = tmp_path / "frontiers.sqlite3"
    bench.capture(path, [{"n": 3, "word": [1, 2]}], max_steps=1000, batch_size=4)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE metadata SET value='false' WHERE key='capture_complete'")
    with pytest.raises(ValueError, match="not marked complete"):
        bench.load_frontiers(path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE metadata SET value='true' WHERE key='capture_complete'")
        db.execute("DELETE FROM workloads")
        db.execute("DELETE FROM frontiers")
    with pytest.raises(ValueError, match="intended capture manifest"):
        bench.load_frontiers(path)
    with pytest.raises(ValueError, match="positive step budget"):
        bench.capture(tmp_path / "zero.sqlite3", [{"n": 3, "word": [1]}], max_steps=0)
