"""Cross-epoch comparison must retain exact evidence and resource refusals."""
import hashlib
import json
from pathlib import Path
import sqlite3
import zlib

import pytest

from experiments import zipper_census_reuse_benchmark as reuse
from experiments import zipper_cold_epoch_benchmark as cold


def database(tmp_path: Path) -> Path:
    path = tmp_path / "epochs.sqlite3"
    raw = cold.encode({"completion": "complete", "terms": []}).encode()
    telemetry = {"stats": {"independently_checked": True}, "input_classification": {"primitive": True}}
    with sqlite3.connect(path) as db:
        reuse.schema(db)
        db.execute("INSERT INTO reuse_bench_cases VALUES(1,'case',2,'[1]','{}')")
        reuse.frozen_cases(db, establish=True)
        db.execute("""CREATE TABLE cold_epoch_runs(run INTEGER,epoch_name TEXT,epoch_path TEXT,
          cohort TEXT,round_index INTEGER,source_manifest TEXT,driver_sha256 TEXT,
          conditions_before TEXT,conditions_after TEXT,process_seconds REAL,process_cpu_seconds REAL)""")
        for run in (1, 2):
            db.execute("INSERT INTO reuse_bench_runs VALUES(?, 'census','fresh','','complete','{}','{}',0,NULL)", (run,))
            db.execute("INSERT INTO cold_epoch_runs VALUES(?,?,?,'cold',0,?,'driver','{}','{}',1,1)",
                       (run, f"epoch{run}", f"/epoch{run}", cold.encode({"source": str(run), "experiments/zipper_census_reuse_benchmark.py": "worker"})))
            db.execute("INSERT INTO reuse_bench_measurements VALUES(?,0,1,'complete',1,0.2,0.5,0,0.1,100,?,NULL)",
                       (run, cold.encode(telemetry)))
            db.execute("INSERT INTO reuse_bench_payloads VALUES(?,0,1,NULL,?,?)",
                       (run, zlib.compress(raw), hashlib.sha256(raw).hexdigest()))
    return path


def test_pair_retains_changed_sources_and_exact_outputs(tmp_path):
    result = cold.pair_report(database(tmp_path), 1, 2)
    assert len(result["matches"]) == 1
    assert result["epochs"][1]["source_manifest"] != result["epochs"][2]["source_manifest"]
    assert result["left_over_right"] == 1


def test_pair_excludes_resource_limited_output(tmp_path):
    path = database(tmp_path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE reuse_bench_measurements SET status='incomplete' WHERE run=2")
    result = cold.pair_report(path, 1, 2)
    assert result["matches"] == []
    assert result["left_seconds"] == 0
    assert result["excluded"][0]["right"] == "incomplete"


@pytest.mark.parametrize("change", ["payload", "checked", "classification", "hash", "config", "driver", "run_failed"])
def test_pair_rejects_invalid_comparisons(tmp_path, change):
    path = database(tmp_path)
    with sqlite3.connect(path) as db:
        if change in ("payload", "hash"):
            raw = cold.encode({"completion": "complete", "terms": [1]}).encode()
            sha = hashlib.sha256(raw).hexdigest() if change == "payload" else "wrong"
            db.execute("UPDATE reuse_bench_payloads SET output=?,output_sha=? WHERE run=2", (zlib.compress(raw), sha))
        elif change in ("checked", "classification"):
            telemetry = json.loads(db.execute("SELECT telemetry FROM reuse_bench_measurements WHERE run=2").fetchone()[0])
            if change == "checked":
                telemetry["stats"]["independently_checked"] = False
            else:
                telemetry["input_classification"] = {}
            db.execute("UPDATE reuse_bench_measurements SET telemetry=? WHERE run=2", (cold.encode(telemetry),))
        elif change == "config":
            db.execute("UPDATE reuse_bench_runs SET config=? WHERE id=2", (cold.encode({"cache_bytes": 1}),))
        elif change == "driver":
            db.execute("UPDATE cold_epoch_runs SET driver_sha256='changed' WHERE run=2")
        else:
            db.execute("UPDATE reuse_bench_runs SET status='failed' WHERE id=2")
    with pytest.raises((ValueError, AssertionError)):
        cold.pair_report(path, 1, 2)


def test_manifest_includes_native_python_and_catalog(tmp_path):
    for filename in ("src/z/example.py", "src/z/catalog.json", "native/engine.cpp", "experiments/zipper_census_reuse_benchmark.py"):
        path = tmp_path / filename
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(filename)
    assert len(cold.source_manifest(tmp_path)) == 4


def test_module_entrypoint_rejects_preloaded_foreign_experiments(tmp_path):
    import os
    import subprocess
    import sys

    cases, config = tmp_path / "cases.json", tmp_path / "config.json"
    cases.write_text("[]")
    config.write_text("{}")
    root = Path(__file__).resolve().parents[1]
    environment = dict(os.environ, PYTHONPATH=str(root))
    result = subprocess.run([sys.executable, "-m", "experiments.zipper_cold_epoch_benchmark",
        str(tmp_path / "result.sqlite3"), "--epoch", str(tmp_path / "foreign"),
        "--epoch-name", "foreign", "--cases", str(cases), "--config", str(config)],
        env=environment, cwd=root, capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert "foreign experiments package" in result.stderr
    assert not (tmp_path / "result.sqlite3").exists()


def test_driver_import_snapshot_rejects_drift(monkeypatch):
    monkeypatch.setattr(cold, "_IMPORTED_DRIVER_SHA256", "changed")
    with pytest.raises(RuntimeError, match="changed after import"):
        cold.assert_driver_source()


def test_pair_suppresses_optimization_ratio_across_other_axes(tmp_path):
    path = database(tmp_path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE reuse_bench_runs SET backend='standalone' WHERE id=2")
    result = cold.pair_report(path, 1, 2)
    assert len(result["matches"]) == 1
    assert result["left_over_right"] is None
    assert not result["optimization_timing_comparable"]


def test_pair_rejects_unknown_measurement_case(tmp_path):
    path = database(tmp_path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE reuse_bench_measurements SET case_id=99 WHERE run=2")
    with pytest.raises(ValueError, match="unknown frozen case"):
        cold.pair_report(path, 1, 2)


@pytest.mark.parametrize("mode,passes", [("retained", 1), ("fresh", 2)])
def test_cold_requires_fresh_single_pass(mode, passes):
    with pytest.raises(ValueError, match="fresh sessions and one pass"):
        cold.validate_cold_cases([], {"passes": passes}, mode, lambda word, n: {})


def test_cold_rejects_duplicate_or_reduced_primitive_requests():
    case = {"label": "input", "word": [1], "n": 2}
    classify = lambda word, n: {"unchanged_primitive": True, "unique_primitive_requests": [{"key": "same"}]}
    with pytest.raises(ValueError, match="repeats a primitive"):
        cold.validate_cold_cases([case, case], {"passes": 1}, "fresh", classify)
    with pytest.raises(ValueError, match="unchanged primitive"):
        cold.validate_cold_cases([case], {"passes": 1}, "fresh", lambda word, n: {
            "unchanged_primitive": False, "unique_primitive_requests": []})
