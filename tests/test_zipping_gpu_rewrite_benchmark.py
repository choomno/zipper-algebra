"""Actual scheduler admissions, exact weights, and complete capture provenance."""

from collections import Counter
import json
import sqlite3
import sys

import pytest

from experiments import zipper_gpu_rewrite_benchmark as bench


@pytest.fixture
def captured(tmp_path):
    path = tmp_path / "rewrites.sqlite3"
    bench.capture(path, [{"n": 3, "word": [1, 2], "name": "branching"},
                         {"n": 3, "word": [-1, -2], "name": "nonoffensive"}],
                  max_steps=1000, batch_size=4)
    return path


def test_capture_covers_actual_catalog_histogram_and_preserves_recursive_rank(captured):
    chunks, meta = bench.load_cases(captured)
    assert meta["rewrite_cases"] == 12 and meta["selected_steps"] == 98
    assert meta["used_cases"] < meta["rewrite_cases"]
    assert {case["n"] for chunk in chunks for case in chunk["cases"]} == {2, 3}
    assert {kind for chunk in chunks for kind in chunk["kinds"]} == {
        "single_catalog_instantiation", "contextual_apply_rule"}
    with sqlite3.connect(captured) as db:
        for root, payload, sha in db.execute("SELECT root,certificate,certificate_sha FROM workloads"):
            cert = bench.decode(payload, sha)
            observed = Counter(case["rule"]["id"] for chunk in chunks if chunk["root"] == root
                               for case in chunk["cases"])
            assert observed == {name: count for name, count in
                                cert["metrics"]["rule_applications"].items()
                                if name.startswith(("single.", "zip."))}
    assert sys.getprofile() is None


def test_reference_restores_huge_actual_weight_without_cross_root_collection(captured):
    chunks, _ = bench.load_cases(captured, {0})
    case = next(case for chunk in chunks for case in chunk["cases"])
    huge = 1 << 512
    case["predecessor"]["coefficient"] = [[-7, huge], [4, -3]]
    other = json.loads(json.dumps(case))
    other["root"] = 1
    output = bench.reference_expand([case, other])
    assert output[0]["root"] == 0 and output[1]["root"] == 1
    assert output[0]["branches"] == output[1]["branches"]
    assert any(abs(c).bit_length() >= 512
               for result in output for branch in result["branches"] for _, c in branch["coefficient"])


def test_single_projection_uses_explicit_exact_branch_permutation():
    first = {"word": [1], "coefficient": [[0, 1]], "control": {"side": "left"}}
    second = {"word": [-1], "coefficient": [[2, -1]], "control": {"side": "right"}}
    assert bench.word_control_permutation([first, second], [second, first]) == [1, 0]
    with pytest.raises(ValueError, match="multiset"):
        bench.word_control_permutation([first], [second])
    with pytest.raises(ValueError, match="unmatched"):
        bench.word_control_permutation([first, second], [first])


def test_legacy_single_instance_hash_is_reinstantiated_before_projection(captured):
    with sqlite3.connect(captured) as db:
        changed = False
        for root, batch, payload, sha in db.execute(
            "SELECT root,batch,payload,payload_sha FROM batches").fetchall():
            saved = bench.decode(payload, sha)
            for case in saved["cases"]:
                if case["capture_kind"] == "single_catalog_instantiation":
                    case.pop("single_projected_to_source_order", None)
                    case["actual_single_instance_sha256"] = "wrong"
                    changed = True
                    break
            if changed:
                payload, sha = bench.encode(saved)
                db.execute("UPDATE batches SET payload=?,payload_sha=? WHERE root=? AND batch=?",
                           (payload, sha, root, batch))
                break
    assert changed
    with pytest.raises(ValueError, match="single source instantiation hash"):
        bench.load_cases(captured)


def test_atomic_refused_batch_is_not_saved(tmp_path):
    path = tmp_path / "prefix.sqlite3"
    manifest = bench.capture(path, [{"n": 3, "word": [1, -2] * 3}],
                             max_steps=10, batch_size=4)
    assert manifest["workloads"][0]["status"] == "incomplete"
    chunks, metadata = bench.load_cases(path)
    with sqlite3.connect(path) as db:
        cert = bench.decode(*db.execute("SELECT certificate,certificate_sha FROM workloads").fetchone())
    assert len(chunks) == cert["metrics"]["batches"]
    assert metadata["selected_steps"] == cert["metrics"]["steps"]


@pytest.mark.parametrize("mutation", ["missing_batch", "orphan_root", "incomplete_marker",
                                     "missing_workload", "corrupt_template"])
def test_capture_loader_rejects_incomplete_or_misbound_data(captured, mutation):
    with sqlite3.connect(captured) as db:
        if mutation == "missing_batch":
            db.execute("DELETE FROM batches WHERE root=0 AND batch=1")
        elif mutation == "orphan_root":
            db.execute("UPDATE batches SET root=99 WHERE root=0")
        elif mutation == "incomplete_marker":
            db.execute("UPDATE metadata SET value='false' WHERE key='complete'")
        elif mutation == "missing_workload":
            db.execute("DELETE FROM workloads WHERE root=0")
        else:
            db.execute("UPDATE templates SET sha='wrong' WHERE rowid=(SELECT min(rowid) FROM templates)")
    with pytest.raises(ValueError):
        bench.load_cases(captured)


def test_dropped_case_with_rehashed_payload_is_rejected_by_rule_counts(captured):
    with sqlite3.connect(captured) as db:
        root, batch, payload, sha = db.execute(
            "SELECT root,batch,payload,payload_sha FROM batches WHERE cases>0 ORDER BY root,batch LIMIT 1"
        ).fetchone()
        saved = bench.decode(payload, sha)
        saved["cases"].pop()
        payload, sha = bench.encode(saved)
        db.execute("UPDATE batches SET payload=?,payload_sha=?,cases=? WHERE root=? AND batch=?",
                   (payload, sha, len(saved["cases"]), root, batch))
    with pytest.raises(ValueError, match="sequence|counts"):
        bench.load_cases(captured)


def test_validly_rehashed_foreign_template_is_rejected_by_catalog_binding(captured):
    with sqlite3.connect(captured) as db:
        old_sha, payload = db.execute("SELECT sha,payload FROM templates LIMIT 1").fetchone()
        rule = bench.decode(payload, old_sha)
        rule["rhs"][0]["coefficient"][0]["integer"] *= -1
        replacement, new_sha = bench.encode(rule)
        db.execute("INSERT INTO templates VALUES(?,?)", (new_sha, replacement))
        for root, batch, payload, sha in db.execute(
            "SELECT root,batch,payload,payload_sha FROM batches").fetchall():
            value = bench.decode(payload, sha)
            for case in value["cases"]:
                if case["rule_sha256"] == old_sha:
                    case["rule_sha256"] = new_sha
                    admitted = {k: case[k] for k in (
                        "root", "id", "batch", "n", "parameters", "words", "offset",
                        "predecessor", "successor_controls", "limits")}
                    admitted["rule"] = rule
                    case["expected"] = bench.reference_expand([admitted])[0]
            payload, sha = bench.encode(value)
            db.execute("UPDATE batches SET payload=?,payload_sha=? WHERE root=? AND batch=?",
                       (payload, sha, root, batch))
    with pytest.raises(ValueError, match="accepted catalog"):
        bench.load_cases(captured)


def test_source_metadata_binding_is_checked(captured):
    with sqlite3.connect(captured) as db:
        root, saved = db.execute("SELECT root,metadata FROM workloads LIMIT 1").fetchone()
        metadata = json.loads(saved)
        metadata["source_sha256"]["catalog.py"] = "wrong"
        db.execute("UPDATE workloads SET metadata=? WHERE root=?", (json.dumps(metadata), root))
    with pytest.raises(ValueError, match="source"):
        bench.load_cases(captured)


def test_strict_atomic_term_limit_does_not_fail_inside_weighted_projection(tmp_path):
    path = tmp_path / "one-term.sqlite3"
    result = bench.capture(path, [{"n": 3, "word": [1, 2]}],
                           max_steps=100, max_terms=1, batch_size=4)
    assert result["workloads"][0]["status"] == "incomplete"
    with sqlite3.connect(path) as db:
        certificate = bench.decode(*db.execute(
            "SELECT certificate,certificate_sha FROM workloads").fetchone())
        assert certificate["metrics"]["batches"] == 0
        assert db.execute("SELECT count(*) FROM batches").fetchone()[0] == 0


def test_unknown_root_and_capture_overwrite_fail_without_mutation(captured):
    original = captured.read_bytes()
    with pytest.raises(ValueError, match="unknown"):
        bench.load_cases(captured, {0, 99})
    with pytest.raises(ValueError, match="new database"):
        bench.capture(captured, [{"n": 2, "word": [1]}])
    assert captured.read_bytes() == original


def test_profiler_is_restored_after_capture_exception(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("capture interruption")
    monkeypatch.setattr(bench.RewriteCaptureSession, "_record", fail)
    path = tmp_path / "failed.sqlite3"
    with pytest.raises(RuntimeError, match="interruption"):
        bench.capture(path, [{"n": 3, "word": [1, 2]}], max_steps=100)
    assert sys.getprofile() is None
    with pytest.raises(ValueError, match="incomplete"):
        bench.load_cases(path)


def test_expansion_benchmark_default_excludes_unused_single_work(captured, monkeypatch):
    from types import SimpleNamespace
    from zipper_algebra.zipping_algorithm import gpu_rewrite

    seen = []

    class FakeExpander:
        metadata = {"backend": "opencl"}
        def expand(self, cases, *, verify):
            assert verify is False
            seen.extend(cases)
            return SimpleNamespace(output=bench.reference_expand(cases),
                                   timings={"packing": 0.0, "h2d": 0.0, "d2h": 0.0},
                                   metadata={"backend": "opencl"})

    monkeypatch.setattr(gpu_rewrite, "expand_cpu", bench.reference_expand)
    cpu = bench.benchmark(captured, backend="cpu", repeats=1)
    gpu = bench.benchmark(captured, backend="opencl", repeats=1, expander=FakeExpander())
    assert gpu["kind"] == cpu["kind"] == "contextual"
    assert gpu["workload"]["selected_rewrite_cases"] < gpu["workload"]["captured_all_rewrite_cases"]
    assert cpu["rows"][0]["output_sha256"] == gpu["rows"][0]["output_sha256"]
    assert gpu["first_use_including_setup_wall"] is None
    assert seen and all(c["rule"]["id"].startswith("zip.") for c in seen)
    for result in [cpu, gpu]:
        timings = result["rows"][0]["timings"]
        assert timings["total"] >= timings["verification"] > 0
        assert timings["captured_output_verification"] > 0


def test_prepared_selection_reuses_integrity_but_repeats_exact_trial_verification(captured, monkeypatch):
    prepared = bench.prepare_benchmark(captured)
    calls = []
    original = bench.reference_expand

    def verify(cases):
        calls.append(len(cases))
        return original(cases)

    def no_reload(*args, **kwargs):
        raise AssertionError("prepared data should not repeat database integrity loading")

    monkeypatch.setattr(bench, "reference_expand", verify)
    monkeypatch.setattr(bench, "load_cases", no_reload)
    result = bench.benchmark(captured, backend="cpu", repeats=2, prepared=prepared)
    assert calls == [len(c) for c in prepared.chunks] * 2
    assert result["preparation_reused"]
    assert result["loading_integrity_and_reference_replay_wall"] == 0
    assert result["preparation_original_wall"] > 0
    assert result["prepared_snapshot_revalidation_wall"] > 0
    assert result["rows"][0]["output_sha256"] == result["rows"][1]["output_sha256"]


def test_prepared_capture_replaced_at_same_path_cannot_receive_old_results(captured, tmp_path):
    prepared = bench.prepare_benchmark(captured)
    replacement = tmp_path / "replacement.sqlite3"
    bench.capture(replacement, [{"n": 3, "word": [1], "name": "different-capture"}],
                  max_steps=1000, batch_size=4)
    replacement.replace(captured)
    with pytest.raises(ValueError, match="destination capture differs"):
        bench.benchmark(captured, backend="cpu", prepared=prepared)
    with sqlite3.connect(captured) as db:
        assert db.execute("SELECT count(*) FROM benchmarks").fetchone()[0] == 0


@pytest.mark.parametrize("table", ["metadata", "workloads", "templates", "batches"])
def test_prepared_capture_mutation_is_rejected_atomically_before_append(captured, table):
    prepared = bench.prepare_benchmark(captured)
    with sqlite3.connect(captured) as db:
        if table == "metadata":
            db.execute("UPDATE metadata SET value='false' WHERE key='complete'")
        elif table == "workloads":
            db.execute("UPDATE workloads SET name='changed' WHERE root=0")
        elif table == "templates":
            db.execute("UPDATE templates SET payload=x'00' WHERE rowid=(SELECT min(rowid) FROM templates)")
        else:
            db.execute("UPDATE batches SET selected_steps=selected_steps+1 WHERE root=0 AND batch=0")
    with pytest.raises(ValueError, match="destination capture differs"):
        bench.benchmark(captured, backend="cpu", prepared=prepared)
    with sqlite3.connect(captured) as db:
        assert db.execute("SELECT count(*) FROM benchmarks").fetchone()[0] == 0


def test_benchmark_appends_do_not_change_capture_identity_and_write_binding_is_locked(captured, monkeypatch):
    prepared = bench.prepare_benchmark(captured)
    trace = []
    original = sqlite3.connect

    def connect(*args, **kwargs):
        connection = original(*args, **kwargs)
        connection.set_trace_callback(trace.append)
        return connection

    monkeypatch.setattr(bench.sqlite3, "connect", connect)
    first = bench.benchmark(captured, backend="cpu", prepared=prepared)
    second = bench.benchmark(captured, backend="cpu", prepared=prepared)
    assert first["workload"]["prepared_capture_tables_sha256"] == second["workload"]["prepared_capture_tables_sha256"]
    assert first["capture_write_identity_validation_wall"] > 0
    assert "BEGIN IMMEDIATE" in trace
    assert next(i for i, statement in enumerate(trace) if statement == "BEGIN IMMEDIATE") < next(
        i for i, statement in enumerate(trace) if statement.startswith("INSERT INTO benchmarks"))
    with original(captured) as db:
        assert db.execute("SELECT count(*) FROM benchmarks").fetchone()[0] == 2
        assert bench.capture_tables_identity(db) == prepared.workload["prepared_capture_tables_sha256"]


def test_write_identity_and_append_hold_a_real_sqlite_writer_fence(captured, monkeypatch):
    prepared = bench.prepare_benchmark(captured)
    original_identity = bench.capture_tables_identity
    competing_writes = []

    def authenticate_under_lock(connection):
        with sqlite3.connect(captured, timeout=0) as other:
            with pytest.raises(sqlite3.OperationalError, match="locked"):
                other.execute("UPDATE metadata SET value='false' WHERE key='complete'")
            competing_writes.append("blocked")
        return original_identity(connection)

    monkeypatch.setattr(bench, "capture_tables_identity", authenticate_under_lock)
    result = bench.benchmark(captured, backend="cpu", prepared=prepared)
    assert competing_writes == ["blocked"]
    assert result["rows"][0]["output_sha256"] == prepared.workload["expected_sha256"]


def test_cpp_constructor_selection_uses_same_external_verifier_without_gpu(captured, monkeypatch):
    from types import SimpleNamespace
    from experiments import zipper_native_rewrite as native

    calls = []

    class FakeNative:
        metadata = {"backend": "cpp"}
        def __init__(self):
            calls.append("constructed")

        def expand(self, cases, *, verify):
            assert verify is False
            calls.append(len(cases))
            return SimpleNamespace(
                output=bench.reference_expand(cases),
                timings={"native_buffer_preparation": 0.0, "word_expansion": 0.0},
                metadata={"backend": "cpp"})

    monkeypatch.setattr(native, "NativeRewriteExpander", FakeNative)
    prepared = bench.prepare_benchmark(captured)
    result = bench.benchmark(captured, backend="cpp", repeats=1, prepared=prepared)
    assert calls == ["constructed", *map(len, prepared.chunks)]
    assert result["backend"] == "cpp" and not result["expander_reused"]
    assert result["rows"][0]["output_sha256"] == prepared.workload["expected_sha256"]
    for chunk in result["rows"][0]["chunks"]:
        assert chunk["metadata"]["verified"]
        assert chunk["timings"]["verification"] > 0
        assert "h2d" not in chunk["timings"] and "d2h" not in chunk["timings"]
    with pytest.raises(ValueError, match="backend differs"):
        bench.benchmark(captured, backend="opencl", prepared=prepared, expander=FakeNative())


@pytest.mark.parametrize("change", ["roots", "mode", "kind", "path", "input", "source"])
def test_prepared_selection_rejects_changed_config_input_or_source(captured, monkeypatch, change):
    prepared = bench.prepare_benchmark(captured)
    arguments = {"database": captured, "backend": "cpu", "prepared": prepared}
    if change == "roots":
        arguments["roots"] = {0}
    elif change == "mode":
        arguments["mode"] = "pooled"
    elif change == "kind":
        arguments["kind"] = "all"
    elif change == "path":
        arguments["database"] = captured.with_name("other.sqlite3")
    elif change == "input":
        prepared.chunks[0][0]["offset"] += 1
    else:
        original = bench.expansion_sources()
        monkeypatch.setattr(bench, "expansion_sources", lambda: {**original, "changed": "source"})
    with pytest.raises(ValueError, match="prepared selection"):
        bench.benchmark(**arguments)
