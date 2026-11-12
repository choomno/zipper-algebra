"""The async comparison requires complete, exact, input-bound evidence."""
import hashlib
import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace
import queue
import zlib

import pytest

from experiments import zipper_async_benchmark as bench
from zipper_algebra.zipping_algorithm.census_reuse import _canonical_digest, _json
from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import completed_output
from zipper_algebra.zipping_algorithm.native_verifier import VerifierError, VerifierIncomplete


def database(path: Path, *, word=None, state="certified") -> Path:
    word = [-1, -2] if word is None else word
    certificate = CollectedSession().reduce(word, n=3)
    raw = _json(certificate).encode()
    with sqlite3.connect(path) as db:
        db.executescript("""
          CREATE TABLE async_jobs(policy TEXT,key TEXT,n INTEGER,word TEXT,state TEXT);
          CREATE TABLE async_candidates(policy TEXT,key TEXT,payload BLOB,digest TEXT,output_digest TEXT,bytes INTEGER);
          CREATE TABLE async_attempts(id INTEGER PRIMARY KEY,policy TEXT,key TEXT,stage TEXT,status TEXT,record TEXT);
        """)
        db.execute("INSERT INTO async_jobs VALUES('p','k',3,?,?)", (_json(word), state))
        blob = zlib.compress(raw)
        db.execute("INSERT INTO async_candidates VALUES('p','k',?,?,?,?)",
                   (blob, hashlib.sha256(raw).hexdigest(), _canonical_digest(completed_output(certificate)), len(blob)))
        db.execute("INSERT INTO async_attempts VALUES(NULL,'p','k','compute','candidate',?)", ('{"seconds":2}',))
        if state == "certified":
            db.execute("INSERT INTO async_attempts VALUES(NULL,'p','k','certify','certified',?)", (json.dumps({
                "digest": hashlib.sha256(raw).hexdigest(), "seconds": 3,
                "verifier": {"abi": "test", "source_sha256": {"test": "test"}, "binary_sha256": "test"}}),))
    return path


def test_exact_outputs_and_store_accounting(tmp_path):
    a, b = database(tmp_path / "a.sqlite3"), database(tmp_path / "b.sqlite3")
    report = bench.compare(a, b)
    assert report["all_inputs_certified_and_equal"]
    assert report["exact_output_matches"] == [[3, [-1, -2]]]
    summary = bench.summarize_store(a)
    assert summary["counts"] == {"certified": 1}
    assert summary["stages"] == {"compute": {"attempts": 1, "seconds": 2},
                                  "certify": {"attempts": 1, "seconds": 3}}
    assert summary["database_bytes"] >= summary["candidate_bytes"] > 0


def test_uncertified_candidates_are_not_compared(tmp_path):
    a = database(tmp_path / "a.sqlite3")
    b = database(tmp_path / "b.sqlite3", state="verification_incomplete")
    report = bench.compare(a, b)
    assert not report["all_inputs_certified_and_equal"]
    assert report["exact_output_matches"] == []
    assert report["unresolved"][0]["right"] == "verification_incomplete"


@pytest.mark.parametrize("field,value,message", [
    ("digest", "corrupt", "attempt digest"),
    ("output_digest", "corrupt", "output digest"),
])
def test_corrupt_digests_rejected(tmp_path, field, value, message):
    path = database(tmp_path / "a.sqlite3")
    with sqlite3.connect(path) as db:
        db.execute(f"UPDATE async_candidates SET {field}=?", (value,))
    with pytest.raises(ValueError, match=message):
        bench.read_outputs(path)


def test_exact_comparison_does_not_trust_matching_hashes(tmp_path, monkeypatch):
    a, b = database(tmp_path / "a.sqlite3"), database(tmp_path / "b.sqlite3")
    original = bench.read_outputs
    def altered(path):
        result = original(path)
        if path == b:
            key = next(iter(result))
            result[key][1]["untrusted_extra"] = True
        return result
    monkeypatch.setattr(bench, "read_outputs", altered)
    with pytest.raises(ValueError, match="full outputs differ"):
        bench.compare(a, b)


def test_inputs_must_match_and_payload_must_bind(tmp_path):
    a, b = database(tmp_path / "a.sqlite3"), database(tmp_path / "b.sqlite3", word=[-2, -1])
    with pytest.raises(ValueError, match="input sets differ"):
        bench.compare(a, b)
    with sqlite3.connect(b) as db:
        db.execute("UPDATE async_jobs SET word='[1,2]'")
    with pytest.raises(ValueError, match="input binding"):
        bench.read_outputs(b)


def test_duplicate_policy_inputs_rejected(tmp_path):
    path = database(tmp_path / "a.sqlite3", state="candidate")
    with sqlite3.connect(path) as db:
        db.execute("INSERT INTO async_jobs SELECT 'other',key,n,word,state FROM async_jobs")
    with pytest.raises(ValueError, match="multiple policies"):
        bench.read_outputs(path)


@pytest.mark.parametrize("change,match", [("missing", "payload is missing"), ("attempt", "no certification"),
                                         ("verifier", "authenticated verifier")])
def test_certified_status_requires_bound_evidence(tmp_path, change, match):
    path = database(tmp_path / "a.sqlite3")
    with sqlite3.connect(path) as db:
        if change == "missing":
            db.execute("DELETE FROM async_candidates")
        elif change == "attempt":
            db.execute("DELETE FROM async_attempts WHERE stage='certify'")
        else:
            record = json.loads(db.execute("SELECT record FROM async_attempts WHERE stage='certify'").fetchone()[0])
            record["verifier"] = {}
            db.execute("UPDATE async_attempts SET record=? WHERE stage='certify'", (json.dumps(record),))
    with pytest.raises(ValueError, match=match):
        bench.read_outputs(path)


def test_receive_detects_worker_death_and_timeout():
    channel = queue.Queue()
    channel.put("ready")
    assert bench.receive(channel, [], 1) == "ready"
    with pytest.raises(RuntimeError, match="exited unsuccessfully"):
        bench.receive(channel, [SimpleNamespace(exitcode=1)], .01)
    with pytest.raises(queue.Empty):
        bench.receive(channel, [SimpleNamespace(exitcode=0)], .001)


@pytest.mark.parametrize("mode", ["reuse", "standalone"])
def test_shared_native_verification_adapter(mode):
    verifier = SimpleNamespace(check=lambda result: result, check_reused=lambda result: result)
    assert bench.verify_result(verifier, True, mode) is True


@pytest.mark.parametrize("report,error,expected", [
    ({"ok": False, "replay_status": "complete"}, VerifierError("mismatch"), False),
    ({"ok": False, "error": "native certificate differs from independent C++ replay"}, VerifierError("mismatch"), False),
    ({"ok": False, "replay_status": "incomplete"}, VerifierIncomplete("limit"), VerifierIncomplete),
    ({}, VerifierError("process died"), VerifierError),
])
def test_shared_adapter_keeps_mismatch_distinct_from_operational_failure(report, error, expected):
    def fail(result):
        raise error
    verifier = SimpleNamespace(check_reused=fail, last_result=report)
    if expected is False:
        assert bench.verify_result(verifier, {}, "reuse") is False
    else:
        with pytest.raises(expected):
            bench.verify_result(verifier, {}, "reuse")


@pytest.mark.parametrize("requests", [[], None, [[2, [1]]], [[3, [True]]], [[3, [3]]], [[3, [-1]], [3, [-1]]]])
def test_request_validation(requests):
    with pytest.raises(ValueError):
        bench.validate_requests(requests)


def test_distinct_requests_and_manifest(tmp_path):
    assert bench.validate_requests([[3, [1]], [4, [1]]]) == [[3, [1]], [4, [1]]]
    for filename in ("src/a.py", "native/a.cpp", "native/vendor/a.hpp",
                     "experiments/zipper_async_benchmark.py", "experiments/zipper_async_census.py",
                     "experiments/zipper_census_reuse_benchmark.py"):
        p = tmp_path / filename
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(filename))
    before = bench.manifest(tmp_path)
    (tmp_path / "native/a.cpp").write_text("changed")
    assert len(before) == 6 and before != bench.manifest(tmp_path)
