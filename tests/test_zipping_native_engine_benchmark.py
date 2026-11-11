"""Stored timings remain bound to verified certificates and frozen sources."""
import hashlib
import json
import sqlite3
import zlib

import pytest

from experiments import zipper_native_engine_benchmark as bench
from zipper_algebra.zipping_algorithm.interpreter import Limits


def test_python_smoke_stores_exact_replay_bound_compressed_certificate(tmp_path, monkeypatch):
    sources = bench.source_identity()
    monkeypatch.setattr(bench, "source_identity", lambda: sources)
    path = tmp_path / "timings.sqlite3"
    campaign = bench.run(path, cases=["smoke"], rounds=1,
                         backends=["python_collected", "python_memoized"],
                         limits=Limits(max_steps=2000), batch_size=64)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT status FROM campaigns WHERE id=?", (campaign,)).fetchone() == ("complete",)
        rows = db.execute("""SELECT backend,certificate_zlib,certificate_sha256,mathematical_sha256,
                           producer_seconds,verification_seconds,verified_total_seconds,verified
                           FROM measurements ORDER BY position""").fetchall()
        assert len(rows) == 2
        mathematical_hashes = []
        for backend, compressed, raw_hash, math_hash, producer, verification, total, verified in rows:
            blob = zlib.decompress(compressed)
            assert hashlib.sha256(blob).hexdigest() == raw_hash
            certificate = json.loads(blob)
            assert bench.digest(bench.reference_form(certificate, backend)) == math_hash
            assert producer > 0 and verification > 0 and total == producer + verification
            assert verified == 1 and certificate["status"] == "complete"
            mathematical_hashes.append(math_hash)
        assert len(set(mathematical_hashes)) == 1


@pytest.mark.parametrize("cases,backends,rounds", [([], ["native"], 1), (["smoke"], [], 1),
    (["smoke"], ["native", "native"], 1), (["smoke"], ["native"], 0),
    (["missing"], ["native"], 1), (["smoke"], ["missing"], 1)])
def test_invalid_benchmark_configuration_fails_before_creating_database(tmp_path, cases, backends, rounds):
    path = tmp_path / "timings.sqlite3"
    with pytest.raises(ValueError):
        bench.run(path, cases=cases, rounds=rounds, backends=backends, limits=Limits(), batch_size=64)
    assert not path.exists()


def test_mathematical_certificate_normalization_does_not_hide_wrong_metrics():
    session = bench.memoized.MemoizedSession()
    cert = session.reduce((), n=1)
    expected = bench.collected.CollectedSession().reduce((), n=1)
    assert bench.reference_form(cert, "python_memoized") == expected
    cert["metrics"]["steps"] += 1
    assert bench.reference_form(cert, "python_memoized") != expected


def test_source_drift_refuses_mixed_epoch_measurement(tmp_path, monkeypatch):
    identities = iter([{"epoch": 1}, {"epoch": 2}])
    monkeypatch.setattr(bench, "source_identity", lambda: next(identities))
    path = tmp_path / "timings.sqlite3"
    with pytest.raises(RuntimeError, match="source changed"):
        bench.run(path, cases=["smoke"], rounds=1, backends=["python_collected"],
                  limits=Limits(), batch_size=64)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT COUNT(*) FROM measurements").fetchone() == (0,)
        assert db.execute("SELECT status FROM campaigns").fetchone() == ("failed",)


def test_unverified_result_never_enters_sqlite(tmp_path, monkeypatch):
    normalize = bench.reference_form
    def corrupt(certificate, backend):
        result = normalize(certificate, backend)
        result["metrics"]["steps"] += 1
        return result
    monkeypatch.setattr(bench, "reference_form", corrupt)
    path = tmp_path / "timings.sqlite3"
    with pytest.raises(AssertionError, match="exact replay mismatch"):
        bench.run(path, cases=["smoke"], rounds=1, backends=["python_collected"],
                  limits=Limits(), batch_size=64)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT COUNT(*) FROM measurements").fetchone() == (0,)
        assert db.execute("SELECT status FROM campaigns").fetchone() == ("failed",)


@pytest.mark.parametrize("alias", [True, 1.0])
def test_equal_valued_boolean_or_float_coefficient_is_rejected(tmp_path, monkeypatch, alias):
    monkeypatch.setitem(bench.CASES, "unit", (1, []))
    normalize = bench.reference_form
    def corrupt_type(certificate, backend):
        result = normalize(certificate, backend)
        assert result["output"][0]["coefficient"] == [[0, 1]]
        result["output"][0]["coefficient"][0][1] = alias
        return result
    monkeypatch.setattr(bench, "reference_form", corrupt_type)
    path = tmp_path / "timings.sqlite3"
    with pytest.raises(AssertionError, match="exact replay mismatch"):
        bench.run(path, cases=["unit"], rounds=1, backends=["python_collected"],
                  limits=Limits(), batch_size=64)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT COUNT(*) FROM measurements").fetchone() == (0,)
