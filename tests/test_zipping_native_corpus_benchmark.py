"""Literal corpus completeness and independently verified SQLite storage."""
import json
import sqlite3

import pytest

from experiments import zipper_native_corpus_benchmark as corpus


def test_default_corpus_is_every_literal_word_with_no_reduction():
    rows = corpus.literal_corpus()
    assert len(rows) == 600
    for n, cutoff, expected in ((3, 4, 341), (4, 3, 259)):
        selected = [row for row in rows if row["n"] == n]
        assert len(selected) == expected
        assert len({tuple(row["word"]) for row in selected}) == expected
        for length in range(cutoff + 1):
            assert sum(row["length"] == length for row in selected) == (2 * (n - 1)) ** length
        assert any(row["word"] == [1, -1] for row in selected)
        assert any(row["word"] == [-1, 1] for row in selected)
        assert any(row["word"] == [] for row in selected)


@pytest.mark.parametrize("bounds", [(), ((1, 2),), ((3, -1),), ((True, 2),), ((3, True),), ((3, 1), (3, 2))])
def test_invalid_corpus_boundaries_fail(bounds):
    with pytest.raises(ValueError):
        corpus.literal_corpus(bounds)


def test_retained_python_smoke_checks_every_input_and_rotates_producer_order(tmp_path, monkeypatch):
    path = tmp_path / "corpus.sqlite3"
    manifest = corpus.literal_corpus(((3, 1),))
    identity = corpus.corpus_source_identity(manifest)
    monkeypatch.setattr(corpus, "corpus_source_identity", lambda rows: identity)
    campaign = corpus.run(path, bounds=((3, 1),), backends=("python_memoized", "python_collected"))
    with sqlite3.connect(path) as db:
        config, status = db.execute("SELECT configuration_json,status FROM campaigns WHERE id=?", (campaign,)).fetchone()
        assert status == "complete" and json.loads(config)["type"] == "retained_literal_corpus"
        rows = db.execute("SELECT case_name,backend,position,mathematical_sha256,verified,telemetry_json FROM measurements ORDER BY id").fetchall()
        assert len(rows) == 10
        for index, item in enumerate(manifest):
            pair = rows[2 * index:2 * index + 2]
            assert all(row[0] == item["case_name"] and row[4] == 1 for row in pair)
            assert pair[0][2] == 0 and pair[1][2] == 1 and pair[0][3] == pair[1][3]
            assert pair[0][1] == ("python_memoized" if index % 2 == 0 else "python_collected")
        infos = [json.loads(row[5])["cache_info"] for row in rows if row[1] == "python_memoized"]
        assert infos[-1]["hits"] > 0
        assert all(b["evaluations"] >= a["evaluations"] for a, b in zip(infos, infos[1:]))


def test_corpus_source_drift_refuses_measurement(tmp_path, monkeypatch):
    identities = iter([{"epoch": 1}, {"epoch": 2}])
    monkeypatch.setattr(corpus, "corpus_source_identity", lambda rows: next(identities))
    path = tmp_path / "corpus.sqlite3"
    with pytest.raises(RuntimeError, match="source changed"):
        corpus.run(path, bounds=((3, 0),), backends=("python_collected",))
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT COUNT(*) FROM measurements").fetchone() == (0,)
        assert db.execute("SELECT status FROM campaigns").fetchone() == ("failed",)


@pytest.mark.parametrize("alias", [True, 1.0])
def test_corpus_rejects_equal_valued_boolean_or_float_coefficient(tmp_path, monkeypatch, alias):
    normalize = corpus.shared.reference_form
    def corrupt_type(certificate, backend):
        result = normalize(certificate, backend)
        assert result["output"][0]["coefficient"] == [[0, 1]]
        result["output"][0]["coefficient"][0][1] = alias
        return result
    monkeypatch.setattr(corpus.shared, "reference_form", corrupt_type)
    path = tmp_path / "corpus.sqlite3"
    with pytest.raises(AssertionError, match="exact replay mismatch"):
        corpus.run(path, bounds=((3, 0),), backends=("python_collected",))
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT COUNT(*) FROM measurements").fetchone() == (0,)
