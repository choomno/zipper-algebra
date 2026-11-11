"""Scaling samples are reproducible, disjoint, and bound to saved outputs."""

import hashlib
import json
import sqlite3
import zlib

import pytest

pytest.importorskip("resource", reason="worker scaling uses Unix resource metrics")
from experiments import zipper_worker_scaling as scaling


def make_database(tmp_path):
    database = tmp_path / "words.sqlite3"
    words = [[-1, -1], [-1, 1], [1, -1], [1, 1]]
    digests = {}
    with sqlite3.connect(database) as db:
        db.execute(
            "CREATE TABLE words(n INTEGER,length INTEGER,ordinal INTEGER,word TEXT,"
            "status TEXT,certificate BLOB,certificate_sha TEXT)"
        )
        for ordinal, word in enumerate(words):
            certificate = {
                "input_word": word,
                "output": [{"word": [1], "coefficient": [[0, 1]]}],
                "metrics": {"steps": ordinal + 3},
            }
            raw = json.dumps(certificate, indent=2).encode()
            digests[ordinal] = hashlib.sha256(raw).hexdigest()
            db.execute("INSERT INTO words VALUES(2,2,?,?,?,?,?)", (
                ordinal, json.dumps(word), "complete", zlib.compress(raw), digests[ordinal],
            ))
    return database, words, digests


def test_seeded_selection_is_disjoint_and_preserves_output_binding(tmp_path):
    database, words, digests = make_database(tmp_path)
    original = database.read_bytes()
    options = {"strata": ((2, 2),), "measured_per_stratum": 2,
               "warmup_per_stratum": 1, "seed": 17}
    first = scaling.prepare(database, tmp_path / "first", **options)
    second = scaling.prepare(database, tmp_path / "second", **options)
    assert first == second
    assert database.read_bytes() == original
    assert len(first["warmup"]) == 1 and len(first["measured"]) == 2
    warmup = {case["ordinal"] for case in first["warmup"]}
    measured = {case["ordinal"] for case in first["measured"]}
    assert warmup.isdisjoint(measured) and len(warmup | measured) == 3
    expected_output_digest = hashlib.sha256(
        b'["{\\"coefficient\\":[[0,1]],\\"word\\":[1]}"]'
    ).hexdigest()
    for case in first["warmup"] + first["measured"]:
        ordinal = case["ordinal"]
        assert case["word"] == words[ordinal]
        assert case["steps"] == ordinal + 3
        assert case["source_certificate_sha"] == digests[ordinal]
        assert case["output_sha"] == expected_output_digest
    saved = json.loads((tmp_path / "first" / "worker-scaling-inputs.json").read_text())
    assert saved == first


def test_incomplete_layer_is_rejected_before_creating_a_manifest(tmp_path):
    database, _, _ = make_database(tmp_path)
    with sqlite3.connect(database) as db:
        db.execute("UPDATE words SET status='incomplete' WHERE ordinal=3")
    output = tmp_path / "sample"
    with pytest.raises(ValueError):
        scaling.prepare(database, output, strata=((2, 2),),
                        measured_per_stratum=2, warmup_per_stratum=1)
    assert not (output / "worker-scaling-inputs.json").exists()


def test_corrupt_source_certificate_hash_is_rejected(tmp_path):
    database, _, _ = make_database(tmp_path)
    with sqlite3.connect(database) as db:
        db.execute("UPDATE words SET certificate_sha='wrong'")
    output = tmp_path / "sample"
    with pytest.raises(ValueError):
        scaling.prepare(database, output, strata=((2, 2),),
                        measured_per_stratum=2, warmup_per_stratum=1)
    assert not (output / "worker-scaling-inputs.json").exists()
