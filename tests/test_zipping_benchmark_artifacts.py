"""Benchmark archives retain exact bytes and detect conflicting/corrupt data."""

import json
import sqlite3
import zlib

import pytest

from experiments import zipper_benchmark_artifacts as artifacts


def make_run(tmp_path):
    run = tmp_path / "variant-a"
    run.mkdir()
    summary = '{"label":"benchmark", "passes":2}\n'
    (run / "summary.json").write_text(summary)
    files = {
        "pass1-case10.json": b'{"coefficient": [[-2, -7], [3, 9]], "case": 10}\n',
        "pass0-case2.json": b'{\n  "case": 2, "word": [1, -1]\n}\n',
        "pass1-case2.json": b'{"case":2,"trace":false}\n',
    }
    for name, raw in files.items():
        (run / name).write_bytes(raw)
    return run, summary, files


def test_lossless_roundtrip_numeric_alias_order_and_verified_removal(tmp_path):
    run, summary, files = make_run(tmp_path)
    report = artifacts.archive(tmp_path)
    database = tmp_path / "benchmark-results.sqlite3"
    assert report == {"database": str(database), "verified_files": 3, "removed_files": 0}
    assert artifacts.aliases(database, "variant-a") == [
        "variant-a/pass0-case2.json", "variant-a/pass1-case2.json", "variant-a/pass1-case10.json",
    ]
    assert artifacts.aliases(database, "variant-a", pass_index=1) == [
        "variant-a/pass1-case2.json", "variant-a/pass1-case10.json",
    ]
    for name, raw in files.items():
        assert artifacts.load(database, f"variant-a/{name}") == raw
        assert (run / name).read_bytes() == raw
    removed = artifacts.archive(tmp_path, remove=True)
    assert removed["verified_files"] == removed["removed_files"] == 3
    assert not list(run.glob("pass*-case*.json"))
    assert (run / "summary.json").read_text() == summary
    with sqlite3.connect(database) as db:
        assert db.execute("SELECT count(*) FROM certificates").fetchone()[0] == 3
        assert db.execute("SELECT summary FROM runs WHERE alias='variant-a'").fetchone()[0] == summary
    for name, raw in files.items():
        assert artifacts.load(database, f"variant-a/{name}") == raw


@pytest.mark.parametrize("mutation", [
    "UPDATE certificates SET sha256='wrong'",
    "UPDATE certificates SET raw_bytes=raw_bytes+1",
])
def test_hash_or_length_corruption_is_rejected(tmp_path, mutation):
    make_run(tmp_path)
    artifacts.archive(tmp_path)
    database = tmp_path / "benchmark-results.sqlite3"
    with sqlite3.connect(database) as db:
        db.execute(mutation)
    with pytest.raises(ValueError, match="hash or size"):
        artifacts.load(database, "variant-a/pass0-case2.json")


def test_changed_compressed_payload_is_rejected(tmp_path):
    make_run(tmp_path)
    artifacts.archive(tmp_path)
    database = tmp_path / "benchmark-results.sqlite3"
    with sqlite3.connect(database) as db:
        db.execute("UPDATE certificates SET payload=?", (zlib.compress(b'{"changed":true}'),))
    with pytest.raises(ValueError, match="hash or size"):
        artifacts.load(database, "variant-a/pass0-case2.json")


def test_conflicting_alias_rolls_back_and_preserves_sources(tmp_path):
    run, _, files = make_run(tmp_path)
    artifacts.archive(tmp_path)
    changed = run / "pass0-case2.json"
    changed.write_text(json.dumps({"case": 2, "different": True}))
    with pytest.raises(ValueError, match="conflicting artifact alias"):
        artifacts.archive(tmp_path, remove=True)
    assert all((run / name).is_file() for name in files)
    assert artifacts.load(tmp_path / "benchmark-results.sqlite3", "variant-a/pass0-case2.json") == files["pass0-case2.json"]


def test_unknown_alias_and_traversal_are_rejected(tmp_path):
    make_run(tmp_path)
    artifacts.archive(tmp_path)
    database = tmp_path / "benchmark-results.sqlite3"
    with pytest.raises(KeyError):
        artifacts.load(database, "variant-a/missing.json")
    assert artifacts.aliases(database, "missing-run") == []
    for alias in ("../variant-a/pass0-case2.json", "/absolute/certificate.json"):
        with pytest.raises(ValueError, match="relative path"):
            artifacts.load(database, alias)
