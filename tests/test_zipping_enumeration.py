"""Literal-word prefix, resumability and persisted evidence regressions."""
from __future__ import annotations

import itertools
import json
import sqlite3

import pytest

from zipper_algebra.zipping_algorithm.enumeration import alphabet, exhaust, iter_words, word_count


def test_literal_words_include_every_inverse_pair():
    assert alphabet(3) == (-2, -1, 1, 2)
    words = list(iter_words(3, 2))
    assert words == list(itertools.product((-2, -1, 1, 2), repeat=2))
    assert (1, -1) in words and (-1, 1) in words
    assert len(words) == word_count(3, 2) == 16
    assert list(iter_words(1, 0)) == [()]
    assert list(iter_words(1, 2)) == []
    assert word_count(1, 0) == 1 and word_count(1, 2) == 0
    assert next(iter_words(10, 100)) == (-9,)*100


def test_complete_layers_and_reloaded_certificates(tmp_path):
    summary = exhaust(tmp_path, ranks=(1, 2, 3), max_length=2)
    assert summary["stop_reason"] == "max_length"
    assert summary["unsupported_ranks"] == [3]
    assert summary["new_or_resumed_words"] == 8
    assert all(layer["layer_complete"] for layer in summary["layers"])
    with sqlite3.connect(tmp_path / "words.sqlite3") as db:
        rows = db.execute("SELECT word,certificate,ordinal FROM words WHERE n=2 AND length=2 ORDER BY ordinal").fetchall()
        assert [json.loads(row[0]) for row in rows] == [[-1,-1],[-1,1],[1,-1],[1,1]]
        assert all(json.loads(row[1])["trace_enabled"] is False for row in rows)
    resumed = exhaust(tmp_path, ranks=(1, 2, 3), max_length=3)
    assert resumed["audited_words"] == 8
    assert resumed["new_or_resumed_words"] == 8
    assert all(layer["layer_complete"] for layer in resumed["layers"])
    assert json.loads((tmp_path / "summary.json").read_text())["stop_reason"] == "max_length"


def test_limit_hit_word_never_skipped_and_larger_limits_resume(tmp_path):
    first = exhaust(tmp_path, ranks=(2,), max_steps=0, max_length=1)
    assert first["stop_reason"] == "word_limit"
    assert first["next_cursor"] == {"length":0,"rank_index":0,"ordinal":0}
    assert first["layers"][0]["layer_complete"] is False
    unchanged = exhaust(tmp_path, ranks=(2,), max_steps=0, max_length=1)
    assert unchanged["stop_reason"] == "word_limit_needs_larger_limits"
    assert unchanged["new_or_resumed_words"] == 0
    completed = exhaust(tmp_path, ranks=(2,), max_steps=30, max_length=1)
    assert completed["new_or_resumed_words"] == 3
    assert all(layer["layer_complete"] for layer in completed["layers"])


def test_soft_wall_budget_preserves_exact_next_ordinal(tmp_path, monkeypatch):
    import zipper_algebra.zipping_algorithm.enumeration as module
    clock = itertools.chain([0, 0, 0, 0], itertools.repeat(1))
    monkeypatch.setattr(module.time, "perf_counter", lambda: next(clock))
    summary = exhaust(tmp_path, ranks=(2,), wall_seconds=0.5)
    assert summary["stop_reason"] == "wall_budget_between_words"
    assert summary["new_or_resumed_words"] == 1
    assert summary["next_cursor"] == {"length":1,"rank_index":0,"ordinal":0}


def test_rank_one_and_unsupported_rank_stop_explicitly(tmp_path):
    one = exhaust(tmp_path / "one", ranks=(1,))
    assert one["stop_reason"] == "all_nonempty_layers_exhausted"
    assert one["layers"] == [{"n":1,"length":0,"expected_words":1,"visited_words":1,"complete_words":1,"layer_complete":True}]
    unsupported = exhaust(tmp_path / "higher", ranks=(3,4))
    assert unsupported["stop_reason"] == "unsupported_ranks"
    assert unsupported["new_or_resumed_words"] == 0
    assert unsupported["layers"] == []


def test_resume_rejects_source_configuration_and_database_tampering(tmp_path):
    exhaust(tmp_path, ranks=(2,), max_length=1)
    with pytest.raises(ValueError, match="configuration changed"):
        exhaust(tmp_path, ranks=(2,), trace=True, max_length=1)
    db = sqlite3.connect(tmp_path / "words.sqlite3")
    db.execute("UPDATE words SET certificate_sha='bad' WHERE length=1 AND ordinal=0")
    db.commit()
    db.close()
    with pytest.raises(ValueError, match="hash mismatch"):
        exhaust(tmp_path, ranks=(2,), max_length=1)


@pytest.mark.parametrize("kwargs", [{"ranks":()}, {"ranks":(2,2)}, {"ranks":(0,)}, {"ranks":(True,)},
    {"wall_seconds":0}, {"wall_seconds":float("inf")}, {"max_length":-1}, {"max_steps":-1},
    {"max_terms":0}, {"trace":1}])
def test_invalid_configuration_rejected_before_writing(tmp_path, kwargs):
    with pytest.raises(ValueError):
        exhaust(tmp_path, **kwargs)
    assert list(tmp_path.iterdir()) == []


def test_invalid_word_enumeration_parameters():
    for n, length in [(0,1),(2,-1),(True,0),(2,False)]:
        with pytest.raises(ValueError):
            list(iter_words(n,length))


@pytest.mark.parametrize("mutation,match", [
    ("DELETE FROM words WHERE length=0", "consecutive"),
    ("DELETE FROM configurations", "configuration is missing"),
    ("UPDATE configurations SET config='{}'", "configuration hash"),
    ("UPDATE words SET metrics='{}' WHERE length=0", "status"),
    ("UPDATE words SET diagnostic='bad' WHERE length=0", "status"),
    ("UPDATE metadata SET value='[0,0,0]' WHERE key='cursor'", "cursor"),
])
def test_resume_rejects_missing_prefix_and_tampered_compact_records(tmp_path, mutation, match):
    exhaust(tmp_path, ranks=(2,), max_length=1)
    with sqlite3.connect(tmp_path / "words.sqlite3") as db:
        db.execute(mutation)
    with pytest.raises(ValueError, match=match):
        exhaust(tmp_path, ranks=(2,), max_length=1)


def test_rejected_producer_is_never_committed(tmp_path, monkeypatch):
    import zipper_algebra.zipping_algorithm.enumeration as module
    monkeypatch.setattr(module.TwoStrandSession, "check", lambda self, value: {"pass":False})
    with pytest.raises(ValueError, match="failed exact checking"):
        exhaust(tmp_path, ranks=(2,), max_length=0)
    with sqlite3.connect(tmp_path / "words.sqlite3") as db:
        assert db.execute("SELECT count(*) FROM words").fetchone()[0] == 0
        assert db.execute("SELECT value FROM metadata WHERE key='cursor'").fetchone()[0] == "[0,0,0]"


def test_exhaustive_rank_two_steps_have_independent_closed_formula(tmp_path):
    exhaust(tmp_path, ranks=(2,), max_length=3)
    with sqlite3.connect(tmp_path / "words.sqlite3") as db:
        for length in range(4):
            metrics = [json.loads(row[0]) for row in db.execute("SELECT metrics FROM words WHERE length=?", (length,))]
            assert sum(row["steps"] for row in metrics) == 4*3**length-3*2**length
