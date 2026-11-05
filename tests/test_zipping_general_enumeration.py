"""General scheduler enumeration uses actual words and preserves stopped prefixes."""
from __future__ import annotations

import itertools
import json
import sqlite3

import pytest

from zipper_algebra.zipping_algorithm.general_enumeration import exhaust_general
from zipper_algebra.zipping_algorithm.interpreter import Limits


def test_actual_general_rank_three_word_layers(tmp_path):
    summary = exhaust_general(tmp_path, ranks=(3,), max_length=2)
    assert summary["stop_reason"] == "max_length"
    assert all(row["layer_complete"] for row in summary["layers"])
    assert summary["new_or_resumed_words"] == 21
    with sqlite3.connect(tmp_path / "words.sqlite3") as db:
        inputs = [json.loads(row[0]) for row in db.execute("SELECT word FROM words WHERE length=2 ORDER BY ordinal")]
        assert inputs == [list(w) for w in itertools.product((-2,-1,1,2), repeat=2)]
        assert all(json.loads(row[0])["schema"] == "zipper-algorithm-run-v1" for row in db.execute("SELECT certificate FROM words"))
    resumed = exhaust_general(tmp_path, ranks=(3,), max_length=2)
    assert resumed["audited_words"] == 21
    assert resumed["new_or_resumed_words"] == 0


def test_bounded_word_is_not_skipped_and_resumes_with_larger_limits(tmp_path):
    first = exhaust_general(tmp_path, ranks=(3,), max_length=0, limits=Limits(max_steps=1))
    assert first["stop_reason"] == "word_limit"
    assert first["layers"][0]["layer_complete"] is False
    assert first["next_cursor"] == {"length":0,"rank_index":0,"ordinal":0}
    second = exhaust_general(tmp_path, ranks=(3,), max_length=0, limits=Limits(max_steps=1))
    assert second["stop_reason"] == "word_limit_needs_larger_limits"
    third = exhaust_general(tmp_path, ranks=(3,), max_length=0, limits=Limits(max_steps=30))
    assert third["stop_reason"] == "max_length"
    assert third["layers"][0]["layer_complete"]


def test_general_driver_wall_cutoff_and_rank_one_termination(tmp_path, monkeypatch):
    import zipper_algebra.zipping_algorithm.general_enumeration as module
    one = exhaust_general(tmp_path / "one", ranks=(1,))
    assert one["stop_reason"] == "all_nonempty_layers_exhausted"
    clock = itertools.chain([0,0,0], itertools.repeat(1))
    monkeypatch.setattr(module.time, "perf_counter", lambda: next(clock))
    stopped = exhaust_general(tmp_path / "wall", ranks=(3,), wall_seconds=0.5)
    assert stopped["stop_reason"] == "wall_budget_between_words"
    assert stopped["new_or_resumed_words"] == 0


def test_general_identity_and_compact_metrics_tampering_rejected(tmp_path):
    exhaust_general(tmp_path, ranks=(3,), max_length=0)
    with pytest.raises(ValueError, match="configuration changed"):
        exhaust_general(tmp_path, ranks=(3,), max_length=0, trace=True)
    with sqlite3.connect(tmp_path / "words.sqlite3") as db:
        db.execute("UPDATE words SET metrics='{}'")
    with pytest.raises(ValueError, match="status"):
        exhaust_general(tmp_path, ranks=(3,), max_length=0)


@pytest.mark.parametrize("kwargs", [{"ranks":()}, {"ranks":(3,3)}, {"ranks":(0,)},
    {"wall_seconds":0}, {"wall_seconds":float("nan")}, {"max_length":-1},
    {"limits":{}}, {"trace":1}])
def test_general_invalid_configuration_rejected_before_writes(tmp_path, kwargs):
    with pytest.raises(ValueError):
        exhaust_general(tmp_path, **kwargs)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("mutation,match", [
    ("UPDATE words SET ordinal=2 WHERE length=0", "consecutive"),
    ("UPDATE words SET certificate_sha='bad'", "hash mismatch"),
    ("DELETE FROM configurations", "configuration is missing"),
    ("UPDATE configurations SET config='{}'", "configuration hash"),
    ("UPDATE metadata SET value='[0,0,0]' WHERE key='cursor'", "cursor"),
])
def test_general_audit_rejects_tampered_prefix_and_configuration(tmp_path, mutation, match):
    exhaust_general(tmp_path, ranks=(3,), max_length=0)
    with sqlite3.connect(tmp_path / "words.sqlite3") as db:
        db.execute(mutation)
    with pytest.raises(ValueError, match=match):
        exhaust_general(tmp_path, ranks=(3,), max_length=0)


def test_general_failed_checker_cannot_commit_a_word(tmp_path, monkeypatch):
    import zipper_algebra.zipping_algorithm.general_enumeration as module
    monkeypatch.setattr(module.AlgorithmSession, "check", lambda self, value: False)
    with pytest.raises(ValueError, match="failed exact checking"):
        exhaust_general(tmp_path, ranks=(3,), max_length=0)
    with sqlite3.connect(tmp_path / "words.sqlite3") as db:
        assert db.execute("SELECT count(*) FROM words").fetchone()[0] == 0
