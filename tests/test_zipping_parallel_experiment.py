"""The parallel experiment retains the exact checked SQLite input prefix."""
from concurrent.futures import ThreadPoolExecutor
import json
import sqlite3
import threading
import time

import pytest

pytest.importorskip("fcntl", reason="parallel experiment requires Unix flock")

from experiments import zipper_exhaustive_parallel as runner
from zipper_algebra.zipping_algorithm.general_enumeration import exhaust_general
from zipper_algebra.zipping_algorithm.interpreter import Limits


def run(path, **kwargs):
    return runner.exhaust_parallel(path, workers=2, _executor_factory=ThreadPoolExecutor,
                                   wall_seconds=30, **kwargs)


def rows(path):
    with sqlite3.connect(path / "words.sqlite3") as db:
        return db.execute("SELECT n,length,ordinal,word,status,certificate FROM words ORDER BY length,rank_order,ordinal").fetchall()


def test_out_of_order_workers_commit_the_original_literal_order(tmp_path, monkeypatch):
    actual = runner._compute
    later_finished = threading.Event()
    order = []
    def reordered(n, word):
        if not word:
            assert later_finished.wait(10)
        result = actual(n, word)
        order.append(word)
        if word:
            later_finished.set()
        return result
    monkeypatch.setattr(runner, "_compute", reordered)
    result = run(tmp_path, ranks=(3,), max_length=1)
    assert result["stop_reason"] == "max_length"
    assert result["new_or_resumed_words"] == 5
    assert order[0] != []
    assert [json.loads(row[3]) for row in rows(tmp_path)] == [[], [-2], [-1], [1], [2]]
    assert result["peak_in_flight"] <= 4
    assert exhaust_general(tmp_path, ranks=(3,), max_length=1)["audited_words"] == 5


def test_first_incomplete_word_is_kept_and_resumed_without_skipping(tmp_path):
    first = run(tmp_path, ranks=(2,3), max_length=1, limits=Limits(max_steps=4))
    assert first["stop_reason"] == "word_limit"
    assert [(r[0], r[1], r[4]) for r in rows(tmp_path)] == [(2,0,"complete"),(3,0,"incomplete")]
    assert first["next_cursor"] == {"length":0,"rank_index":1,"ordinal":0}
    repeated = run(tmp_path, ranks=(2,3), max_length=1, limits=Limits(max_steps=4))
    assert repeated["stop_reason"] == "word_limit_needs_larger_limits"
    resumed = run(tmp_path, ranks=(2,3), max_length=1, limits=Limits(max_steps=100))
    assert resumed["stop_reason"] == "max_length"
    assert len(rows(tmp_path)) == 8
    assert all(r[4] == "complete" for r in rows(tmp_path))


def test_sequential_prefix_preserved_and_runner_source_bound(tmp_path, monkeypatch):
    exhaust_general(tmp_path, ranks=(3,), max_length=0)
    original = rows(tmp_path)[0]
    result = run(tmp_path, ranks=(3,), max_length=1)
    assert result["audited_words"] == 1
    assert rows(tmp_path)[0] == original
    monkeypatch.setattr(runner, "_runner_hash", lambda: "changed")
    with pytest.raises(ValueError, match="runner source changed"):
        run(tmp_path, ranks=(3,), max_length=1)
    assert len(rows(tmp_path)) == 5


def test_worker_exception_does_not_commit_later_words(tmp_path, monkeypatch):
    actual = runner._compute
    def fail(n, word):
        if word == [-2]:
            time.sleep(.02)
            raise RuntimeError("worker failed")
        return actual(n, word)
    monkeypatch.setattr(runner, "_compute", fail)
    with pytest.raises(RuntimeError, match="worker failed"):
        run(tmp_path, ranks=(3,), max_length=1)
    assert [(r[1],json.loads(r[3])) for r in rows(tmp_path)] == [(0,[])]
    monkeypatch.setattr(runner, "_compute", actual)
    assert run(tmp_path, ranks=(3,), max_length=1)["stop_reason"] == "max_length"
    assert len(rows(tmp_path)) == 5


def test_total_row_cap_and_single_strand_termination(tmp_path):
    capped = run(tmp_path / "cap", ranks=(3,), max_words=3)
    assert capped["stop_reason"] == "max_words"
    assert len(rows(tmp_path / "cap")) == 3
    assert run(tmp_path / "cap", ranks=(3,), max_words=5)["new_or_resumed_words"] == 2
    assert run(tmp_path / "one", ranks=(1,))["stop_reason"] == "all_nonempty_layers_exhausted"


def test_hot_source_change_cannot_commit_a_certificate(tmp_path, monkeypatch):
    actual = runner._assert_source
    def changed(identity, runner_hash):
        raise ValueError("enumeration source or catalog changed during the calculation")
    monkeypatch.setattr(runner, "_assert_source", changed)
    with pytest.raises(ValueError, match="source or catalog changed"):
        run(tmp_path, ranks=(3,), max_length=0)
    assert rows(tmp_path) == []
    monkeypatch.setattr(runner, "_assert_source", actual)
    assert run(tmp_path, ranks=(3,), max_length=0)["stop_reason"] == "max_length"


def test_process_worker_smoke(tmp_path):
    result = runner.exhaust_parallel(tmp_path, ranks=(2,), workers=2, max_length=1, wall_seconds=30)
    assert result["stop_reason"] == "max_length"
    assert len(rows(tmp_path)) == 3


def test_parallel_writer_lock_rejects_overlap_and_releases(tmp_path):
    import fcntl
    with (tmp_path / ".parallel.lock").open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(ValueError, match="another parallel runner"):
            run(tmp_path, ranks=(2,), max_length=0)
        assert not (tmp_path / "words.sqlite3").exists()
    assert run(tmp_path, ranks=(2,), max_length=0)["stop_reason"] == "max_length"
