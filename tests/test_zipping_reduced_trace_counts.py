"""Cross-check candidate counts by independently reducing every small raw word."""
from itertools import product
from pathlib import Path
import sqlite3

import pytest

import experiments.zipper_reduced_trace_counts as counting
from experiments.zipper_reduced_trace_counts import count_layers, export_counts, write_counts


def reduced_trace_normal(word):
    remaining = list(word)
    while True:
        pair = next(((i, j) for i in range(len(remaining))
                     for j in range(i + 1, len(remaining))
                     if remaining[j] == -remaining[i]
                     and all(abs(abs(a) - abs(remaining[i])) >= 2
                             for a in remaining[i + 1:j])), None)
        if pair is None:
            break
        i, j = pair
        del remaining[j]
        del remaining[i]
    canonical = []
    while remaining:
        available = [(a, i) for i, a in enumerate(remaining)
                     if all(abs(abs(a) - abs(b)) >= 2 for b in remaining[:i])]
        a, i = min(available)
        canonical.append(a)
        del remaining[i]
    return tuple(canonical)


def has_under(word, n):
    labels = list(range(1, n + 1))
    under = False
    for a in reversed(word):
        i = abs(a) - 1
        left, right = labels[i:i + 2]
        over = left if a > 0 else right
        under |= n in (left, right) and over != n
        labels[i], labels[i + 1] = right, left
    return under


@pytest.mark.parametrize("n,cutoff", [(2, 6), (3, 5), (4, 4), (5, 3)])
def test_counts_against_independent_reduction(n, cutoff):
    letters = tuple(range(1 - n, 0)) + tuple(range(1, n))
    by_length = [set() for _ in range(cutoff + 1)]
    for length in range(cutoff + 1):
        for word in product(letters, repeat=length):
            normal = reduced_trace_normal(word[::-1])[::-1]
            if {abs(a) for a in normal} == set(range(1, n)):
                by_length[len(normal)].add(normal)
    expected = [(len(words), sum(has_under(w, n) for w in words)) for words in by_length]
    assert count_layers(n, cutoff) == expected


def test_cancellation_can_cross_commuting_letters():
    assert reduced_trace_normal((1, 3, -1)) == (3,)
    assert reduced_trace_normal((1, 2, -1)) == (1, 2, -1)
    assert count_layers(2, 3) == [(0, 0), (2, 1), (2, 2), (2, 2)]


def test_sqlite_provenance_idempotence_and_conflict(tmp_path: Path):
    database = tmp_path / "candidate.sqlite3"
    run_id = write_counts(database, maximum_rank=3, maximum_length=3)
    assert write_counts(database, maximum_rank=3, maximum_length=3) == run_id
    with sqlite3.connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM counts").fetchone() == (8,)
        assert connection.execute("SELECT reduced_trace_full_support FROM counts WHERE n=3 AND length=3").fetchone() == ("32",)
        connection.execute("UPDATE counts SET reduced_trace_full_support='999' WHERE n=3 AND length=3")
    with pytest.raises(ValueError, match="disagrees"):
        write_counts(database, maximum_rank=3, maximum_length=3)


@pytest.mark.parametrize("n,length", [(1, 2), (True, 2), (3, -1), (3, False)])
def test_reject_invalid_domain(n, length):
    with pytest.raises(ValueError):
        count_layers(n, length)


def test_source_drift_rolls_back(tmp_path, monkeypatch):
    snapshots = iter([{"script_sha256": "before"}, {"script_sha256": "after"}])
    monkeypatch.setattr(counting, "_source_identity", lambda: next(snapshots))
    database = tmp_path / "drift.sqlite3"
    with pytest.raises(ValueError, match="source changed"):
        write_counts(database, maximum_rank=2, maximum_length=1)
    with sqlite3.connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM runs").fetchone() == (0,)


def test_export_selects_one_run(tmp_path):
    import json
    database = tmp_path / "counts.sqlite3"
    small = write_counts(database, maximum_rank=2, maximum_length=1)
    write_counts(database, maximum_rank=2, maximum_length=2)
    output = tmp_path / "counts.json"
    export_counts(database, small, output)
    result = json.loads(output.read_text())
    assert len(result["rows"]) == 2
    assert result["totals_excluding_identity"]["reduced_trace_full_support"] == 2
    assert result["totals_excluding_identity"]["current_nonpeelable"] == 1


@pytest.mark.parametrize("table", ["counts", "current_comparison"])
@pytest.mark.parametrize("damage", ["missing", "extra"])
def test_export_rejects_incomplete_or_extra_grid(tmp_path, table, damage):
    database = tmp_path / "counts.sqlite3"
    run_id = write_counts(database, maximum_rank=2, maximum_length=1)
    with sqlite3.connect(database) as connection:
        if damage == "missing":
            connection.execute(f"DELETE FROM {table} WHERE run_id=? AND length=1", (run_id,))
        else:
            connection.execute(f"UPDATE {table} SET length=2 WHERE run_id=? AND length=1", (run_id,))
    output = tmp_path / "counts.json"
    output.write_text("preserved")
    with pytest.raises(ValueError, match="complete declared grid"):
        export_counts(database, run_id, output)
    assert output.read_text() == "preserved"
