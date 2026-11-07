"""Independent exact actions and brute-force checks of the revised policy."""
from itertools import product
import json
from pathlib import Path
import sqlite3

import pytest

from experiments.zipper_artin_baseline_counts import (
    count_layers, decreasing_rules, decreasing_triples, export_counts, representative, write_counts,
)
from experiments.zipper_reduced_trace_counts import count_layers as reduced_counts
from zipper_algebra.zipping_algorithm.extraction import _artin_action


def action(word, n):
    return _artin_action(tuple(word), n, max_image_letters=100000,
                         max_operations=1000000, operations=[0])


def has_under(word, n):
    labels = list(range(1, n + 1))
    result = False
    for a in reversed(word):
        i = abs(a) - 1
        left, right = labels[i:i + 2]
        result |= n in (left, right) and (left if a > 0 else right) != n
        labels[i], labels[i + 1] = right, left
    return result


def test_signed_artin_relations_exactly_six():
    rules = decreasing_triples(5)
    assert len(rules) == 18
    for left, right in rules.items():
        assert left > right
        assert action(left, 5) == action(right, 5)
    for signs in product((-1, 1), repeat=3):
        e, f, g = signs
        left, right = (e, 2 * f, g), (2 * g, f, 2 * e)
        assert (action(left, 3) == action(right, 3)) == (not (e == g and e != f))


@pytest.mark.parametrize("length_four", [False, True])
@pytest.mark.parametrize("n,cutoff", [(2, 5), (3, 5), (4, 4), (5, 3)])
def test_count_fixed_words_by_brute_force(n, cutoff, length_four):
    letters = tuple(range(1 - n, 0)) + tuple(range(1, n))
    expected = []
    for length in range(cutoff + 1):
        fixed = []
        for word in product(letters, repeat=length):
            normal = representative(word, n=n, length_four=length_four)
            assert representative(normal, n=n, length_four=length_four) == normal
            assert len(normal) <= len(word)
            if length <= 3:
                assert action(normal, n) == action(word, n)
            if normal == word and {abs(a) for a in word} == set(range(1, n)):
                fixed.append(word)
        expected.append((len(fixed), sum(has_under(w, n) for w in fixed)))
    assert count_layers(n, cutoff, length_four=length_four) == expected
    assert all(a <= c and b <= d for (a, b), (c, d) in zip(expected, reduced_counts(n, cutoff)))


def test_sqlite_bound_idempotent_conflict(tmp_path: Path):
    path = tmp_path / 'artin.sqlite3'
    run = write_counts(path, maximum_rank=3, maximum_length=3)
    assert write_counts(path, maximum_rank=3, maximum_length=3) == run
    with sqlite3.connect(path) as db:
        assert db.execute('SELECT COUNT(*) FROM counts').fetchone() == (8,)
        db.execute("UPDATE counts SET normal_words='999' WHERE n=3 AND length=3")
    with pytest.raises(ValueError, match='disagrees'):
        write_counts(path, maximum_rank=3, maximum_length=3)


@pytest.mark.parametrize('n,length', [(1, 2), (True, 2), (3, -1), (3, False)])
def test_invalid(n, length):
    with pytest.raises(ValueError):
        count_layers(n, length)


def test_export_requires_complete_grid(tmp_path):
    database = tmp_path / 'artin.sqlite3'
    output = tmp_path / 'counts.json'
    run = write_counts(database, maximum_rank=3, maximum_length=2)
    export_counts(database, run, output)
    assert output.is_file()
    with sqlite3.connect(database) as c:
        c.execute('DELETE FROM counts WHERE n=3 AND length=2')
    with pytest.raises(ValueError, match='incomplete'):
        export_counts(database, run, output)
    with pytest.raises(ValueError, match='unknown'):
        export_counts(database, 'missing', output)


def test_length_four_rules_are_exact():
    rules = decreasing_rules(5, length_four=True)
    assert any(len(left) == 4 for left in rules)
    for left, right in rules.items():
        assert (len(right), right) < (len(left), left)
        assert action(left, 5) == action(right, 5)
    assert all(a <= c and b <= d for (a, b), (c, d) in zip(
        count_layers(4, 8, length_four=True), count_layers(4, 8)))
    with pytest.raises(ValueError):
        decreasing_rules(3, length_four=1)


def test_export_rejects_tampered_definition_without_replacing_output(tmp_path):
    database = tmp_path / 'artin.sqlite3'
    output = tmp_path / 'counts.json'
    run = write_counts(database, maximum_rank=3, maximum_length=2)
    export_counts(database, run, output)
    original = output.read_bytes()
    with sqlite3.connect(database) as c:
        definition = json.loads(c.execute('SELECT definition FROM runs WHERE run_id=?', (run,)).fetchone()[0])
        definition['policy'] = 'tampered'
        c.execute('UPDATE runs SET definition=? WHERE run_id=?', (json.dumps(definition), run))
    with pytest.raises(ValueError, match='definition hash mismatch'):
        export_counts(database, run, output)
    assert output.read_bytes() == original
