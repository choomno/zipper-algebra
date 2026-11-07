"""Essential representatives have distinct coverage, exact replay and resumability."""
from concurrent.futures import ThreadPoolExecutor
import json
import sqlite3
import zlib

import pytest

pytest.importorskip('fcntl', reason='experiment requires Unix flock')
from experiments import zipper_exhaustive_essential as runner
from zipper_algebra.zipping_algorithm.interpreter import Limits


@pytest.fixture(autouse=True)
def close_sessions():
    yield
    if runner.engine._verifier is not None:
        runner.engine._verifier.close()
    if runner.engine._session is not None:
        runner.engine._session.close()


def run(path, **kwargs):
    return runner.exhaust(path, workers=1, wall_seconds=30,
                          _executor_factory=ThreadPoolExecutor, **kwargs)


def rows(path):
    with sqlite3.connect(path / 'essential.sqlite3') as db:
        return db.execute('SELECT * FROM representatives ORDER BY n,length,ordinal').fetchall()


def test_resume_real_cpp_checks_and_distinct_legacy_store(tmp_path):
    legacy = tmp_path / 'words.sqlite3'
    legacy.write_bytes(b'legacy sentinel')
    first = run(tmp_path, targets={1: 0, 2: 2, 3: 2}, max_new_words=4)
    assert first['new_or_resumed_representatives'] == 4
    assert first['stop_reason'] == 'word_budget'
    previous = rows(tmp_path)
    result = run(tmp_path, targets={1: 0, 2: 2, 3: 2})
    assert result['status'] == 'completed'
    assert result['audit']['representatives'] == 4
    assert len(rows(tmp_path)) == 15
    assert legacy.read_bytes() == b'legacy sentinel'
    assert set(previous).issubset(rows(tmp_path))
    assert result['peak_in_flight'] <= 2
    assert result['coverage_scope'] == runner.LITERAL_COVERAGE
    assert result['identity']['commutations'] is False
    for row in rows(tmp_path):
        certificate = json.loads(zlib.decompress(row[6]))
        assert certificate['input_word'] == json.loads(row[3])
        assert row[13] == runner.engine.CPP_CHECK
        assert set(map(abs, certificate['input_word'])) == set(range(1, row[0]))


def test_real_spawn_and_far_commutation_counts(tmp_path):
    result = runner.exhaust(tmp_path, targets={4: 3}, workers=2, max_new_words=6,
                            wall_seconds=30, commutations=True)
    assert result['coverage_scope'] == runner.TRACE_COVERAGE
    assert result['new_or_resumed_representatives'] == 6
    assert result['ranks'][0]['target_representatives'] == 32
    assert result['peak_in_flight'] <= 4
    for row in rows(tmp_path):
        word = json.loads(row[3])
        assert runner._enumerator(4, 3, True).ordinal(word) == row[2]


def test_optional_trace_is_verified_and_bound(tmp_path):
    result = run(tmp_path, targets={2: 1}, trace=True)
    assert result['status'] == 'completed'
    assert all(json.loads(zlib.decompress(row[6]))['trace_enabled'] for row in rows(tmp_path))
    with pytest.raises(ValueError, match='identity changed'):
        run(tmp_path, targets={2: 1}, trace=False)


def test_incomplete_retry_archives_original_and_keeps_prefix(tmp_path):
    first = run(tmp_path, targets={2: 1, 3: 2}, limits=Limits(max_steps=1))
    assert first['status'] == 'paused'
    saved = rows(tmp_path)
    assert saved and all(row[5] == 'incomplete' for row in saved)
    repeated = run(tmp_path, targets={2: 1, 3: 2}, limits=Limits(max_steps=1))
    assert repeated['new_or_resumed_representatives'] == 0
    finished = run(tmp_path, targets={2: 1, 3: 2}, limits=Limits(max_steps=10000))
    assert finished['status'] == 'completed'
    with sqlite3.connect(tmp_path / 'essential.sqlite3') as db:
        archived = [row[:-1] for row in db.execute('SELECT * FROM attempts ORDER BY n,length,ordinal')]
        assert archived == saved
        assert runner.audit(db, finished['identity'])['archived_attempts'] == len(saved)


@pytest.mark.parametrize('mutation', [
    "UPDATE representatives SET word='[99]'",
    "UPDATE representatives SET certificate_sha='bad'",
    "UPDATE representatives SET certificate=x'00'",
    "UPDATE representatives SET certificate_bytes=0",
    "UPDATE representatives SET status='incomplete'",
    "UPDATE representatives SET metrics='{}'",
    "UPDATE representatives SET production_check='unchecked'",
    "DELETE FROM representatives WHERE ordinal=0",
    "UPDATE ranks SET ordinal=1",
    "UPDATE configurations SET config='{}'",
])
def test_audit_rejects_corruption(tmp_path, mutation):
    run(tmp_path, targets={2: 1})
    with sqlite3.connect(tmp_path / 'essential.sqlite3') as db:
        db.execute(mutation)
    with pytest.raises((ValueError, KeyError)):
        run(tmp_path, targets={2: 1})


def test_zero_budget_and_unpopulated_layers(tmp_path):
    result = run(tmp_path, targets={8: 6}, max_new_words=0)
    assert result['status'] == 'completed'
    assert result['new_or_resumed_representatives'] == 0
    assert result['ranks'][0]['target_representatives'] == 0
    result = run(tmp_path, targets={8: 7}, max_new_words=0)
    assert result['stop_reason'] == 'word_budget'
    assert result['status'] == 'paused'


def test_source_change_and_mode_change_rejected(tmp_path):
    run(tmp_path, targets={2: 1})
    with pytest.raises(ValueError, match='identity changed'):
        run(tmp_path, targets={2: 1}, commutations=True)
    with sqlite3.connect(tmp_path / 'essential.sqlite3') as db:
        expected = json.loads(db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0])
    expected['enumeration_sha256'] = 'bad'
    with pytest.raises(ValueError, match='source changed'):
        runner._assert_source(expected)


@pytest.mark.parametrize('kwargs', [
    {'workers': 0}, {'trace': 1}, {'commutations': 'yes'},
    {'wall_seconds': 0}, {'cache_bytes': -1}, {'max_new_words': -1},
    {'targets': {1: 1}}, {'targets': {0: 1}}, {'limits': {}},
])
def test_invalid_options_leave_no_database(tmp_path, kwargs):
    with pytest.raises(ValueError):
        runner.exhaust(tmp_path, **kwargs)
    assert not (tmp_path / 'essential.sqlite3').exists()


def test_failed_worker_records_failure_without_advancing_cursor(tmp_path, monkeypatch):
    actual = runner._compute
    def fail(n, word):
        if n == 2:
            raise RuntimeError('worker failure')
        return actual(n, word)
    monkeypatch.setattr(runner, '_compute', fail)
    result = run(tmp_path, targets={2: 1, 3: 2})
    assert result['status'] == 'paused'
    assert result['new_failures'] == 1
    assert result['failure_records'][0]['diagnostic']['message'] == 'worker failure'
    with sqlite3.connect(tmp_path / 'essential.sqlite3') as db:
        assert db.execute('SELECT * FROM ranks WHERE n=2').fetchall() == [(2, 1, 0)]
        failure = json.loads(db.execute('SELECT failure FROM failures').fetchone()[0])
        assert failure['word'] == [-1] and failure['stage'] == 'worker'
    assert all(row[0] == 3 for row in rows(tmp_path))
    assert len(rows(tmp_path)) == 8
    monkeypatch.setattr(runner, '_compute', actual)
    assert run(tmp_path, targets={2: 1})['new_or_resumed_representatives'] == 0
    assert run(tmp_path, targets={2: 1}, retry_failures=True)['status'] == 'completed'
    assert runner.audit_store(tmp_path)['audit']['failures'] == 1


def test_read_only_audit_leaves_database_unchanged(tmp_path):
    run(tmp_path, targets={2: 1})
    path = tmp_path / 'essential.sqlite3'
    before = path.read_bytes()
    result = runner.audit_store(tmp_path)
    assert result['audit']['representatives'] == 2
    assert result['audit']['replayed'] == 0
    assert path.read_bytes() == before


def test_replay_audit_rejects_rehashed_wrong_output(tmp_path):
    run(tmp_path, targets={2: 1})
    assert runner.audit_store(tmp_path, replay=True)['audit']['replayed'] == 2
    with sqlite3.connect(tmp_path / 'essential.sqlite3') as db:
        row = db.execute('SELECT certificate FROM representatives WHERE ordinal=0').fetchone()
        certificate = json.loads(zlib.decompress(row[0]))
        certificate['output'][0]['coefficient'][0][1] += 1
        encoded = runner.engine._json(certificate)
        raw = encoded.encode()
        db.execute('UPDATE representatives SET certificate=?,certificate_sha=?,certificate_bytes=? WHERE ordinal=0',
                   (zlib.compress(raw), runner.engine._sha(encoded), len(raw)))
    assert runner.audit_store(tmp_path)['audit']['replayed'] == 0
    with pytest.raises(ValueError, match='independent replay rejected'):
        runner.audit_store(tmp_path, replay=True)


def test_verifier_rejection_retains_certificate_and_input(tmp_path, monkeypatch):
    actual = runner._initialize
    def initialize(*args):
        actual(*args)
        runner.engine._verifier.check = lambda certificate: False
    monkeypatch.setattr(runner, '_initialize', initialize)
    result = run(tmp_path, targets={2: 1})
    assert result['new_failures'] == 1
    assert result['new_or_resumed_representatives'] == 0
    assert not rows(tmp_path)
    with sqlite3.connect(tmp_path / 'essential.sqlite3') as db:
        row = db.execute('SELECT * FROM failures').fetchone()
        failure = json.loads(row[5])
        certificate = json.loads(zlib.decompress(row[7]))
        assert failure['stage'] == 'verifier'
        assert certificate['input_word'] == failure['word'] == [-1]
        assert certificate['n'] == failure['n'] == 2
        assert db.execute('PRAGMA journal_mode').fetchone()[0] == 'wal'
    assert runner.audit_store(tmp_path)['audit']['failures'] == 1


def test_summary_counts_only_current_target_and_replay_archives(tmp_path):
    run(tmp_path, targets={2: 2}, limits=Limits(max_steps=1))
    result = run(tmp_path, targets={2: 2}, limits=Limits(max_steps=10000))
    assert result['ranks'][0]['complete_representatives'] == 6
    assert runner.audit_store(tmp_path, replay=True)['audit']['replayed'] == 7
    narrowed = run(tmp_path, targets={2: 1})
    assert narrowed['ranks'][0]['complete_representatives'] == 2
    assert narrowed['ranks'][0]['target_representatives'] == 2


def test_failure_diagnostic_and_certificate_tampering_rejected(tmp_path, monkeypatch):
    actual = runner._initialize
    def initialize(*args):
        actual(*args)
        runner.engine._verifier.check = lambda certificate: False
    monkeypatch.setattr(runner, '_initialize', initialize)
    run(tmp_path, targets={2: 1})
    with sqlite3.connect(tmp_path / 'essential.sqlite3') as db:
        db.execute("UPDATE failures SET failure='{}'")
    with pytest.raises(ValueError, match='failure diagnostic hash'):
        runner.audit_store(tmp_path)


@pytest.mark.parametrize('replay', [False, True])
def test_read_only_audit_pins_wal_snapshot_and_closes_connection(tmp_path, monkeypatch, replay):
    run(tmp_path, targets={2: 1})
    actual = runner.audit
    connections = []
    def concurrent_write(db, expected, **kwargs):
        assert db.in_transaction
        connections.append(db)
        # Commit after audit_store read the identity. The audit must continue
        # seeing the same database version, even while independently replaying.
        with sqlite3.connect(tmp_path / 'essential.sqlite3') as writer:
            writer.execute("UPDATE representatives SET word='[99]' WHERE ordinal=0")
        return actual(db, expected, **kwargs)
    monkeypatch.setattr(runner, 'audit', concurrent_write)
    result = runner.audit_store(tmp_path, replay=replay)
    assert result['audit']['representatives'] == 2
    assert result['audit']['replayed'] == (2 if replay else 0)
    with pytest.raises(sqlite3.ProgrammingError, match='closed'):
        connections[0].execute('SELECT 1')
    with sqlite3.connect(tmp_path / 'essential.sqlite3') as writer:
        assert writer.execute('SELECT word FROM representatives WHERE ordinal=0').fetchone()[0] == '[99]'


def test_connect_failure_releases_writer_lock(tmp_path, monkeypatch):
    actual = runner.sqlite3.connect
    def fail(*args, **kwargs):
        raise sqlite3.OperationalError('cannot open database')
    monkeypatch.setattr(runner.sqlite3, 'connect', fail)
    with pytest.raises(sqlite3.OperationalError, match='cannot open database'):
        run(tmp_path, targets={2: 1})
    monkeypatch.setattr(runner.sqlite3, 'connect', actual)
    assert run(tmp_path, targets={2: 1})['status'] == 'completed'


def test_default_keeps_every_connected_literal_order(tmp_path):
    result = run(tmp_path, targets={4: 3}, max_new_words=0)
    assert result['ranks'][0]['target_representatives'] == 48
    assert result['identity']['commutations'] is False
    words = runner._enumerator(4, 3, False)
    assert words.word_at(words.ordinal((1, 3, 2))) == (1, 3, 2)
    assert words.word_at(words.ordinal((3, 1, 2))) == (3, 1, 2)


@pytest.mark.parametrize(('flags', 'commutations'), [
    ([], False), (['--retain-commuting-orders'], False), (['--commutation-classes'], True),
])
def test_cli_commutation_quotient_requires_explicit_opt_in(tmp_path, monkeypatch, flags, commutations):
    captured = {}
    def exhaust(output, **kwargs):
        captured.update(kwargs)
        return {}
    monkeypatch.setattr(runner, 'exhaust', exhaust)
    monkeypatch.setattr('sys.argv', ['zipper_exhaustive_essential', str(tmp_path), *flags])
    runner.main()
    assert captured['commutations'] is commutations
