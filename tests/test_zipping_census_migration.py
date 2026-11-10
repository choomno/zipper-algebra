from copy import deepcopy
from dataclasses import asdict
import sqlite3

import pytest

from experiments.zipper_census_migrate import history_digest, migrate, validate_coverage
from zipper_algebra.zipping_algorithm.census_queue import CensusQueue, SCHEMA
from zipper_algebra.zipping_algorithm.census_reuse import reconstruction_sources
from zipper_algebra.zipping_algorithm.interpreter import Limits


def seed(tmp_path):
    old = {'schema': SCHEMA, 'producer': {'test': 'old'}, 'producer_metadata': {'test': True},
           'limits': asdict(Limits()), 'mode': 'reuse', 'package_sources': reconstruction_sources()}
    new = deepcopy(old)
    new['producer'] = {'test': 'new'}
    queue = CensusQueue(tmp_path / 'q.sqlite3', old)
    for word, state in [([-1], 'certified'), ([-2], 'candidate'), ([1], 'computing'), ([2], 'queued')]:
        key = queue.enqueue(3, word)
        queue.db.execute('UPDATE async_jobs SET state=? WHERE key=?', (state, key))
        queue.db.execute('INSERT INTO async_attempts(policy,key,stage,status,record,created) VALUES (?,?,?,?,?,?)',
                         (queue.policy_key, key, 'compute', state, '{}', 1))
        if state == 'candidate':
            queue.db.execute('INSERT INTO async_candidates VALUES (?,?,?,?,?,?)', (queue.policy_key, key, b'\x00\xff', 'digest', 'output', 2))
    queue.db.execute('INSERT INTO async_cursors VALUES (?,?,?)', (queue.policy_key, 'layer', '4'))
    return queue, new


def test_preserves_history_partitions_and_is_idempotent(tmp_path):
    queue, new = seed(tmp_path)
    before = history_digest(queue.db, queue.policy_key)
    result = migrate(queue.path, queue.policy_key, new)
    assert result['source_history'] == before == history_digest(queue.db, queue.policy_key)
    assert result['retained_certified'] == 1
    assert result['queued_for_recomputation'] == 3
    assert tuple(queue.db.execute('SELECT state,count(*) FROM async_jobs WHERE policy=? GROUP BY state', (result['target'],)).fetchone()) == ('queued', 3)
    assert queue.db.execute('SELECT count(*) FROM async_candidates WHERE policy=?', (result['target'],)).fetchone()[0] == 0
    assert migrate(queue.path, queue.policy_key, new) == result
    queue.db.execute("UPDATE async_jobs SET state='certified' WHERE policy=?", (result['target'],))
    assert migrate(queue.path, queue.policy_key, new) == result
    queue.close()


def test_failure_rolls_back_policy_jobs_and_linkage(tmp_path):
    queue, new = seed(tmp_path)
    queue.db.execute("CREATE TRIGGER fail_new BEFORE INSERT ON async_jobs WHEN NEW.policy != '" + queue.policy_key + "' BEGIN SELECT RAISE(ABORT, 'injected'); END")
    before = history_digest(queue.db, queue.policy_key)
    with pytest.raises(Exception, match='injected'):
        migrate(queue.path, queue.policy_key, new)
    assert history_digest(queue.db, queue.policy_key) == before
    assert queue.db.execute('SELECT count(*) FROM async_policies').fetchone()[0] == 1
    assert queue.db.execute("SELECT count(*) FROM sqlite_master WHERE name='async_migrations'").fetchone()[0] == 0
    queue.close()


def test_detects_source_mutation_and_coverage_damage(tmp_path):
    queue, new = seed(tmp_path)
    result = migrate(queue.path, queue.policy_key, new)
    queue.db.execute('DELETE FROM async_jobs WHERE policy=?', (result['target'],))
    with pytest.raises(ValueError, match='coverage'):
        validate_coverage(queue.db, queue.policy_key, result['target'], result)
    with pytest.raises(sqlite3.IntegrityError, match='retired census'):
        queue.db.execute("UPDATE async_jobs SET updated=updated+1 WHERE policy=?", (queue.policy_key,))
    queue.db.execute(f'DROP TRIGGER retired_{queue.policy_key}_async_jobs_update')
    queue.db.execute("UPDATE async_jobs SET updated=updated+1 WHERE policy=?", (queue.policy_key,))
    with pytest.raises(ValueError, match='source epoch changed'):
        migrate(queue.path, queue.policy_key, new)
    queue.close()
