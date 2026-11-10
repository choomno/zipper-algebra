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


def test_refuses_forked_target_and_new_admissions(tmp_path):
    queue, new = seed(tmp_path)
    result = migrate(queue.path, queue.policy_key, new)
    fork = deepcopy(new)
    fork['producer'] = {'test': 'other'}
    with pytest.raises(ValueError, match='another target'):
        migrate(queue.path, queue.policy_key, fork)
    target = CensusQueue(queue.path, new)
    with pytest.raises(sqlite3.IntegrityError, match='original uncertified'):
        target.enqueue(3, [1, 2])
    with pytest.raises(sqlite3.IntegrityError, match='original uncertified'):
        target.enqueue(3, [-1])
    assert target.db.execute('SELECT count(*) FROM async_cursors WHERE policy=?', (result['target'],)).fetchone()[0] == 1
    target.close()
    queue.close()


def test_rejects_incomplete_declared_frontier_and_count(tmp_path):
    queue, new = seed(tmp_path)
    with pytest.raises(ValueError, match='frontier is incomplete'):
        migrate(queue.path, queue.policy_key, new, expected_count=4, completed_layers=((3, 9),))
    with pytest.raises(ValueError, match='literal count'):
        migrate(queue.path, queue.policy_key, new, expected_count=5)
    assert queue.db.execute('SELECT count(*) FROM async_policies').fetchone()[0] == 1
    queue.close()


def chain_seed(tmp_path):
    from zipper_algebra.zipping_algorithm.peeling_enumeration import NonPeelableWords
    old = {'schema': SCHEMA, 'producer': {'test': 'root'}, 'producer_metadata': {'test': True},
           'limits': asdict(Limits()), 'mode': 'reuse', 'package_sources': reconstruction_sources()}
    root = CensusQueue(tmp_path / 'chain.sqlite3', old)
    layer = NonPeelableWords(3, 2)
    root.enqueue_interval('nonpeelable-v1:3:2', layer.count,
                          lambda i: (3, list(layer.word_at(i))), batch=layer.count)
    root.db.execute("UPDATE async_jobs SET state='certified' WHERE policy=? AND rowid=(SELECT min(rowid) FROM async_jobs)", (root.policy_key,))
    next_policy = deepcopy(old)
    next_policy['producer'] = {'test': 'one'}
    first = migrate(root.path, root.policy_key, next_policy, expected_count=5, completed_layers=((3, 2),))
    return root, next_policy, first


def advance_tail(db, policy):
    key = db.execute('SELECT key FROM async_jobs WHERE policy=? ORDER BY rowid LIMIT 1', (policy,)).fetchone()[0]
    db.execute("UPDATE async_jobs SET state='certified' WHERE policy=? AND key=?", (policy, key))
    db.execute('INSERT INTO async_attempts(policy,key,stage,status,record,created) VALUES (?,?,?,?,?,?)',
               (policy, key, 'certify', 'certified', '{"provenance":"preserve"}', 1))


def test_multi_hop_preserves_all_history_and_disjoint_scope(tmp_path):
    import json
    root, one, first = chain_seed(tmp_path)
    advance_tail(root.db, first['target'])
    histories = {key: history_digest(root.db, key) for key in (root.policy_key, first['target'])}
    two = deepcopy(one)
    two['producer'] = {'test': 'two'}
    second = migrate(root.path, first['target'], two, expected_count=4, inherit_closed_scope=True)
    assert all(history_digest(root.db, key) == value for key, value in histories.items())
    inherited = second['inherited_closed_scope']
    assert inherited['origin_source'] == root.policy_key
    assert inherited['origin_literal_requests'] == 5
    assert inherited['ancestor_retained_certified'] == 1
    assert second['retained_certified'] == 1
    assert second['queued_for_recomputation'] == 3
    advance_tail(root.db, second['target'])
    assert migrate(root.path, first['target'], two, expected_count=4, inherit_closed_scope=True) == second
    three = deepcopy(two)
    three['producer'] = {'test': 'three'}
    third = migrate(root.path, second['target'], three, expected_count=3, inherit_closed_scope=True)
    assert third['inherited_closed_scope']['ancestor_retained_certified'] == 2
    assert len(third['inherited_closed_scope']['edges']) == 2
    assert third['queued_for_recomputation'] == 2
    assert third['inherited_closed_scope']['ancestor_retained_certified'] + third['retained_certified'] + third['queued_for_recomputation'] == 5
    assert histories[root.policy_key] == history_digest(root.db, root.policy_key)
    assert histories[first['target']] == history_digest(root.db, first['target'])
    assert json.loads(root.db.execute('SELECT record FROM async_attempts WHERE policy=?', (first['target'],)).fetchone()[0]) == {'provenance': 'preserve'}
    root.close()


@pytest.mark.parametrize('damage', ['missing_parent', 'coverage', 'guard', 'manifest', 'overlap', 'omission', 'literal', 'frontier', 'ancestor_history'])
def test_rejects_corrupt_inherited_scope_without_changes(tmp_path, damage):
    import json
    from experiments.zipper_census_migrate import closure_sql
    root, one, first = chain_seed(tmp_path)
    target = first['target']
    if damage == 'missing_parent':
        root.db.execute('PRAGMA foreign_keys=OFF')
        root.db.execute('DELETE FROM async_migrations')
    elif damage == 'coverage':
        root.db.execute('DELETE FROM async_migration_coverage WHERE rowid=(SELECT min(rowid) FROM async_migration_coverage)')
    elif damage == 'guard':
        root.db.execute(f'DROP TRIGGER retired_{root.policy_key}_async_jobs_update')
    elif damage == 'manifest':
        manifest = deepcopy(first)
        manifest['retained_certified'] += 1
        root.db.execute('UPDATE async_migrations SET manifest=?', (json.dumps(manifest),))
    elif damage == 'overlap':
        root.db.execute(f'DROP TRIGGER closed_migration_{target}')
        root.db.execute("INSERT INTO async_jobs SELECT ?,key,n,word,'queued',NULL,NULL,created,updated FROM async_jobs WHERE policy=? AND state='certified'", (target, root.policy_key))
        root.db.execute(closure_sql(root.policy_key, target))
    elif damage == 'omission':
        root.db.execute('DELETE FROM async_jobs WHERE policy=? AND key=(SELECT key FROM async_jobs WHERE policy=? LIMIT 1)', (target, target))
    elif damage == 'literal':
        root.db.execute("UPDATE async_jobs SET word='[]' WHERE policy=?", (target,))
    elif damage == 'frontier':
        root.db.execute(f'DROP TRIGGER retired_{root.policy_key}_async_cursors_update')
        root.db.execute("UPDATE async_cursors SET ordinal='4' WHERE policy=?", (root.policy_key,))
    elif damage == 'ancestor_history':
        root.db.execute(f'DROP TRIGGER retired_{root.policy_key}_async_jobs_update')
        root.db.execute('UPDATE async_jobs SET updated=0 WHERE policy=?', (root.policy_key,))
    before = history_digest(root.db, target)
    two = deepcopy(one)
    two['producer'] = {'test': 'two'}
    with pytest.raises(ValueError):
        migrate(root.path, target, two, inherit_closed_scope=True)
    assert before == history_digest(root.db, target)
    assert root.db.execute('SELECT count(*) FROM async_policies').fetchone()[0] == 2
    root.close()


def test_inherited_scope_rejects_wrong_mode_limits_count_and_bypass(tmp_path):
    root, one, first = chain_seed(tmp_path)
    two = deepcopy(one)
    two['producer'] = {'test': 'two'}
    with pytest.raises(ValueError, match='literal count'):
        migrate(root.path, first['target'], two, expected_count=5, inherit_closed_scope=True)
    with pytest.raises(ValueError, match='requires inherited'):
        migrate(root.path, first['target'], two)
    for field, value in [('mode', 'standalone'), ('limits', {'max_steps': 1})]:
        invalid = deepcopy(two)
        invalid[field] = value
        with pytest.raises(ValueError, match='limits and execution mode'):
            migrate(root.path, first['target'], invalid, inherit_closed_scope=True)
    with pytest.raises(ValueError, match='mutually exclusive'):
        migrate(root.path, first['target'], two, completed_layers=((3, 2),), inherit_closed_scope=True)
    root.close()


def test_inherited_migration_rolls_back_every_change(tmp_path):
    root, one, first = chain_seed(tmp_path)
    advance_tail(root.db, first['target'])
    two = deepcopy(one)
    two['producer'] = {'test': 'two'}
    from zipper_algebra.zipping_algorithm.census_reuse import _sha
    target = _sha(two)
    histories = {key: history_digest(root.db, key) for key in (root.policy_key, first['target'])}
    root.db.execute(f"CREATE TRIGGER fail_chain BEFORE INSERT ON async_jobs WHEN NEW.policy='{target}' BEGIN SELECT RAISE(ABORT, 'injected chained failure'); END")
    with pytest.raises(sqlite3.IntegrityError, match='injected chained'):
        migrate(root.path, first['target'], two, expected_count=4, inherit_closed_scope=True)
    assert all(history_digest(root.db, key) == digest for key, digest in histories.items())
    assert root.db.execute('SELECT count(*) FROM async_policies').fetchone()[0] == 2
    assert root.db.execute('SELECT count(*) FROM async_migrations').fetchone()[0] == 1
    assert root.db.execute("SELECT count(*) FROM sqlite_master WHERE type='trigger' AND name=?", (f'retired_{first["target"]}_async_jobs_update',)).fetchone()[0] == 0
    root.close()


def test_cli_scope_flags_are_mutually_exclusive(monkeypatch):
    from experiments.zipper_census_migrate import main
    monkeypatch.setattr('sys.argv', ['migrate', 'unused.sqlite3', '--source-policy', 'source',
                                    '--expected-literal-count', '4', '--expected-target', 'prepared', '--completed-layer', '3:2',
                                    '--inherit-closed-scope'])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2


def test_inherited_fully_certified_tail_has_no_recomputation(tmp_path):
    root, one, first = chain_seed(tmp_path)
    root.db.execute("UPDATE async_jobs SET state='certified' WHERE policy=?", (first['target'],))
    two = deepcopy(one)
    two['producer'] = {'test': 'two'}
    second = migrate(root.path, first['target'], two, expected_count=4, inherit_closed_scope=True)
    assert second['queued_for_recomputation'] == 0
    assert second['inherited_closed_scope']['ancestor_retained_certified'] + second['retained_certified'] == 5
    assert migrate(root.path, first['target'], two, expected_count=4, inherit_closed_scope=True) == second
    root.close()


def test_idempotence_refuses_missing_guards_without_row_changes(tmp_path):
    root, one, first = chain_seed(tmp_path)
    before = history_digest(root.db, root.policy_key)
    root.db.execute(f'DROP TRIGGER closed_migration_{first["target"]}')
    with pytest.raises(ValueError, match='guard missing or changed'):
        migrate(root.path, root.policy_key, one, expected_count=5, completed_layers=((3, 2),))
    assert before == history_digest(root.db, root.policy_key)
    root.close()


def test_expected_target_and_external_ancestry_anchors(tmp_path):
    import hashlib
    from zipper_algebra.zipping_algorithm.census_reuse import _json, _sha
    root, one, first = chain_seed(tmp_path)
    two = deepcopy(one)
    two['producer'] = {'test': 'two'}
    digest = hashlib.sha256(_json(first).encode()).hexdigest()
    for kwargs, message in [({'expected_target': 'wrong'}, 'target policy differs'),
                            ({'expected_parent_manifest_sha256': 'wrong'}, 'parent manifest differs'),
                            ({'expected_origin': 'wrong'}, 'origin differs')]:
        with pytest.raises(ValueError, match=message):
            migrate(root.path, first['target'], two, inherit_closed_scope=True, **kwargs)
    second = migrate(root.path, first['target'], two, expected_count=4, inherit_closed_scope=True,
                     expected_target=_sha(two), expected_parent_manifest_sha256=digest,
                     expected_origin=root.policy_key)
    assert second['inherited_closed_scope']['edges'][0]['manifest_sha256'] == digest
    root.close()


def test_origin_membership_detects_count_preserving_foreign_literal(tmp_path):
    from experiments.zipper_census_migrate import validate_origin_membership
    root, one, first = chain_seed(tmp_path)
    root.db.execute(f'DROP TRIGGER retired_{root.policy_key}_async_jobs_update')
    root.db.execute("UPDATE async_jobs SET word='[]' WHERE policy=? AND key=(SELECT key FROM async_jobs WHERE policy=? LIMIT 1)", (root.policy_key, root.policy_key))
    with pytest.raises(ValueError, match='literal membership differs'):
        validate_origin_membership(root.db, root.policy_key, first['completed_layers'])
    root.close()


def test_idempotence_authenticates_target_payload(tmp_path):
    root, one, first = chain_seed(tmp_path)
    root.db.execute("UPDATE async_policies SET payload='{}' WHERE key=?", (first['target'],))
    with pytest.raises(ValueError, match='target policy binding'):
        migrate(root.path, root.policy_key, one, expected_count=5, completed_layers=((3, 2),))
    root.close()


@pytest.mark.parametrize('field,value', [('schema', 'wrong'), ('source', 'wrong'), ('target', 'wrong'),
                                         ('source_retired', False), ('source_states', {}), ('literal_requests', 6)])
def test_idempotence_authenticates_existing_manifest(tmp_path, field, value):
    import json
    root, one, first = chain_seed(tmp_path)
    corrupted = deepcopy(first)
    corrupted[field] = value
    root.db.execute('UPDATE async_migrations SET manifest=?', (json.dumps(corrupted),))
    with pytest.raises(ValueError, match='manifest identity or source counts'):
        migrate(root.path, root.policy_key, one, expected_count=5, completed_layers=((3, 2),))
    root.close()


def test_origin_membership_detects_real_count_preserving_substitution(tmp_path):
    from experiments.zipper_census_migrate import validate_origin_membership
    from zipper_algebra.zipping_algorithm.peeling_enumeration import NonPeelableWords
    old = {'schema': SCHEMA, 'producer': {'test': 'root'}, 'producer_metadata': {'test': True},
           'limits': asdict(Limits()), 'mode': 'reuse', 'package_sources': reconstruction_sources()}
    root = CensusQueue(tmp_path / 'substitution.sqlite3', old)
    layer = NonPeelableWords(3, 2)
    for i in range(layer.count):
        root.enqueue(3, list(layer.word_at(i)))
    root.db.execute('DELETE FROM async_jobs WHERE rowid=(SELECT min(rowid) FROM async_jobs)')
    root.enqueue(3, [])
    assert root.db.execute('SELECT count(*) FROM async_jobs').fetchone()[0] == layer.count
    with pytest.raises(ValueError, match='literal membership differs'):
        validate_origin_membership(root.db, root.policy_key, [{'n': 3, 'length': 2, 'ordinal': layer.count}])
    root.close()


def test_production_origin_count_is_pinned():
    from zipper_algebra.zipping_algorithm.peeling_enumeration import NonPeelableWords
    assert NonPeelableWords(3, 9).count == 252497


def test_cli_inherited_scope_requires_external_anchors(monkeypatch):
    from experiments.zipper_census_migrate import main
    monkeypatch.setattr('sys.argv', ['migrate', 'unused.sqlite3', '--source-policy', 'source',
                                    '--expected-literal-count', '4', '--expected-target', 'prepared',
                                    '--inherit-closed-scope'])
    with pytest.raises(SystemExit) as error:
        main()
    assert error.value.code == 2
