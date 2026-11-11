"""Preserve an async census epoch and queue its uncertified literal requests.

Stop and authenticate all source workers before invoking this command. Migration
holds one SQLite write transaction; it does not interrupt or retire workers.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import closing
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import time
from typing import Any

from experiments.zipper_async_census import initialize, read_policy
from zipper_algebra.zipping_algorithm.census_queue import SCHEMA
from zipper_algebra.zipping_algorithm.census_reuse import _json, _sha, reconstruction_sources
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.peeling_enumeration import NonPeelableWords

TABLES = ('async_policies', 'async_jobs', 'async_candidates', 'async_attempts', 'async_cursors')


def history_digest(db: sqlite3.Connection, policy: str) -> dict[str, Any]:
    """Hash every old row, including binary candidate payloads, deterministically."""
    result = {}
    for table in TABLES:
        column = 'key' if table == 'async_policies' else 'policy'
        digest = hashlib.sha256()
        count = 0
        for row in db.execute(f'SELECT * FROM {table} WHERE {column}=? ORDER BY rowid', (policy,)):
            for value in row:
                data = value if isinstance(value, bytes) else _json(value).encode()
                digest.update(type(value).__name__.encode() + b'\x00')
                digest.update(len(data).to_bytes(8, 'big'))
                digest.update(data)
            count += 1
        result[table] = {'rows': count, 'sha256': digest.hexdigest()}
    return result


def migrate(path: Path, source: str, policy: dict[str, Any], *, expected_count: int | None = None, completed_layers: tuple[tuple[int, int], ...] = (), inherit_closed_scope: bool = False, expected_target: str | None = None, expected_parent_manifest_sha256: str | None = None, expected_origin: str | None = None) -> dict[str, Any]:
    if inherit_closed_scope and completed_layers:
        raise ValueError('inherited scope and completed layers are mutually exclusive')
    target = _sha(policy)
    if expected_target is not None and target != expected_target:
        raise ValueError('target policy differs from independently prepared target')
    if source == target or policy.get('schema') != SCHEMA or policy.get('package_sources') != reconstruction_sources() or not policy.get('producer') or not policy.get('producer_metadata') or not policy.get('limits'):
        raise ValueError('a distinct policy bound to current package sources is required')
    with closing(sqlite3.connect(path.resolve().as_uri() + "?mode=rw", uri=True, timeout=60, isolation_level=None)) as db:
        db.execute('PRAGMA foreign_keys=ON')
        db.execute('BEGIN IMMEDIATE')
        try:
            stored = db.execute('SELECT payload FROM async_policies WHERE key=?', (source,)).fetchone()
            if stored is None or _sha(json.loads(stored[0])) != source:
                raise ValueError('source policy missing or digest mismatch')
            old_policy = json.loads(stored[0])
            if _json(policy['limits']) != _json(old_policy['limits']) or _json(policy['mode']) != _json(old_policy['mode']):
                raise ValueError('epoch switch must preserve limits and execution mode')
            migration_table = db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='async_migrations'").fetchone()
            parent = db.execute('SELECT 1 FROM async_migrations WHERE target=?', (source,)).fetchone() if migration_table else None
            if parent and not inherit_closed_scope:
                raise ValueError('migrated source requires inherited closed scope')
            inherited = inherited_scope(db, source) if inherit_closed_scope else None
            if inherited is not None:
                if expected_origin is not None and inherited['origin_source'] != expected_origin:
                    raise ValueError('origin differs from external anchor')
                if expected_parent_manifest_sha256 is not None and inherited['edges'][-1]['manifest_sha256'] != expected_parent_manifest_sha256:
                    raise ValueError('parent manifest differs from external anchor')
            frontier = []
            for n, length in completed_layers:
                layer = NonPeelableWords(n, length)
                cursor = db.execute('SELECT ordinal FROM async_cursors WHERE policy=? AND layer=?', (source, f'nonpeelable-v1:{n}:{length}')).fetchone()
                if cursor is None or int(cursor[0]) != layer.count:
                    raise ValueError('source enumeration frontier is incomplete')
                frontier.append({'n': n, 'length': length, 'ordinal': layer.count})
            source_count = db.execute('SELECT count(*) FROM async_jobs WHERE policy=?', (source,)).fetchone()[0]
            if expected_count is not None and source_count != expected_count:
                raise ValueError('source literal count differs from declared closed census')
            if frontier and sum(item['ordinal'] for item in frontier) != source_count:
                raise ValueError('completed layers do not cover the declared closed census')
            if frontier:
                validate_origin_membership(db, source, frontier)
            db.execute('CREATE TABLE IF NOT EXISTS async_migrations(source TEXT NOT NULL UNIQUE REFERENCES async_policies(key), target TEXT NOT NULL UNIQUE REFERENCES async_policies(key), manifest TEXT NOT NULL, PRIMARY KEY(source,target))')
            db.execute('CREATE TABLE IF NOT EXISTS async_migration_coverage(source TEXT NOT NULL,target TEXT NOT NULL,key TEXT NOT NULL,old_state TEXT NOT NULL,disposition TEXT NOT NULL,PRIMARY KEY(source,target,key),FOREIGN KEY(source,target) REFERENCES async_migrations(source,target),FOREIGN KEY(source,key) REFERENCES async_jobs(policy,key))')
            existing = db.execute('SELECT manifest FROM async_migrations WHERE source=? AND target=?', (source, target)).fetchone()
            before = history_digest(db, source)
            if existing:
                target_row = db.execute('SELECT payload FROM async_policies WHERE key=?', (target,)).fetchone()
                if target_row is None or target_row[0] != _json(policy):
                    raise ValueError('stored target policy binding mismatch')
                manifest = json.loads(existing[0])
                if manifest['source_history'] != before:
                    raise ValueError('source epoch changed after migration')
                validate_coverage(db, source, target, manifest)
                verify_guards(db, source, target)
                if manifest.get('completed_layers') != frontier:
                    raise ValueError('migration completed frontier differs')
                if manifest.get('inherited_closed_scope') != inherited:
                    raise ValueError('migration inherited scope differs')
                db.execute('COMMIT')
                return manifest
            if db.execute('SELECT 1 FROM async_migrations WHERE source=?', (source,)).fetchone():
                raise ValueError('source epoch already migrated to another target')
            if db.execute('SELECT 1 FROM async_policies WHERE key=?', (target,)).fetchone():
                raise ValueError('target policy already exists without migration linkage')
            now = time.time()
            counts = dict(db.execute('SELECT state,count(*) FROM async_jobs WHERE policy=? GROUP BY state', (source,)))
            manifest = {'schema': 'zipper-async-migration-v1', 'source': source, 'target': target,
                        'created': now, 'source_history': before, 'source_states': counts,
                        'scope': 'closed source admitted literal set; no further admission',
                        'completed_layers': frontier, 'inherited_closed_scope': inherited, 'source_retired': True, 'source_guard': 'SQLite triggers reject all old-policy row mutations',
                        'literal_requests': sum(counts.values()), 'retained_certified': counts.get('certified', 0),
                        'queued_for_recomputation': sum(v for k, v in counts.items() if k != 'certified'),
                        'coverage': 'Each source literal is retained as source certified evidence or queued once under target; provisional candidates are recomputed.'}
            db.execute('INSERT INTO async_policies VALUES (?,?)', (target, _json(policy)))
            db.execute('INSERT INTO async_migrations VALUES (?,?,?)', (source, target, _json(manifest)))
            db.execute("INSERT INTO async_migration_coverage SELECT policy,?,key,state,CASE WHEN state='certified' THEN 'retained_source_certificate' ELSE 'recompute_target' END FROM async_jobs WHERE policy=?", (target, source))
            db.execute("INSERT INTO async_jobs SELECT ?,key,n,word,'queued',NULL,NULL,?,? FROM async_jobs WHERE policy=? AND state!='certified' ORDER BY created,rowid", (target, now, now, source))
            db.execute('INSERT INTO async_cursors SELECT ?,layer,ordinal FROM async_cursors WHERE policy=?', (target, source))
            validate_coverage(db, source, target, manifest)
            if inherited is not None and inherited != inherited_scope(db, source):
                raise ValueError('ancestor history changed during migration')
            if before != history_digest(db, source):
                raise ValueError('source history changed during migration')
            retire_source(db, source)
            db.execute(closure_sql(source, target))
            db.execute('COMMIT')
            return manifest
        except BaseException:
            if db.in_transaction:
                db.execute('ROLLBACK')
            raise


def retirement_sql(source: str) -> list[str]:
    statements = []
    for table in TABLES:
        column = 'key' if table == 'async_policies' else 'policy'
        for operation in ('INSERT', 'UPDATE', 'DELETE'):
            refs = ('OLD', 'NEW') if operation == 'UPDATE' else (('OLD',) if operation == 'DELETE' else ('NEW',))
            condition = ' OR '.join(f"{ref}.{column}='{source}'" for ref in refs)
            name = f'retired_{source}_{table}_{operation.lower()}'
            statements.append(f"CREATE TRIGGER {name} BEFORE {operation} ON {table} WHEN {condition} BEGIN SELECT RAISE(ABORT, 'retired census source policy is immutable'); END")
    return statements


def closure_sql(source: str, target: str) -> str:
    return f"CREATE TRIGGER closed_migration_{target} BEFORE INSERT ON async_jobs WHEN NEW.policy='{target}' AND NOT EXISTS (SELECT 1 FROM async_migration_coverage WHERE source='{source}' AND target=NEW.policy AND key=NEW.key AND disposition='recompute_target') BEGIN SELECT RAISE(ABORT, 'migration target accepts only original uncertified literals'); END"


def retire_source(db: sqlite3.Connection, source: str) -> None:
    """Reject accidental restart or modification of any retired policy row."""
    for statement in retirement_sql(source):
        db.execute(statement)


def verify_guards(db: sqlite3.Connection, source: str, target: str) -> None:
    for statement in retirement_sql(source) + [closure_sql(source, target)]:
        name = statement.split()[2]
        row = db.execute("SELECT sql FROM sqlite_master WHERE type='trigger' AND name=?", (name,)).fetchone()
        if row is None or row[0] != statement:
            raise ValueError('migration retirement or closure guard missing or changed')


def validate_origin_membership(db: sqlite3.Connection, origin: str, layers: list[dict[str, int]]) -> None:
    """Reconstruct the declared enumeration and match each exact literal row."""
    if len({(item['n'], item['length']) for item in layers}) != len(layers):
        raise ValueError('origin completed layers are duplicated')
    remaining = {key: (n, word) for key, n, word in db.execute('SELECT key,n,word FROM async_jobs WHERE policy=?', (origin,))}
    for item in layers:
        layer = NonPeelableWords(item['n'], item['length'])
        for ordinal in range(layer.count):
            word = list(layer.word_at(ordinal))
            key = _sha({'n': layer.n, 'word': word})
            if remaining.pop(key, None) != (layer.n, _json(word)):
                raise ValueError('origin literal membership differs from enumeration')
    if remaining:
        raise ValueError('origin contains literals outside declared enumeration')


def inherited_scope(db: sqlite3.Connection, source: str) -> dict[str, Any]:
    """Authenticate every historical partition back to an exhausted origin."""
    tables = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if not {'async_migrations', 'async_migration_coverage'} <= tables:
        raise ValueError('inherited closed scope requires a parent migration')
    edges = []
    seen = {source}
    cursor = source
    while True:
        rows = db.execute('SELECT source,target,manifest FROM async_migrations WHERE target=?', (cursor,)).fetchall()
        if not rows:
            break
        if len(rows) != 1:
            raise ValueError('ambiguous parent migration')
        parent, target, payload = rows[0]
        if parent in seen:
            raise ValueError('migration ancestry cycle')
        seen.add(parent)
        manifest = json.loads(payload)
        if manifest.get('schema') != 'zipper-async-migration-v1' or manifest.get('source') != parent or manifest.get('target') != target or manifest.get('source_retired') is not True:
            raise ValueError('parent migration identity or retirement mismatch')
        policies = []
        for key in (parent, target):
            row = db.execute('SELECT payload FROM async_policies WHERE key=?', (key,)).fetchone()
            policy = json.loads(row[0]) if row else {}
            if _sha(policy) != key or policy.get('schema') != SCHEMA:
                raise ValueError('parent policy missing or digest mismatch')
            policies.append(policy)
        if any(_json(policies[0].get(field)) != _json(policies[1].get(field)) for field in ('limits', 'mode')):
            raise ValueError('parent migration limits or mode mismatch')
        if history_digest(db, parent) != manifest.get('source_history'):
            raise ValueError('ancestor source history changed')
        states = dict(db.execute('SELECT state,count(*) FROM async_jobs WHERE policy=? GROUP BY state', (parent,)))
        if states != manifest.get('source_states') or sum(states.values()) != manifest.get('literal_requests'):
            raise ValueError('parent source counts differ')
        validate_coverage(db, parent, target, manifest)
        verify_guards(db, parent, target)
        edges.append((parent, target, manifest, payload))
        cursor = parent
    if not edges:
        raise ValueError('inherited closed scope requires a parent migration')
    scope: dict[str, Any] = {}
    for parent, target, manifest, payload in reversed(edges):
        if not scope:
            layers = manifest.get('completed_layers')
            if not layers or manifest.get('inherited_closed_scope'):
                raise ValueError('origin has no completed frontier proof')
            for layer in layers:
                n, length = layer['n'], layer['length']
                count = NonPeelableWords(n, length).count
                row = db.execute('SELECT ordinal FROM async_cursors WHERE policy=? AND layer=?', (parent, f'nonpeelable-v1:{n}:{length}')).fetchone()
                if layer.get('ordinal') != count or row is None or int(row[0]) != count:
                    raise ValueError('origin frontier proof differs')
            if sum(item['ordinal'] for item in layers) != manifest['literal_requests']:
                raise ValueError('origin completed layers do not cover its literal set')
            validate_origin_membership(db, parent, layers)
            scope = {'origin_source': parent, 'origin_literal_requests': manifest['literal_requests'],
                     'origin_completed_layers': layers, 'ancestor_retained_certified': 0,
                     'edges': [], 'ancestor_history': {}}
        elif manifest.get('inherited_closed_scope') != scope:
            raise ValueError('parent inherited scope linkage differs')
        scope['ancestor_retained_certified'] += manifest['retained_certified']
        scope['edges'].append({'source': parent, 'target': target, 'manifest_sha256': hashlib.sha256(payload.encode()).hexdigest()})
        scope['ancestor_history'][parent] = manifest['source_history']
    count = db.execute('SELECT count(*) FROM async_jobs WHERE policy=?', (source,)).fetchone()[0]
    if scope['origin_literal_requests'] != scope['ancestor_retained_certified'] + count:
        raise ValueError('inherited closed scope has overlaps or omissions')
    return scope


def validate_coverage(db: sqlite3.Connection, source: str, target: str, manifest: dict[str, Any]) -> None:
    states = dict(db.execute('SELECT state,count(*) FROM async_jobs WHERE policy=? GROUP BY state', (source,)))
    if (manifest.get('schema') != 'zipper-async-migration-v1' or manifest.get('source') != source
            or manifest.get('target') != target or manifest.get('source_retired') is not True
            or manifest.get('source_states') != states or manifest.get('literal_requests') != sum(states.values())
            or manifest.get('retained_certified') != states.get('certified', 0)
            or manifest.get('queued_for_recomputation') != sum(v for k, v in states.items() if k != 'certified')):
        raise ValueError('migration manifest identity or source counts differ')
    counts = Counter(dict(db.execute('SELECT disposition,count(*) FROM async_migration_coverage WHERE source=? AND target=? GROUP BY disposition', (source, target))))
    if counts != Counter({'retained_source_certificate': manifest['retained_certified'], 'recompute_target': manifest['queued_for_recomputation']}):
        raise ValueError('migration coverage counts differ')
    missing = db.execute("SELECT count(*) FROM async_jobs s LEFT JOIN async_migration_coverage c ON c.source=s.policy AND c.target=? AND c.key=s.key WHERE s.policy=? AND (c.key IS NULL OR c.old_state!=s.state OR c.disposition!=CASE WHEN s.state='certified' THEN 'retained_source_certificate' ELSE 'recompute_target' END)", (target, source)).fetchone()[0]
    wrong = db.execute("SELECT count(*) FROM async_jobs s LEFT JOIN async_jobs t ON t.policy=? AND t.key=s.key WHERE s.policy=? AND ((s.state='certified' AND t.key IS NOT NULL) OR (s.state!='certified' AND (t.key IS NULL OR t.n!=s.n OR t.word!=s.word)))", (target, source)).fetchone()[0]
    total = db.execute('SELECT count(*) FROM async_jobs WHERE policy=?', (target,)).fetchone()[0]
    if missing or wrong or total != manifest['queued_for_recomputation']:
        raise ValueError('migration literal coverage has overlaps or omissions')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database', type=Path)
    parser.add_argument('--source-policy', required=True)
    parser.add_argument('--expected-literal-count', required=True, type=int)
    parser.add_argument('--expected-target', required=True, help='independently prepared current production policy key')
    parser.add_argument('--expected-parent-manifest-sha256', help='SHA256 of canonical externally archived immediate parent manifest')
    parser.add_argument('--expected-origin', help='externally recorded origin policy key')
    scope = parser.add_mutually_exclusive_group(required=True)
    scope.add_argument('--completed-layer', action='append', help='fully admitted rank:length layer')
    scope.add_argument('--inherit-closed-scope', action='store_true', help='authenticate the complete parent migration chain')
    args = parser.parse_args()
    if args.inherit_closed_scope and (not args.expected_parent_manifest_sha256 or not args.expected_origin):
        parser.error('inherited scope requires external parent manifest and origin anchors')
    old = read_policy(args.database, args.source_policy)
    with tempfile.TemporaryDirectory(prefix='zipper-policy-') as directory:
        path = Path(directory) / 'policy.sqlite3'
        key = initialize(path, Limits(**old['limits']), old['mode'], 256 * 2**20)
        policy = read_policy(path, key)
    layers = tuple(tuple(map(int, item.split(':'))) for item in (args.completed_layer or []))
    if any(len(item) != 2 for item in layers):
        parser.error('completed layers must be rank:length')
    print(json.dumps(migrate(args.database, args.source_policy, policy, expected_count=args.expected_literal_count, completed_layers=layers, inherit_closed_scope=args.inherit_closed_scope, expected_target=args.expected_target, expected_parent_manifest_sha256=args.expected_parent_manifest_sha256, expected_origin=args.expected_origin), sort_keys=True))


if __name__ == '__main__':
    main()
