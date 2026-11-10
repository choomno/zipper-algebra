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
                digest.update(len(data).to_bytes(8, 'big'))
                digest.update(data)
            count += 1
        result[table] = {'rows': count, 'sha256': digest.hexdigest()}
    return result


def migrate(path: Path, source: str, policy: dict[str, Any]) -> dict[str, Any]:
    target = _sha(policy)
    if source == target or policy.get('schema') != SCHEMA or policy.get('package_sources') != reconstruction_sources():
        raise ValueError('a distinct policy bound to current package sources is required')
    with closing(sqlite3.connect(path, timeout=60, isolation_level=None)) as db:
        db.execute('PRAGMA foreign_keys=ON')
        db.execute('BEGIN IMMEDIATE')
        try:
            stored = db.execute('SELECT payload FROM async_policies WHERE key=?', (source,)).fetchone()
            if stored is None or _sha(json.loads(stored[0])) != source:
                raise ValueError('source policy missing or digest mismatch')
            old_policy = json.loads(stored[0])
            if policy['limits'] != old_policy['limits'] or policy['mode'] != old_policy['mode']:
                raise ValueError('epoch switch must preserve limits and execution mode')
            db.execute('CREATE TABLE IF NOT EXISTS async_migrations(source TEXT NOT NULL REFERENCES async_policies(key), target TEXT NOT NULL UNIQUE REFERENCES async_policies(key), manifest TEXT NOT NULL, PRIMARY KEY(source,target))')
            db.execute('CREATE TABLE IF NOT EXISTS async_migration_coverage(source TEXT NOT NULL,target TEXT NOT NULL,key TEXT NOT NULL,old_state TEXT NOT NULL,disposition TEXT NOT NULL,PRIMARY KEY(source,target,key),FOREIGN KEY(source,target) REFERENCES async_migrations(source,target),FOREIGN KEY(source,key) REFERENCES async_jobs(policy,key))')
            existing = db.execute('SELECT manifest FROM async_migrations WHERE source=? AND target=?', (source, target)).fetchone()
            before = history_digest(db, source)
            if existing:
                manifest = json.loads(existing[0])
                if manifest['source_history'] != before:
                    raise ValueError('source epoch changed after migration')
                validate_coverage(db, source, target, manifest)
                db.execute('COMMIT')
                return manifest
            if db.execute('SELECT 1 FROM async_policies WHERE key=?', (target,)).fetchone():
                raise ValueError('target policy already exists without migration linkage')
            now = time.time()
            counts = dict(db.execute('SELECT state,count(*) FROM async_jobs WHERE policy=? GROUP BY state', (source,)))
            manifest = {'schema': 'zipper-async-migration-v1', 'source': source, 'target': target,
                        'created': now, 'source_history': before, 'source_states': counts,
                        'source_retired': True, 'source_guard': 'SQLite triggers reject all old-policy row mutations',
                        'literal_requests': sum(counts.values()), 'retained_certified': counts.get('certified', 0),
                        'queued_for_recomputation': sum(v for k, v in counts.items() if k != 'certified'),
                        'coverage': 'Each source literal is retained as source certified evidence or queued once under target; provisional candidates are recomputed.'}
            db.execute('INSERT INTO async_policies VALUES (?,?)', (target, _json(policy)))
            db.execute('INSERT INTO async_migrations VALUES (?,?,?)', (source, target, _json(manifest)))
            db.execute("INSERT INTO async_migration_coverage SELECT policy,?,key,state,CASE WHEN state='certified' THEN 'retained_source_certificate' ELSE 'recompute_target' END FROM async_jobs WHERE policy=?", (target, source))
            db.execute("INSERT INTO async_jobs SELECT ?,key,n,word,'queued',NULL,NULL,?,? FROM async_jobs WHERE policy=? AND state!='certified'", (target, now, now, source))
            validate_coverage(db, source, target, manifest)
            if before != history_digest(db, source):
                raise ValueError('source history changed during migration')
            retire_source(db, source)
            db.execute('COMMIT')
            return manifest
        except BaseException:
            db.execute('ROLLBACK')
            raise


def retire_source(db: sqlite3.Connection, source: str) -> None:
    """Reject accidental restart or modification of any retired policy row."""
    for table in TABLES:
        column = 'key' if table == 'async_policies' else 'policy'
        for operation in ('INSERT', 'UPDATE', 'DELETE'):
            refs = ('OLD', 'NEW') if operation == 'UPDATE' else (('OLD',) if operation == 'DELETE' else ('NEW',))
            condition = ' OR '.join(f"{ref}.{column}='{source}'" for ref in refs)
            name = f'retired_{source}_{table}_{operation.lower()}'
            db.execute(f"CREATE TRIGGER {name} BEFORE {operation} ON {table} WHEN {condition} BEGIN SELECT RAISE(ABORT, 'retired census source policy is immutable'); END")


def validate_coverage(db: sqlite3.Connection, source: str, target: str, manifest: dict[str, Any]) -> None:
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
    args = parser.parse_args()
    old = read_policy(args.database, args.source_policy)
    with tempfile.TemporaryDirectory(prefix='zipper-policy-') as directory:
        path = Path(directory) / 'policy.sqlite3'
        key = initialize(path, Limits(**old['limits']), old['mode'], 256 * 2**20)
        policy = read_policy(path, key)
    print(json.dumps(migrate(args.database, args.source_policy, policy), sort_keys=True))


if __name__ == '__main__':
    main()
