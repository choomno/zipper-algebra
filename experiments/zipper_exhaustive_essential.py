"""Census of full-support literal signed words, stored independently.

Each completed input has an independently replayed zipper certificate.
Connected-block output separation permits reducing redundant ambient ranks and
disconnected support; this does not execute every ambient word or measure its
cost. Optional far-commutation representatives are an experimental input domain:
words in a class need not have identical zipper outputs. The legacy literal
census is never opened or modified by this runner.
"""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import FIRST_COMPLETED, wait
from contextlib import closing
from dataclasses import asdict
from datetime import datetime, timezone
from functools import lru_cache
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import time
import zlib

from experiments import zipper_exhaustive_collected as engine
from zipper_algebra.zipping_algorithm import essential_enumeration as enumeration
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import _strict_equal

SCHEMA = 'zipper-essential-representative-census-v1'
LITERAL_COVERAGE = ('connected-full-support-literal-words; '
                    'proved-connected-block-output-separation; '
                    'ambient-literal-executions-and-costs-not-covered')
TRACE_COVERAGE = ('experimental-connected-full-support-trace-representatives; '
                  'no-claim-of-identical-output-for-noncanonical-commutation-representatives; '
                  'literal-input-factorization-does-not-certify-ambient-execution')
DEFAULT_TARGETS = {1: 0, **{n: 10 for n in range(2, 9)}}
_worker_trace = False


@lru_cache(maxsize=16)
def _enumerator(n, length, commutations):
    return enumeration.EssentialWords(n, length, commute=commutations)


def _count(n, length, commutations):
    return _enumerator(n, length, commutations).count


def _word(n, length, ordinal, commutations):
    return list(_enumerator(n, length, commutations).word_at(ordinal))


def _first(n):
    return [0 if n == 1 else n - 1, 0]


def _advance(n, cursor, commutations):
    length, ordinal = cursor
    return [length + 1, 0] if ordinal + 1 == _count(n, length, commutations) else [length, ordinal + 1]


def identity(batch_size, commutations=False, trace=False, verifier_allocator=None):
    return {'schema': SCHEMA, 'coverage_scope': TRACE_COVERAGE if commutations else LITERAL_COVERAGE,
            'commutations': commutations, 'trace': trace,
            'runner_sha256': engine._hash(__file__),
            'enumeration_sha256': engine._hash(enumeration.__file__),
            'engine': engine.identity(batch_size, 'native', 'cpp', verifier_allocator)}


def _assert_source(expected):
    engine._assert_source(expected['engine'])
    if (engine._hash(__file__) != expected['runner_sha256']
            or engine._hash(enumeration.__file__) != expected['enumeration_sha256']):
        raise ValueError('essential enumeration source changed')


def _initialize(limits, expected, cache_bytes, cache_entries, verifier_cache_bytes):
    global _worker_trace
    _worker_trace = expected['trace']
    engine._initialize(limits, expected['engine'], cache_bytes, cache_entries, 0, 0,
                       verifier_cache_bytes, cache_entries)


def _certificate_identity(certificate, n, word, limits, expected):
    if not _strict_equal(certificate.get('trace_enabled'), expected['trace']):
        raise ValueError('certificate trace setting differs')
    # Reuse all existing strict build/input/catalog/limit checks. The legacy
    # campaign has trace=False fixed; only its trace assertion is adapted here.
    shape = dict(certificate, trace_enabled=False)
    engine._certificate_identity(shape, n, word, limits, expected['engine'])


def _failure(n, word, stage, error, encoded=None, producer_ns=0, checker_ns=0):
    return {'schema': 'zipper-essential-worker-failure-v1', 'n': n, 'word': word,
            'stage': stage, 'diagnostic': {'exception': type(error).__name__, 'message': str(error)},
            'certificate': encoded, 'producer_ns': producer_ns, 'checker_ns': checker_ns}


def _compute(n, word):
    encoded = None
    producer_ns = checker_ns = 0
    stage = 'producer'
    start = time.perf_counter_ns()
    try:
        certificate = engine._session.reduce(word, n=n, limits=engine._worker_limits,
            trace=_worker_trace, workers=1, batch_size=engine._worker_identity['batch_size'])
        producer_ns = time.perf_counter_ns() - start
        stage = 'encoding'
        encoded = engine._json(certificate)
        stage = 'verifier'
        start = time.perf_counter_ns()
        reloaded = engine._loads(encoded)
        if engine._verifier.check(reloaded) is not True:
            raise ValueError('representative certificate failed independent C++ replay')
        checker_ns = time.perf_counter_ns() - start
        return encoded, producer_ns, checker_ns
    except Exception as error:
        elapsed = time.perf_counter_ns() - start
        if stage == 'producer':
            producer_ns = elapsed
        elif stage == 'verifier':
            checker_ns = elapsed
        return _failure(n, word, stage, error, encoded, producer_ns, checker_ns)


def _schema(db):
    db.executescript('''
        CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS configurations(hash TEXT PRIMARY KEY,config TEXT NOT NULL);
        CREATE TABLE IF NOT EXISTS ranks(n INTEGER PRIMARY KEY,length INTEGER,ordinal INTEGER);
        CREATE TABLE IF NOT EXISTS representatives(
          n INTEGER,length INTEGER,ordinal INTEGER,word TEXT,input_sha TEXT,status TEXT,
          certificate BLOB,certificate_sha TEXT,certificate_bytes INTEGER,config_sha TEXT,
          producer_ns INTEGER,checker_ns INTEGER,metrics TEXT,production_check TEXT,
          PRIMARY KEY(n,length,ordinal));
        CREATE TABLE IF NOT EXISTS attempts(
          n INTEGER,length INTEGER,ordinal INTEGER,word TEXT,input_sha TEXT,status TEXT,
          certificate BLOB,certificate_sha TEXT,certificate_bytes INTEGER,config_sha TEXT,
          producer_ns INTEGER,checker_ns INTEGER,metrics TEXT,production_check TEXT,
          archived_utc TEXT);
        CREATE TABLE IF NOT EXISTS failures(
          id INTEGER PRIMARY KEY,n INTEGER,length INTEGER,ordinal INTEGER,config_sha TEXT,
          failure TEXT,failure_sha TEXT,certificate BLOB,certificate_sha TEXT,
          certificate_bytes INTEGER,recorded_utc TEXT);
    ''')


def _validate_row(row, configs, expected):
    (n, length, ordinal, word_json, input_sha, status, blob, digest, size,
     config_sha, producer, checker, metrics, check) = row
    if (type(n) is not int or n < 1 or type(length) is not int or length < 0
            or type(ordinal) is not int or not 0 <= ordinal < _count(n, length, expected['commutations'])):
        raise ValueError('invalid representative ordinal')
    word = _word(n, length, ordinal, expected['commutations'])
    if word_json != engine._json(word) or input_sha != engine._sha(engine._json({'n': n, 'word': word})):
        raise ValueError('representative input binding differs')
    config = configs.get(config_sha)
    if config is None or str(n) not in config['targets'] or length > config['targets'][str(n)]:
        raise ValueError('representative configuration binding differs')
    if check != expected['engine']['production_check']:
        raise ValueError('independent verification provenance differs')
    try:
        raw = zlib.decompress(blob)
        if len(raw) != size or hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError('certificate hash or length differs')
        certificate = engine._loads(raw)
    except (zlib.error, UnicodeError, json.JSONDecodeError, TypeError) as error:
        raise ValueError('invalid compressed certificate') from error
    _certificate_identity(certificate, n, word, config['limits'], expected)
    if (status != certificate['status'] or metrics != engine._json(certificate['metrics'])
            or type(producer) is not int or producer < 0 or type(checker) is not int or checker < 0):
        raise ValueError('representative columns differ from certificate')
    return certificate


def audit(db, expected, *, verifier=None):
    metadata = dict(db.execute('SELECT key,value FROM metadata'))
    if metadata.get('identity') != engine._json(expected):
        raise ValueError('essential campaign identity changed; use a separate campaign directory')
    configs = {}
    for digest, encoded in db.execute('SELECT hash,config FROM configurations'):
        config = engine._loads(encoded)
        if engine._sha(encoded) != digest or config['identity_sha256'] != engine._sha(engine._json(expected)):
            raise ValueError('configuration hash or identity differs')
        configs[digest] = config
    cursors = {n: _first(n) for n, _, _ in db.execute('SELECT * FROM ranks')}
    stopped = set()
    total = replayed = 0
    complete_by_layer = {}
    def replay(certificate):
        nonlocal replayed
        if verifier is not None:
            try:
                accepted = verifier.check(certificate)
            except Exception as error:
                raise ValueError('independent replay rejected saved certificate') from error
            if accepted is not True:
                raise ValueError('independent replay rejected saved certificate')
            replayed += 1
    for row in db.execute('SELECT * FROM representatives ORDER BY n,length,ordinal'):
        n, length, ordinal = row[:3]
        if n not in cursors or n in stopped or [length, ordinal] != cursors[n]:
            raise ValueError('representatives are not consecutive rank prefixes')
        certificate = _validate_row(row, configs, expected)
        replay(certificate)
        total += 1
        if certificate['status'] == 'complete':
            complete_by_layer[n, length] = complete_by_layer.get((n, length), 0) + 1
            cursors[n] = _advance(n, cursors[n], expected['commutations'])
        else:
            stopped.add(n)
    for n, length, ordinal in db.execute('SELECT * FROM ranks'):
        if cursors[n] != [length, ordinal]:
            raise ValueError('saved representative cursor differs')
    archived = 0
    for row in db.execute('SELECT * FROM attempts'):
        replay(_validate_row(row[:-1], configs, expected))
        if row[0] not in cursors or list(row[1:3]) > cursors[row[0]]:
            raise ValueError('archived attempt lies beyond its rank cursor')
        if row[5] == 'complete':
            raise ValueError('complete representative cannot be archived')
        archived += 1
    failures = 0
    for row in db.execute('SELECT * FROM failures'):
        _validate_failure(row, configs, expected, cursors)
        failures += 1
    return {'representatives': total, 'archived_attempts': archived, 'failures': failures,
            'complete_by_layer': [[n, length, count] for (n, length), count in sorted(complete_by_layer.items())],
            'hash_and_binding_checked': total + archived, 'replayed': replayed}


def _validate_failure(row, configs, expected, cursors):
    _, n, length, ordinal, config_sha, encoded, digest, blob, cert_sha, size, _ = row
    config = configs.get(config_sha)
    if (n not in cursors or [length, ordinal] > cursors[n]
            or type(ordinal) is not int or not 0 <= ordinal < _count(n, length, expected['commutations'])
            or config is None or str(n) not in config['targets'] or length > config['targets'][str(n)]):
        raise ValueError('failure input/configuration/cursor binding differs')
    if engine._sha(encoded) != digest:
        raise ValueError('failure diagnostic hash differs')
    failure = engine._loads(encoded)
    if (not isinstance(failure, dict) or failure.get('schema') != 'zipper-essential-worker-failure-v1'
            or not _strict_equal(failure.get('n'), n)
            or not _strict_equal(failure.get('word'), _word(n, length, ordinal, expected['commutations']))
            or failure.get('stage') not in ('producer', 'encoding', 'verifier', 'worker')
            or not isinstance(failure.get('diagnostic'), dict)
            or any(type(failure.get(key)) is not int or failure[key] < 0 for key in ('producer_ns', 'checker_ns'))):
        raise ValueError('invalid worker failure record')
    if blob is None:
        if cert_sha is not None or size is not None:
            raise ValueError('missing failure certificate has metadata')
    else:
        try:
            raw = zlib.decompress(blob)
        except (zlib.error, TypeError) as error:
            raise ValueError('invalid failure certificate compression') from error
        if hashlib.sha256(raw).hexdigest() != cert_sha or len(raw) != size:
            raise ValueError('failure certificate hash or length differs')
    return failure


def _commit_failure(db, n, cursor, failure, config_sha, config, expected):
    payload = dict(failure)
    certificate = payload.pop('certificate', None)
    raw = None if certificate is None else certificate.encode()
    encoded = engine._json(payload)
    row = (None, n, *cursor, config_sha, encoded, engine._sha(encoded),
           None if raw is None else zlib.compress(raw),
           None if raw is None else hashlib.sha256(raw).hexdigest(),
           None if raw is None else len(raw), datetime.now(timezone.utc).isoformat())
    _validate_failure(row, {config_sha: config}, expected, {n: cursor})
    with db:
        result = db.execute('INSERT INTO failures VALUES(' + ','.join('?' * 11) + ')', row)
    return result.lastrowid


def _commit(db, n, cursor, word, result, config_sha, config, expected):
    encoded, producer, checker = result
    certificate = engine._loads(encoded)
    _certificate_identity(certificate, n, word, config['limits'], expected)
    raw = encoded.encode()
    row = (n, *cursor, engine._json(word), engine._sha(engine._json({'n': n, 'word': word})),
           certificate['status'], zlib.compress(raw), engine._sha(encoded), len(raw), config_sha,
           producer, checker, engine._json(certificate['metrics']), expected['engine']['production_check'])
    _validate_row(row, {config_sha: config}, expected)
    following = _advance(n, cursor, expected['commutations']) if certificate['status'] == 'complete' else cursor
    with db:
        previous = db.execute('SELECT * FROM representatives WHERE n=? AND length=? AND ordinal=?', (n, *cursor)).fetchone()
        if previous:
            if previous[5] == 'complete':
                raise ValueError('completed representatives are immutable')
            db.execute('INSERT INTO attempts VALUES(' + ','.join('?' * 15) + ')',
                       (*previous, datetime.now(timezone.utc).isoformat()))
        db.execute('INSERT OR REPLACE INTO representatives VALUES(' + ','.join('?' * 14) + ')', row)
        db.execute('UPDATE ranks SET length=?,ordinal=? WHERE n=?', (*following, n))
    return certificate['status'], following


def exhaust(output, *, targets=None, workers=8, batch_size=64, wall_seconds=3600,
            limits=None, max_new_words=None, commutations=False, trace=False,
            cache_bytes=256 * 1024**2, cache_entries=262144,
            verifier_cache_bytes=128 * 1024**2, verifier_allocator=None, retry_failures=False,
            _executor_factory=engine._process_pool):
    targets = DEFAULT_TARGETS.copy() if targets is None else dict(targets)
    if (not targets or any(type(n) is not int or n < 1 or type(length) is not int or length < 0
            or (n == 1 and length != 0) for n, length in targets.items())):
        raise ValueError('positive integer ranks and nonnegative lengths required; rank one has only length zero')
    if any(type(v) is not int or v < 1 for v in (workers, batch_size)):
        raise ValueError('positive workers and batch_size required')
    if any(type(v) is not int or v < 0 for v in (cache_bytes, cache_entries, verifier_cache_bytes)):
        raise ValueError('nonnegative integer cache budgets required')
    if type(wall_seconds) not in (int, float) or not math.isfinite(wall_seconds) or wall_seconds <= 0:
        raise ValueError('positive finite wall_seconds required')
    if max_new_words is not None and (type(max_new_words) is not int or max_new_words < 0):
        raise ValueError('nonnegative integer max_new_words required')
    if type(commutations) is not bool or type(trace) is not bool or type(retry_failures) is not bool:
        raise ValueError('commutations, trace and retry_failures must be booleans')
    limits = Limits(max_steps=1000000, max_terms=100000, expansion_budget=1000000) if limits is None else limits
    if not isinstance(limits, Limits):
        raise ValueError('limits must be a Limits record')
    expected = identity(batch_size, commutations, trace, verifier_allocator)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    lock = (output / '.essential.lock').open('a+')
    try:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        lock.close()
        raise ValueError('another writer owns this essential campaign') from error
    try:
        db = sqlite3.connect(output / 'essential.sqlite3')
    except BaseException:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()
        raise
    pool = None
    queues = {n: deque() for n in targets}
    discarded = set()
    drain = [False]
    previous_handler = engine._install_drain_handler(drain)
    try:
        db.execute('PRAGMA journal_mode=WAL')
        db.execute('PRAGMA synchronous=NORMAL')
        _schema(db)
        saved = db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()
        if saved is None:
            if db.execute('SELECT count(*) FROM representatives').fetchone()[0]:
                raise ValueError('representatives lack campaign identity')
            with db:
                db.execute("INSERT INTO metadata VALUES('identity',?)", (engine._json(expected),))
        audited = audit(db, expected)
        config = {'targets': {str(n): length for n, length in targets.items()}, 'workers': workers,
                  'limits': asdict(limits), 'identity_sha256': engine._sha(engine._json(expected)),
                  'wall_seconds': wall_seconds, 'max_new_words': max_new_words, 'retry_failures': retry_failures,
                  'cache_bytes': cache_bytes, 'cache_entries': cache_entries,
                  'verifier_cache_bytes': verifier_cache_bytes,
                  'created_utc': datetime.now(timezone.utc).isoformat()}
        config_sha = engine._sha(engine._json(config))
        with db:
            db.execute('INSERT INTO configurations VALUES(?,?)', (config_sha, engine._json(config)))
            for n in targets:
                db.execute('INSERT OR IGNORE INTO ranks VALUES(?,?,?)', (n, *_first(n)))
        cursors = {n: [length, ordinal] for n, length, ordinal in db.execute('SELECT * FROM ranks') if n in targets}
        next_submit = {n: cursor[:] for n, cursor in cursors.items()}
        active = {n for n, cursor in cursors.items() if cursor[0] <= targets[n]}
        reasons = {n: 'target_length' for n in targets if n not in active}
        for n in list(active):
            if not retry_failures and db.execute('SELECT 1 FROM failures WHERE n=? AND length=? AND ordinal=? LIMIT 1', (n, *cursors[n])).fetchone():
                active.remove(n)
                reasons[n] = 'worker_failure_needs_explicit_retry'
                continue
            prior = db.execute('SELECT certificate,status FROM representatives WHERE n=? AND length=? AND ordinal=?', (n, *cursors[n])).fetchone()
            if prior:
                old_limits = engine._loads(zlib.decompress(prior[0]))['limits']
                if prior[1] == 'control_obligation' or old_limits == asdict(limits) or any(asdict(limits)[key] < value for key, value in old_limits.items()):
                    active.remove(n)
                    reasons[n] = 'control_obligation' if prior[1] == 'control_obligation' else 'word_limit_needs_larger_limits'
        target_counts = {n: sum(_count(n, length, commutations)
                                for length in range(targets[n] + 1)) for n in targets}
        complete_counts = {n: sum(count for rank, length, count in audited['complete_by_layer']
                                  if rank == n and length <= targets[n]) for n in targets}
        started = time.perf_counter()
        new = peak = failed = 0
        failure_records = []
        rotation = deque(sorted(targets))
        admitting = bool(active) and max_new_words != 0
        stop = 'running' if admitting else ('word_budget' if max_new_words == 0 else 'rank_limits')
        peak_rank = {n: 0 for n in targets}
        def snapshot(state):
            data = {'schema': SCHEMA, 'coverage_scope': expected['coverage_scope'], 'identity': expected,
                    'database': str((output / 'essential.sqlite3').resolve()), 'status': state,
                    'stop_reason': stop, 'audit': audited, 'config_sha256': config_sha,
                    'new_or_resumed_representatives': new, 'new_failures': failed,
                    'failure_records': failure_records, 'new_work_seconds': time.perf_counter() - started,
                    'peak_in_flight': peak, 'peak_per_rank_in_flight': peak_rank,
                    'ranks': [{'n': n, 'next_cursor': cursors[n], 'stop_reason': reasons.get(n, stop),
                               'target_representatives': target_counts[n],
                               'complete_representatives': complete_counts[n]}
                              for n in sorted(targets)]}
            with db:
                db.execute("INSERT OR REPLACE INTO metadata VALUES('run_state',?)", (engine._json({'status': state, 'pid': os.getpid(), 'recorded_utc': datetime.now(timezone.utc).isoformat(), 'stop_reason': stop}),))
            temporary = output / 'essential-summary.json.tmp'
            temporary.write_text(engine._json(data) + '\n')
            temporary.replace(output / 'essential-summary.json')
            return data
        snapshot('running')
        if admitting:
            pool = _executor_factory(max_workers=workers, initializer=_initialize,
                initargs=(asdict(limits), expected, cache_bytes, cache_entries, verifier_cache_bytes))
        last_snapshot = started
        while (active and admitting) or any(queues.values()) or discarded:
            discarded = {future for future in discarded if not future.done()}
            if admitting and (drain[0] or time.perf_counter() - started >= wall_seconds):
                admitting = False
                stop = 'signal_drain' if drain[0] else 'wall_budget'
            occupancy = sum(map(len, queues.values())) + len(discarded)
            per_rank_cap = max(1, math.ceil(2 * workers / max(1, len(active))))
            if admitting:
                for _ in range(2 * workers):
                    if occupancy >= 2 * workers or (max_new_words is not None and new + occupancy >= max_new_words):
                        break
                    chosen = None
                    for _ in range(len(rotation)):
                        n = rotation[0]
                        rotation.rotate(-1)
                        if n in active and next_submit[n][0] <= targets[n] and len(queues[n]) < per_rank_cap:
                            chosen = n
                            break
                    if chosen is None:
                        break
                    n = chosen
                    cursor = next_submit[n][:]
                    word = _word(n, *cursor, commutations)
                    queues[n].append((cursor, word, pool.submit(_compute, n, word)))
                    next_submit[n] = _advance(n, cursor, commutations)
                    occupancy += 1
                    peak = max(peak, occupancy)
                    peak_rank[n] = max(peak_rank[n], len(queues[n]))
            ready = [n for n in rotation if queues[n] and queues[n][0][2].done()]
            if not ready:
                futures = [job[2] for queue in queues.values() for job in queue if not job[2].done()] + list(discarded)
                if not futures:
                    if max_new_words is not None and new >= max_new_words:
                        stop = 'word_budget'
                    break
                wait(futures, timeout=.2, return_when=FIRST_COMPLETED)
            else:
                _assert_source(expected)
                for n in ready:
                    cursor, word, future = queues[n].popleft()
                    if cursor != cursors[n]:
                        raise ValueError('worker is not the next representative')
                    try:
                        result = future.result()
                    except Exception as error:
                        result = _failure(n, word, 'worker', error)
                    if isinstance(result, dict):
                        failure_id = _commit_failure(db, n, cursor, result, config_sha, config, expected)
                        failure_records.append({'id': failure_id, 'n': n, 'length': cursor[0],
                                                'ordinal': cursor[1], 'word': word, 'stage': result['stage'],
                                                'diagnostic': result['diagnostic']})
                        status, following = 'worker_failure', cursor
                        failed += 1
                    else:
                        status, following = _commit(db, n, cursor, word, result, config_sha, config, expected)
                        new += 1
                        if status == 'complete':
                            complete_counts[n] += 1
                    cursors[n] = following
                    if status != 'complete' or following[0] > targets[n]:
                        active.discard(n)
                        reasons[n] = 'target_length' if status == 'complete' else status
                        for _, _, later in queues[n]:
                            later.cancel()
                            if not later.done():
                                discarded.add(later)
                        queues[n].clear()
            if time.perf_counter() - last_snapshot >= 5:
                snapshot('running' if admitting else 'draining')
                last_snapshot = time.perf_counter()
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
            pool = None
        _assert_source(expected)
        complete = all(cursors[n][0] > targets[n] for n in targets)
        if complete:
            stop = 'targets_completed'
        elif stop == 'running':
            stop = 'rank_limits'
        return snapshot('completed' if complete else 'paused')
    except BaseException:
        if 'snapshot' in locals():
            try:
                db.rollback()
                stop = 'failed'
                snapshot('failed')
            except Exception:
                pass
        raise
    finally:
        try:
            for queue in queues.values():
                for _, _, future in queue:
                    future.cancel()
            if pool is not None:
                pool.shutdown(wait=True, cancel_futures=True)
        finally:
            db.close()
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
            lock.close()
            engine._restore_drain_handler(previous_handler)



def audit_store(output, *, replay=False):
    """Validate source-bound saved prefixes without opening a writable connection."""
    path = (Path(output) / 'essential.sqlite3').resolve()
    with closing(sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)) as db:
        # A deferred transaction pins one WAL snapshot at the first SELECT and
        # keeps it through the complete hash/binding audit and optional replay.
        db.execute('BEGIN')
        saved = db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()
        if saved is None:
            raise ValueError('essential campaign identity is missing')
        expected = engine._loads(saved[0])
        _assert_source(expected)
        if not replay:
            checked = audit(db, expected)
        else:
            from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier
            session = engine._new_session(expected['engine'], cache_bytes=0, cache_entries=0)
            try:
                allocator = expected['engine'].get('verifier_allocator')
                with CppVerifier(producer_metadata=session.metadata,
                        allocator_library=None if allocator is None else Path(allocator['path'])) as verifier:
                    if (engine._verifier_build(verifier.metadata) != expected['engine']['verifier_build']
                            or verifier.metadata.get('allocator_library') != engine._allocator_runtime(allocator)):
                        raise ValueError('replay verifier build or allocator differs')
                    checked = audit(db, expected, verifier=verifier)
            finally:
                session.close()
        return {'schema': SCHEMA, 'coverage_scope': expected['coverage_scope'],
                'database': str(path), 'audit': checked}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--targets', type=json.loads, default=DEFAULT_TARGETS)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--wall-seconds', type=float, default=3600)
    parser.add_argument('--batch-size', type=int, default=64)
    parser.add_argument('--cache-mib', type=int, default=256)
    parser.add_argument('--verifier-cache-mib', type=int, default=128)
    parser.add_argument('--cache-entries', type=int, default=262144)
    parser.add_argument('--max-new-words', type=int)
    parser.add_argument('--max-steps', type=int, default=1000000)
    parser.add_argument('--max-terms', type=int, default=100000)
    parser.add_argument('--expansion-budget', type=int, default=1000000)
    parser.add_argument('--trace', action='store_true')
    parser.add_argument('--audit-only', action='store_true')
    parser.add_argument('--replay-audit', action='store_true')
    parser.add_argument('--retry-failures', action='store_true')
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--commutation-classes', dest='commutations', action='store_true',
                       help='Experimental representatives; outputs need not agree within a class')
    modes.add_argument('--retain-commuting-orders', dest='commutations', action='store_false',
                       help='Keep every literal order within connected support (default)')
    parser.set_defaults(commutations=False)
    parser.add_argument('--verifier-allocator', type=Path)
    args = parser.parse_args()
    if args.audit_only:
        print(engine._json(audit_store(args.output, replay=args.replay_audit)))
        return
    if args.replay_audit:
        parser.error('--replay-audit requires --audit-only')
    result = exhaust(args.output, targets={int(n): length for n, length in args.targets.items()},
                     workers=args.workers, wall_seconds=args.wall_seconds, batch_size=args.batch_size,
                     cache_bytes=args.cache_mib * 1024**2, cache_entries=args.cache_entries,
                     verifier_cache_bytes=args.verifier_cache_mib * 1024**2,
                     max_new_words=args.max_new_words, trace=args.trace,
                     commutations=args.commutations,
                     verifier_allocator=args.verifier_allocator, retry_failures=args.retry_failures,
                     limits=Limits(max_steps=args.max_steps, max_terms=args.max_terms,
                                   expansion_budget=args.expansion_budget))
    print(engine._json(result))


if __name__ == '__main__':
    main()
