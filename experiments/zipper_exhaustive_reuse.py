"""Fair current-procedure census with SQLite shards and verified recursive reuse.

The original literal census remains a separate entry point. Each worker owns a
persistent producer, independent verifier and sidecar. A coordinator ledger
stores coverage and resumable enumeration cursors; no legacy database is opened.
"""
from __future__ import annotations

import argparse
import atexit
from collections import OrderedDict, deque
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
from concurrent.futures.process import BrokenProcessPool
from contextlib import ExitStack, closing
from dataclasses import asdict, replace
import fcntl
from functools import lru_cache
import hashlib
import json
import multiprocessing
from pathlib import Path
import sqlite3
import zlib
import uuid
from time import perf_counter
from typing import Any

from zipper_algebra.zipping_algorithm.census_reuse import (
    ReuseCensus, reconstruction_sources, verifier_binding, _canonical_digest, _IMPORT_SOURCES,
)
from zipper_algebra.zipping_algorithm.input_normal_form import input_normal_form, leaf_requests
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native import NativeSession, _decode_coefficient_object
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import completed_output, execution_cache_key, execution_identity, _identity
from zipper_algebra.zipping_algorithm.peeling_enumeration import NonPeelableWords

SCHEMA = 'zipper-current-reuse-campaign-v1'
RUNNER_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
PACKAGE_SOURCE_SHA256 = dict(_IMPORT_SOURCES)
_store: ReuseCensus | None = None
_sessions: tuple[Any, ...] = ()
_work: dict[str, Any] = {}


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _key(value):
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _assert_runner_source():
    if hashlib.sha256(Path(__file__).read_bytes()).hexdigest()!=RUNNER_SHA256:
        raise ValueError('runner source changed since import')


def _assert_package_sources():
    if reconstruction_sources()!=PACKAGE_SOURCE_SHA256:
        raise ValueError('reconstruction package sources changed since import')


def _close():
    global _store, _sessions
    store, _store = _store, None
    sessions, _sessions = _sessions, ()
    with ExitStack() as cleanup:
        for session in sessions:
            cleanup.callback(session.close)
        if store is not None:
            cleanup.callback(store.close)


def _initialize(path, settings):
    global _store, _sessions, _work
    cleanup = ExitStack()
    try:
        producer = NativeSession(cache_bytes=settings['cache_bytes'], cache_entries=settings['cache_entries'])
        cleanup.callback(producer.close)
        verifier = CppVerifier(producer_metadata=producer.metadata,
                               cache_bytes=settings['cache_bytes'], cache_entries=settings['cache_entries'],
                               allocator_library=None if settings['allocator'] is None else Path(settings['allocator']))
        cleanup.callback(verifier.close)
        _initialize_store(path, settings, producer, verifier)
    except BaseException:
        cleanup.close()
        raise
    cleanup.pop_all()
    atexit.register(_close)


def _initialize_store(path, settings, producer, verifier):
    global _store, _sessions, _work
    limits = Limits(**settings['limits'])
    if settings['mode'] == 'reuse':
        reduce = producer.reduce_reused
        check = verifier.check_reused
    else:
        reduce = producer.reduce
        check = verifier.check
    # Bootstrap only the source/execution identity with a completed rank-one
    # identity. Actual requests and the persisted policy retain the exact
    # requested limits, including budgets too small to finish this seed.
    seed_limits = replace(limits, max_steps=max(2, limits.max_steps),
                          max_terms=max(2, limits.max_terms))
    seed = reduce([], n=1, limits=seed_limits)
    if check(seed) is not True:
        raise ValueError('independent verifier rejected initialization result')
    _work = {side: {'calls': 0, 'dag_node_expansions': 0, 'evaluated_rule_applications': {}}
             for side in ('producer', 'verifier')}
    def record(side, value):
        state = _work[side]
        state['calls'] += 1
        metrics = value.get('metrics', {})
        state['dag_node_expansions'] += metrics.get('transition_evaluations', metrics.get('expanded_nodes', 0))
        for rule,count in metrics.get('evaluated_rule_applications', metrics.get('rule_applications', {})).items():
            state['evaluated_rule_applications'][rule] = state['evaluated_rule_applications'].get(rule, 0) + count
    def produce(n, word):
        result = reduce(word, n=n, limits=limits)
        record('producer', result)
        return result
    def verify(value):
        previous = verifier.last_result
        try:
            return check(value)
        finally:
            if verifier.last_result is not previous:
                record('verifier', verifier.last_result)
    _store = ReuseCensus(path, produce=produce,
                         verify=verify, producer_identity=execution_identity(seed),
                         verifier_identity=verifier_binding(verifier.metadata), limits=asdict(limits))
    _sessions = producer, verifier


def _ready():
    assert _store is not None
    _assert_runner_source()
    _assert_package_sources()
    return {'policy': _store.policy_key, 'identity': _store.policy,
            'runner_sha256':RUNNER_SHA256}


def _evaluate(n, word):
    assert _store is not None
    result = _store.evaluate(word, n=n)
    result.pop('output', None)  # Full exact output is materializable from the shard.
    if _sessions:
        result['worker_telemetry'] = {
            'scope': 'cumulative-actual-work-in-this-worker-process; excludes-initialization',
            'dag_metrics_scope': 'reuse leaf backend only; standalone transition work is in cache counters',
            'work': _work,
            'producer_cache': _sessions[0].cache_info(),
            'verifier_cache': _sessions[1].cache_info()}
    return result


def _snapshot_counts():
    assert _store is not None
    return _store.counts()


def _probe_ledger(path):
    """Identify an existing file before any writable connection or recovery."""
    if path.exists():
        with closing(sqlite3.connect(path.resolve().as_uri()+'?mode=ro',uri=True)) as probe:
            names={row[0] for row in probe.execute("SELECT name FROM sqlite_master WHERE type='table'")}
            if names and 'reuse_campaigns' not in names:
                raise ValueError('coordinator requires its own reuse campaign database')


def _pending_jobs(db, campaign, retry_incomplete):
    """Bounded keyset scan of retries present at the start of this invocation."""
    ceiling=db.execute('SELECT coalesce(max(rowid),0) FROM reuse_jobs WHERE campaign=?',(campaign,)).fetchone()[0]
    after=0
    statuses="('interrupted','incomplete','failed')" if retry_incomplete else "('interrupted')"
    while after<ceiling:
        rows=db.execute('SELECT rowid,key,n,word FROM reuse_jobs WHERE campaign=? '
                        'AND rowid>? AND rowid<=? AND status IN '+statuses+' ORDER BY rowid LIMIT 256',
                        (campaign,after,ceiling)).fetchall()
        if not rows:
            break
        for rowid,key,n,word in rows:
            after=rowid
            yield key,n,word


@lru_cache(maxsize=32)
def _layer(n, length):
    return NonPeelableWords(n, length)


def _small_layer(n, length):
    if n == 1:
        return [()] if length == 0 else []
    words = [(1,) * length, (-1,) * length] if length else []
    return [word for word in words if leaf_requests(input_normal_form(word, n=2))]


def _next(n, cursor, target):
    """Direct ordinal enumeration; rank two uses proved exponent representatives."""
    length, ordinal = cursor
    while length <= target:
        count = len(_small_layer(n, length)) if n <= 2 else _layer(n, length).count
        if ordinal < count:
            word = (_small_layer(n, length)[ordinal] if n <= 2
                    else _layer(n, length).word_at(ordinal))
            return list(word), [length, ordinal + 1], [length, ordinal]
        length, ordinal = length + 1, 0
    return None


def _ledger(db):
    names = {r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if names and 'reuse_campaigns' not in names:
        raise ValueError('coordinator requires its own reuse campaign database')
    db.execute('PRAGMA journal_mode=WAL')
    db.executescript('''
      CREATE TABLE IF NOT EXISTS reuse_campaigns(key TEXT PRIMARY KEY, identity TEXT NOT NULL);
      CREATE TABLE IF NOT EXISTS reuse_cursors(campaign TEXT NOT NULL, n INTEGER NOT NULL,
          length INTEGER NOT NULL, ordinal TEXT NOT NULL, PRIMARY KEY(campaign,n));
      CREATE TABLE IF NOT EXISTS reuse_worker_telemetry(campaign TEXT NOT NULL, run TEXT NOT NULL,
          worker INTEGER NOT NULL, payload TEXT NOT NULL, PRIMARY KEY(campaign,run,worker));
      CREATE TABLE IF NOT EXISTS reuse_round_robin(campaign TEXT PRIMARY KEY, last_rank INTEGER NOT NULL);
      CREATE TABLE IF NOT EXISTS reuse_jobs(campaign TEXT NOT NULL, key TEXT NOT NULL,
          n INTEGER NOT NULL, word TEXT NOT NULL, status TEXT NOT NULL,
          shard INTEGER, result TEXT, PRIMARY KEY(campaign,key));
      CREATE TABLE IF NOT EXISTS reuse_job_attempts(id INTEGER PRIMARY KEY,
          campaign TEXT NOT NULL, job TEXT NOT NULL, status TEXT NOT NULL, result TEXT NOT NULL);
    ''')


def _audit_leaf(row, request, producer, limits):
    """Check integrity, binding, and structure of one complete stored leaf exactly.

    JSON parsing owns these fresh dictionaries, so coefficient decoding can
    happen in its object hook. The terminal sum is canonicalized once; the
    digest then uses precisely the same bytes as output_fingerprint on that
    validated canonical sum. This remains an integrity audit, not rewrite replay.
    """
    if row[0] != _json(request):
        raise ValueError('verified leaf request mismatch')
    raw=zlib.decompress(row[1])
    if hashlib.sha256(raw).hexdigest()!=row[2]:
        raise ValueError('verified leaf content or binding mismatch')
    value=json.loads(raw, object_hook=_decode_coefficient_object)
    if type(value) is not dict:
        raise ValueError('verified leaf content or binding mismatch')
    metadata=_identity({'schema':value.get('schema'),
                        'source_sha256':value.get('source_sha256'),
                        'catalog_sha256':value.get('catalog_sha256'),
                        'execution':value.get('execution',{})})
    if (_json(metadata)!=_json(producer)
            or _json(value.get('n'))!=_json(request['n'])
            or _json(value.get('input_word'))!=_json(request['word'])
            or _json(value.get('limits'))!=_json(limits)):
        raise ValueError('verified leaf content or binding mismatch')
    normal=completed_output(value)
    if normal['completion']!='complete' or _canonical_digest(normal)!=row[3]:
        raise ValueError('verified leaf content or binding mismatch')


def audit_store(output):
    """Read-only integrity and dependency audit; does not claim fresh replay.

    One coordinator transaction fixes its campaigns/cursors/jobs snapshot. Each
    later shard transaction contains every completed root referenced there:
    the runner commits immutable shard rows before marking a job complete.
    This relies on that append-only contract, not arbitrary external mutations
    being globally atomic across separate databases. No writer lock is needed.
    """
    output = Path(output)
    path = output/'census.sqlite3'
    if not path.exists():
        return {'jobs':0, 'complete':0}
    _assert_runner_source()
    _assert_package_sources()
    db=sqlite3.connect(path.resolve().as_uri()+'?mode=ro',uri=True)
    shards={}
    counts={'jobs':0,'complete':0}
    checked=OrderedDict()
    try:
        db.execute('BEGIN')
        identities={}
        for key,payload in db.execute('SELECT key,identity FROM reuse_campaigns'):
            identity=json.loads(payload)
            if _key(identity)!=key:
                raise ValueError('campaign identity digest mismatch')
            worker_identity=identity.get('worker',{})
            if (type(worker_identity) is not dict or 'identity' not in worker_identity
                    or _key(worker_identity['identity'])!=worker_identity.get('policy')):
                raise ValueError('embedded worker policy identity digest mismatch')
            identities[key]=identity
        cursors={}
        for campaign,n,length,ordinal in db.execute('SELECT * FROM reuse_cursors'):
            if (campaign not in identities or type(n) is not int or n<1
                    or type(length) is not int or length<0
                    or type(ordinal) is not str or not ordinal.isascii() or not ordinal.isdigit()
                    or str(int(ordinal))!=ordinal):
                raise ValueError('invalid persisted census cursor')
            cursors[campaign,n]=(length,int(ordinal))
        admitted={}
        for campaign,key,n,encoded,status,worker,payload in db.execute('SELECT * FROM reuse_jobs'):
            word=json.loads(encoded)
            form=input_normal_form(word,n=n)
            if (_key({'n':n,'word':word})!=key or campaign not in identities
                    or status not in ('running','interrupted','complete','incomplete','failed')):
                raise ValueError('invalid persisted job identity or status')
            if identities[campaign]['domain']=='current-primitive-baseline-v1':
                cursor=cursors.get((campaign,n))
                length=len(word)
                ordinal=(_small_layer(n,length).index(tuple(word)) if n<=2
                         else _layer(n,length).ordinal(word))
                if cursor is None or (length,ordinal)>=cursor:
                    raise ValueError('admitted word lies outside its enumeration cursor')
                count,maximum=admitted.get((campaign,n),(0,0))
                admitted[campaign,n]=(count+1,max(maximum,length))
            counts['jobs']+=1
            if status!='complete':
                continue
            result=json.loads(payload)
            identity=identities[campaign]['worker']
            if (result.get('status')!='complete' or result.get('root')!=key
                    or result.get('policy')!=identity['policy']
                    or type(worker) is not int or worker<0):
                raise ValueError('completed job result binding mismatch')
            if worker not in shards:
                shard_path=output/'shards'/f'worker-{worker}.sqlite3'
                shards[worker]=sqlite3.connect(shard_path.resolve().as_uri()+'?mode=ro',uri=True)
                shards[worker].execute('BEGIN')
            shard=shards[worker]
            policy=identity['policy']
            row=shard.execute('SELECT payload FROM reuse_policies WHERE key=?',(policy,)).fetchone()
            if row is None or row[0]!=_json(identity['identity']) or _key(json.loads(row[0]))!=policy:
                raise ValueError('shard policy identity mismatch')
            row=shard.execute('SELECT witness,output_digest FROM reuse_roots WHERE policy=? AND key=?',
                              (policy,key)).fetchone()
            if row is None or row[0]!=_json(form) or row[1]!=result.get('output_digest'):
                raise ValueError('completed root missing or corrupted in shard')
            producer=identity['identity']['producer']
            requests={execution_cache_key(r,producer):r for r in leaf_requests(form)}
            dependencies={r[0] for r in shard.execute('SELECT leaf FROM reuse_dependencies WHERE policy=? AND root=?',
                                                     (policy,key))}
            if set(requests)!=dependencies:
                raise ValueError('completed root dependency mismatch')
            for leaf,request in requests.items():
                row=shard.execute('SELECT request,payload,digest,output_digest FROM reuse_leaves WHERE policy=? AND key=?',
                                  (policy,leaf)).fetchone()
                if row is None:
                    raise ValueError('missing verified leaf')
                if row[0]!=_json(request):
                    raise ValueError('verified leaf request mismatch')
                checked_key=(worker,policy,leaf,hashlib.sha256(row[1]).digest(),row[2],row[3])
                if checked_key in checked:
                    checked.move_to_end(checked_key)
                    continue
                _audit_leaf(row,request,producer,identity['identity']['limits'])
                checked[checked_key]=None
                if len(checked)>65536:
                    checked.popitem(last=False)
            counts['complete']+=1
        for (campaign,n),(length,ordinal) in cursors.items():
            if identities[campaign]['domain']!='current-primitive-baseline-v1':
                raise ValueError('explicit-input campaign cannot have enumeration cursors')
            count,maximum=admitted.get((campaign,n),(0,-1))
            if maximum!=length:
                raise ValueError('census cursor length disagrees with admitted words')
            sizes=[len(_small_layer(n,size)) if n<=2 else _layer(n,size).count
                   for size in range(length+1)]
            if ordinal>sizes[-1] or sum(sizes[:-1])+ordinal!=count:
                raise ValueError('census cursor skips or repeats admitted words')
        _assert_runner_source()
        _assert_package_sources()
        return counts
    finally:
        for shard in shards.values():
            shard.close()
        db.close()


def run(output, *, targets=None, words=None, workers=1, wall_seconds=3600,
        max_new_words=None, mode='reuse', limits=None, cache_bytes=256*1024**2,
        cache_entries=262144, allocator=None, retry_incomplete=False):
    """Run a bounded admission window and drain admitted work before returning.

    Explicit words are (rank, literal-word) pairs. Without them, enumerate the
    proved primitive baseline fairly by rank. Cache reuse is shard-local.
    """
    _assert_runner_source()
    _assert_package_sources()
    if (type(workers) is not int or workers < 1 or wall_seconds <= 0
            or mode not in ('reuse', 'standalone')
            or (max_new_words is not None and (type(max_new_words) is not int or max_new_words < 0))):
        raise ValueError('invalid worker, mode, duration, or admission limit')
    targets = dict(targets if targets is not None else {n:10 for n in range(2,9)})
    if any(type(n) is not int or n<1 or type(length) is not int or length<0 for n,length in targets.items()):
        raise ValueError('targets require positive ranks and nonnegative lengths')
    if words is not None:
        words = [(n, list(word)) for n,word in words]
        for n,word in words:
            input_normal_form(word,n=n)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    lock = (output/'coordinator.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BaseException:
        lock.close(); raise
    pools = []
    db = None
    started = perf_counter()
    run_id = uuid.uuid4().hex
    try:
        _probe_ledger(output/'census.sqlite3')
        db = sqlite3.connect(output/'census.sqlite3')
        _ledger(db)
        audit_started = perf_counter()
        audit_store(output)
        audit_seconds = perf_counter() - audit_started
        _assert_package_sources()
        settings = {'mode':mode, 'limits':asdict(limits or Limits()),
                    'cache_bytes':cache_bytes, 'cache_entries':cache_entries,
                    'allocator':None if allocator is None else str(allocator)}
        for worker in range(workers):
            pools.append(ProcessPoolExecutor(max_workers=1,
                         mp_context=multiprocessing.get_context('spawn'),
                         initializer=_initialize,
                         initargs=(str(output/'shards'/f'worker-{worker}.sqlite3'),settings)))
        identities = [pool.submit(_ready).result() for pool in pools]
        _assert_runner_source()
        _assert_package_sources()
        if any(item['runner_sha256']!=RUNNER_SHA256 for item in identities):
            raise ValueError('worker runner source differs from coordinator')
        if any(_json(item['identity'].get('reconstruction_sources'))!=_json(PACKAGE_SOURCE_SHA256) for item in identities):
            raise ValueError('worker reconstruction sources differ from coordinator import')
        if len({item['policy'] for item in identities}) != 1:
            raise ValueError('worker execution identities disagree')
        identity = {'schema':SCHEMA, 'worker':identities[0],
                    'runner_sha256':RUNNER_SHA256,
                    'domain':'explicit-literal-requests' if words is not None else 'current-primitive-baseline-v1'}
        campaign = _key(identity)
        with db:
            db.execute('INSERT OR IGNORE INTO reuse_campaigns VALUES (?,?)',(campaign,_json(identity)))
            if db.execute('SELECT identity FROM reuse_campaigns WHERE key=?',(campaign,)).fetchone()[0] != _json(identity):
                raise ValueError('stored campaign identity mismatch')
            db.execute("UPDATE reuse_jobs SET status='interrupted' WHERE campaign=? AND status='running'",(campaign,))
        pending = _pending_jobs(db,campaign,retry_incomplete)
        explicit = deque(words or [])
        ranks = deque(sorted(targets))
        last = db.execute('SELECT last_rank FROM reuse_round_robin WHERE campaign=?',(campaign,)).fetchone()
        if last and ranks:
            order=sorted(targets)
            ranks=deque([n for n in order if n>last[0]]+[n for n in order if n<=last[0]])
        cursors = {n: [0 if n==1 else max(1,n-1),0] for n in targets}
        for n,length,ordinal in db.execute('SELECT n,length,ordinal FROM reuse_cursors WHERE campaign=?',(campaign,)):
            if n in cursors:
                cursors[n] = [length,int(ordinal)]
        active = {}
        idle = deque(range(workers))
        admitted = 0
        snapshot_counts = {}
        failed_workers = set()
        stop_reason = None
        source_snapshot = PACKAGE_SOURCE_SHA256
        admission_started = perf_counter()
        setup_seconds = admission_started - started
        while active or idle:
            while idle and stop_reason is None and perf_counter()-admission_started < wall_seconds and (max_new_words is None or admitted<max_new_words):
                item = None
                retry = next(pending,None)
                if retry is not None:
                    key,n,encoded = retry
                    item = (key,n,json.loads(encoded),None)
                elif words is not None:
                    while explicit:
                        n,word = explicit.popleft()
                        key = _key({'n':n,'word':word})
                        if db.execute('SELECT 1 FROM reuse_jobs WHERE campaign=? AND key=?',(campaign,key)).fetchone() is None:
                            item = key,n,word,None
                            break
                else:
                    while ranks:
                        n = ranks.popleft()
                        candidate = _next(n,cursors[n],targets[n])
                        if candidate is None:
                            continue
                        word,next_cursor,_ = candidate
                        ranks.append(n)
                        cursors[n] = next_cursor
                        key = _key({'n':n,'word':word})
                        item = key,n,word,next_cursor
                        break
                if item is None:
                    break
                key,n,word,cursor = item
                worker = idle.popleft()
                with db:
                    db.execute('INSERT INTO reuse_jobs VALUES (?,?,?,?,?,?,?) '
                               'ON CONFLICT(campaign,key) DO UPDATE SET status=excluded.status,shard=excluded.shard',
                               (campaign,key,n,_json(word),'running',worker,None))
                    if cursor is not None:
                        db.execute('INSERT INTO reuse_round_robin VALUES (?,?) '
                                   'ON CONFLICT(campaign) DO UPDATE SET last_rank=excluded.last_rank',(campaign,n))
                        db.execute('INSERT INTO reuse_cursors VALUES (?,?,?,?) '
                                   'ON CONFLICT(campaign,n) DO UPDATE SET length=excluded.length,ordinal=excluded.ordinal',
                                   (campaign,n,cursor[0],str(cursor[1])))
                admitted += 1
                try:
                    future = pools[worker].submit(_evaluate,n,word)
                except BrokenProcessPool as error:
                    # Admission/cursor are already committed. Leave this exact
                    # job retryable and stop consuming the enumeration.
                    result={'error':type(error).__name__,'message':str(error),'stage':'submission'}
                    with db:
                        db.execute("UPDATE reuse_jobs SET status='interrupted',result=? WHERE campaign=? AND key=?",
                                   (_json(result),campaign,key))
                        db.execute("INSERT INTO reuse_job_attempts(campaign,job,status,result) VALUES (?,?,'interrupted',?)",
                                   (campaign,key,_json(result)))
                    failed_workers.add(worker)
                    stop_reason='worker_failure'
                    break
                active[future] = worker,key
            if not active:
                break
            ready,_ = wait(active,timeout=1,return_when=FIRST_COMPLETED)
            for future in ready:
                worker,key = active.pop(future)
                try:
                    result = future.result()
                    if result['policy'] != identities[worker]['policy'] or result.get('root') != key:
                        raise ValueError('worker result policy mismatch')
                    telemetry = result.pop('worker_telemetry', None)
                    state = result['status']
                    if state not in ('complete','incomplete'):
                        raise ValueError('invalid worker completion status')
                except BrokenProcessPool as error:
                    state = 'interrupted'
                    telemetry = None
                    result = {'error':type(error).__name__,'message':str(error)}
                    failed_workers.add(worker)
                    stop_reason='worker_failure'
                except Exception as error:
                    state = 'failed'
                    telemetry = None
                    result = {'error':type(error).__name__,'message':str(error)}
                    stop_reason='worker_result_failure'
                with db:
                    if telemetry is not None:
                        db.execute('INSERT INTO reuse_worker_telemetry VALUES (?,?,?,?) '
                                   'ON CONFLICT(campaign,run,worker) DO UPDATE SET payload=excluded.payload',
                                   (campaign,run_id,worker,_json(telemetry)))
                    db.execute('UPDATE reuse_jobs SET status=?,result=? WHERE campaign=? AND key=?',
                               (state,_json(result),campaign,key))
                    db.execute('INSERT INTO reuse_job_attempts(campaign,job,status,result) VALUES (?,?,?,?)',
                               (campaign,key,state,_json(result)))
                if worker not in failed_workers:
                    idle.append(worker)
            if reconstruction_sources()!=source_snapshot:
                raise ValueError('sources changed during census campaign')
            _assert_runner_source()
        # Exact shard totals are scanned once at the end, not for every word.
        for worker,pool in enumerate(pools):
            if worker not in failed_workers:
                try:
                    snapshot_counts[worker]=pool.submit(_snapshot_counts).result()
                except BrokenProcessPool:
                    failed_workers.add(worker)
                    stop_reason='worker_failure'
        counts = dict(db.execute('SELECT status,count(*) FROM reuse_jobs WHERE campaign=? GROUP BY status',(campaign,)))
        _assert_runner_source()
        _assert_package_sources()
        finished = perf_counter()
        summary = {'schema':SCHEMA,'campaign':campaign,'run':run_id,'admitted':admitted,'jobs':counts,
                   'worker_counts':snapshot_counts,'wall_seconds':finished-started,
                   'setup_seconds':setup_seconds,'audit_seconds':audit_seconds,
                   'admission_and_drain_seconds':finished-admission_started,
                   'admission_window_seconds':wall_seconds,
                   'stop_reason':stop_reason or 'admission_finished',
                   'unavailable_workers':sorted(failed_workers),
                   'coverage':'completed requested outputs via checked reconstruction; no full ambient traces',
                   'reuse_scope':'shard-local verified leaves and native completion caches',
                   'count_scope':'worker totals are shard-local and may overlap across workers'}
        (output/'summary.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
        return summary
    finally:
        for pool in pools:
            pool.shutdown(wait=True,cancel_futures=True)
        if db is not None:
            db.close()
        lock.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--targets',type=json.loads)
    parser.add_argument('--words',type=Path,help='JSON list of [rank, word] pairs')
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--wall-seconds',type=float,default=3600,
                        help='admission window after full audit and worker setup; admitted jobs then drain')
    parser.add_argument('--max-new-words',type=int)
    parser.add_argument('--mode',choices=('reuse','standalone'),default='reuse')
    parser.add_argument('--limits',type=json.loads,default={})
    parser.add_argument('--cache-mib',type=int,default=256)
    parser.add_argument('--verifier-allocator',type=Path)
    parser.add_argument('--retry-incomplete',action='store_true')
    parser.add_argument('--audit-only',action='store_true')
    args=parser.parse_args()
    if args.audit_only:
        print(json.dumps(audit_store(args.output),indent=2))
        return
    print(json.dumps(run(args.output,targets=None if args.targets is None else {int(n):l for n,l in args.targets.items()},
                         words=None if args.words is None else json.loads(args.words.read_text()),
                         workers=args.workers,wall_seconds=args.wall_seconds,max_new_words=args.max_new_words,
                         mode=args.mode,limits=Limits(**args.limits),cache_bytes=args.cache_mib*1024**2,
                         allocator=args.verifier_allocator,retry_incomplete=args.retry_incomplete),indent=2))


if __name__=='__main__':
    main()
