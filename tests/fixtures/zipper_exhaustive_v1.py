"""Exhaust literal words with collected arithmetic and independent rank prefixes.

One parent writes one SQLite database. Spawned workers share a CollectedSession
within each process, reduce one literal word, JSON-reload it and serially replay
its certificate before returning it. Certificates are stored as zlib-compressed
UTF-8 JSON, with their uncompressed SHA-256 and byte length. No word is removed
by cancellation or isotopy. Ambient strand count and output slider size are
reported separately.

Each rank commits only its consecutive prefix; a limited word pauses that rank
while other ranks continue. At most twice the worker count is outstanding,
including completed results waiting behind an earlier word. Resume audits all
hashes, literal inputs, configurations and cursors. Recorded production replay
is not repeated unless --replay-audit is requested. Larger limits can retry the
first incomplete word. The wall budget stops admissions/commits between words;
shutdown waits for bounded speculative work and discards uncommitted results.

Example:
  python experiments/zipper_exhaustive_collected.py results --workers 8 \
    --targets '{"2":10,"3":10,"4":10,"5":10,"6":10,"7":10,"8":10}' \
    --wall-seconds 1200
"""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import asdict
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import sqlite3
import sys
import time
import zlib
from typing import Any

from zipper_algebra.zipping_algorithm import collected, enumeration
from zipper_algebra.zipping_algorithm.catalog import Catalog, _data, SCHEMA_SHA256
from zipper_algebra.zipping_algorithm.enumeration import _json, _sha, _word_at, word_count
from zipper_algebra.zipping_algorithm.interpreter import Limits

CODEC = 'zlib-utf8-json-v1'
CHECK = 'json-roundtrip-collected-serial-replay-v1'
DEFAULT_TARGETS = {n:10 for n in range(2,9)}
_session = None
_worker_limits = None
_worker_identity = None


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(batch_size):
    return {'schema':'zipper-collected-literal-enumeration-v1',
            'algorithm_source_sha256':collected._source(),
            'catalog_sha256':Catalog.bundled().sha256,'catalog_schema_sha256':SCHEMA_SHA256,
            'runner_sha256':_hash(__file__),'enumeration_sha256':_hash(enumeration.__file__),
            'batch_size':batch_size,'coefficient_ring':'Z[q,q^-1]',
            'codec':CODEC,'production_check':CHECK,'trace':False,'term_workers':1}


def _assert_source(expected):
    if (collected._source()!=expected['algorithm_source_sha256']
            or _hash(__file__)!=expected['runner_sha256']
            or _hash(enumeration.__file__)!=expected['enumeration_sha256']
            or hashlib.sha256(_data('full-catalog.json')).hexdigest()!=expected['catalog_sha256']
            or hashlib.sha256(_data('schema.json')).hexdigest()!=expected['catalog_schema_sha256']):
        raise ValueError('enumeration source or catalog changed')


def _initialize(limits, expected):
    global _session, _worker_limits, _worker_identity
    _session=collected.CollectedSession()
    _worker_limits=Limits(**limits)
    _worker_identity=expected


def _certificate_identity(certificate,n,word,limits,expected):
    if (certificate['schema']!='zipper-collected-run-v1'
            or certificate['n']!=n or certificate['input_word']!=word
            or certificate['source_sha256']!=expected['algorithm_source_sha256']
            or certificate['catalog_sha256']!=expected['catalog_sha256']
            or certificate['limits']!=limits or certificate['trace_enabled'] is not False
            or certificate['coefficient_ring']!=expected['coefficient_ring']
            or certificate['execution']['batch_size']!=expected['batch_size']
            or certificate['status'] not in ('complete','incomplete','control_obligation')):
        raise ValueError('certificate identity does not match literal input/configuration')


def _slider_size(certificate):
    return max((component['size'] for term in certificate['output']
                for component in term['components']),default=0)


def _compute(n,word):
    assert _session is not None and _worker_limits is not None and _worker_identity is not None
    start=time.perf_counter_ns()
    certificate=_session.reduce(word,n=n,limits=_worker_limits,trace=False,
                                workers=1,batch_size=_worker_identity['batch_size'])
    producer=time.perf_counter_ns()-start
    encoded=_json(certificate)
    start=time.perf_counter_ns()
    if _session.check(json.loads(encoded)) is not True:
        raise ValueError('certificate failed exact replay')
    checker=time.perf_counter_ns()-start
    _certificate_identity(certificate,n,word,asdict(_worker_limits),_worker_identity)
    return encoded,producer,checker


def _process_pool(**kwargs):
    return ProcessPoolExecutor(mp_context=multiprocessing.get_context('spawn'),**kwargs)


def _advance(n,cursor):
    length,ordinal=cursor
    return [length+1,0] if ordinal+1==word_count(n,length) else [length,ordinal+1]


def _schema(db):
    db.execute('PRAGMA journal_mode=WAL')
    db.execute('CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL)')
    db.execute('CREATE TABLE IF NOT EXISTS configurations(hash TEXT PRIMARY KEY,config TEXT NOT NULL)')
    db.execute('CREATE TABLE IF NOT EXISTS ranks(n INTEGER PRIMARY KEY,length INTEGER NOT NULL,ordinal INTEGER NOT NULL)')
    db.execute('''CREATE TABLE IF NOT EXISTS words(
        n INTEGER,length INTEGER,ordinal INTEGER,word TEXT,input_sha TEXT,status TEXT,
        diagnostic TEXT,metrics TEXT,certificate BLOB,certificate_sha TEXT,certificate_bytes INTEGER,
        codec TEXT,config_sha TEXT,producer_ns INTEGER,checker_ns INTEGER,output_slider_size INTEGER,
        production_check TEXT,PRIMARY KEY(n,length,ordinal))''')


def audit(db,expected,*,replay=False):
    """Validate every saved binding; optionally repeat its exact arithmetic."""
    cursors={n:[length,ordinal] for n,length,ordinal in db.execute('SELECT n,length,ordinal FROM ranks')}
    if any(type(n) is not int or n<2 or any(type(v) is not int or v<0 for v in cursor) for n,cursor in cursors.items()):
        raise ValueError("invalid saved rank cursor")
    next_cursor={n:[0,0] for n in cursors};blocked=set();count=0
    configs={};session=collected.CollectedSession() if replay else None
    for row in db.execute('SELECT * FROM words ORDER BY n,length,ordinal'):
        n,length,ordinal,word_json,input_sha,status,diagnostic,metrics,blob,digest,size,codec,config_sha,producer,checker,slider,check=row
        if n not in cursors or n in blocked or [length,ordinal]!=next_cursor[n]:
            raise ValueError('saved rows are not consecutive literal rank prefixes')
        word=list(_word_at(n,length,ordinal))
        if word_json!=_json(word) or input_sha!=_sha(_json({'n':n,'word':word})):
            raise ValueError('saved literal input binding mismatch')
        if codec!=CODEC or check!=CHECK or type(size) is not int or size<0:
            raise ValueError('saved certificate codec/check provenance mismatch')
        try:
            raw=zlib.decompress(blob)
            if len(raw)!=size or hashlib.sha256(raw).hexdigest()!=digest:
                raise ValueError('saved certificate length/hash mismatch')
            certificate=json.loads(raw)
        except (zlib.error,UnicodeError,json.JSONDecodeError,TypeError) as error:
            raise ValueError('saved compressed certificate is invalid') from error
        if config_sha not in configs:
            record=db.execute('SELECT config FROM configurations WHERE hash=?',(config_sha,)).fetchone()
            if not record or _sha(record[0])!=config_sha:
                raise ValueError('saved configuration hash mismatch')
            configs[config_sha]=json.loads(record[0])
        config=configs[config_sha]
        if (config['identity_sha256']!=_sha(_json(expected))
                or str(n) not in config['targets'] or length>config['targets'][str(n)]):
            raise ValueError('saved configuration source binding mismatch')
        _certificate_identity(certificate,n,word,config['limits'],expected)
        if (status!=certificate['status'] or diagnostic!=_json(certificate['diagnostic'])
                or metrics!=_json(certificate['metrics']) or slider!=_slider_size(certificate)
                or type(producer) is not int or producer<0 or type(checker) is not int or checker<0):
            raise ValueError('saved certificate columns differ from its contents')
        if replay:
            assert session is not None
            if session.check(certificate) is not True:
                raise ValueError('saved certificate failed exact replay')
        count+=1
        if status=='complete':next_cursor[n]=_advance(n,next_cursor[n])
        else:blocked.add(n)
    if next_cursor!=cursors:
        raise ValueError('saved rank cursor differs from its literal prefix')
    return {'words':count,'hash_and_binding_checked':count,'replayed':count if replay else 0}


def summarize(db,targets,reasons):
    result=[]
    for n,length,ordinal in db.execute('SELECT n,length,ordinal FROM ranks ORDER BY n'):
        layers=[];through=-1
        for L,count,complete in db.execute("SELECT length,count(*),sum(status='complete') FROM words WHERE n=? GROUP BY length ORDER BY length",(n,)):
            expected=word_count(n,L);finished=count==complete==expected
            if finished and L==through+1:through=L
            layers.append({'length':L,'expected_words':expected,'visited_words':count,
                           'complete_words':complete,'layer_complete':finished})
        totals=db.execute("SELECT count(*),sum(status='complete'),max(output_slider_size),sum(length(certificate)),sum(certificate_bytes) FROM words WHERE n=?",(n,)).fetchone()
        unresolved=db.execute("SELECT word,status,diagnostic FROM words WHERE n=? AND length=? AND ordinal=?",(n,length,ordinal)).fetchone()
        result.append({'n':n,'target_length':targets.get(n),'complete_through_length':through,
            'visited_words':totals[0],'complete_words':totals[1] or 0,'maximum_output_slider_size':totals[2] or 0,
            'compressed_certificate_bytes':totals[3] or 0,'raw_certificate_bytes':totals[4] or 0,
            'next_cursor':{'length':length,'ordinal':ordinal},'stop_reason':reasons.get(n),
            'first_unresolved':None if not unresolved else {'word':json.loads(unresolved[0]),'status':unresolved[1],'diagnostic':json.loads(unresolved[2])},'layers':layers})
    return result


def exhaust(output,*,targets=None,workers=8,batch_size=64,wall_seconds=86400,limits=None,
            max_new_words=None,replay_audit=False,audit_only=False,_executor_factory=_process_pool):
    targets=DEFAULT_TARGETS.copy() if targets is None else dict(targets)
    if not targets or any(type(n) is not int or n<2 or type(L) is not int or L<0 for n,L in targets.items()):
        raise ValueError('integer rank >=2 and nonnegative target length required')
    if type(workers) is not int or workers<1 or type(batch_size) is not int or batch_size<1:
        raise ValueError('positive integer workers and batch_size required')
    if type(wall_seconds) not in (int,float) or not math.isfinite(wall_seconds) or wall_seconds<=0:
        raise ValueError('positive finite wall_seconds required')
    if max_new_words is not None and (type(max_new_words) is not int or max_new_words<0):
        raise ValueError('nonnegative max_new_words required')
    limits=Limits(max_steps=200000,max_terms=100000,expansion_budget=1000000) if limits is None else limits
    if not isinstance(limits,Limits):raise ValueError('limits must be a Limits record')
    output=Path(output);output.mkdir(parents=True,exist_ok=True);database=output/'words.sqlite3'
    expected=identity(batch_size)
    config={'targets':targets,'workers':workers,'max_in_flight':2*workers,'limits':asdict(limits),
            'wall_seconds':wall_seconds,'identity_sha256':_sha(_json(expected)),
            'python':sys.version,'created_utc':datetime.now(timezone.utc).isoformat()}
    config_sha=_sha(_json(config));lock=(output/'.parallel.lock').open('a+')
    try:fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError as error:
        lock.close();raise ValueError('another writer owns this campaign') from error
    db=sqlite3.connect(database);pool=None;new=0;peak=0
    queues: dict[int, deque[Any]] = {}
    discarded: set[Any] = set()
    try:
        _schema(db)
        old=db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()
        if old and old[0]!=_json(expected):raise ValueError('enumeration identity changed; use a new database')
        with db:
            db.execute("INSERT OR IGNORE INTO metadata VALUES('identity',?)",(_json(expected),))
        audit_start=time.perf_counter();audited=audit(db,expected,replay=replay_audit)
        audit_seconds=time.perf_counter()-audit_start
        with db:
            db.execute('INSERT OR IGNORE INTO configurations VALUES(?,?)',(config_sha,_json(config)))
            db.execute("INSERT OR REPLACE INTO metadata VALUES('last_config',?)",(_json(config),))
            for n in targets:db.execute('INSERT OR IGNORE INTO ranks VALUES(?,0,0)',(n,))
        cursors={n:[L,o] for n,L,o in db.execute('SELECT n,length,ordinal FROM ranks') if n in targets}
        next_submit={n:cursor[:] for n,cursor in cursors.items()};queues={n:deque() for n in targets}
        reasons={};active=set(targets)
        for n,cursor in cursors.items():
            old=db.execute('SELECT status,certificate FROM words WHERE n=? AND length=? AND ordinal=?',(n,*cursor)).fetchone()
            if old:
                before=json.loads(zlib.decompress(old[1]))['limits'];after=asdict(limits)
                if old[0]=='control_obligation':reasons[n]='control_obligation'
                elif before==after or any(after[k]<v for k,v in before.items()):reasons[n]='word_limit_needs_larger_limits'
            if cursor[0]>targets[n]:reasons[n]='target_length'
            if n in reasons:active.discard(n)
        started=time.perf_counter();last_snapshot=started;rotation=deque(sorted(targets));stop='finished'
        def snapshot():
            data={'schema':expected['schema'],'identity':expected,'database':str(database.resolve()),
                  'config':config,'audit':audited,'audit_seconds':audit_seconds,
                  'new_or_resumed_words':new,'new_work_seconds':time.perf_counter()-started,
                  'peak_in_flight':peak,'stop_reason':stop,'ranks':summarize(db,targets,reasons)}
            temporary=output/'summary.json.tmp';temporary.write_text(_json(data)+'\n');temporary.replace(output/'summary.json')
            return data
        if active and not audit_only and max_new_words!=0:
            pool=_executor_factory(max_workers=workers,initializer=_initialize,initargs=(asdict(limits),expected))
        while active and not audit_only and max_new_words!=0:
            assert pool is not None
            if time.perf_counter()-started>=wall_seconds:stop='wall_budget';break
            if max_new_words is not None and new>=max_new_words:stop='word_budget';break
            discarded={future for future in discarded if not future.done()}
            occupancy=sum(len(q) for q in queues.values())+len(discarded)
            for _ in range(2*workers):
                if occupancy>=2*workers or (max_new_words is not None and new+occupancy>=max_new_words):break
                chosen=None
                for _ in range(len(rotation)):
                    n=rotation[0];rotation.rotate(-1)
                    if n in active and next_submit[n][0]<=targets[n]:chosen=n;break
                if chosen is None:break
                n=chosen;cursor=next_submit[n][:];word=list(_word_at(n,*cursor))
                queues[n].append((cursor,word,pool.submit(_compute,n,word)))
                next_submit[n]=_advance(n,cursor);occupancy+=1;peak=max(peak,occupancy)
            ready=[n for n in rotation if queues[n] and queues[n][0][2].done()]
            if not ready:
                futures=[job[2] for q in queues.values() for job in q if not job[2].done()]+list(discarded)
                if not futures:raise RuntimeError('active rank has no pending literal input')
                wait(futures,timeout=.2,return_when=FIRST_COMPLETED);continue
            _assert_source(expected)
            for n in ready:
                cursor,word,future=queues[n].popleft()
                if cursor!=cursors[n]:raise ValueError('worker is not the next literal input')
                encoded,producer,checker=future.result();certificate=json.loads(encoded)
                _certificate_identity(certificate,n,word,asdict(limits),expected)
                status=certificate['status'];following=_advance(n,cursor) if status=='complete' else cursor
                raw=encoded.encode();blob=zlib.compress(raw)
                with db:
                    db.execute('INSERT OR REPLACE INTO words VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                        (n,*cursor,_json(word),_sha(_json({'n':n,'word':word})),status,
                         _json(certificate['diagnostic']),_json(certificate['metrics']),blob,_sha(encoded),len(raw),CODEC,
                         config_sha,producer,checker,_slider_size(certificate),CHECK))
                    db.execute('UPDATE ranks SET length=?,ordinal=? WHERE n=?',(*following,n))
                cursors[n]=following;new+=1
                if status!='complete' or following[0]>targets[n]:
                    reasons[n]='target_length' if status=='complete' else status
                    active.discard(n)
                    for _,_,later in queues[n]:
                        later.cancel()
                        if not later.done():discarded.add(later)
                    queues[n].clear()
            if time.perf_counter()-last_snapshot>=5:
                snapshot();last_snapshot=time.perf_counter()
        if audit_only:stop='audit_only'
        elif max_new_words==0:stop='word_budget'
        for q in queues.values():
            for _,_,future in q:future.cancel()
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True);pool=None
        _assert_source(expected)
        return snapshot()
    finally:
        for q in queues.values():
            for _,_,future in q:future.cancel()
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True)
        db.close();fcntl.flock(lock.fileno(),fcntl.LOCK_UN);lock.close()


def status(output):
    """Read one consistent SQLite snapshot without claiming the process is alive."""
    database=Path(output).resolve()/"words.sqlite3"
    with sqlite3.connect(database.as_uri()+"?mode=ro",uri=True) as db:
        db.execute("BEGIN")
        metadata=dict(db.execute("SELECT key,value FROM metadata"))
        expected=json.loads(metadata["identity"]);config=json.loads(metadata["last_config"])
        targets={int(n):L for n,L in config["targets"].items()}
        return {"schema":expected["schema"],"identity":expected,"database":str(database),
                "snapshot_utc":datetime.now(timezone.utc).isoformat(),"config":config,
                "ranks":summarize(db,targets,{})}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output');parser.add_argument('--targets',type=json.loads)
    parser.add_argument('--workers',type=int,default=8);parser.add_argument('--batch-size',type=int,default=64)
    parser.add_argument('--wall-seconds',type=float,default=86400);parser.add_argument('--limits',type=json.loads)
    parser.add_argument('--max-new-words',type=int);parser.add_argument('--replay-audit',action='store_true')
    parser.add_argument('--audit-only',action='store_true');parser.add_argument('--status',action='store_true')
    args=vars(parser.parse_args())
    if args.pop('status'):
        print(json.dumps(status(args['output']),indent=2),flush=True);return
    if args['targets'] is not None:args['targets']={int(n):L for n,L in args['targets'].items()}
    if args['limits'] is not None:args['limits']=Limits(**args['limits'])
    print(json.dumps(exhaust(**args),indent=2),flush=True)


if __name__=='__main__':main()
