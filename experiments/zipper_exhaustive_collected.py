"""Exhaust literal words with collected arithmetic and independent rank prefixes.

One parent writes one SQLite database. Each spawned worker retains a C++
NativeSession and a separate Python or C++ verification session. It reduces one literal
word, JSON-reloads it and independently replays the complete calculation before
returning it. Python production modes remain available for historical replay.
Certificates are stored as zlib-compressed
UTF-8 JSON, with their uncompressed SHA-256 and byte length. No word is removed
by cancellation or isotopy. Ambient strand count and output slider size are
reported separately.

Each rank commits only its consecutive prefix; a limited word pauses that rank
while other ranks continue. At most twice the worker count is outstanding,
including completed results waiting behind an earlier word. Resume audits all
hashes, literal inputs, configurations and cursors. Recorded production replay
is not repeated unless --replay-audit is requested. Larger limits can retry the
first incomplete word. The wall budget stops admissions, then drains checked
rank prefixes already admitted. Per-rank windows prevent a slow rank from
occupying the entire worker queue. This bounds queued work per rank; it does
not reserve a physical worker for each rank. Draining may exceed the admission
deadline by the cost of all already admitted computations. SIGUSR1 requests the
same graceful admission stop and drain; the prior signal handler is restored.
An optional absolute UTC admission deadline also applies during startup audits,
so migration time cannot extend a previously scheduled admission window.

Explicit --migrate audits and backs up the existing database before registering
a new source epoch. Legacy identities and configurations remain immutable;
completed rows are retained, and retried checkpoints are archived verbatim.
Audit-only uses a read-only connection. Replay of historical certificates needs
the corresponding source checkout. Worker cache/memory telemetry is separate
from certificates; logical rule histograms describe the current saved prefixes.

Example:
  python experiments/zipper_exhaustive_collected.py results --workers 8 \
    --targets '{"2":10,"3":10,"4":10,"5":10,"6":10,"7":10,"8":10}' \
    --wall-seconds 1200
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import asdict
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import multiprocessing
import os
import resource
import signal
import threading
from multiprocessing.util import Finalize
from pathlib import Path
import sqlite3
import sys
import time
import zlib
from typing import Any

from zipper_algebra.zipping_algorithm import collected, enumeration
from zipper_algebra.zipping_algorithm.catalog import Catalog, _data, SCHEMA_SHA256
from zipper_algebra.zipping_algorithm.enumeration import _sha, _word_at, word_count
from zipper_algebra.zipping_algorithm.interpreter import Limits

CODEC = 'zlib-utf8-json-v1'
CHECK = 'json-roundtrip-exact-serial-replay-v2'
NATIVE_CHECK = 'json-roundtrip-independent-python-cached-replay-v1'
CPP_CHECK = 'json-roundtrip-independent-cpp-cached-replay-v1'
VERIFIER_BUILD_FIELDS = ('abi','binary_sha256','compiler','compiler_version','flags',
                         'wrapper_source_sha256','link_libraries','libraries')
NATIVE_BUILD_FIELDS = ('abi','binary_sha256','compiler','compiler_version','flags',
                       'wrapper_source_sha256','link_libraries')
DEFAULT_TARGETS = {n:10 for n in range(2,9)}
_session = None
_verifier = None
_worker_limits = None
_worker_identity = None
_worker_sequence = 0


def _integer_text(value):
    try:return str(value)
    except ValueError:
        from zipper_algebra.zipping_algorithm.native import _integer_decimal
        return _integer_decimal(value)


def _json(value):
    """Keep canonical JSON numeric integers exact beyond Python's digit limit."""
    try:return enumeration._json(value)
    except ValueError as error:
        if 'Exceeds the limit' not in str(error):raise
    active=set()
    def encode(item):
        if item is None or isinstance(item,(str,bool,float)):
            return json.dumps(item,allow_nan=False)
        if isinstance(item,int):return _integer_text(item)
        if isinstance(item,(list,tuple,dict)):
            marker=id(item)
            if marker in active:raise ValueError('Circular reference detected')
            active.add(marker)
            try:
                if isinstance(item,dict):
                    entries=[]
                    for key,child in sorted(item.items()):
                        if isinstance(key,str):label=key
                        elif key is None:label='null'
                        elif isinstance(key,bool):label='true' if key else 'false'
                        elif isinstance(key,int):label=_integer_text(key)
                        elif isinstance(key,float):label=json.dumps(key,allow_nan=False)
                        else:raise TypeError('JSON object keys must be primitive values')
                        entries.append(json.dumps(label)+':'+encode(child))
                    return '{'+','.join(entries)+'}'
                return '['+','.join(encode(child) for child in item)+']'
            finally:active.remove(marker)
        raise TypeError('Object is not JSON serializable')
    return encode(value)


def _loads(value):
    def integer(text):
        try:return int(text)
        except ValueError:
            from zipper_algebra.zipping_algorithm.native import _decimal_integer
            return _decimal_integer(text)
    return json.loads(value,parse_int=integer)


def _hash(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def _engine(execution):
    if execution == 'native':
        from zipper_algebra.zipping_algorithm import native
        return native
    if execution == 'memoized':
        from zipper_algebra.zipping_algorithm import memoized
        return memoized
    return collected


def _source(execution):
    engine = _engine(execution)
    return engine.source_identity() if execution == 'native' else engine._source()


def _verifier_source(backend):
    if backend == 'cpp':
        from zipper_algebra.zipping_algorithm import native_verifier
        return {**native_verifier.source_identity(), 'native_verifier.py':_hash(native_verifier.__file__)}
    from zipper_algebra.zipping_algorithm import independent_replay
    return independent_replay.source_identity()


def _verifier_build(metadata):
    return {key:metadata[key] for key in VERIFIER_BUILD_FIELDS}


def _verifier_allocator(path):
    if path is None:return None
    from zipper_algebra.zipping_algorithm.native_verifier import _allocator_identity
    return _allocator_identity(Path(path))


def _allocator_runtime(identity):
    return None if identity is None else dict(identity,mapped_verified=True)


def identity(batch_size, execution='native', verifier='python', verifier_allocator=None):
    engine = _engine(execution)
    extra = {}
    if execution == 'native':
        _, build, _ = engine.build_engine()
        build['wrapper_source_sha256'] = _hash(engine.__file__)
        extra = {'native_build':{key:build[key] for key in NATIVE_BUILD_FIELDS},
                 'python_reference_source_sha256':collected._source(),
                 'verifier_source_sha256':_verifier_source(verifier)}
        if verifier == 'cpp':
            from zipper_algebra.zipping_algorithm import native_verifier
            _, verifier_build, _ = native_verifier.build_verifier()
            verifier_build['wrapper_source_sha256']=_hash(native_verifier.__file__)
            extra.update(verifier_backend='cpp',verifier_build=_verifier_build(verifier_build))
            if verifier_allocator is not None:
                extra['verifier_allocator']=_verifier_allocator(verifier_allocator)
    return {'schema':'zipper-versioned-literal-enumeration-v2',
            'algorithm_source_sha256':_source(execution),
            'certificate_schema':{'native':'zipper-native-run-v1','memoized':'zipper-memoized-run-v1',
                                  'collected':'zipper-collected-run-v1'}[execution],
            'execution':execution,
            'catalog_sha256':Catalog.bundled().sha256,'catalog_schema_sha256':SCHEMA_SHA256,
            'runner_sha256':_hash(__file__),'enumeration_sha256':_hash(enumeration.__file__),
            'batch_size':batch_size,'coefficient_ring':'Z[q,q^-1]',
            'codec':CODEC,'production_check':(CPP_CHECK if verifier=='cpp' else NATIVE_CHECK) if execution=='native' else CHECK,
            'trace':False,'term_workers':1,**extra}


def _assert_source(expected):
    if (_source(expected.get('execution','collected'))!=expected['algorithm_source_sha256']
            or _hash(__file__)!=expected['runner_sha256']
            or _hash(enumeration.__file__)!=expected['enumeration_sha256']
            or hashlib.sha256(_data('full-catalog.json')).hexdigest()!=expected['catalog_sha256']
            or hashlib.sha256(_data('schema.json')).hexdigest()!=expected['catalog_schema_sha256']):
        raise ValueError('enumeration source or catalog changed')
    if expected.get('execution') == 'native':
        from zipper_algebra.zipping_algorithm import native
        if (collected._source()!=expected['python_reference_source_sha256']
                or _verifier_source(expected.get('verifier_backend','python'))!=expected['verifier_source_sha256']
                or _hash(native.__file__)!=expected['native_build']['wrapper_source_sha256']):
            raise ValueError('native wrapper or independent verifier source changed')
        allocator=expected.get('verifier_allocator')
        if allocator is not None and _verifier_allocator(allocator['path'])!=allocator:
            raise ValueError('independent verifier allocator identity changed')


def _new_session(expected, cache_bytes=64*1024*1024, cache_entries=4096, rule_cache_bytes=4*1024*1024, rule_cache_entries=256):
    engine = _engine(expected.get('execution','collected'))
    if expected.get('execution') == 'native':
        session=engine.NativeSession(cache_bytes=cache_bytes,cache_entries=cache_entries)
        if {key:session.metadata[key] for key in NATIVE_BUILD_FIELDS} != expected['native_build']:
            session.close()
            raise ValueError('native executable build differs from campaign identity')
        return session
    if expected.get('execution') == 'memoized':
        return engine.MemoizedSession(cache_bytes=cache_bytes,cache_entries=cache_entries,rule_cache_bytes=rule_cache_bytes,rule_cache_entries=rule_cache_entries)
    return collected.CollectedSession()


def _initialize(limits, expected, cache_bytes, cache_entries, rule_cache_bytes, rule_cache_entries,
                verifier_cache_bytes=256*1024*1024, verifier_cache_entries=262144):
    global _session, _verifier, _worker_limits, _worker_identity, _worker_sequence
    _worker_sequence=0
    _session=_new_session(expected,cache_bytes,cache_entries,rule_cache_bytes,rule_cache_entries)
    _verifier=None
    if expected.get('execution')=='native':
        try:
            if expected.get('verifier_backend','python') == 'cpp':
                from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier
                _verifier=CppVerifier(producer_metadata=_session.metadata,
                                      cache_bytes=verifier_cache_bytes,cache_entries=verifier_cache_entries,
                                      allocator_library=(Path(expected['verifier_allocator']['path']) if expected.get('verifier_allocator') else None))
                if _verifier_build(_verifier.metadata)!=expected['verifier_build']:
                    raise ValueError('independent verifier build differs from campaign identity')
                if _verifier.metadata.get('allocator_library')!=_allocator_runtime(expected.get('verifier_allocator')):
                    raise ValueError('independent verifier allocator differs from campaign identity')
                Finalize(_verifier,_verifier.close,exitpriority=20)
            else:
                from zipper_algebra.zipping_algorithm.independent_replay import IndependentVerifier
                _verifier=IndependentVerifier(cache_bytes=verifier_cache_bytes,cache_entries=verifier_cache_entries,
                                              rule_cache_bytes=rule_cache_bytes,rule_cache_entries=rule_cache_entries)
            Finalize(_session,_session.close,exitpriority=10)
        except BaseException:
            try:
                if _verifier is not None and hasattr(_verifier,'close'):_verifier.close()
            finally:_session.close()
            raise
    _worker_limits=Limits(**limits)
    _worker_identity=expected


def _certificate_identity(certificate,n,word,limits,expected):
    if (certificate['schema']!=expected.get('certificate_schema','zipper-collected-run-v1')
            or certificate['n']!=n or certificate['input_word']!=word
            or certificate['source_sha256']!=expected['algorithm_source_sha256']
            or certificate['catalog_sha256']!=expected['catalog_sha256']
            or certificate['limits']!=limits or certificate['trace_enabled'] is not False
            or certificate['coefficient_ring']!=expected['coefficient_ring']
            or certificate['execution']['batch_size']!=expected['batch_size']
            or certificate['status'] not in ('complete','incomplete','control_obligation')):
        raise ValueError('certificate identity does not match literal input/configuration')
    if expected.get('execution')=='native':
        from zipper_algebra.zipping_algorithm.native import _strict_equal
        if (not _strict_equal(certificate.get('n'),n)
                or not _strict_equal(certificate.get('input_word'),word)
                or not _strict_equal(certificate.get('limits'),limits)
                or not _strict_equal(certificate['execution'].get('batch_size'),expected['batch_size'])
                or not _strict_equal(certificate.get('source_sha256'),expected['algorithm_source_sha256'])
                or not _strict_equal(certificate.get('native_build'),expected['native_build'])
                or not _strict_equal(certificate.get('python_reference_source_sha256'),
                                     expected['python_reference_source_sha256'])):
            raise ValueError('native certificate build or reference/input identity differs')


def _slider_size(certificate):
    return max((component['size'] for term in certificate['output']
                for component in term['components']),default=0)


def _cache_info(session):
    return session.cache_info() if hasattr(session,'cache_info') else {}


def _numeric_delta(before, after):
    result = {}
    for key,value in after.items():
        previous=before.get(key,{} if isinstance(value,dict) else 0)
        if isinstance(value,dict) and isinstance(previous,dict):
            result[key]=_numeric_delta(previous,value)
        elif type(value) in (int,float) and type(previous) in (int,float):
            result[key]=value-previous
    return result


def _memory_info():
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak=int(peak if sys.platform=='darwin' else peak*1024)
    try:
        current=int(Path('/proc/self/statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
    except (OSError,ValueError,IndexError):
        current=None
    return {'process_peak_rss_bytes':peak,'process_current_rss_bytes':current}


def _native_memory_info(pid):
    try:
        fields=dict(line.split(':',1) for line in Path(f'/proc/{pid}/status').read_text().splitlines() if ':' in line)
        return {name:int(fields[key].split()[0])*1024 for key,name in
                (('VmRSS','native_process_current_rss_bytes'),('VmHWM','native_process_peak_rss_bytes'))}
    except (OSError,ValueError,KeyError,IndexError):
        return {'native_process_current_rss_bytes':None,'native_process_peak_rss_bytes':None}


def _compute(n,word):
    global _worker_sequence
    _worker_sequence+=1
    sequence=_worker_sequence
    assert _session is not None and _worker_limits is not None and _worker_identity is not None
    before=_cache_info(_session)
    verifier_before=_cache_info(_verifier)
    start=time.perf_counter_ns()
    certificate=_session.reduce(word,n=n,limits=_worker_limits,trace=False,
                                workers=1,batch_size=_worker_identity['batch_size'])
    producer=time.perf_counter_ns()-start
    after_reduce=_cache_info(_session)
    encoded=_json(certificate)
    start=time.perf_counter_ns()
    reloaded=_loads(encoded)
    if _worker_identity.get('verifier_backend')=='cpp':accepted=_verifier.check(reloaded)
    else:accepted=(_verifier.check_native(_session,reloaded) if _verifier is not None else _session.check(reloaded))
    if accepted is not True:
        raise ValueError('certificate failed exact replay')
    checker=time.perf_counter_ns()-start
    after_check=_cache_info(_session)
    _certificate_identity(reloaded,n,word,asdict(_worker_limits),_worker_identity)
    telemetry={'schema':'zipper-worker-telemetry-v1','worker_pid':os.getpid(),'worker_sequence':sequence,'before':before,'after_reduce':after_reduce,
               'after_check':after_check,'producer_delta':_numeric_delta(before,after_reduce),
               'checker_producer_cache_delta':_numeric_delta(after_reduce,after_check),
               'checker_evaluation_counters':None,
               'checker_mode':'independent-cold-replay' if _worker_identity.get('execution')=='memoized' else 'collected-serial-replay',
               'evaluation_scope':'producer-attempted-transitions-including-rolled-back-batches',**_memory_info()}
    if _verifier is not None:
        verifier_after=_cache_info(_verifier)
        telemetry.update(checker_mode=('independent-cpp-cached-replay' if _worker_identity.get('verifier_backend')=='cpp' else 'independent-python-cached-replay'),
                         verifier_before=verifier_before,verifier_after=verifier_after,
                         verifier_delta=_numeric_delta(verifier_before,verifier_after),
                         verifier_source_sha256=_worker_identity['verifier_source_sha256'],
                         native_pid=_session.pid,native_runtime=_session.last_runtime_info,
                         **_native_memory_info(_session.pid))
        if _worker_identity.get('verifier_backend')=='cpp':
            telemetry.update(verifier_pid=_verifier.pid,verifier_build=_verifier_build(_verifier.metadata),
                             verifier_allocator=_verifier.metadata.get('allocator_library'),
                             verifier_process_memory=_native_memory_info(_verifier.pid))
    return encoded,producer,checker,telemetry


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


def _has_table(db,name):
    return db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",(name,)).fetchone() is not None


def _load_versions(db):
    metadata=dict(db.execute('SELECT key,value FROM metadata'))
    legacy_raw=metadata['identity'];legacy=_loads(legacy_raw)
    if _json(legacy)!=legacy_raw:
        raise ValueError('legacy identity is not canonical JSON')
    legacy_sha=_sha(legacy_raw)
    if _has_table(db,'versions'):
        versions={}
        for digest,encoded,_ in db.execute('SELECT hash,identity,registered_utc FROM versions'):
            if _sha(encoded)!=digest:raise ValueError('version identity hash mismatch')
            versions[digest]=_loads(encoded)
        if versions.get(legacy_sha)!=legacy or metadata.get('active_version') not in versions:
            raise ValueError('version history does not retain its legacy/active identity')
        if _has_table(db,'migration_events'):
            for old,new,snapshot in db.execute('SELECT old_version,new_version,metadata_snapshot FROM migration_events'):
                saved=_loads(snapshot)
                if old not in versions or new not in versions or saved.get('identity')!=legacy_raw or saved.get('active_version',legacy_sha)!=old:
                    raise ValueError('migration history version/metadata binding mismatch')
        return versions,metadata['active_version']
    return {legacy_sha:legacy},legacy_sha


def _extended_schema(db):
    db.execute('CREATE TABLE IF NOT EXISTS versions(hash TEXT PRIMARY KEY,identity TEXT NOT NULL,registered_utc TEXT NOT NULL)')
    db.execute('CREATE TABLE IF NOT EXISTS migration_events(id INTEGER PRIMARY KEY,old_version TEXT,new_version TEXT,backup_path TEXT,backup_sha TEXT,metadata_snapshot TEXT,audited_words INTEGER,registered_utc TEXT)')
    db.execute('CREATE TABLE IF NOT EXISTS word_attempts AS SELECT *,CAST(NULL AS TEXT) AS runtime,CAST(NULL AS TEXT) AS runtime_sha,CAST(NULL AS TEXT) AS archived_utc FROM words WHERE 0')
    db.execute('CREATE UNIQUE INDEX IF NOT EXISTS word_attempt_identity ON word_attempts(n,length,ordinal,certificate_sha,config_sha)')
    db.execute('CREATE TABLE IF NOT EXISTS word_runtime(n INTEGER,length INTEGER,ordinal INTEGER,config_sha TEXT,runtime TEXT,runtime_sha TEXT,PRIMARY KEY(n,length,ordinal))')
    db.execute('CREATE TABLE IF NOT EXISTS worker_state(worker_pid INTEGER,config_sha TEXT,runtime TEXT,runtime_sha TEXT,updated_utc TEXT,PRIMARY KEY(worker_pid,config_sha))')
    db.execute('CREATE TABLE IF NOT EXISTS prefix_rule_totals(version TEXT,n INTEGER,rule TEXT,count INTEGER,PRIMARY KEY(version,n,rule))')
    db.execute('CREATE TABLE IF NOT EXISTS prefix_evaluation_totals(version TEXT,n INTEGER,kind TEXT,rule TEXT,count INTEGER,PRIMARY KEY(version,n,kind,rule))')


def _configuration_map(db,versions):
    configs={}
    for digest,encoded in db.execute('SELECT hash,config FROM configurations'):
        config=_loads(encoded)
        if not isinstance(config,dict) or _sha(encoded)!=digest or config.get('identity_sha256') not in versions:
            raise ValueError('saved configuration hash/version binding mismatch')
        expected=versions[config['identity_sha256']]
        if expected.get('execution')=='native' and (config.get('execution')!='native'
                or config.get('verifier_source_sha256')!=expected['verifier_source_sha256']
                or config.get('verifier_backend','python')!=expected.get('verifier_backend','python')
                or config.get('verifier_build')!=expected.get('verifier_build')
                or config.get('verifier_allocator')!=expected.get('verifier_allocator')):
            raise ValueError('native configuration verifier identity differs from epoch')
        configs[digest]=config
    return configs


def _validate_record(row,configs,versions):
    n,length,ordinal,word_json,input_sha,status,diagnostic,metrics,blob,digest,size,codec,config_sha,producer,checker,slider,check=row
    if (type(n) is not int or n<2 or type(length) is not int or length<0
            or type(ordinal) is not int or ordinal<0 or ordinal>=word_count(n,length)):
        raise ValueError('invalid literal row index')
    word=list(_word_at(n,length,ordinal))
    if word_json!=_json(word) or input_sha!=_sha(_json({'n':n,'word':word})):
        raise ValueError('saved literal input binding mismatch')
    config=configs.get(config_sha)
    if config is None:raise ValueError('saved configuration is missing')
    version=config['identity_sha256'];expected=versions[version]
    if str(n) not in config['targets'] or length>config['targets'][str(n)]:
        raise ValueError('saved configuration rank/length binding mismatch')
    if codec!=CODEC or codec!=expected['codec'] or check!=expected['production_check'] or type(size) is not int or size<0:
        raise ValueError('saved certificate codec/check provenance mismatch')
    try:
        raw=zlib.decompress(blob)
        if len(raw)!=size or hashlib.sha256(raw).hexdigest()!=digest:
            raise ValueError('saved certificate length/hash mismatch')
        certificate=_loads(raw)
    except (zlib.error,UnicodeError,json.JSONDecodeError,TypeError) as error:
        raise ValueError('saved compressed certificate is invalid') from error
    _certificate_identity(certificate,n,word,config['limits'],expected)
    if (status!=certificate['status'] or diagnostic!=_json(certificate['diagnostic'])
            or metrics!=_json(certificate['metrics']) or slider!=_slider_size(certificate)
            or type(producer) is not int or producer<0 or type(checker) is not int or checker<0):
        raise ValueError('saved certificate columns differ from its contents')
    return certificate,version


def _validate_runtime(encoded,digest,config=None):
    if not isinstance(encoded,str) or not isinstance(digest,str) or _sha(encoded)!=digest:raise ValueError('worker telemetry hash mismatch')
    value=_loads(encoded)
    if type(value.get('worker_pid')) is not int or value['worker_pid']<1:
        raise ValueError('worker telemetry PID is invalid')
    if value.get('schema')=='zipper-worker-telemetry-v1':
        if value['producer_delta']!=_numeric_delta(value['before'],value['after_reduce']):
            raise ValueError('worker producer telemetry delta mismatch')
        if value['checker_producer_cache_delta']!=_numeric_delta(value['after_reduce'],value['after_check']):
            raise ValueError('worker checker telemetry delta mismatch')
    if config and config.get('execution')=='memoized':
        for stage in ('before','after_reduce','after_check'):
            info=value[stage];rule=info['rule_cache']
            if (info['max_bytes']!=config['cache_bytes'] or info['max_entries']!=config['cache_entries']
                    or info['bytes']>config['cache_bytes'] or info['entries']>config['cache_entries']
                    or rule['max_bytes']!=config['rule_cache_bytes'] or rule['max_entries']!=config['rule_cache_entries']
                    or rule['retained_bytes']>config['rule_cache_bytes'] or rule['entries']>config['rule_cache_entries']):
                raise ValueError('worker cache telemetry exceeds or differs from configured budgets')
    if config and config.get('execution')=='native':
        for stage in ('before','after_reduce','after_check'):
            info=value[stage]
            if (info['max_bytes']!=config['cache_bytes'] or info['max_entries']!=config['cache_entries']
                    or info['bytes']>config['cache_bytes'] or info['entries']>config['cache_entries']):
                raise ValueError('native cache telemetry differs from configured budgets')
        backend=config.get('verifier_backend','python')
        mode='independent-cpp-cached-replay' if backend=='cpp' else 'independent-python-cached-replay'
        if (backend not in ('python','cpp') or value.get('checker_mode')!=mode
                or value.get('verifier_source_sha256')!=config['verifier_source_sha256']
                or value['verifier_delta']!=_numeric_delta(value['verifier_before'],value['verifier_after'])
                or type(value.get('native_pid')) is not int or value['native_pid']<1):
            raise ValueError('independent verifier telemetry identity or delta differs')
        if backend=='cpp' and (type(value.get('verifier_pid')) is not int or value['verifier_pid']<1
                or value.get('verifier_build')!=config['verifier_build']):
            raise ValueError('independent C++ verifier build or PID differs')
        if backend=='cpp':
            allocator=value.get('verifier_allocator')
            if (allocator!=_allocator_runtime(config.get('verifier_allocator'))
                    or allocator is not None and allocator.get('mapped_verified') is not True):
                raise ValueError('independent C++ verifier allocator identity differs')
        for stage in ('verifier_before','verifier_after'):
            info=value[stage]
            if (info['max_bytes']!=config['verifier_cache_bytes']
                    or info['max_entries']!=config['verifier_cache_entries']
                    or info['bytes']>config['verifier_cache_bytes']
                    or info['entries']>config['verifier_cache_entries']):
                raise ValueError('independent verifier cache differs from configured budgets')
            if backend=='python':
                rule=info['rule_cache']
                if (rule['max_bytes']!=config['rule_cache_bytes']
                    or rule['max_entries']!=config['rule_cache_entries']
                    or rule['retained_bytes']>config['rule_cache_bytes']
                    or rule['entries']>config['rule_cache_entries']):
                    raise ValueError('independent verifier rule cache differs from configured budgets')
    return value


def _evaluation_counts(telemetry):
    for kind in ('actual_unit_evaluations','unit_lru_hits','unit_batch_reuses'):
        for rule,count in telemetry.get('producer_delta',{}).get(kind,{}).items():
            if type(count) is not int or count<0:raise ValueError('invalid actual evaluation counter')
            if count:yield kind,rule,count


def audit(db,expected=None,*,replay=False):
    """Validate all source epochs; exact replay requires matching loaded source."""
    versions,active=_load_versions(db)
    configs=_configuration_map(db,versions)
    cursors={n:[length,ordinal] for n,length,ordinal in db.execute('SELECT n,length,ordinal FROM ranks')}
    if any(type(n) is not int or n<2 or any(type(v) is not int or v<0 for v in cursor) for n,cursor in cursors.items()):
        raise ValueError('invalid saved rank cursor')
    next_cursor={n:[0,0] for n in cursors};blocked=set();count=0;history_count=0
    histograms: Counter[tuple[str,int,str]] = Counter()
    evaluations: Counter[tuple[str,int,str,str]] = Counter()
    sessions={}
    def replay_one(certificate,version):
        wanted=versions[version]
        if expected is None or wanted['algorithm_source_sha256']!=expected['algorithm_source_sha256']:
            raise ValueError('historical replay requires its matching source checkout; hash audit remains available')
        if version not in sessions:sessions[version]=_new_session(wanted)
        if sessions[version].check(certificate) is not True:raise ValueError('saved certificate failed exact replay')
    for row in db.execute('SELECT * FROM words ORDER BY n,length,ordinal'):
        n,length,ordinal=row[:3]
        if n not in cursors or n in blocked or [length,ordinal]!=next_cursor[n]:
            raise ValueError('saved rows are not consecutive literal rank prefixes')
        certificate,version=_validate_record(row,configs,versions)
        for rule,value in certificate['metrics']['rule_applications'].items():
            if type(value) is not int or value<0:raise ValueError('invalid logical rule count')
            if value:histograms[version,n,rule]+=value
        runtime=db.execute('SELECT config_sha,runtime,runtime_sha FROM word_runtime WHERE n=? AND length=? AND ordinal=?',(n,length,ordinal)).fetchone() if _has_table(db,'word_runtime') else None
        if runtime:
            if runtime[0]!=row[12]:raise ValueError('worker telemetry configuration mismatch')
            telemetry=_validate_runtime(runtime[1],runtime[2],configs[row[12]])
            for kind,rule,value in _evaluation_counts(telemetry):evaluations[version,n,kind,rule]+=value
        elif versions[version].get('schema')=='zipper-versioned-literal-enumeration-v2':
            raise ValueError('versioned word is missing its worker telemetry')
        if replay:replay_one(certificate,version)
        count+=1
        if certificate['status']=='complete':next_cursor[n]=_advance(n,next_cursor[n])
        else:blocked.add(n)
    if next_cursor!=cursors:raise ValueError('saved rank cursor differs from its literal prefix')
    if _has_table(db,'word_attempts'):
        for archived in db.execute('SELECT * FROM word_attempts'):
            certificate,version=_validate_record(archived[:17],configs,versions)
            if certificate['status']=='complete':raise ValueError('completed prefix rows must not be archived/replaced')
            if archived[17] is not None:_validate_runtime(archived[17],archived[18],configs[archived[12]])
            elif versions[version].get('schema')=='zipper-versioned-literal-enumeration-v2':
                raise ValueError('archived attempt is missing its telemetry')
            if replay:replay_one(certificate,version)
            history_count+=1
    if _has_table(db,'prefix_rule_totals'):
        stored={(v,n,rule):value for v,n,rule,value in db.execute('SELECT * FROM prefix_rule_totals')}
        if stored!=dict(histograms):raise ValueError('logical rule histogram differs from current prefix certificates')
    if _has_table(db,'prefix_evaluation_totals'):
        stored_evaluations={(v,n,kind,rule):value for v,n,kind,rule,value in db.execute('SELECT * FROM prefix_evaluation_totals')}
        if stored_evaluations!=dict(evaluations):raise ValueError('actual evaluation histogram differs from saved producer telemetry')
    if _has_table(db,'worker_state'):
        for pid,config_sha,encoded,digest,_ in db.execute('SELECT * FROM worker_state'):
            if config_sha not in configs or _validate_runtime(encoded,digest,configs[config_sha])['worker_pid']!=pid:
                raise ValueError('worker state configuration/PID mismatch')
    return {'words':count,'archived_attempts':history_count,
            'hash_and_binding_checked':count+history_count,'replayed':count+history_count if replay else 0,
            'active_version':active,'versions':list(versions),
            'evaluation_counts_by_version':[{'version':v,'n':n,'kind':kind,'rule':rule,'count':value}
                                            for (v,n,kind,rule),value in sorted(evaluations.items())],
            'rule_counts_by_version':[{'version':v,'n':n,'rule':rule,'count':value}
                                      for (v,n,rule),value in sorted(histograms.items())]}


def _file_sha(path):
    return _hash(path)


def register_version(db,expected,output,*,migrate=False,backup_path=None):
    """Audit and back up before atomically adding an epoch; never rewrite identities."""
    versions,active=_load_versions(db);wanted=_sha(_json(expected))
    if wanted==active:return None
    if not migrate:raise ValueError('enumeration identity changed; explicit --migrate and backup are required')
    for field in ('catalog_sha256','catalog_schema_sha256','coefficient_ring','codec'):
        if versions[active][field]!=expected[field]:raise ValueError('migration cannot change the catalog, algebra or storage codec')
    verified=audit(db)
    metadata=dict(db.execute('SELECT key,value FROM metadata'))
    now=datetime.now(timezone.utc)
    backup=Path(backup_path) if backup_path else Path(output)/'backups'/f"before-{active[:12]}-{now.strftime('%Y%m%dT%H%M%S%f')}.sqlite3"
    backup.parent.mkdir(parents=True,exist_ok=True)
    with backup.open('xb'):
        pass
    destination=sqlite3.connect(backup)
    try:
        db.backup(destination)
        if destination.execute('PRAGMA integrity_check').fetchall()!=[('ok',)]:raise ValueError('migration backup integrity failed')
        destination.execute('PRAGMA wal_checkpoint(TRUNCATE)')
    finally:destination.close()
    backup_sha=_file_sha(backup);stamp=now.isoformat()
    _assert_source(expected)
    with db:
        db.execute('BEGIN IMMEDIATE')
        _extended_schema(db)
        for key,value in versions.items():
            db.execute('INSERT OR IGNORE INTO versions VALUES(?,?,?)',(key,_json(value),stamp))
        db.execute('INSERT OR IGNORE INTO versions VALUES(?,?,?)',(wanted,_json(expected),stamp))
        db.execute("INSERT OR REPLACE INTO metadata VALUES('active_version',?)",(wanted,))
        db.execute('DELETE FROM prefix_evaluation_totals')
        db.executemany('INSERT INTO prefix_evaluation_totals VALUES(?,?,?,?,?)',
                       [(x['version'],x['n'],x['kind'],x['rule'],x['count']) for x in verified['evaluation_counts_by_version']])
        db.execute('DELETE FROM prefix_rule_totals')
        db.executemany('INSERT INTO prefix_rule_totals VALUES(?,?,?,?)',
                       [(x['version'],x['n'],x['rule'],x['count']) for x in verified['rule_counts_by_version']])
        db.execute('INSERT INTO migration_events(old_version,new_version,backup_path,backup_sha,metadata_snapshot,audited_words,registered_utc) VALUES(?,?,?,?,?,?,?)',
                   (active,wanted,str(backup.resolve()),backup_sha,_json(metadata),verified['words'],stamp))
    return {'old_version':active,'new_version':wanted,'backup':str(backup.resolve()),
            'backup_sha256':backup_sha,'audited_words':verified['words'],'audit':verified}


def _update_histogram(db,version,n,certificate,sign):
    for rule,value in certificate['metrics']['rule_applications'].items():
        db.execute('INSERT INTO prefix_rule_totals VALUES(?,?,?,?) ON CONFLICT(version,n,rule) DO UPDATE SET count=count+excluded.count',
                   (version,n,rule,sign*value))
    db.execute('DELETE FROM prefix_rule_totals WHERE count=0')


def _update_evaluations(db,version,n,telemetry,sign):
    for kind,rule,value in _evaluation_counts(telemetry):
        db.execute('INSERT INTO prefix_evaluation_totals VALUES(?,?,?,?,?) ON CONFLICT(version,n,kind,rule) DO UPDATE SET count=count+excluded.count',
                   (version,n,kind,rule,sign*value))
    db.execute('DELETE FROM prefix_evaluation_totals WHERE count=0')


def _commit_word(db,n,cursor,word,encoded,producer,checker,telemetry,config_sha,expected):
    config=_loads(db.execute('SELECT config FROM configurations WHERE hash=?',(config_sha,)).fetchone()[0])
    certificate=_loads(encoded);_certificate_identity(certificate,n,word,config['limits'],expected)
    status=certificate['status'];following=_advance(n,cursor) if status=='complete' else cursor
    raw=encoded.encode();runtime=_json(telemetry);runtime_sha=_sha(runtime);version=_sha(_json(expected));stamp=datetime.now(timezone.utc).isoformat()
    _validate_runtime(runtime,runtime_sha,config)
    with db:
        previous=db.execute('SELECT * FROM words WHERE n=? AND length=? AND ordinal=?',(n,*cursor)).fetchone()
        if previous:
            if previous[5]=='complete':raise ValueError('completed prefix words are immutable')
            previous_runtime=db.execute('SELECT runtime,runtime_sha FROM word_runtime WHERE n=? AND length=? AND ordinal=?',(n,*cursor)).fetchone()
            db.execute('INSERT OR IGNORE INTO word_attempts VALUES('+','.join('?' for _ in range(20))+')',
                       (*previous,*(previous_runtime or (None,None)),stamp))
            old_certificate=_loads(zlib.decompress(previous[8]))
            old_version=_loads(db.execute('SELECT config FROM configurations WHERE hash=?',(previous[12],)).fetchone()[0])['identity_sha256']
            _update_histogram(db,old_version,n,old_certificate,-1)
            if previous_runtime:_update_evaluations(db,old_version,n,_loads(previous_runtime[0]),-1)
        db.execute('INSERT OR REPLACE INTO words VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                   (n,*cursor,_json(word),_sha(_json({'n':n,'word':word})),status,
                    _json(certificate['diagnostic']),_json(certificate['metrics']),zlib.compress(raw),_sha(encoded),len(raw),CODEC,
                    config_sha,producer,checker,_slider_size(certificate),expected['production_check']))
        db.execute('INSERT OR REPLACE INTO word_runtime VALUES(?,?,?,?,?,?)',(n,*cursor,config_sha,runtime,runtime_sha))
        prior=db.execute('SELECT runtime FROM worker_state WHERE worker_pid=? AND config_sha=?',(telemetry['worker_pid'],config_sha)).fetchone()
        if prior is None or _loads(prior[0]).get('worker_sequence',0)<telemetry['worker_sequence']:
            db.execute('INSERT OR REPLACE INTO worker_state VALUES(?,?,?,?,?)',(telemetry['worker_pid'],config_sha,runtime,runtime_sha,stamp))
        _update_histogram(db,version,n,certificate,1)
        _update_evaluations(db,version,n,telemetry,1)
        db.execute('UPDATE ranks SET length=?,ordinal=? WHERE n=?',(*following,n))
    return status,following


def summarize(db,targets,reasons):
    result=[];active_version=_load_versions(db)[1]
    for n,length,ordinal in db.execute('SELECT n,length,ordinal FROM ranks ORDER BY n'):
        layers=[];through=-1
        for L,count,complete in db.execute("SELECT length,count(*),sum(status='complete') FROM words WHERE n=? GROUP BY length ORDER BY length",(n,)):
            expected=word_count(n,L);finished=count==complete==expected
            if finished and L==through+1:through=L
            layers.append({'length':L,'expected_words':expected,'visited_words':count,
                           'complete_words':complete,'layer_complete':finished})
        totals=db.execute("SELECT count(*),sum(status='complete'),max(output_slider_size),sum(length(certificate)),sum(certificate_bytes) FROM words WHERE n=?",(n,)).fetchone()
        unresolved=db.execute("SELECT word,status,diagnostic FROM words WHERE n=? AND length=? AND ordinal=?",(n,length,ordinal)).fetchone()
        version_counts=[{'version':v,'words':count,'complete_words':complete} for v,count,complete in db.execute(
            "SELECT json_extract(c.config,'$.identity_sha256'),count(*),sum(w.status='complete') FROM words w JOIN configurations c ON c.hash=w.config_sha WHERE w.n=? GROUP BY json_extract(c.config,'$.identity_sha256')",(n,))]
        active_counts=next((x for x in version_counts if x['version']==active_version),
                           {'version':active_version,'words':0,'complete_words':0})
        result.append({'n':n,'target_length':targets.get(n),'complete_through_length':through,
            'coverage_scope':'all_registered_versions','active_version_counts':active_counts,
            'visited_words':totals[0],'complete_words':totals[1] or 0,'maximum_output_slider_size':totals[2] or 0,
            'compressed_certificate_bytes':totals[3] or 0,'raw_certificate_bytes':totals[4] or 0,
            'next_cursor':{'length':length,'ordinal':ordinal},'stop_reason':reasons.get(n),
            'first_unresolved':None if not unresolved else {'word':_loads(unresolved[0]),'status':unresolved[1],'diagnostic':_loads(unresolved[2])},'layers':layers,
            'version_counts':version_counts,
            'actual_evaluation_counts_by_version':[{'version':v,'kind':kind,'rule':rule,'count':count} for v,kind,rule,count in db.execute(
                'SELECT version,kind,rule,count FROM prefix_evaluation_totals WHERE n=? ORDER BY version,kind,rule',(n,))] if _has_table(db,'prefix_evaluation_totals') else None,
            'rule_counts_by_version':[{'version':v,'rule':rule,'count':count} for v,rule,count in db.execute(
                'SELECT version,rule,count FROM prefix_rule_totals WHERE n=? ORDER BY version,rule',(n,))] if _has_table(db,'prefix_rule_totals') else None})
    return result


def _read_only_audit(output,*,expected=None,replay=False):
    database=Path(output).resolve()/'words.sqlite3'
    if not database.is_file():raise ValueError('audit-only requires an existing database')
    with sqlite3.connect(database.as_uri()+'?mode=ro',uri=True) as db:
        db.execute('BEGIN')
        if replay:_preflight_replay(db,expected)
        result=audit(db,expected,replay=replay)
        return {'schema':'zipper-versioned-prefix-audit-v2','database':str(database),
                'status':'audit_only','audit':result,'mutated_database':False}


def _preflight_replay(db,expected):
    versions,_=_load_versions(db)
    if expected is None or any(value['algorithm_source_sha256']!=expected['algorithm_source_sha256'] for value in versions.values()):
        raise ValueError('historical replay requires matching source checkouts; no migration was performed')


def _empty_bootstrap(db):
    """Only a recognized empty initialization may acquire its first identity."""
    known={'metadata','configurations','ranks','words','versions','migration_events',
           'word_attempts','word_runtime','worker_state','prefix_rule_totals','prefix_evaluation_totals'}
    tables=[row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")]
    if any(name not in known for name in tables):raise ValueError('unbound database has an unknown schema')
    for table in tables:
        if db.execute(f'SELECT count(*) FROM "{table}"').fetchone()[0]:
            raise ValueError('unbound database contains saved records; refusing to assign an identity')
    columns={'metadata':['key','value'],'configurations':['hash','config'],'ranks':['n','length','ordinal'],
             'words':['n','length','ordinal','word','input_sha','status','diagnostic','metrics','certificate',
                      'certificate_sha','certificate_bytes','codec','config_sha','producer_ns','checker_ns',
                      'output_slider_size','production_check']}
    for table,wanted in columns.items():
        if table in tables and [row[1] for row in db.execute(f'PRAGMA table_info("{table}")')]!=wanted:
            raise ValueError('unbound database has an incompatible table layout')


def _start_store(db,expected,output,*,migrate=False,backup_path=None):
    has_identity=_has_table(db,'metadata') and db.execute("SELECT 1 FROM metadata WHERE key='identity'").fetchone()
    if not has_identity:
        _empty_bootstrap(db);_schema(db)
        stamp=datetime.now(timezone.utc).isoformat();version=_sha(_json(expected))
        with db:
            db.execute('BEGIN IMMEDIATE');_extended_schema(db)
            db.execute("INSERT INTO metadata VALUES('identity',?)",(_json(expected),))
            db.execute('INSERT INTO versions VALUES(?,?,?)',(version,_json(expected),stamp))
            db.execute("INSERT INTO metadata VALUES('active_version',?)",(version,))
        return None
    return register_version(db,expected,output,migrate=migrate,backup_path=backup_path)


def _goals_complete(db,targets):
    cursors={n:L for n,L in db.execute('SELECT n,length FROM ranks')}
    return all(n in cursors and cursors[n]>target for n,target in targets.items())


def _install_drain_handler(requested):
    """SIGUSR1 stops admissions and drains checked prefixes; restore on exit."""
    if threading.current_thread() is not threading.main_thread():return None
    previous=signal.getsignal(signal.SIGUSR1)
    def request_drain(signum,frame):
        requested[0]=True
    signal.signal(signal.SIGUSR1,request_drain)
    return previous


def _restore_drain_handler(previous):
    if previous is not None:signal.signal(signal.SIGUSR1,previous)


def exhaust(output,*,targets=None,workers=8,batch_size=64,wall_seconds=86400,limits=None,
            max_new_words=None,replay_audit=False,audit_only=False,migrate=False,
            migration_backup=None,execution='native',verifier='python',verifier_allocator=None,cache_bytes=256*1024*1024,
            cache_entries=262144,rule_cache_bytes=4*1024*1024,rule_cache_entries=256,
            verifier_cache_bytes=256*1024*1024,verifier_cache_entries=262144,
            admission_deadline_utc=None,_executor_factory=_process_pool):
    targets=DEFAULT_TARGETS.copy() if targets is None else dict(targets)
    if not targets or any(type(n) is not int or n<2 or type(L) is not int or L<0 for n,L in targets.items()):
        raise ValueError('integer rank >=2 and nonnegative target length required')
    if type(workers) is not int or workers<1 or type(batch_size) is not int or batch_size<1:
        raise ValueError('positive integer workers and batch_size required')
    if execution not in ('native','memoized','collected'):raise ValueError('execution must be native, memoized or collected')
    if any(type(v) is not int or v<0 for v in (cache_bytes,cache_entries,rule_cache_bytes,rule_cache_entries,verifier_cache_bytes,verifier_cache_entries)):
        raise ValueError('cache budgets must be nonnegative integers')
    if type(wall_seconds) not in (int,float) or not math.isfinite(wall_seconds) or wall_seconds<=0:
        raise ValueError('positive finite wall_seconds required')
    if max_new_words is not None and (type(max_new_words) is not int or max_new_words<0):
        raise ValueError('nonnegative max_new_words required')
    if verifier not in ('python','cpp') or (execution!='native' and verifier!='python'):
        raise ValueError('C++ verification requires native execution; verifier must be python or cpp')
    if verifier_allocator is not None and (execution!='native' or verifier!='cpp'):
        raise ValueError('verifier allocator requires native execution with C++ verification')
    deadline=None
    if admission_deadline_utc is not None:
        try:
            if not isinstance(admission_deadline_utc,str):raise ValueError
            deadline=datetime.fromisoformat(admission_deadline_utc)
            if deadline.utcoffset() is None or deadline.utcoffset().total_seconds()!=0:
                raise ValueError
        except (TypeError,ValueError) as error:
            raise ValueError('admission deadline must be a timezone-aware UTC timestamp') from error
        admission_deadline_utc=deadline.isoformat()
    limits=Limits(max_steps=200000,max_terms=100000,expansion_budget=1000000) if limits is None else limits
    if not isinstance(limits,Limits):raise ValueError('limits must be a Limits record')
    if audit_only:
        if migrate:raise ValueError('audit-only cannot migrate a database')
        expected=(identity(batch_size,execution) if verifier=='python' else identity(batch_size,execution,verifier,verifier_allocator)) if replay_audit else None
        return _read_only_audit(output,expected=expected,replay=replay_audit)
    expected=(identity(batch_size,execution) if verifier=='python' else identity(batch_size,execution,verifier,verifier_allocator))
    output=Path(output);output.mkdir(parents=True,exist_ok=True);database=output/'words.sqlite3'
    config={'targets':targets,'workers':workers,'max_in_flight':2*workers,
            'per_rank_window':'ceil_global_window_over_active_ranks','limits':asdict(limits),
            'wall_seconds':wall_seconds,'admission_deadline_utc':admission_deadline_utc,
            'max_new_words':max_new_words,'identity_sha256':_sha(_json(expected)),
            'execution':execution,'cache_bytes':cache_bytes,'cache_entries':cache_entries,
            'rule_cache_bytes':rule_cache_bytes,'rule_cache_entries':rule_cache_entries,
            'verifier_cache_bytes':verifier_cache_bytes,'verifier_cache_entries':verifier_cache_entries,
            'python':sys.version,'created_utc':datetime.now(timezone.utc).isoformat()}
    if execution=='native':
        config['verifier_source_sha256']=expected['verifier_source_sha256']
        if verifier=='cpp':
            config.update(verifier_backend='cpp',verifier_build=expected['verifier_build'])
            if 'verifier_allocator' in expected:config['verifier_allocator']=expected['verifier_allocator']
    config_sha=_sha(_json(config));lock=(output/'.parallel.lock').open('a+')
    try:fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError as error:
        lock.close();raise ValueError('another writer owns this campaign') from error
    db=sqlite3.connect(database);pool=None;new=0;peak=0;run_started=False
    queues: dict[int, deque[Any]] = {};discarded: set[Any] = set()
    drain_requested=[False];previous_handler=_install_drain_handler(drain_requested)
    try:
        audit_start=time.perf_counter();pre_replayed=None
        if replay_audit and _has_table(db,'metadata') and db.execute("SELECT 1 FROM metadata WHERE key='identity'").fetchone():
            _preflight_replay(db,expected)
            pre_replayed=audit(db,expected,replay=True)
        migration=_start_store(db,expected,output,migrate=migrate,backup_path=migration_backup)
        if pre_replayed is not None:
            audited=pre_replayed.copy()
            audited['active_version']=_sha(_json(expected))
            audited['versions']=list(_load_versions(db)[0])
        elif migration:
            audited=migration['audit'].copy()
            audited['active_version']=_sha(_json(expected))
            audited['versions']=list(_load_versions(db)[0])
        else:audited=audit(db,expected,replay=replay_audit)
        audit_seconds=time.perf_counter()-audit_start
        with db:
            db.execute('INSERT OR IGNORE INTO configurations VALUES(?,?)',(config_sha,_json(config)))
            db.execute("INSERT OR REPLACE INTO metadata VALUES('last_config',?)",(_json(config),))
            for n in targets:db.execute('INSERT OR IGNORE INTO ranks VALUES(?,0,0)',(n,))
        cursors={n:[L,o] for n,L,o in db.execute('SELECT n,length,ordinal FROM ranks') if n in targets}
        next_submit={n:cursor[:] for n,cursor in cursors.items()};queues={n:deque() for n in targets}
        reasons={};active=set(targets)
        for n,cursor in cursors.items():
            old=db.execute('SELECT status,certificate,config_sha FROM words WHERE n=? AND length=? AND ordinal=?',(n,*cursor)).fetchone()
            if old:
                before=_loads(zlib.decompress(old[1]))['limits'];after=asdict(limits)
                old_version=_loads(db.execute('SELECT config FROM configurations WHERE hash=?',(old[2],)).fetchone()[0])['identity_sha256']
                # A new source epoch may retry the old checkpoint even with equal limits.
                if old_version==_sha(_json(expected)):
                    if old[0]=='control_obligation':reasons[n]='control_obligation'
                    elif before==after or any(after[k]<v for k,v in before.items()):reasons[n]='word_limit_needs_larger_limits'
            if cursor[0]>targets[n]:reasons[n]='target_length'
            if n in reasons:active.discard(n)
        started=time.perf_counter();last_snapshot=started;rotation=deque(sorted(targets))
        deadline_reached=deadline is not None and datetime.now(timezone.utc)>=deadline
        admitting=bool(active) and max_new_words!=0 and not drain_requested[0] and not deadline_reached
        stop=('signal_drain' if drain_requested[0] else
              'admission_deadline' if deadline_reached else
              'running' if admitting else ('word_budget' if max_new_words==0 else 'rank_limits'))
        peak_per_rank={n:0 for n in targets}
        def snapshot(state='running',phase='computing'):
            stamp=datetime.now(timezone.utc).isoformat()
            run_state={'status':state,'phase':phase,'stop_reason':stop,'pid':os.getpid(),
                       'config_sha':config_sha,'recorded_utc':stamp}
            with db:db.execute("INSERT OR REPLACE INTO metadata VALUES('run_state',?)",(_json(run_state),))
            versions,active_version=_load_versions(db)
            states=[_loads(row[0]) for row in db.execute('SELECT runtime FROM worker_state WHERE config_sha=?',(config_sha,))]
            data={'schema':expected['schema'],'identity':expected,'active_version':active_version,
                  'versions':versions,'legacy_identity':_loads(db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0]),
                  'identity_scope':'active_version','coverage_scope':'all_registered_versions',
                  'database':str(database.resolve()),'status':state,'phase':phase,
                  'config':config,'audit':audited,'audit_seconds':audit_seconds,'migration':migration,
                  'new_or_resumed_words':new,'new_work_seconds':time.perf_counter()-started,
                  'peak_in_flight':peak,'peak_per_rank_in_flight':peak_per_rank,'stop_reason':stop,
                  'worker_memory_and_cache':states,'ranks':summarize(db,targets,reasons)}
            temporary=output/'summary.json.tmp';temporary.write_text(_json(data)+'\n');temporary.replace(output/'summary.json')
            return data
        run_started=True
        snapshot(phase='starting')
        if admitting:
            pool=_executor_factory(max_workers=workers,initializer=_initialize,
                                   initargs=(asdict(limits),expected,cache_bytes,cache_entries,rule_cache_bytes,rule_cache_entries,
                                             verifier_cache_bytes,verifier_cache_entries))
        while (active and admitting) or any(queues.values()) or discarded:
            assert pool is not None
            discarded={future for future in discarded if not future.done()}
            if admitting and drain_requested[0]:
                admitting=False;stop='signal_drain'
            if admitting and deadline is not None and datetime.now(timezone.utc)>=deadline:
                admitting=False;stop='admission_deadline'
            if admitting and time.perf_counter()-started>=wall_seconds:
                admitting=False;stop='wall_budget'
            if admitting and max_new_words is not None and new>=max_new_words:
                admitting=False;stop='word_budget'
            occupancy=sum(len(q) for q in queues.values())+len(discarded)
            per_rank_cap=max(1,math.ceil(2*workers/max(1,len(active))))
            if admitting:
                for _ in range(2*workers):
                    if deadline is not None and datetime.now(timezone.utc)>=deadline:
                        admitting=False;stop='admission_deadline';break
                    if occupancy>=2*workers or (max_new_words is not None and new+occupancy>=max_new_words):break
                    chosen=None
                    for _ in range(len(rotation)):
                        n=rotation[0];rotation.rotate(-1)
                        if (n in active and next_submit[n][0]<=targets[n]
                                and len(queues[n])<per_rank_cap):chosen=n;break
                    if chosen is None:break
                    n=chosen;cursor=next_submit[n][:];word=list(_word_at(n,*cursor))
                    queues[n].append((cursor,word,pool.submit(_compute,n,word)))
                    next_submit[n]=_advance(n,cursor);occupancy+=1;peak=max(peak,occupancy)
                    peak_per_rank[n]=max(peak_per_rank[n],len(queues[n]))
            ready=[n for n in rotation if queues[n] and queues[n][0][2].done()]
            if not ready:
                futures=[job[2] for q in queues.values() for job in q if not job[2].done()]+list(discarded)
                if not futures:
                    if not admitting or not active:break
                    raise RuntimeError('active rank has no pending literal input')
                wait(futures,timeout=.2,return_when=FIRST_COMPLETED)
                if time.perf_counter()-last_snapshot>=5:
                    snapshot(phase='computing' if admitting else 'draining');last_snapshot=time.perf_counter()
                continue
            _assert_source(expected)
            for n in ready:
                cursor,word,future=queues[n].popleft()
                if cursor!=cursors[n]:raise ValueError('worker is not the next literal input')
                encoded,producer,checker,telemetry=future.result()
                status,following=_commit_word(db,n,cursor,word,encoded,producer,checker,telemetry,config_sha,expected)
                cursors[n]=following;new+=1
                if status!='complete' or following[0]>targets[n]:
                    reasons[n]='target_length' if status=='complete' else status
                    active.discard(n)
                    for _,_,later in queues[n]:
                        later.cancel()
                        if not later.done():discarded.add(later)
                    queues[n].clear()
            if time.perf_counter()-last_snapshot>=5:
                snapshot(phase='computing' if admitting else 'draining');last_snapshot=time.perf_counter()
        if pool is not None:pool.shutdown(wait=True,cancel_futures=True);pool=None
        _assert_source(expected)
        complete=_goals_complete(db,targets)
        if complete:stop='targets_completed'
        elif stop=='running':stop='rank_limits'
        return snapshot(state='completed' if complete else 'paused',phase='idle')
    except BaseException as error:
        if run_started:
            try:
                db.rollback()
                interrupted=isinstance(error,(KeyboardInterrupt,SystemExit))
                stop='interrupted' if interrupted else 'failed'
                snapshot(state='paused' if interrupted else 'failed',phase='error')
            except Exception:
                pass  # The original failure remains primary if the database is unusable.
        raise
    finally:
        try:
            for q in queues.values():
                for _,_,future in q:future.cancel()
            if pool is not None:pool.shutdown(wait=True,cancel_futures=True)
        finally:
            try:
                db.close();fcntl.flock(lock.fileno(),fcntl.LOCK_UN);lock.close()
            finally:
                _restore_drain_handler(previous_handler)


def status(output):
    """Read a consistent, unaudited snapshot; recorded process status may be stale."""
    database=Path(output).resolve()/'words.sqlite3'
    if not database.is_file():raise ValueError('status requires an existing database')
    db=sqlite3.connect(database.as_uri()+'?mode=ro',uri=True)
    try:
        db.execute('BEGIN')
        metadata=dict(db.execute('SELECT key,value FROM metadata'))
        versions,active_version=_load_versions(db);expected=versions[active_version]
        config=_loads(metadata['last_config']);config_sha=_sha(metadata['last_config'])
        targets={int(n):L for n,L in config['targets'].items()}
        states=[_loads(row[0]) for row in db.execute('SELECT runtime FROM worker_state WHERE config_sha=?',(config_sha,))] if _has_table(db,'worker_state') else []
        return {'schema':expected['schema'],'identity':expected,'database':str(database),
                'active_version':active_version,'versions':versions,
                'legacy_identity':_loads(metadata['identity']),'identity_scope':'active_version',
                'coverage_scope':'all_registered_versions',
                'snapshot_utc':datetime.now(timezone.utc).isoformat(),'config':config,
                'recorded_run_state':_loads(metadata['run_state']) if 'run_state' in metadata else None,
                'prefix_audited_by_this_read':False,'process_liveness_verified':False,
                'worker_memory_and_cache':states,'ranks':summarize(db,targets,{})}
    finally:db.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output');parser.add_argument('--targets',type=json.loads)
    parser.add_argument('--workers',type=int,default=8);parser.add_argument('--batch-size',type=int,default=64)
    parser.add_argument('--wall-seconds',type=float,default=86400);parser.add_argument('--limits',type=json.loads)
    parser.add_argument('--admission-deadline-utc')
    parser.add_argument('--max-new-words',type=int);parser.add_argument('--replay-audit',action='store_true')
    parser.add_argument('--audit-only',action='store_true');parser.add_argument('--status',action='store_true')
    parser.add_argument('--migrate',action='store_true');parser.add_argument('--migration-backup')
    parser.add_argument('--execution',choices=('native','memoized','collected'),default='native')
    parser.add_argument('--verifier',choices=('python','cpp'),default='python')
    parser.add_argument('--verifier-allocator')
    parser.add_argument('--cache-bytes',type=int,default=256*1024*1024)
    parser.add_argument('--cache-entries',type=int,default=262144)
    parser.add_argument('--verifier-cache-bytes',type=int,default=256*1024*1024)
    parser.add_argument('--verifier-cache-entries',type=int,default=262144)
    parser.add_argument('--rule-cache-bytes',type=int,default=4*1024*1024)
    parser.add_argument('--rule-cache-entries',type=int,default=256)
    args=vars(parser.parse_args())
    if args.pop('status'):
        print(_json(status(args['output'])),flush=True);return
    if args['targets'] is not None:args['targets']={int(n):L for n,L in args['targets'].items()}
    if args['limits'] is not None:args['limits']=Limits(**args['limits'])
    print(_json(exhaust(**args)),flush=True)


if __name__=='__main__':main()
