"""Native census migration conserves old rows and independently checks new rows."""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import json
import os
import signal
import sqlite3
import zlib

import pytest

pytest.importorskip('fcntl',reason='experiment requires Unix flock')
from experiments import zipper_exhaustive_collected as runner
from zipper_algebra.zipping_algorithm.interpreter import Limits


@pytest.mark.timeout(300)
def test_native_default_migration_preserves_completed_rows_and_audits(tmp_path):
    runner.exhaust(tmp_path,targets={2:1},workers=1,wall_seconds=30,
                   execution='collected',_executor_factory=ThreadPoolExecutor)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        original=db.execute('SELECT * FROM words ORDER BY n,length,ordinal').fetchall()
        old_identity=db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0]
    result=runner.exhaust(tmp_path,targets={2:2,3:1},workers=1,wall_seconds=30,migrate=True,
                          cache_bytes=1024*1024,cache_entries=128,
                          verifier_cache_bytes=2*1024*1024,verifier_cache_entries=256,
                          rule_cache_bytes=65536,rule_cache_entries=32)
    assert result['status']=='completed' and result['identity']['execution']=='native'
    assert result['identity']['production_check']==runner.NATIVE_CHECK
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        assert db.execute('SELECT * FROM words WHERE n=2 AND length<=1 ORDER BY n,length,ordinal').fetchall()==original
        assert db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0]==old_identity
        assert db.execute('SELECT count(*) FROM versions').fetchone()[0]==2
        assert db.execute('SELECT count(*) FROM words').fetchone()[0]==12
        saved=db.execute('SELECT certificate,runtime FROM words JOIN word_runtime USING(n,length,ordinal) WHERE n=3').fetchall()
        for blob,encoded in saved:
            certificate=json.loads(zlib.decompress(blob));info=json.loads(encoded)
            assert certificate['schema']=='zipper-native-run-v1'
            assert certificate['native_build']==result['identity']['native_build']
            assert info['checker_mode']=='independent-python-cached-replay'
            assert info['after_reduce']['max_bytes']==1024*1024
            assert info['verifier_after']['max_bytes']==2*1024*1024
            assert info['native_pid']!=info['worker_pid']
        assert any(json.loads(encoded)['verifier_delta']['hits'] for _,encoded in saved)
    backup=sqlite3.connect(result['migration']['backup'])
    try:assert backup.execute('SELECT * FROM words ORDER BY n,length,ordinal').fetchall()==original
    finally:backup.close()
    assert runner.exhaust(tmp_path,audit_only=True)['audit']['words']==12
    # Native build and independent verifier bindings survive a complete database audit.
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        n,length,ordinal,config_sha,encoded,digest=db.execute(
            'SELECT * FROM word_runtime WHERE n=3 LIMIT 1').fetchone()
        info=json.loads(encoded);info['verifier_after']['max_bytes']+=1
        encoded=runner._json(info)
        db.execute('UPDATE word_runtime SET runtime=?,runtime_sha=? WHERE n=? AND length=? AND ordinal=?',
                   (encoded,runner._sha(encoded),n,length,ordinal))
    with pytest.raises(ValueError,match='verifier.*(delta|budgets)'):
        runner.exhaust(tmp_path,audit_only=True)


def test_native_certificate_build_and_reference_bindings():
    expected={'execution':'native','certificate_schema':'zipper-native-run-v1',
              'algorithm_source_sha256':{'a':'b'},'catalog_sha256':'cat','coefficient_ring':'Z[q,q^-1]',
              'batch_size':64,'native_build':{'binary_sha256':'binary'},
              'python_reference_source_sha256':{'python':'source'}}
    certificate={'schema':'zipper-native-run-v1','n':3,'input_word':[1],
                 'source_sha256':{'a':'b'},'catalog_sha256':'cat','limits':{},'trace_enabled':False,
                 'coefficient_ring':'Z[q,q^-1]','execution':{'batch_size':64},'status':'complete',
                 'native_build':{'binary_sha256':'binary'},'python_reference_source_sha256':{'python':'source'}}
    runner._certificate_identity(certificate,3,[1],{},expected)
    for key in ('native_build','python_reference_source_sha256'):
        corrupt=deepcopy(certificate);corrupt[key]={}
        with pytest.raises(ValueError,match='build or reference'):
            runner._certificate_identity(corrupt,3,[1],{},expected)
    for key,value in [('n',3.0),('input_word',[1.0]),('input_word',[True]),('execution',{'batch_size':64.0})]:
        corrupt=deepcopy(certificate);corrupt[key]=value
        with pytest.raises(ValueError,match='identity differs'):
            runner._certificate_identity(corrupt,3,[1],{},expected)


def test_graceful_signal_drain_commits_admitted_prefix_and_restores_handler(tmp_path,monkeypatch):
    previous=signal.getsignal(signal.SIGUSR1)
    actual=runner._compute
    def request(n,word):
        os.kill(os.getpid(),signal.SIGUSR1)
        return actual(n,word)
    monkeypatch.setattr(runner,'_compute',request)
    result=runner.exhaust(tmp_path,targets={2:10},workers=1,wall_seconds=30,
                          execution='collected',_executor_factory=ThreadPoolExecutor)
    assert result['status']=='paused' and result['stop_reason']=='signal_drain'
    assert 1<=result['new_or_resumed_words']<=2
    assert signal.getsignal(signal.SIGUSR1)==previous
    assert runner.exhaust(tmp_path,audit_only=True)['audit']['words']==result['new_or_resumed_words']


def test_drain_handler_restored_after_worker_failure(tmp_path,monkeypatch):
    previous=signal.getsignal(signal.SIGUSR1)
    def fail(n,word):raise RuntimeError('worker stopped')
    monkeypatch.setattr(runner,'_compute',fail)
    with pytest.raises(RuntimeError,match='worker stopped'):
        runner.exhaust(tmp_path,targets={2:1},workers=1,wall_seconds=30,
                       execution='collected',_executor_factory=ThreadPoolExecutor)
    assert signal.getsignal(signal.SIGUSR1)==previous


def test_native_verifier_cache_budgets_are_validated_before_database_creation(tmp_path):
    with pytest.raises(ValueError,match='cache budgets'):
        runner.exhaust(tmp_path,verifier_cache_entries=-1)
    assert not (tmp_path/'words.sqlite3').exists()


@pytest.mark.timeout(300)
def test_native_checkpoint_retry_archives_old_attempt(tmp_path):
    first=runner.exhaust(tmp_path,targets={3:0},workers=1,wall_seconds=30,
                         execution='collected',limits=Limits(max_steps=1),
                         _executor_factory=ThreadPoolExecutor)
    assert first['status']=='paused'
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        old=db.execute('SELECT * FROM words').fetchone()
    result=runner.exhaust(tmp_path,targets={3:0},workers=1,wall_seconds=30,migrate=True,
                          limits=Limits(max_steps=100))
    assert result['status']=='completed'
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        archived=db.execute('SELECT * FROM word_attempts').fetchone()
        assert archived[:17]==old
    audit=runner.exhaust(tmp_path,audit_only=True)['audit']
    assert audit['archived_attempts']==1 and audit['words']==1


def test_exact_json_wide_integers_preserves_standard_codec_and_global_limit():
    import sys
    from zipper_algebra.zipping_algorithm import enumeration
    from zipper_algebra.zipping_algorithm.native import _integer_decimal
    limit=sys.get_int_max_str_digits()
    ordinary={'integer':123,'coefficient':[[0,-127],[1,2]],'unicode':'α',
              'nested':[True,False,None,1.5,{'a':[1,2]}]}
    assert runner._json(ordinary)==enumeration._json(ordinary)
    huge=2**32770+17
    value={'coefficient':[[0,huge],[1,-huge]],'ordinary':ordinary}
    encoded=runner._json(value)
    assert _integer_decimal(huge) in encoded
    assert runner._loads(encoded)==value
    assert sys.get_int_max_str_digits()==limit
    # Standard failures remain failures; the fallback must not admit NaN.
    with pytest.raises(ValueError):runner._json({'huge':huge,'nan':float('nan')})


def test_file_hashes_stream_without_read_bytes(tmp_path,monkeypatch):
    import hashlib
    from pathlib import Path
    payload=bytes(range(256))*8193
    path=tmp_path/'backup.sqlite3';path.write_bytes(payload)
    expected=hashlib.sha256(payload).hexdigest()
    def forbidden(*args,**kwargs):raise AssertionError('whole-file read forbidden')
    monkeypatch.setattr(Path,'read_bytes',forbidden)
    assert runner._hash(path)==expected
    assert runner._file_sha(path)==expected


def test_expired_absolute_deadline_admits_no_work_and_restores_handler(tmp_path):
    previous=signal.getsignal(signal.SIGUSR1)
    def forbidden(**kwargs):raise AssertionError('expired deadline started workers')
    result=runner.exhaust(tmp_path,targets={2:1},execution='collected',
                          admission_deadline_utc='2000-01-01T00:00:00Z',_executor_factory=forbidden)
    assert result['stop_reason']=='admission_deadline'
    assert result['status']=='paused' and result['new_or_resumed_words']==0
    assert result['config']['admission_deadline_utc']=='2000-01-01T00:00:00+00:00'
    assert signal.getsignal(signal.SIGUSR1)==previous


def test_absolute_deadline_includes_startup_audit_time(tmp_path,monkeypatch):
    from datetime import datetime, timezone
    before=datetime(2030,1,1,tzinfo=timezone.utc)
    after=datetime(2030,1,3,tzinfo=timezone.utc)
    audited=[False]
    class Clock(datetime):
        @classmethod
        def now(cls,tz=None):return after if audited[0] else before
    original=runner.audit
    def audit(*args,**kwargs):
        result=original(*args,**kwargs);audited[0]=True;return result
    def forbidden(**kwargs):raise AssertionError('audit time extended deadline')
    monkeypatch.setattr(runner,'datetime',Clock)
    monkeypatch.setattr(runner,'audit',audit)
    result=runner.exhaust(tmp_path,targets={2:1},execution='collected',wall_seconds=3600,
                          admission_deadline_utc='2030-01-02T00:00:00Z',_executor_factory=forbidden)
    assert audited[0] and result['new_or_resumed_words']==0
    assert result['stop_reason']=='admission_deadline'


@pytest.mark.parametrize('deadline',['2030-01-01','2030-01-01T12:00:00',
                                    '2030-01-01T12:00:00+01:00','invalid',123])
def test_absolute_deadline_requires_aware_utc_before_creating_store(tmp_path,deadline):
    with pytest.raises(ValueError,match='timezone-aware UTC'):
        runner.exhaust(tmp_path,admission_deadline_utc=deadline)
    assert not (tmp_path/'words.sqlite3').exists()


def test_absolute_deadline_stops_new_admissions_and_drains_checked_prefix(tmp_path,monkeypatch):
    from datetime import datetime, timezone
    submitted=[False]
    class Clock(datetime):
        @classmethod
        def now(cls,tz=None):
            return datetime(2030,1,3 if submitted[0] else 1,tzinfo=timezone.utc)
    class Pool(ThreadPoolExecutor):
        def submit(self,*args,**kwargs):
            result=super().submit(*args,**kwargs);submitted[0]=True;return result
    monkeypatch.setattr(runner,'datetime',Clock)
    result=runner.exhaust(tmp_path,targets={2:10},execution='collected',workers=1,
                          admission_deadline_utc='2030-01-02T00:00:00Z',_executor_factory=Pool)
    assert result['stop_reason']=='admission_deadline'
    assert result['new_or_resumed_words']==1
    assert result['ranks'][0]['complete_words']==1
    assert runner.exhaust(tmp_path,audit_only=True)['audit']['words']==1
