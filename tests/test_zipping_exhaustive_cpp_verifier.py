"""C++ census verification preserves source epochs, certificates, and allocator provenance."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sqlite3

import pytest

pytest.importorskip('fcntl',reason='experiment requires Unix flock')

from experiments import zipper_exhaustive_collected as runner
from zipper_algebra.zipping_algorithm.interpreter import Limits

def stop_sessions(module):
    if module._verifier is not None and hasattr(module._verifier,'close'):module._verifier.close()
    if module._session is not None:module._session.close()


def test_migration_conserves_old_native_rows_and_runtime(tmp_path):
    try:
        runner.exhaust(tmp_path,verifier='python',targets={2:1},workers=1,wall_seconds=30,_executor_factory=ThreadPoolExecutor)
        with sqlite3.connect(tmp_path/'words.sqlite3') as db:
            saved=db.execute('SELECT * FROM words ORDER BY n,length,ordinal').fetchall()
            runtime=db.execute('SELECT * FROM word_runtime ORDER BY n,length,ordinal').fetchall()
            oldmeta=dict(db.execute('SELECT * FROM metadata'))
        result=runner.exhaust(tmp_path,targets={2:2,3:1},workers=1,wall_seconds=30,verifier='cpp',migrate=True,
                              _executor_factory=ThreadPoolExecutor)
        assert result['status']=='completed'
        assert result['identity']['production_check']==runner.CPP_CHECK
        assert result['identity']['verifier_backend']=='cpp'
        assert result['identity']['verifier_build']['libraries']['gmp']['sha256']
        with sqlite3.connect(tmp_path/'words.sqlite3') as db:
            assert db.execute('SELECT * FROM words WHERE n=2 AND length<=1 ORDER BY n,length,ordinal').fetchall()==saved
            assert db.execute('SELECT * FROM word_runtime WHERE n=2 AND length<=1 ORDER BY n,length,ordinal').fetchall()==runtime
            assert db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0]==oldmeta['identity']
            assert db.execute('SELECT count(*) FROM versions').fetchone()[0]==2
            audit=runner.audit(db)
            assert audit['words']==12 and audit['replayed']==0
            row=db.execute('SELECT * FROM word_runtime WHERE n=3 LIMIT 1').fetchone()
            info=runner._loads(row[4])
            assert info['checker_mode']=='independent-cpp-cached-replay'
            assert info['verifier_pid']!=info['native_pid']!=info['worker_pid']
            assert info['verifier_process_memory']['native_process_current_rss_bytes']>0
            assert 'rule_cache' not in info['verifier_after']
            info['verifier_build']['binary_sha256']='tampered'
            encoded=runner._json(info)
            db.execute('UPDATE word_runtime SET runtime=?,runtime_sha=? WHERE n=? AND length=? AND ordinal=?',
                       (encoded,runner._sha(encoded),*row[:3]))
            with pytest.raises(ValueError,match='verifier build'):runner.audit(db)
        with sqlite3.connect(result['migration']['backup']) as db:
            assert db.execute('SELECT * FROM words ORDER BY n,length,ordinal').fetchall()==saved
    finally:
        stop_sessions(runner)


def test_migration_keeps_expired_absolute_admission_deadline(tmp_path):
    runner.exhaust(tmp_path,verifier='python',targets={2:0},workers=1,wall_seconds=30,_executor_factory=ThreadPoolExecutor)
    try:
        result=runner.exhaust(tmp_path,targets={2:2},workers=1,verifier='cpp',migrate=True,
                              admission_deadline_utc='2000-01-01T00:00:00Z',
                              _executor_factory=lambda **kwargs:pytest.fail('expired deadline admitted workers'))
        assert result['stop_reason']=='admission_deadline'
        assert result['new_or_resumed_words']==0
        assert result['config']['admission_deadline_utc']=='2000-01-01T00:00:00+00:00'
        assert result['migration']['audited_words']==1
    finally:stop_sessions(runner)


def test_cpp_checkpoint_retry_archives_exact_old_attempt(tmp_path):
    try:
        runner.exhaust(tmp_path,verifier='python',targets={3:0},workers=1,limits=Limits(max_steps=1),_executor_factory=ThreadPoolExecutor)
        with sqlite3.connect(tmp_path/'words.sqlite3') as db: saved=db.execute('SELECT * FROM words').fetchone()
        result=runner.exhaust(tmp_path,targets={3:0},workers=1,verifier='cpp',migrate=True,
                              limits=Limits(max_steps=100),_executor_factory=ThreadPoolExecutor)
        assert result['status']=='completed'
        with sqlite3.connect(tmp_path/'words.sqlite3') as db:
            assert db.execute('SELECT * FROM word_attempts').fetchone()[:17]==saved
            assert runner.audit(db)['archived_attempts']==1
    finally:stop_sessions(runner)


def test_cpp_source_and_config_binding(tmp_path,monkeypatch):
    identity=runner.identity(64,'native','cpp')
    runner._assert_source(identity)
    monkeypatch.setattr(runner,'_verifier_source',lambda backend: {})
    with pytest.raises(ValueError,match='verifier source changed'):runner._assert_source(identity)
    with pytest.raises(ValueError,match=r'C\+\+ verification'):runner.exhaust(tmp_path,execution='collected',verifier='cpp')
    assert not (tmp_path/'words.sqlite3').exists()


def test_allocator_migration_binds_runtime_and_preserves_default_cpp_epoch(tmp_path):
    library=Path('/lib64/libtcmalloc_minimal.so.4')
    if not library.is_file():pytest.skip('tcmalloc minimal is unavailable')
    try:
        runner.exhaust(tmp_path,targets={2:1},workers=1,verifier='cpp',_executor_factory=ThreadPoolExecutor)
        with sqlite3.connect(tmp_path/'words.sqlite3') as db:
            saved=db.execute('SELECT * FROM words ORDER BY n,length,ordinal').fetchall()
        result=runner.exhaust(tmp_path,targets={2:2},workers=1,verifier='cpp',verifier_allocator=library,
                              migrate=True,_executor_factory=ThreadPoolExecutor)
        assert result['status']=='completed'
        allocator=result['identity']['verifier_allocator']
        assert allocator==runner._verifier_allocator(library)
        assert result['config']['verifier_allocator']==allocator
        runner._assert_source(result['identity'])
        with sqlite3.connect(tmp_path/'words.sqlite3') as db:
            assert db.execute('SELECT * FROM words WHERE length<=1 ORDER BY n,length,ordinal').fetchall()==saved
            assert db.execute('SELECT count(*) FROM versions').fetchone()[0]==2
            assert runner.audit(db)['words']==7
            row=db.execute('SELECT * FROM word_runtime WHERE length=2 LIMIT 1').fetchone()
            config=runner._loads(db.execute('SELECT config FROM configurations WHERE hash=?',(row[3],)).fetchone()[0])
            info=runner._loads(row[4])
            assert info['verifier_allocator']==dict(allocator,mapped_verified=True)
            for changed in [dict(allocator,mapped_verified=False),dict(allocator,mapped_verified=1),
                            dict(allocator,sha256='changed',mapped_verified=True),None]:
                info['verifier_allocator']=changed;encoded=runner._json(info)
                with pytest.raises(ValueError,match='allocator identity'):
                    runner._validate_runtime(encoded,runner._sha(encoded),config)
    finally:stop_sessions(runner)


def test_allocator_source_guard_detects_changed_library(tmp_path):
    library=Path('/lib64/libtcmalloc_minimal.so.4')
    if not library.is_file():pytest.skip('tcmalloc minimal is unavailable')
    copy=tmp_path/'allocator.so';copy.write_bytes(library.read_bytes())
    identity=runner.identity(64,'native','cpp',copy)
    runner._assert_source(identity)
    with copy.open('ab') as output:output.write(b'changed')
    with pytest.raises(ValueError,match='allocator identity changed'):runner._assert_source(identity)


@pytest.mark.parametrize('execution,verifier',[('native','python'),('collected','python'),('memoized','python')])
def test_allocator_requires_cpp_before_database_creation(tmp_path,execution,verifier):
    with pytest.raises(ValueError,match='allocator requires'):
        runner.exhaust(tmp_path,execution=execution,verifier=verifier,verifier_allocator='/lib64/libtcmalloc_minimal.so.4')
    assert not (tmp_path/'words.sqlite3').exists()



def test_failed_verifier_cleanup_still_closes_producer(monkeypatch):
    from zipper_algebra.zipping_algorithm import native_verifier
    closed=[]
    class Producer:
        metadata={}
        def close(self):closed.append('producer')
    class Verifier:
        metadata={}
        def __init__(self,**kwargs):pass
        def close(self):
            closed.append('verifier')
            raise RuntimeError('verifier cleanup interrupted')
    monkeypatch.setattr(runner,'_session',None)
    monkeypatch.setattr(runner,'_verifier',None)
    monkeypatch.setattr(runner,'_new_session',lambda *args:Producer())
    monkeypatch.setattr(runner,'_verifier_build',lambda metadata:{'build':'different'})
    monkeypatch.setattr(native_verifier,'CppVerifier',Verifier)
    expected={'execution':'native','verifier_backend':'cpp','verifier_build':{'build':'expected'}}
    with pytest.raises(RuntimeError,match='cleanup interrupted'):
        runner._initialize({},expected,1024,1,1024,1)
    assert closed==['verifier','producer']
