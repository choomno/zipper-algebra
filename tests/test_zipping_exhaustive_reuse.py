from concurrent.futures import Future
from concurrent.futures.process import BrokenProcessPool
from dataclasses import asdict
import json
import sqlite3

import pytest

from experiments import zipper_exhaustive_reuse as runner
from zipper_algebra.zipping_algorithm.census_reuse import ReuseCensus
from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import execution_identity


@pytest.fixture
def fake_pool(monkeypatch):
    class Pool:
        def __init__(self, *, initargs, **kwargs):
            path,settings=initargs
            session=CollectedSession()
            limits=Limits(**settings['limits'])
            seed=session.reduce([],n=1,limits=limits)
            self.store=ReuseCensus(path,produce=lambda n,w:session.reduce(w,n=n,limits=limits),
                                  verify=session.check,producer_identity=execution_identity(seed),
                                  verifier_identity={'name':'test-reference'},limits=asdict(limits))
        def submit(self, fn, *args):
            runner._store=self.store
            result=Future()
            try:
                result.set_result(fn(*args))
            except BaseException as e:
                result.set_exception(e)
            return result
        def shutdown(self, **kwargs):
            self.store.close()
            runner._store=None
    monkeypatch.setattr(runner,'ProcessPoolExecutor',Pool)


def test_round_robin_resume_and_audit(tmp_path,fake_pool):
    for i in range(4):
        result=runner.run(tmp_path,targets={2:3,3:2},max_new_words=1)
        assert result['jobs']=={'complete':i+1}
    db=sqlite3.connect(tmp_path/'census.sqlite3')
    assert [r[0] for r in db.execute('SELECT n FROM reuse_jobs ORDER BY rowid')]==[2,3,2,3]
    assert runner.audit_store(tmp_path)=={'jobs':4,'complete':4}
    db.execute("UPDATE reuse_jobs SET result='{}' WHERE rowid=1")
    db.commit();db.close()
    with pytest.raises(ValueError,match='result binding'):
        runner.audit_store(tmp_path)


def test_explicit_dedup_and_missing_shard(tmp_path,fake_pool):
    result=runner.run(tmp_path,words=[(4,[1,3]),(4,[1,3]),(2,[1])],workers=2)
    assert result['admitted']==2 and result['jobs']=={'complete':2}
    assert runner.run(tmp_path,words=[(4,[1,3])])['admitted']==0
    (tmp_path/'shards/worker-1.sqlite3').unlink()
    with pytest.raises(sqlite3.OperationalError):
        runner.audit_store(tmp_path)


def test_incomplete_attempts_and_retry(tmp_path,fake_pool):
    result=runner.run(tmp_path,words=[(3,[1,2,1])],limits=Limits(max_steps=1))
    assert result['jobs']=={'incomplete':1}
    assert runner.run(tmp_path,words=[(3,[1,2,1])],limits=Limits(max_steps=1))['admitted']==0
    assert runner.run(tmp_path,words=[],limits=Limits(max_steps=1),retry_incomplete=True)['admitted']==1
    assert runner.audit_store(tmp_path)=={'jobs':1,'complete':0}


def test_interrupted_retry_and_preserved_cursor(tmp_path,fake_pool):
    runner.run(tmp_path,targets={2:1,3:2},max_new_words=1)
    db=sqlite3.connect(tmp_path/'census.sqlite3')
    db.execute("UPDATE reuse_jobs SET status='running',result=NULL")
    db.commit();db.close()
    result=runner.run(tmp_path,targets={2:1,3:2},max_new_words=1)
    assert result['admitted']==1 and result['jobs']=={'complete':1}
    # Retrying the dispatched root does not advance the cursor twice.
    result=runner.run(tmp_path,targets={2:1,3:2},max_new_words=1)
    assert result['jobs']=={'complete':2}


def test_domain_counts_and_validation(tmp_path,fake_pool):
    assert runner._small_layer(1,0)==[()]
    assert runner._small_layer(1,1)==[]
    assert runner._small_layer(2,1)==[(1,)]
    assert runner._small_layer(2,2)==[(1,1),(-1,-1)]
    assert sum(len(runner._small_layer(2,l)) for l in range(1,11))==19
    assert runner._next(3,[2,0],2)[0] == [1,-2]
    assert runner._next(3,[2,1000],2) is None
    assert runner.audit_store(tmp_path)=={'jobs':0,'complete':0}
    result=runner.run(tmp_path,targets={1:0,2:1},workers=2)
    assert result['jobs']=={'complete':2}
    for kw in ({'workers':0},{'wall_seconds':0},{'mode':'bad'},{'max_new_words':-1},{'targets':{0:1}}):
        with pytest.raises(ValueError):
            runner.run(tmp_path,**kw)
    db=sqlite3.connect(tmp_path/'legacy.sqlite3');db.execute('CREATE TABLE words(x)')
    with pytest.raises(ValueError,match='own reuse'):
        runner._ledger(db)
    db.close()


def test_failure_recorded_and_root_validation(tmp_path,fake_pool,monkeypatch):
    original=runner._evaluate
    monkeypatch.setattr(runner,'_evaluate',lambda n,w:{'policy':'bad','root':'bad'})
    result=runner.run(tmp_path,words=[(2,[1])])
    assert result['jobs']=={'failed':1}
    monkeypatch.setattr(runner,'_evaluate',original)
    assert runner.run(tmp_path,words=[],retry_incomplete=True)['jobs']=={'complete':1}


def test_cli_audit(tmp_path,monkeypatch,capsys):
    monkeypatch.setattr('sys.argv',['zipper_exhaustive_reuse',str(tmp_path),'--audit-only'])
    runner.main()
    assert json.loads(capsys.readouterr().out)=={'jobs':0,'complete':0}


@pytest.mark.parametrize('when', ['submit', 'result'])
def test_broken_worker_stops_admission_and_retries_exact_job(tmp_path,fake_pool,monkeypatch,when):
    pool_type=runner.ProcessPoolExecutor
    original=pool_type.submit
    def broken(self,fn,*args):
        if fn is runner._evaluate:
            if when=='submit':
                raise BrokenProcessPool('terminated worker')
            result=Future()
            result.set_exception(BrokenProcessPool('terminated worker'))
            return result
        return original(self,fn,*args)
    monkeypatch.setattr(pool_type,'submit',broken)
    result=runner.run(tmp_path,targets={2:4,3:3},workers=1,max_new_words=8)
    assert result['admitted']==1
    assert result['jobs']=={'interrupted':1}
    assert result['stop_reason']=='worker_failure'
    monkeypatch.setattr(pool_type,'submit',original)
    result=runner.run(tmp_path,targets={2:4,3:3},max_new_words=1)
    assert result['jobs']=={'complete':1}
    assert runner.audit_store(tmp_path)=={'jobs':1,'complete':1}


def test_counts_scanned_only_at_summary(tmp_path,fake_pool,monkeypatch):
    original=ReuseCensus.counts
    calls=[]
    def count(self):
        calls.append(self.path)
        return original(self)
    monkeypatch.setattr(ReuseCensus,'counts',count)
    result=runner.run(tmp_path,targets={2:4},max_new_words=5)
    assert result['admitted']==5 and len(calls)==1


def test_pending_keyset_scan_is_bounded_and_does_not_retry_new_failures():
    db=sqlite3.connect(':memory:')
    runner._ledger(db)
    db.executemany('INSERT INTO reuse_jobs VALUES (?,?,?,?,?,?,?)',
                   [('policy',str(i),2,'[1]','incomplete',None,None) for i in range(600)])
    pending=runner._pending_jobs(db,'policy',True)
    assert next(pending)[0]=='0'
    db.execute('INSERT INTO reuse_jobs VALUES (?,?,?,?,?,?,?)',
               ('policy','later',2,'[1]','failed',None,None))
    rest=list(pending)
    assert len(rest)==599 and rest[-1][0]=='599'
    db.close()


def test_readonly_coordinator_probe_preserves_legacy_wal(tmp_path,fake_pool):
    path=tmp_path/'census.sqlite3'
    original=sqlite3.connect(path)
    original.execute('PRAGMA journal_mode=WAL')
    original.execute('PRAGMA wal_autocheckpoint=0')
    original.execute('CREATE TABLE words(value TEXT)')
    original.execute("INSERT INTO words VALUES ('preserved')")
    original.commit()
    before=(path.read_bytes(),path.with_name(path.name+'-wal').read_bytes())
    with pytest.raises(ValueError,match='own reuse'):
        runner.run(tmp_path,max_new_words=1)
    assert (path.read_bytes(),path.with_name(path.name+'-wal').read_bytes())==before
    original.close()


def test_worker_seed_must_be_independently_accepted(monkeypatch):
    closed=[]
    class Producer:
        metadata={}
        def __init__(self,**kwargs):
            pass
        def reduce_reused(self,*args,**kwargs):
            return {}
        def close(self):
            closed.append('producer')
    class Verifier:
        def __init__(self,**kwargs):
            pass
        def check_reused(self,result):
            return False
        def close(self):
            closed.append('verifier')
    monkeypatch.setattr(runner,'NativeSession',Producer)
    monkeypatch.setattr(runner,'CppVerifier',Verifier)
    with pytest.raises(ValueError,match='initialization'):
        runner._initialize('unused',{'cache_bytes':0,'cache_entries':0,'allocator':None,
                                     'mode':'reuse','limits':asdict(Limits())})
    assert closed==['verifier','producer']


@pytest.mark.parametrize('change', ["UPDATE reuse_cursors SET ordinal='999'",
                                  "UPDATE reuse_cursors SET length=1000000000",
                                  'DELETE FROM reuse_cursors'])
def test_corrupted_cursor_cannot_skip_inputs(tmp_path,fake_pool,change):
    runner.run(tmp_path,targets={2:3},max_new_words=1)
    with sqlite3.connect(tmp_path/'census.sqlite3') as db:
        db.execute(change)
    with pytest.raises(ValueError,match='cursor'):
        runner.audit_store(tmp_path)


def test_runner_import_binding_rejects_drift_before_store_open(tmp_path,monkeypatch):
    monkeypatch.setattr(runner,'RUNNER_SHA256','0'*64)
    with pytest.raises(ValueError,match='since import'):
        runner.run(tmp_path/'unused',max_new_words=0)
    assert not (tmp_path/'unused').exists()


def test_worker_runner_identity_checked_before_admission(tmp_path,fake_pool,monkeypatch):
    ready=runner._ready
    monkeypatch.setattr(runner,'_ready',lambda:dict(ready(),runner_sha256='0'*64))
    with pytest.raises(ValueError,match='worker runner source'):
        runner.run(tmp_path,targets={2:3},max_new_words=1)
    with sqlite3.connect(tmp_path/'census.sqlite3') as db:
        assert db.execute('SELECT count(*) FROM reuse_jobs').fetchone()[0]==0


@pytest.mark.parametrize('bounds',[{'max_steps':1},{'max_terms':1},
    dict(max_steps=1,max_terms=1,max_word_size=1,max_coefficient_support=1,max_integer_bits=1,expansion_budget=1)])
def test_native_low_limits_bootstrap_and_record_incomplete(tmp_path,bounds):
    limits=Limits(**bounds)
    runner._close()
    try:
        runner._initialize(str(tmp_path/'worker.sqlite3'),
                           {'cache_bytes':0,'cache_entries':0,'allocator':None,
                            'mode':'reuse','limits':asdict(limits)})
        assert runner._store is not None
        assert runner._store.policy['limits']==asdict(limits)
        sessions=runner._sessions
        result=runner._evaluate(2,[1])
        assert result['status']=='incomplete'
        assert runner._snapshot_counts()=={'leaves':0,'roots':0,'attempts':1}
        assert result['diagnostic']['leaf_status']=='incomplete'
        assert result['worker_telemetry']['work']['producer']['calls']==1
    finally:
        runner._close()
    assert runner._store is None and runner._sessions==()
    assert all(session.pid is None for session in sessions)
    runner._close()
