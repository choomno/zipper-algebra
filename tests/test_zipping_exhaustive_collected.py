"""Independent literal prefixes preserve all words and replay provenance."""
from concurrent.futures import ThreadPoolExecutor
import json
import sqlite3
import threading
import zlib

import pytest
pytest.importorskip('fcntl',reason='experiment requires Unix flock')
from experiments import zipper_exhaustive_collected as runner
from zipper_algebra.zipping_algorithm.interpreter import Limits


def run(path,**kwargs):
    return runner.exhaust(path,workers=2,wall_seconds=30,execution='collected',_executor_factory=ThreadPoolExecutor,**kwargs)


def rows(path):
    with sqlite3.connect(path/'words.sqlite3') as db:
        return db.execute('SELECT n,length,ordinal,word,status,certificate FROM words ORDER BY n,length,ordinal').fetchall()


def test_literal_layers_compressed_certificates_resume_and_explicit_replay(tmp_path):
    first=run(tmp_path,targets={2:2,3:1},max_new_words=4)
    assert first['new_or_resumed_words']==4
    result=run(tmp_path,targets={2:2,3:1})
    assert len(rows(tmp_path))==12
    assert result['audit']['words']==4 and result['audit']['replayed']==0
    assert result['peak_in_flight']<=4
    assert [(r['n'],r['complete_through_length']) for r in result['ranks']]==[(2,2),(3,1)]
    assert [json.loads(r[3]) for r in rows(tmp_path) if r[0]==2 and r[1]==2]==[[-1,-1],[-1,1],[1,-1],[1,1]]
    for row in rows(tmp_path):
        assert json.loads(zlib.decompress(row[5]))['input_word']==json.loads(row[3])
    audited=run(tmp_path,targets={2:2,3:1},audit_only=True,replay_audit=True)
    assert audited['audit']['words']==audited['audit']['hash_and_binding_checked']==audited['audit']['replayed']==12


def test_incomplete_rank_pauses_without_blocking_other_rank(tmp_path):
    run(tmp_path,targets={2:1,3:1},limits=Limits(max_steps=5))
    assert [(r[0],r[1],r[4]) for r in rows(tmp_path) if r[0]==3]==[(3,0,'incomplete')]
    assert any(r[0]==2 and json.loads(r[3])==[-1] and r[4]=='complete' for r in rows(tmp_path))
    repeated=run(tmp_path,targets={2:1,3:1},limits=Limits(max_steps=5))
    assert repeated['new_or_resumed_words']==0
    assert next(r for r in repeated['ranks'] if r['n']==3)['stop_reason']=='word_limit_needs_larger_limits'
    resumed=run(tmp_path,targets={2:1,3:1},limits=Limits(max_steps=100))
    assert all(r['complete_through_length']==1 for r in resumed['ranks'])
    assert len(rows(tmp_path))==8


def test_slow_rank_does_not_hold_other_rank_commits(tmp_path,monkeypatch):
    actual=runner._compute;seen=threading.Event()
    def delayed(n,word):
        if n==2 and not word:
            assert seen.wait(5)
        result=actual(n,word)
        if n==3:seen.set()
        return result
    monkeypatch.setattr(runner,'_compute',delayed)
    assert run(tmp_path,targets={2:1,3:1})['new_or_resumed_words']==8


def test_worker_failure_keeps_next_input_and_other_prefixes_valid(tmp_path,monkeypatch):
    actual=runner._compute
    def fail(n,word):
        if n==2 and word==[-1]:raise RuntimeError('worker failure')
        return actual(n,word)
    monkeypatch.setattr(runner,'_compute',fail)
    with pytest.raises(RuntimeError,match='worker failure'):run(tmp_path,targets={2:1,3:1})
    assert [(r[1],r[2]) for r in rows(tmp_path) if r[0]==2]==[(0,0)]
    monkeypatch.setattr(runner,'_compute',actual)
    assert run(tmp_path,targets={2:1,3:1})['status']=='completed'
    assert len(rows(tmp_path))==8


@pytest.mark.parametrize('mutation',[
    "UPDATE words SET word='[99]' WHERE n=2 AND length=0",
    "UPDATE words SET certificate_sha='bad' WHERE n=2 AND length=0",
    "UPDATE words SET certificate=x'00' WHERE n=2 AND length=0",
    "UPDATE words SET certificate_bytes=0 WHERE n=2 AND length=0",
    "UPDATE words SET status='incomplete' WHERE n=2 AND length=0",
    "UPDATE words SET metrics='{}' WHERE n=2 AND length=0",
    "UPDATE words SET production_check='unchecked' WHERE n=2 AND length=0",
    "DELETE FROM words WHERE n=2 AND length=0",
    "UPDATE ranks SET ordinal=1 WHERE n=2",
    "UPDATE configurations SET config='{}'",
])
def test_audit_detects_corruption(tmp_path,mutation):
    run(tmp_path,targets={2:1})
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:db.execute(mutation)
    with pytest.raises(ValueError):run(tmp_path,targets={2:1},audit_only=True)


def test_identity_binding_and_adding_new_rank(tmp_path,monkeypatch):
    run(tmp_path,targets={2:1})
    extended=run(tmp_path,targets={3:1})
    assert extended['audit']['words']==3 and len(rows(tmp_path))==8
    actual=runner.identity
    def changed(batch_size,execution='collected'):
        result=actual(batch_size,execution);result['runner_sha256']='changed';return result
    monkeypatch.setattr(runner,'identity',changed)
    with pytest.raises(ValueError,match='identity changed'):run(tmp_path,targets={2:1})


def test_real_spawn_workers(tmp_path):
    result=runner.exhaust(tmp_path,targets={2:1,3:1},workers=2,wall_seconds=30,execution='collected')
    assert result['new_or_resumed_words']==8
    assert all(r['complete_through_length']==1 for r in result['ranks'])


def test_first_seen_nested_cache_counters_are_counted():
    before={'hits':0,'phase_misses':{},'actual_unit_evaluations':{}}
    after={'hits':3,'phase_misses':{'single':5},'actual_unit_evaluations':{'zip.rule':7},'label':'cache'}
    assert runner._numeric_delta(before,after)=={'hits':3,'phase_misses':{'single':5},'actual_unit_evaluations':{'zip.rule':7}}


def test_slow_rank_cannot_fill_the_entire_outstanding_window(tmp_path,monkeypatch):
    actual=runner._compute;other_rank_finished=threading.Event()
    def delayed(n,word):
        if n==2 and not word:
            assert other_rank_finished.wait(10), 'slow rank monopolized all admission slots'
        result=actual(n,word)
        if n==3 and word==[2,2]:other_rank_finished.set()
        return result
    monkeypatch.setattr(runner,'_compute',delayed)
    result=run(tmp_path,targets={2:1,3:2})
    assert result['status']=='completed' and len(rows(tmp_path))==24
    assert result['peak_in_flight']<=4


def test_deadline_stops_admissions_but_commits_already_admitted_prefix(tmp_path,monkeypatch):
    import time
    actual=runner._compute
    def delayed(n,word):
        time.sleep(.2)
        return actual(n,word)
    monkeypatch.setattr(runner,'_compute',delayed)
    result=runner.exhaust(tmp_path,targets={2:3},workers=2,wall_seconds=.1,
                          execution='collected',_executor_factory=ThreadPoolExecutor)
    assert result['status']=='paused' and result['stop_reason']=='wall_budget'
    assert result['new_or_resumed_words']==4
    assert result['new_work_seconds']>.1
    assert len(rows(tmp_path))==4
    assert run(tmp_path,targets={2:3},audit_only=True)['audit']['words']==4


def test_periodic_status_running_and_final_pause_are_explicit(tmp_path,monkeypatch):
    actual=runner._compute;seen=[]
    def observe(n,word):
        with sqlite3.connect(tmp_path/'words.sqlite3') as db:
            seen.append(json.loads(db.execute("SELECT value FROM metadata WHERE key='run_state'").fetchone()[0])['status'])
        return actual(n,word)
    monkeypatch.setattr(runner,'_compute',observe)
    result=run(tmp_path,targets={3:1},limits=Limits(max_steps=1))
    assert seen and set(seen)=={'running'}
    assert result['status']=='paused' and result['stop_reason']=='rank_limits'
    status=runner.status(tmp_path)
    assert status['recorded_run_state']['status']=='paused'
    assert status['prefix_audited_by_this_read'] is False
    assert status['process_liveness_verified'] is False


def test_audit_only_never_creates_or_updates_store(tmp_path):
    nonexistent=tmp_path/'absent'
    with pytest.raises(ValueError,match='existing database'):
        run(nonexistent,targets={2:1},audit_only=True)
    assert not nonexistent.exists()
    run(tmp_path,targets={2:1})
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        before={table:db.execute(f'SELECT * FROM {table}').fetchall() for table in ('metadata','configurations','ranks','versions')}
    result=run(tmp_path,targets={9:10},audit_only=True)
    assert result['mutated_database'] is False
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        assert before=={table:db.execute(f'SELECT * FROM {table}').fetchall() for table in before}


def _legacy_store(path):
    """Exercise the actual committed v1 writer against a small temporary store."""
    import importlib.util
    import hashlib
    from pathlib import Path
    filename=Path(__file__).parent/'fixtures/zipper_exhaustive_v1.py'
    assert hashlib.sha256(filename.read_bytes()).hexdigest()=='7c340660db34953f29467bdef4116f97611c8d47c635624d0058114b0f48bdda'
    spec=importlib.util.spec_from_file_location('legacy_test_writer',filename)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.exhaust(path,targets={2:1,3:1},workers=2,wall_seconds=30,
                   limits=Limits(max_steps=5),_executor_factory=ThreadPoolExecutor)
    with sqlite3.connect(path/'words.sqlite3') as db:
        return {table:db.execute(f'SELECT * FROM {table}').fetchall() for table in ('metadata','configurations','words','ranks')}


def test_legacy_migration_retains_identity_rows_backup_and_archives_checkpoint(tmp_path):
    old=_legacy_store(tmp_path)
    with pytest.raises(ValueError,match='explicit --migrate'):
        run(tmp_path,targets={2:1,3:1},limits=Limits(max_steps=100))
    result=run(tmp_path,targets={2:1,3:1},limits=Limits(max_steps=100),migrate=True)
    assert result['status']=='completed'
    assert result['migration']['audited_words']==len(old['words'])
    with sqlite3.connect(result['migration']['backup']) as backup:
        assert backup.execute('SELECT * FROM words').fetchall()==old['words']
        assert backup.execute('SELECT * FROM metadata').fetchall()==old['metadata']
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        assert dict(db.execute('SELECT key,value FROM metadata'))['identity']==dict(old['metadata'])['identity']
        assert set(old['configurations']).issubset(set(db.execute('SELECT * FROM configurations')))
        assert db.execute('SELECT count(*) FROM versions').fetchone()[0]==2
        archived=[row[:17] for row in db.execute('SELECT * FROM word_attempts')]
        assert set(archived)=={row for row in old['words'] if row[5]!='complete'}
        for row in old['words']:
            if row[5]=='complete':
                assert db.execute('SELECT * FROM words WHERE n=? AND length=? AND ordinal=?',row[:3]).fetchone()==row
    audited=run(tmp_path,targets={2:1,3:1},audit_only=True)
    assert audited['audit']['words']==8 and audited['audit']['archived_attempts']==2
    assert audited['audit']['hash_and_binding_checked']==10


def test_failed_version_registration_rolls_back_without_touching_old_results(tmp_path):
    old=_legacy_store(tmp_path);expected=runner.identity(64,'collected')
    class Fault:
        def __init__(self,db):self.db=db
        def __getattr__(self,name):return getattr(self.db,name)
        def __enter__(self):self.db.__enter__();return self
        def __exit__(self,*args):return self.db.__exit__(*args)
        def execute(self,sql,parameters=()):
            if sql.startswith('INSERT OR IGNORE INTO versions') and parameters[0]==runner._sha(runner._json(expected)):
                raise RuntimeError('simulated migration failure')
            return self.db.execute(sql,parameters)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        with pytest.raises(RuntimeError,match='simulated migration failure'):
            runner.register_version(Fault(db),expected,tmp_path,migrate=True)
        assert not runner._has_table(db,'versions')
        for table,records in old.items():assert db.execute(f'SELECT * FROM {table}').fetchall()==records


def test_histogram_and_archived_attempt_tampering_are_detected(tmp_path):
    _legacy_store(tmp_path)
    run(tmp_path,targets={2:1,3:1},limits=Limits(max_steps=100),migrate=True)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        db.execute('UPDATE prefix_rule_totals SET count=count+1 WHERE rowid=(SELECT min(rowid) FROM prefix_rule_totals)')
    with pytest.raises(ValueError,match='logical rule histogram'):
        run(tmp_path,targets={2:1,3:1},audit_only=True)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        db.execute('UPDATE prefix_rule_totals SET count=count-1 WHERE rowid=(SELECT min(rowid) FROM prefix_rule_totals)')
        db.execute("UPDATE word_attempts SET certificate_sha='tampered'")
    with pytest.raises(ValueError,match='hash mismatch'):
        run(tmp_path,targets={2:1,3:1},audit_only=True)


def test_actual_memoized_worker_telemetry_and_rule_histograms(tmp_path):
    result=runner.exhaust(tmp_path,targets={3:1},workers=1,wall_seconds=30,
                          execution='memoized',cache_bytes=1024*1024,cache_entries=128,
                          rule_cache_bytes=65536,rule_cache_entries=32)
    assert result['status']=='completed'
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        states=[json.loads(row[0]) for row in db.execute('SELECT runtime FROM word_runtime')]
        assert len(states)==5
        assert all(s['checker_mode']=='independent-cold-replay' for s in states)
        assert all(s['after_reduce']['max_bytes']==1024*1024 for s in states)
        assert all(s['after_reduce']['rule_cache']['max_bytes']==65536 for s in states)
        assert any(s['producer_delta']['hits']>0 for s in states)
        assert db.execute("SELECT sum(count) FROM prefix_evaluation_totals WHERE kind='actual_unit_evaluations'").fetchone()[0]>0
    assert runner.exhaust(tmp_path,audit_only=True)['audit']['words']==5


def test_incompatible_migration_replay_is_rejected_before_any_database_change(tmp_path):
    old=_legacy_store(tmp_path)
    with pytest.raises(ValueError,match='matching source checkouts'):
        runner.exhaust(tmp_path,targets={2:1,3:1},execution='memoized',migrate=True,
                       replay_audit=True,max_new_words=0)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        assert not runner._has_table(db,'versions')
        for table,records in old.items():assert db.execute(f'SELECT * FROM {table}').fetchall()==records
    assert not (tmp_path/'backups').exists()


def test_all_epoch_coverage_and_active_only_counts_are_explicit(tmp_path):
    _legacy_store(tmp_path)
    result=run(tmp_path,targets={2:1,3:1},migrate=True,max_new_words=0)
    assert result['coverage_scope']=='all_registered_versions'
    assert result['identity_scope']=='active_version'
    assert result['legacy_identity']!=result['identity']
    assert sum(rank['visited_words'] for rank in result['ranks'])==4
    assert all(rank['active_version_counts']['words']==0 for rank in result['ranks'])
    status=runner.status(tmp_path)
    assert status['coverage_scope']=='all_registered_versions'
    assert all(rank['active_version_counts']['complete_words']==0 for rank in status['ranks'])


@pytest.mark.parametrize('failure,state,reason',[(RuntimeError,'failed','failed'),(KeyboardInterrupt,'paused','interrupted')])
def test_worker_exception_records_nonrunning_state(tmp_path,monkeypatch,failure,state,reason):
    def fail(n,word):raise failure('simulated interruption')
    monkeypatch.setattr(runner,'_compute',fail)
    with pytest.raises(failure):run(tmp_path,targets={2:1})
    status=runner.status(tmp_path)
    assert status['recorded_run_state']['status']==state
    assert status['recorded_run_state']['stop_reason']==reason
    assert status['ranks'][0]['visited_words']==0


def test_empty_crash_after_schema_creation_can_restart(tmp_path):
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:runner._schema(db)
    assert run(tmp_path,targets={2:1})['status']=='completed'
    assert len(rows(tmp_path))==3


def test_nonempty_store_without_identity_is_never_adopted(tmp_path):
    run(tmp_path,targets={2:1})
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        db.execute("DELETE FROM metadata WHERE key='identity'")
        before=db.execute('SELECT * FROM words').fetchall()
    with pytest.raises(ValueError,match='unbound database contains saved records'):
        run(tmp_path,targets={2:1},migrate=True)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        assert db.execute('SELECT * FROM words').fetchall()==before
        assert not db.execute("SELECT 1 FROM metadata WHERE key='identity'").fetchone()
