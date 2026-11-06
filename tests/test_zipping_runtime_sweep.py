"""Runtime sampling keeps source versions and censored observations separate."""
import json
import sqlite3
import time

import pytest
pytest.importorskip('fcntl',reason='runtime sweep requires Unix flock')
from experiments import zipper_runtime_sweep as runtime
from zipper_algebra.zipping_algorithm.interpreter import Limits


def sleepy(sender,case,limits,identity):
    time.sleep(5)


def test_sampler_reproducible_full_support_and_permutation_words():
    cases=runtime.sample_cases([2,3,4,5,6,7,8],[4,10,12])
    assert len(cases)==81
    assert cases==runtime.sample_cases([2,3,4,5,6,7,8],[4,10,12])
    assert [c['n'] for c in cases[:7]]==list(range(2,9))
    assert all(c['length']==10 for c in cases[:7])
    assert all(len(c['word'])==c['length'] for c in cases)
    for case in cases:
        if case['family']=='negative_permutation':
            ports=list(range(case['n']));seen=set()
            for letter in reversed(case['word']):
                at=abs(letter)-1;pair=frozenset(ports[at:at+2])
                assert pair not in seen;seen.add(pair)
                ports[at],ports[at+1]=ports[at+1],ports[at]
            assert case['length']<=case['n']*(case['n']-1)//2
    assert all(c['family']!='alternating_cycle' for c in cases if c['n']==2)
    assert not any(c['n']==5 and c['length']==12 and c['family']=='negative_permutation' for c in cases)


def test_checked_child_has_cpu_wall_timings_and_timeout_has_no_certificate():
    expected=runtime.identity();case=runtime.sample_cases([2],[10],random_per_cell=0)[0]
    done=runtime.run_case(case,Limits(),expected,10)
    assert done['status']=='complete'
    assert json.loads(done['certificate'])['input_word']==case['word']
    assert done['detail']['reduce_cpu']>=0 and done['detail']['check_cpu']>=0
    assert done['attempt_wall']>=done['detail']['reduce_wall']
    timeout=runtime.run_case(case,Limits(),expected,.15,_target=sleepy)
    assert timeout['status']=='timeout' and 'certificate' not in timeout


def test_runtime_campaign_resume_and_independent_source_versions(tmp_path,monkeypatch):
    db=tmp_path/'words.sqlite3'
    with sqlite3.connect(db) as connection:
        connection.execute('CREATE TABLE metadata (key TEXT PRIMARY KEY,value TEXT)')
        connection.execute("INSERT INTO metadata VALUES ('identity','old exhaustive source untouched')")
    first=runtime.run_sweep(db,ranks=(2,),lengths=(4,),random_per_cell=0,campaign='first',max_new_cases=1)
    assert first['new_attempts']==1
    second=runtime.run_sweep(db,ranks=(2,),lengths=(4,),random_per_cell=0,campaign='first')
    assert second['new_attempts']==1
    assert all(s['completion_fraction']==1 for s in second['strata'])
    actual=runtime.identity
    def changed():
        result=actual();result['runner_sha256']='new-version';return result
    monkeypatch.setattr(runtime,'identity',changed)
    with pytest.raises(ValueError,match='source or configuration changed'):
        runtime.run_sweep(db,ranks=(2,),lengths=(4,),random_per_cell=0,campaign='first')
    third=runtime.run_sweep(db,ranks=(2,),lengths=(4,),random_per_cell=0,campaign='second',max_new_cases=0)
    assert third['new_attempts']==0
    with sqlite3.connect(db) as connection:
        assert connection.execute('SELECT count(*) FROM runtime_campaigns').fetchone()[0]==2
        assert connection.execute('SELECT value FROM metadata').fetchone()[0]=='old exhaustive source untouched'


def test_censored_results_not_complete_timing_samples(tmp_path):
    def capped(case,limits,expected,wall):
        return {'status':'timeout','attempt_wall':wall,'detail':{'phase':'reduce'}}
    result=runtime.run_sweep(tmp_path/'words.sqlite3',ranks=(3,),lengths=(10,),random_per_cell=0,_run_case=capped)
    for row in result['strata']:
        assert row['completion_fraction']==0
        assert row['complete_reduce_cpu_mean'] is None
        assert row['statuses']=={'timeout':1}


def test_saved_certificate_tamper_rejected(tmp_path):
    db=tmp_path/'words.sqlite3'
    runtime.run_sweep(db,ranks=(2,),lengths=(2,),random_per_cell=0,campaign='test')
    with sqlite3.connect(db) as connection:
        connection.execute("UPDATE runtime_cases SET certificate_sha='bad'")
    with pytest.raises(ValueError,match='hash or identity mismatch'):
        runtime.run_sweep(db,ranks=(2,),lengths=(2,),random_per_cell=0,campaign='test')


@pytest.mark.parametrize('column,value',[('n','9'),('length','99'),('family',"'forged'"),('word',"'[1,1]'"),('input_sha',"'bad'"),('full_support','3'),('ordinal','99')])
def test_saved_case_columns_are_bound(tmp_path,column,value):
    db=tmp_path/'words.sqlite3'
    runtime.run_sweep(db,ranks=(2,),lengths=(2,),random_per_cell=0,campaign='bound',max_new_cases=0)
    with sqlite3.connect(db) as connection:
        connection.execute(f'UPDATE runtime_cases SET {column}={value} WHERE ordinal=0')
    with pytest.raises(ValueError,match='deterministic sample'):
        runtime.run_sweep(db,ranks=(2,),lengths=(2,),random_per_cell=0,campaign='bound',max_new_cases=0)


def test_bounded_parallelism_and_check_timeout_classification(tmp_path):
    import threading
    lock=threading.Lock();active=0;peak=0
    def fake(case,limits,expected,wall):
        nonlocal active,peak
        with lock:active+=1;peak=max(peak,active)
        time.sleep(.02)
        with lock:active-=1
        return {'status':'timeout','attempt_wall':wall,'detail':{'phase':'check','reduce_status':'complete','reduce_wall':.1,'reduce_cpu':.1}}
    result=runtime.run_sweep(tmp_path/'words.sqlite3',ranks=(3,),lengths=(10,),workers=2,_run_case=fake)
    assert peak==2
    assert all(row['timeout_phases']=={'check':1} for row in result['strata'])
    assert all(row['complete_reduce_cpu_mean'] is None for row in result['strata'])
