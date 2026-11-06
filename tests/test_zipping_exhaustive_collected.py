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
    return runner.exhaust(path,workers=2,wall_seconds=30,_executor_factory=ThreadPoolExecutor,**kwargs)


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
    assert audited['audit']=={'words':12,'hash_and_binding_checked':12,'replayed':12}


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
    assert run(tmp_path,targets={2:1,3:1})['stop_reason']=='finished'
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
    def changed(batch_size):
        result=actual(batch_size);result['runner_sha256']='changed';return result
    monkeypatch.setattr(runner,'identity',changed)
    with pytest.raises(ValueError,match='identity changed'):run(tmp_path,targets={2:1})


def test_real_spawn_workers(tmp_path):
    result=runner.exhaust(tmp_path,targets={2:1,3:1},workers=2,wall_seconds=30)
    assert result['new_or_resumed_words']==8
    assert all(r['complete_through_length']==1 for r in result['ranks'])
