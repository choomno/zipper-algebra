"""Selected probe import leaves exhaustive-prefix claims separate and intact."""
import json
import sqlite3

import pytest

pytest.importorskip("fcntl", reason="probe importer requires Unix flock")
from experiments.zipper_probe_store import import_probes
from zipper_algebra.zipping_algorithm.algorithm import AlgorithmSession, _source
from zipper_algebra.zipping_algorithm.general_enumeration import exhaust_general


def report(tmp_path):
    exhaust_general(tmp_path, ranks=(2,), max_length=0)
    cert = AlgorithmSession().reduce([1,-1],n=2)
    value = {"schema":"zipper-selected-long-word-probes-v1","exhaustive":False,"source_sha256":_source(),
        "uncheckpointed_attempts":[{"reason":"wall cap","word":[1,1,1]}],
        "probes":[{"name":"inverse_pair","n":2,"length":2,"word":[1,-1],"status":cert["status"],
                   "metrics":cert["metrics"],"limits":cert["limits"],"certificate":cert,
                   "reduce_seconds":.01,"check_seconds":.01,"json_reload_replay_checked":True}]}
    path=tmp_path/'probes.json';path.write_text(json.dumps(value));return path,value


def test_probe_import_replays_and_keeps_nonexhaustive_records_separate(tmp_path):
    path,_=report(tmp_path)
    result=import_probes(tmp_path/'words.sqlite3',path)
    assert result['probe_count']==1 and result['replay_method']=='replayed_on_import'
    import_probes(tmp_path/'words.sqlite3',path)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        assert db.execute('SELECT count(*) FROM words').fetchone()[0]==1
        assert db.execute('SELECT count(*) FROM long_word_probes').fetchone()[0]==1
        metadata=json.loads(db.execute('SELECT metadata FROM probe_imports').fetchone()[0])
        assert len(metadata['uncheckpointed_attempts'])==1
    assert exhaust_general(tmp_path,ranks=(2,),max_length=0)['audited_words']==1


@pytest.mark.parametrize('mutation',['source','certificate','declaration'])
def test_bad_source_or_certificate_is_rejected_before_any_probe_rows(tmp_path,mutation):
    path,value=report(tmp_path)
    if mutation=='source':value['source_sha256']['algorithm.py']='bad'
    elif mutation=='certificate':value['probes'][0]['certificate']['n']=3
    else:value['probes'][0]['json_reload_replay_checked']=False
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError):import_probes(tmp_path/'words.sqlite3',path,recorded_replay=True)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        assert not db.execute("SELECT name FROM sqlite_master WHERE name='long_word_probes'").fetchall()


def test_invalid_second_probe_cannot_partially_import_valid_first(tmp_path):
    from copy import deepcopy
    path,value=report(tmp_path)
    invalid=deepcopy(value['probes'][0]);invalid['name']='invalid_second'
    invalid['certificate']['limits']['max_steps']+=1
    # Separate summary limits from the certificate before changing its limit.
    invalid['limits']=deepcopy(value['probes'][0]['limits'])
    value['probes'].append(invalid);path.write_text(json.dumps(value))
    with pytest.raises(ValueError,match='identity or summary mismatch'):
        import_probes(tmp_path/'words.sqlite3',path,recorded_replay=True)
    with sqlite3.connect(tmp_path/'words.sqlite3') as db:
        assert not db.execute("SELECT name FROM sqlite_master WHERE name='long_word_probes'").fetchall()
        assert db.execute('SELECT count(*) FROM words').fetchone()[0]==1
