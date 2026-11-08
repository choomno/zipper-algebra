from copy import deepcopy
from dataclasses import asdict
import hashlib
import sqlite3
import zlib

import pytest

from zipper_algebra.zipping_algorithm import census_reuse as cr
from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.input_normal_form import input_normal_form, leaf_requests
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import (
    completed_output, execution_identity, reconstruct_output,
)
from zipper_algebra.zipping_algorithm.output_normal_form import from_certificate
from zipper_algebra.zipping_algorithm.reconstruction_check import check_reconstruction


@pytest.fixture
def backend():
    session = CollectedSession()
    limits = Limits()
    identity = execution_identity(session.reduce([], n=1, limits=limits))
    calls = []
    def produce(n, word):
        calls.append((n, tuple(word)))
        return session.reduce(word, n=n, limits=limits)
    checks = []
    def verify(cert):
        checks.append(cert['input_word'])
        return session.check(cert)
    return dict(produce=produce, verify=verify, producer_identity=identity,
                verifier_identity={'name': 'independent-python-test'}, limits=asdict(limits)), calls, checks, session


@pytest.mark.parametrize('word,n', [([], 1), ([], 4), ([1,-1],2), ([-1,-1],2),
    ([1,3],4), ([1,2,1,5],6), ([1,-3,1,3,1,-3,2,1],4)])
def test_exact_reconstruction_and_reuse(tmp_path, backend, word, n):
    settings, calls, checks, session = backend
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3', **settings) as store:
        result = store.evaluate(word,n=n)
        assert result['output'] == from_certificate(session.reduce(word,n=n))
        count = len(calls)
        repeated = store.evaluate(word,n=n)
        assert repeated['output'] == result['output'] and len(calls) == count
        assert len(checks) == count
        assert store.counts() == {'roots':1,'leaves':count,'attempts':0}
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3', **settings) as store:
        assert store.evaluate(word,n=n)['output'] == result['output']
        assert len(calls) == count and len(checks) == 2*count


def test_complete_only_and_rejection(tmp_path, backend):
    settings, _, _, _ = backend
    original = settings['produce']
    def incomplete(n, word):
        value = original(n,word)
        value['status'] = 'incomplete'
        return value
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3', **dict(settings,produce=incomplete)) as store:
        assert store.evaluate([1],n=2)['status'] == 'incomplete'
        assert store.counts() == {'roots':0,'leaves':0,'attempts':1}
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3', **dict(settings, verify=lambda _: False)) as store:
        with pytest.raises(ValueError,match='verifier rejected'):
            store.evaluate([1],n=2)
        assert store.counts()['leaves'] == 0


def test_corrupt_payload_and_dependencies(tmp_path, backend):
    settings, _, _, _ = backend
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3', **settings) as store:
        store.evaluate([1,3],n=4)
        row = store.db.execute('SELECT payload FROM reuse_leaves').fetchone()[0]
        store.db.execute("UPDATE reuse_leaves SET digest='incorrect'")
        store.db.commit()
        with pytest.raises(ValueError,match='digest mismatch'):
            store.evaluate([1,3],n=4)
        store.db.execute('UPDATE reuse_leaves SET digest=?', (hashlib.sha256(zlib.decompress(row)).hexdigest(),))
        store.db.execute('DELETE FROM reuse_dependencies')
        store.db.commit()
        with pytest.raises(ValueError,match='dependency mismatch'):
            store.evaluate([1,3],n=4)


def test_reconstruction_checker_mutations(backend):
    settings, _, _, _ = backend
    form = input_normal_form([1,3],n=4)
    leaves = {r['key']:settings['produce'](r['n'],r['word']) for r in leaf_requests(form)}
    output = reconstruct_output(form,leaves,source_identity=settings['producer_identity'])
    assert check_reconstruction(form,leaves,output,source_identity=settings['producer_identity'])
    bad = deepcopy(output)
    bad['terms'][0]['coefficient'][0][1] = '100'
    assert not check_reconstruction(form,leaves,bad,source_identity=settings['producer_identity'])
    assert not check_reconstruction({},leaves,output,source_identity=settings['producer_identity'])
    assert not check_reconstruction(form,{},output,source_identity=settings['producer_identity'])
    assert not check_reconstruction(form,dict(leaves,extra={}),output,source_identity=settings['producer_identity'])
    bad = deepcopy(leaves)
    next(iter(bad.values()))['input_word'] = [-1]
    assert not check_reconstruction(form,bad,output,source_identity=settings['producer_identity'])
    bad = deepcopy(leaves)
    next(iter(bad.values()))['status'] = 'incomplete'
    assert not check_reconstruction(form,bad,output,source_identity=settings['producer_identity'])


def test_policy_isolation_and_legacy_untouched(tmp_path, backend):
    settings, _, _, _ = backend
    path=tmp_path/'legacy.sqlite3'
    db=sqlite3.connect(path)
    db.execute('CREATE TABLE words(value TEXT)')
    db.execute("INSERT INTO words VALUES ('original')")
    db.commit(); db.close()
    before=path.read_bytes()
    with pytest.raises(ValueError,match='sidecar'):
        cr.ReuseCensus(path,**settings)
    assert path.read_bytes() == before
    for identity in ({'verifier':'a'},{'verifier':'b'}):
        with cr.ReuseCensus(tmp_path/'reuse.sqlite3',**dict(settings,verifier_identity=identity)) as store:
            assert store.evaluate([1],n=2)['stats']['leaf_misses'] == 1
    for kwargs in ({'cache_entries':-1},{'verifier_identity':{}},{'limits':{}}):
        with pytest.raises(ValueError):
            cr.ReuseCensus(tmp_path/'bad.sqlite3',**dict(settings,**kwargs))


def test_lru_and_source_drift(tmp_path,backend,monkeypatch):
    settings, _, checks, _ = backend
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3',**settings,cache_entries=0) as store:
        store.evaluate([1],n=2); store.evaluate([1],n=2)
        assert len(checks)==2 and not store._checked
        monkeypatch.setattr(cr,'reconstruction_sources',lambda:{})
        with pytest.raises(ValueError,match='sources changed'):
            store.evaluate([1],n=2)


def test_native_reuse_adapter(backend):
    settings, _, _, _ = backend
    cert=settings['produce'](2,[1])
    expected=from_certificate(cert)
    cert.update(schema='zipper-native-reuse-v1',
                scope='literal_word_native_reused_spanning_procedure',trace_enabled=False,
                execution={'schedule':'intrinsic-completion-dag-v1'})
    assert completed_output(cert)==expected
    cert['trace_enabled']=True
    with pytest.raises(ValueError,match='native reuse'):
        completed_output(cert)


def test_more_corruptions(tmp_path,backend):
    settings, _, _, _ = backend
    path=tmp_path/'reuse.sqlite3'
    with cr.ReuseCensus(path,**settings) as store:
        store.evaluate([1],n=2)
        original=store.db.execute('SELECT output_digest FROM reuse_leaves').fetchone()[0]
        store.db.execute("UPDATE reuse_leaves SET output_digest='bad'");store.db.commit()
        with pytest.raises(ValueError,match='output digest'):
            store.evaluate([1],n=2)
        store.db.execute('UPDATE reuse_leaves SET output_digest=?',(original,))
        store.db.execute("UPDATE reuse_roots SET output_digest='bad'");store.db.commit()
        with pytest.raises(ValueError,match='immutable reconstruction'):
            store.evaluate([1],n=2)
        store.db.execute("UPDATE reuse_policies SET payload='{}'");store.db.commit()
    with pytest.raises(ValueError,match='policy identity'):
        cr.ReuseCensus(path,**settings)


def test_verifier_resource_refusal_and_stored_rejection(tmp_path, backend):
    from zipper_algebra.zipping_algorithm.native_verifier import VerifierIncomplete
    settings, _, _, _ = backend
    path=tmp_path/'reuse.sqlite3'
    with cr.ReuseCensus(path,**settings) as store:
        store.evaluate([1],n=2)
    def limited(_):
        raise VerifierIncomplete('budget',replay_diagnostic={'limit':'steps'})
    with cr.ReuseCensus(path,**dict(settings,verify=limited)) as store:
        result=store.evaluate([1],n=2)
        assert result['status']=='incomplete' and result['diagnostic']['stage']=='verification'
        assert store.counts()=={'leaves':1,'roots':1,'attempts':1}
    with cr.ReuseCensus(path,**dict(settings,verify=lambda _:False)) as store:
        with pytest.raises(ValueError,match='stored primitive'):
            store.evaluate([1],n=2)


def test_reconstruction_bounds_and_check_failure(tmp_path,backend,monkeypatch):
    settings, _, _, _ = backend
    form=input_normal_form([1,3],n=4)
    leaves={r['key']:settings['produce'](r['n'],r['word']) for r in leaf_requests(form)}
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3',**settings) as store:
        for key,value,message in [('max_terms',1,'Cartesian'),('max_coefficient_support',1,'support'),
                                  ('max_word_size',1,'word budget'),('max_integer_bits',1,'integer budget')]:
            original=store.limits[key]
            store.limits[key]=value
            assert message in store._reconstruction_budget(form,leaves)
            store.limits[key]=original
        monkeypatch.setattr(store,'_reconstruction_budget',lambda *args:'budget')
        result=store.evaluate([1,3],n=4)
        assert result['status']=='incomplete' and store.counts()=={'leaves':1,'roots':0,'attempts':1}
        monkeypatch.setattr(store,'_reconstruction_budget',lambda *args:None)
        monkeypatch.setattr(cr,'check_reconstruction',lambda *args,**kwargs:False)
        with pytest.raises(ValueError,match='reconstruction check'):
            store.evaluate([1,3],n=4)
        original=cr.reconstruct_output
        store.limits['max_integer_bits']=1
        # Test actual output check with a leaf-free identity and a faulty coefficient.
        monkeypatch.setattr(cr,'reconstruct_output',lambda *args,**kwargs:
                            {'schema':'zipper-formal-output-normal-form-v1','coefficient_ring':'Z[q,q^-1]',
                             'n':1,'completion':'complete','terms':[
                                 dict(original(input_normal_form([],n=1),{},source_identity=settings['producer_identity'])['terms'][0],
                                      coefficient=[['0','100000']])]})
        assert store.evaluate([],n=1)['status']=='incomplete'


def test_invalid_binding_and_verifier_identity(tmp_path,backend):
    settings, _, _, _ = backend
    binding=cr.verifier_binding({'abi':'1','source_sha256':{'x':'a'},'binary_sha256':'b','initialization_seconds':99})
    assert binding=={'abi':'1','source_sha256':{'x':'a'},'binary_sha256':'b'}
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3',**settings) as store:
        request=leaf_requests(input_normal_form([1],n=2))[0]
        result=settings['produce'](2,[1])
        for changes in ({'n':True},{'input_word':[True]},{'input_word':[-1]},
                        {'limits':{}},{'status':'incomplete'}):
            with pytest.raises(ValueError,match='binding'):
                store._validate_leaf(dict(result,**changes),request)


@pytest.mark.parametrize('conflict',[False,True])
def test_concurrent_leaf_insertion_is_exact(tmp_path,backend,conflict):
    from zipper_algebra.zipping_algorithm.normal_form_reconstruction import execution_cache_key
    from zipper_algebra.zipping_algorithm.output_normal_form import output_fingerprint
    settings, _, _, _ = backend
    request=leaf_requests(input_normal_form([1],n=2))[0]
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3',**settings) as store:
        def racing_verifier(result):
            assert settings['verify'](result)
            encoded=cr._json(result).encode()
            key=execution_cache_key(request,settings['producer_identity'])
            with store.db:
                store.db.execute('INSERT INTO reuse_leaves VALUES (?,?,?,?,?,?)',
                                 (store.policy_key,key,cr._json(request),zlib.compress(encoded),
                                  'wrong' if conflict else hashlib.sha256(encoded).hexdigest(),
                                  output_fingerprint(completed_output(result))))
            return True
        store.verify=racing_verifier
        if conflict:
            with pytest.raises(ValueError,match='immutable primitive'):
                store.evaluate([1],n=2)
            assert store.counts()['roots']==0
        else:
            assert store.evaluate([1],n=2)['status']=='complete'
            assert store.counts()['leaves']==1


def test_identity_embedding_fastpath_and_tampering(tmp_path, backend, monkeypatch):
    from zipper_algebra.zipping_algorithm.reconstruction_check import identity_leaf_request
    settings, _, _, session = backend
    form=input_normal_form([1,2,1],n=3)
    requests=leaf_requests(form)
    assert identity_leaf_request(form,requests)==requests[0]['key']
    for bad in ({},dict(form,root=dict(form['root'],top_labels=[2,1,3])),
                dict(form,root=dict(form['root'],bottom_labels=[2,1,3]))):
        assert identity_leaf_request(bad,requests) is None
    assert identity_leaf_request(form,[]) is None
    assert identity_leaf_request(form,[dict(requests[0],n=True)]) is None
    split=input_normal_form([1,3],n=4)
    assert identity_leaf_request(split,leaf_requests(split)) is None
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3',**settings) as store:
        def forbidden(*args,**kwargs):
            raise AssertionError('identity embedding should not assemble an output')
        monkeypatch.setattr(cr,'reconstruct_output',forbidden)
        monkeypatch.setattr(cr,'check_reconstruction',forbidden)
        result=store.evaluate([1,2,1],n=3)
        assert result['stats']['identity_reconstruction'] is True
        assert result['output']==from_certificate(session.reduce([1,2,1],n=3))


def test_primitive_normalizes_once_and_reuses_digest(tmp_path,backend,monkeypatch):
    from zipper_algebra.zipping_algorithm.output_normal_form import output_fingerprint
    settings, _, _, _=backend
    original=cr.completed_output
    calls=[]
    def count(value):
        calls.append(value['input_word'])
        return original(value)
    monkeypatch.setattr(cr,'completed_output',count)
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3',**settings) as store:
        first=store.evaluate([1,2,1],n=3)
        assert calls==[[1,2,1]]
        assert first['output_digest']==output_fingerprint(first['output'])
        first['output']['terms'].clear()  # Caller mutation cannot poison saved data.
        second=store.evaluate([1,2,1],n=3)
        assert len(calls)==2 and second['output']['terms']
        assert second['output_digest']==output_fingerprint(second['output'])


def test_source_request_binding_and_store_lock(tmp_path,backend,monkeypatch):
    settings, _, _, _=backend
    path=tmp_path/'reuse.sqlite3'
    with cr.ReuseCensus(path,**settings) as store:
        with pytest.raises(BlockingIOError):
            cr.ReuseCensus(path,**settings)
        different=input_normal_form([1,2,1],n=3)
        monkeypatch.setattr(cr,'input_normal_form',lambda *args,**kwargs:different)
        with pytest.raises(ValueError,match='requested literal'):
            store.evaluate([2,1,2],n=3)
        assert store.counts()=={'leaves':0,'roots':0,'attempts':0}
    with pytest.raises(ValueError,match='authenticated verifier'):
        cr.verifier_binding({'abi':'x'})
    monkeypatch.setattr(cr,'reconstruction_sources',lambda:{})
    with pytest.raises(ValueError,match='since import'):
        cr.ReuseCensus(path,**settings)


def test_legacy_uncheckpointed_wal_is_read_only(tmp_path,backend):
    import subprocess
    import sys
    settings, _, _, _=backend
    path=tmp_path/'legacy.sqlite3'
    script="""import os,sqlite3,sys
c=sqlite3.connect(sys.argv[1]);c.execute('PRAGMA journal_mode=WAL')
c.execute('CREATE TABLE words(value)');c.execute('INSERT INTO words VALUES (17)');c.commit();os._exit(0)
"""
    subprocess.run([sys.executable,'-c',script,str(path)],check=True)
    files=[path,path.with_name(path.name+'-wal')]
    before=[p.read_bytes() for p in files]
    with pytest.raises(ValueError,match='sidecar'):
        cr.ReuseCensus(path,**settings)
    assert [p.read_bytes() for p in files]==before


def test_rehashed_wrong_leaf_still_requires_independent_replay(tmp_path,backend):
    import json
    from zipper_algebra.zipping_algorithm.native import _decode_coefficients
    from zipper_algebra.zipping_algorithm.output_normal_form import output_fingerprint
    settings, _, checks, _=backend
    with cr.ReuseCensus(tmp_path/'reuse.sqlite3',**settings) as store:
        store.evaluate([1],n=2)
        row=store.db.execute('SELECT payload FROM reuse_leaves').fetchone()[0]
        value=_decode_coefficients(json.loads(zlib.decompress(row)))
        value['output'][0]['coefficient']=[[p,2*c] for p,c in value['output'][0]['coefficient']]
        encoded=cr._json(value).encode()
        with store.db:
            store.db.execute('UPDATE reuse_leaves SET payload=?,digest=?,output_digest=?',
                             (zlib.compress(encoded),hashlib.sha256(encoded).hexdigest(),
                              output_fingerprint(completed_output(value))))
        with pytest.raises(ValueError,match='deterministic replay'):
            store.evaluate([1],n=2)
        assert len(checks)==2


@pytest.mark.parametrize('failure', ['connect', 'legacy-race'])
def test_initialization_failure_releases_lock_and_preserves_legacy(tmp_path, backend,
                                                                  monkeypatch, failure):
    """A failed open or file created after preflight cannot strand the store lock."""
    import fcntl
    settings, _, _, _ = backend
    path = tmp_path / 'racing.sqlite3'
    connect = sqlite3.connect
    opened = []

    def injected(*args, **kwargs):
        if failure == 'connect':
            raise sqlite3.OperationalError('test connection failure')
        connection = connect(*args, **kwargs)
        connection.execute('CREATE TABLE words(value)')
        connection.execute('INSERT INTO words VALUES (31)')
        connection.commit()
        opened.append(connection)
        return connection

    monkeypatch.setattr(cr.sqlite3, 'connect', injected)
    error = sqlite3.OperationalError if failure == 'connect' else ValueError
    with pytest.raises(error, match='connection failure|sidecar'):
        cr.ReuseCensus(path, **settings)
    with path.with_name(path.name + '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if opened:
        with pytest.raises(sqlite3.ProgrammingError, match='closed'):
            opened[0].execute('SELECT 1')
        with connect(path) as db:
            assert db.execute('SELECT value FROM words').fetchall() == [(31,)]
            assert db.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall() == [('words',)]
