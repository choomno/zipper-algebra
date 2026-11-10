"""Schedule comparisons reject mixed provenance and retain censored inputs."""
import hashlib
import json
from pathlib import Path
import sqlite3
import zlib

import pytest

from experiments import zipper_census_reuse_benchmark as reuse
from experiments.zipper_schedule_report import report
from zipper_algebra.zipping_algorithm.output_normal_form import canonical_bytes, canonical_output


def database(tmp_path: Path) -> Path:
    path = tmp_path / 'schedule.sqlite3'
    identity = {'producer': {'binary_sha256': 'p', 'source_sha256': {'a': 'p'}, 'wrapper_source_sha256': 'pw'},
                'verifier': {'binary_sha256': 'v', 'source_sha256': {'b': 'v'}, 'wrapper_source_sha256': 'vw'},
                'catalog_sha256': 'c', 'comparison_sources': {'wrapper.py': 'w'}, 'harness_sha256': 'h'}
    telemetry = json.dumps({'identity': identity, 'stats': {'independently_checked': True},
                            'session_generation': 1})
    raw = canonical_bytes(canonical_output([], n=3, complete=True))
    with sqlite3.connect(path) as db:
        reuse.schema(db)
        db.execute("INSERT INTO reuse_bench_cases VALUES(1,'case',3,'[1,2]','{}')")
        reuse.frozen_cases(db, establish=True)
        db.execute('''CREATE TABLE cold_epoch_runs(run INTEGER,epoch_name TEXT,epoch_path TEXT,
          cohort TEXT,round_index INTEGER,source_manifest TEXT,driver_sha256 TEXT,
          conditions_before TEXT,conditions_after TEXT,process_seconds REAL,process_cpu_seconds REAL)''')
        for run, backend in ((1, 'dag'), (2, 'standalone')):
            db.execute("INSERT INTO reuse_bench_runs VALUES(?,?,'retained','','complete',?,?,0,NULL)",
                       (run, backend, json.dumps({'passes': 1, 'case_seconds': 20, 'limits': {}, 'cache_bytes': 256, 'cache_entries': 100, 'allocator': None, 'cpu': 6, 'process_address_space_bytes': 4096}), json.dumps(identity)))
            db.execute("INSERT INTO cold_epoch_runs VALUES(?,'epoch','/epoch','retained',0,?,'driver','{}','{}',1,1)",
                       (run, json.dumps({'src.py': 'a'})))
            db.execute("INSERT INTO reuse_bench_measurements VALUES(?,0,1,'complete',1,.2,.5,0,.1,100,?,NULL)",
                       (run, telemetry))
            db.execute('INSERT INTO reuse_bench_payloads VALUES(?,0,1,?,?,?)',
                       (run, zlib.compress(b'{"status":"complete"}'), zlib.compress(raw), hashlib.sha256(raw).hexdigest()))
    return path


def test_report_orients_by_schedule_and_keeps_full_output_parity(tmp_path):
    result = report(database(tmp_path), 1, 2)
    assert [a['backend'] for a in result['arms']] == ['standalone', 'dag']
    assert result['exact_matching_complete'] == 1
    assert result['observed_frontier_over_dag_call_seconds'] == 1
    assert result['status_cross_table'] == {'complete/complete': 1}


@pytest.mark.parametrize('change', ['failed', 'running', 'manifest', 'driver', 'identity', 'cpu', 'checked', 'hash', 'payload'])
def test_report_rejects_unbound_or_corrupted_pairs(tmp_path, change):
    path = database(tmp_path)
    with sqlite3.connect(path) as db:
        if change in ('failed', 'running'):
            db.execute('UPDATE reuse_bench_runs SET status=? WHERE id=1', (change,))
        elif change == 'manifest':
            db.execute("UPDATE cold_epoch_runs SET source_manifest='{}' WHERE run=1")
        elif change == 'driver':
            db.execute("UPDATE cold_epoch_runs SET driver_sha256='other' WHERE run=1")
        elif change == 'identity':
            db.execute("UPDATE reuse_bench_runs SET identity='{}' WHERE id=1")
        elif change == 'cpu':
            db.execute("UPDATE reuse_bench_runs SET config=? WHERE id=1", (json.dumps({'passes': 1, 'case_seconds': 20, 'cpu': 7}),))
        elif change == 'checked':
            t = json.loads(db.execute('SELECT telemetry FROM reuse_bench_measurements WHERE run=1').fetchone()[0])
            t['stats']['independently_checked'] = False
            db.execute('UPDATE reuse_bench_measurements SET telemetry=? WHERE run=1', (json.dumps(t),))
        elif change == 'hash':
            db.execute("UPDATE reuse_bench_payloads SET output_sha='bad' WHERE run=1")
        else:
            raw = canonical_bytes(canonical_output([], n=4, complete=True))
            db.execute('UPDATE reuse_bench_payloads SET output=?,output_sha=? WHERE run=1',
                       (zlib.compress(raw), hashlib.sha256(raw).hexdigest()))
    with pytest.raises((ValueError, AssertionError)):
        report(path, 1, 2)


def test_report_retains_timeout_cost_without_claiming_equal_completion(tmp_path):
    path = database(tmp_path)
    with sqlite3.connect(path) as db:
        db.execute("UPDATE reuse_bench_measurements SET status='timeout',seconds=20,producer_seconds=NULL,checker_seconds=NULL WHERE run=1")
        db.execute('DELETE FROM reuse_bench_payloads WHERE run=1')
    result = report(path, 1, 2)
    assert result['exact_matching_complete'] == 0
    assert result['observed_frontier_over_dag_call_seconds'] is None
    assert result['status_cross_table'] == {'complete/timeout': 1}
    assert result['outcomes'][0]['dag']['seconds'] == 20
    assert result['outcomes'][0]['dag']['timeout_phase'] == 'unknown'
