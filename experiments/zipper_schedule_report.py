"""Strict same-source frontier/DAG reports with complete and censored outcomes.

Reports describe an observed traversal. They do not estimate a causal speedup
or production throughput. Native startup is separate from stored-result time.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import statistics
from typing import Any
import zlib

from experiments.zipper_census_reuse_benchmark import encoded, frozen_cases
from zipper_algebra.zipping_algorithm.output_normal_form import canonical_bytes


def file_sha256(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def report(database: Path, left_run: int, right_run: int, *, cases_file: Path | None = None) -> dict[str, Any]:
    """Validate a same-epoch pair, oriented as frontier then DAG."""
    wal = Path(str(database) + '-wal')
    if wal.exists() and wal.stat().st_size:
        raise ValueError('report requires a consolidated database without pending WAL pages')
    database_sha = file_sha256(database)
    reporter_path = Path(__file__)
    reporter_sha = hashlib.sha256(reporter_path.read_bytes()).hexdigest()
    tests_path = reporter_path.parents[1] / 'tests/test_zipping_schedule_report.py'
    tests_sha = hashlib.sha256(tests_path.read_bytes()).hexdigest() if tests_path.exists() else None
    with sqlite3.connect(database.resolve().as_uri() + '?mode=ro', uri=True) as db:
        db.row_factory = sqlite3.Row
        db.execute('BEGIN')
        frozen = frozen_cases(db)
        cases_sha = None
        if cases_file is not None:
            raw_cases = cases_file.read_bytes()
            if json.loads(raw_cases) != frozen:
                raise ValueError('cases file differs from frozen database requests')
            cases_sha = hashlib.sha256(raw_cases).hexdigest()
        runs = [db.execute('SELECT * FROM reuse_bench_runs WHERE id=?', (r,)).fetchone()
                for r in (left_run, right_run)]
        if any(r is None or r['status'] != 'complete' for r in runs):
            raise ValueError('both schedule runs must be complete')
        runs.sort(key=lambda r: r['backend'] != 'standalone')
        if [r['backend'] for r in runs] != ['standalone', 'dag']:
            raise ValueError('one frontier and one DAG run required')
        configs = [json.loads(r['config']) for r in runs]
        required = ('passes', 'case_seconds', 'limits', 'cache_bytes', 'cache_entries', 'allocator', 'cpu', 'process_address_space_bytes')
        if any(any(k not in c for k in required) for c in configs):
            raise ValueError('resource and cache configuration is incomplete')
        if configs[0] != configs[1] or configs[0]['passes'] != 1:
            raise ValueError('identical one-pass configurations required')
        if runs[0]['cache_mode'] != runs[1]['cache_mode']:
            raise ValueError('cache regimes differ')
        epochs = [db.execute('SELECT * FROM cold_epoch_runs WHERE run=?', (r['id'],)).fetchone()
                  for r in runs]
        if any(e is None or e['conditions_after'] is None or e['process_seconds'] is None
               for e in epochs):
            raise ValueError('finished epoch provenance required')
        for key in ('epoch_name', 'epoch_path', 'cohort', 'round_index', 'source_manifest', 'driver_sha256'):
            if epochs[0][key] != epochs[1][key]:
                raise ValueError('schedule epoch binding differs: ' + key)
        if not json.loads(epochs[0]['source_manifest']) or not epochs[0]['driver_sha256']:
            raise ValueError('full source and driver hashes required')
        identities = [json.loads(r['identity']) for r in runs]
        for identity in identities:
            if not all(identity.get(k) for k in ('producer', 'verifier', 'catalog_sha256', 'comparison_sources', 'harness_sha256')):
                raise ValueError('complete native and interface provenance required')
            for implementation in ('producer', 'verifier'):
                if not all(identity[implementation].get(k) for k in ('binary_sha256', 'source_sha256', 'wrapper_source_sha256')):
                    raise ValueError('native binary and source hashes required')
        if identities[0] != identities[1]:
            raise ValueError('native build, wrapper, catalog or allocator identities differ')
        cases = {r['id']: dict(r) for r in db.execute('SELECT * FROM reuse_bench_cases')}
        labels = configs[0].get('labels')
        expected = {k for k, c in cases.items() if not labels or c['label'] in labels}
        arms = []
        for run, epoch, identity in zip(runs, epochs, identities):
            rows = {r['case_id']: dict(r) for r in db.execute(
                'SELECT * FROM reuse_bench_measurements WHERE run=? AND pass=0', (run['id'],))}
            total = db.execute('SELECT count(*) FROM reuse_bench_measurements WHERE run=?', (run['id'],)).fetchone()[0]
            if set(rows) != expected or total != len(rows):
                raise ValueError('schedule measurements omit or add selected inputs')
            for case_id, row in rows.items():
                telemetry = json.loads(row['telemetry'])
                if telemetry.get('identity') != identity:
                    raise ValueError('per-measurement identity differs from run identity')
                row['session_generation'] = telemetry.get('session_generation')
                row['cache'] = telemetry.get('cache')
                row['independently_checked'] = telemetry.get('stats', {}).get('independently_checked', False)
                row['output_sha256'] = None
                payload = db.execute('SELECT * FROM reuse_bench_payloads WHERE run=? AND pass=0 AND case_id=?',
                                     (run['id'], case_id)).fetchone()
                row['producer_diagnostic'] = None
                if payload is not None and payload['result'] is not None:
                    result = json.loads(zlib.decompress(payload['result']))
                    row['producer_diagnostic'] = result.get('producer_result', result).get('diagnostic')
                if row['status'] == 'complete':
                    for field in ('seconds', 'producer_seconds', 'checker_seconds', 'store_seconds'):
                        value = row[field]
                        if value is None or not math.isfinite(value) or value < 0:
                            raise ValueError('completed output lacks a valid stage duration')
                    if telemetry.get('stats', {}).get('independently_checked') is not True:
                        raise ValueError('completed output lacks independent acceptance')
                    if payload is None or payload['output'] is None:
                        raise ValueError('completed output absent')
                    raw = zlib.decompress(payload['output'])
                    normal = json.loads(raw)
                    if (hashlib.sha256(raw).hexdigest() != payload['output_sha']
                            or canonical_bytes(normal) != raw or normal.get('completion') != 'complete'):
                        raise ValueError('completed output integrity or canonicalization differs')
                    row['output_bytes'] = raw
                    row['output_sha256'] = payload['output_sha']
            verification_scope = ('independent full standalone certificate replay, including execution fields'
                                  if run['backend'] == 'standalone' else
                                  'independent full formal-output replay; standalone schedule, trace and producer counters are not certified')
            arm = dict(run=int(run['id']), backend=run['backend'], started_utc=run['started_utc'],
                       verification_scope=verification_scope,
                       counts=dict(Counter(r['status'] for r in rows.values())),
                       measured_seconds=sum(r['seconds'] for r in rows.values()),
                       verified_call_seconds=(sum(r['producer_seconds'] + r['checker_seconds'] for r in rows.values())
                                              if all(r['producer_seconds'] is not None and r['checker_seconds'] is not None
                                                     for r in rows.values()) else None),
                       recorded_stage_seconds=sum((r['producer_seconds'] or 0) + (r['checker_seconds'] or 0)
                                                  for r in rows.values()),
                       calls_with_unknown_durations=sum(r['producer_seconds'] is None or r['checker_seconds'] is None
                                                        for r in rows.values()),
                       store_seconds=sum(r['store_seconds'] or 0 for r in rows.values()),
                       setup_seconds=run['setup_seconds'], process_seconds=epoch['process_seconds'],
                       process_cpu_seconds=epoch['process_cpu_seconds'],
                       sampled_peak_tree_rss_bytes=max(r['peak_tree_rss_bytes'] or 0 for r in rows.values()),
                       session_generations=sorted({r['session_generation'] for r in rows.values()}),
                       conditions_before=json.loads(epoch['conditions_before']),
                       conditions_after=json.loads(epoch['conditions_after']), rows=rows)
            arms.append(arm)
        outcomes = []
        for case_id in sorted(expected):
            l, r = (a['rows'][case_id] for a in arms)
            if l['status'] == r['status'] == 'complete' and l['output_bytes'] != r['output_bytes']:
                raise AssertionError('full canonical output differs for case ' + str(case_id))
            c = cases[case_id]
            outcome = dict(case_id=case_id, label=c['label'], n=c['n'], length=len(json.loads(c['word'])))
            for name, row in zip(('frontier', 'dag'), (l, r)):
                outcome[name] = {k: row[k] for k in ('status', 'seconds', 'producer_seconds', 'checker_seconds',
                    'store_seconds', 'session_generation', 'producer_diagnostic', 'output_sha256', 'cache', 'independently_checked')}
                outcome[name]['verified_call_seconds'] = (row['producer_seconds'] + row['checker_seconds']
                                                         if row['producer_seconds'] is not None and row['checker_seconds'] is not None
                                                         else None)
                outcome[name]['timeout_phase'] = 'unknown' if row['status'] == 'timeout' else None
                outcome[name]['wall_deadline_seconds'] = configs[0]['case_seconds']
            outcomes.append(outcome)
        common = [o for o in outcomes if o['frontier']['status'] == o['dag']['status'] == 'complete']
        matched = {name: {key: sum(o[name][key] or 0 for o in common)
                          for key in ('seconds', 'producer_seconds', 'checker_seconds', 'verified_call_seconds', 'store_seconds')}
                   for name in ('frontier', 'dag')}
        for outcome in common:
            denominator = outcome['dag']['verified_call_seconds']
            outcome['observed_frontier_over_dag_call_seconds'] = (outcome['frontier']['verified_call_seconds'] / denominator
                                                                 if denominator > 0 else None)
        ratios = sorted(o['observed_frontier_over_dag_call_seconds'] for o in common
                        if o['observed_frontier_over_dag_call_seconds'] is not None)
        distribution = (dict(count=len(ratios), minimum=ratios[0], median=statistics.median(ratios), maximum=ratios[-1],
                             p10=ratios[max(0, math.ceil(.1 * len(ratios)) - 1)],
                             p90=ratios[max(0, math.ceil(.9 * len(ratios)) - 1)],
                             quantile_method='nearest rank; descriptive empirical distribution') if ratios else None)
        blocks = []
        for start in range(0, len(outcomes), 32):
            block = outcomes[start:start + 32]
            completed = [o for o in block if o['frontier']['status'] == o['dag']['status'] == 'complete']
            times = {name: sum(o[name]['verified_call_seconds'] for o in completed) for name in ('frontier', 'dag')}
            blocks.append(dict(first_case_id=block[0]['case_id'], last_case_id=block[-1]['case_id'],
                               common_complete=len(completed), producer_checker_seconds=times,
                               observed_frontier_over_dag_call_seconds=(times['frontier'] / times['dag'] if times['dag'] > 0 else None),
                               cache_end={name: block[-1][name]['cache'] for name in ('frontier', 'dag')}))
        output_manifest = [{'case_id': o['case_id'], 'output_sha256': o['frontier']['output_sha256']} for o in common]
        if file_sha256(database) != database_sha:
            raise ValueError('database changed during report generation')
        if hashlib.sha256(reporter_path.read_bytes()).hexdigest() != reporter_sha:
            raise ValueError('report generator changed during execution')
        denominator = matched['dag']['verified_call_seconds']
        ratio = matched['frontier']['verified_call_seconds'] / denominator if common and denominator > 0 else None
        for arm in arms:
            del arm['rows']
        return dict(scope='single observed traversal under concurrent host load; no replicated causal estimate',
                    database_name=database.name, database_file_sha256=database_sha,
                    reporter_sha256=reporter_sha, report_tests_sha256=tests_sha, cases_file_sha256=cases_sha,
                    frozen_cases_sha256=hashlib.sha256(encoded(frozen)).hexdigest(),
                    complete_output_manifest_sha256=hashlib.sha256(encoded(output_manifest)).hexdigest(),
                    other_runs=[dict(r) for r in db.execute('SELECT id,backend,cache_mode,status FROM reuse_bench_runs WHERE id NOT IN (?,?)', (left_run, right_run))],
                    cache_mode=runs[0]['cache_mode'], config=configs[0],
                    source_manifest=json.loads(epochs[0]['source_manifest']),
                    driver_sha256=epochs[0]['driver_sha256'], identity=identities[0],
                    exact_matching_complete=len(common), arms=arms, matched=matched,
                    blocks=blocks, per_case_call_ratio_distribution=distribution,
                    observed_frontier_over_dag_call_seconds=ratio,
                    status_cross_table=dict(Counter(o['frontier']['status'] + '/' + o['dag']['status'] for o in outcomes)),
                    outcomes=outcomes)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database', type=Path)
    parser.add_argument('runs', type=int, nargs=2)
    parser.add_argument('--cases-file', type=Path)
    args = parser.parse_args()
    print(encoded(report(args.database, *args.runs, cases_file=args.cases_file)).decode())


if __name__ == '__main__':
    main()
