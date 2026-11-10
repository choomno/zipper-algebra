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
from typing import Any
import zlib

from experiments.zipper_census_reuse_benchmark import encoded, frozen_cases
from zipper_algebra.zipping_algorithm.output_normal_form import canonical_bytes


def report(database: Path, left_run: int, right_run: int) -> dict[str, Any]:
    """Validate a same-epoch pair, oriented as frontier then DAG."""
    with sqlite3.connect(database) as db:
        db.row_factory = sqlite3.Row
        frozen_cases(db)
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
            arm = dict(run=int(run['id']), backend=run['backend'],
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
                    'store_seconds', 'session_generation', 'producer_diagnostic')}
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
        ratio = matched['frontier']['verified_call_seconds'] / matched['dag']['verified_call_seconds'] if common else None
        for arm in arms:
            del arm['rows']
        return dict(scope='single observed traversal under concurrent host load; no replicated causal estimate',
                    cache_mode=runs[0]['cache_mode'], config=configs[0],
                    source_manifest=json.loads(epochs[0]['source_manifest']),
                    driver_sha256=epochs[0]['driver_sha256'], identity=identities[0],
                    exact_matching_complete=len(common), arms=arms, matched=matched,
                    observed_frontier_over_dag_call_seconds=ratio,
                    status_cross_table=dict(Counter(o['frontier']['status'] + '/' + o['dag']['status'] for o in outcomes)),
                    outcomes=outcomes)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database', type=Path)
    parser.add_argument('runs', type=int, nargs=2)
    args = parser.parse_args()
    print(encoded(report(args.database, *args.runs)).decode())


if __name__ == '__main__':
    main()
