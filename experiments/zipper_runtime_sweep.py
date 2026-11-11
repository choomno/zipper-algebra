"""Selected, source-bound runtime sweeps in one SQLite database (Unix).

These samples are not exhaustive. Each word runs in a fresh spawned process;
its wall cap includes startup, reduction, JSON reload/replay and transfer.
Bounded parent threads supervise independent children; CPU times remain useful
when wall times reflect concurrent load. Worker count is campaign metadata.
The campaign wall budget stops new admissions; active cases finish or reach
their individual caps, so total time can exceed it by a case cap plus cleanup.
Reduction/check timings exclude child startup and catalog setup, recorded
separately. CPU times use process_time for reference execution. Collected execution adds
reaped term-worker CPU time; its default is one outer word worker, and
--term-workers controls parallel arithmetic inside that word. Timed-out words
terminate their entire process group, including nested workers.
Incomplete and timeout observations
are censored, never substituted for complete-word timings.

New source versions receive independent runtime campaigns without changing
any existing exhaustive-prefix tables. A single parent holds the campaign
flock and writes SQLite; workers never open the database. Resume verifies
saved hashes and identities and runs pending cases only. Every stored
certificate was JSON-reloaded and deterministically replay-checked before
its original commit; resume does not repeat all old computations.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from dataclasses import asdict
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import multiprocessing
import os
import resource
import signal
from pathlib import Path
import platform
import random
import sqlite3
import statistics
import subprocess
import sys
import time

from zipper_algebra.zipping_algorithm.algorithm import AlgorithmSession, _source
from zipper_algebra.zipping_algorithm.catalog import Catalog, SCHEMA_SHA256, _data
from zipper_algebra.zipping_algorithm.enumeration import _json, _sha
from zipper_algebra.zipping_algorithm.interpreter import Limits


def _algorithm_source(execution):
    if execution == "collected":
        from zipper_algebra.zipping_algorithm.collected import _source as collected_source
        return collected_source()
    return _source()


def identity(execution="reference"):
    return {"algorithm_source_sha256": _algorithm_source(execution), "catalog_sha256": Catalog.bundled().sha256,
            "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "execution": execution,
            "catalog_schema_sha256": SCHEMA_SHA256,
            "coefficient_ring": "Z[q,q^-1]", "python": sys.version, "platform": platform.platform()}


def _source_matches(expected):
    return (_algorithm_source(expected.get("execution", "reference")) == expected["algorithm_source_sha256"]
            and hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == expected["runner_sha256"]
            and hashlib.sha256(_data("full-catalog.json")).hexdigest() == expected["catalog_sha256"]
            and hashlib.sha256(_data("schema.json")).hexdigest() == expected["catalog_schema_sha256"])


def sample_cases(ranks, lengths, *, seed=20261106, random_per_cell=1):
    """Interleave ranks within strata, visiting L=10 first when requested."""
    if (not ranks or not lengths or len(set(ranks)) != len(ranks) or len(set(lengths)) != len(lengths)
            or any(type(n) is not int or n < 2 for n in ranks)
            or any(type(length) is not int or length < 1 for length in lengths)
            or type(seed) is not int or type(random_per_cell) is not int or random_per_cell < 0):
        raise ValueError("distinct ranks >=2, positive lengths, integer seed and nonnegative sample count required")
    result = []
    for length in sorted(lengths, key=lambda value: (value != 10, value)):
        cells = {}
        for n in ranks:
            positive = [1 + p % (n - 1) for p in range(length)]
            delta = [-x for width in range(n - 1, 0, -1) for x in range(1, width + 1)]
            negative = ("negative_permutation", delta[:length]) if length <= len(delta) else (
                "negative_cycle", [-x for x in positive])
            families = [negative, ("positive_cycle", positive)]
            if n > 2:
                families.append(("alternating_cycle", [x if p % 2 == 0 else -x for p,x in enumerate(positive)]))
            rng = random.Random(_sha(_json([seed,n,length])))
            alphabet = list(range(1-n,0)) + list(range(1,n))
            families.extend((f"random_signed_{i}", [rng.choice(alphabet) for _ in range(length)])
                            for i in range(random_per_cell))
            cells[n] = families
        for slot in range(max(map(len,cells.values()))):
            for n in ranks:
                if slot >= len(cells[n]):
                    continue
                family, word = cells[n][slot]
                case = {"n":n,"length":length,"family":family,"word":word,
                        "full_support":set(map(abs,word)) == set(range(1,n)),
                        "adjacent_inverse_pairs":sum(a == -b for a,b in zip(word,word[1:]))}
                case["case_id"] = _sha(_json(case))
                result.append(case)
    return result


def _execution_matches(certificate, expected):
    if expected.get("execution", "reference") == "collected":
        return (certificate.get("schema") == "zipper-collected-run-v1"
                and certificate.get("execution", {}).get("batch_size") == expected.get("batch_size",64))
    return certificate.get("schema") != "zipper-collected-run-v1"


def _cpu_time(execution):
    value = time.process_time()
    if execution == "collected":
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        value += children.ru_utime + children.ru_stime
    return value


def _process_entry(target, sender, case, limits, expected):
    # A separate session lets the supervising parent stop nested term pools.
    os.setsid()
    target(sender, case, limits, expected)


def _worker(sender, case, limits, expected):
    try:
        execution = expected.get("execution", "reference")
        wall, cpu = time.perf_counter(), _cpu_time(execution)
        if execution == "collected":
            from zipper_algebra.zipping_algorithm import collected
            session = collected.CollectedSession()
            reduce = lambda: session.reduce(case["word"], n=case["n"],
                workers=expected.get("term_workers",1), batch_size=expected.get("batch_size",64),
                limits=Limits(**limits), trace=False)
            check = session.check
        else:
            session = AlgorithmSession()
            reduce = lambda: session.reduce(case["word"],n=case["n"],limits=Limits(**limits),trace=False)
            check = session.check
        setup = {"setup_wall":time.perf_counter()-wall,"setup_cpu":_cpu_time(execution)-cpu,
                 "cpu_scope":"process_and_reaped_children" if execution == "collected" else "process"}
        sender.send({"kind":"phase","phase":"reduce",**setup})
        wall, cpu = time.perf_counter(), _cpu_time(execution)
        certificate = reduce()
        timing = {**setup,"reduce_wall":time.perf_counter()-wall,"reduce_cpu":_cpu_time(execution)-cpu,
                  "reduce_status":certificate["status"],"reduce_steps":certificate["metrics"]["steps"]}
        sender.send({"kind":"phase","phase":"check",**timing})
        encoded = _json(certificate)
        wall, cpu = time.perf_counter(), _cpu_time(execution)
        if check(json.loads(encoded)) is not True:
            raise ValueError("certificate failed JSON reload/replay")
        timing.update(check_wall=time.perf_counter()-wall,check_cpu=_cpu_time(execution)-cpu)
        if (certificate["source_sha256"] != expected["algorithm_source_sha256"]
                or certificate["catalog_sha256"] != expected["catalog_sha256"]):
            raise ValueError("source changed during computation")
        payload = encoded.encode()
        sender.send({"kind":"certificate","phase":"transfer","length":len(payload),
                     "certificate_sha":_sha(encoded),**timing})
        for start in range(0,len(payload),65536):
            sender.send_bytes(payload[start:start+65536])
    except BaseException as error:
        try:
            sender.send({"kind":"error","error":type(error).__name__,"message":str(error)})
        except (BrokenPipeError,OSError):
            pass
    finally:
        sender.close()


def run_case(case, limits, expected, wall_seconds, *, _target=_worker):
    """Return a checked result or a censored observation with no certificate."""
    ctx = multiprocessing.get_context("spawn")
    receiver, sender = ctx.Pipe(duplex=False)
    process = ctx.Process(target=_process_entry,args=(_target,sender,case,asdict(limits),expected))
    started = time.perf_counter()
    process.start()
    sender.close()
    detail = {"phase":"startup"}
    header, chunks, received = None, [], 0
    result = None
    try:
        while time.perf_counter()-started < wall_seconds:
            if not receiver.poll(min(.05,max(0,wall_seconds-(time.perf_counter()-started)))):
                if not process.is_alive():
                    break
                continue
            try:
                if header is not None:
                    chunk = receiver.recv_bytes(65536)
                    chunks.append(chunk); received += len(chunk)
                    if received > header["length"]:
                        raise ValueError("certificate transfer exceeds declared size")
                    if received == header["length"]:
                        encoded = b"".join(chunks).decode()
                        if _sha(encoded) != header["certificate_sha"]:
                            raise ValueError("certificate transfer hash mismatch")
                        certificate = json.loads(encoded)
                        if (certificate["n"] != case["n"] or certificate["input_word"] != case["word"]
                                or certificate["source_sha256"] != expected["algorithm_source_sha256"]
                                or certificate["catalog_sha256"] != expected["catalog_sha256"]
                                or certificate["limits"] != asdict(limits)
                                or certificate["trace_enabled"] is not False
                                or not _execution_matches(certificate,expected)):
                            raise ValueError("checked certificate does not match its case")
                        result = {"status":certificate["status"],"certificate":encoded,
                                  "certificate_sha":header["certificate_sha"],"metrics":certificate["metrics"],
                                  "detail":header}
                        break
                else:
                    message = receiver.recv()
                    if message["kind"] == "error":
                        result = {"status":"error","detail":message};break
                    detail = message
                    if message["kind"] == "certificate":
                        header = message
            except EOFError:
                break
        elapsed = time.perf_counter()-started
        if elapsed >= wall_seconds:
            result = None
        if result is None:
            result = {"status":"timeout" if elapsed >= wall_seconds else "worker_exit",
                      "detail":{**detail,"exitcode":process.exitcode}}
        result["attempt_wall"] = elapsed
        return result
    finally:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        if process.is_alive():
            process.terminate()
        process.join(2)
        if process.is_alive():
            process.kill();process.join()
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        receiver.close()


def summarize(db, campaign_id):
    groups = defaultdict(list)
    for n,length,family,status,wall,cpu,word,detail in db.execute(
            "SELECT n,length,family,status,reduce_wall,reduce_cpu,word,detail FROM runtime_cases WHERE campaign_id=?",
            (campaign_id,)):
        groups[n,length,family].append((status,wall,cpu,word,json.loads(detail)))
    result = []
    for (n,length,family), rows in sorted(groups.items()):
        counts = {status:sum(row[0] == status for row in rows) for status in sorted({row[0] for row in rows})}
        completed = [row for row in rows if row[0] == "complete"]
        attempted = len(rows)-counts.get("pending",0)
        timeout_phases = {}
        for row in rows:
            if row[0] == "timeout":
                phase = row[4].get("phase","unknown")
                timeout_phases[phase] = timeout_phases.get(phase,0)+1
        result.append({"n":n,"length":length,"family":family,"planned":len(rows),"attempted":attempted,
            "unique_words":len({row[3] for row in rows}),"statuses":counts,"timeout_phases":timeout_phases,
            "completion_fraction":len(completed)/attempted if attempted else None,
            "complete_reduce_wall_mean":statistics.mean(row[1] for row in completed) if completed else None,
            "complete_reduce_cpu_mean":statistics.mean(row[2] for row in completed) if completed else None,
            "complete_reduce_cpu_max":max((row[2] for row in completed),default=None)})
    return result


def run_sweep(database, *, ranks=(2,3,4,5,6,7,8), lengths=(4,10,12), seed=20261106,
              random_per_cell=1, limits=None, per_word_wall=30, wall_seconds=1800,
              campaign=None, max_new_cases=None, workers=None, execution="reference",
              term_workers=1, batch_size=64, _run_case=run_case):
    if execution not in ("reference", "collected"):
        raise ValueError("execution must be reference or collected")
    if type(term_workers) is not int or term_workers < 1 or type(batch_size) is not int or batch_size < 1:
        raise ValueError("positive integer term_workers and batch_size required")
    if execution == "reference" and term_workers != 1:
        raise ValueError("term_workers requires collected execution")
    if workers is None:
        workers = 1 if execution == "collected" else 4
    if type(workers) is not int or workers < 1:
        raise ValueError("positive worker count required")
    cases = sample_cases(ranks,lengths,seed=seed,random_per_cell=random_per_cell)
    limits = Limits(max_steps=30000,max_terms=100000,expansion_budget=1000000) if limits is None else limits
    if not isinstance(limits,Limits):
        raise ValueError("limits must be a Limits record")
    for value in (per_word_wall,wall_seconds):
        if type(value) not in (int,float) or not math.isfinite(value) or value <= 0:
            raise ValueError("positive finite wall budgets required")
    if max_new_cases is not None and (type(max_new_cases) is not int or max_new_cases < 0):
        raise ValueError("max_new_cases must be nonnegative")
    expected = identity() if execution == "reference" else identity(execution)
    expected.update(execution=execution, term_workers=term_workers, batch_size=batch_size)
    config = {"ranks":list(ranks),"lengths":list(lengths),"seed":seed,"random_per_cell":random_per_cell,
              "limits":asdict(limits),"per_word_wall":per_word_wall,"trace":False,
              "cases_sha256":_sha(_json(cases)),"concurrency":workers,
              "execution":execution,"term_workers":term_workers,"batch_size":batch_size}
    campaign = campaign or "runtime-" + _sha(_json([expected,config]))[:20]
    if not isinstance(campaign,str) or not campaign:
        raise ValueError("nonempty campaign identifier required")
    database = Path(database);database.parent.mkdir(parents=True,exist_ok=True)
    with (database.parent/'.parallel.lock').open('a+') as lock:
        try:fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as error:raise ValueError("another writer owns this database") from error
        db = sqlite3.connect(database)
        try:
            db.execute('PRAGMA journal_mode=WAL')
            db.execute('CREATE TABLE IF NOT EXISTS runtime_campaigns (campaign_id TEXT PRIMARY KEY, identity TEXT, config TEXT, created_utc TEXT, git_commit TEXT)')
            db.execute('CREATE TABLE IF NOT EXISTS runtime_cases (campaign_id TEXT, case_id TEXT, ordinal INTEGER, n INTEGER, length INTEGER, family TEXT, word TEXT, input_sha TEXT, full_support INTEGER, status TEXT, attempt_wall REAL, reduce_wall REAL, reduce_cpu REAL, check_wall REAL, check_cpu REAL, metrics TEXT, certificate TEXT, certificate_sha TEXT, detail TEXT, PRIMARY KEY(campaign_id,case_id), UNIQUE(campaign_id,ordinal))')
            old = db.execute('SELECT identity,config FROM runtime_campaigns WHERE campaign_id=?',(campaign,)).fetchone()
            if old and old != (_json(expected),_json(config)):
                raise ValueError('runtime campaign source or configuration changed; use a new campaign')
            if not old:
                commit = subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True).stdout.strip()
                with db:
                    db.execute('INSERT INTO runtime_campaigns VALUES (?,?,?,?,?)',(campaign,_json(expected),_json(config),datetime.now(timezone.utc).isoformat(),commit))
                    for ordinal,case in enumerate(cases):
                        db.execute('INSERT INTO runtime_cases VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',
                            (campaign,case['case_id'],ordinal,case['n'],case['length'],case['family'],_json(case['word']),
                             _sha(_json({'n':case['n'],'word':case['word']})),int(case['full_support']),'pending',
                             None,None,None,None,None,None,None,None,_json({'adjacent_inverse_pairs':case['adjacent_inverse_pairs']})))
            rows = db.execute('SELECT case_id,status,certificate,certificate_sha,n,length,family,word,input_sha,full_support,ordinal FROM runtime_cases WHERE campaign_id=? ORDER BY ordinal',(campaign,)).fetchall()
            if [row[0] for row in rows] != [case['case_id'] for case in cases]:
                raise ValueError('stored runtime cases do not match the deterministic sample')
            for ordinal,(case,row) in enumerate(zip(cases,rows)):
                _,status,certificate,digest,*stored = row
                required = [case["n"],case["length"],case["family"],_json(case["word"]),
                            _sha(_json({"n":case["n"],"word":case["word"]})),int(case["full_support"]),ordinal]
                if stored != required:
                    raise ValueError("stored runtime case columns differ from the deterministic sample")
                if status not in ("pending","complete","incomplete","control_obligation","timeout","error","worker_exit"):
                    raise ValueError("unknown stored runtime status")
                if certificate is None and status in ("complete","incomplete","control_obligation"):
                    raise ValueError("checked runtime status requires a certificate")
                if certificate is not None:
                    value = json.loads(certificate)
                    if (_sha(certificate) != digest or value['source_sha256'] != expected['algorithm_source_sha256']
                            or value['catalog_sha256'] != expected['catalog_sha256'] or value['status'] != status
                            or value['n'] != case['n'] or value['input_word'] != case['word']
                            or value['limits'] != asdict(limits) or value['trace_enabled'] is not False
                            or not _execution_matches(value,expected)):
                        raise ValueError('stored runtime certificate hash or identity mismatch')
            started = time.perf_counter();new = 0;submitted = 0;reason = 'finished'
            remaining = [case for case,row in zip(cases,rows) if row[1] == 'pending']
            pending = {}
            with ThreadPoolExecutor(max_workers=workers) as pool:
                while remaining or pending:
                    while remaining and len(pending) < workers:
                        if time.perf_counter()-started >= wall_seconds:
                            reason='campaign_wall_budget';break
                        if max_new_cases is not None and submitted >= max_new_cases:
                            reason='case_budget';break
                        case = remaining.pop(0)
                        pending[pool.submit(_run_case,case,limits,expected,per_word_wall)] = case
                        submitted += 1
                    if not pending:break
                    finished,_ = wait(pending,return_when=FIRST_COMPLETED)
                    for future in finished:
                        case = pending.pop(future)
                        result = future.result()
                        if not _source_matches(expected):raise ValueError('source changed during runtime sweep')
                        detail=result.get('detail',{})
                        with db:
                            db.execute('UPDATE runtime_cases SET status=?,attempt_wall=?,reduce_wall=?,reduce_cpu=?,check_wall=?,check_cpu=?,metrics=?,certificate=?,certificate_sha=?,detail=? WHERE campaign_id=? AND case_id=?',
                                (result['status'],result['attempt_wall'],detail.get('reduce_wall'),detail.get('reduce_cpu'),
                                 detail.get('check_wall'),detail.get('check_cpu'),_json(result.get('metrics')),
                                 result.get('certificate'),result.get('certificate_sha'),_json(detail),campaign,case['case_id']))
                        new += 1
            return {'schema':'zipper-selected-runtime-sweep-v1','campaign_id':campaign,'database':str(database.resolve()),
                    'identity':expected,'config':config,'new_attempts':new,'stop_reason':reason,
                    'elapsed_wall':time.perf_counter()-started,'strata':summarize(db,campaign)}
        finally:db.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('database');parser.add_argument('--campaign')
    parser.add_argument('--ranks',nargs='+',type=int,default=[2,3,4,5,6,7,8])
    parser.add_argument('--lengths',nargs='+',type=int,default=[4,10,12])
    parser.add_argument('--seed',type=int,default=20261106)
    parser.add_argument('--random-per-cell',type=int,default=1)
    parser.add_argument('--limits',type=json.loads)
    parser.add_argument('--per-word-wall',type=float,default=30)
    parser.add_argument('--wall-seconds',type=float,default=1800)
    parser.add_argument('--max-new-cases',type=int)
    parser.add_argument('--workers',type=int)
    parser.add_argument('--execution',choices=('reference','collected'),default='reference')
    parser.add_argument('--term-workers',type=int,default=1)
    parser.add_argument('--batch-size',type=int,default=64)
    args=vars(parser.parse_args());args['limits']=Limits(**args['limits']) if args['limits'] else None
    print(json.dumps(run_sweep(**args),indent=2),flush=True)


if __name__=='__main__':main()
