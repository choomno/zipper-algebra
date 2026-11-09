from copy import deepcopy
from dataclasses import asdict, replace
import json
import sqlite3
import time

import pytest

from zipper_algebra.zipping_algorithm.census_queue import CensusQueue, Claim, LostLease, SCHEMA
from zipper_algebra.zipping_algorithm.census_reuse import reconstruction_sources
from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native_verifier import VerifierIncomplete, VerifierError
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import execution_identity
from experiments.zipper_async_census import admit, read_policy

IDENTITY = {"abi": "test", "source_sha256": {"a": "b"}, "binary_sha256": "test"}


@pytest.fixture
def setup(tmp_path):
    session = CollectedSession()
    limits = Limits()
    seed = session.reduce([], n=1, limits=limits)
    policy = {"schema": SCHEMA, "producer": execution_identity(seed), "producer_metadata": {"test": True},
              "limits": asdict(limits), "package_sources": reconstruction_sources()}
    queue = CensusQueue(tmp_path / "queue.sqlite3", policy)
    yield queue, lambda n, word: session.reduce(word, n=n, limits=limits), session.check
    queue.close()


def test_work_presence_ignores_certified_history(setup):
    queue, produce, verify = setup
    assert not any(queue.work_presence().values())
    queue.enqueue(3, [-1])
    assert queue.work_presence()["queued"]
    queue.compute_one(produce)
    assert queue.work_presence()["candidate"]
    queue.certify_one(verify, IDENTITY)
    statements = []
    queue.db.set_trace_callback(statements.append)
    assert not any(queue.work_presence().values())
    queue.db.set_trace_callback(None)
    assert len(statements) == 1
    assert "COUNT(" not in statements[0].upper()
    assert "ASYNC_CANDIDATES" not in statements[0].upper()
    plans = queue.db.execute("EXPLAIN QUERY PLAN " + statements[0]).fetchall()
    assert all("SCAN async_jobs" not in row[3] for row in plans)


def test_durable_separate_compute_and_certify(setup):
    queue, produce, verify = setup
    key = queue.enqueue(3, [-1, -2])
    assert queue.enqueue(3, [-1, -2]) == key
    assert queue.compute_one(produce)
    assert queue.counts()["candidate"] == 1
    assert "certified" not in queue.counts()
    other = CensusQueue(queue.path, queue.policy)
    try:
        assert other.certify_one(verify, IDENTITY)
        assert other.counts()["certified"] == 1
        assert not other.certify_one(verify, IDENTITY)
        records = other.db.execute("SELECT stage,status,record FROM async_attempts").fetchall()
        assert [r[1] for r in records] == ["candidate", "certified"]
        assert json.loads(records[0][2])["digest"] == json.loads(records[1][2])["digest"]
    finally:
        other.close()


def test_claim_recovery_prevents_stale_commit(setup):
    queue, produce, _ = setup
    queue.enqueue(3, [-1])
    old = queue.claim("compute")
    assert old
    queue.db.execute("UPDATE async_jobs SET lease=0")
    new = queue.claim("compute")
    assert new and new.key == old.key and new.token != old.token
    assert not queue.heartbeat(old, 300)
    with pytest.raises(LostLease):
        queue.publish(old, produce(3, [-1]), seconds=1)
    queue.publish(new, produce(3, [-1]), seconds=1)
    check = queue.claim("certify")
    assert check
    queue.db.execute("UPDATE async_jobs SET lease=0")
    recovered = queue.claim("certify")
    assert recovered
    with pytest.raises(LostLease):
        queue.finish(check, "certified", {})
    assert queue.heartbeat(recovered, 30)
    assert queue.db.execute("SELECT count(*) FROM async_attempts WHERE status='lease_expired'").fetchone()[0] == 2


def test_heartbeat_and_finalize_race(setup):
    queue, produce, _ = setup
    queue.enqueue(3, [-1])
    claim = queue.claim("compute", lease_seconds=.1)
    assert claim
    with queue.keepalive(claim, .1):
        time.sleep(.15)
        assert queue.claim("compute") is None
        queue.publish(claim, produce(3, [-1]), seconds=1)
        time.sleep(.1)  # heartbeat sees token cleared by successful completion
    assert queue.counts()["candidate"] == 1


def test_backpressure_and_fifo(setup):
    queue, produce, verify = setup
    first = queue.enqueue(3, [-1])
    second = queue.enqueue(3, [-2])
    claim = queue.claim("compute", max_pending=1)
    assert claim and claim.key == first
    assert queue.claim("compute", max_pending=1) is None
    queue.publish(claim, produce(3, [-1]), seconds=0)
    assert queue.claim("compute", max_pending_bytes=1) is None
    assert queue.certify_one(verify, IDENTITY)
    claim = queue.claim("compute", max_pending=1)
    assert claim and claim.key == second


@pytest.mark.parametrize("failure,status,reason", [
    (False, "rejected", "replay_mismatch"),
    (VerifierIncomplete("budget"), "verification_incomplete", "resource_limit"),
    (VerifierIncomplete("control", replay_status="control_obligation",
                        replay_diagnostic={"rule": "pending"}),
     "verification_incomplete", "control_obligation"),
    (VerifierIncomplete("unsupported", replay_status="unsupported"),
     "verification_incomplete", "unsupported"),
    (VerifierError("source changed"), "verification_error", "operational_error"),
    (RuntimeError("process died"), "verification_error", "operational_error"),
])
def test_verification_outcomes_preserve_candidate(setup, failure, status, reason):
    queue, produce, _ = setup
    queue.enqueue(3, [-1])
    queue.compute_one(produce)
    blob = queue.db.execute("SELECT payload FROM async_candidates").fetchone()[0]
    def check(_):
        if isinstance(failure, Exception):
            raise failure
        return failure
    assert queue.certify_one(check, IDENTITY)
    assert queue.counts()[status] == 1
    assert queue.db.execute("SELECT payload FROM async_candidates").fetchone()[0] == blob
    record = json.loads(queue.db.execute("SELECT record FROM async_attempts ORDER BY id DESC").fetchone()[0])
    assert record["reason"] == reason
    if isinstance(failure, VerifierIncomplete):
        assert record["replay_status"] == failure.replay_status
        assert record["replay_diagnostic"] == failure.replay_diagnostic
    if status == "rejected":
        with pytest.raises(ValueError):
            queue.retry(status)
    else:
        assert queue.retry(status) == 1
        assert queue.counts()["candidate"] == 1


@pytest.mark.parametrize("column,value", [("digest", "bad"), ("output_digest", "bad"), ("payload", b"broken")])
def test_corruption_is_not_certified(setup, column, value):
    queue, produce, verify = setup
    queue.enqueue(3, [-1])
    queue.compute_one(produce)
    queue.db.execute(f"UPDATE async_candidates SET {column}=?", (value,))
    assert queue.certify_one(verify, IDENTITY)
    assert queue.counts()["rejected"] == 1
    record = json.loads(queue.db.execute("SELECT record FROM async_attempts ORDER BY id DESC").fetchone()[0])
    assert record["reason"] == "candidate_integrity"


def test_atomic_publish_and_binding(setup):
    queue, produce, _ = setup
    queue.enqueue(3, [-1])
    claim = queue.claim("compute")
    assert claim
    value = produce(3, [-1])
    wrong = deepcopy(value)
    wrong["input_word"] = [-2]
    with pytest.raises(ValueError, match="binding"):
        queue.publish(claim, wrong, seconds=0)
    queue.db.execute("CREATE TRIGGER fail_attempt BEFORE INSERT ON async_attempts BEGIN SELECT RAISE(ABORT,'test'); END")
    with pytest.raises(sqlite3.IntegrityError):
        queue.publish(claim, value, seconds=0)
    assert queue.db.execute("SELECT count(*) FROM async_candidates").fetchone()[0] == 0
    assert queue.counts()["computing"] == 1
    queue.db.execute("DROP TRIGGER fail_attempt")
    queue.publish(claim, value, seconds=0)


def test_incomplete_compute_and_retry(setup):
    queue, _, _ = setup
    queue.enqueue(3, [-1])
    assert queue.compute_one(lambda *_: {"status": "incomplete", "diagnostic": "limit"})
    assert queue.counts()["compute_incomplete"] == 1
    assert queue.retry("compute_incomplete") == 1
    def fail(*_):
        raise RuntimeError("failure")
    with pytest.raises(RuntimeError):
        queue.compute_one(fail)
    assert queue.counts()["compute_error"] == 1


def test_policy_isolation_sources_and_legacy(setup, tmp_path):
    queue, _, _ = setup
    queue.enqueue(3, [-1])
    policy = deepcopy(queue.policy)
    policy["limits"]["max_steps"] += 1
    other = CensusQueue(queue.path, policy)
    try:
        assert other.counts() == {"stored_candidate_bytes": 0}
        assert read_policy(queue.path, queue.policy_key) == queue.policy
        with pytest.raises(ValueError):
            read_policy(queue.path, None)
    finally:
        other.close()
    queue.policy["package_sources"] = {"wrong": "source"}
    with pytest.raises(ValueError, match="source"):
        queue.compute_one(lambda *_: {})
    old = tmp_path / "legacy.sqlite3"
    with sqlite3.connect(old) as db:
        db.execute("CREATE TABLE old_records(a TEXT)")
    before = old.read_bytes()
    with pytest.raises(ValueError, match="sidecar"):
        CensusQueue(old, policy)
    assert old.read_bytes() == before


def test_cursor_atomicity_and_primitive_admission(setup):
    queue, _, _ = setup
    def bad(i):
        if i == 1:
            raise RuntimeError("boom")
        return 3, [-1]
    with pytest.raises(RuntimeError):
        queue.enqueue_interval("test", 2, bad)
    assert queue.counts() == {"stored_candidate_bytes": 0}
    assert admit(queue, {"3": 3, "4": 3}, maximum=17) == 17
    assert admit(queue, {"3": 3, "4": 3}, maximum=10000) > 0
    assert admit(queue, {"3": 3, "4": 3}, maximum=10000) == 0
    assert queue.db.execute("SELECT min(n) FROM async_jobs").fetchone()[0] == 3
    with pytest.raises(ValueError, match="rank"):
        admit(queue, {"2": 10}, maximum=10)


@pytest.mark.parametrize("n,word", [(0, []), (3, [0]), (3, [3]), (3, [True]), (True, [])])
def test_bad_requests(setup, n, word):
    with pytest.raises(ValueError):
        setup[0].enqueue(n, word)


def test_bad_claims(setup):
    queue, produce, _ = setup
    for options in ({"stage": "bad"}, {"stage": "compute", "lease_seconds": 0},
                    {"stage": "compute", "max_pending": 0}):
        with pytest.raises(ValueError):
            queue.claim(**options)
    queue.enqueue(3, [-1])
    claim = queue.claim("compute")
    assert claim
    with pytest.raises(ValueError):
        queue.finish(claim, "certified", {})
    with pytest.raises(LostLease):
        queue.publish(replace(claim, word=[-2]), produce(3, [-2]), seconds=0)
    with pytest.raises(ValueError):
        queue.publish(Claim(claim.key, 3, [-1], claim.token, "certify"), produce(3, [-1]), seconds=0)
    with pytest.raises(ValueError):
        queue.certify_one(lambda _: True, {})


def test_native_build_binding_and_exact_output(tmp_path):
    from zipper_algebra.zipping_algorithm.native import NativeSession
    with NativeSession() as producer:
        limits = Limits(max_steps=500)
        result = producer.reduce_reused([-1, -2], n=3, limits=limits)
        policy = {"schema": SCHEMA, "producer": execution_identity(result),
                  "producer_metadata": producer.metadata, "limits": asdict(limits),
                  "package_sources": reconstruction_sources()}
        queue = CensusQueue(tmp_path / "native.sqlite3", policy)
        try:
            queue.enqueue(3, [-1, -2])
            claim = queue.claim("compute")
            assert claim
            wrong = deepcopy(result)
            wrong["native_build"]["binary_sha256"] = "bad"
            with pytest.raises(ValueError, match="native build"):
                queue.publish(claim, wrong, seconds=0)
            queue.publish(claim, result, seconds=0)
            check = queue.claim("certify")
            assert check
            loaded, _ = queue.candidate(check)
            assert loaded == result
        finally:
            queue.close()


def test_missing_payload_lost_leases_and_no_work(setup, monkeypatch):
    queue, produce, verify = setup
    assert not queue.compute_one(produce)
    queue.enqueue(3, [-1])
    def lost(*_):
        raise LostLease("expired")
    with pytest.raises(LostLease):
        queue.compute_one(lost)
    queue.db.execute("UPDATE async_jobs SET lease=0")
    queue.compute_one(produce)
    with monkeypatch.context() as patch:
        patch.setattr(queue, "candidate", lost)
        with pytest.raises(LostLease):
            queue.certify_one(verify, IDENTITY)
    queue.db.execute("UPDATE async_jobs SET lease=0")
    claim = queue.claim("certify")
    assert claim
    queue.db.execute("DELETE FROM async_candidates")
    with pytest.raises(ValueError, match="missing"):
        queue.candidate(claim)


def test_heartbeat_connection_failure_is_recoverable(setup, monkeypatch):
    import threading
    import zipper_algebra.zipping_algorithm.census_queue as module
    queue, produce, _ = setup
    queue.enqueue(3, [-1])
    claim = queue.claim("compute")
    assert claim
    called = threading.Event()
    def fail(*args, **kwargs):
        called.set()
        raise sqlite3.OperationalError("offline")
    with monkeypatch.context() as patch:
        patch.setattr(module.sqlite3, "connect", fail)
        with queue.keepalive(claim, .1):
            assert called.wait(1)
            queue.publish(claim, produce(3, [-1]), seconds=0)
    assert queue.counts()["candidate"] == 1


def test_invalid_policy_and_cursor_bounds(setup, tmp_path):
    queue, _, _ = setup
    with pytest.raises(ValueError, match="explicit"):
        CensusQueue(tmp_path / "bad.sqlite3", {})
    wrong = deepcopy(queue.policy)
    wrong["package_sources"] = {"different": "source"}
    with pytest.raises(ValueError, match="imported"):
        CensusQueue(tmp_path / "bad.sqlite3", wrong)
    queue.db.execute("UPDATE async_policies SET payload='{}'")
    with pytest.raises(ValueError, match="policy binding"):
        CensusQueue(queue.path, queue.policy)
    for count, batch in ((-1, 3), (2, 0)):
        with pytest.raises(ValueError):
            queue.enqueue_interval("bad", count, lambda i: (3, [1]), batch=batch)


def test_retry_is_logged_and_admission_rejects_source_drift(setup, monkeypatch):
    import zipper_algebra.zipping_algorithm.census_queue as module
    queue, produce, _ = setup
    queue.enqueue(3, [-1])
    queue.compute_one(lambda *_: {"status": "incomplete"})
    assert queue.retry("compute_incomplete") == 1
    stage, status, record = queue.db.execute("SELECT stage,status,record FROM async_attempts ORDER BY id DESC").fetchone()
    assert (stage, status) == ("retry", "queued")
    assert json.loads(record) == {"previous_state": "compute_incomplete"}
    monkeypatch.setattr(module, "_IMPORT_SOURCES", {"altered": "snapshot"})
    with pytest.raises(ValueError, match="source"):
        queue.enqueue(3, [-2])
    with pytest.raises(ValueError, match="source"):
        queue.enqueue_interval("bad", 2, lambda _: (3, [-2]))
    assert queue.db.execute("SELECT count(*) FROM async_cursors").fetchone()[0] == 0


def test_source_change_during_load_is_operational(setup, monkeypatch):
    from zipper_algebra.zipping_algorithm.census_queue import SourceChanged
    queue, produce, verify = setup
    queue.enqueue(3, [-1])
    queue.compute_one(produce)
    def drift(_):
        raise SourceChanged("changed during loading")
    monkeypatch.setattr(queue, "candidate", drift)
    with pytest.raises(SourceChanged):
        queue.certify_one(verify, IDENTITY)
    assert queue.counts()["certifying"] == 1
    assert "rejected" not in queue.counts()


def test_short_admissions_rotate_durably(setup):
    queue, _, _ = setup
    assert admit(queue, {"3": 4, "4": 3}, maximum=1) == 1
    assert admit(queue, {"3": 4, "4": 3}, maximum=1) == 1
    assert admit(queue, {"3": 4, "4": 3}, maximum=1) == 1
    assert admit(queue, {"3": 4, "4": 3}, maximum=1) == 1
    rows = queue.db.execute("SELECT n,word FROM async_jobs ORDER BY rowid").fetchall()
    assert [(r[0], len(json.loads(r[1]))) for r in rows] == [(3, 2), (3, 3), (3, 4), (4, 3)]


def test_readonly_historical_status(setup, capsys):
    from zipper_algebra.zipping_algorithm.census_reuse import _sha, _json
    from experiments.zipper_async_census import main
    queue, _, _ = setup
    old = deepcopy(queue.policy)
    old["package_sources"] = {"old": "source"}
    key = _sha(old)
    queue.db.execute("INSERT INTO async_policies VALUES (?,?)", (key, _json(old)))
    before = queue.db.total_changes
    assert main(["status", str(queue.path), "--policy", key]) == 0
    answer = json.loads(capsys.readouterr().out)
    assert answer["counts"] == {"stored_candidate_bytes": 0}
    assert queue.db.total_changes == before
    old["schema"] = "foreign"
    queue.db.execute("UPDATE async_policies SET payload=? WHERE key=?", (_json(old), key))
    with pytest.raises(ValueError, match="schema or digest"):
        main(["status", str(queue.path), "--policy", key])
    missing = queue.path.with_name("does-not-exist.sqlite3")
    with pytest.raises(sqlite3.OperationalError):
        main(["status", str(missing)])
    assert not missing.exists()


def _worker_policy(queue):
    import hashlib
    from pathlib import Path
    import experiments.zipper_async_census as cli
    policy = deepcopy(queue.policy)
    policy.update(mode="standalone", runner_sha256=hashlib.sha256(Path(cli.__file__).read_bytes()).hexdigest())
    return policy


def _worker_settings():
    return dict(cache_bytes=0, allocator=None, wall_seconds=3, once=False, drain=True,
                max_pending=3, max_pending_bytes=1000000)


def test_dead_verifier_stops_after_first_candidate(setup, monkeypatch):
    import experiments.zipper_async_census as cli
    queue, produce, _ = setup
    policy = _worker_policy(queue)
    work = CensusQueue(queue.path, policy)
    try:
        for word in ([-1], [-2]):
            work.enqueue(3, word)
            work.compute_one(produce)
        class DeadVerifier:
            metadata = IDENTITY
            last_result = {}
            pid = 1
            def __init__(self, **kwargs):
                pass
            def check(self, result):
                self.pid = None
                raise VerifierError("Verifier exited without a response")
            def close(self):
                pass
        monkeypatch.setattr(cli, "CppVerifier", DeadVerifier)
        with pytest.raises(RuntimeError, match="remaining candidates preserved"):
            cli.worker(queue.path, policy, "certify", _worker_settings())
        assert work.counts()["verification_error"] == 1
        assert work.counts()["candidate"] == 1
    finally:
        work.close()


def test_worker_lost_lease_retries_normally(setup, monkeypatch):
    import experiments.zipper_async_census as cli
    queue, _, _ = setup
    policy = _worker_policy(queue)
    calls = []
    class Verifier:
        metadata = IDENTITY
        pid = 1
        def __init__(self, **kwargs):
            pass
        def check(self, result):
            return True
        def close(self):
            pass
    def certify(self, *args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise LostLease("reclaimed")
        return False
    monkeypatch.setattr(cli, "CppVerifier", Verifier)
    monkeypatch.setattr(cli.CensusQueue, "certify_one", certify)
    cli.worker(queue.path, policy, "certify", _worker_settings())
    assert len(calls) == 2


def test_interval_without_rotation(setup):
    queue, _, _ = setup
    assert queue.enqueue_interval("plain", 1, lambda _: (3, [-1])) == 1
    assert queue.enqueue_interval("plain", 1, lambda _: (3, [-1])) == 0
    assert queue.db.execute("SELECT ordinal FROM async_cursors WHERE layer='plain'").fetchone()[0] == "1"


@pytest.mark.parametrize("once,dead,expected_calls", [(True, False, 1), (False, False, 2), (True, True, 1)])
def test_lost_lease_respects_once_and_dead_verifier(setup, monkeypatch, once, dead, expected_calls):
    import experiments.zipper_async_census as cli
    queue, _, _ = setup
    calls = []
    class Verifier:
        metadata = IDENTITY
        pid = None if dead else 1
        def __init__(self, **kwargs):
            pass
        def close(self):
            pass
    def certify(self, *args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise LostLease("reclaimed")
        return False
    monkeypatch.setattr(cli, "CppVerifier", Verifier)
    monkeypatch.setattr(cli.CensusQueue, "certify_one", certify)
    settings = dict(_worker_settings(), once=once)
    if dead:
        with pytest.raises(RuntimeError, match="remaining candidates preserved"):
            cli.worker(queue.path, _worker_policy(queue), "certify", settings)
    else:
        cli.worker(queue.path, _worker_policy(queue), "certify", settings)
    assert len(calls) == expected_calls


def test_minimum_length_admission_and_cli(setup, capsys):
    import experiments.zipper_async_census as cli
    from zipper_algebra.zipping_algorithm.peeling_enumeration import NonPeelableWords
    queue, _, _ = setup
    assert NonPeelableWords(3, 9).count == 252497
    policy = _worker_policy(queue)
    work = CensusQueue(queue.path, policy)
    try:
        assert cli.main(["enqueue", str(queue.path), "--policy", work.policy_key,
                         "--targets", '{"3":9}', "--min-length", "9", "--admit-max", "3"]) == 0
        assert json.loads(capsys.readouterr().out)["counts"]["queued"] == 3
        rows = work.db.execute("SELECT n,word FROM async_jobs WHERE policy=?", (work.policy_key,)).fetchall()
        assert all(n == 3 and len(json.loads(word)) == 9 for n, word in rows)
        assert admit(work, {"3": 8}, minimum_length=9, maximum=2) == 0
        for value in (-1, True, 1.5):
            with pytest.raises(ValueError, match="minimum length"):
                admit(work, {"3": 9}, minimum_length=value, maximum=1)
        with pytest.raises(SystemExit) as error:
            cli.main(["enqueue", str(queue.path), "--min-length", "-1"])
        assert error.value.code == 2
    finally:
        work.close()
