"""Captured benchmarks preserve exact census rows without writing the census."""
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import sqlite3
import zlib

import pytest

from experiments.zipper_cpp_verifier_benchmark import capture, encoded, load_certificate


def test_capture_preserves_census_and_freezes_case_selection(tmp_path):
    source, target = tmp_path / "census.sqlite3", tmp_path / "benchmark.sqlite3"
    version = encoded({"engine": "test"})
    version_hash = hashlib.sha256(version.encode()).hexdigest()
    with sqlite3.connect(source) as db:
        db.executescript('''
            CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);
            CREATE TABLE versions(hash TEXT PRIMARY KEY,identity TEXT);
            CREATE TABLE words(n,length,ordinal,word,input_sha,status,metrics,
              certificate,certificate_sha,certificate_bytes,config_sha,checker_ns);
        ''')
        db.executemany("INSERT INTO metadata VALUES(?,?)", [
            ("last_config", '{"execution":"native"}'), ("run_state", '{"config_sha":"test"}'),
            ("active_version", version_hash)])
        db.execute("INSERT INTO versions VALUES(?,?)", (version_hash, version))
        for n in range(3, 9):
            for ordinal in range(10):
                raw = encoded({"schema": "zipper-native-run-v1", "n": n, "input_word": [1]}).encode()
                db.execute("INSERT INTO words VALUES(?,?,?,?,?,?,?,?,?,?,?,?)", (
                    n, 1, ordinal, "[1]", "input-hash", "complete", encoded({"steps": ordinal}),
                    zlib.compress(raw), hashlib.sha256(raw).hexdigest(), len(raw), "test", ordinal))
    before = source.read_bytes()
    assert capture(source, target, recent=10) == 18
    assert source.read_bytes() == before
    with sqlite3.connect(target) as db:
        assert db.execute("SELECT ordinal FROM cases WHERE n=3 ORDER BY id").fetchall() == [(4,), (8,), (9,)]
        row = db.execute("SELECT certificate,certificate_sha,certificate_bytes FROM cases LIMIT 1").fetchone()
        assert load_certificate(*row)["input_word"] == [1]
    assert capture(source, target, recent=3) == 18
    assert source.read_bytes() == before
    other = tmp_path / "other-census.sqlite3"
    shutil.copyfile(source, other)
    with pytest.raises(ValueError, match="different census source"):
        capture(other, target)
    # Later progress/epochs may retain the frozen historical source version.
    with sqlite3.connect(source) as db:
        db.execute("UPDATE metadata SET value='new-version' WHERE key='active_version'")
    assert capture(source, target) == 18
    with sqlite3.connect(source) as db:
        db.execute("DELETE FROM versions")
    with pytest.raises(ValueError, match="version is absent"):
        capture(source, target)
    with sqlite3.connect(source) as db:
        db.execute("INSERT INTO metadata VALUES('identity','changed')")
    with pytest.raises(ValueError, match="identity differs"):
        capture(source, target)


def test_capture_refuses_census_as_destination(tmp_path):
    source = tmp_path / "census.sqlite3"
    with pytest.raises(ValueError, match="must differ"):
        capture(source, source)
    assert not source.exists()


def test_certificate_integrity_and_large_integer_decoding():
    raw = b'{"coefficient":[[0,' + b'1' * 5000 + b']]}'
    assert load_certificate(zlib.compress(raw), hashlib.sha256(raw).hexdigest(), len(raw))["coefficient"][0][1] > 10**4999
    with pytest.raises(ValueError, match="SHA256"):
        load_certificate(zlib.compress(raw), "incorrect", len(raw))
    with pytest.raises(ValueError, match="size"):
        load_certificate(zlib.compress(raw), hashlib.sha256(raw).hexdigest(), len(raw) + 1)
    with pytest.raises(ValueError, match="object"):
        load_certificate(zlib.compress(b"[]"), hashlib.sha256(b"[]").hexdigest(), 2)


def test_worker_releases_case_read_before_another_writer(tmp_path, monkeypatch):
    from experiments.zipper_cpp_verifier_benchmark import worker
    from zipper_algebra.zipping_algorithm import independent_replay, native

    database = tmp_path / "benchmark.sqlite3"
    raw = b'{"n":3,"input_word":[1]}'
    with sqlite3.connect(database) as db:
        db.executescript('''
          CREATE TABLE cases(id INTEGER PRIMARY KEY,n,length,ordinal,certificate,certificate_sha,certificate_bytes);
          CREATE TABLE runs(id INTEGER PRIMARY KEY,backend,started_utc,status,cache_bytes,cache_entries,
                            benchmark_sha,setup_seconds,identity,error);
          CREATE TABLE measurements(run,pass,case_id,seconds,cache);
          CREATE TABLE external_writes(value INTEGER);
        ''')
        for index in (1, 2):
            db.execute("INSERT INTO cases VALUES(?,?,?,?,?,?,?)", (
                index, 3, 1, index, zlib.compress(raw), hashlib.sha256(raw).hexdigest(), len(raw)))

    class Producer:
        metadata = {}

        def __init__(self, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    class Checker:
        def __init__(self, **kwargs):
            pass

        def check_native(self, producer, certificate):
            # A real second connection commits between the case read and the
            # measurement write. A retained SELECT cursor would leave the first
            # connection on an old WAL snapshot and make its write fail.
            with sqlite3.connect(database, timeout=0.1) as other:
                other.execute("INSERT INTO external_writes VALUES(1)")
            return True

        def cache_info(self):
            return {}

    monkeypatch.setattr(native, "NativeSession", Producer)
    monkeypatch.setattr(independent_replay, "IndependentVerifier", Checker)
    monkeypatch.setattr(independent_replay, "source_identity", lambda: {})
    assert worker(database, "python", 1, 0, 0) == 0
    with sqlite3.connect(database) as db:
        assert db.execute("SELECT COUNT(*) FROM external_writes").fetchone()[0] == 2
        assert db.execute("SELECT COUNT(*) FROM measurements").fetchone()[0] == 2
        assert db.execute("SELECT status FROM runs").fetchone()[0] == "complete"


@pytest.mark.parametrize("interrupted,initial_status,double_interrupt", [
    (True, "running", False), (False, "running", False),
    (True, "complete", False), (True, "failed", False), (True, "running", True),
])
def test_exit_terminates_stubborn_owned_worker_group(tmp_path, monkeypatch, interrupted,
                                                     initial_status, double_interrupt):
    from experiments import zipper_cpp_verifier_benchmark as benchmark

    database = tmp_path / "benchmark.sqlite3"
    with sqlite3.connect(database) as db:
        db.execute("CREATE TABLE runs(id INTEGER PRIMARY KEY,backend,status,error,peak_tree_rss_bytes,peak_process_rss)")
    child_pid_file = tmp_path / "child.pid"
    ready_file = tmp_path / "child.ready"
    child_code = ("import signal,time,pathlib; signal.signal(signal.SIGTERM,signal.SIG_IGN); "
                  f"pathlib.Path({str(ready_file)!r}).write_text('ready'); time.sleep(60)")
    code = (
        "import subprocess,sys,time,pathlib,sqlite3,signal; "
        + ("signal.signal(signal.SIGTERM,signal.SIG_IGN); " if double_interrupt else "")
        +
        f"db=sqlite3.connect({str(database)!r}); "
        f"db.execute(\"INSERT INTO runs(backend,status,error) VALUES(?,?,?)\",('cpp',{initial_status!r},{'saved failure' if initial_status == 'failed' else None!r})); db.commit(); db.close(); "
        f"child=subprocess.Popen([sys.executable,'-c',{child_code!r}]); "
        f"ready=pathlib.Path({str(ready_file)!r}); "
        "exec('while not ready.exists(): time.sleep(0.01)'); "
        f"pathlib.Path({str(child_pid_file)!r}).write_text(str(child.pid)); "
        + ("time.sleep(60)" if interrupted else "sys.exit(3)")
    )
    real_popen = subprocess.Popen
    processes = []

    def launch(command, **kwargs):
        assert kwargs["start_new_session"] is True
        process = real_popen([sys.executable, "-c", code], **kwargs)
        processes.append(process)
        deadline = time.monotonic() + 5
        while not child_pid_file.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        assert child_pid_file.exists()
        assert os.getpgid(int(child_pid_file.read_text())) == process.pid
        if double_interrupt:
            original_wait = process.wait
            def wait(timeout=None):
                if timeout == 5:
                    raise KeyboardInterrupt("second interrupt")
                return original_wait(timeout=timeout)
            process.wait = wait
        return process

    def interrupt(pid):
        raise KeyboardInterrupt

    monkeypatch.setattr(benchmark.subprocess, "Popen", launch)
    if interrupted:
        monkeypatch.setattr(benchmark, "tree_rss", interrupt)
    try:
        with pytest.raises(KeyboardInterrupt if interrupted else RuntimeError):
            benchmark.run_backend(database, "cpp", 1, 0, 0)
        assert processes[0].poll() is not None
        with sqlite3.connect(database) as db:
            status, error = db.execute("SELECT status,error FROM runs").fetchone()
            expected = initial_status if initial_status in ("complete", "failed") else ("interrupted" if interrupted else "failed")
            assert status == expected
            if initial_status == "failed":
                assert error == "saved failure"
        child_pid = int(child_pid_file.read_text())
        deadline = time.monotonic() + 5
        while time.monotonic() < deadline:
            try:
                state = Path(f"/proc/{child_pid}/stat").read_text().split()[2]
            except FileNotFoundError:
                break
            if state == "Z":
                break  # exited, awaiting its new parent's reap
            time.sleep(0.01)
        else:
            pytest.fail("Owned verifier child remained running after interruption")
    finally:
        for process in processes:
            try:
                os.killpg(process.pid, 9)
            except ProcessLookupError:
                pass
            process.wait()
