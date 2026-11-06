"""Historical executable provenance and paired benchmark acceptance checks."""
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import tarfile
from types import SimpleNamespace

import pytest

from experiments import zipper_native_epoch_benchmark as bench
from zipper_algebra.zipping_algorithm.interpreter import Limits


def archive_epoch(tmp_path, *, duplicate=False):
    directory = tmp_path / "build"
    directory.mkdir()
    binary = directory / "zipper-engine"
    binary.write_bytes(b"test executable")
    source = b"native test source"
    name = "core.cpp"
    digest = hashlib.sha256(source).hexdigest()
    archive = tmp_path / "source.tar.gz"
    with tarfile.open(archive, "w:gz") as saved:
        for _ in range(2 if duplicate else 1):
            member = tarfile.TarInfo(bench.PREFIX + name)
            member.size = len(source)
            saved.addfile(member, io.BytesIO(source))
    identity = tmp_path / "source-identity.json"
    identity.write_text(json.dumps({"archive": archive.name, "archive_sha256": bench.sha(archive),
                                    "sources": {"native": {bench.PREFIX + name: digest},
                                                "python": bench.memoized._source()}}))
    binary.with_suffix(".json").write_text(json.dumps({"abi": bench.ABI,
                                                      "binary_sha256": bench.sha(binary),
                                                      "source_sha256": {name: digest}}))
    return {"label": "initial", "binary": str(binary), "source_identity": str(identity)}


def test_historical_archive_matches_each_source_and_binary(tmp_path):
    spec = archive_epoch(tmp_path)
    result = bench.validate_epoch(spec)
    assert result["build"]["source_sha256"] == {
        "core.cpp": hashlib.sha256(b"native test source").hexdigest()}
    assert result["source_identity_sha256"] == bench.sha(Path(spec["source_identity"]))


@pytest.mark.parametrize("target", ["binary", "archive", "source_record", "python"])
def test_tampered_epoch_is_refused(tmp_path, target):
    spec = archive_epoch(tmp_path)
    if target == "binary":
        Path(spec["binary"]).write_bytes(b"changed")
    elif target == "archive":
        (tmp_path / "source.tar.gz").write_bytes(b"changed")
    else:
        path = Path(spec["source_identity"])
        value = json.loads(path.read_text())
        value["sources"]["native" if target == "source_record" else "python"] = {}
        path.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        bench.validate_epoch(spec)


def test_duplicate_archive_source_is_refused(tmp_path):
    with pytest.raises(ValueError, match="duplicate"):
        bench.validate_epoch(archive_epoch(tmp_path, duplicate=True))


def test_cpu_stat_parser_handles_spaces_and_parentheses(monkeypatch):
    fields = ["S"] + ["0"] * 49
    fields[11], fields[12], fields[19] = "17", "4", "999"
    monkeypatch.setattr(Path, "read_text", lambda self: "42 (a (word) name) " + " ".join(fields))
    assert bench.cpu_snapshot(42) == (21, 999)


def test_exact_decode_and_strict_comparison():
    huge = "9" * 10000
    result = bench.decode({"coefficient": [[-4, huge]], "child": [{"coefficient": [[0, "-1"]]}]})
    assert result["coefficient"][0][1] == 10**10000 - 1
    assert result["child"][0]["coefficient"] == [[0, -1]]
    assert not bench.equal({"flag": False}, {"flag": 0})
    assert not bench.equal([1], [1.0])
    assert bench.arrays({"ports": (1, 2)}) == {"ports": [1, 2]}
    for value in (True, "", "-", "1a", "+1"):
        with pytest.raises(ValueError):
            bench.decimal_integer(value)


def test_failed_request_permanently_closes_session():
    class Writer:
        def write(self, value):
            raise KeyboardInterrupt
    session = object.__new__(bench.EpochSession)
    session.process = SimpleNamespace(stdin=Writer(), stdout=io.StringIO())
    closed = []
    session.close = lambda: closed.append(True)
    with pytest.raises(KeyboardInterrupt):
        session.request({"op": "reduce"})
    assert closed == [True]


def mocked_sessions(monkeypatch, *, wrong=False):
    expected = {"result": [1, True, "exact"], "coefficient": [[0, 1]]}
    references = []
    class Oracle:
        def __init__(self):
            self._catalog = SimpleNamespace(configure_rule_cache=lambda **kwargs: None)
        def reduce(self, *args, **kwargs):
            references.append(kwargs)
            return deepcopy(expected)
    class Session:
        def __init__(self, spec, *, cpu=None):
            self.affinity = [cpu] if cpu is not None else [0]
            self.identity = {"specification": spec}
            self.startup_seconds = 0.01
            self.pid = 42
        def request(self, request):
            return {}
        def reduce(self, *args):
            actual = deepcopy(expected)
            if wrong:
                actual["result"][1] = 1
            return actual, {"seconds": {"total_seconds": 0.1}}
        def close(self):
            pass
    monkeypatch.setattr(bench.collected, "CollectedSession", Oracle)
    monkeypatch.setattr(bench, "EpochSession", Session)
    monkeypatch.setattr(bench, "validate_epoch", lambda spec: {"specification": spec})
    ticks = iter(range(0, 1000, 2))
    monkeypatch.setattr(bench, "cpu_snapshot", lambda pid: (next(ticks), 99))
    monkeypatch.setattr(bench, "memory_snapshot", lambda pid: {"VmRSS_kib": 10})
    return references


def test_three_epochs_each_occupy_each_position_and_share_one_reference(tmp_path, monkeypatch):
    references = mocked_sessions(monkeypatch)
    specs = [{"label": label} for label in "abc"]
    database = tmp_path / "results.sqlite3"
    campaign = bench.run(database, specs, cases=["smoke"], rounds=3, limits=Limits())
    assert len(references) == 1
    with sqlite3.connect(database) as db:
        assert db.execute("select status from native_epoch_campaigns").fetchone() == ("complete",)
        rows = db.execute("select epoch,position,count(*) from native_epoch_measurements group by epoch,position").fetchall()
        assert len(rows) == 9 and all(row[2] == 1 for row in rows)
        assert db.execute("select count(distinct certificate_sha256) from native_epoch_measurements").fetchone() == (1,)
        assert db.execute("select min(native_cpu_ticks) from native_epoch_measurements").fetchone() == (2,)
        assert db.execute("select count(*) from native_epoch_references where campaign_id=?", (campaign,)).fetchone() == (1,)


def test_wrong_scalar_type_stops_before_measurement_insert(tmp_path, monkeypatch):
    mocked_sessions(monkeypatch, wrong=True)
    database = tmp_path / "results.sqlite3"
    with pytest.raises(AssertionError, match="exact reference"):
        bench.run(database, [{"label": "wrong"}], cases=["smoke"], rounds=1, limits=Limits())
    with sqlite3.connect(database) as db:
        assert db.execute("select count(*) from native_epoch_measurements").fetchone() == (0,)
        assert db.execute("select status from native_epoch_campaigns").fetchone() == ("failed",)


def test_harness_source_drift_stops_before_measurement_insert(tmp_path, monkeypatch):
    mocked_sessions(monkeypatch)
    original = bench.EpochSession.reduce
    changed = False
    def reduce_and_change(self, *args):
        nonlocal changed
        result = original(self, *args)
        changed = True
        return result
    monkeypatch.setattr(bench.EpochSession, "reduce", reduce_and_change)
    monkeypatch.setattr(bench, "sha", lambda path: "new" if changed else "old")
    database = tmp_path / "results.sqlite3"
    with pytest.raises(RuntimeError, match="source changed during measurement"):
        bench.run(database, [{"label": "changed"}], cases=["smoke"], rounds=1, limits=Limits())
    with sqlite3.connect(database) as db:
        assert db.execute("select count(*) from native_epoch_measurements").fetchone() == (0,)
        assert db.execute("select status from native_epoch_campaigns").fetchone() == ("failed",)


def test_requested_affinity_is_recorded_and_invalid_cpu_rejected(tmp_path, monkeypatch):
    mocked_sessions(monkeypatch)
    monkeypatch.setattr(bench.os, "sched_getaffinity", lambda pid: {3, 7})
    database = tmp_path / "results.sqlite3"
    bench.run(database, [{"label": "pinned"}], cases=["smoke"], rounds=1, limits=Limits(), cpu=7)
    with sqlite3.connect(database) as db:
        config, setup = db.execute("select config_json,setup_json from native_epoch_campaigns").fetchone()
        assert json.loads(config)["native_child_cpu"] == 7
        assert json.loads(setup)["pinned"]["native_affinity"] == [7]
    for invalid in (2, True):
        with pytest.raises(ValueError, match="allowed integer CPU"):
            bench.run(database, [{"label": "invalid"}], cases=["smoke"], rounds=1, limits=Limits(), cpu=invalid)
