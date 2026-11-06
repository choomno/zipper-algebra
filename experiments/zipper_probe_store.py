"""Import selected long-word certificates into an existing exhaustive SQLite store.

The selected probes are kept separate from the exhaustive prefix. By default
certificates are replayed on import. --trust-recorded-replay accepts an explicitly
prechecked report without repeating its computation and records that provenance;
it is intended for reports produced and checked in the same trusted session.
Uncheckpointed attempts are retained as metadata, never as checked probe rows.
This Unix experiment uses the parallel runner's campaign flock. Do not run the
sequential driver concurrently on the same database.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import sqlite3

from zipper_algebra.zipping_algorithm import algorithm
from zipper_algebra.zipping_algorithm.catalog import Catalog
from zipper_algebra.zipping_algorithm.enumeration import _json, _sha


def import_probes(database, report_path, *, recorded_replay=False):
    database, report_path = Path(database), Path(report_path)
    raw = report_path.read_bytes()
    report = json.loads(raw)
    if report.get("schema") != "zipper-selected-long-word-probes-v1" or report.get("exhaustive") is not False:
        raise ValueError("an explicitly nonexhaustive selected-probe report is required")
    report_sha = hashlib.sha256(raw).hexdigest()
    importer_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    source_dir = Path(algorithm.__file__).parent
    for name, digest in report["source_sha256"].items():
        if (name != Path(name).name or not name.endswith(".py")
                or hashlib.sha256((source_dir / name).read_bytes()).hexdigest() != digest):
            raise ValueError("probe report source binding does not match the checkout")
    method = "recorded_replay_before_export" if recorded_replay else "replayed_on_import"
    session = None if recorded_replay else algorithm.AlgorithmSession()
    with (database.parent / ".parallel.lock").open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError("another runner owns this campaign") from error
        with sqlite3.connect(database.resolve().as_uri() + "?mode=rw", uri=True) as db:
            identity = json.loads(db.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()[0])
            expected = identity["algorithm_source_sha256"]
            if algorithm._source() != expected or Catalog.bundled().sha256 != identity["catalog_sha256"]:
                raise ValueError("main campaign source or catalog differs from the checkout")
            if any(report["source_sha256"].get(name) != digest for name, digest in expected.items()):
                raise ValueError("probe report differs from the main campaign source")
            records = []
            keys = set()
            for probe in report["probes"]:
                certificate = probe["certificate"]
                encoded = _json(certificate)
                if (certificate["source_sha256"] != expected
                        or certificate["catalog_sha256"] != identity["catalog_sha256"]
                        or certificate["n"] != probe["n"]
                        or certificate["input_word"] != probe["word"]
                        or certificate["status"] != probe["status"]
                        or certificate["metrics"] != probe["metrics"]
                        or certificate["limits"] != probe["limits"]
                        or probe["length"] != len(probe["word"])
                        or certificate["coefficient_ring"] != identity["coefficient_ring"]
                        or certificate["trace_enabled"] != report.get("trace", False)
                        or certificate["status"] not in ("complete", "incomplete", "control_obligation")):
                    raise ValueError("probe certificate identity or summary mismatch")
                if recorded_replay:
                    if probe.get("json_reload_replay_checked") is not True:
                        raise ValueError("recorded replay must be explicitly declared for every probe")
                elif session.check(json.loads(encoded)) is not True:
                    raise ValueError("probe certificate failed exact replay")
                key = (probe["name"], probe["n"], tuple(probe["word"]))
                if key in keys:
                    raise ValueError("duplicate selected probe in one report")
                keys.add(key)
                word = _json(probe["word"])
                records.append((report_sha, probe["name"], probe["n"], probe["length"], word,
                    _sha(_json({"n":probe["n"],"word":probe["word"]})), probe["status"], encoded,
                    _sha(encoded), _json(expected), identity["catalog_sha256"], _json(probe["limits"]),
                    _json(probe["metrics"]), probe["reduce_seconds"], probe["check_seconds"], method))
            if algorithm._source() != expected:
                raise ValueError("source changed while importing probes")
            metadata = {key:value for key,value in report.items() if key != "probes"}
            with db:
                db.execute("CREATE TABLE IF NOT EXISTS probe_imports (report_sha TEXT PRIMARY KEY, importer_sha TEXT NOT NULL, report_path TEXT NOT NULL, replay_method TEXT NOT NULL, metadata TEXT NOT NULL)")
                db.execute("CREATE TABLE IF NOT EXISTS long_word_probes (report_sha TEXT, name TEXT, n INTEGER, length INTEGER, word TEXT, input_sha TEXT, status TEXT, certificate TEXT, certificate_sha TEXT, algorithm_source TEXT, catalog_sha TEXT, limits TEXT, metrics TEXT, producer_seconds REAL, checker_seconds REAL, replay_method TEXT, PRIMARY KEY(report_sha,name,n,input_sha))")
                db.execute("INSERT OR IGNORE INTO probe_imports VALUES (?,?,?,?,?)",
                           (report_sha,importer_sha,str(report_path.resolve()),method,_json(metadata)))
                db.executemany("INSERT OR IGNORE INTO long_word_probes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", records)
            return {"report_sha256":report_sha,"importer_sha256":importer_sha,"probe_count":len(records),
                    "replay_method":method,"uncheckpointed_attempt_count":len(report.get("uncheckpointed_attempts",[])),
                    "database":str(database.resolve())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database")
    parser.add_argument("report")
    parser.add_argument("--trust-recorded-replay", dest="recorded_replay", action="store_true")
    args = parser.parse_args()
    print(json.dumps(import_probes(args.database,args.report,recorded_replay=args.recorded_replay),indent=2))


if __name__ == "__main__":
    main()
