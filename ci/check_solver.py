"""Check native solver integration on a small, explicitly zero-slider presentation."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
TARGET = "X_2"
PRESENTATION = ("--n", "2", "--untwist", "none", "--zero-from", "2")
REVISION = "cb37d6570a42cf2f786cc96180b18456ae6d1bcf"


def isolated(*arguments, directory, timeout=120):
    env = dict(os.environ, PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1")
    env.pop("PYTHONPATH", None)
    result = subprocess.run([sys.executable, "-I", *arguments], cwd=directory,
                            env=env, text=True, capture_output=True, timeout=timeout)
    if result.returncode:
        raise AssertionError(f"{arguments}: exit {result.returncode}\n"
                             f"{result.stdout[-5000:]}\n{result.stderr[-2000:]}")
    return result.stdout


def check(binary, revision, directory, solver_io="files"):
    binary = Path(binary).resolve()
    directory = Path(directory).resolve()
    output, cache = directory / "production", directory / "cache"
    package = Path(isolated("-c", "import pathlib, zipper_algebra; "
                            "print(pathlib.Path(zipper_algebra.__file__).resolve())",
                            directory=directory).strip())
    assert not package.is_relative_to(ROOT), "the solver check requires an installed wheel"
    description = json.loads(isolated("-m", "zipper_algebra", "describe", *PRESENTATION,
                                     "--json", directory=directory))
    relation_count = len(description["relations"])
    # X_2=0 is explicitly part of this small test presentation. The test
    # exercises native I/O, proof validation, coefficient reconstruction and
    # both certificate checks; it is not a research computation.
    text = isolated("-m", "zipper_algebra", "discover", *PRESENTATION,
                    "--target", TARGET, "--solver", str(binary), "--output", str(output),
                    "--cache", str(cache), "--anchor", "5,7", "--start", "48",
                    "--step", "12", "--max", "96", "--jobs", "1",
                    "--timeout", "180", "--max-iterations", "100000",
                    "--solver-revision", revision, "--solver-io", solver_io, "--json",
                    directory=directory, timeout=900)
    report = json.loads(text)
    assert report == json.loads((output / "result.json").read_text())
    assert report["status"] == "verified" and report["verified"] is True, report
    assert report["basis_complete"] is False
    assert report['producer']['origin'] == 'build'
    assert report['producer']['changes_since_build'] == []
    assert report["exact"]["pass"] is True
    assert report["target_name"] == TARGET and report["runs"]
    assert len(report["presentation"]["relation_labels"]) == relation_count
    certificate = Path(report["certificate"])
    verified = json.loads(isolated("-m", "zipper_algebra", "verify", str(certificate),
                                    *PRESENTATION, "--target", TARGET, "--pit", "--json",
                                    directory=directory))
    assert verified["pass"] and verified["exact"]["pass"] and verified["pit"]["pass"]
    for run in report["runs"]:
        sample = json.loads(isolated(
            "-c", "import json,sys; from zipper_algebra.execution.store import read; "
            "files,record=read(sys.argv[1],sys.argv[2],generation=sys.argv[3]); "
            "print(json.dumps({'record':record,'input':files['input.ms'].decode('utf-8')}))",
            run['store'], run['key'], run['generation'], directory=directory))
        record = sample['record']
        assert record["status"] == "proof-validated" and record["basis_complete"] is False
        request = record["request"]
        assert request['producer'] == report['producer']
        assert request['solver']['io_mode'] == solver_io
        assert request["solver"]["revision"] == revision
        assert request["solver"]["prefix"] == []
        assert request["point"][:2] == run["point"]
        source = sample['input'].encode('utf-8')
        assert hashlib.sha256(source).hexdigest() == request["input_sha256"]
        assert len(source.splitlines()) - 2 == relation_count
        expected = directory / "expected-input.ms"
        q, t = run["point"]
        isolated("-m", "zipper_algebra", "emit", *PRESENTATION, "--q", str(q),
                 "--t", str(t), "--out", str(expected), "--json", directory=directory)
        assert source == expected.read_bytes(), "a solver input changed presentation"
    summary = {"status": "verified", "solver_revision": revision, "solver_io": solver_io,
               "installed_package": str(package), "target": TARGET,
               "relation_count": relation_count, "sample_points": report["points"],
               "cofactor_terms": report["tuples"], "basis_complete": False,
               "exact": verified["exact"], "pit": verified["pit"]}
    (directory / "check-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--revision", default=REVISION)
    parser.add_argument("--solver-io", choices=("files", "sqlite"), default="files")
    parser.add_argument("--output", type=Path, help="new directory retaining the complete run")
    args = parser.parse_args(argv)
    if args.output is None:
        with tempfile.TemporaryDirectory(prefix="zipper-real-solver-") as directory:
            check(args.binary, args.revision, directory, args.solver_io)
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        check(args.binary, args.revision, args.output, args.solver_io)


if __name__ == "__main__":
    main()
