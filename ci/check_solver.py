"""Exercise fresh discovery with the pinned real solver and an installed wheel."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
REVISION = "789504d7c7da6af6946faecb73421b037766a5b8"


def isolated(*arguments, directory, timeout=120):
    env = dict(os.environ, PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1")
    env.pop("PYTHONPATH", None)
    result = subprocess.run([sys.executable, "-I", *arguments], cwd=directory,
                            env=env, text=True, capture_output=True, timeout=timeout)
    if result.returncode:
        raise AssertionError(f"{arguments}: exit {result.returncode}\n"
                             f"{result.stdout[-5000:]}\n{result.stderr[-2000:]}")
    return result.stdout


def check(binary, revision, directory, *, n=3):
    target = "X_3" if n == 4 else "(Theta_3-q^2)X_3"
    start, step, maximum = (100, 40, 260) if n == 4 else (12, 12, 96)
    strands = str(n)
    binary = Path(binary).resolve()
    directory = Path(directory).resolve()
    output, cache = directory / "production", directory / "cache"
    package = Path(isolated("-c", "import pathlib, zipper_algebra; "
                            "print(pathlib.Path(zipper_algebra.__file__).resolve())",
                            directory=directory).strip())
    assert not package.is_relative_to(ROOT), "the solver check requires an installed wheel"
    description = json.loads(isolated("-m", "zipper_algebra", "describe", "--n", strands,
                                     "--json", directory=directory))
    relation_count = len(description["relations"])
    # Neither target is added to the defining presentation; no saved lift or
    # certificate is supplied. The n=4 case uses the original 31 relations.
    if n == 4:
        assert relation_count == 31
        assert description["presentation_hash"] == "091bbb92271ab870"
    started = time.monotonic()
    text = isolated("-m", "zipper_algebra", "discover", "--n", strands,
                    "--target", target, "--solver", str(binary), "--output", str(output),
                    "--cache", str(cache), "--anchor", "5,7", "--start", str(start),
                    "--step", str(step), "--max", str(maximum), "--jobs", "1",
                    "--timeout", "180", "--max-iterations", "100000",
                    "--solver-revision", revision, "--json",
                    directory=directory, timeout=1800 if n == 4 else 900)
    discovery_seconds = time.monotonic() - started
    report = json.loads(text)
    assert report == json.loads((output / "result.json").read_text())
    assert report["status"] == "verified" and report["verified"] is True, report
    assert report["basis_complete"] is False
    assert report["exact"]["pass"] is True
    assert report["target_name"] == target and report["runs"]
    assert len(report["presentation"]["relation_labels"]) == relation_count
    certificate = Path(report["certificate"])
    verified = json.loads(isolated("-m", "zipper_algebra", "verify", str(certificate),
                                    "--n", strands, "--target", target, "--pit", "--json",
                                    directory=directory, timeout=600))
    assert verified["pass"] and verified["exact"]["pass"] and verified["pit"]["pass"]
    for run in report["runs"]:
        record_path = Path(run["record"])
        record = json.loads(record_path.read_text())
        assert record["status"] == "proof-validated" and record["basis_complete"] is False
        request = record["request"]
        assert request["solver"]["revision"] == revision
        assert request["solver"]["prefix"] == []
        assert request["point"][:2] == run["point"]
        source = record_path.parent / record["generation"] / "input.ms"
        assert hashlib.sha256(source.read_bytes()).hexdigest() == request["input_sha256"]
        assert len(source.read_text().splitlines()) - 2 == relation_count
        expected = directory / "expected-input.ms"
        q, t = run["point"]
        isolated("-m", "zipper_algebra", "emit", "--n", strands, "--q", str(q),
                 "--t", str(t), "--out", str(expected), "--json", directory=directory)
        assert source.read_bytes() == expected.read_bytes(), "a solver input changed presentation"
    summary = {"status": "verified", "solver_revision": revision,
               "installed_package": str(package), "n": n, "target": target,
               "discovery_seconds": discovery_seconds,
               "pipeline_seconds": report["seconds"],
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
    parser.add_argument("--n", type=int, choices=(3, 4), default=3,
                        help="3: small integration case; 4: fresh X_3 discovery")
    parser.add_argument("--output", type=Path, help="new directory retaining the complete run")
    args = parser.parse_args(argv)
    if args.output is None:
        with tempfile.TemporaryDirectory(prefix="zipper-real-solver-") as directory:
            check(args.binary, args.revision, directory, n=args.n)
    else:
        args.output.mkdir(parents=True, exist_ok=False)
        check(args.binary, args.revision, args.output, n=args.n)


if __name__ == "__main__":
    main()
