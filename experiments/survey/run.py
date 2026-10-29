"""Record trace-matrix observations and optional basis-coefficient candidates.

Each invocation recomputes its observations. Existing summaries are historical
records, not completion certificates or a reason to silently skip a new run.
"""
from __future__ import annotations
import argparse
from experiments.output import (prepare_output_dir, capture_experiment_provenance,
                                save_pickle_output)
from zipper_algebra.provenance import write_output_provenance
import json
import os
from pathlib import Path
import time
from .spec import Spec, UNTWIST
from . import semisimple as ss

SURVEY_SCHEMA = "zipper-experimental-survey-v2"


def run_spec(spec: Spec, scratch: str, outdir: str, do_basis: bool = True,
             char: int = ss.DEFAULT_CHAR, grid=None, log=print):
    prepare_output_dir(outdir)
    family = Path(__file__).resolve().parents[1]
    producer = capture_experiment_provenance(
        __file__, family / "research/theta_min_poly.py",
        family / "variants/square_zero_regular_rep_qt.py")
    summary_path = os.path.join(outdir, f"{spec.name}.json")
    result = {"schema": SURVEY_SCHEMA, "status": "observed", "basis_complete": False,
              "spec": spec.name, "untwist": spec.untwist, "n": spec.n,
              "max_untwist_level": spec.max_untwist_level}
    started = time.time()
    trace = ss.semisimplicity(spec, scratch, char=char)
    result["trace_observations"] = trace
    log(f"[{spec.name}] observed {len(trace['points'])} trace matrices over "
        f"{trace['coefficient_domain']}")

    if do_basis:
        import sympy as sp
        from . import reconstruct as rec
        candidate = rec.reconstruct_gb(spec, scratch, grid=grid, log=log)
        passed, message = rec.validate(spec, scratch, candidate)
        artifact = os.path.join(outdir, f"{spec.name}.candidate-basis.pkl")
        result["basis_candidate"] = {
            "status": "candidate", "basis_complete": False,
            "atoms": candidate["atoms"], "size": len(candidate["LMS"]),
            "sample_check_passed": passed, "sample_check_message": message,
            "artifact": artifact,
        }
        save_pickle_output(artifact,
                           {"status": "candidate", "basis_complete": False,
                            "spec": spec.name, "n": spec.n, "LMS": candidate["LMS"],
                            "atoms": candidate["atoms"], "D": sp.srepr(candidate["D"]),
                            "GEN": {"|".join(lm): {"*".join(mon): sp.srepr(c)
                                                    for mon, c in terms.items()}
                                    for lm, terms in candidate["GEN"].items()}},
                           producer=producer)
        log(f"[{spec.name}] candidate coefficient sample check={passed}")

    result["seconds"] = round(time.time() - started, 1)
    with open(summary_path, "w") as stream:
        json.dump(result, stream, indent=2)
    write_output_provenance(summary_path, producer=producer)
    log(f"[{spec.name}] saved observations -> {summary_path}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--untwist", required=True, choices=list(UNTWIST))
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--max-untwist-level", type=int, default=2)
    parser.add_argument("--outdir", default="runs/survey")
    parser.add_argument("--scratch", default="/tmp")
    parser.add_argument("--char", type=int, default=ss.DEFAULT_CHAR,
                        help="trace-matrix field characteristic: a prime, or 0 for Q")
    parser.add_argument("--no-basis", action="store_true", help="trace observations only")
    args = parser.parse_args()
    spec = Spec(args.untwist, args.n, args.max_untwist_level)
    run_spec(spec, args.scratch, args.outdir, do_basis=not args.no_basis, char=args.char)


if __name__ == "__main__":
    main()
