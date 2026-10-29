"""Run specialized survey observations and aggregate explicitly qualified records.

Run one preset with --index I, every preset with --all, or read existing records
with --aggregate. An invocation recomputes the requested observations.
"""
from __future__ import annotations
import argparse, glob, json, os
from .spec import Spec
from . import run as _run
from . import semisimple as _ss

# The survey grid. Extend untwistings / levels here; add higher-untwisting variants later.
UNTWISTS = ["linear", "square-zero", "framing"]
LEVELS = [3, 4, 5]


def enumerate_specs():
    return [Spec(u, n) for u in UNTWISTS for n in LEVELS]


def aggregate(outdir):
    rows = []
    for p in sorted(glob.glob(os.path.join(outdir, "*.json"))):
        rows.append(json.load(open(p)))
    return rows


def _fmt(rows):
    header = f"{'spec':22s} {'point':>12s} {'field':>12s} {'words':>7s} {'trace-rank':>10s} {'status':>20s}"
    output = [header, "-" * len(header)]
    for record in sorted(rows, key=lambda row: row.get("spec", "")):
        if record.get("schema") != _run.SURVEY_SCHEMA:
            output.append(f"{record.get('spec', '?'):22s} {'?':>12s} {'?':>12s} "
                          f"{'?':>7s} {'?':>10s} {'legacy-unqualified':>20s}")
            continue
        trace = record["trace_observations"]
        for point in trace["points"]:
            pair = f"({point['q']},{point['t']})"
            output.append(f"{record['spec']:22s} {pair:>12s} {trace['coefficient_domain']:>12s} "
                          f"{point['standard_word_count']:>7} {point['trace_rank']:>10} {'observed':>20s}")
    return "\n".join(output)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--index", type=int, default=None, help="run only the i-th spec (SLURM array task)")
    ap.add_argument("--all", action="store_true", help="recompute every spec locally")
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--outdir", default="survey_results")
    ap.add_argument("--scratch", default="/tmp")
    ap.add_argument("--char", type=int, default=_ss.DEFAULT_CHAR)
    ap.add_argument("--no-basis", action="store_true")
    args = ap.parse_args()
    specs = enumerate_specs()

    if args.list:
        for i, s in enumerate(specs):
            print(f"{i:3d}  {s.name}")
        return
    if args.aggregate:
        print(_fmt(aggregate(args.outdir)))
        return
    todo = [specs[args.index]] if args.index is not None else (specs if args.all else [])
    if not todo:
        ap.error("choose --list, --index I, --all, or --aggregate")
    for s in todo:
        try:
            _run.run_spec(s, args.scratch, args.outdir, do_basis=not args.no_basis, char=args.char)
        except Exception as e:
            print(f"[{s.name}] FAILED: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
