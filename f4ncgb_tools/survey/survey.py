"""Mass-survey fan-out over a list of specs. Each spec is independent -> one cluster task each
(job array), or run locally in sequence. Resumable: run_spec skips specs whose result JSON exists.

CLI:
  python3 -m f4ncgb_tools.survey.survey --list                 # print the spec list (index: name)
  python3 -m f4ncgb_tools.survey.survey --index 5 --outdir R   # run one spec (for SLURM array)
  python3 -m f4ncgb_tools.survey.survey --all --outdir R       # run all locally (resumable)
  python3 -m f4ncgb_tools.survey.survey --aggregate --outdir R # collate results into a table
"""
from __future__ import annotations
import argparse, glob, json, os
from .spec import Spec, UNTWIST
from . import run as _run
from . import semisimple as _ss

# The survey grid. Extend untwistings / levels here; add higher-untwisting variants later.
UNTWISTS = ["linear", "jordan", "framing"]
LEVELS = [3, 4, 5]


def enumerate_specs():
    return [Spec(u, n) for u in UNTWISTS for n in LEVELS]


def aggregate(outdir):
    rows = []
    for p in sorted(glob.glob(os.path.join(outdir, "*.json"))):
        rows.append(json.load(open(p)))
    return rows


def _fmt(rows):
    hdr = f"{'spec':22s} {'dim':>5s} {'rank':>5s} {'radical':>7s} {'semisimple':>10s} {'atoms':>24s}"
    out = [hdr, "-" * len(hdr)]
    for r in sorted(rows, key=lambda r: (r.get("untwist"), r.get("n"))):
        out.append(f"{r['spec']:22s} {r.get('dim','?'):>5} {r.get('rank','?'):>5} "
                   f"{r.get('radical_dim','?'):>7} {str(r.get('semisimple','?')):>10s} "
                   f"{','.join(r.get('atoms',[])):>24s}")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--index", type=int, default=None, help="run only the i-th spec (SLURM array task)")
    ap.add_argument("--all", action="store_true", help="run every spec locally, resumable")
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
