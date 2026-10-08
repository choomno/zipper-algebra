"""Driver for one survey point: reconstruct the full generic Q(q,t) basis (atoms) and decide
(non-)semisimplicity (trace form), then save a structured result. Resumable (skips if the result
file already exists). CLI:

    python3 -m tools.survey.run --untwist linear --n 4 --outdir results [--no-basis]

Outputs <outdir>/<spec>.json (summary) and, unless --no-basis, <outdir>/<spec>.gb.pkl (full GB).
"""
from __future__ import annotations
import argparse, json, os, pickle, time
import sympy as sp
from .spec import Spec, UNTWIST
from . import semisimple as ss
from . import reconstruct as rec


def run_spec(spec: Spec, scratch: str, outdir: str, do_basis: bool = True,
             char: int = ss.DEFAULT_CHAR, grid=None, log=print):
    os.makedirs(outdir, exist_ok=True)
    summ_path = os.path.join(outdir, f"{spec.name}.json")
    if os.path.exists(summ_path):
        log(f"[{spec.name}] already done -> {summ_path}")
        return json.load(open(summ_path))

    result = {"spec": spec.name, "untwist": spec.untwist, "n": spec.n,
              "max_untwist_level": spec.max_untwist_level}
    t0 = time.time()

    # --- semisimplicity (fast, mod-p; the survey's primary invariant) ---
    S = ss.semisimplicity(spec, scratch, char=char)
    result.update({"dim": S["dim"], "rank": S["rank"], "radical_dim": S["radical_dim"],
                   "semisimple": S["semisimple"], "trace_form_symmetric": S["symmetric"],
                   "ss_points": S["points"], "ss_char": char})
    log(f"[{spec.name}] dim={S['dim']} rank={S['rank']} semisimple={S['semisimple']} "
        f"radical={S['radical_dim']}")

    # --- full basis reconstruction over Q(q,t) via atoms ---
    if do_basis:
        R = rec.reconstruct_gb(spec, scratch, grid=grid, log=log)
        ok, msg = rec.validate(spec, scratch, R)
        result.update({"atoms": R["atoms"], "gb_size": len(R["LMS"]),
                       "basis_validated": ok, "basis_validation_msg": msg})
        with open(os.path.join(outdir, f"{spec.name}.gb.pkl"), "wb") as f:
            pickle.dump({"spec": spec.name, "n": spec.n,
                         "LMS": R["LMS"], "atoms": R["atoms"], "D": sp.srepr(R["D"]),
                         "GEN": {"|".join(lm): {"*".join(mon): sp.srepr(c) for mon, c in d.items()}
                                 for lm, d in R["GEN"].items()}}, f)
        log(f"[{spec.name}] atoms={R['atoms']} basis_validated={ok}")

    result["seconds"] = round(time.time() - t0, 1)
    json.dump(result, open(summ_path, "w"), indent=2)
    log(f"[{spec.name}] saved {summ_path}  [{result['seconds']}s]")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--untwist", required=True, choices=list(UNTWIST))
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--max-untwist-level", type=int, default=2)
    ap.add_argument("--outdir", default="survey_results")
    ap.add_argument("--scratch", default="/tmp")
    ap.add_argument("--char", type=int, default=ss.DEFAULT_CHAR,
                    help="trace-form char: large prime (fast, default) or 0 (exact-Q verify)")
    ap.add_argument("--no-basis", action="store_true", help="semisimplicity only, skip full-basis reconstruction")
    args = ap.parse_args()
    spec = Spec(args.untwist, args.n, args.max_untwist_level)
    run_spec(spec, args.scratch, args.outdir, do_basis=not args.no_basis, char=args.char)


if __name__ == "__main__":
    main()
