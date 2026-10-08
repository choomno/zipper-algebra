#!/usr/bin/env python3
"""Closed loop: sample -> lift -> compose -> symbolic check, growing the grid
until the check passes.

**The pass criterion is the symbolic identity check**, not any proxy: the composed
cofactor must satisfy

    sum_j c_j(q,t) * L_j * R_{i_j} * R'_j  ==  X_3

as an identity in the free algebra Q(q,t)<x1,x2,x3,y1,y2,y3>, with q,t
transcendental. Nothing is shipped that has not passed it.

A failure has TWO possible causes and the loop must not confuse them:

  (a) UNDER-SAMPLED. The held-out check fails: fits agree on every point shown and
      mispredict unseen ones. More points help; the loop returns to sampling.

  (b) WRONG HYPOTHESIS CLASS. The held-out check PASSES and the symbolic identity
      still fails. Then the fits predict points they never saw and are *still*
      wrong in the free algebra -- so the true coefficient lies outside the model:
      a denominator beyond the six atoms {q, t, q-1, q+1, t-q, t+q^2}, or a
      bidegree above fit_rational's cap of 5, or a fault in composition rather
      than the lift. **More points will never fix this.** One can go on sampling
      points that agree with the fit so far and never obtain a passing cofactor.
      The loop stops and says which knob to widen.

      Since 2026-09-17 the ATOM half of the class is no longer a knob: the lifter
      infers the atoms from the sampled denominators (tools/cofactor/atoms.py) and
      reports `unexplained > 0` when a denominator has a factor outside its search
      box. What remains for (b) is the numerator bidegree cap, the search box, and
      composition.

Why this exists rather than a fixed grid: `run_lift_grid.py` carries a hard-coded
37-point grid, and the abandoned "unified" cofactor was fit on exactly those 37
points and disagreed with the exact answer on all 695 slots. Fitting is not
verification -- with fewer points than unknowns, many polynomials pass through
every point shown. The number of points needed is not known in advance, so it is
discovered.

This is orchestration the pipeline did not have; the phases themselves are the
existing tools (`extract_cofactor_qt.py`, `symbolic_lift_x3.py`,
`compose_x3_cofactor_qt.py`, `symbolic_identity_check.py`) and are not
reimplemented here.

Usage:
    lift_until_verified.py                      # grow from what is on disk
    lift_until_verified.py --start 60 --step 40 --max 260
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import glob
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent                      # zipper-algebra/tools
ZA = TOOLS.parent                        # zipper-algebra
REPO = ZA.parent
sys.path.insert(0, str(REPO))
LOGS = ZA / "logs_and_traces" / "f4ncgb"

PROOF_RE = re.compile(r"proof-n4-Q-q(\d+)-t(\d+)-noex\.txt$")


def have_points() -> set[tuple[int, int]]:
    out = set()
    for f in glob.glob(str(LOGS / "proof-n4-Q-q*-t*-noex.txt")):
        m = PROOF_RE.search(f)
        if m:
            out.add((int(m.group(1)), int(m.group(2))))
    return out


def candidate_points():
    """(q,t) candidates, in the order the pipeline's own grid prefers.

    Filters documented in run_lift_grid.py: q not in {2,4} (coefficient collapse
    on that hypersurface) and t != q (t-q is a denominator atom, so the diagonal
    is degenerate).
    """
    qs = [3, 5, 7, 8, 9, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    ts = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    for t in ts:
        for q in qs:
            if q in (2, 4) or t == q:
                continue
            yield (q, t)


def generate(points, jobs: int) -> int:
    """Run the existing per-point extractor at each (q,t)."""
    def one(qt):
        q, t = qt
        r = subprocess.run([sys.executable, str(HERE / "extract_cofactor_qt.py"),
                            "--q", str(q), "--t", str(t)],
                           capture_output=True, text=True)
        return qt, r.returncode
    ok = 0
    # Threads, not processes: each task is a subprocess call, so the GIL is not the
    # constraint, and a nested function cannot be pickled for a process pool.
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        for qt, rc in ex.map(one, points):
            if rc == 0:
                ok += 1
            else:
                print(f"    point {qt} failed (rc={rc})")
    return ok


def run(cmd, **kw):
    return subprocess.run(cmd, cwd=str(HERE), capture_output=True, text=True, **kw)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=60, help="points to reach before the first attempt")
    ap.add_argument("--step", type=int, default=40, help="points added per failed attempt")
    ap.add_argument("--max", type=int, default=260, help="give up beyond this many points")
    ap.add_argument("--jobs", type=int, default=None)
    a = ap.parse_args()

    try:
        from repo_settings import settings
        jobs = a.jobs or int(settings().limit("jobs"))
    except Exception:
        jobs = a.jobs or 6

    target = a.start
    while True:
        have = have_points()
        if len(have) < target:
            want = [p for p in candidate_points() if p not in have][: target - len(have)]
            print(f"[sample] have {len(have)}, generating {len(want)} more "
                  f"(jobs={jobs}, ~12.4 s each) ...", flush=True)
            generate(want, jobs)
            have = have_points()
        print(f"[sample] {len(have)} grid points on disk")

        print("[lift]   symbolic_lift_x3.py ...", flush=True)
        r = run([sys.executable, str(HERE / "symbolic_lift_x3.py")])
        if r.returncode == 3:
            last = [l for l in r.stdout.splitlines() if "NEED_MORE_SAMPLES" in l]
            print(f"[lift]   {last[-1] if last else 'needs more samples'}")
            target += a.step
            if target > a.max:
                print(f"[stop]   exceeded --max {a.max} points without verifying"); return 1
            continue
        if r.returncode == 4:
            print("[lift]   WRONG HYPOTHESIS CLASS -- more points will not help.")
            print("   " + "\n   ".join(r.stdout.strip().splitlines()[-4:]))
            return 4
        if r.returncode != 0:
            print("[lift]   FAILED:"); print(r.stdout[-1500:]); print(r.stderr[-800:]); return 2

        print("[compose] compose_x3_cofactor_qt.py ...", flush=True)
        r = run([sys.executable, str(HERE / "compose_x3_cofactor_qt.py")])
        if r.returncode != 0:
            print("[compose] FAILED:"); print(r.stdout[-1500:]); print(r.stderr[-800:]); return 2

        print("[verify]  symbolic_identity_check.py  <-- THE pass criterion", flush=True)
        r = run([sys.executable, str(HERE / "symbolic_identity_check.py")])
        if r.returncode == 0 and "PASS" in r.stdout:
            print("\n" + "\n".join(r.stdout.strip().splitlines()[-8:]))
            print(f"\n[done] verified with {len(have_points())} grid points.")
            return 0

        # TWO DISTINCT FAILURE MODES -- they need different escalations.
        #
        # (a) held-out verification failed  -> under-sampled. More points help.
        # (b) held-out verification PASSED and the symbolic identity still fails
        #     -> the fits predict unseen points and are still wrong, so the
        #     hypothesis class is wrong, not the sample count: the true
        #     coefficient is outside {six denominator atoms} x {bidegree <= 5},
        #     or the composition is at fault. MORE POINTS WILL NEVER FIX THIS --
        #     one can keep sampling points that agree with the fit so far and
        #     never obtain a passing cofactor. Stop and say so.
        print("[verify]  FAILED")
        print("   " + "\n   ".join(r.stdout.strip().splitlines()[-6:]))
        print()
        print("   Held-out verification PASSED but the symbolic identity did not.")
        print("   That is NOT an under-sampling failure: these fits predict points they")
        print("   never saw and are still wrong in the free algebra. Sampling more points")
        print("   cannot fix a wrong hypothesis class. Escalate the CLASS instead:")
        print("     - fit_rational(max_pq, max_pt) is capped at 5 -- raise it;")
        print("     - DENOMS is built from the six atoms {q, t, q-1, q+1, t-q, t+q^2};")
        print("       a true denominator outside that set can never be fit;")
        print("     - or the fault is in composition, not the lift.")
        print("   Re-run with a widened class rather than a larger grid.")
        return 4


if __name__ == "__main__":
    raise SystemExit(main())
