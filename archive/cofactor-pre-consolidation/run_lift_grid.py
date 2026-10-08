#!/usr/bin/env python3
"""Run extract_cofactor_qt.py over a 30-point (q,t) grid, save outputs to
logs_and_traces/f4ncgb/lift_samples/, then build a consolidated JSON of
all (left, idx, right, coef) tuples ready for Q(q,t) interpolation.
"""
from __future__ import annotations
import concurrent.futures as cf
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
# cofactor/ -> tools/ -> zipper-algebra/ ; written when this lived at the repo root
REPO = os.path.dirname(os.path.dirname(HERE))
LOGS = os.path.join(REPO, "logs_and_traces", "f4ncgb")
SAMPLES = os.path.join(LOGS, "lift_samples")
EXTRACT = os.path.join(HERE, "extract_cofactor_qt.py")
SRC_PREFIX = os.path.join(LOGS, "cofactor-X3-n4-Q-")

os.makedirs(SAMPLES, exist_ok=True)

# 30-point grid. Filters:
#   q in {3,5,7,8,9,11,13}   (skip q in {2,4}: observed coefficient collapse)
#   t in {3,5,7,11,13,17,19}
#   t != q, t != q^2 (q^2 never falls in t-list so this is automatic)
#   Empirically: q in {2,4} hit a hypersurface where some coefficients vanish
#     (679 / 674 tuples) — kept out of the interpolation grid; documented below.
GRID = [
    (3, 5), (3, 7), (3, 11), (3, 13), (3, 17), (3, 19),
    (5, 3), (5, 7), (5, 11), (5, 13), (5, 17), (5, 19),
    (7, 3), (7, 5), (7, 11), (7, 13), (7, 17), (7, 19),
    (8, 3), (8, 5), (8, 11), (8, 13), (8, 17),
    (9, 3), (9, 5), (9, 7), (9, 11), (9, 13), (9, 17), (9, 19),
    (11, 3), (11, 5), (11, 7), (11, 19),
    (13, 3), (13, 7), (13, 19),
]
assert len(GRID) == 37, len(GRID)

# The 37-point grid above is the ORIGINAL one, and it is too small: the abandoned
# "unified" cofactor was fit from exactly these 37 samples and disagreed with the
# exact answer on all 695 slots. The production lift used 212 points, recorded in
# lift_samples/coefficients.json. --grid-from reads that list back.
def grid_from_coefficients(path):
    import json
    with open(path) as fh:
        d = json.load(fh)
    pts = sorted(tuple(int(x) for x in k.split(",")) for k in d["per_point_sizes"])
    return pts



def run_one(qt):
    q, t = qt
    out_default = os.path.join(LOGS, f"cofactor-X3-n4-Q-q{q}-t{t}.txt")
    out_dest = os.path.join(SAMPLES, f"cofactor-X3-n4-Q-q{q}-t{t}.txt")
    if os.path.exists(out_dest):
        return qt, "cached", out_dest, 0.0
    t0 = time.time()
    res = subprocess.run(
        [sys.executable, EXTRACT, "--q", str(q), "--t", str(t)],
        capture_output=True, text=True, cwd=REPO,
    )
    elapsed = time.time() - t0
    if res.returncode != 0:
        return qt, f"FAIL rc={res.returncode}: {res.stderr[-300:]}", None, elapsed
    # extract_cofactor_qt.py writes to LOGS; copy to lift_samples
    if os.path.exists(out_default):
        shutil.copy2(out_default, out_dest)
    # parse the verify line + tuple count
    m = re.search(r"verify:\s*(\S+);\s*tuples:\s*(\d+)", res.stdout)
    status = f"OK {m.group(2)} tuples" if m else "OK (no parse)"
    return qt, status, out_dest, elapsed


def parse_cofactor(path):
    """Return list of (left, idx, right, coef_str)."""
    tuples = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("*")]
            coef = parts[0]
            left = ""
            right = ""
            idx = None
            for i, p in enumerate(parts[1:], 1):
                if p.startswith("[") and p.endswith("]"):
                    idx = int(p[1:-1])
                    if i > 1:
                        left = "*".join(parts[1:i])
                    if i < len(parts) - 1:
                        right = "*".join(parts[i + 1:])
                    break
            assert idx is not None, f"unparsed: {line}"
            # canonicalise coef: drop leading '+'
            if coef.startswith("+"):
                coef = coef[1:]
            tuples.append((left, idx, right, coef))
    return tuples


def main():
    t_start = time.time()
    import argparse as _ap
    _p = _ap.ArgumentParser()
    _p.add_argument("--grid-from", default=None,
                    help="read the (q,t) list from a coefficients.json (the 212-point production grid)")
    _p.add_argument("--jobs", type=int, default=None)
    _a, _ = _p.parse_known_args()
    global GRID
    if _a.grid_from:
        GRID = grid_from_coefficients(_a.grid_from)
        print(f"Grid from {_a.grid_from}: {len(GRID)} points")
    print(f"Running grid of {len(GRID)} points with parallel workers...")
    results = {}
    timings = {}
    with cf.ProcessPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(run_one, qt): qt for qt in GRID}
        for fut in cf.as_completed(futs):
            qt, status, path, elapsed = fut.result()
            results[qt] = (status, path)
            timings[qt] = elapsed
            print(f"  ({qt[0]:2d},{qt[1]:2d}) {status}  [{elapsed:.1f}s]")

    failures = [qt for qt, (s, _) in results.items() if not s.startswith("OK") and s != "cached"]
    if failures:
        print(f"FAILURES at: {failures}")
        return 1

    # Build consolidated JSON
    print("\nBuilding consolidated JSON ...")
    samples = []
    skeleton_per_point = {}
    for qt in GRID:
        _, path = results[qt]
        tups = parse_cofactor(path)
        q, t = qt
        samples.append({
            "q": q, "t": t,
            "n_tuples": len(tups),
            "tuples": [
                {"left": l, "idx": i, "right": r, "coef": c}
                for (l, i, r, c) in tups
            ],
        })
        skeleton_per_point[qt] = set((l, i, r) for (l, i, r, _) in tups)

    # Universal skeleton = union of all skeletons present
    all_skels = set()
    for s in skeleton_per_point.values():
        all_skels |= s
    skel_index = sorted(all_skels, key=lambda lir: (lir[1], len(lir[0]), lir[0], len(lir[2]), lir[2]))

    # Per-point cardinality check
    sizes = {qt: len(s) for qt, s in skeleton_per_point.items()}
    full_set = max(skeleton_per_point.values(), key=len)
    full_size = len(full_set)
    mismatches = {qt: full_size - sz for qt, sz in sizes.items() if sz != full_size}

    # Verify all sub-skeletons are subsets of the universal
    bad = []
    for qt, sk in skeleton_per_point.items():
        if not (sk <= full_set):
            bad.append(qt)
    universal_consistent = (not bad)

    out = {
        "n_samples": len(samples),
        "n_skeleton_slots": len(skel_index),
        "skeleton_index": [[l, i, r] for (l, i, r) in skel_index],
        "per_point_sizes": {f"{qt[0]},{qt[1]}": sz for qt, sz in sizes.items()},
        "size_mismatches": {f"{qt[0]},{qt[1]}": missing for qt, missing in mismatches.items()},
        "universal_consistent": universal_consistent,
        "samples": samples,
    }
    out_path = os.path.join(SAMPLES, "coefficients.json")
    with open(out_path, "w") as f:
        json.dump(out, f)
    size = os.path.getsize(out_path)
    print(f"Wrote {out_path}  ({size/1024:.1f} KB)")

    # Summary timing
    total_wall = time.time() - t_start
    sum_individual = sum(timings.values())
    print(f"\nTotal wall time: {total_wall:.1f}s; sum of per-job times: {sum_individual:.1f}s")
    print(f"per_point_sizes (full = {full_size}):")
    for qt in GRID:
        marker = "" if sizes[qt] == full_size else f"  <-- missing {full_size - sizes[qt]}"
        print(f"  ({qt[0]:2d},{qt[1]:2d}): {sizes[qt]}{marker}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
