#!/usr/bin/env python3
"""Re-run the Python cofactor extraction pipeline at all 212 cached (q,t)
sample points using the canonicalized extract_cofactor.py, write outputs to
logs_and_traces/f4ncgb/lift_samples_canonical/, and consolidate to
coefficients.json.

This script does NOT regenerate basis/proof via f4ncgb -- it reads the
already-cached basis-n4-Q-q{Q}-t{T}.txt and proof-n4-Q-q{Q}-t{T}-noex.txt
files. f4ncgb's reduction order is the documented upstream source of the 73
stubborn slots; the Python side is verified deterministic. The
canonicalisation patch in extract_cofactor.py (explicit lex+index tiebreak
on lt_basis) is a no-op for the n=4 basis (no LM collisions across all 137
polys) but is applied here for documentation/future-proofing.
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py lives
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import concurrent.futures as cf
import pathlib as _pl
sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[3]))
import json
import time
REPO = "/home/choom/braid-research/zipper-algebra"
TOOLS = os.path.join(REPO, "tools")
LOGS = os.path.join(REPO, "logs_and_traces", "f4ncgb")
SAMPLES_OLD = os.path.join(LOGS, "lift_samples")
SAMPLES = os.path.join(LOGS, "lift_samples_canonical")
os.makedirs(SAMPLES, exist_ok=True)

sys.path.insert(0, TOOLS)


def run_one(qt):
    q, t = qt
    # Import inside worker (process pool) so each child has clean state.
    from extract_cofactor import (
        parse_basis, parse_noex_proof, expand_proofs,
        reduce_against_basis, compose_cofactor,
        verify_cofactor, VAR_NAMES, _mon_str, _frac_str,
    )
    from zipper_in_ms import build_relations, Poly, Coef

    char = 0
    n = 4

    basis_path = os.path.join(LOGS, f"basis-n4-Q-q{q}-t{t}.txt")
    proof_path = os.path.join(LOGS, f"proof-n4-Q-q{q}-t{t}-noex.txt")
    out_path = os.path.join(SAMPLES, f"cofactor-X3-n4-Q-q{q}-t{t}.txt")

    if not (os.path.exists(basis_path) and os.path.exists(proof_path)):
        return qt, f"FAIL: missing cached basis/proof for ({q},{t})", None, 0.0

    # Skip if already done.
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return qt, "CACHED", out_path, 0.0

    t0 = time.time()
    basis = parse_basis(basis_path, char)
    noex = parse_noex_proof(proof_path)
    if len(basis) != len(noex):
        return qt, f"FAIL: len mismatch basis={len(basis)} noex={len(noex)}", None, time.time() - t0

    expanded = expand_proofs(noex)

    var_names, relations, S = build_relations(n, q, t, char, "bigelow")
    if len(relations) != 31:
        return qt, f"FAIL: relations count {len(relations)}", None, time.time() - t0
    x3 = S[(3, 1)]

    # Sanity-check expanded reconstructs basis
    bad = 0
    for k, exp in enumerate(expanded):
        acc = Poly(char)
        for c, u, idx, v in exp:
            lp = Poly(char); lp.add_term(u, Coef.one(char))
            rp = Poly(char); rp.add_term(v, Coef.one(char))
            prod = (lp * relations[idx] * rp).scalar_mul(Coef(c, char))
            acc = acc + prod
        if (acc - basis[k]).terms:
            bad += 1
    if bad:
        return qt, f"FAIL: {bad} bad expansions", None, time.time() - t0

    remainder, steps = reduce_against_basis(x3, basis, char)
    if remainder.terms:
        return qt, "FAIL: nonzero remainder", None, time.time() - t0

    cofactor_terms = compose_cofactor(steps, expanded, char)
    ok = verify_cofactor(cofactor_terms, relations, x3, char)
    if not ok:
        return qt, "FAIL: verify failed", None, time.time() - t0

    # Write cofactor file
    terms = sorted(cofactor_terms,
                   key=lambda tt: (tt[2], len(tt[1]), tt[1], len(tt[3]), tt[3]))
    used = sorted({tt[2] for tt in terms})
    with open(out_path, "w") as f:
        f.write(f"# Cofactor of X_3 at (q,t)=({q},{t}) over Q (canonicalised).\n")
        f.write(f"# Number of terms: {len(terms)}\n")
        f.write(f"# Distinct inputs used: {used}\n")
        f.write("# --- Original input relations ---\n")
        for i, g in enumerate(relations):
            f.write(f"# [{i}] = {g.to_msolve_string(VAR_NAMES)}\n")
        f.write("# --- End ---\n")
        for coef, u, idx, v in terms:
            cs = _frac_str(coef)
            left = _mon_str(u); right = _mon_str(v)
            pieces = [("+" + cs) if coef > 0 else cs]
            if left:
                pieces.append(left)
            pieces.append(f"[{idx}]")
            if right:
                pieces.append(right)
            f.write(" * ".join(pieces) + "\n")

    return qt, f"OK {len(cofactor_terms)} tuples ({len(steps)} steps)", out_path, time.time() - t0


def parse_cofactor(path):
    tuples = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split("*")]
            coef = parts[0]
            left = ""; right = ""; idx = None
            for i, p in enumerate(parts[1:], 1):
                if p.startswith("[") and p.endswith("]"):
                    idx = int(p[1:-1])
                    if i > 1:
                        left = "*".join(parts[1:i])
                    if i < len(parts) - 1:
                        right = "*".join(parts[i + 1:])
                    break
            assert idx is not None, f"unparsed: {line}"
            if coef.startswith("+"):
                coef = coef[1:]
            tuples.append((left, idx, right, coef))
    return tuples


def main():
    # Discover all 212 sample points from existing lift_samples coefficients.json
    with open(os.path.join(SAMPLES_OLD, "coefficients.json")) as f:
        old = json.load(f)
    GRID = sorted({(s["q"], s["t"]) for s in old["samples"]})
    print(f"Discovered {len(GRID)} sample points from lift_samples/coefficients.json")

    t_start = time.time()
    results = {}
    timings = {}
    # parallelism is a machine fact -> [limits].jobs in settings.toml
    try:
        from repo_settings import settings as _settings
        _jobs = int(_settings().limit('jobs'))
    except Exception:
        _jobs = 6
    with cf.ProcessPoolExecutor(max_workers=_jobs) as ex:
        futs = {ex.submit(run_one, qt): qt for qt in GRID}
        done = 0
        for fut in cf.as_completed(futs):
            qt, status, path, elapsed = fut.result()
            results[qt] = (status, path)
            timings[qt] = elapsed
            done += 1
            if done % 20 == 0 or not status.startswith("OK"):
                print(f"  [{done:3d}/{len(GRID)}] ({qt[0]:2d},{qt[1]:2d}) {status}  [{elapsed:.1f}s]")

    failures = [qt for qt, (s, _) in results.items()
                if not (s.startswith("OK") or s == "CACHED")]
    if failures:
        print(f"FAILURES at {len(failures)} points: {failures[:10]}")
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
            "tuples": [{"left": l, "idx": i, "right": r, "coef": c}
                       for (l, i, r, c) in tups],
        })
        skeleton_per_point[qt] = set((l, i, r) for (l, i, r, _) in tups)

    all_skels = set()
    for s in skeleton_per_point.values():
        all_skels |= s
    skel_index = sorted(all_skels,
                        key=lambda lir: (lir[1], len(lir[0]), lir[0],
                                         len(lir[2]), lir[2]))
    sizes = {qt: len(s) for qt, s in skeleton_per_point.items()}
    full_size = max(sizes.values())
    mismatches = {qt: full_size - sz for qt, sz in sizes.items() if sz != full_size}
    bad = [qt for qt, sk in skeleton_per_point.items() if not (sk <= all_skels)]
    universal_consistent = (not bad)

    out = {
        "n_samples": len(samples),
        "n_skeleton_slots": len(skel_index),
        "skeleton_index": [[l, i, r] for (l, i, r) in skel_index],
        "per_point_sizes": {f"{qt[0]},{qt[1]}": sz for qt, sz in sizes.items()},
        "size_mismatches": {f"{qt[0]},{qt[1]}": missing
                            for qt, missing in mismatches.items()},
        "universal_consistent": universal_consistent,
        "samples": samples,
    }
    out_path = os.path.join(SAMPLES, "coefficients.json")
    with open(out_path, "w") as f:
        json.dump(out, f)
    size = os.path.getsize(out_path)
    print(f"Wrote {out_path}  ({size/1024:.1f} KB)")

    total_wall = time.time() - t_start
    sum_individual = sum(timings.values())
    print(f"\nTotal wall time: {total_wall:.1f}s; sum of per-job times: {sum_individual:.1f}s")
    print(f"per_point_sizes full = {full_size}; n_skeleton_slots = {len(skel_index)}")
    print(f"size mismatches: {len(mismatches)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
