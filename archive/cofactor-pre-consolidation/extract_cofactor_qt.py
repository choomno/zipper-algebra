#!/usr/bin/env python3
"""Run extract_cofactor pipeline at a custom (q, t) over Q.

Generates basis + non-expanded proof via f4ncgb, then performs the same
reduction + composition as tools/extract_cofactor.py but
parameterised over (q, t).
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py lives
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import argparse
import subprocess
import time

REPO = os.path.dirname(PARENT)   # tools/ -> zipper-algebra/
TOOLS = os.path.join(REPO, "tools")
LOGS = os.path.join(REPO, "logs_and_traces", "f4ncgb")
F4 = "/home/choom/src/f4ncgb/build-release/f4ncgb"

sys.path.insert(0, TOOLS)
from extract_cofactor import (  # noqa: E402
    parse_basis, parse_noex_proof, expand_proofs,
    reduce_against_basis, compose_cofactor,
    verify_cofactor, write_cofactor, VAR_NAMES, _mon_str, _frac_str,
)
from zipper_in_ms import build_relations, Coef, Poly  # noqa: E402
def write_cofactor_qt(path, terms, originals, q_int, t_int):
    terms = sorted(terms, key=lambda t: (t[2], len(t[1]), t[1], len(t[3]), t[3]))
    used = sorted({t[2] for t in terms})
    with open(path, "w") as f:
        f.write(f"# Cofactor of X_3 at (q,t)=({q_int},{t_int}) over Q.\n")
        f.write(f"# Number of terms: {len(terms)}\n")
        f.write(f"# Distinct inputs used: {used}\n")
        f.write("# --- Original input relations ---\n")
        for i, g in enumerate(originals):
            f.write(f"# [{i}] = {g.to_msolve_string(VAR_NAMES)}\n")
        f.write("# --- End ---\n")
        for coef, u, idx, v in terms:
            cs = _frac_str(coef)
            left = _mon_str(u); right = _mon_str(v)
            pieces = []
            pieces.append(("+" + cs) if coef > 0 else cs)
            if left: pieces.append(left)
            pieces.append(f"[{idx}]")
            if right: pieces.append(right)
            f.write(" * ".join(pieces) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, required=True)
    ap.add_argument("--t", type=int, required=True)
    args = ap.parse_args()
    q_int, t_int = args.q, args.t
    char = 0
    n = 4

    ms_path = os.path.join(REPO, "ms_inputs", "qt_sweep", f"zipper-n4-Q-q{q_int}-t{t_int}.ms")
    basis_path = os.path.join(LOGS, f"basis-n4-Q-q{q_int}-t{t_int}.txt")
    proof_path = os.path.join(LOGS, f"proof-n4-Q-q{q_int}-t{t_int}-noex.txt")
    cofactor_path = os.path.join(LOGS, f"cofactor-X3-n4-Q-q{q_int}-t{t_int}.txt")

    # 1) Generate .ms file (no --target, we want 31 inputs)
    print(f"=== (q,t)=({q_int},{t_int}) ===")
    t0 = time.time()
    # --untwist passed EXPLICITLY: L97 builds the Python-side relations with
    # "bigelow", and relying on the emitter's default to match would let a change
    # of default silently compute the cofactor against a different ideal than the
    # basis f4ncgb produced.
    subprocess.run([sys.executable, os.path.join(TOOLS, "zipper_in_ms.py"),
                    "--n", str(n), "--q", str(q_int), "--t", str(t_int),
                    "--char", "0", "--untwist", "bigelow",
                    "--out", ms_path], check=True)
    print(f"ms gen: {time.time()-t0:.2f}s")

    # 2) Run f4ncgb with -p (proof on, NOT expanded) and -o for basis
    t0 = time.time()
    res = subprocess.run([F4, "-i", ms_path, "-o", basis_path,
                          "-p", proof_path, "-e", "off"],
                         capture_output=True, text=True)
    print(f"f4ncgb: {time.time()-t0:.2f}s; rc={res.returncode}")
    if res.returncode != 0:
        print("STDOUT:", res.stdout[-500:])
        print("STDERR:", res.stderr[-500:])
        return 1

    # 3) Extract cofactor
    basis = parse_basis(basis_path, char)
    noex = parse_noex_proof(proof_path)
    print(f"basis polys: {len(basis)}; proof lines: {len(noex)}")
    assert len(basis) == len(noex)

    expanded = expand_proofs(noex)

    var_names, relations, S = build_relations(n, q_int, t_int, char, "bigelow")
    assert len(relations) == 31
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
    print(f"expanded sanity: {bad} bad / {len(basis)}")
    if bad:
        return 1

    remainder, steps = reduce_against_basis(x3, basis, char)
    print(f"reduction steps: {len(steps)}; remainder terms: {len(remainder.terms)}")
    if remainder.terms:
        return 1

    cofactor_terms = compose_cofactor(steps, expanded, char)
    ok = verify_cofactor(cofactor_terms, relations, x3, char)
    print(f"verify: {'OK' if ok else 'FAIL'}; tuples: {len(cofactor_terms)}")
    if not ok:
        return 1

    write_cofactor_qt(cofactor_path, cofactor_terms, relations, q_int, t_int)
    print(f"wrote {cofactor_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
