#!/usr/bin/env python3
r"""Dimension of the cocenter ``A / [A, A]`` over F_p, via online row reduction.

Same target as :mod:`trace_space` but over F_p, with the post-hoc cubic
Gaussian-elimination step replaced by an online reduction interleaved with
commutator generation. This makes large ``N`` (e.g. ``N = 186`` for our
``n = 4`` variants) tractable in pure Python.

Strategy:

1. Build the AC automaton, enumerate standard monomials, build
   ``gen_L[v][j] = reduce(v · m_j)`` for each generator.
2. For each ``i``, precompute ``L_{m_i}`` as a sparse column map. Then
   ``m_i · m_j`` reduced is just ``L_{m_i}[j]``.
3. For each pair ``(i, j)`` with ``i < j``, form
   ``v_{ij} = L_{m_i}[j] − L_{m_j}[i]`` and reduce it against the running
   echelon form ``E`` of ``[A, A]``; if a non-zero residue remains, add
   it as a new pivot.
4. ``dim A / [A, A] = N − |E|``.

The post-mortem rank is just ``|E|``; no separate Gaussian-elimination
pass is needed.


Experimental interpretation: dimension, minimality and quotient-structure
conclusions require a certified complete basis or a separate exact witness.
A successful solver exit and agreement at several specializations do not supply
that evidence. State the corresponding evidence when interpreting these observations.
"""
from __future__ import annotations

import argparse
import os
import sys
import time


from zipper_algebra import core, Presentation
from zipper_algebra.linear import validate_field_characteristic


def precompute_left_mult_columns(
    standard, gen_L, p: int,
) -> list[dict[int, dict[int, int]]]:
    """``mi_action[i][j] = reduce(m_i · m_j)`` as a sparse vector."""
    N = len(standard)
    mi_action: list[dict[int, dict[int, int]]] = []
    for m_i in standard:
        col_map: dict[int, dict[int, int]] = {}
        if not m_i:
            for j in range(N):
                col_map[j] = {j: 1}
        else:
            for j in range(N):
                vec = {j: 1}
                for v in reversed(m_i):
                    nxt: dict[int, int] = {}
                    col = gen_L[v]
                    for jj, c in vec.items():
                        cj = col.get(jj)
                        if cj is None:
                            continue
                        for k, ck in cj.items():
                            new = (nxt.get(k, 0) + c * ck) % p
                            if new == 0 and k in nxt:
                                del nxt[k]
                            else:
                                nxt[k] = new
                    vec = nxt
                    if not vec:
                        break
                if vec:
                    col_map[j] = vec
        mi_action.append(col_map)
    return mi_action


def cocenter_rank_Fp(
    mi_action: list[dict[int, dict[int, int]]], N: int, p: int,
) -> int:
    """Online row reduction of ``{m_i m_j − m_j m_i}_{i<j}`` modulo p.

    ``echelon[pivot] = {idx: coef}`` is a sparse row with the pivot
    coefficient normalised to 1. The pivot is the smallest non-zero index
    in the row.
    """
    echelon: dict[int, dict[int, int]] = {}
    for i in range(N):
        for j in range(i + 1, N):
            v1 = mi_action[i].get(j, {})
            v2 = mi_action[j].get(i, {})
            diff: dict[int, int] = dict(v1)
            for k, c in v2.items():
                new = (diff.get(k, 0) - c) % p
                if new == 0 and k in diff:
                    del diff[k]
                else:
                    diff[k] = new
            if not diff:
                continue
            while True:
                pivot_in_diff = None
                for idx in diff:
                    if idx in echelon and (
                        pivot_in_diff is None or idx < pivot_in_diff
                    ):
                        pivot_in_diff = idx
                if pivot_in_diff is None:
                    break
                c = diff[pivot_in_diff]
                row = echelon[pivot_in_diff]
                for k, rc in row.items():
                    new = (diff.get(k, 0) - c * rc) % p
                    if new == 0 and k in diff:
                        del diff[k]
                    else:
                        diff[k] = new
            if not diff:
                continue
            piv = min(diff.keys())
            piv_inv = pow(diff[piv], p - 2, p)
            echelon[piv] = {k: (c * piv_inv) % p for k, c in diff.items()}
    return len(echelon)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="f4ncgb F_p .basis file")
    ap.add_argument("--char", type=int, default=2147483647)
    ap.add_argument("--max-dim", type=int, default=2500)
    ap.add_argument("--n", type=int, required=True, help="declared strand count for the full generator alphabet")
    args = ap.parse_args()
    p = validate_field_characteristic(args.char)
    if not p:
        ap.error("this recipe requires a positive prime characteristic")

    t0 = time.time()
    polys = core.parse_basis_Fp(args.basis)
    variables = Presentation(args.n).var_names()
    print(f"[parse] {len(polys)} polys, {len(variables)} vars: {variables}",
          flush=True)

    lms = [lm for lm, _ in polys]
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        delta, is_match, variables, max_total=args.max_dim,
    )
    if standard is None:
        print(f"[abort] dim > {args.max_dim}")
        return 1
    N = len(standard)
    print(f"[dim] N = observed standard-word count = {N}", flush=True)

    mon_idx = {m: i for i, m in enumerate(standard)}
    lm_to_tail = core.lm_to_tail_Fp(polys, p)

    print(f"[gen-L] building L_v for {len(variables)} generators...", flush=True)
    gen_L: dict[str, dict[int, dict[int, int]]] = {}
    for v in variables:
        col: dict[int, dict[int, int]] = {}
        for j, m_j in enumerate(standard):
            red = core.reduce_to_basis_Fp(
                [(1, (v,) + m_j)],
                mon_idx, lm_to_tail, delta, is_match, match_pat, p,
            )
            if red:
                col[j] = red
        gen_L[v] = col
    print(f"[gen-L] done in {time.time()-t0:.1f}s.", flush=True)

    print(f"[L-action] precomputing L_{{m_i}} columns...", flush=True)
    t_la = time.time()
    mi_action = precompute_left_mult_columns(standard, gen_L, p)
    print(f"[L-action] done in {time.time()-t_la:.1f}s.", flush=True)

    print(f"[echelon] online reduction of {N*(N-1)//2} commutator vectors...",
          flush=True)
    t_ech = time.time()
    rank_comm = cocenter_rank_Fp(mi_action, N, p)
    dim_trace = N - rank_comm
    print(f"[echelon] done in {time.time()-t_ech:.1f}s.", flush=True)

    total = time.time() - t0
    print()
    print(f"[result] dim A         = {N}")
    print(f"[result] dim [A, A]    = {rank_comm}")
    print(f"[result] dim A/[A,A]   = {dim_trace}")
    print(f"[result] total time    = {total:.1f}s")
    print()
    print(f"[interpretation] Trace space dimension = {dim_trace}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
