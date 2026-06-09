#!/usr/bin/env python3
"""Trace-form rank and determinant for a noncommutative quotient, over F_p.

Computes the Gram matrix ``T_{ij} = tr(L_{m_i m_j})`` of the trace form and
its rank. The non-vanishing-det direction of the criterion is
characteristic-free: ``J(A) ⊆ rad T`` always holds, so ``rank T = N``
implies ``A`` is semisimple. The reverse direction (``rank T < N ⇒
non-semisimple``) requires char 0 (Dieudonné); see ``semisimple_Q.py`` for
that path.

The bulk of the work is the ``N × N`` matrix fill. We avoid the naive
quadratic blow-up by:

(1) precomputing each *generator's* left-multiplication matrix once
    (``gen_L[v][j] = reduce(v · m_j)``, sparse), then
(2) computing ``L_{m_l} · e_j`` by walking the letters of ``m_l`` through
    these generator matrices, then
(3) using the structure-constant identity
    ``T_{ij} = sum_l c_l(i,j) · tr(L_{m_l})`` where ``m_i · m_j = sum_l c_l·m_l``.

The final step is Gaussian elimination over F_p.
"""
from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import core


def build_generator_matrices(
    variables: list[str],
    standard: list[tuple[str, ...]],
    mon_idx: dict[tuple[str, ...], int],
    lm_to_tail,
    delta,
    is_match,
    match_pat,
    p: int,
) -> dict[str, dict[int, dict[int, int]]]:
    """``gen_L[v][j] = reduce(v · m_j)`` as a sparse vector."""
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
    return gen_L


def apply_generator(
    gen_L: dict[str, dict[int, dict[int, int]]],
    v: str,
    vec: dict[int, int],
    p: int,
) -> dict[int, int]:
    """Multiply the sparse vector ``vec`` (over the standard basis) on the
    left by the generator ``v``."""
    out: dict[int, int] = {}
    col = gen_L[v]
    for j, c in vec.items():
        cj = col.get(j)
        if cj is None:
            continue
        for k, ck in cj.items():
            new = (out.get(k, 0) + c * ck) % p
            if new == 0 and k in out:
                del out[k]
            else:
                out[k] = new
    return out


def left_mult_by_word(
    gen_L: dict[str, dict[int, dict[int, int]]],
    word: tuple[str, ...],
    vec: dict[int, int],
    p: int,
) -> dict[int, int]:
    """``L_{word}(vec)``: apply ``L_{v_1} L_{v_2} … L_{v_d}`` to ``vec``,
    i.e. walk the generators in reverse and left-multiply each."""
    for v in reversed(word):
        vec = apply_generator(gen_L, v, vec, p)
        if not vec:
            break
    return vec


def trace_form_rank_Fp(
    polys, standard, mon_idx, p, *, log=lambda *a: None,
) -> tuple[int, int, int]:
    """Compute ``(N, rank(T), det(T) mod p)`` for the trace form ``T``."""
    N = len(standard)
    lms = [lm for lm, _ in polys]
    variables = core.collect_variables(polys)
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    lm_to_tail = core.lm_to_tail_Fp(polys, p)

    log("building generator matrices")
    gen_L = build_generator_matrices(
        variables, standard, mon_idx, lm_to_tail,
        delta, is_match, match_pat, p,
    )

    log("computing tr(L_{m_l}) for each standard monomial")
    traces = [0] * N
    for l, m_l in enumerate(standard):
        if not m_l:
            traces[l] = N % p
            continue
        s = 0
        for j in range(N):
            vec = left_mult_by_word(gen_L, m_l, {j: 1}, p)
            s = (s + vec.get(j, 0)) % p
        traces[l] = s

    log(f"filling {N}×{N} trace form")
    T = [[0] * N for _ in range(N)]
    for i, m_i in enumerate(standard):
        for j in range(N):
            vec = left_mult_by_word(gen_L, m_i, {j: 1}, p)
            T[i][j] = sum(c * traces[l] for l, c in vec.items()) % p

    log("Gaussian elimination over F_p")
    det, rank = gauss_rank_det_Fp(T, p)
    return N, rank, det


def gauss_rank_det_Fp(T: list[list[int]], p: int) -> tuple[int, int]:
    """Reduce ``T`` in place; return ``(det mod p, rank)``."""
    N = len(T)
    det = 1
    rank = 0
    for col in range(N):
        pivot = -1
        for r in range(rank, N):
            if T[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            det = 0
            continue
        if pivot != rank:
            T[rank], T[pivot] = T[pivot], T[rank]
            det = (-det) % p
        det = (det * T[rank][col]) % p
        inv = pow(T[rank][col], p - 2, p)
        for r in range(N):
            if r == rank or T[r][col] == 0:
                continue
            fac = (T[r][col] * inv) % p
            for c2 in range(col, N):
                T[r][c2] = (T[r][c2] - fac * T[rank][c2]) % p
        rank += 1
    return det, rank


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis")
    ap.add_argument("--char", type=int, default=2147483647)
    ap.add_argument("--max-dim", type=int, default=2500)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    p = args.char

    def log(msg: str) -> None:
        if args.verbose:
            print(f"[{msg}]", file=sys.stderr, flush=True)

    polys = core.parse_basis_Fp(args.basis)
    variables = core.collect_variables(polys)
    log(f"parsed {len(polys)} polys, {len(variables)} vars")

    lms = [lm for lm, _ in polys]
    delta, is_match, _, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        delta, is_match, variables, max_total=args.max_dim,
    )
    if standard is None:
        print(f"[abort] dim > {args.max_dim}", file=sys.stderr)
        return 1
    N = len(standard)
    log(f"N = {N}")

    mon_idx = {m: i for i, m in enumerate(standard)}
    N, rank, det = trace_form_rank_Fp(polys, standard, mon_idx, p, log=log)
    defect = N - rank
    print(f"[result] N = {N}, rank(T) = {rank}, "
          f"deficit = {defect}, det(T) mod p = {det}")
    if defect == 0:
        print("[conclusion] SEMISIMPLE at this F_p specialization.")
        print("             (Generic semisimplicity over Q(q,t) follows by lift.)")
    else:
        print("[conclusion] NOT SEMISIMPLE at this F_p specialization.")
        print(f"             rad(T) has dim {defect}; J(A) ⊆ rad(T) so dim J ≤ {defect}.")
        print("             Char-0 verdict requires semisimple_Q.py.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
