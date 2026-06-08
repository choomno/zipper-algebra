#!/usr/bin/env python3
"""
Fast semisimplicity check via the trace form, using only generator
left-multiplication matrices.

  STRATEGY:
  1. Build generator left-mult: for each generator v in {x_1, y_1, ...},
     gen_L[v] is a dict j -> {k: coef} encoding column j of L_v.
  2. To compute L_{m_l} as a matrix (sparse): for each column j, walk the
     basis vector e_j through the generator chain m_l = v_1...v_d, applying
     L_{v_d}, ..., L_{v_1} in turn.
  3. Trace of L_{m_l} = sum over j of (column-j-result)[j].
  4. For each (i, j) compute reduce(m_i * m_j) by walking standard[i] through
     the same chain applied to e_j. Then trace form T[i, j] = sum_k coef * t[k].

Total cost:
  - Build gen_L: |gens| * N reductions of short words.
  - Compute trace[l] for all l: N * N * avg_len * avg_sparsity = O(N^2 * d * s).
  - Compute T: O(N^2 * d * s).

The constant factor is much smaller than the from-scratch reduce-to-basis
approach because every step is a sparse matrix-vector multiply on a small
domain (avg_sparsity << N).
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from collections import deque

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from semisimple_check import (
    parse_basis_full, collect_variables, build_aho_corasick,
    enumerate_standard_monomials, find_leading_match,
)


def reduce_to_basis(word_terms, mon_idx, lm_to_tail, ac_delta, ac_is_match,
                    ac_match_pat, p):
    out = {}
    worklist = list(word_terms)
    while worklist:
        coef, w = worklist.pop()
        if coef == 0:
            continue
        if w in mon_idx:
            k = mon_idx[w]
            new = (out.get(k, 0) + coef) % p
            if new == 0 and k in out:
                del out[k]
            else:
                out[k] = new
            continue
        m = find_leading_match(w, ac_delta, ac_is_match, ac_match_pat)
        if m is None:
            raise RuntimeError(f"Standard word {w} not enumerated")
        start, end, lm = m
        u, v = w[:start], w[end:]
        for c, t_mon in lm_to_tail[lm]:
            worklist.append(((coef * c) % p, u + t_mon + v))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis")
    ap.add_argument("--char", type=int, default=2147483647)
    ap.add_argument("--max-dim", type=int, default=2500)
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    p = args.char

    polys = parse_basis_full(args.basis)
    variables = collect_variables(polys)
    print(f"[parse] {len(polys)} polys, {len(variables)} vars: {variables}",
          file=sys.stderr, flush=True)

    lms = [lm for lm, _ in polys]
    ac_delta, ac_is_match, ac_match_pat, _ = build_aho_corasick(lms, variables)
    standard = enumerate_standard_monomials(ac_delta, ac_is_match, None,
                                            variables, max_total=args.max_dim)
    if standard is None:
        print(f"[abort] dim > {args.max_dim}", file=sys.stderr)
        return 1
    N = len(standard)
    print(f"[dim] N = {N}", file=sys.stderr, flush=True)

    lm_to_tail = {}
    for lm, terms in polys:
        lead = terms[0][0]
        inv = pow(lead, p - 2, p)
        lm_to_tail[lm] = [((-c * inv) % p, mon) for c, mon in terms[1:]]
    mon_idx = {mon: i for i, mon in enumerate(standard)}

    # Build gen_L: for each generator v, gen_L[v][j] = dict {k: coef}
    # representing reduce(v * m_j) in standard basis.
    print(f"[gen-L] building L_v for {len(variables)} generators...",
          file=sys.stderr, flush=True)
    gen_L = {}
    for v in variables:
        col = {}
        for j, m_j in enumerate(standard):
            red = reduce_to_basis([(1, (v,) + m_j)],
                                  mon_idx, lm_to_tail,
                                  ac_delta, ac_is_match, ac_match_pat, p)
            if red:
                col[j] = red
        gen_L[v] = col
    print(f"[gen-L] done.", file=sys.stderr, flush=True)

    def apply_gen(v, vec):
        """Apply L_v to a sparse vector vec = dict {j: coef}; returns
        new sparse vector dict."""
        out = {}
        col = gen_L[v]
        for j, c in vec.items():
            if j not in col:
                continue
            for k, ck in col[j].items():
                new = (out.get(k, 0) + c * ck) % p
                if new == 0 and k in out:
                    del out[k]
                else:
                    out[k] = new
        return out

    # For each standard monomial m_l = v_1...v_d, compute its trace t[l] =
    # trace(L_{m_l}) by applying L_{v_d}, ..., L_{v_1} to each basis vector e_j
    # and reading the j-th coordinate.
    print(f"[traces] computing tr(L_{{m_l}}) for l=0..{N-1}...",
          file=sys.stderr, flush=True)
    traces = [0] * N
    for l, m_l in enumerate(standard):
        if not m_l:
            # m_l is the empty monomial = 1; L_1 = identity; trace = N
            traces[l] = N % p
            continue
        s = 0
        for j in range(N):
            vec = {j: 1}
            # apply L_{v_d}, ..., L_{v_1} to e_j: actually L_{m_l} = L_{v_1} L_{v_2} ... L_{v_d}
            # acting on m_j gives v_1*v_2*...*v_d*m_j; we apply from RIGHT to LEFT in word.
            # If word = v_1 v_2 ... v_d (left to right), then L_{m_l} = L_{v_1} after L_{v_2} after ...
            # We want L_{m_l} e_j = L_{v_1}(L_{v_2}(...L_{v_d}(e_j))).
            # So we apply v_d, then v_{d-1}, ..., then v_1 — i.e. reversed.
            for v in reversed(m_l):
                vec = apply_gen(v, vec)
                if not vec:
                    break
            s = (s + vec.get(j, 0)) % p
        traces[l] = s
        if args.verbose and (l + 1) % 25 == 0:
            print(f"   t[{l+1}/{N}]", file=sys.stderr, flush=True)
    print(f"[traces] done.", file=sys.stderr, flush=True)

    # Now build trace form: T[i, j] = sum_l (m_i * m_j reduced coef of m_l) * traces[l]
    # We compute m_i * m_j by walking e_j through letters of m_i (reversed).
    print(f"[T] filling N x N trace form...", file=sys.stderr, flush=True)
    T = [[0] * N for _ in range(N)]
    for i in range(N):
        m_i = standard[i]
        for j in range(N):
            vec = {j: 1}
            for v in reversed(m_i):
                vec = apply_gen(v, vec)
                if not vec:
                    break
            # vec is now reduce(m_i * m_j) as a sparse dict {l: coef}
            s = 0
            for l, c in vec.items():
                s = (s + c * traces[l]) % p
            T[i][j] = s
        if args.verbose and (i + 1) % 25 == 0:
            print(f"   T row {i+1}/{N}", file=sys.stderr, flush=True)
    print(f"[T] done.", file=sys.stderr, flush=True)

    # Compute rank (and det) via Gaussian elimination over GF(p)
    print(f"[rank] computing rank(T) by Gaussian elimination...",
          file=sys.stderr, flush=True)
    det = 1
    rank = 0
    M = [row[:] for row in T]
    for col in range(N):
        pivot = -1
        for r in range(rank, N):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            det = 0
            continue
        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]
            det = (-det) % p
        det = (det * M[rank][col]) % p
        inv = pow(M[rank][col], p - 2, p)
        for r in range(N):
            if r == rank or M[r][col] == 0:
                continue
            fac = (M[r][col] * inv) % p
            for c2 in range(col, N):
                M[r][c2] = (M[r][c2] - fac * M[rank][c2]) % p
        rank += 1
    defect = N - rank
    print(f"[result] N = {N}, rank(T) = {rank}, deficit = {defect}, det(T) mod p = {det}")
    if defect == 0:
        print(f"[conclusion] SEMISIMPLE at (q,t)=specialization in F_p.")
        print(f"             (By the lift argument: also generically semisimple over Q(q,t).)")
    else:
        print(f"[conclusion] NOT SEMISIMPLE at this specialization in F_p.")
        print(f"             Radical of trace form has dim {defect};")
        print(f"             J(A) is contained in it but may be smaller.")
        print(f"             Generic semisimplicity over Q(q,t) is NOT ruled out;")
        print(f"             retry at independent (q,t) to distinguish.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
