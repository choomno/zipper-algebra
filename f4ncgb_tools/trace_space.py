#!/usr/bin/env python3
"""
Compute dim A / [A, A] for a quotient algebra given by a Groebner basis.

A trace tau : A -> k with tau(ab) = tau(ba) factors through A / [A, A]
(the abelianization), so the space of traces is the dual to A / [A, A].
This is the cocenter, dim of which equals:
  - the number of conjugacy classes for a group algebra k[G],
  - the number of irreducible modules for a semisimple algebra over alg closed k.

Implementation: build the subspace [A, A] = span{ab - ba : a, b basis}
inside A, then dim A/[A,A] = N - dim [A, A].

In practice: for each pair (i, j) of basis monomials, reduce m_i m_j - m_j m_i
to a sparse vector in the basis, collect all of them, and compute the rank
of the resulting matrix.

The dim of the trace space tells us the maximum possible number of
independent traces on A. For Markov-trace purposes we then need:
  - a tower inclusion Z_{n-1} -> Z_n,
  - one trace per algebra in the tower compatibly with the inclusion,
  - the Markov property: tau_{n+1}(a sigma_n^{+/-1}) = mu^{+/-1} tau_n(a)
    for a in Z_n.
The trace space dimension is the dimension of the "playground" of
candidate traces; the Markov property is a strong linear constraint
that picks out a subspace (or single line) within it.
"""
from __future__ import annotations
import argparse
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from semisimple_Q import (
    parse_basis_Q, collect_variables, build_aho_corasick,
    enumerate_standard_monomials, reduce_to_basis_Q,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="char-0 .basis file from f4ncgb")
    ap.add_argument("--max-dim", type=int, default=500)
    args = ap.parse_args()

    polys = parse_basis_Q(args.basis)
    variables = collect_variables(polys)
    print(f"[parse] {len(polys)} polys over Q, variables = {variables}", flush=True)

    lms = [lm for lm, _ in polys]
    ac_d, ac_m, ac_p = build_aho_corasick(lms, variables)
    standard = enumerate_standard_monomials(ac_d, ac_m, variables, max_total=args.max_dim)
    if standard is None:
        print(f"[abort] dim > {args.max_dim}")
        return 1
    N = len(standard)
    print(f"[dim A] N = {N}", flush=True)

    lm_to_tail = {}
    for lm, terms in polys:
        lead = terms[0][0]
        inv = Fraction(1) / lead
        lm_to_tail[lm] = [(-c * inv, mon) for c, mon in terms[1:]]
    mon_idx = {mon: i for i, mon in enumerate(standard)}

    # For each pair (i, j), compute m_i*m_j - m_j*m_i as a vector in basis.
    # Collect all these vectors as ROWS of an Nxsomething matrix.
    print(f"[commutators] computing m_i*m_j - m_j*m_i for all pairs (i, j)...", flush=True)
    rows = []
    for i in range(N):
        for j in range(i + 1, N):
            v1 = reduce_to_basis_Q([(Fraction(1), standard[i] + standard[j])],
                                   mon_idx, lm_to_tail, ac_d, ac_m, ac_p)
            v2 = reduce_to_basis_Q([(Fraction(1), standard[j] + standard[i])],
                                   mon_idx, lm_to_tail, ac_d, ac_m, ac_p)
            # diff = v1 - v2
            diff = dict(v1)
            for k, c in v2.items():
                new = diff.get(k, Fraction(0)) - c
                if new == 0:
                    if k in diff:
                        del diff[k]
                else:
                    diff[k] = new
            if diff:
                row = [Fraction(0)] * N
                for k, c in diff.items():
                    row[k] = c
                rows.append(row)
        if (i + 1) % 25 == 0:
            print(f"   row {i+1}/{N} done", file=sys.stderr, flush=True)
    print(f"[commutators] collected {len(rows)} commutator rows; reducing...", flush=True)

    # rank of the rows = dim [A, A] inside A
    # Gaussian elimination over Q
    M = [row[:] for row in rows]
    rank = 0
    for col in range(N):
        pivot = -1
        for r in range(rank, len(M)):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]
        inv = Fraction(1) / M[rank][col]
        for r in range(len(M)):
            if r != rank and M[r][col] != 0:
                fac = M[r][col] * inv
                for c2 in range(col, N):
                    M[r][c2] -= fac * M[rank][c2]
        rank += 1
    dim_comm = rank
    dim_trace = N - dim_comm
    print()
    print(f"[result] dim A         = {N}")
    print(f"[result] dim [A, A]    = {dim_comm}")
    print(f"[result] dim A/[A,A]   = {dim_trace}")
    print()
    print(f"[interpretation] The space of trace functionals tau : A -> Q")
    print(f"                 with tau(ab) = tau(ba) has dimension {dim_trace}.")
    print(f"                 For a semisimple algebra over alg closed k, this")
    print(f"                 equals the number of irreducible summands.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
