#!/usr/bin/env python3
r"""Dimension of the cocenter ``A / [A, A]`` over Q.

A trace functional ``τ: A → Q`` satisfies ``τ(ab) = τ(ba)``, equivalently
vanishes on the subspace ``[A, A] := span{ab − ba}``. So the space of
traces is dual to ``A / [A, A]``, and ``dim(A/[A,A])`` is the number of
linearly independent traces. For semisimple ``A`` over an algebraically
closed field, this equals the number of irreducible summands in the
Wedderburn decomposition.

Algorithm: build ``[A, A]`` row-by-row (commutators of basis-monomial
pairs reduced to standard form), then ``dim A / [A, A] = N − rank``. The
exact-rational arithmetic makes this slow at large ``N``; use
:mod:`trace_space_fast` for the F_p analogue.


Experimental interpretation: dimension, minimality and quotient-structure
conclusions require a certified complete basis or a separate exact witness.
A successful solver exit and agreement at several specializations do not supply
that evidence. State the corresponding evidence when interpreting these observations.
"""
from __future__ import annotations

import argparse
from fractions import Fraction


from zipper_algebra import core, Presentation
from zipper_algebra.linear import matrix_rank


def build_commutator_matrix_Q(
    standard, mon_idx, lm_to_tail, delta, is_match, match_pat, N,
) -> list[list[Fraction]]:
    """Each row is ``reduce(m_i m_j − m_j m_i)`` in the standard basis."""
    rows: list[list[Fraction]] = []
    for i in range(N):
        for j in range(i + 1, N):
            a = core.reduce_to_basis_Q(
                [(Fraction(1), standard[i] + standard[j])],
                mon_idx, lm_to_tail, delta, is_match, match_pat,
            )
            b = core.reduce_to_basis_Q(
                [(Fraction(1), standard[j] + standard[i])],
                mon_idx, lm_to_tail, delta, is_match, match_pat,
            )
            diff = dict(a)
            for k, c in b.items():
                new = diff.get(k, Fraction(0)) - c
                if new == 0 and k in diff:
                    del diff[k]
                else:
                    diff[k] = new
            if diff:
                row = [Fraction(0)] * N
                for k, c in diff.items():
                    row[k] = c
                rows.append(row)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="char-0 .basis file from f4ncgb")
    ap.add_argument("--max-dim", type=int, default=500)
    ap.add_argument("--n", type=int, required=True, help="declared strand count for the full generator alphabet")
    args = ap.parse_args()

    polys = core.parse_basis_Q(args.basis)
    variables = Presentation(args.n).var_names()
    print(f"[parse] {len(polys)} polys over Q, variables = {variables}",
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
    print(f"[dim A] N = {N}", flush=True)

    mon_idx = {m: i for i, m in enumerate(standard)}
    lm_to_tail = core.lm_to_tail_Q(polys)

    print("[commutators] building [A, A] matrix...", flush=True)
    rows = build_commutator_matrix_Q(
        standard, mon_idx, lm_to_tail, delta, is_match, match_pat, N,
    )
    print(f"[commutators] {len(rows)} non-zero commutator rows", flush=True)

    print("[rank] Gaussian elimination over Q...", flush=True)
    dim_comm = matrix_rank(rows)
    dim_trace = N - dim_comm

    print()
    print(f"[result] dim A         = {N}")
    print(f"[result] dim [A, A]    = {dim_comm}")
    print(f"[result] dim A/[A,A]   = {dim_trace}")
    print()
    print(f"[interpretation] The space of trace functionals τ: A → Q has")
    print(f"                 dimension {dim_trace}. For a semisimple algebra")
    print(f"                 over an algebraically closed field this equals")
    print(f"                 the number of irreducible summands.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
