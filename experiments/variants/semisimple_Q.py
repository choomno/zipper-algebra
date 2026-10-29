#!/usr/bin/env python3
r"""Trace-form rank and determinant for a noncommutative quotient, over Q.

In characteristic zero the trace-form criterion (Pierce, *Associative
Algebras*, §9.4) is an iff: ``A`` is semisimple iff ``det T ≠ 0``, and
when ``det T = 0`` the Jacobson radical has dimension ``N − rank T``
exactly. This script runs the criterion over the rationals (using exact
:class:`fractions.Fraction` arithmetic).

The algorithm is the structure-constant decomposition of
:mod:`semisimple_v3`, repeated over ``Q`` -- slower than F_p because of
fraction reduction in inner loops, but it produces a theorem at the
``(q, t)`` specialisation rather than the lift-only statement available in
``F_p``.


Experimental interpretation: dimension, minimality and quotient-structure
conclusions require a certified complete basis or a separate exact witness.
A successful solver exit and agreement at several specializations do not supply
that evidence. State the corresponding evidence when interpreting these observations.
"""
from __future__ import annotations

import argparse
import sys
from fractions import Fraction


from zipper_algebra import core, Presentation


def build_generator_matrices_Q(
    variables, standard, mon_idx, lm_to_tail, delta, is_match, match_pat,
) -> dict[str, dict[int, dict[int, Fraction]]]:
    gen_L: dict[str, dict[int, dict[int, Fraction]]] = {}
    for v in variables:
        col: dict[int, dict[int, Fraction]] = {}
        for j, m_j in enumerate(standard):
            red = core.reduce_to_basis_Q(
                [(Fraction(1), (v,) + m_j)],
                mon_idx, lm_to_tail, delta, is_match, match_pat,
            )
            if red:
                col[j] = red
        gen_L[v] = col
    return gen_L


def apply_generator_Q(gen_L, v, vec):
    out: dict[int, Fraction] = {}
    col = gen_L[v]
    for j, c in vec.items():
        cj = col.get(j)
        if cj is None:
            continue
        for k, ck in cj.items():
            new = out.get(k, Fraction(0)) + c * ck
            if new == 0 and k in out:
                del out[k]
            else:
                out[k] = new
    return out


def left_mult_by_word_Q(gen_L, word, vec):
    for v in reversed(word):
        vec = apply_generator_Q(gen_L, v, vec)
        if not vec:
            break
    return vec


def gauss_rank_det_Q(T: list[list[Fraction]]) -> tuple[Fraction, int]:
    """Reduce ``T`` in place; return ``(det, rank)`` over Q."""
    N = len(T)
    det = Fraction(1)
    rank = 0
    sign = 1
    for col in range(N):
        pivot = -1
        for r in range(rank, N):
            if T[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            det = Fraction(0)
            continue
        if pivot != rank:
            T[rank], T[pivot] = T[pivot], T[rank]
            sign = -sign
        det *= T[rank][col]
        inv = Fraction(1) / T[rank][col]
        for r in range(N):
            if r == rank or T[r][col] == 0:
                continue
            fac = T[r][col] * inv
            for c2 in range(col, N):
                T[r][c2] -= fac * T[rank][c2]
        rank += 1
    return (det * sign if det != 0 else Fraction(0)), rank


def trace_form_rank_Q(
    polys, standard, mon_idx, *, log=lambda *a: None, variables=None,
) -> tuple[int, int, Fraction]:
    N = len(standard)
    lms = [lm for lm, _ in polys]
    if variables is None:
        variables = core.collect_variables(polys)
        log("observed-alphabet reduction diagnostic; the full declared alphabet was not supplied")
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    lm_to_tail = core.lm_to_tail_Q(polys)

    log("building generator matrices over Q")
    gen_L = build_generator_matrices_Q(
        variables, standard, mon_idx, lm_to_tail, delta, is_match, match_pat,
    )

    log("computing traces")
    traces: list[Fraction] = [Fraction(0)] * N
    for l, m_l in enumerate(standard):
        if not m_l:
            traces[l] = Fraction(N)
            continue
        s = Fraction(0)
        for j in range(N):
            vec = left_mult_by_word_Q(gen_L, m_l, {j: Fraction(1)})
            s += vec.get(j, Fraction(0))
        traces[l] = s

    log(f"filling {N}×{N} trace form")
    T = [[Fraction(0)] * N for _ in range(N)]
    for i, m_i in enumerate(standard):
        for j in range(N):
            vec = left_mult_by_word_Q(gen_L, m_i, {j: Fraction(1)})
            T[i][j] = sum((c * traces[l] for l, c in vec.items()), Fraction(0))

    log("Gaussian elimination over Q")
    det, rank = gauss_rank_det_Q(T)
    return N, rank, det


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="path to f4ncgb basis file (char 0)")
    ap.add_argument("--max-dim", type=int, default=500)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--n", type=int, required=True, help="declared strand count for the full generator alphabet")
    args = ap.parse_args()

    def log(msg: str) -> None:
        if args.verbose:
            print(f"[{msg}]", file=sys.stderr, flush=True)

    polys = core.parse_basis_Q(args.basis)
    variables = Presentation(args.n).var_names()
    log(f"parsed {len(polys)} polys over Q, vars = {variables}")

    lms = [lm for lm, _ in polys]
    delta, is_match, _, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        delta, is_match, variables, max_total=args.max_dim,
    )
    if standard is None:
        print(f"[abort] dim > {args.max_dim}")
        return 1
    N = len(standard)
    log(f"N = {N}")

    mon_idx = {m: i for i, m in enumerate(standard)}
    N, rank, det = trace_form_rank_Q(polys, standard, mon_idx, log=log, variables=variables)
    defect = N - rank

    print()
    print("[CHAR-0 RESULT]")
    print(f"  N         = {N}")
    print(f"  rank(T)   = {rank}    (over Q)")
    print(f"  deficit   = {defect}  (radical interpretation requires a certified finite algebra)")
    print(f"  det T     = {det}")
    print()
    if defect == 0:
        print("[CONDITIONAL: certified regular representation required] det T != 0 over Q  ⇒  Z_n^variant is SEMISIMPLE over Q.")
    else:
        print(f"[CONDITIONAL: certified regular representation required] rank deficit {defect} over Q + char-0 trace-form criterion")
        print(f"          ⇒ dim J(Z_n^variant) = {defect} over Q.")
        print("          ⇒ Z_n^variant is GENUINELY NON-SEMISIMPLE over Q.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
