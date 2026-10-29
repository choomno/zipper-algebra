#!/usr/bin/env python3
r"""Verify that ``Z_{n−1}^{var} ↪ Z_n^{var}`` is an algebra embedding for the
framing and square-zero variants of the zipper algebra.

The map is induced by the natural braid-group inclusion ``B_{n−1} ↪ B_n``
(send ``σ_j`` to ``σ_j`` for ``1 ≤ j ≤ n−2``). At the variable level, the
small-algebra alphabet is a subset of the big-algebra alphabet, so each
small-basis monomial is also a word in the big algebra; reducing it
modulo the larger Gröbner basis gives its image in ``span(B_n)``. The map
is injective iff the resulting ``N_small × N_big`` matrix has rank
``N_small``.

The relations of ``Z_{n−1}`` are automatically in the ideal of ``Z_n``:
they involve only the small-alphabet letters and are by construction
among the generators of the bigger ideal (the relation list at ``n`` is a
superset of the relation list at ``n − 1``, modulo the new generator).


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


def check_tower(
    basis_small: str, basis_big: str, label: str, max_dim: int = 500,
    *, n_small: int = 3, n_big: int = 4,
) -> dict:
    print(f"\n=== {label} ===", flush=True)

    polys_small = core.parse_basis_Q(basis_small)
    polys_big = core.parse_basis_Q(basis_big)
    vars_small = Presentation(n_small).var_names()
    vars_big = Presentation(n_big).var_names()
    print(f"[parse] small: {len(polys_small)} polys, vars = {vars_small}",
          flush=True)
    print(f"[parse] big:   {len(polys_big)} polys, vars = {vars_big}",
          flush=True)

    missing = [v for v in vars_small if v not in vars_big]
    if missing:
        raise RuntimeError(
            f"Small-algebra variables not all present in big algebra: {missing}"
        )

    lms_small = [lm for lm, _ in polys_small]
    delta_s, im_s, mp_s, _ = core.build_aho_corasick(lms_small, vars_small)
    standard_small = core.enumerate_standard_monomials(
        delta_s, im_s, vars_small, max_total=max_dim,
    )
    if standard_small is None:
        raise RuntimeError(f"Small basis dim > {max_dim}")
    N_small = len(standard_small)
    print(f"[dim] small algebra: N = {N_small}", flush=True)

    lms_big = [lm for lm, _ in polys_big]
    delta_b, im_b, mp_b, _ = core.build_aho_corasick(lms_big, vars_big)
    standard_big = core.enumerate_standard_monomials(
        delta_b, im_b, vars_big, max_total=max_dim,
    )
    if standard_big is None:
        raise RuntimeError(f"Big basis dim > {max_dim}")
    N_big = len(standard_big)
    print(f"[dim] big algebra:   N = {N_big}", flush=True)

    mon_idx_big = {m: i for i, m in enumerate(standard_big)}
    lm_to_tail_big = core.lm_to_tail_Q(polys_big)

    print(f"[reduce] reducing {N_small} small-basis monomials in big algebra...",
          flush=True)
    rows: list[list[Fraction]] = []
    for m in standard_small:
        red = core.reduce_to_basis_Q(
            [(Fraction(1), m)],
            mon_idx_big, lm_to_tail_big, delta_b, im_b, mp_b,
        )
        row = [Fraction(0)] * N_big
        for idx, c in red.items():
            row[idx] = c
        rows.append(row)
    print("[reduce] done. Computing rank...", flush=True)

    rank = matrix_rank(rows)
    is_tower = rank == N_small
    print(f"[result] N_small = {N_small}, N_big = {N_big}, "
          f"rank of image = {rank}, is_tower = {is_tower}")
    return {
        "label": label,
        "N_small": N_small,
        "N_big": N_big,
        "rank": rank,
        "is_tower": is_tower,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--framing-small", required=True)
    ap.add_argument("--framing-big", required=True)
    ap.add_argument("--square-zero-small", required=True)
    ap.add_argument("--square-zero-big", required=True)
    ap.add_argument("--max-dim", type=int, default=500)
    args = ap.parse_args()

    results = [
        check_tower(args.framing_small, args.framing_big,
                    "FRAMING n=3 -> n=4 (a=2)", args.max_dim),
        check_tower(args.square_zero_small, args.square_zero_big,
                    "SQUARE-ZERO n=3 -> n=4", args.max_dim),
    ]

    print("\n=== SUMMARY ===")
    for r in results:
        print(f"  {r['label']}: N_small={r['N_small']}, N_big={r['N_big']}, "
              f"rank={r['rank']}, is_tower={r['is_tower']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
