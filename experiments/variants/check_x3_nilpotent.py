#!/usr/bin/env python3
r"""Probe whether the slider X_3 = S[3][1] is nilpotent in a variant zipper
quotient.

A non-zero nilpotent ``w`` in a finite-dimensional algebra ``A`` is one
necessary witness of non-semisimplicity, but not sufficient (e.g. E_12 in
M_n(k) is nilpotent yet M_n(k) is simple). This script reports the support
of X_3 and of X_3^k modulo the Gröbner basis; the verdict it prints flags
the case ``X_3 ≠ 0 yet X_3^k = 0`` for some small k as a useful
preliminary signal, with the rigorous semisimplicity verdict deferred to
the trace-form computation in ``semisimple_v3.py`` / ``semisimple_Q.py``.


Experimental interpretation: dimension, minimality and quotient-structure
conclusions require a certified complete basis or a separate exact witness.
A successful solver exit and agreement at several specializations do not supply
that evidence. State the corresponding evidence when interpreting these observations.
"""
from __future__ import annotations

import argparse


from zipper_algebra import core, Presentation
from zipper_algebra.linear import validate_field_characteristic


def build_X3_Fp(q: int, p: int) -> list[tuple[int, tuple[str, ...]]]:
    """Build X_3 = S[3][1] as a polynomial in the free algebra at the (q,t)
    specialisation, over F_p.

    Mirrors the slider recursion in ``zipper_algebra.relations``::

        S[2][1] = q*y_1 + (1-q) - x_1
        S[3][1] = (q^2*y_1*y_2 - x_1*x_2) * S[2][1].
    """
    S2 = [
        (q % p, ("y1",)),
        ((1 - q) % p, ()),
        ((-1) % p, ("x1",)),
    ]
    head = [
        (q * q % p, ("y1", "y2")),
        ((-1) % p, ("x1", "x2")),
    ]
    return [
        ((c1 * c2) % p, m1 + m2)
        for c1, m1 in head
        for c2, m2 in S2
    ]


def multiply_vec_Fp(
    v: dict[int, int],
    w: dict[int, int],
    standard: list[tuple[str, ...]],
    mon_idx: dict[tuple[str, ...], int],
    lm_to_tail,
    delta,
    is_match,
    match_pat,
    p: int,
) -> dict[int, int]:
    """Multiply two sparse vectors and reduce the result to standard form."""
    result: dict[int, int] = {}
    for i, ci in v.items():
        for j, cj in w.items():
            red = core.reduce_to_basis_Fp(
                [((ci * cj) % p, standard[i] + standard[j])],
                mon_idx, lm_to_tail, delta, is_match, match_pat, p,
            )
            for k, c in red.items():
                new = (result.get(k, 0) + c) % p
                if new == 0 and k in result:
                    del result[k]
                else:
                    result[k] = new
    return result


def vec_str(v, standard, max_terms=6):
    if not v:
        return "0"
    items = sorted(v.items(), key=lambda kc: (len(standard[kc[0]]), standard[kc[0]]))
    parts = [
        f"{c}*{'*'.join(standard[k]) if standard[k] else '1'}"
        for k, c in items[:max_terms]
    ]
    if len(items) > max_terms:
        parts.append(f"... ({len(items) - max_terms} more)")
    return " + ".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="path to f4ncgb basis file")
    ap.add_argument("--n", type=int, required=True, help="strand count")
    ap.add_argument("--q", type=int, default=5)
    ap.add_argument("--t", type=int, default=7)
    ap.add_argument("--char", type=int, default=2147483647)
    ap.add_argument("--max-power", type=int, default=8)
    ap.add_argument("--max-dim", type=int, default=5000)
    args = ap.parse_args()
    p = validate_field_characteristic(args.char)
    if not p:
        ap.error("this recipe requires a positive prime characteristic")

    polys = core.parse_basis_Fp(args.basis)
    print(f"[parse] {len(polys)} basis polynomials.")
    variables = Presentation(args.n).var_names()
    print(f"[vars] {variables}")

    lms = [lm for lm, _ in polys]
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        delta, is_match, variables, max_total=args.max_dim,
    )
    if standard is None:
        print(f"[abort] dim > {args.max_dim}")
        return 1
    print(f"[dim] N = {len(standard)}")

    mon_idx = {mon: i for i, mon in enumerate(standard)}
    lm_to_tail = core.lm_to_tail_Fp(polys, p)

    X3_poly = build_X3_Fp(args.q, p)
    print(f"[X_3] {len(X3_poly)} free-algebra terms before reduction.")

    v = core.reduce_to_basis_Fp(
        X3_poly, mon_idx, lm_to_tail, delta, is_match, match_pat, p)
    print(f"\nX_3       = {vec_str(v, standard)}  [supp={len(v)}]")

    cur = v
    for k in range(2, args.max_power + 1):
        cur = multiply_vec_Fp(
            cur, v, standard, mon_idx, lm_to_tail,
            delta, is_match, match_pat, p,
        )
        print(f"X_3^{k:<3}   = {vec_str(cur, standard)}  [supp={len(cur)}]")
        if not cur:
            print(f"\n[VERDICT] X_3 is nilpotent of index <= {k}.")
            if not v:
                print("          (But X_3 = 0 already, so this is the canonical case.)")
            else:
                print(f"          X_3 != 0 yet X_3^{k} = 0: candidate radical element.")
                print("          (Run semisimple_v3 / semisimple_Q for a rigorous verdict.)")
            return 0

    if not v:
        print("\n[VERDICT] X_3 = 0 (canonical BMW collapse).")
    else:
        print(
            f"\n[VERDICT] X_3^{args.max_power} != 0; X_3 not nilpotent within "
            "the power bound."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
