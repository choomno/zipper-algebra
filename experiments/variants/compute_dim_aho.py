#!/usr/bin/env python3
r"""Hilbert series and dimension of a noncommutative quotient algebra, from a
Gröbner basis.

The standard monomials of ``k⟨x⟩/I`` are the free words not containing any
leading monomial of ``I`` as a contiguous substring. The set of such words
is a regular language: the complement of ``Σ* · L · Σ*`` where ``L`` is
the leading-monomial set, recognised by the Aho--Corasick automaton on
``L`` with the match states pruned. The Hilbert series is then the
generating function of length-d walks in the pruned subautomaton, computed
by an O(|states| · |alphabet|) recurrence per degree.

Usage::

    python3 compute_dim_aho.py path/to/basis.txt [--max-deg 40] [--vars "x1,..."]


Experimental interpretation: dimension, minimality and quotient-structure
conclusions require a certified complete basis or a separate exact witness.
A successful solver exit and agreement at several specializations do not supply
that evidence. State the corresponding evidence when interpreting these observations.
"""
from __future__ import annotations

import argparse
import os
import sys


from zipper_algebra import core, Presentation


def hilbert_series(
    delta: dict[tuple[int, str], int],
    is_match: list[bool],
    nstates: int,
    variables: list[str],
    max_deg: int,
) -> list[int]:
    """Return ``[h_0, h_1, ..., h_{max_deg}]``, where ``h_d`` is the number of
    standard monomials of degree ``d``.

    Iterates a sparse state-multiplicity vector ``f[s]`` (number of length-d
    walks from the root ending at state ``s``, avoiding match states); the
    Hilbert coefficient at degree ``d`` is ``sum(f)``. Stops early if
    ``f`` becomes identically zero.
    """
    counts = [0] * (max_deg + 1)
    if is_match[0]:
        return counts
    counts[0] = 1
    f = [0] * nstates
    f[0] = 1
    for d in range(1, max_deg + 1):
        nxt = [0] * nstates
        for s in range(nstates):
            if f[s] == 0 or is_match[s]:
                continue
            for v in variables:
                t = delta[(s, v)]
                if not is_match[t]:
                    nxt[t] += f[s]
        f = nxt
        counts[d] = sum(f)
        if counts[d] == 0 and not any(f):
            break
    return counts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="path to f4ncgb basis output")
    ap.add_argument("--max-deg", type=int, default=40)
    alphabet = ap.add_mutually_exclusive_group(required=True)
    alphabet.add_argument("--vars", help="complete ordered generator alphabet, comma-separated")
    alphabet.add_argument("--n", type=int, help="declared zipper strand count")
    args = ap.parse_args()

    if args.n is not None and args.n < 2:
        ap.error("--n must be at least two")
    lms = core.parse_leading_monomials(args.basis)
    print(f"Read {len(lms)} leading monomials.")
    variables = (
        args.vars.split(",") if args.vars is not None else Presentation(args.n).var_names()
    )
    if not variables or any(not name for name in variables) or len(set(variables)) != len(variables):
        ap.error("the ordered alphabet must contain distinct nonempty names")
    missing = {name for word in lms for name in word} - set(variables)
    if missing:
        ap.error(f"basis variables are absent from the declared alphabet: {sorted(missing)}")
    print(f"Variables (count {len(variables)}): {variables}")

    delta, is_match, _, nstates = core.build_aho_corasick(lms, variables)
    print(f"AC automaton: {nstates} states, "
          f"{sum(is_match)} match states excluded.")

    counts = hilbert_series(delta, is_match, nstates, variables, args.max_deg)
    total = sum(counts)
    last = max((d for d, c in enumerate(counts) if c != 0), default=0)
    print(f"Standard-word count (through degree {args.max_deg}): {total}")
    print(f"Leading-word avoidance counts (deg 0..{args.max_deg}): {counts}")
    print(f"  (last nonzero: deg {last}, count {counts[last]})")
    if last == args.max_deg and counts[last] > 0:
        print("  WARNING: last nonzero is at the truncation bound;")
        print("           algebra may be infinite-dim. Try higher --max-deg.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
