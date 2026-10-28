#!/usr/bin/env python3
"""Phase 5: prove the cofactor identity over Q(q,t), with q and t transcendental.

This is the certificate.  Everything upstream -- the sampled Groebner bases, the fitted
coefficients, the inferred atom set -- is search, and search can be wrong in ways more sampling
never reveals.  What makes a result a result is this check:

    sum_j  c_j(q,t) * L_j * rel_{i_j} * R_j  ==  target       in Q(q,t)<x_1..y_{n-1}>,

an identity between two free-algebra polynomials, verified exactly.  Denominators are cleared by
the common denominator M of the c_j, so the identity becomes one in Q[q,t]<X>, and both sides
are compared term by term. M must be a nonzero polynomial. The conclusion is over Q(q,t);
specializations of the displayed rational identity require its denominators to remain nonzero.

Ported from symbolic_identity_check.py, which was the trusted verifier but only for X_3 at n=4:
it carried a 130-line hand copy of the presentation (now `Presentation.symbolic_ncpolys`, the
single builder) and parsed denominators as a 4-tuple over four hard-coded atoms `q, t, q-1, t+q^2`
(now the cofactor's own AtomSet, whatever it contains).  The free-algebra arithmetic is the same
code, moved to `ncpoly.py` unchanged.

`pit_check` is a second, independent verifier along the lines of `tools/survey/certify.py`: read
the bidegree off the cleared identity a posteriori, evaluate both sides on an integer product grid
one larger in each variable, and conclude by interpolation.  It shares no arithmetic with `check`,
so agreement between them is real evidence rather than a repeated computation.
"""
from __future__ import annotations

import os
import time
from fractions import Fraction

from .cofactor_io import Cofactor, symbolic_target
from .ncpoly import (NCPoly, qt_add, qt_const, qt_mul, qt_pow, nc_add, nc_is_zero,
                     nc_mul, nc_scalar, nc_size, nc_sub, nc_zero)
from ..presentation import Presentation
from ..zipper_in_ms import NumericRing, SymbolicRing, Poly


def _atom_qt(atom) -> dict:
    """An Atom as a QTPoly."""
    return {(i, j): Fraction(c) for (i, j), c in atom.coefs.items()}


def _atom_powers(atoms, max_exps) -> list:
    """[{e: atom^e as QTPoly} for each atom], up to the LCM exponents."""
    out = []
    for atom, top in zip(atoms.atoms, max_exps):
        a = _atom_qt(atom)
        powers = {0: qt_const(1)}
        for e in range(1, top + 1):
            powers[e] = qt_mul(powers[e - 1], a)
        out.append(powers)
    return out


def check(cof: Cofactor, pres: Presentation | None = None, target=None, verbose=True) -> dict:
    """Verify the identity exactly over Q(q,t).  Returns a report; report['pass'] is the verdict."""
    t_start = time.time()
    if pres is None and not cof.presentation_hash:
        raise ValueError("a certificate without a presentation hash requires explicit pres=")
    pres = pres or Presentation(cof.n)
    # A cofactor is a certificate FOR one algebra. Relation LABELS coincide across presentations
    # that differ only in the untwisting (T_2 is T_2 in both bigelow and square-zero), so checking
    # labels is not enough -- it would "verify" a bigelow certificate against square-zero and fail
    # confusingly, or worse, pass. The presentation hash is the identity; a legacy file carries
    # none and is taken on trust, which is why new files always carry one.
    cof.validate(pres)
    relations, sliders = pres.symbolic_ncpolys()
    tgt_name = target if target is not None else cof.target
    tgt = symbolic_target(pres, tgt_name)
    labels = pres.labels()
    index_of = {lab: i for i, lab in enumerate(labels)}
    t_built = time.time()
    if verbose:
        print(f"[build] {len(relations)} input relations and {tgt_name} built in "
              f"{t_built - t_start:.2f}s.  {tgt_name} has {nc_size(tgt)} NC monomials.")
        print(f"[cof]   {len(cof.terms)} cofactor tuples.")

    atoms = cof.atoms
    lcm = cof.lcm_denominator()
    powers = _atom_powers(atoms, lcm)
    M_poly = qt_const(1)
    for pw, e in zip(powers, lcm):
        M_poly = qt_mul(M_poly, pw[e])
    # Validation gives nonzero atoms, so their product is nonzero in Q[q,t].
    if verbose:
        print(f"[lcm]   M = {atoms.denom_str(lcm)}")
        print(f"[lcm]   M has {len(M_poly)} terms when expanded.")

    total: NCPoly = nc_zero()
    used = set()
    t_pre = time.time()
    for n_idx, (left, lab, right, num, den) in enumerate(cof.terms):
        idx = index_of[lab]
        used.add(idx)
        # (M / D_j) as a QTPoly, then the coefficient N_j * (M / D_j)
        scale = qt_const(1)
        for pw, top, e in zip(powers, lcm, den):
            scale = qt_mul(scale, pw[top - e])
        coef_poly = qt_mul({k: Fraction(v) for k, v in num.items()}, scale)
        if not coef_poly:
            continue
        term = nc_mul(nc_mul({tuple(left): qt_const(1)}, relations[idx]),
                      {tuple(right): qt_const(1)})
        total = nc_add(total, nc_scalar(term, coef_poly))
        if verbose and (n_idx + 1) % 100 == 0:
            print(f"[acc]   {n_idx + 1}/{len(cof.terms)} tuples accumulated "
                  f"({time.time() - t_pre:.1f}s).  total has {nc_size(total)} NC monomials.")

    if verbose:
        print(f"[sum]   LHS built in {time.time() - t_built:.2f}s with "
              f"{nc_size(total)} NC monomials.")
    rhs = nc_scalar(tgt, M_poly)
    if verbose:
        print(f"[rhs]   M * {tgt_name} has {nc_size(rhs)} NC monomials.")
    diff = nc_sub(total, rhs)
    ok = nc_is_zero(diff)
    report = {"pass": ok, "tuples": len(cof.terms), "relations_used": sorted(used),
              "labels_used": sorted({cof.terms[i][1] for i in range(len(cof.terms))}),
              "lcm": atoms.denom_str(lcm), "lcm_exps": list(lcm),
              "residual_monomials": 0 if ok else nc_size(diff),
              "seconds": time.time() - t_start, "target": tgt_name,
              "presentation": pres.slug(), "presentation_hash": pres.content_hash()}
    if verbose:
        names = pres.var_names()
        print("\n" + "=" * 72)
        if ok:
            print(f"PASS: sum_j c_j(q,t) * L_j * R_{{i_j}} * R'_j  ==  {tgt_name}")
            print(f"      identity verified in the FREE algebra")
            print(f"      Q(q,t)<{','.join(names)}>.")
            print(f"      Tuples used:           {len(cof.terms)}")
            print(f"      Distinct R_i used:     {sorted(used)}")
            print(f"      LCM cleared:           M = {atoms.denom_str(lcm)}")
            print(f"      Total wall time:       {report['seconds']:.2f}s")
        else:
            print(f"FAIL: residual has {nc_size(diff)} non-zero NC monomial(s).")
            for mon, coef in sorted(diff.items(), key=lambda kv: (len(kv[0]), kv[0]))[:20]:
                mstr = "*".join(names[i - 1] for i in mon) if mon else "1"
                print(f"   {mstr}: {len(coef)} q,t-terms")
        print("=" * 72)
    return report


def _symbolic_target(pres: Presentation, sliders: dict, name) -> NCPoly:
    """Compatibility wrapper around the shared target resolver."""
    return symbolic_target(pres, name)


# ---------------------------------------------------------------------------
# Independent second verifier: deterministic PIT on an integer product grid
# ---------------------------------------------------------------------------

def pit_check(cof: Cofactor, pres: Presentation | None = None, target=None,
              verbose=True) -> dict:
    """Verify the cleared identity by evaluation on an integer grid, a posteriori.

    Both sides of `M * target == sum_j N_j * (M/D_j) * L_j * rel * R_j` are free-algebra
    polynomials whose word-coefficients lie in Q[q,t] of a bidegree that can be READ OFF the
    certificate.  A polynomial of bidegree (Dq, Dt) vanishing on a (Dq+1) x (Dt+1) integer
    product grid is identically zero, so agreement on such a grid is a proof, not a sample.
    Shares no arithmetic with `check` -- it only ever evaluates at integers.
    """
    t0 = time.time()
    if pres is None and not cof.presentation_hash:
        raise ValueError("a certificate without a presentation hash requires explicit pres=")
    pres = pres or Presentation(cof.n)
    cof.validate(pres)
    tgt_name = target if target is not None else cof.target
    tgt_sym = symbolic_target(pres, tgt_name)
    labels = pres.labels()
    index_of = {lab: i for i, lab in enumerate(labels)}
    atoms = cof.atoms
    lcm = cof.lcm_denominator()

    # Bound both the cleared cofactor sum and the cleared target separately.
    def bidegree(p: dict):
        return (max((i for i, _ in p), default=0), max((j for _, j in p), default=0))

    dq = dt = 0
    for a, e in zip(atoms.atoms, lcm):
        bi, bj = a.bidegree()
        dq += bi * e
        dt += bj * e
    nq = nt = 0
    for _l, _lab, _r, num, _d in cof.terms:
        bi, bj = bidegree(num)
        nq, nt = max(nq, bi), max(nt, bj)
    # relation coefficients: built from q, t and q^{i-1}; bound generously from the presentation
    _, rels_sym, _, _ = _symbolic_build(pres)
    rq = rt = 0
    for r in rels_sym:
        for coef in r.values():
            bi, bj = bidegree(coef)
            rq, rt = max(rq, bi), max(rt, bj)
    tq = tt = 0
    for coef in tgt_sym.values():
        bi, bj = bidegree(coef)
        tq, tt = max(tq, bi), max(tt, bj)
    Dq, Dt = max(dq + nq + rq, dq + tq), max(dt + nt + rt, dt + tt)
    grid_q = _grid_points(Dq + 1)
    grid_t = _grid_points(Dt + 1)
    if verbose:
        print(f"[pit]   bidegree bound ({Dq},{Dt}) read off the certificate; "
              f"grid {len(grid_q)}x{len(grid_t)} = {len(grid_q) * len(grid_t)} points")

    bad = []
    for q in grid_q:
        for t in grid_t:
            ring = NumericRing(q, t, 0)
            tgt = _evaluate_target(tgt_sym, ring, q, t)
            rels = pres.relations(ring)
            M = Fraction(atoms.eval_denom(lcm, q, t))
            acc = None
            for left, lab, right, num, den in cof.terms:
                c = sum(Fraction(v) * Fraction(q) ** i * Fraction(t) ** j
                        for (i, j), v in num.items())
                if not c:
                    continue
                # Evaluate the polynomial M/D directly. Dividing evaluations would be 0/0
                # at a grid point on a denominator's zero set, although the cleared
                # identity is a polynomial identity and is defined there.
                complement = tuple(top - e for top, e in zip(lcm, den))
                scale = Fraction(atoms.eval_denom(complement, q, t))
                lp, rp = _word_poly(left, ring), _word_poly(right, ring)
                piece = (lp * rels[index_of[lab]] * rp).scalar_mul(ring.scalar(c * scale))
                acc = piece if acc is None else acc + piece
            rhs = tgt.scalar_mul(ring.scalar(M))
            diff = (acc - rhs) if acc is not None else rhs
            if diff.terms:
                bad.append((q, t, len(diff.terms)))
    ok = not bad
    if verbose:
        print("[pit]   " + ("PASS: identity holds on the full grid -- exact by interpolation"
                            if ok else f"FAIL at {len(bad)} point(s), e.g. {bad[:3]}"))
    return {"pass": ok, "bidegree": [Dq, Dt], "grid": [len(grid_q), len(grid_t)],
            "failures": bad[:10], "seconds": time.time() - t0}


def _symbolic_build(pres: Presentation):
    names, rels, S, tags = pres.build(SymbolicRing())
    return names, [{m: dict(c.v) for m, c in r.terms.items()} for r in rels], S, tags


def _word_poly(word, ring):
    p = Poly(ring.char)
    p.add_term(tuple(word), ring.one)
    return p


def _evaluate_target(target, ring, q, t):
    p = Poly(0)
    for word, num in target.items():
        c = sum(Fraction(v) * Fraction(q) ** i * Fraction(t) ** j
                for (i, j), v in num.items())
        p.add_term(word, ring.scalar(c))
    return p


def _grid_points(n: int) -> list:
    """n distinct small integers. Cleared polynomials can also be evaluated at poles."""
    out, k = [], 2
    while len(out) < n:
        out.append(k)
        k += 1
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="verify a cofactor certificate over Q(q,t)")
    ap.add_argument("cofactor")
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--untwist", default="bigelow")
    ap.add_argument("--max-untwist-level", type=int, default=None)
    ap.add_argument("--target", default=None)
    ap.add_argument("--pit", action="store_true", help="also run the independent grid check")
    a = ap.parse_args()

    # This module entry point uses n=4 unless --n is supplied; certificate metadata
    # is checked against that explicitly constructed presentation.
    pres = Presentation(a.n or 4, a.untwist, max_untwist_level=a.max_untwist_level)
    cof = Cofactor.load(a.cofactor, pres=pres)
    tgt = a.target or cof.target
    if not tgt:
        raise SystemExit(f"{os.path.basename(a.cofactor)} carries no target in its header "
                         f"(a legacy file); pass --target, e.g. --target 'S[3][1]'")
    rep = check(cof, pres, tgt)
    if a.pit and rep["pass"]:
        rep2 = pit_check(cof, pres, tgt)
        rep["pass"] = rep["pass"] and rep2["pass"]
    raise SystemExit(0 if rep["pass"] else 1)
