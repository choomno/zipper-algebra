#!/usr/bin/env python3
"""Phase 5: prove the cofactor identity over Q(q,t), with q and t transcendental.

This is the certificate.  Everything upstream -- the sampled Groebner bases, the fitted
coefficients, the inferred atom set -- is search, and search can be wrong in ways more sampling
never reveals.  What makes a result a result is this check:

    sum_j  c_j(q,t) * L_j * rel_{i_j} * R_j  ==  target       in Q(q,t)<x_1..y_{n-1}>,

an identity between two free-algebra polynomials, verified exactly.  Denominators are cleared by
the LCM M of the c_j, so the identity becomes one in Z[q,t]<X>, and both sides are compared term by
term.  No specialisation, no genericity argument, no bad locus.

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
import sys
import time
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))
_TOOLS = os.path.dirname(_HERE)
for _p in (_HERE, _TOOLS):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cofactor_io import Cofactor, resolve_target                      # noqa: E402
from ncpoly import (NCPoly, qt_add, qt_const, qt_mul, qt_pow,          # noqa: E402
                    nc_add, nc_is_zero, nc_mul, nc_scalar, nc_size,
                    nc_sub, nc_zero)
from presentation import Presentation                                  # noqa: E402
from zipper_in_ms import NumericRing                                   # noqa: E402


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
    pres = pres or Presentation(cof.n)
    relations, sliders = pres.symbolic_ncpolys()
    tgt_name = target if target is not None else cof.target
    tgt = _symbolic_target(pres, sliders, tgt_name)
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
    if verbose:
        print(f"[lcm]   M = {atoms.denom_str(lcm)}")
        print(f"[lcm]   M has {len(M_poly)} terms when expanded.")

    total: NCPoly = nc_zero()
    used = set()
    t_pre = time.time()
    for n_idx, (left, lab, right, num, den) in enumerate(cof.terms):
        if lab not in index_of:
            raise KeyError(f"relation label {lab!r} is not in presentation {pres.slug()}; "
                           f"the cofactor was made for a different algebra")
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
    """The target as an NCPoly over Q[q,t].  Accepts 'S[i][j]', 'X_k', (i, j), or an NCPoly."""
    if isinstance(name, dict):
        return name
    if isinstance(name, (tuple, list)) and len(name) == 2:
        key = (int(name[0]), int(name[1]))
    else:
        s = str(name).strip()
        if s.startswith("S[") and s.count("[") == 2:
            i, j = s[2:].replace("]", "").split("[")
            key = (int(i), int(j))
        elif s.upper().startswith("X_"):
            key = (int(s[2:]), 1)
        else:
            raise ValueError(f"cannot resolve target {name!r} symbolically")
    if key not in sliders:
        raise ValueError(f"target S{list(key)} does not exist at n={pres.n}")
    return sliders[key]


# ---------------------------------------------------------------------------
# Independent second verifier: deterministic PIT on an integer product grid
# ---------------------------------------------------------------------------

def pit_check(cof: Cofactor, pres: Presentation | None = None, target=None,
              verbose=True) -> dict:
    """Verify the cleared identity by evaluation on an integer grid, a posteriori.

    Both sides of `M * target == sum_j N_j * (M/D_j) * L_j * rel * R_j` are free-algebra
    polynomials whose word-coefficients lie in Z[q,t] of a bidegree that can be READ OFF the
    certificate.  A polynomial of bidegree (Dq, Dt) vanishing on a (Dq+1) x (Dt+1) integer
    product grid is identically zero, so agreement on such a grid is a proof, not a sample.
    Shares no arithmetic with `check` -- it only ever evaluates at integers.
    """
    t0 = time.time()
    pres = pres or Presentation(cof.n)
    tgt_name = target if target is not None else cof.target
    labels = pres.labels()
    index_of = {lab: i for i, lab in enumerate(labels)}
    atoms = cof.atoms
    lcm = cof.lcm_denominator()

    # a posteriori bidegree bound: deg(M) + max_j deg(N_j) + deg of the relations in (q,t)
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
    Dq, Dt = dq + nq + rq, dt + nt + rt
    grid_q = _grid_points(Dq + 1)
    grid_t = _grid_points(Dt + 1)
    if verbose:
        print(f"[pit]   bidegree bound ({Dq},{Dt}) read off the certificate; "
              f"grid {len(grid_q)}x{len(grid_t)} = {len(grid_q) * len(grid_t)} points")

    bad = []
    for q in grid_q:
        for t in grid_t:
            ring = NumericRing(q, t, 0)
            tgt = resolve_target(pres, tgt_name, ring)
            rels = pres.relations(ring)
            M = Fraction(atoms.eval_denom(lcm, q, t))
            acc = None
            for left, lab, right, num, den in cof.terms:
                c = sum(Fraction(v) * Fraction(q) ** i * Fraction(t) ** j
                        for (i, j), v in num.items())
                if not c:
                    continue
                scale = M / Fraction(atoms.eval_denom(den, q, t))
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
    from zipper_in_ms import SymbolicRing
    names, rels, S, tags = pres.build(SymbolicRing())
    return names, [{m: dict(c.v) for m, c in r.terms.items()} for r in rels], S, tags


def _word_poly(word, ring):
    from zipper_in_ms import Poly
    p = Poly(ring.char)
    p.add_term(tuple(word), ring.one)
    return p


def _grid_points(n: int) -> list:
    """n distinct small integers avoiding 0, 1 and -1 (they collapse atoms like q-1, q+1)."""
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

    cof = Cofactor.load(a.cofactor)
    pres = Presentation(a.n or cof.n, a.untwist, max_untwist_level=a.max_untwist_level)
    rep = check(cof, pres, a.target)
    if a.pit and rep["pass"]:
        rep2 = pit_check(cof, pres, a.target)
        rep["pass"] = rep["pass"] and rep2["pass"]
    raise SystemExit(0 if rep["pass"] else 1)
