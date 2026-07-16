"""Systematic atom-based reconstruction of the generic Q(q,t) Groebner basis.

Every tail coefficient c(q,t) of the generic GB is num(q,t)/D(q,t) where D factors into a short
list of irreducible ATOMS (each atom is a reduction pivot; reductions are few, so the set is
small). We DISCOVER atoms (bounded rational fit on a probe sample -> factor denominators), form the
shared denominator D = lcm, VERIFY by fitting every coefficient as num/D (exact polynomial
interpolation with all-point check), and SELF-HEAL (enlarge D from any failing entry's own
denominator). Result: the generic GB + the atom set, validated at fresh points.

Reuses the reconstruction core in variants/jordan_regular_rep_qt.py (fit functions are generic).
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction as Fr
import sympy as sp
from f4ncgb_tools.variants import jordan_regular_rep_qt as jr
from .pointalg import point_gb
from .spec import Spec

q, t = jr.q, jr.t


def _parse(polys):
    return {lm: {mon: c for c, mon in terms[1:]} for lm, terms in polys}


def _atoms_of(D):
    """Irreducible factors of a sympy polynomial D (as a set of sympy exprs)."""
    return {sp.factor(f) for f, _ in sp.factor_list(D)[1]}


def wide_fit_denominator(fit_pts, fit_vals, all_pts, all_vals, maxnum=10, maxden=8):
    """Bivariate rational fit with WIDE numerator/denominator boxes (the small-shape probe in jr
    fails at n>=4 where numerators are high-degree). Returns the factored sympy denominator or None."""
    for total in range(0, maxnum + maxden + 1):
        for dq_n in range(min(total, maxnum) + 1):
            dt_n = total - dq_n
            if dt_n > maxnum:
                continue
            for ddt in range(min(total, maxden) + 1):
                for dq_d in range(ddt + 1):
                    dt_d = ddt - dq_d
                    if dq_n + dt_n + dq_d + dt_d == 0:
                        continue
                    res = jr.fit_rational(fit_pts, fit_vals, all_pts, all_vals, dq_n, dt_n, dq_d, dt_d)
                    if res is not None:
                        _, _, den_c, den_m = res
                        den = sum(sp.Rational(b.numerator, b.denominator) * q ** i * t ** j
                                  for b, (i, j) in zip(den_c, den_m))
                        return sp.factor(den)
    return None


def _Dvals(D, pts):
    f = sp.lambdify((q, t), D, "math") if False else None
    out = []
    for (qv, tv) in pts:
        val = D.subs({q: qv, t: tv})
        r = sp.nsimplify(val); out.append(Fr(int(sp.numer(r)), int(sp.denom(r))))
    return out


def reconstruct_gb(spec: Spec, scratch: str, grid=None, n_probe: int = 16, log=print):
    if grid is None:
        grid = jr.grid_points()
    gbs = {}
    for (qv, tv) in grid:
        try:
            gbs[(qv, tv)] = _parse(point_gb(spec, qv, tv, scratch, char=0))
        except Exception:
            pass
    skel = Counter(frozenset(g) for g in gbs.values())
    gen = skel.most_common(1)[0][0]
    all_pts = [pt for pt, g in gbs.items() if frozenset(g) == gen]
    LMS = sorted(gen)
    log(f"[{spec.name}] {len(gbs)}/{len(grid)} solved, {len(all_pts)} generic, |GB|={len(LMS)}, "
        f"skeleton multiplicity {len(skel)}")
    nfit = min(len(all_pts), max(30, len(all_pts) * 2 // 3))
    fit_pts = all_pts[:nfit]

    slots = {}
    for lm in LMS:
        mons = set()
        for pt in all_pts:
            mons |= set(gbs[pt][lm].keys())
        for mon in mons:
            slots[(lm, mon)] = [gbs[pt][lm].get(mon, Fr(0)) for pt in all_pts]

    # --- discover atoms from the richest (most-varying) slots via a WIDE rational fit ---
    rich = sorted(slots.items(), key=lambda kv: -len(set(kv[1])))[:n_probe]
    D = sp.Integer(1)
    for (lm, mon), vals in rich:
        if all(v == vals[0] for v in vals):
            continue
        d = wide_fit_denominator(fit_pts, vals[:nfit], all_pts, vals)
        if d is not None:
            D = sp.lcm(D, d)
    atoms = _atoms_of(D)
    log(f"[{spec.name}] discovered atoms: {sorted(map(str, atoms))}; shared D = {sp.factor(D)}")

    # --- fit every slot as num/D, self-heal on failure ---
    GEN = {lm: {} for lm in LMS}
    Dfit = _Dvals(D, fit_pts); Dall = _Dvals(D, all_pts)
    healed = 0
    for (lm, mon), vals in slots.items():
        if all(v == vals[0] for v in vals):
            if vals[0] != 0:
                GEN[lm][mon] = sp.Rational(vals[0].numerator, vals[0].denominator)
            continue
        expr = jr.fit_entry_with_denominator(fit_pts, vals[:nfit], all_pts, vals, D, Dfit, Dall,
                                             max_dq=14, max_dt=14)
        if expr is None:
            d = wide_fit_denominator(fit_pts, vals[:nfit], all_pts, vals)
            if d is not None:
                D = sp.lcm(D, d); atoms = _atoms_of(D)
                Dfit = _Dvals(D, fit_pts); Dall = _Dvals(D, all_pts); healed += 1
            expr = jr.fit_entry_with_denominator(fit_pts, vals[:nfit], all_pts, vals, D, Dfit, Dall,
                                                 max_dq=14, max_dt=14)
        if expr is None:
            raise RuntimeError(f"reconstruction failed for {lm} :: {mon}")
        GEN[lm][mon] = sp.cancel(expr)
    if healed:
        log(f"[{spec.name}] self-healed {healed} times; final atoms: {sorted(map(str, atoms))}")

    return {"spec": spec.name, "n": spec.n, "LMS": [list(m) for m in LMS],
            "GEN": GEN, "atoms": sorted(map(str, atoms)), "D": D, "all_pts": all_pts}


def validate(spec: Spec, scratch: str, recon: dict, fresh=((19, 23), (25, 29), (31, 37))):
    """Check the reconstructed GB against fresh f4ncgb points off the grid."""
    LMS = [tuple(m) for m in recon["LMS"]]; GEN = recon["GEN"]
    bad = 0
    for (qv, tv) in fresh:
        g = _parse(point_gb(spec, qv, tv, scratch, char=0))
        if frozenset(g.keys()) != frozenset(LMS):
            return False, f"skeleton differs at ({qv},{tv})"
        for lm in LMS:
            for mon in set(g[lm]) | set(GEN[lm]):
                got = GEN[lm].get(mon, sp.Integer(0))
                gv = got.subs({q: qv, t: tv}) if hasattr(got, "subs") else got
                gr = sp.nsimplify(gv); gf = Fr(int(sp.numer(gr)), int(sp.denom(gr)))
                if gf != g[lm].get(mon, Fr(0)):
                    bad += 1
    return bad == 0, f"{bad} mismatches over {len(fresh)} fresh points"
