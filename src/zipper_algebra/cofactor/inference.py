"""Heuristic denominator-factor inference from exact parameter samples.

Line reconstruction and factor lifting propose an ordered denominator family.
They do not establish a generic identity. Optional Sympy factorization is
imported only by the residual-factor operation that needs it.
"""
from __future__ import annotations

import math
from fractions import Fraction
from itertools import product

from .atoms import Atom, AtomSet, poly_str
from .univariate import _poly_divmod, _specialise, _solve, reconstruct_univariate


Rational = int | Fraction


def _candidate_polys(max_bideg=(2, 1), coef_bound=2):
    """Primitive polynomials of bidegree <= max_bideg with |coefficients| <= coef_bound,
    non-constant, sign-normalised, in ascending (total degree, coefficient size)."""
    mons = [(i, j) for i in range(max_bideg[0] + 1) for j in range(max_bideg[1] + 1)]
    mons.sort(reverse=True)
    rng = range(-coef_bound, coef_bound + 1)
    out = []
    for cs in product(rng, repeat=len(mons)):
        coefs = {m: c for m, c in zip(mons, cs) if c}
        if not coefs or list(coefs) == [(0, 0)]:
            continue
        if coefs[max(coefs)] < 0:
            continue
        if math.gcd(*[abs(c) for c in coefs.values()]) != 1:
            continue
        out.append((max(i + j for i, j in coefs), sum(abs(c) for c in coefs.values()), coefs))
    out.sort(key=lambda x: (x[0], x[1]))
    return [c for _, _, c in out]


def _primitive(coefs: list) -> tuple:
    """A rational coefficient list -> the primitive integer one, sign-normalised on the lead."""
    from fractions import Fraction as F
    cs = [F(c) for c in coefs]
    den = 1
    for c in cs:
        den = den * c.denominator // math.gcd(den, c.denominator)
    ints = [int(c * den) for c in cs]
    g = 0
    for v in ints:
        g = math.gcd(g, abs(v))
    if g:
        ints = [v // g for v in ints]
    if ints and ints[-1] < 0:
        ints = [-v for v in ints]
    return tuple(ints)


def residual_factors(lines, verbose=True):
    """Factor every leftover line denominator over Q and report the irreducible factors.

    `lines` is the [(fixed_var, x0, D)] left by infer() after stripping the atoms it accepted.
    The brute-force candidate search cannot be widened far -- it enumerates all coefficient
    vectors in a bidegree box, which is exponential -- so when it fails the useful thing is not a
    bigger box but the ACTUAL factors, which sympy gets by factoring each univariate residual.
    Returns {(fixed_var, primitive coefficient tuple): [x0, ...]}.
    """
    import sympy as sp
    x = sp.Symbol("x")
    out: dict[tuple[str, tuple[int, ...]], list[Rational]] = {}
    for fixed, x0, D in lines:
        if len(D) < 2:
            continue
        poly = sp.Poly(list(reversed([sp.Rational(c) for c in D])), x)
        for f, _e in sp.factor_list(poly)[1]:
            factor_coefs = list(reversed(f.all_coeffs()))
            out.setdefault((fixed, _primitive(factor_coefs)), []).append(x0)
    if verbose and out:
        print(f"[atoms] residual irreducible factors ({len(out)} distinct):")
        for (fixed, cs), xs in sorted(out.items(), key=lambda kv: -len(kv[1]))[:8]:
            var = "q" if fixed == "t" else "t"
            print(f"[atoms]   {_fmt_univariate(cs, var)}   on {len(xs)} lines "
                  f"({'t' if fixed == 't' else 'q'} = {sorted(set(xs))[:5]}...)")
    return out


def _fmt_univariate(cs, var) -> str:
    parts = []
    for i, c in enumerate(cs):
        if not c:
            continue
        mon = "" if i == 0 else (var if i == 1 else f"{var}^{i}")
        body = f"{abs(c)}" if not mon else (mon if abs(c) == 1 else f"{abs(c)}*{mon}")
        parts.append(("-" if c < 0 else "+") + body)
    return "".join(parts).lstrip("+") or "0"


def atoms_from_residual_factors(factors, max_other_deg=4, min_lines=3, verbose=True):
    """Lift line factors to bivariate atoms.

    A factor of the residual on the line t = t0 is the restriction a(q, t0) of an atom a(q,t).
    Group the factors by their degree in q, and for each coefficient position interpolate its
    values across the lines as a polynomial in t (exactly, verified on held-out lines).  This
    is performed one variable at a time up to max_other_deg, with held-out lines.
    At least max(min_lines + 2, max_other_deg + 1) distinct usable coordinates
    are required. Ambiguous same-degree factor matches are left unresolved.
    """
    from fractions import Fraction as F
    by_deg: dict[tuple[str, int], list[tuple[tuple[int, ...], list[Rational]]]] = {}
    for (fixed, cs), xs in factors.items():
        by_deg.setdefault((fixed, len(cs)), []).append((cs, xs))
    found = []
    for (fixed, ln), group in sorted(by_deg.items()):
        pts: dict[Rational, set[tuple[int, ...]]] = {}
        for cs, xs in group:
            for x0 in xs:
                pts.setdefault(x0, set()).add(cs)
        usable = {x0: next(iter(v)) for x0, v in pts.items() if len(v) == 1}
        if len(usable) < max(min_lines + 2, max_other_deg + 1):
            continue
        xs = sorted(usable)
        coefs = {}
        ok = True
        for pos in range(ln):
            data = [(x0, F(usable[x0][pos])) for x0 in xs]
            fit = _fit_univariate_exact(data, max_other_deg)
            if fit is None:
                ok = False
                break
            coefs[pos] = fit
        if not ok:
            continue
        # assemble: coefficient of q^pos is a polynomial in t (or vice versa)
        bi: dict[tuple[int, int], Fraction] = {}
        for pos, poly in coefs.items():
            for d, c in poly.items():
                key = (pos, d) if fixed == "t" else (d, pos)
                bi[key] = bi.get(key, F(0)) + c
        bi = {k: v for k, v in bi.items() if v}
        if not bi or all(k == (0, 0) for k in bi):
            continue
        den = 1
        for c in bi.values():
            den = den * F(c).denominator // math.gcd(den, F(c).denominator)
        primitive_bi = {k: int(F(v) * den) for k, v in bi.items()}
        a = Atom(poly_str(primitive_bi), primitive_bi)
        if a not in found:
            found.append(a)
            if verbose:
                print(f"[atoms]   lifted to a bivariate atom: {a.name}")
    return found


def _fit_univariate_exact(data, max_deg):
    """Fit {x: value} by a polynomial in x of degree <= max_deg, exactly, with two held out."""
    from fractions import Fraction as F
    data = sorted(data)
    if len(data) < 3:
        return None
    fit = data[:-2]
    for deg in range(max_deg + 1):
        if len(fit) < deg + 1:
            break
        rows = [[F(x) ** i for i in range(deg + 1)] + [F(v)] for x, v in fit[: deg + 1]]
        sol = _solve(rows, deg + 1)
        if sol is None:
            continue
        good = all(sum(sol[i] * F(x) ** i for i in range(deg + 1)) == F(v)
                   for x, v in data)
        if good:
            return {i: sol[i] for i in range(deg + 1) if sol[i]}
    return None


def infer(values_by_slot: dict, seed: AtomSet | None = None, max_bideg=(2, 1), coef_bound=2,
          min_lines=3, verbose=True):
    """Infer the atom set from sampled coefficient values.

    Exact univariate reconstruction supplies candidate line denominators. Factoring
    those lines first proposes repeated one-variable factors and then lifts
    parameter-dependent factors. Bounded polynomial enumeration is a fallback
    when factorization is unavailable or leaves ambiguous factors. Its cost grows
    exponentially in the bidegree and coefficient bounds.

    values_by_slot : {slot: [((q, t), Fraction), ...]} -- for every coefficient the lift must
                     fit, its value at every grid point.  Grid points sharing a q or a t form
                     the lines along which denominators are reconstructed.
    seed           : atoms assumed present.  None = start empty (the honest test).

    Returns (atomset, report).  Seed positions are unchanged; discovered atoms are appended.
    report['unexplained'] > 0 means some reconstructed line factors remain unresolved.
    More samples, less ambiguous factor matching or a revised search class may help.
    Zero unexplained factors does not prove generic completeness; every reconstructed
    coefficient and final certificate still needs its separate verification.
    report['residual_factors'] records a diagnostic snapshot before the latest
    factor-lifting pass; accepted factors can therefore remain in that list.
    """
    seed = seed if seed is not None else AtomSet([], [])
    # 1. reconstruct denominators along lines: list of (fixed_var, x0, D_univariate)
    lines = []
    n_rec = n_try = 0
    for slot, pts in values_by_slot.items():
        by_t: dict[Rational, list[tuple[Rational, Rational]]] = {}
        by_q: dict[Rational, list[tuple[Rational, Rational]]] = {}
        for (q, t), v in pts:
            by_t.setdefault(t, []).append((q, v))
            by_q.setdefault(q, []).append((t, v))
        for fixed, groups in (("t", by_t), ("q", by_q)):
            for x0, line in groups.items():
                if len(line) < 6:
                    continue
                n_try += 1
                rec = reconstruct_univariate(line)
                if rec is None:
                    continue
                n_rec += 1
                D = rec[1]
                if len(D) > 1:
                    lines.append((fixed, x0, D))
    if verbose:
        print(f"[atoms] reconstructed {n_rec}/{n_try} line denominators over "
              f"{len(values_by_slot)} slots; {len(lines)} non-trivial")
    report = {"seed": seed.names(), "lines_reconstructed": n_rec, "lines_tried": n_try,
              "found": [], "multiplicity": {}, "unexplained": 0, "residual_factors": [],
              "search": {"max_bideg": list(max_bideg), "coef_bound": coef_bound}}

    # 2. strip the seed atoms (polynomial division, recording multiplicities)
    mult = {a.name: 0 for a in seed.atoms}

    def strip(atom):
        m_max = 0
        for idx, (fixed, x0, D) in enumerate(lines):
            a1 = _specialise(atom.coefs, fixed, x0)
            if len(a1) < 2:
                continue
            m = 0
            while len(D) >= len(a1):
                qq, r = _poly_divmod(D, a1)
                if r:
                    break
                D, m = qq, m + 1
            lines[idx] = (fixed, x0, D)
            m_max = max(m_max, m)
        return m_max

    for a in seed.atoms:
        mult[a.name] = strip(a)

    # 3. Prefer factors read from the reconstructed denominators. Identical
    # univariate factors on several distinct parallel lines are cheap proposals;
    # stripping them removes ambiguity before parameter-dependent factor lifting.
    # These remain sampled hypotheses, subject to the later identity check.
    found = []
    known = set(seed.atoms)
    max_other_deg = 4
    if lines:
        phase = "constant"
        while True:
            try:
                factors = residual_factors(lines, verbose=False)
                if not factors:
                    break
                report["residual_factors"] = [
                    _fmt_univariate(cs, "q" if fixed == "t" else "t")
                    for fixed, cs in factors]
                if phase == "constant":
                    proposals = []
                    for (fixed, cs), xs in factors.items():
                        if len(set(xs)) >= max(max_other_deg + 1, min_lines):
                            polynomial = {(degree, 0) if fixed == "t" else (0, degree): c
                                          for degree, c in enumerate(cs) if c}
                            proposals.append(Atom(poly_str(polynomial), polynomial))
                else:
                    proposals = atoms_from_residual_factors(
                        factors, max_other_deg=max_other_deg, min_lines=min_lines, verbose=False)
            except Exception as error:
                # Optional factorization may fail; arithmetic failures while
                # stripping must propagate rather than conceal partial changes.
                if verbose:
                    print(f"[atoms]   (factor-first proposals unavailable: {error})")
                break
            progress = False
            for atom in proposals:
                if atom not in known:
                    multiplicity = strip(atom)
                    if multiplicity:
                        found.append(atom)
                        known.add(atom)
                        mult[atom.name] = multiplicity
                        progress = True
                        if verbose:
                            print(f"[atoms]   found {atom.name} from factored sample lines")
            if phase == "varying" and not progress:
                break
            phase = "varying"
    cands = (_candidate_polys(max_bideg, coef_bound)
             if any(len(D) >= 2 for _, _, D in lines) else [])
    while True:
        best: tuple[Atom, tuple[int]] | None = None
        for coefs in cands:
            a = Atom(poly_str(coefs), coefs)
            if a in known:
                continue
            hits = 0
            tvals: set[Rational] = set()
            qvals: set[Rational] = set()
            for fixed, x0, D in lines:
                if len(D) < 2:
                    continue
                a1 = _specialise(coefs, fixed, x0)
                if len(a1) < 2 or len(a1) > len(D):
                    continue
                if not _poly_divmod(D, a1)[1]:
                    hits += 1
                    (tvals if fixed == "t" else qvals).add(x0)
            dep_t = any(j for _, j in coefs)
            dep_q = any(i for i, _ in coefs)
            if hits < min_lines:
                continue
            if dep_t and len(tvals) < 2 and len(qvals) < 2:
                continue
            if dep_q and len(qvals) < 2 and len(tvals) < 2:
                continue
            key = (hits,)
            if best is None or key > best[1]:
                best = a, key
        if best is None:
            break
        best_atom, best_key = best
        found.append(best_atom); known.add(best_atom)
        mult[best_atom.name] = strip(best_atom)
        if verbose:
            print(f"[atoms]   found {best_atom.name}  (multiplicity <= {mult[best_atom.name]}, "
                  f"{best_key[0]} lines)")
        if all(len(D) < 2 for _, _, D in lines):
            break

    left = sum(1 for _, _, D in lines if len(D) >= 2)
    report["found"] = [a.name for a in found]
    report["multiplicity"] = mult
    report["unexplained"] = left
    while left:
        # The box can remove ambiguities that blocked the first lifting pass.
        # Repeat residual lifting while newly accepted factors expose others.
        try:
            rf = residual_factors(lines, verbose=verbose)
            lifted = atoms_from_residual_factors(
                rf, max_other_deg=max_other_deg, min_lines=min_lines, verbose=verbose)
        except Exception as e:                      # sympy missing, or a factorisation failure
            rf, lifted = {}, []
            if verbose:
                print(f"[atoms]   (residual factorisation unavailable: {e})")
        report["residual_factors"] = [_fmt_univariate(cs, "q" if f == "t" else "t")
                                      for (f, cs) in rf]
        progress = False
        for atom in lifted:
            if atom not in known:
                multiplicity = strip(atom)
                if multiplicity:
                    found.append(atom)
                    known.add(atom)
                    mult[atom.name] = multiplicity
                    progress = True
        left = sum(1 for _, _, D in lines if len(D) >= 2)
        report["found"] = [a.name for a in found]
        report["multiplicity"] = mult
        report["unexplained"] = left
        if not progress:
            break
    if verbose:
        if left:
            ex = next(D for _, _, D in lines if len(D) >= 2)
            print(f"[atoms] UNRESOLVED: {left} line denominators keep an unexplained factor "
                  f"at bidegree {max_bideg} / |coef| <= {coef_bound}, e.g. {[str(c) for c in ex]}")
        else:
            print(f"[atoms] every reconstructed denominator is a product of "
                  f"{seed.names() + report['found']}")
    out = seed.extended(found, provenance=f"inferred over {len(values_by_slot)} slots; seed {seed.names()}")
    # exponent box from observed multiplicities (+1 slack), never below the seed's box
    out.max_exp = [max(m0, mult.get(a.name, 0) + 1) for a, m0 in zip(out.atoms, out.max_exp)]
    return out, report
