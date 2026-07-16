"""PIT certification core (rigor backbone).

A MEMBERSHIP certificate for `target ∈ ideal(relations)` is a free-algebra identity
    sum_j  c_j(q,t) · L_j · rel_{i_j} · R_j  =  target          (in Z[q,t]<x_1..y_{n-1}>)
with c_j = M_j(q,t)/Den(q,t), M_j,Den in Z[q,t], L_j,R_j words, rel_i the FULL presentation
relations (braid + inverse + mixed + zipper + untwisting — all explicit, so NO group/Garside
normalization is needed: it is a pure free-algebra polynomial identity).

VERIFICATION = deterministic PIT. Both sides are free-algebra polynomials whose word-coefficients
are in Z[q,t]. Read the bidegree (Dq,Dt) off the certificate (a posteriori — no a-priori bound).
A bivariate polynomial of bidegree (Dq,Dt) that vanishes on a (Dq+1)x(Dt+1) integer product grid
is identically zero (iterated 1-var interpolation / finite differences) — EXACT, no probability, no
genericity, no bad locus. So: evaluate the cleared identity  sum_j M_j·L_j rel R_j = Den·target
at each grid point (integers), expand each side as a free-algebra poly (dict word->int), compare.
Agreement on the full grid == proof.

Relations are obtained (per integer point) from the tested emitter `zipper_in_ms.py --char 0`,
so we reuse validated code rather than re-deriving the presentation.
"""
from __future__ import annotations
import os, subprocess
from fractions import Fraction as Fr
import sympy as sp
from .spec import Spec

_HERE = os.path.dirname(os.path.abspath(__file__))
EMIT = os.path.normpath(f"{_HERE}/../zipper_in_ms.py")
q, t = sp.symbols("q t")

# ---------------------------------------------------------------------------
# free-algebra polynomials over a coefficient ring: dict {word(tuple[str]): coeff}
# ---------------------------------------------------------------------------
def fp_add(a, b, out=None):
    out = dict(a) if out is None else out
    for w, c in b.items():
        out[w] = out.get(w, 0) + c
        if out[w] == 0: del out[w]
    return out

def fp_scale(a, s):
    return {w: c * s for w, c in a.items() if c * s != 0}

def fp_mul(a, b):
    out = {}
    for wa, ca in a.items():
        for wb, cb in b.items():
            w = wa + wb; out[w] = out.get(w, 0) + ca * cb
            if out[w] == 0: del out[w]
    return out

def word(*vs):
    return {tuple(vs): 1}

# ---------------------------------------------------------------------------
# presentation relations at an integer point (q0,t0), as free-algebra polys over Z
# ---------------------------------------------------------------------------
def _parse_ms(path):
    lines = [l.strip().rstrip(",") for l in open(path) if l.strip()]
    rels = []
    for line in lines[2:]:  # skip variable header + modulus
        poly = {}
        for term in line.replace("-", "+-").split("+"):
            term = term.strip()
            if not term: continue
            factors = term.split("*")
            coeff = 1; vs = []
            for f in factors:
                f = f.strip()
                if f.lstrip("-").isdigit(): coeff *= int(f)
                else: vs.append(f)
            w = tuple(vs)
            poly[w] = poly.get(w, 0) + coeff
        rels.append({w: c for w, c in poly.items() if c})
    return rels

def presentation_at(spec: Spec, q0: int, t0: int, scratch: str):
    """List of relations (free-algebra polys over Z) of the FULL presentation at (q0,t0)."""
    ms = os.path.join(scratch, f"cert-{spec.name}-{q0}-{t0}.ms")
    subprocess.run(["python3", EMIT, "--n", str(spec.n), "--untwist", spec.emit_untwist,
                    "--framing-power", str(spec.framing_power), "--q", str(q0), "--t", str(t0),
                    "--char", "0", "--max-untwist-level", str(spec.max_untwist_level),
                    "--out", ms], check=True, capture_output=True)
    return _parse_ms(ms)

# ---------------------------------------------------------------------------
# certificate verification (deterministic grid PIT)
# ---------------------------------------------------------------------------
def _int_at(expr, q0, t0):
    if isinstance(expr, (int, Fr)):
        v = Fr(expr)
    else:
        v = sp.nsimplify(sp.sympify(expr).subs({q: q0, t: t0}))
        v = Fr(int(sp.numer(v)), int(sp.denom(v)))
    assert v.denominator == 1, f"cleared coeff not integer at ({q0},{t0}): {v}"
    return v.numerator

def verify_certificate(spec: Spec, cofactor, target_builder, scratch,
                       Dq: int, Dt: int, margin: int = 2, log=print):
    """cofactor: list of (M_expr, Den_expr, L_word, rel_index, R_word), the CLEARED identity being
        sum_j M_j · L_j rel_{i_j} R_j == Den · target,   Den common (pass as each tuple's Den_expr,
        must be identical). target_builder(q0,t0)->free-algebra poly over Z (e.g. X_k or g(Theta)X_k).
    Verify at a (Dq+1+margin)x(Dt+1+margin) integer grid avoiding q0 in {0,1,-1}, t0 in {0,1}.
    Returns (certified: bool, details)."""
    Dens = {sp.simplify(sp.sympify(c[1])) for c in cofactor}
    assert len(Dens) == 1, "certificate must have a single common denominator"
    Den = Dens.pop()
    qs = _pick(Dq + 1 + margin, skip={0, 1, -1})
    ts = _pick(Dt + 1 + margin, skip={0, 1})
    npts = 0
    for q0 in qs:
        for t0 in ts:
            rels = presentation_at(spec, q0, t0, scratch)
            lhs = {}
            for (M, _Den, L, i, R) in cofactor:
                m = _int_at(M, q0, t0)
                if m == 0: continue
                term = fp_mul(fp_mul(word(*L), rels[i]), word(*R))
                fp_add(lhs, fp_scale(term, m), out=lhs)
            den = _int_at(Den, q0, t0)
            rhs = fp_scale(target_builder(q0, t0), den)
            diff = fp_add(dict(lhs), fp_scale(rhs, -1))
            if any(v != 0 for v in diff.values()):
                return False, f"FAILED at (q,t)=({q0},{t0}); residual has {len([v for v in diff.values() if v])} words"
            npts += 1
    log(f"[certify {spec.name}] PIT PASSED on {npts} integer points (bidegree bound ({Dq},{Dt}))")
    return True, f"certified on {npts} points"

def _pick(n, skip):
    out = []; k = 2
    while len(out) < n:
        if k not in skip: out.append(k)
        k += 1
    return out
