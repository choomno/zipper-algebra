#!/usr/bin/env python3
r"""Candidate regular-representation matrices for the n=3 square-zero quotient.

The matrices have entries in Q(q,t), reconstructed from samples and checked
through exact matrix identities. Identifying them with the regular
representation additionally needs an independently proved upper dimension bound.

This is the regular-rep analogue of the 3-dimensional symbolic verification in
``variants/explicit_rep_n3.py``.  Where that script exhibits one small (3-dim)
representation and checks the relations on it, here we build candidate 24-by-24
matrices and verify every defining relation as an exact matrix identity.
The cyclic-vector check examines every listed word, including products of
several generators; checking generator images alone would be insufficient.

------------------------------------------------------------------------------
THE ALGEBRA
------------------------------------------------------------------------------
Generators x1, x2 (= sigma_1, sigma_2) and their inverses y1, y2.
Sliders   S2_j = q*y_j + (1-q) - x_j,   X_2 = S2_1.

Relations (square-zero variant; mirrors zipper_in_ms.build_relations(n=3,'square-zero')
and singular_run/square-zero-n3-symbolic.sing):

  R1: x_j y_j = y_j x_j = 1                       (j = 1, 2)
  R2: x1 x2 x1 = x2 x1 x2 ;  y1 y2 y1 = y2 y1 y2
  R4: x1 x2 y1 = y2 x1 x2 ;  y1 y2 x1 = x2 y1 y2 ;
      x1 y2 y1 = y2 y1 x2 ;  y1 x2 x1 = x2 x1 y2
  R5 (level-2 zipper, j=1):
      (q y2 + (1-q) - x2 - q y1 y2 + x1 x2) * S2_1 = 0
  R7 (square-zero, j=1,2):  (x_j - t)^2 * S2_j = 0
  R8 (level-3 untwist): (x2 x1 - q^2 y2 y1) * S3_1 = 0,
      with S3_1 = (q^2 y1 y2 - x1 x2) * S2_1
  (R3, R6 vacuous at n=3.)

------------------------------------------------------------------------------
CONVENTION  (homomorphism, not anti-homomorphism)
------------------------------------------------------------------------------
The regular module is A acting on itself by LEFT multiplication.  Define
L_a(m) = a*m for m in A.  Then
        L_a L_b (m) = a*(b*m) = (a b)*m = L_{ab}(m),
so  a |-> L_a  is an ALGEBRA HOMOMORPHISM:  L_{a b} = L_a L_b,  L_1 = I.
In the standard-monomial basis (e_l <-> m_l), column l of L_g is the
coefficient vector of the reduced product g*m_l:
        L_g[k][l] = c^{g}_{l,k},   where  g*m_l = sum_k c^{g}_{l,k} m_k.

A defining relation must act as the zero matrix. These exact checks establish
a representation when they all pass; inferring equality in the presented
algebra from equality of matrices would additionally require faithfulness.

------------------------------------------------------------------------------
PIPELINE
------------------------------------------------------------------------------
1. STANDARD BASIS.  Build the square-zero-3 .ms at reference (q,t)=(5,7) over Q,
   run f4ncgb (char 0), parse, enumerate the 24 standard monomials.  This is
   a reference candidate, not a certified generic basis. Sample points are
   selected when their listed standard words agree with this reference.

2. SAMPLE.  Over a grid of integer (q0,t0) points (avoiding q0 in {0,1,-1},
   t0 in {0,1}) build the GB over Q at each point, enumerate standard
   monomials, CHECK the 24-monomial set matches the reference exactly (else
   discard the point as bad-locus), and read off the rational coefficients
   c^{g}_{l,k}(q0,t0) of each generator action.

3. RATIONAL RECONSTRUCT.  For each (g,l,k) reconstruct c^{g}_{l,k}(q,t) as an
   element of Q(q,t) from the grid samples (dense bivariate rational
   interpolation, numerator/denominator of increasing degree until stable and
   held-out points reproduce exactly).  Assemble the four 24x24 sympy matrices
   L_{x1}, L_{x2}, L_{y1}, L_{y2}.

4. VERIFY (the rigorous core).  With q, t symbolic, check every defining
   relation as an EXACT zero-matrix identity over Q(q,t) using cancel/together.

5. CYCLIC VECTORS.  Check L_word e_unit = e_word for EVERY listed standard
   word, multiplying generator matrices in the declared order. If the defining
   relations and these checks all pass, the resulting cyclic representation
   gives the lower bound dim A >= N. A matching generic upper bound is separate;
   specialization and matching sample words do not establish it.

Run this producer to generate its pickle payload for the symbolic trace-form
recipe. The default output directory is runs/variants;
ZIPPER_EXPERIMENT_OUTPUT overrides it.

Run with the sympy venv:
    python -m experiments.variants.square_zero_regular_rep_qt


Experimental interpretation: dimension, minimality and quotient-structure
conclusions require a certified complete basis or a separate exact witness.
A successful solver exit and agreement at several specializations do not supply
that evidence. State the corresponding evidence when interpreting these observations.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from fractions import Fraction

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

from experiments.output import (prepare_output_dir, capture_experiment_provenance,
                                save_pickle_output, tag_solver_input, tag_solver_basis)
from zipper_algebra import core                            # noqa: E402
from zipper_algebra import Presentation

F4NCGB = os.environ.get("F4NCGB", "f4ncgb")
OUT_DIR = os.environ.get("ZIPPER_EXPERIMENT_OUTPUT",
                         os.path.join(os.path.dirname(os.path.dirname(HERE)), "runs", "variants"))
OUT_PKL = os.path.join(OUT_DIR, "square-zero-n3-regrep-qt.pkl")

GENERATORS = ["x1", "x2", "y1", "y2"]
MAXITER = 500
REF_Q, REF_T = 5, 7

q, t = sp.symbols("q t")


# ===========================================================================
# f4ncgb driver + standard-monomial machinery at a single integer (q0, t0)
# ===========================================================================
def write_square_zero_ms(path: str, q0: int, t0: int, *, provenance=True) -> None:
    """Emit the current square-zero presentation at the requested exact point."""
    Presentation(3, untwist="square-zero").emit_ms(path, q0, t0, 0, provenance=provenance)


def run_f4ncgb(ms_path: str, basis_path: str) -> None:
    producer = capture_experiment_provenance(__file__)
    inputs = tag_solver_input(ms_path, producer=producer)
    res = subprocess.run(
        [F4NCGB, "-i", ms_path, "-o", basis_path, "--maxiter", str(MAXITER),
         "-v", "0"],
        capture_output=True, text=True,
    )
    if res.returncode != 0 or not os.path.exists(basis_path):
        raise RuntimeError(
            f"f4ncgb failed (rc={res.returncode}) on {ms_path}\n"
            f"stderr:\n{res.stderr}\nstdout:\n{res.stdout}"
        )
    tag_solver_basis(basis_path, producer=producer, inputs=inputs)


class Quotient:
    """Standard-monomial model of the algebra at a fixed integer (q0, t0)."""

    def __init__(self, q0: int, t0: int, workdir: str):
        self.q0, self.t0 = q0, t0
        ms = os.path.join(workdir, f"square-zero-n3-Q-{q0}-{t0}.ms")
        basis = os.path.join(workdir, f"square-zero-n3-Q-{q0}-{t0}.basis")
        write_square_zero_ms(ms, q0, t0, provenance=False)
        run_f4ncgb(ms, basis)
        self.basis_path = basis
        self.polys = core.parse_basis_Q(basis)
        self.variables = list(GENERATORS)
        lms = [lm for lm, _ in self.polys]
        self.delta, self.is_match, self.match_pat, _ = core.build_aho_corasick(
            lms, self.variables
        )
        self.standard = core.enumerate_standard_monomials(
            self.delta, self.is_match, self.variables, max_total=500
        )
        if self.standard is None:
            self.dim = None
            return
        self.dim = len(self.standard)
        self.mon_idx = {m: i for i, m in enumerate(self.standard)}
        self.lm_to_tail = core.lm_to_tail_Q(self.polys)

    def gen_matrix_Q(self, g: str) -> list[list[Fraction]]:
        """24x24 matrix of left-multiplication by generator g over Q at (q0,t0).
        Column l = coefficient vector of the reduced product g*m_l."""
        N = self.dim
        M = [[Fraction(0)] * N for _ in range(N)]
        for l, m_l in enumerate(self.standard):
            red = core.reduce_to_basis_Q(
                [(Fraction(1), (g,) + m_l)], self.mon_idx, self.lm_to_tail,
                self.delta, self.is_match, self.match_pat,
            )
            for k, c in red.items():
                M[k][l] = c
        return M


# ===========================================================================
# Bivariate rational reconstruction of a single coefficient c(q, t) in Q(q,t)
# ===========================================================================
def _linsolve_exact(M, rhs):
    """Solve M x = rhs over Q (Fractions); return x or None if inconsistent."""
    n = len(M)
    if n == 0:
        return []
    cols = len(M[0])
    A = [list(r) + [rhs[i]] for i, r in enumerate(M)]
    pivot_row = 0
    pivot_cols = []
    for col in range(cols):
        if pivot_row >= n:
            break
        pr = None
        for r in range(pivot_row, n):
            if A[r][col] != 0:
                pr = r
                break
        if pr is None:
            continue
        A[pivot_row], A[pr] = A[pr], A[pivot_row]
        piv = A[pivot_row][col]
        A[pivot_row] = [v / piv for v in A[pivot_row]]
        for r in range(n):
            if r != pivot_row and A[r][col] != 0:
                f = A[r][col]
                A[r] = [A[r][c] - f * A[pivot_row][c] for c in range(cols + 1)]
        pivot_cols.append(col)
        pivot_row += 1
    for r in range(pivot_row, n):
        if A[r][cols] != 0:
            return None
    x = [Fraction(0)] * cols
    for i, col in enumerate(pivot_cols):
        x[col] = A[i][cols]
    return x


def _num_mons(dq, dt):
    return [(i, j) for i in range(dq + 1) for j in range(dt + 1)]


def _den_mons(dq, dt):
    """Denominator monomial exponents -- the FULL box q^i t^j (i<=dq, j<=dt),
    INCLUDING the constant term.  We do NOT force constant term = 1 (that
    would forbid denominators like q, t, q*t that vanish at the origin);
    instead the rational fit is solved as a homogeneous null-space problem and
    normalised afterwards."""
    return [(i, j) for i in range(dq + 1) for j in range(dt + 1)]


def _nullspace_vector(M):
    """Given a (consistent) homogeneous system M x = 0 with M a list of rows
    over Q, return a single non-trivial null vector if the null space is
    exactly 1-dimensional, else None (0-dim => only trivial; >1-dim =>
    ambiguous, treat as no clean fit).  Exact Fraction Gaussian elimination.
    """
    if not M:
        return None
    n_rows = len(M)
    n_cols = len(M[0])
    A = [list(r) for r in M]
    pivot_cols = []
    pr = 0
    for col in range(n_cols):
        piv = None
        for r in range(pr, n_rows):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        A[pr], A[piv] = A[piv], A[pr]
        inv = Fraction(1) / A[pr][col]
        A[pr] = [v * inv for v in A[pr]]
        for r in range(n_rows):
            if r != pr and A[r][col] != 0:
                f = A[r][col]
                A[r] = [A[r][c] - f * A[pr][c] for c in range(n_cols)]
        pivot_cols.append(col)
        pr += 1
        if pr == n_rows:
            break
    free_cols = [c for c in range(n_cols) if c not in pivot_cols]
    if len(free_cols) != 1:
        return None  # need exactly one free variable for a unique ratio
    free = free_cols[0]
    x = [Fraction(0)] * n_cols
    x[free] = Fraction(1)
    for i, pc in enumerate(pivot_cols):
        x[pc] = -A[i][free]
    return x


def fit_rational(fit_pts, fit_vals, all_pts, all_vals, dq_n, dt_n, dq_d, dt_d):
    """Fit value(q,t) = N(q,t)/D(q,t) exactly.  Homogeneous formulation:
    for each point, value*D(q,t) - N(q,t) = 0, with unknowns the coefficients
    of N (box dq_n x dt_n) and D (box dq_d x dt_d, constant term allowed).
    Solve the null space on ``fit_pts``; require a 1-dim ratio; verify on ALL
    points.  Returns (num_coeffs, num_mons, den_coeffs, den_mons) or None.
    """
    num_mons = _num_mons(dq_n, dt_n)
    den_mons = _den_mons(dq_d, dt_d)
    n_unk = len(num_mons) + len(den_mons)
    # homogeneous: need at least n_unk-1 independent equations to pin the ratio
    if len(fit_pts) < n_unk - 1:
        return None
    M = []
    for (qv, tv), c in zip(fit_pts, fit_vals):
        qv = Fraction(qv)
        tv = Fraction(tv)
        c = Fraction(c)
        row = [-(qv ** i * tv ** j) for (i, j) in num_mons]
        row += [c * qv ** i * tv ** j for (i, j) in den_mons]
        M.append(row)
    vec = _nullspace_vector(M)
    if vec is None:
        return None
    num_c = vec[:len(num_mons)]
    den_c = vec[len(num_mons):]
    # denominator must not be identically zero
    if all(b == 0 for b in den_c):
        return None
    # verify exact on ALL points (fit + held-out)
    for (qv, tv), c in zip(all_pts, all_vals):
        qv = Fraction(qv)
        tv = Fraction(tv)
        c = Fraction(c)
        num = sum(a * qv ** i * tv ** j for a, (i, j) in zip(num_c, num_mons))
        den = sum(b * qv ** i * tv ** j for b, (i, j) in zip(den_c, den_mons))
        if den == 0 or num != c * den:
            return None
    return num_c, num_mons, den_c, den_mons


# --- Shared-denominator reconstruction -------------------------------------
# All 24x24 entries of L_g are rational functions whose denominators divide a
# common polynomial D_g(q,t) (the product of the pivot denominators appearing
# in the standard-monomial reduction, a determinant-like quantity).  We
#   (a) discover D_g once from a handful of probe entries via a small bounded
#       rational fit, taking the lcm of the discovered denominators, then
#   (b) fit EVERY entry as numerator/D_g, where numerator = c*D_g is a pure
#       polynomial -- fast exact polynomial interpolation, no per-entry search.
# A final exact re-check on all sample points (train+holdout) certifies the
# fit; if any entry fails, D_g is enlarged and the pass repeats.

# Small shapes for the probe denominator discovery (numerator dq_n,dt_n ;
# denominator dq_d,dt_d), in increasing-complexity order.  Kept modest: the
# regular-rep entries have low-degree denominators.
_PROBE_SHAPES = []
for total in range(0, 7):
    for dq_n in range(0, total + 1):
        dt_n = total - dq_n
        for dd_total in range(0, 4):
            for dq_d in range(0, dd_total + 1):
                dt_d = dd_total - dq_d
                _PROBE_SHAPES.append((dq_n, dt_n, dq_d, dt_d))
_PROBE_SHAPES.sort(key=lambda s: ((s[0] + 1) * (s[1] + 1)
                                  + (s[2] + 1) * (s[3] + 1) - 1,
                                  s[2] + s[3], s[0] + s[1]))

# Wider shapes for the self-healing per-entry fallback fit.  Numerator is a
# full box [0..dq_n] x [0..dt_n]; denominator box similarly.  Ordered by total
# number of unknowns then by denominator size (prefer simpler denominators).
_ENTRY_SHAPES = []
for dq_n in range(0, 6):
    for dt_n in range(0, 6):
        for dq_d in range(0, 4):
            for dt_d in range(0, 4):
                _ENTRY_SHAPES.append((dq_n, dt_n, dq_d, dt_d))
_ENTRY_SHAPES.sort(key=lambda s: ((s[0] + 1) * (s[1] + 1)
                                  + (s[2] + 1) * (s[3] + 1) - 1,
                                  (s[2] + 1) * (s[3] + 1),
                                  s[0] + s[1]))


def _coeffs_to_expr(num_c, num_mons, den_c, den_mons):
    """Build sympy cancel(N/D) from null-space coefficient vectors."""
    num_e = sum(
        (sp.Rational(a.numerator, a.denominator) * q ** i * t ** j
         for a, (i, j) in zip(num_c, num_mons)), sp.Integer(0))
    den_e = sum(
        (sp.Rational(b.numerator, b.denominator) * q ** i * t ** j
         for b, (i, j) in zip(den_c, den_mons)), sp.Integer(0))
    return sp.cancel(num_e / den_e)


def probe_denominator(fit_pts, fit_vals, all_pts, all_vals):
    """Reconstruct one entry c(q,t) by a small bounded rational fit (fit on
    ``fit_pts`` subset, verify on all) and return its denominator as a sympy
    polynomial (factored).  None if the entry is constant."""
    if all(v == fit_vals[0] for v in fit_vals):
        return None
    for (dq_n, dt_n, dq_d, dt_d) in _PROBE_SHAPES:
        res = fit_rational(fit_pts, fit_vals, all_pts, all_vals,
                           dq_n, dt_n, dq_d, dt_d)
        if res is None:
            continue
        num_c, num_mons, den_c, den_mons = res
        expr = _coeffs_to_expr(num_c, num_mons, den_c, den_mons)
        _n, d = sp.fraction(sp.together(expr))
        d = sp.factor(d)
        return d
    raise RuntimeError("probe: no small rational shape fit this entry")


def _poly_mons(dq, dt):
    return [(i, j) for i in range(dq + 1) for j in range(dt + 1)]


def interp_polynomial(fit_pts, fit_vals, all_pts, all_vals, dq, dt):
    """Solve for poly(q,t) of per-variable degree (dq,dt) on the small
    ``fit_pts`` subset, then verify exactly on ALL points.  Returns sympy
    poly or None if no exact fit of that shape on every point."""
    mons = _poly_mons(dq, dt)
    if len(fit_pts) < len(mons):
        return None
    M, rhs = [], []
    for (qv, tv), c in zip(fit_pts, fit_vals):
        qv = Fraction(qv)
        tv = Fraction(tv)
        M.append([qv ** i * tv ** j for (i, j) in mons])
        rhs.append(Fraction(c))
    sol = _linsolve_exact(M, rhs)
    if sol is None:
        return None
    # exact verification on EVERY point (fit + held-out)
    for (qv, tv), c in zip(all_pts, all_vals):
        qv = Fraction(qv)
        tv = Fraction(tv)
        s = sum(a * qv ** i * tv ** j for a, (i, j) in zip(sol, mons))
        if s != Fraction(c):
            return None
    return sum(
        (sp.Rational(a.numerator, a.denominator) * q ** i * t ** j
         for a, (i, j) in zip(sol, mons)),
        sp.Integer(0),
    )


def fit_entry_with_denominator(fit_pts, fit_vals, all_pts, all_vals,
                               D, fit_Dvals, all_Dvals,
                               max_dq=10, max_dt=10):
    """Given common denominator D(q,t) (sympy), fit c = numer/D where
    numer = c*D is a polynomial.  Fit numerator coefficients on the small
    ``fit_pts`` subset and verify on ALL points.  Returns sympy expr or None.
    """
    fit_num = [Fraction(c) * dv for c, dv in zip(fit_vals, fit_Dvals)]
    all_num = [Fraction(c) * dv for c, dv in zip(all_vals, all_Dvals)]
    for total in range(0, max_dq + max_dt + 1):
        for dq in range(0, min(total, max_dq) + 1):
            dt = total - dq
            if dt > max_dt:
                continue
            poly = interp_polynomial(fit_pts, fit_num, all_pts, all_num, dq, dt)
            if poly is not None:
                return sp.cancel(poly / D)
    return None


# ===========================================================================
# Build the symbolic 24x24 generator matrices
# ===========================================================================
def grid_points():
    """13x13 grid of small integers, q0 in [2..] avoiding {-1,0,1},
    t0 avoiding {0,1}.  Mix of signs to pin down q- and t- atoms."""
    qs = [2, 3, 4, 6, 8, 9, 11, 13, -2, -3, -5, -7, 12]
    ts = [2, 3, 4, 5, 6, 8, 9, 10, -2, -3, -4, -6, 11]
    pts = [(qq, tt) for qq in qs for tt in ts]
    return pts


def sample_actions(log):
    """Run f4ncgb on the grid; return (ref_standard, samples) where samples is
    a list of (q0, t0, {g: 24x24 Fraction matrix}) for points that reproduce
    the reference 24-monomial basis."""
    with tempfile.TemporaryDirectory() as wd:
        log("building reference quotient at (q,t)=(5,7)")
        ref = Quotient(REF_Q, REF_T, wd)
        if ref.dim != 24:
            raise RuntimeError(f"reference dim = {ref.dim} != 24")
        ref_standard = list(ref.standard)
        ref_set = set(ref_standard)

        # record leading-coefficient info (always 1 after f4ncgb normalisation,
        # but confirm none is degenerate / the basis is monic) and GB size.
        log(f"reference GB has {len(ref.polys)} elements, dim 24")

        # include reference point itself as a sample
        samples = []
        all_pts = [(REF_Q, REF_T)] + [
            p for p in grid_points() if p != (REF_Q, REF_T)
        ]
        bad = []
        for (q0, t0) in all_pts:
            if q0 in (0, 1, -1) or t0 in (0, 1):
                continue
            try:
                quo = Quotient(q0, t0, wd)
            except RuntimeError as e:
                bad.append((q0, t0, f"f4ncgb error: {e}"))
                continue
            if quo.dim != 24 or list(quo.standard) != ref_standard:
                bad.append((q0, t0,
                            f"dim={quo.dim}, basis "
                            f"{'matches' if quo.dim == 24 and set(quo.standard) == ref_set else 'DIFFERS'}"))
                continue
            mats = {g: quo.gen_matrix_Q(g) for g in GENERATORS}
            samples.append((q0, t0, mats))
        log(f"usable sample points: {len(samples)}; "
            f"bad-locus/skipped: {len(bad)}")
        for (q0, t0, why) in bad:
            log(f"  bad point (q={q0}, t={t0}): {why}")
    return ref_standard, samples


def _eval_frac(expr, q0, t0):
    """Evaluate a sympy expr in q,t at integers; return a Fraction (or None
    if the denominator vanishes there).  Done exactly: substitute into the
    numerator and denominator polynomials separately."""
    num, den = sp.fraction(sp.together(expr))
    nv = num.subs({q: sp.Integer(q0), t: sp.Integer(t0)})
    dv = den.subs({q: sp.Integer(q0), t: sp.Integer(t0)})
    nv = sp.Rational(nv)
    dv = sp.Rational(dv)
    if dv == 0:
        return None
    r = nv / dv
    return Fraction(int(sp.numer(r)), int(sp.denom(r)))


def _split_fit_all(samples, fit_bound=7):
    """Return (fit_idx, all_idx, pts) where fit_idx selects small-magnitude
    points (|q|,|t| <= fit_bound) used to SOLVE the linear systems (keeps
    integer powers small / fast), and ALL points are used to VERIFY.
    """
    pts = [(q0, t0) for (q0, t0, _m) in samples]
    fit_idx = [i for i, (q0, t0) in enumerate(pts)
               if abs(q0) <= fit_bound and abs(t0) <= fit_bound]
    # ensure enough fit points; relax bound if needed
    bound = fit_bound
    while len(fit_idx) < 60 and bound < 50:
        bound += 2
        fit_idx = [i for i, (q0, t0) in enumerate(pts)
                   if abs(q0) <= bound and abs(t0) <= bound]
    all_idx = list(range(len(pts)))
    return fit_idx, all_idx, pts


def _entry_denominator(fit_pts, fit_vals, all_pts, all_vals):
    """Direct per-entry rational fit; return the entry's denominator (factored
    sympy poly) and the entry expression.  Wider shape search than the probe,
    used as the self-healing fallback when the running common denominator does
    not yet cover an entry.  Raises if nothing fits."""
    for (dq_n, dt_n, dq_d, dt_d) in _ENTRY_SHAPES:
        res = fit_rational(fit_pts, fit_vals, all_pts, all_vals,
                           dq_n, dt_n, dq_d, dt_d)
        if res is None:
            continue
        num_c, num_mons, den_c, den_mons = res
        expr = _coeffs_to_expr(num_c, num_mons, den_c, den_mons)
        _n, d = sp.fraction(sp.together(expr))
        return sp.factor(d), expr
    raise RuntimeError("entry: no rational shape fit")


def discover_common_denominator(g, samples, fit_idx, log):
    """Seed a common denominator D_g(q,t) by probing a spread of non-constant
    entries and lcm-ing their denominators.  The main build self-heals D by
    enlarging it whenever an entry needs an extra factor, so this is only a
    fast warm start."""
    N = 24
    pts = [(q0, t0) for (q0, t0, _m) in samples]
    fit_pts = [pts[i] for i in fit_idx]
    D = sp.Integer(1)
    n_probed = 0
    probe_cells = [(k, (k * 7 + 3) % N) for k in range(N)] + \
                  [((l * 5 + 2) % N, l) for l in range(N)]
    seen = set()
    for (k, l) in probe_cells:
        if (k, l) in seen:
            continue
        seen.add((k, l))
        all_vals = [mats[g][k][l] for (_q, _t, mats) in samples]
        if all(v == all_vals[0] for v in all_vals):
            continue
        fit_vals = [all_vals[i] for i in fit_idx]
        try:
            d = probe_denominator(fit_pts, fit_vals, pts, all_vals)
        except RuntimeError:
            continue
        n_probed += 1
        if d is not None and not d.is_number:
            D = sp.lcm(D, sp.factor(d))
        if n_probed >= 10:
            break
    D = sp.factor(sp.cancel(D))
    log(f"  L_{g}: warm-start common denominator D_{g} = {D}")
    return D


def _dvals_at(D, pts):
    out = []
    for (q0, t0) in pts:
        dv = _eval_frac(D, q0, t0)
        if dv is None or dv == 0:
            return None
        out.append(dv)
    return out


def build_symbolic_matrices(ref_standard, samples, log):
    """Reconstruct each L_g as a 24x24 sympy matrix over Q(q,t).  A common
    denominator D_g is grown on demand: each entry is fit as numerator/D_g,
    and if that fails the entry's own denominator is computed and multiplied
    into D_g (self-healing).  Every entry is verified EXACTLY on all points."""
    N = len(ref_standard)
    fit_idx, all_idx, pts = _split_fit_all(samples)
    fit_pts = [pts[i] for i in fit_idx]
    log(f"reconstruction: solving on {len(fit_pts)} small-magnitude points, "
        f"verifying on all {len(pts)}")

    L = {}
    for g in GENERATORS:
        D = discover_common_denominator(g, samples, fit_idx, log)
        all_Dvals = _dvals_at(D, pts)
        if all_Dvals is None:
            D = sp.Integer(1)
            all_Dvals = [Fraction(1)] * len(pts)
        fit_Dvals = [all_Dvals[i] for i in fit_idx]

        Mexpr = sp.zeros(N, N)
        n_const = n_rat = n_heal = 0
        for k in range(N):
            for l in range(N):
                all_vals = [mats[g][k][l] for (_q, _t, mats) in samples]
                if all(v == all_vals[0] for v in all_vals):
                    v0 = all_vals[0]
                    Mexpr[k, l] = sp.Rational(v0.numerator, v0.denominator)
                    n_const += 1
                    continue
                fit_vals = [all_vals[i] for i in fit_idx]
                expr = fit_entry_with_denominator(
                    fit_pts, fit_vals, pts, all_vals,
                    D, fit_Dvals, all_Dvals)
                if expr is None:
                    # self-heal: discover this entry's denominator and grow D
                    d_entry, e_entry = _entry_denominator(
                        fit_pts, fit_vals, pts, all_vals)
                    Dnew = sp.factor(sp.cancel(sp.lcm(D, d_entry)))
                    dvn = _dvals_at(Dnew, pts)
                    if dvn is None:
                        # Dnew vanishes somewhere on the grid; keep the direct
                        # per-entry expression and don't fold into D.
                        Mexpr[k, l] = e_entry
                        n_heal += 1
                        n_rat += 1
                        continue
                    D, all_Dvals = Dnew, dvn
                    fit_Dvals = [all_Dvals[i] for i in fit_idx]
                    expr = fit_entry_with_denominator(
                        fit_pts, fit_vals, pts, all_vals,
                        D, fit_Dvals, all_Dvals)
                    if expr is None:
                        expr = e_entry
                    n_heal += 1
                Mexpr[k, l] = expr
                n_rat += 1
        log(f"  L_{g}: grown common denominator D_{g} = {D}  "
            f"({n_heal} self-heal steps)")
        # EXACT validation of the whole matrix on every sample point
        _validate_matrix(g, Mexpr, samples, log)
        L[g] = Mexpr
        log(f"  L_{g}: reconstructed and validated "
            f"({n_const} constant, {n_rat} rational entries)")
    return L


def _validate_matrix(g, Mexpr, samples, log):
    """Re-evaluate the symbolic matrix at every sample point and confirm it
    equals the directly computed Fraction matrix.  Raises on mismatch."""
    N = Mexpr.shape[0]
    # cache per-entry sympy lambda-free eval via subs (entries are simple)
    for (q0, t0, mats) in samples:
        sub = {q: sp.Integer(q0), t: sp.Integer(t0)}
        for k in range(N):
            for l in range(N):
                e = Mexpr[k, l]
                if e.is_number:
                    got = Fraction(int(sp.numer(e)), int(sp.denom(e)))
                else:
                    fv = _eval_frac(e, q0, t0)
                    if fv is None:
                        raise RuntimeError(
                            f"L_{g}[{k},{l}] denominator vanishes at "
                            f"(q={q0}, t={t0}) during validation")
                    got = fv
                want = mats[g][k][l]
                if got != want:
                    raise RuntimeError(
                        f"VALIDATION FAIL L_{g}[{k},{l}] at (q={q0},t={t0}): "
                        f"reconstructed {got} != sampled {want}")


# ===========================================================================
# Symbolic relation checks  (homomorphism convention: L_{ab} = L_a L_b)
# ===========================================================================
def is_zero_matrix(M):
    Mc = M.applyfunc(lambda e: sp.cancel(sp.together(e)))
    return Mc == sp.zeros(*M.shape)


def check_cyclic_basis(matrices, standard):
    """Check every named basis vector as L_word applied to the unit vector.

    Columns encode left actions, so L_(g1...gk) = L_g1 ... L_gk acts on a
    vector by applying the rightmost generator first. This check concerns
    cyclicity only; defining-relation checks are separate.
    """
    words = [tuple(word) for word in standard]
    if () not in words:
        raise ValueError("the listed standard words must include the empty word")
    n = len(words)
    if any(matrix.shape != (n, n) for matrix in matrices.values()):
        raise ValueError("generator matrix size does not match the listed words")
    unit = sp.zeros(n, 1)
    unit[words.index(()), 0] = 1
    for index, word in enumerate(words):
        column = unit
        for generator in reversed(word):
            column = matrices[generator] * column
        expected = sp.zeros(n, 1)
        expected[index, 0] = 1
        if not is_zero_matrix(column - expected):
            return False
    return True


def relation_checks(L):
    """Return ordered list of (name, matrix-that-must-be-zero).

    Words are mapped to operators via the homomorphism  L_{g1 g2 ... gr}
    = L_{g1} L_{g2} ... L_{gr}.  Inverses come from the relation x_j y_j = 1,
    but for the operator algebra we use L_{y_j} directly (the reconstructed
    matrix); R1 below certifies L_{x_j} L_{y_j} = I so the two notions agree.
    """
    N = L["x1"].shape[0]
    I = sp.eye(N)
    X1, X2, Y1, Y2 = L["x1"], L["x2"], L["y1"], L["y2"]

    # Operator form of sliders.  S2_j as an element is q y_j + (1-q) - x_j,
    # so L_{S2_j} = q L_{y_j} + (1-q) I - L_{x_j}.
    S2_1 = q * Y1 + (1 - q) * I - X1
    S2_2 = q * Y2 + (1 - q) * I - X2
    # S3_1 = (q^2 y1 y2 - x1 x2) * S2_1  (an element);  L of a product of
    # elements = product of their L's (homomorphism):
    S3_1 = (q**2 * (Y1 * Y2) - (X1 * X2)) * S2_1

    rels = [
        # R1: inverses (also certifies invertibility of L_{x_j})
        ("R1  L_x1 L_y1 - I",            X1 * Y1 - I),
        ("R1  L_y1 L_x1 - I",            Y1 * X1 - I),
        ("R1  L_x2 L_y2 - I",            X2 * Y2 - I),
        ("R1  L_y2 L_x2 - I",            Y2 * X2 - I),
        # R2: braid
        ("R2  x1x2x1 - x2x1x2",          X1 * X2 * X1 - X2 * X1 * X2),
        ("R2  y1y2y1 - y2y1y2",          Y1 * Y2 * Y1 - Y2 * Y1 * Y2),
        # R4: mixed length-3 identities
        ("R4a x1x2y1 - y2x1x2",          X1 * X2 * Y1 - Y2 * X1 * X2),
        ("R4b y1y2x1 - x2y1y2",          Y1 * Y2 * X1 - X2 * Y1 * Y2),
        ("R4c x1y2y1 - y2y1x2",          X1 * Y2 * Y1 - Y2 * Y1 * X2),
        ("R4d y1x2x1 - x2x1y2",          Y1 * X2 * X1 - X2 * X1 * Y2),
        # R5: level-2 zipper.  Element:
        #   (q y2 + (1-q) - x2 - q y1 y2 + x1 x2) * S2_1 = 0
        ("R5  level-2 zipper (j=1)",
         (q * Y2 + (1 - q) * I - X2 - q * (Y1 * Y2) + (X1 * X2)) * S2_1),
        # R7: square-zero.  Element (x_j - t)^2 * S2_j = 0
        ("R7  (x1 - t)^2 S2_1",          (X1 - t * I) * (X1 - t * I) * S2_1),
        ("R7  (x2 - t)^2 S2_2",          (X2 - t * I) * (X2 - t * I) * S2_2),
        # R8: level-3 untwist.  Element (x2 x1 - q^2 y2 y1) * S3_1 = 0
        ("R8  (x2x1 - q^2 y2y1) S3_1",
         (X2 * X1 - q**2 * (Y2 * Y1)) * S3_1),
    ]
    return rels, S2_1, S2_2, S3_1


# ===========================================================================
# Main
# ===========================================================================
def main():
    producer = capture_experiment_provenance(__file__)
    def log(msg):
        print(f"[{msg}]", file=sys.stderr, flush=True)

    print("=" * 74)
    print("CANDIDATE REGULAR-REPRESENTATION MATRICES over Q(q,t)  (24x24)")
    print("=" * 74)

    # ---- step 1+2: reference basis and grid sampling -------------------
    ref_standard, samples = sample_actions(log)
    N = len(ref_standard)
    print(f"\n[STEP 1] STANDARD BASIS: dim = {N}   "
          f"({'PASS' if N == 24 else 'FAIL'})")
    print("  reference standard words (candidate basis, e_1 = empty word):")
    for i, m in enumerate(ref_standard):
        print(f"    e{i:2d}: {'*'.join(m) if m else '<1>'}")
    print(f"\n[STEP 2] SAMPLED generator actions at {len(samples)} grid "
          f"points reproducing the reference 24-monomial basis.")

    if N != 24:
        print("\n*** ABORT: reference dimension is not 24 ***")
        return 1

    # ---- step 3: rational reconstruction -------------------------------
    print("\n[STEP 3] RATIONAL RECONSTRUCTION over Q(q,t) ...")
    L = build_symbolic_matrices(ref_standard, samples, log)
    print("  reconstructed L_x1, L_x2, L_y1, L_y2 as 24x24 matrices "
          "over Q(q,t).")

    # collect denominator atoms appearing across all entries
    atoms = {}
    for g in GENERATORS:
        for k in range(N):
            for l in range(N):
                e = sp.cancel(L[g][k, l])
                _num, den = sp.fraction(sp.together(e))
                den = sp.factor(den)
                if den.is_number:
                    continue
                for fac in (den.args if den.func == sp.Mul else [den]):
                    base = fac
                    if base.func == sp.Pow:
                        base = base.args[0]
                    if base.is_number:
                        continue
                    key = sp.sstr(sp.factor(base))
                    atoms[key] = atoms.get(key, 0) + 1

    # ---- step 4: symbolic verification ---------------------------------
    print("\n[STEP 4] SYMBOLIC VERIFICATION of every relation over Q(q,t)")
    print("  convention: L_{a b} = L_a L_b (homomorphism); each 'expr = 0'")
    print("  relation is checked as the matrix identity  L_expr = 0_(24x24).")
    rels, S2_1, S2_2, S3_1 = relation_checks(L)
    all_ok = True
    results = []
    for name, mat in rels:
        ok = is_zero_matrix(mat)
        results.append((name, ok))
        all_ok &= ok
        print(f"  {name:34s} {'PASS' if ok else '*** FAIL ***'}")
        if not ok:
            res = mat.applyfunc(lambda e: sp.cancel(sp.together(e)))
            nz = [(i, j, res[i, j]) for i in range(res.shape[0])
                  for j in range(res.shape[1]) if res[i, j] != 0]
            print(f"      {len(nz)} nonzero residual entries; first few:")
            for (i, j, v) in nz[:5]:
                print(f"        [{i},{j}] = {v}")

    # sanity: X_2 = S2_1 operator is nonzero (the slider does not collapse)
    X2op = S2_1.applyfunc(lambda e: sp.cancel(sp.together(e)))
    x2_nonzero = X2op != sp.zeros(N, N)
    print(f"  {'X_2 = L_{S2_1} nonzero (slider survives)':34s} "
          f"{'PASS' if x2_nonzero else '*** FAIL ***'}")
    all_ok &= x2_nonzero

    # ---- step 5: all-word cyclic-vector check -------------------------
    print("\n[STEP 5] CYCLIC-VECTOR CHECK")
    cyclic_ok = check_cyclic_basis(L, ref_standard)
    print(f"  L_word e_unit = e_word for all {N} listed words: "
          f"{'PASS' if cyclic_ok else 'FAIL'}")
    all_ok &= cyclic_ok
    if cyclic_ok and all(ok for _name, ok in results):
        print(f"  Verified defining relations and {N} independent orbit vectors give dim A >= {N}.")
    else:
        print("  These checks do not establish a lower dimension bound.")
    print("  A generic upper dimension bound still requires independent evidence;")
    print("  cyclicity alone does not identify this as the regular representation.")

    # ---- denominator atoms report -------------------------------------
    print("\n[DENOMINATOR ATOMS] (irreducible factors in entry denominators)")
    if atoms:
        for key in sorted(atoms, key=lambda s: (-atoms[s], s)):
            print(f"  {key:28s}  appears in {atoms[key]} entries")
    else:
        print("  (none -- all entries are polynomial in q, t)")
    print("  These define the bad locus where the regular rep degenerates.")

    # ---- step 6: save artifacts ---------------------------------------
    print("\n[STEP 6] SAVING ARTIFACTS")
    prepare_output_dir(OUT_DIR)
    payload = {
        "description": "Candidate regular-representation matrices for Z_3^square-zero over Q(q,t)",
        "basis_complete": False,
        "convention": "L_{ab} = L_a L_b; L_g[k][l] = coef of m_k in g*m_l",
        "standard_monomials": ref_standard,
        "symbols": {"q": "q", "t": "t"},
        "matrices": {g: L[g] for g in GENERATORS},
        "denominator_atoms": atoms,
        "relation_results": results,
        "cyclic_basis_verified": bool(cyclic_ok),
        "all_pass": bool(all_ok),
    }
    save_pickle_output(OUT_PKL, payload, producer=producer)
    print(f"  pickle -> {OUT_PKL}")

    # ---- final verdict -------------------------------------------------
    print("\n" + "=" * 74)
    n_pass = sum(1 for _n, ok in results if ok)
    print(f"RELATIONS: {n_pass}/{len(results)} passed symbolically over Q(q,t)")
    print(f"OVERALL: {'ALL CHECKS PASS' if all_ok else '*** SOME CHECK FAILED ***'}")
    print("=" * 74)
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
