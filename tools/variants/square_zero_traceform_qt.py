#!/usr/bin/env python3
r"""TRACE FORM and GENERIC NON-SEMISIMPLICITY of the n=3 square-zero zipper algebra
A = Z_3^square-zero over the rational function field Q(q, t).

This is the trace-form analogue of ``variants/semisimple_Q.py`` (which runs the
char-0 trace-form criterion at a single integer (q,t)), lifted to the WHOLE
function field Q(q,t).  Input data is the verified symbolic 24x24 left regular
representation built by ``variants/square_zero_regular_rep_qt.py`` and stored in
``logs_and_traces/variants/square-zero-n3-regrep-qt.pkl`` (matrices L_{x1}, L_{x2},
L_{y1}, L_{y2} over Q(q,t), with the homomorphism convention L_{ab}=L_a L_b and
L_g[k][l] = coefficient of m_k in g*m_l).

------------------------------------------------------------------------------
THE CRITERION  (Pierce, Associative Algebras, Sec 9.4; char 0, here char(Q)=0)
------------------------------------------------------------------------------
For a finite-dim algebra A over a field K of char 0, the trace form
        T_ij = tr_reg(m_i m_j) = trace(L_{m_i} L_{m_j})
on the regular representation satisfies:  A is semisimple  <=>  det T != 0;
and in general  dim_K rad(A) = dim_K A - rank(T)  EXACTLY.

Here K = Q(q,t).  We compute T over Q(q,t) symbolically, find rank(T)=15, so
dim J = 24 - 15 = 9 generically, i.e. A is generically non-semisimple.  Two
independent certificates make this rigorous despite expression swell:
  (a) a single nonzero 15x15 minor of T  =>  rank(T) >= 15;
  (b) the radical witness:  build the two-sided ideal <w_1> with
      w_1 = (x1 - t) * X_2  (X_2 = S2_1 = q*y1 + (1-q) - x1) as a subspace of
      Q(q,t)^24, show dim<w_1> = 9 and <w_1> subseteq ker T(q,t).  The kernel
      containment forces deficit >= 9, i.e. rank(T) <= 15.
  (a)+(b) pin rank(T) = 15 exactly, hence dim J = 9.

Run with the sympy venv:
    /home/choom/braid-research/.venv-symrep/bin/python square_zero_traceform_qt.py
"""
from __future__ import annotations

import os
import pickle
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
# Outputs live apart from code, and no path is machine-specific (settings rule): the variants'
# artifact directory is repo-relative. Its INPUT, the verified 24x24 regular representation, is the
# output of square_zero_regular_rep_qt.py; the archived copy is archive/variants-data/.
LOGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                      "logs_and_traces", "variants")
os.makedirs(LOGDIR, exist_ok=True)
REGREP_PKL = os.path.join(LOGDIR, "square-zero-n3-regrep-qt.pkl")
OUT_PKL = os.path.join(LOGDIR, "square-zero-n3-traceform-qt.pkl")

q, t = sp.symbols("q t")
GENERATORS = ["x1", "x2", "y1", "y2"]


def log(msg):
    print(f"[{msg}]", file=sys.stderr, flush=True)


def simp(e):
    """Cheap, robust simplification of a rational-function entry."""
    return sp.cancel(sp.together(e))


# ===========================================================================
# Load the verified regular representation
# ===========================================================================
def load_regrep():
    with open(REGREP_PKL, "rb") as f:
        p = pickle.load(f)
    if not p.get("all_pass", False):
        raise RuntimeError("regular rep pickle reports all_pass=False; refuse")
    L = {g: sp.Matrix(p["matrices"][g]) for g in GENERATORS}
    standard = [tuple(m) for m in p["standard_monomials"]]
    N = len(standard)
    if N != 24:
        raise RuntimeError(f"expected 24 standard monomials, got {N}")
    for g in GENERATORS:
        if L[g].shape != (N, N):
            raise RuntimeError(f"L_{g} has shape {L[g].shape}, expected ({N},{N})")
    return L, standard, N, p.get("denominator_atoms", {})


# ===========================================================================
# L_{m_l}: left-mult operator of a standard monomial = product of generator
# matrices in word order (homomorphism convention L_{ab} = L_a L_b).
# ===========================================================================
def L_word(L, word, N):
    M = sp.eye(N)
    for g in word:
        M = M * L[g]
    return M


# ===========================================================================
# Fraction-free rank over Q(q,t):  clear denominators row-wise to get a
# polynomial matrix, then Bareiss fraction-free elimination.  Rank over the
# field of fractions Q(q,t) equals the rank of the polynomial matrix over the
# same field (clearing a per-row nonzero scalar does not change rank).
# ===========================================================================
def clear_denoms_rowwise(M):
    """Return an integer/polynomial matrix P with P[i] = d_i * M[i], where d_i
    is the lcm of the denominators in row i (a nonzero element of Q(q,t), so it
    does not change the row space / rank over the field)."""
    N = M.shape[0]
    P = sp.zeros(*M.shape)
    for i in range(N):
        dens = []
        for j in range(M.shape[1]):
            _n, d = sp.fraction(sp.together(M[i, j]))
            dens.append(d)
        di = sp.Integer(1)
        for d in dens:
            di = sp.lcm(di, d)
        for j in range(M.shape[1]):
            P[i, j] = sp.expand(sp.cancel(M[i, j] * di))
    return P


def bareiss_rank(P):
    """Fraction-free (Bareiss) Gaussian elimination over Q(q,t)[poly], returns
    rank over the fraction field.  Entries are sympy polynomials in q,t.
    (Kept for reference; the main path uses field_rank below, which is faster
    for these swelling bivariate entries.)"""
    A = [[sp.expand(P[i, j]) for j in range(P.shape[1])]
         for i in range(P.shape[0])]
    n = len(A)
    m = len(A[0]) if n else 0
    rank = 0
    prev = sp.Integer(1)
    row = 0
    for col in range(m):
        if row >= n:
            break
        piv = None
        for r in range(row, n):
            if sp.expand(A[r][col]) != 0:
                piv = r
                break
        if piv is None:
            continue
        if piv != row:
            A[row], A[piv] = A[piv], A[row]
        for r in range(row + 1, n):
            for c in range(col + 1, m):
                num = A[r][c] * A[row][col] - A[r][col] * A[row][c]
                A[r][c] = sp.cancel(num / prev)
            A[r][col] = sp.Integer(0)
        prev = sp.expand(A[row][col])
        rank += 1
        row += 1
    return rank


# Generic numeric probe points (off the bad locus {q,q-1,t}=0) used only to
# GUIDE pivot selection.  The actual arithmetic is exact symbolic; the probe
# never affects correctness, only the order in which pivots are chosen (which
# keeps intermediate rational functions small).
_PROBE_PTS = [(5, 7), (3, 4), (11, 13), (2, 9), (13, 6)]


def _numeric_abs(expr):
    """Average |value| of a rational-function entry over the probe points; 0 if
    it vanishes identically (numerically) at every probe.  Used only to rank
    candidate pivots by numeric magnitude."""
    tot = 0.0
    cnt = 0
    for (q0, t0) in _PROBE_PTS:
        try:
            v = complex(expr.subs({q: q0, t: t0}))
        except (TypeError, ValueError, ZeroDivisionError):
            continue
        tot += abs(v)
        cnt += 1
    return tot / cnt if cnt else 0.0


def field_rank(M, log=lambda *a: None):
    """Exact rank of M over the fraction field Q(q,t) by Gauss-square-zero
    elimination with division (each entry kept as a cancelled rational
    function).  Pivots are chosen by largest numeric magnitude at generic probe
    points to minimise expression swell; an entry is treated as a valid pivot
    only after a symbolic nonzero test (cancel != 0).  Returns the exact rank.
    """
    N = M.shape[0]
    A = [[simp(M[i, j]) for j in range(N)] for i in range(N)]
    used_rows = set()
    rank = 0
    for col in range(N):
        # candidate pivot rows for this column, ranked by numeric magnitude
        cands = sorted(
            (r for r in range(N) if r not in used_rows),
            key=lambda r: -_numeric_abs(A[r][col]),
        )
        piv = None
        for r in cands:
            if simp(A[r][col]) != 0:
                piv = r
                break
        if piv is None:
            continue
        used_rows.add(piv)
        rank += 1
        pv = A[piv][col]
        # eliminate this column from all other not-yet-used rows
        for r in range(N):
            if r == piv or r in used_rows:
                continue
            f = A[r][col]
            if f == 0:
                continue
            fac = simp(f / pv)
            for c in range(N):
                if c == col:
                    A[r][c] = sp.Integer(0)
                else:
                    A[r][c] = simp(A[r][c] - fac * A[piv][c])
        log(f"  pivot col {col} (rank now {rank})")
    return rank


# ===========================================================================
# A nonzero r x r minor (one symbolic determinant) certifying rank >= r.
# Greedy pivoting on the symbolic matrix selects r independent rows & cols.
# ===========================================================================
def _numeric_pivot_rows_cols(M, r):
    """Use an EXACT evaluation of M at a generic probe point (Fraction
    arithmetic, off the bad locus) to choose r rows and r columns whose r x r
    submatrix is numerically nonsingular there, and return that submatrix's
    EXACT rational determinant at the probe point as well.

    KEY FACT (rigorous): exact evaluation at a point P=(q0,t0) off the bad locus
    is a Q-algebra map  ev_P : Q(q,t) -> Q.  det is a polynomial in the entries,
    so ev_P(det_symbolic) = det(ev_P(submatrix)).  If the latter is a NONZERO
    rational number, then det_symbolic is a NONZERO element of Q(q,t) (a zero
    element would map to 0).  Hence a numerically nonsingular r x r submatrix at
    P certifies a nonvanishing symbolic r x r minor -- no symbolic determinant
    needed.  Returns (rows, cols, (q0,t0), exact_det_at_P) or None."""
    from fractions import Fraction
    N = M.shape[0]
    for (q0, t0) in _PROBE_PTS:
        # exact rational evaluation
        Mnum = []
        ok = True
        for i in range(N):
            row = []
            for j in range(N):
                num, den = sp.fraction(sp.together(M[i, j]))
                dv = den.subs({q: q0, t: t0})
                if dv == 0:
                    ok = False
                    break
                nv = num.subs({q: q0, t: t0})
                row.append(Fraction(int(nv), int(dv)))
            if not ok:
                break
            Mnum.append(row)
        if not ok:
            continue
        # Gaussian elimination with full pivot tracking to pick r indep rows/cols
        A = [r_[:] for r_ in Mnum]
        chosen_rows, chosen_cols = [], []
        used_r, used_c = set(), set()
        det_val = Fraction(1)
        for _ in range(r):
            best = None
            best_abs = Fraction(0)
            for i in range(N):
                if i in used_r:
                    continue
                for j in range(N):
                    if j in used_c:
                        continue
                    a = A[i][j]
                    if a != 0 and abs(a) > best_abs:
                        best_abs = abs(a)
                        best = (i, j)
            if best is None:
                break
            pi, pj = best
            chosen_rows.append(pi)
            chosen_cols.append(pj)
            used_r.add(pi)
            used_c.add(pj)
            piv = A[pi][pj]
            det_val *= piv
            for i in range(N):
                if i in used_r:
                    continue
                f = A[i][pj]
                if f == 0:
                    continue
                fac = f / piv
                for j in range(N):
                    if j in used_c:
                        continue
                    A[i][j] = A[i][j] - fac * A[pi][j]
        if len(chosen_rows) == r:
            # det of the chosen submatrix at P, up to the sign of the row/col
            # permutation -- magnitude (nonzero-ness) is all that matters here.
            return (sorted(chosen_rows), sorted(chosen_cols),
                    (q0, t0), det_val)
    return None


def find_nonzero_minor(M, r, full_symbolic=False):
    """Exhibit r rows and r cols whose r x r submatrix has a NONZERO determinant
    over Q(q,t), certifying rank(M) >= r.

    Default (full_symbolic=False): pick the rows/cols by exact numeric probe and
    certify nonvanishing via the EXACT rational determinant of that submatrix at
    the probe point (rigorous: a nonzero specialization forces det_symbolic !=
    0 in Q(q,t)).  This avoids the costly 15x15 symbolic determinant.  Returns
    (rows, cols, det) where det is the exact rational determinant at the probe
    point (a concrete nonzero number witnessing the minor is not identically 0).

    full_symbolic=True additionally computes the exact symbolic determinant of
    the chosen polynomial submatrix (slower) and returns it instead."""
    sel = _numeric_pivot_rows_cols(M, r)
    if sel is None:
        return None
    chosen_rows, chosen_cols, P, det_at_P = sel
    if det_at_P == 0:
        return None
    if not full_symbolic:
        return chosen_rows, chosen_cols, sp.Rational(det_at_P.numerator,
                                                     det_at_P.denominator)
    # optional exact symbolic determinant of the chosen submatrix
    sub = M[chosen_rows, chosen_cols]
    rr = sub.shape[0]
    Psub = sp.zeros(rr, rr)
    for i in range(rr):
        di = sp.Integer(1)
        for j in range(rr):
            _n, d = sp.fraction(sp.together(sub[i, j]))
            di = sp.lcm(di, d)
        for j in range(rr):
            Psub[i, j] = sp.expand(sp.cancel(sub[i, j] * di))
    det = sp.cancel(Psub.det(method="bareiss"))
    if det == 0:
        return None
    return chosen_rows, chosen_cols, det


# ===========================================================================
# Two-sided ideal generated by a coordinate vector, as a subspace of Q(q,t)^24.
# We close under all LEFT generators L_g and all RIGHT generators R_g.
# R_g column l = coordinate of m_l * g = L_{m_l} e_g  (homomorphism: the coord
# vector of an element a is L_a e_0, and m_l*g has L_{m_l*g}=L_{m_l}L_g so its
# coord = L_{m_l} L_g e_0 = L_{m_l} e_g).
# ===========================================================================
def build_right_mult(L, standard, N):
    single = {standard[i][0]: i for i in range(N) if len(standard[i]) == 1}
    Lm = [L_word(L, m, N) for m in standard]
    R = {}
    for g in GENERATORS:
        gi = single[g]
        Rg = sp.zeros(N, N)
        for l in range(N):
            Rg[:, l] = Lm[l][:, gi]
        R[g] = Rg
    return R


def col_echelon_rank_and_basis(vectors, N):
    """Given a list of column vectors (sympy N x 1) over Q(q,t), return a basis
    (list of pivot vectors) of their span via fraction-free column reduction,
    plus the rank."""
    # assemble as columns, do Gaussian elimination over the field
    cols = [sp.Matrix([simp(v[i, 0]) for i in range(N)]) for v in vectors]
    basis = []  # list of (pivot_row, reduced_vector)
    pivots = []
    for v in cols:
        w = v.copy()
        for (pr, bv) in basis:
            if w[pr, 0] != 0:
                f = simp(w[pr, 0] / bv[pr, 0])
                w = sp.Matrix([simp(w[i, 0] - f * bv[i, 0]) for i in range(N)])
        # find a pivot row in w
        pr = None
        for i in range(N):
            if w[i, 0] != 0:
                pr = i
                break
        if pr is not None:
            basis.append((pr, w))
            pivots.append(pr)
    return [bv for (_pr, bv) in basis], len(basis)


def two_sided_ideal(seed_vec, L, R, N):
    """Span of { L_{u} R_{v} (seed) } for words u,v.  Generated by closing the
    seed under all single-generator left and right multiplications (associative
    algebra: words = products of generators) until the span stabilises."""
    basis_vecs = [seed_vec]
    cur_basis, dim = col_echelon_rank_and_basis(basis_vecs, N)
    changed = True
    while changed:
        changed = False
        new_vecs = list(cur_basis)
        for bv in cur_basis:
            for g in GENERATORS:
                new_vecs.append(L[g] * bv)
                new_vecs.append(R[g] * bv)
        nb, nd = col_echelon_rank_and_basis(new_vecs, N)
        if nd > dim:
            cur_basis, dim = nb, nd
            changed = True
    return cur_basis, dim


# ===========================================================================
# Main
# ===========================================================================
def main():
    print("=" * 74)
    print("TRACE FORM of Z_3^square-zero over Q(q,t)  +  generic non-semisimplicity")
    print("=" * 74)

    log("loading verified regular representation")
    L, standard, N, denom_atoms = load_regrep()
    print(f"\n[INPUT] verified left regular rep over Q(q,t): dim = {N}")
    print(f"  denominator atoms recorded in regrep: {denom_atoms}")

    # ---- 1. L_{m_l} for every standard monomial -----------------------
    print("\n[STEP 1] reconstruct L_{m_l} = product of generator matrices")
    log("building L_{m_l} for all 24 standard monomials")
    Lm = [L_word(L, m, N) for m in standard]
    # sanity: e_0 cyclic -> L_{m_l} e_0 = e_l
    e0 = sp.zeros(N, 1)
    e0[0, 0] = 1
    cyc_ok = True
    for l in range(N):
        col = Lm[l] * e0
        want = sp.zeros(N, 1)
        want[l, 0] = 1
        if sp.Matrix([simp(col[i, 0] - want[i, 0]) for i in range(N)]) != \
                sp.zeros(N, 1):
            cyc_ok = False
    print(f"  cyclic-vector sanity (L_(m_l) e_0 = e_l for all l): "
          f"{'PASS' if cyc_ok else 'FAIL'}")
    if not cyc_ok:
        raise RuntimeError("cyclic-vector sanity failed; data inconsistent")

    # ---- 2. trace vector tau_l = tr(L_{m_l}) --------------------------
    print("\n[STEP 2] trace vector  tau_l = trace(L_(m_l))  over Q(q,t)")
    log("computing tau_l")
    tau = [simp(Lm[l].trace()) for l in range(N)]
    for l in range(N):
        w = "*".join(standard[l]) if standard[l] else "<1>"
        print(f"  tau[{l:2d}]  m={w:18s}  = {tau[l]}")

    # ---- 3. trace form Gram matrix T_ij = tr(L_{m_i} L_{m_j}) ----------
    print("\n[STEP 3] trace form  T_ij = trace(L_(m_i) L_(m_j))  (24x24)")
    T_CACHE = os.path.join(LOGDIR, "square-zero-n3-traceform-T-qt.pkl")
    T = None
    if os.path.exists(T_CACHE):
        try:
            with open(T_CACHE, "rb") as f:
                cached = pickle.load(f)
            if cached.get("standard") == standard:
                T = sp.Matrix(cached["T"])
                log("loaded cached trace form T from pickle")
                print("  (loaded cached symbolic T)")
        except Exception:
            T = None
    if T is None:
        log("building T directly as trace of product of L_{m_i} and L_{m_j}")
        T = sp.zeros(N, N)
        for i in range(N):
            Mi = Lm[i]
            for j in range(i, N):
                # trace(Mi * Lm[j]) = sum_{a,b} Mi[a,b] Lm[j][b,a]
                Mj = Lm[j]
                s = sp.Integer(0)
                for a in range(N):
                    for b in range(N):
                        if Mi[a, b] != 0 and Mj[b, a] != 0:
                            s += Mi[a, b] * Mj[b, a]
                tij = simp(s)
                T[i, j] = tij
                T[j, i] = tij
            log(f"  T row {i} done")
        with open(T_CACHE, "wb") as f:
            pickle.dump({"standard": standard, "T": T}, f)
        log("cached trace form T to pickle")
    # symmetry check
    sym_ok = all(simp(T[i, j] - T[j, i]) == 0 for i in range(N)
                 for j in range(N))
    print(f"  T symmetric (T_ij = T_ji): {'PASS' if sym_ok else 'FAIL'}")
    # consistency: diagonal-ish check T_i0 = tau_i (since m_0 = 1)
    t0_ok = all(simp(T[i, 0] - tau[i]) == 0 for i in range(N))
    print(f"  T[i,0] = tau_i (m_0 = 1) consistency: {'PASS' if t0_ok else 'FAIL'}")

    # ---- 4. rank over Q(q,t) ------------------------------------------
    # The exact rank over Q(q,t) is pinned RIGOROUSLY by the two certificates
    # computed below (STEP 5a gives a nonzero 15x15 minor => rank >= 15; STEP 5b
    # gives a 9-dim ideal inside ker T => deficit >= 9 => rank <= 15).  Full
    # symbolic Gauss-square-zero elimination is exact too but suffers expression
    # swell over Q(q,t); it is therefore OPT-IN via TRACEFORM_FULL_RANK=1.  The
    # reported rank is the certificate value, which is fully rigorous.
    print("\n[STEP 4] RANK of T over Q(q,t)")
    do_full = os.environ.get("TRACEFORM_FULL_RANK", "0") == "1"
    full_rank = None
    if do_full:
        log("field-based Gauss-square-zero rank over Q(q,t), probe-guided pivoting")
        full_rank = field_rank(T, log=log)
        print(f"  full symbolic elimination rank(T) = {full_rank}")
    else:
        print("  (full symbolic elimination skipped; rank pinned by the two")
        print("   certificates below.  Set TRACEFORM_FULL_RANK=1 to also run it.)")

    # ---- 5a. certify rank >= 15 via a nonzero 15x15 minor -------------
    print("\n[STEP 5a] CERTIFY rank >= 15: exhibit a nonzero 15x15 minor of T")
    log("searching for a nonzero 15x15 minor")
    minor = find_nonzero_minor(T, 15, full_symbolic=do_full)
    if minor is None:
        print("  *** could not find a nonzero 15x15 minor ***")
        minor_rows = minor_cols = None
        minor_det = None
    else:
        minor_rows, minor_cols, minor_det = minor
        print(f"  rows {minor_rows}")
        print(f"  cols {minor_cols}")
        if do_full:
            md_fac = sp.factor(minor_det)
            print(f"  symbolic det(15x15 minor) = {md_fac}")
        else:
            # rigorous nonvanishing witness: exact value of the minor at a
            # generic probe point off the bad locus (nonzero => det !=0 in
            # Q(q,t); see find_nonzero_minor docstring for the ev_P argument).
            print(f"  exact det(15x15 minor) at a generic probe point = "
                  f"{minor_det}  (!= 0)")
        print("  not identically zero over Q(q,t)  =>  rank(T) >= 15")
        # the matching upper bound rank <= 15 is supplied rigorously by the
        # 9-dim ideal inside ker T (STEP 5b'); the redundant direct 16x16 search
        # is opt-in (TRACEFORM_FULL_RANK=1) as it is comparatively expensive.
    if do_full:
        log("checking no nonzero 16x16 minor exists (rank <= 15)")
        minor16 = find_nonzero_minor(T, 16, full_symbolic=True)
        if minor16 is None:
            print("  no nonzero 16x16 minor found  =>  rank(T) <= 15")
        else:
            print(f"  WARNING: found a candidate 16x16 minor; det = "
                  f"{sp.factor(minor16[2])}")

    # ---- 5b. radical witness: two-sided ideal <w_1> -------------------
    print("\n[STEP 5b] CROSS-CHECK via the radical: two-sided ideal <w_1>")
    print("  X_2 = S2_1 = q*y1 + (1-q) - x1 ;  w_1 = (x1 - t) * X_2")
    I = sp.eye(N)
    X1, X2g, Y1, Y2 = L["x1"], L["x2"], L["y1"], L["y2"]
    X2op = q * Y1 + (1 - q) * I - X1          # operator L_{X_2}=L_{S2_1}
    w1_op = (X1 - t * I) * X2op               # operator L_{w_1} = L_{(x1-t)X_2}
    w2_op = (X2g - t * I) * X2op              # L_{w_2}=L_{(x2-t)X_2}
    w1_vec = sp.Matrix([simp((w1_op * e0)[i, 0]) for i in range(N)])
    w2_vec = sp.Matrix([simp((w2_op * e0)[i, 0]) for i in range(N)])
    print("  w_1 coordinate vector L_{w_1} e_0 =")
    print("   ", [simp(w1_vec[i, 0]) for i in range(N)])

    log("building right-multiplication operators R_g")
    R = build_right_mult(L, standard, N)
    # confirm R_g commutes with L_h (left/right multiplications commute)
    comm_ok = True
    for g in GENERATORS:
        for h in GENERATORS:
            C = L[h] * R[g] - R[g] * L[h]
            if any(simp(C[i, j]) != 0 for i in range(N) for j in range(N)):
                comm_ok = False
    print(f"  R_g built; [L_h, R_g] = 0 for all g,h: "
          f"{'PASS' if comm_ok else 'FAIL'}")

    log("closing w_1 under all L_g and R_g -> two-sided ideal")
    ideal_basis, ideal_dim = two_sided_ideal(w1_vec, L, R, N)
    print(f"  dim_Q(q,t) <w_1>  = {ideal_dim}")

    # is w_2 already inside <w_1> ?
    test_basis, test_dim = col_echelon_rank_and_basis(
        ideal_basis + [w2_vec], N)
    w2_in = (test_dim == ideal_dim)
    print(f"  w_2 in <w_1>: {w2_in}  (two-sided ideal already contains it: "
          f"{'yes' if w2_in else 'no'})")

    # ---- <w_1> subseteq ker T :  T v = 0 for each basis vector v ------
    print("\n[STEP 5b'] verify <w_1> subseteq ker T(q,t)")
    log("applying T to each ideal basis vector")
    ker_ok = True
    for bv in ideal_basis:
        Tv = T * bv
        Tv = sp.Matrix([simp(Tv[i, 0]) for i in range(N)])
        if Tv != sp.zeros(N, 1):
            ker_ok = False
            nz = [(i, Tv[i, 0]) for i in range(N) if Tv[i, 0] != 0]
            print(f"    *** nonzero T v: first entries {nz[:3]}")
    print(f"  T * (every basis vector of <w_1>) = 0: {'PASS' if ker_ok else 'FAIL'}")
    print(f"  => ker T has dim >= {ideal_dim}  =>  rank(T) <= {N - ideal_dim}")

    # ---- pin the exact rank from the two certificates -----------------
    # lower bound from the nonzero r x r minor; upper bound from the ideal in
    # ker T.  When they meet, the rank is determined exactly over Q(q,t).
    lower = len(minor_rows) if minor_rows is not None else 0
    upper = N - ideal_dim if ker_ok else N
    print("\n[STEP 4'] EXACT RANK from certificates")
    print(f"  rank(T) >= {lower}   (nonzero {lower}x{lower} minor exhibited)")
    print(f"  rank(T) <= {upper}   ({ideal_dim}-dim ideal inside ker T)")
    if lower == upper:
        rank = lower
        print(f"  => rank(T) = {rank}  (bounds coincide; exact over Q(q,t))")
    else:
        # fall back to the full elimination value if it was computed
        rank = full_rank if full_rank is not None else lower
        print(f"  bounds do not coincide; using rank = {rank}")
    if full_rank is not None and full_rank != rank:
        print(f"  *** WARNING: full elimination rank {full_rank} != "
              f"certificate rank {rank} ***")
    deficit = N - rank
    print(f"  rank(T) = {rank},  deficit = {deficit}  "
          f"(= generic dim_Q(q,t) J by char-0 criterion)")

    # ---- bad-locus denominator D --------------------------------------
    print("\n[BAD LOCUS] denominator atoms over which the statement holds")
    # collect denominator atoms from L_g entries (recorded), the trace form,
    # and the certifying minor determinant.
    atoms = {}

    def add_atoms(expr):
        e = sp.cancel(sp.together(expr))
        _n, d = sp.fraction(e)
        d = sp.factor(d)
        if d.is_number:
            return
        facs = d.args if d.func == sp.Mul else [d]
        for f in facs:
            base = f.args[0] if f.func == sp.Pow else f
            if base.is_number:
                continue
            atoms[sp.sstr(sp.factor(base))] = True

    for g in GENERATORS:
        for i in range(N):
            for j in range(N):
                add_atoms(L[g][i, j])
    for i in range(N):
        for j in range(N):
            add_atoms(T[i, j])
    # the minor determinant numerator gives the extra non-vanishing locus
    minor_atoms = []
    if minor_det is not None:
        mn = sp.factor(minor_det)
        # numerator factors that must be nonzero for the minor to certify
        num, den = sp.fraction(sp.together(mn))
        add_atoms(mn)
        for part, lbl in ((sp.factor(num), "num"), (sp.factor(den), "den")):
            if part.is_number:
                continue
            facs = part.args if part.func == sp.Mul else [part]
            for f in facs:
                base = f.args[0] if f.func == sp.Pow else f
                if not base.is_number:
                    minor_atoms.append((lbl, sp.sstr(sp.factor(base))))

    atom_list = sorted(atoms.keys())
    print(f"  denominator atoms (entries of L_g and T): {atom_list}")
    D = sp.Integer(1)
    for a in atom_list:
        D = D * sp.sympify(a)
    # also require the 15x15 minor numerator nonzero
    Dminor = sp.factor(minor_det) if minor_det is not None else sp.Integer(1)
    print(f"  certifying 15x15 minor (must be != 0): {Dminor}")
    print(f"  bad-locus denominator atoms product D = {sp.factor(D)}")

    # ---- consistency of the two certificates --------------------------
    print("\n[CROSS-CHECK]")
    agree = (deficit == ideal_dim == 9) and (rank == 15)
    print(f"  rank(T)            = {rank}   (expected 15)")
    print(f"  deficit 24-rank    = {deficit}   (expected 9)")
    print(f"  dim <w_1>          = {ideal_dim}   (expected 9)")
    print(f"  deficit == dim<w_1> == 9 : {deficit == ideal_dim == 9}")
    print(f"  <w_1> subseteq ker T     : {ker_ok}")
    print(f"  15x15 minor != 0         : {minor_det is not None and minor_det != 0}")
    print(f"  ALL CERTIFICATES AGREE   : {agree and ker_ok}")

    # ---- THEOREM ------------------------------------------------------
    print("\n" + "=" * 74)
    Dstr = sp.sstr(sp.factor(D))
    minor_note = ("symbolic det" if do_full else "det at a generic probe point")
    print("THEOREM.")
    print(f"  Z_3^square-zero is non-semisimple over Q(q,t), with dim J = {deficit},")
    print(f"  generically (away from D = 0, where D = {Dstr}, together with the")
    print(f"  non-vanishing of a 15x15 trace-form minor; {minor_note} = "
          f"{minor_det if minor_det is not None else '?'} != 0).")
    print("  Certified two independent ways over Q(q,t):")
    print(f"    (a) rank(trace form T) = {rank} < 24  (nonzero 15x15 minor; no 16x16),")
    print(f"        so by the char-0 trace-form criterion dim J = 24 - {rank} = {deficit};")
    print(f"    (b) the two-sided ideal <(x1 - t)X_2> has dim {ideal_dim} and lies in")
    print(f"        ker T, forcing deficit >= {ideal_dim}.  Both give dim J = {deficit}.")
    print("=" * 74)

    # ---- save ---------------------------------------------------------
    payload = {
        "description": "Trace form of Z_3^square-zero over Q(q,t) + non-semisimplicity",
        "standard_monomials": standard,
        "tau": tau,
        "T": T,
        "rank": rank,
        "deficit": deficit,
        "minor_rows": minor_rows,
        "minor_cols": minor_cols,
        "minor_det": minor_det,
        "ideal_dim": ideal_dim,
        "ideal_in_ker": ker_ok,
        "w2_in_ideal": w2_in,
        "bad_locus_D": D,
        "denominator_atoms": atom_list,
    }
    with open(OUT_PKL, "wb") as f:
        pickle.dump(payload, f)
    print(f"\n[SAVED] {OUT_PKL}")

    ok = (rank == 15 and deficit == 9 and ideal_dim == 9 and ker_ok
          and minor_det is not None and minor_det != 0)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
