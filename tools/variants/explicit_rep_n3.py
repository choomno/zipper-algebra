#!/usr/bin/env python3
"""Explicit 3-dimensional representation of the n=3 zipper-variant algebra
over the rational function field Q(q, lambda).

This certifies, self-containedly, that in the canonical lambda-untwist
quotient C_lambda = Z_3 / (canonical tower at framing lambda):

  * the slider X_2 = q*y_1 + (1-q) - x_1 is NONZERO (rank 1 here), and
  * sigma_1 has eigenvalues exactly {1, -q, lambda}.

Relation set (cross-checked against build_relations in
tools/zipper_in_ms.py and singular_run/square-zero-n3-symbolic.sing;
notation x_j = sigma_j, y_j = sigma_j^{-1}, S2_j = q*y_j + (1-q) - x_j):

  R1 : x_j y_j = y_j x_j = 1                          (j = 1, 2)
  R2 : x1 x2 x1 = x2 x1 x2  and the y-version
  R4 : x1 x2 y1 = y2 x1 x2,  y1 y2 x1 = x2 y1 y2,
       x1 y2 y1 = y2 y1 x2,  y1 x2 x1 = x2 x1 y2      (group identities)
  R5 : (q y2 + (1-q) - x2 - q y1 y2 + x1 x2) * S2_1 = 0
  R7 : (x_j - lambda) * S2_j = 0                       (j = 1, 2; canonical)
  R8 : (x2 x1 - q^2 y2 y1) * S3_1 = 0,
       where S3_1 = (q^2 y1 y2 - x1 x2) * S2_1

NOTE on framing powers: verifying the canonical untwist
(x_j - lambda) * S2_j = 0 suffices for every framing power a
simultaneously, because x_j^a - t = (telescoping sum) * (x_j - lambda)
when t = lambda^a.  Hence we work over Q(q, lambda) with two independent
transcendentals and never introduce t or a.

The representation (Tuba-Wenzl triangular gauge; cf. Tuba & Wenzl,
Pacific J. Math 197 (2001) -- the unique 3-dim B_3 pair with spectrum
{1, -q, lambda}, found here directly from the braid equations):

      sigma_1 = [ 1  -(lam+q^2)   q  ]      sigma_2 = [ lam      0        0 ]
                [ 0     -q        1  ]                [  1      -q        0 ]
                [ 0      0       lam ]                [  q  -(lam+q^2)    1 ]

A second gauge with sigma_1 = diag(1, -q, lam) (and sigma_2 = P D P^{-1},
P explicit) is also constructed and fully re-verified; there
X_2 = diag(0, 0, s) with s = q/lam + 1 - q - lam = -(lam-1)(lam+q)/lam.

Usage:
  python3 explicit_rep_n3.py              # exact symbolic verification
  python3 explicit_rep_n3.py --numeric    # float sanity check, q=5, lam=7^(1/4)
  python3 explicit_rep_n3.py --all        # both

Exit status 0 iff every check passes.
"""
from __future__ import annotations

import argparse
import sys

import sympy as sp

q, lam = sp.symbols('q lam', nonzero=True)
I3 = sp.eye(3)


# --------------------------------------------------------------------------
# The representation
# --------------------------------------------------------------------------
def tw_matrices():
    """sigma_1, sigma_2 in the Tuba-Wenzl triangular gauge."""
    M1 = sp.Matrix([[1, -(lam + q**2), q],
                    [0, -q,            1],
                    [0,  0,          lam]])
    M2 = sp.Matrix([[lam,  0,            0],
                    [1,   -q,            0],
                    [q,   -(lam + q**2), 1]])
    return M1, M2


def diag_gauge_matrices():
    """Same rep conjugated so sigma_1 = diag(1, -q, lam); returns
    (M1d, M2d, U) with M1d = U^{-1} M1 U, M2d = U^{-1} M2 U."""
    M1, M2 = tw_matrices()
    # Columns of U: eigenvectors of the triangular M1 for 1, -q, lam.
    U = sp.Matrix([[1, lam + q**2, lam*(1 - q)/(1 - lam)],
                   [0, 1 + q,      1],
                   [0, 0,          q + lam]])
    Ui = U.inv()
    M1d = (Ui * M1 * U).applyfunc(sp.cancel)
    M2d = (Ui * M2 * U).applyfunc(lambda e: sp.cancel(sp.together(e)))
    return M1d, M2d, U


# --------------------------------------------------------------------------
# Relation checks (symbolic)
# --------------------------------------------------------------------------
def is_zero_matrix(M):
    return M.applyfunc(lambda e: sp.cancel(sp.together(e))) == sp.zeros(*M.shape)


def relation_checks(M1, M2, qv=q, lamv=lam):
    """Return ordered list of (name, matrix-that-should-be-zero).

    ``qv``/``lamv`` default to the symbols; pass numbers for numeric mode.
    """
    Y1 = M1.inv().applyfunc(sp.cancel)
    Y2 = M2.inv().applyfunc(sp.cancel)
    S2_1 = (qv*Y1 + (1 - qv)*I3 - M1).applyfunc(sp.cancel)    # the slider X_2
    S2_2 = (qv*Y2 + (1 - qv)*I3 - M2).applyfunc(sp.cancel)
    S3_1 = ((qv**2*Y1*Y2 - M1*M2) * S2_1).applyfunc(sp.cancel)

    rels = [
        ('R1  x1*y1 - 1',                M1*Y1 - I3),
        ('R1  y1*x1 - 1',                Y1*M1 - I3),
        ('R1  x2*y2 - 1',                M2*Y2 - I3),
        ('R1  y2*x2 - 1',                Y2*M2 - I3),
        ('R2  x1x2x1 - x2x1x2',          M1*M2*M1 - M2*M1*M2),
        ('R2  y1y2y1 - y2y1y2',          Y1*Y2*Y1 - Y2*Y1*Y2),
        ('R4a x1x2y1 - y2x1x2',          M1*M2*Y1 - Y2*M1*M2),
        ('R4b y1y2x1 - x2y1y2',          Y1*Y2*M1 - M2*Y1*Y2),
        ('R4c x1y2y1 - y2y1x2',          M1*Y2*Y1 - Y2*Y1*M2),
        ('R4d y1x2x1 - x2x1y2',          Y1*M2*M1 - M2*M1*Y2),
        ('R5  zipper level-2 (j=1)',
         (qv*Y2 + (1 - qv)*I3 - M2 - qv*Y1*Y2 + M1*M2) * S2_1),
        ('R7  (x1 - lam)*S2_1',          (M1 - lamv*I3) * S2_1),
        ('R7  (x2 - lam)*S2_2',          (M2 - lamv*I3) * S2_2),
        ('R8  (x2x1 - q^2 y2y1)*S3_1',   (M2*M1 - qv**2*Y2*Y1) * S3_1),
        ('cubic minpoly of x1: (x1-1)(x1+q)(x1-lam)',
         (M1 - I3) * (M1 + qv*I3) * (M1 - lamv*I3)),
        ('cubic minpoly of x2: (x2-1)(x2+q)(x2-lam)',
         (M2 - I3) * (M2 + qv*I3) * (M2 - lamv*I3)),
    ]
    return rels, S2_1, S2_2


def run_symbolic():
    all_ok = True
    s = sp.cancel(q/lam + 1 - q - lam)   # scalar S2(lam); the X_2 "weight"

    M1d, M2d, _U = diag_gauge_matrices()
    for gauge_name, M1, M2 in [
            ('Tuba-Wenzl triangular gauge', *tw_matrices()),
            ('diagonal-sigma_1 gauge', M1d, M2d)]:
        print(f"\n=== {gauge_name} ===")
        print("sigma_1 =")
        sp.pprint(M1)
        print("sigma_2 =")
        sp.pprint(M2.applyfunc(sp.factor))

        # determinants (R1 solvability: invertibility over Q(q, lam))
        d1, d2 = sp.factor(M1.det()), sp.factor(M2.det())
        print(f"det sigma_1 = {d1},  det sigma_2 = {d2}  "
              "(nonzero in Q(q,lam): generators invertible)")

        rels, X2, S2_2 = relation_checks(M1, M2)
        for name, mat in rels:
            ok = is_zero_matrix(mat)
            all_ok &= ok
            print(f"  {name:46s} {'PASS' if ok else 'FAIL'}")
            if not ok:
                print("    residual:")
                sp.pprint(mat.applyfunc(sp.cancel))

        # eigenvalues
        ev1, ev2 = M1.eigenvals(), M2.eigenvals()
        want = {sp.Integer(1): 1, -q: 1, lam: 1}
        ok_ev = (ev1 == want and ev2 == want)
        all_ok &= ok_ev
        print(f"  {'eigenvalues of sigma_1, sigma_2 = {1, -q, lam}':46s} "
              f"{'PASS' if ok_ev else 'FAIL'}")
        print(f"    eig(sigma_1) = {ev1}")
        print(f"    eig(sigma_2) = {ev2}")
        # distinct for generic q, lam + cubic minpoly above => diagonalizable,
        # i.e. conjugate to diag(1, -q, lam).

        # slider X_2 nonzero, rank 1
        X2c = X2.applyfunc(sp.cancel)
        rank = X2c.rank()
        ok_x2 = (X2c != sp.zeros(3, 3)) and rank == 1
        all_ok &= ok_x2
        print(f"  {'X_2 = q*y1 + (1-q) - x1  nonzero, rank 1':46s} "
              f"{'PASS' if ok_x2 else 'FAIL'}")
        print(f"    rank X_2 = {rank};  X_2 =")
        sp.pprint(X2c.applyfunc(sp.factor))

    print(f"\nscalar S2(mu) = q/mu + 1 - q - mu  at mu = lam:  "
          f"s = {sp.factor(s)}  (nonzero in Q(q,lam))")
    print("   S2(1) = 0, S2(-q) = 0  => X_2 is rank 1, supported on the "
          "lam-eigenline of sigma_1.")
    return all_ok


# --------------------------------------------------------------------------
# Numeric sanity check
# --------------------------------------------------------------------------
def run_numeric(tol=1e-9):
    qv = sp.Float(5, 30)
    lv = sp.Float(7, 30)**sp.Rational(1, 4)   # 7^(1/4) = 1.62657656169...
    subs = {q: qv, lam: lv}
    M1, M2 = tw_matrices()
    M1n = M1.subs(subs).evalf(30)
    M2n = M2.subs(subs).evalf(30)
    rels, X2, _ = relation_checks(M1n, M2n, qv=qv, lamv=lv)

    all_ok = True
    print("\n=== numeric sanity check:  q = 5, lambda = 7**(1/4) ===")
    for name, mat in rels:
        res = max(abs(complex(e)) for e in mat.evalf(30))
        ok = res < tol
        all_ok &= ok
        print(f"  {name:46s} residual {res:.3e}  {'PASS' if ok else 'FAIL'}")

    for tag, M in [('sigma_1', M1n), ('sigma_2', M2n)]:
        cp = M.charpoly(sp.Symbol('mu'))
        roots = sorted((complex(r) for r in sp.nroots(cp.as_expr())),
                       key=lambda zz: zz.real)
        want = sorted([complex(-5.0), complex(1.0), complex(7**0.25)],
                      key=lambda zz: zz.real)
        err = max(abs(rr - ww) for rr, ww in zip(roots, want))
        ok = err < tol
        all_ok &= ok
        print(f"  eig({tag}) = "
              f"{[f'{rr.real:.6f}{rr.imag:+.1e}j' for rr in roots]}  "
              f"(want ~ [-5, 1, 1.626577])  err {err:.3e}  "
              f"{'PASS' if ok else 'FAIL'}")

    x2max = max(abs(complex(e)) for e in X2.evalf(30))
    ok = x2max > 1e-6
    all_ok &= ok
    print(f"  {'X_2 nonzero numerically':46s} max|entry| {x2max:.6f}  "
          f"{'PASS' if ok else 'FAIL'}")
    return all_ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--numeric', action='store_true',
                    help='run only the float sanity check')
    ap.add_argument('--all', action='store_true',
                    help='run symbolic and numeric checks')
    args = ap.parse_args()

    ok = True
    if args.all or not args.numeric:
        ok &= run_symbolic()
    if args.all or args.numeric:
        ok &= run_numeric()

    print('\nOVERALL:', 'ALL CHECKS PASS' if ok else '*** SOME CHECK FAILED ***')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
