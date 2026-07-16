#!/usr/bin/env python3
"""Derived k-strand untwisting = min poly of the lasso twist Theta_k on the slider
X_k in the two-strand-only quotient Z_n / (all zipper + p_2(sigma_j) X_2 = 0), for
general n. This is the finitization-program readout (see project-research-phase).

Findings (2026-07, cases p_2 in {sigma-t, (sigma-t)^2, sigma^2-t}, p* = reflection):
  n=3:  min poly of Theta_3 on X_3 = p*                       (sharpening drops all q-powers)
  n=4:  min poly of Theta_4 on X_4 = p*(L) * q^{2d} p*(L/q^2) (one q^2-dilate survives)
  tower conjecture:  g_k(L) = prod_{j=0}^{k-3} q^{2 j d} p*(L / q^{2 j}),  d = deg p*
                     spectrum(Theta_k on X_k) = U_{j=0}^{k-3} q^{2j} * roots(p*)

Method: per-point f4ncgb Groebner basis (mod p, or --char 0 for exact Q), then the min
poly of Theta_k on X_k by an incremental Krylov with a MEMOIZED normal-form reducer
(each distinct word reduced once; Theta applied to already-reduced vectors so words stay
short). The stock core.reduce_to_basis_* has no memoization and blows up past n=3.

Usage:
  python3 -m f4ncgb_tools.research.theta_min_poly --n 4 --cases linear jordan framing
  python3 -m f4ncgb_tools.research.theta_min_poly --n 4 --char 0 --grid 5,7 3,10
"""
import argparse, os, subprocess, sys
from fractions import Fraction as Fr
from f4ncgb_tools import core

F4 = os.environ.get("F4NCGB", "/home/choom/src/f4ncgb/build-release/f4ncgb")
HERE = os.path.dirname(os.path.abspath(__file__))
EMIT = os.path.normpath(f"{HERE}/../zipper_in_ms.py")
P31 = 2 ** 31 - 1
UNTW = {"linear": ("bigelow", "1"), "jordan": ("jordan", "1"), "framing": ("bigelow", "2")}
PSTAR = {  # p*(L) as ascending coeff list, symbolic in t (filled per point)
    "linear":  lambda t: [t, 1],
    "jordan":  lambda t: [t * t, 2 * t, 1],
    "framing": lambda t: [-t, 0, 1],
}


def slider(k, q, one, mul_scalar):
    """X_k terms as [(coef, word)] over the chosen ring. one/mul_scalar adapt Q vs F_p."""
    def qp(e):  # q**e in the ring
        r = one
        for _ in range(e): r = mul_scalar(r, q)
        return r
    X = [(q, ("y1",)), (mul_scalar(one, 1) if False else (one - q if isinstance(q, Fr) else (1 - q)), ()),
         (mul_scalar(one, -1), ("x1",))]
    for j in range(3, k + 1):
        ys = tuple(f"y{i}" for i in range(1, j)); xs = tuple(f"x{i}" for i in range(1, j))
        fac = [(qp(j - 1), ys), (mul_scalar(one, -1), xs)]
        X = ncmul(fac, X)
    return X


def ncmul(A, B):
    out = {}
    for ca, wa in A:
        for cb, wb in B:
            w = wa + wb; out[w] = out.get(w, 0) + ca * cb
    return [(c, w) for w, c in out.items() if c != 0]


def theta_word(k):
    pre = tuple(f"x{i}" for i in range(1, k))
    return pre + pre[::-1]


def make_reducer(polys, ring):
    """ring: 'Fp' (mod P31) or 'Q'. Returns (dim, nf, standard)."""
    lms = [lm for lm, _ in polys]; variables = core.collect_variables(polys)
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(delta, is_match, variables, max_total=2_000_000)
    mon_idx = {m: i for i, m in enumerate(standard)}
    if ring == "Fp":
        lm_to_tail = core.lm_to_tail_Fp(polys, P31); zero = 0
        norm = lambda c: c % P31
    else:
        lm_to_tail = core.lm_to_tail_Q(polys); zero = Fr(0)
        norm = lambda c: c
    memo = {}
    def nf(word):
        stack = [word]
        while stack:
            w = stack[-1]
            if w in memo: stack.pop(); continue
            idx = mon_idx.get(w)
            if idx is not None:
                memo[w] = {idx: (1 if ring == "Fp" else Fr(1))}; stack.pop(); continue
            m = core.find_leading_match(w, delta, is_match, match_pat)
            if m is None: raise RuntimeError(f"standard word not in basis: {w!r}")
            st, en, lm = m; u, v = w[:st], w[en:]
            kids = [u + tm + v for _, tm in lm_to_tail[lm]]
            todo = [c for c in kids if c not in memo]
            if todo: stack.extend(todo); continue
            res = {}
            for c, tm in lm_to_tail[lm]:
                for i, cc in memo[u + tm + v].items():
                    res[i] = norm(res.get(i, zero) + c * cc)
                    if res[i] == 0: del res[i]
            memo[w] = res; stack.pop()
        return memo[word]
    return len(standard), nf, standard, zero, norm


def reduce_poly(nf, P, zero, norm):
    out = {}
    for c, w in P:
        for i, cc in nf(w).items():
            out[i] = norm(out.get(i, zero) + c * cc)
            if out[i] == 0: del out[i]
    return out


def min_poly(nf, standard, Xterms, theta, ring, zero, norm, kmax=None):
    if kmax is None: kmax = 3 * len(standard)  # never binding in practice
    inv = (lambda a: pow(a, P31 - 2, P31)) if ring == "Fp" else (lambda a: Fr(1) / a)
    one = 1 if ring == "Fp" else Fr(1)
    colcache = {}
    def apply_theta(vec):
        out = {}
        for idx, c in vec.items():
            col = colcache.get(idx)
            if col is None: col = colcache[idx] = nf(tuple(theta) + standard[idx])
            for i, cc in col.items():
                out[i] = norm(out.get(i, zero) + c * cc)
                if out[i] == 0: del out[i]
        return out
    rows = []; ri = reduce_poly(nf, Xterms, zero, norm)
    i = 0
    while i <= kmax:
        cur = dict(ri); comb = [zero] * i + [one]
        for (rv, cv) in rows:
            piv = next(iter(rv))
            if cur.get(piv):
                f = norm(cur[piv] * inv(rv[piv]))
                for mm, cval in rv.items():
                    cur[mm] = norm(cur.get(mm, zero) - f * cval)
                    if cur[mm] == 0: del cur[mm]
                cv = cv + [zero] * (len(comb) - len(cv)); comb = comb + [zero] * (len(cv) - len(comb))
                comb = [norm(a - f * b) for a, b in zip(comb, cv)]
        if not cur: return [norm(c) for c in comb]
        piv = min(cur); ip = inv(cur[piv])
        rows.append(({m: norm(c * ip) for m, c in cur.items()}, [norm(c * ip) for c in comb]))
        ri = apply_theta(ri); i += 1
    return None


def tower_prediction(case, q, t, ymax):
    """g_k = prod_{j=0}^{ymax} q^{2 j d} p*(L / q^{2 j})  as ascending coeff list (over Q)."""
    q = Fr(q); t = Fr(t)
    def mul(a, b):
        r = [Fr(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): r[i + j] += x * y
        return r
    ps = [Fr(c) for c in PSTAR[case](t)]; d = len(ps) - 1
    acc = [Fr(1)]
    for j in range(ymax + 1):
        dil = [ps[i] * q ** (2 * j * (d - i)) for i in range(len(ps))]  # q^{2 j d} p*(L/q^{2j})
        acc = mul(acc, dil)
    return acc


def build_basis(case, n, q, t, char, scratch):
    untw, fp = UNTW[case]
    ms = f"{scratch}/z{n}-{case}-{q}-{t}-c{char}.ms"
    bs = ms[:-3] + ".basis"
    subprocess.run(["python3", EMIT, "--n", str(n), "--untwist", untw, "--framing-power", fp,
                    "--q", str(q), "--t", str(t), "--char", str(char), "--max-untwist-level", "2",
                    "--out", ms], check=True, capture_output=True)
    subprocess.run([F4, "-m", "200000", "-i", ms, "-o", bs], check=True, capture_output=True)
    return bs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--cases", nargs="+", default=["linear", "jordan", "framing"])
    ap.add_argument("--char", type=int, default=P31, help="P31 default (mod p); 0 for exact Q")
    ap.add_argument("--grid", nargs="+", default=["5,7", "11,13", "3,10"])
    ap.add_argument("--scratch", default="/tmp")
    args = ap.parse_args()
    ring = "Q" if args.char == 0 else "Fp"
    parse = core.parse_basis_Q if ring == "Q" else (lambda f: core.parse_basis_Fp(f))
    grid = [tuple(int(x) for x in g.split(",")) for g in args.grid]
    for case in args.cases:
        print(f"\n=== {case}  (n={args.n}, {'Q' if ring=='Q' else 'mod 2^31-1'}) ===", flush=True)
        for (q, t) in grid:
            bs = build_basis(case, args.n, q, t, args.char, args.scratch)
            polys = parse(bs)
            N, nf, standard, zero, norm = make_reducer(polys, ring)
            one = 1 if ring == "Fp" else Fr(1)
            line = f"  (q,t)=({q},{t}) dim={N}: "
            for k in range(3, args.n + 1):
                Xk = slider(k, (q % P31) if ring == "Fp" else Fr(q), one,
                            (lambda a, b: (a * b) % P31) if ring == "Fp" else (lambda a, b: a * b))
                mp = min_poly(nf, standard, Xk, theta_word(k), ring, zero, norm)
                deg = len(mp) - 1 if mp else None
                tag = ""
                if ring == "Q":
                    pred = tower_prediction(case, q, t, k - 3)
                    tag = "  == tower?" + (" YES" if [Fr(c) for c in mp] == pred else " NO")
                line += f"deg(Theta_{k})={deg}{tag}   "
            print(line, flush=True)


if __name__ == "__main__":
    main()
