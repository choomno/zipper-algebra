#!/usr/bin/env python3
r"""Compare candidate lasso-twist annihilators at specified parameter points.

The predicted formulas below preserve historical research hypotheses. Agreement
at sampled parameters is not a proof of a generic formula or its minimality.
Solver output has no basis-completeness certificate in this experimental driver.
The Krylov calculation concerns the supplied reduction model.
"""
import argparse, os, subprocess, sys
from fractions import Fraction as Fr
from zipper_algebra import core, Presentation
from zipper_algebra.linear import validate_field_characteristic
from experiments.output import (capture_experiment_provenance, tag_solver_input,
                                tag_solver_basis)

F4 = os.environ.get("F4NCGB", "f4ncgb")
P31 = 2 ** 31 - 1
UNTW = {"linear": ("bigelow", "1"), "square-zero": ("square-zero", "1"), "framing": ("bigelow", "2")}
PSTAR = {  # p*(L) as ascending coeff list, symbolic in t (filled per point)
    "linear":  lambda t: [t, 1],
    "square-zero":  lambda t: [t * t, 2 * t, 1],
    "framing": lambda t: [-t, 0, 1],
}


def slider(k, q, one, mul_scalar):
    """X_k terms as [(coef, word)] over the chosen ring. one/mul_scalar adapt Q vs F_p."""
    def qp(e):  # q**e in the ring
        r = one
        for _ in range(e): r = mul_scalar(r, q)
        return r
    X = [(q, ("y1",)), (one - q, ()),
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


def _characteristic(ring, char):
    if ring == "Q":
        if char is not None:
            validate_field_characteristic(char)
            if char:
                raise ValueError("Q reduction requires characteristic zero")
        return 0
    if ring != "Fp":
        raise ValueError("ring must be 'Q' or 'Fp'")
    char = P31 if char is None else char
    validate_field_characteristic(char)
    if not char:
        raise ValueError("Fp reduction requires a positive prime characteristic")
    return char


def make_reducer(polys, ring, *, char=None, variables=None, max_total=2_000_000):
    """Build a reduction model over Q or the specified prime field.

    The optional ordered alphabet includes generators absent from the supplied
    polynomials. Finite enumeration and solver exit do not certify completeness.
    """
    char = _characteristic(ring, char)
    lms = [lm for lm, _ in polys]
    declared = core.collect_variables(polys)
    variables = declared if variables is None else list(variables)
    if len(set(variables)) != len(variables) or set(declared) - set(variables):
        raise ValueError("declared alphabet must be distinct and contain every basis variable")
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(delta, is_match, variables, max_total=max_total)
    if standard is None:
        raise RuntimeError("standard-word enumeration exceeded its budget")
    mon_idx = {m: i for i, m in enumerate(standard)}
    if ring == "Fp":
        lm_to_tail = core.lm_to_tail_Fp(polys, char); zero = 0
        norm = lambda c: c % char
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


def min_poly(nf, standard, Xterms, theta, ring, zero, norm, kmax=None, *, char=None):
    """Krylov dependence for the supplied operator and reduction model."""
    char = _characteristic(ring, char)
    if kmax is None: kmax = 3 * len(standard)  # never binding in practice
    inv = (lambda a: pow(a, -1, char)) if ring == "Fp" else (lambda a: Fr(1) / a)
    one = 1 if ring == "Fp" else Fr(1)
    # apply Theta = word w_1..w_m as a COMPOSITION of single-generator left-mults:
    # L_w = M_{w_1} .. M_{w_m}, so apply letters right-to-left; each M_g needs only nf(g*mon)
    # (short words) -> bounded memory, scales to n>=5. gcol[(g,idx)] cached lazily.
    gcol = {}
    def Lg(g, vec):
        out = {}
        for idx, c in vec.items():
            key = (g, idx); col = gcol.get(key)
            if col is None: col = gcol[key] = nf((g,) + standard[idx])
            for i, cc in col.items():
                out[i] = norm(out.get(i, zero) + c * cc)
                if out[i] == 0: del out[i]
        return out
    def apply_theta(vec):
        for g in reversed(tuple(theta)): vec = Lg(g, vec)
        return vec
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


def predicted_min_poly(case, q, t, k):
    """Historical candidate pattern as ascending coefficients; not a generic theorem."""
    q = Fr(q); t = Fr(t)
    def mul(a, b):
        r = [Fr(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): r[i + j] += x * y
        return r
    ps = [Fr(c) for c in PSTAR[case](t)]; d = len(ps) - 1
    if k <= 3: return ps
    e = k - 2
    dil = [ps[i] * q ** (e * (d - i)) for i in range(len(ps))]  # q^{e d} p*(L / q^e)
    return mul(ps, dil)


def build_basis(case, n, q, t, char, scratch):
    validate_field_characteristic(char)
    producer = capture_experiment_provenance(__file__)
    untw, fp = UNTW[case]
    ms = f"{scratch}/z{n}-{case}-{q}-{t}-c{char}.ms"
    bs = ms[:-3] + ".basis"
    Presentation(n, untw, framing_power=int(fp), max_untwist_level=2).emit_ms(ms, q, t, char, provenance=False)
    inputs = tag_solver_input(ms, producer=producer)
    subprocess.run([F4, "-m", "200000", "-i", ms, "-o", bs], check=True, capture_output=True)
    tag_solver_basis(bs, producer=producer, inputs=inputs)
    return bs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--cases", nargs="+", default=["linear", "square-zero", "framing"])
    ap.add_argument("--char", type=int, default=P31, help="P31 default (mod p); 0 for exact Q")
    ap.add_argument("--grid", nargs="+", default=["5,7", "11,13", "3,10"])
    ap.add_argument("--scratch", default="/tmp")
    args = ap.parse_args()
    validate_field_characteristic(args.char)
    ring = "Q" if args.char == 0 else "Fp"
    parse = core.parse_basis_Q if ring == "Q" else (lambda f: core.parse_basis_Fp(f))
    grid = [tuple(int(x) for x in g.split(",")) for g in args.grid]
    for case in args.cases:
        print(f"\n=== {case}  (n={args.n}, {'Q' if ring=='Q' else f'mod {args.char}'}) ===", flush=True)
        for (q, t) in grid:
            bs = build_basis(case, args.n, q, t, args.char, args.scratch)
            polys = parse(bs)
            N, nf, standard, zero, norm = make_reducer(polys, ring, char=args.char,
                                                   variables=Presentation(args.n).var_names())
            one = 1 if ring == "Fp" else Fr(1)
            line = f"  (q,t)=({q},{t}) standard_words={N}: "
            for k in range(3, args.n + 1):
                Xk = slider(k, (q % args.char) if ring == "Fp" else Fr(q), one,
                            (lambda a, b: (a * b) % args.char) if ring == "Fp" else (lambda a, b: a * b))
                mp = min_poly(nf, standard, Xk, theta_word(k), ring, zero, norm, char=args.char)
                deg = len(mp) - 1 if mp else None
                tag = ""
                if ring == "Q":
                    pred = predicted_min_poly(case, q, t, k)
                    tag = "  == pattern?" + (" YES" if mp is not None and [Fr(c) for c in mp] == pred else " NO")
                line += f"observed Krylov degree(Theta_{k})={deg}{tag}   "
            print(line, flush=True)


if __name__ == "__main__":
    main()
