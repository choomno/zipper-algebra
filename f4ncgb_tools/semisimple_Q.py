#!/usr/bin/env python3
"""
Semisimplicity check OVER THE RATIONALS (Q), not F_p.

Combined with f4ncgb --char 0, this gives a characteristic-zero verification
of (non-)semisimplicity via Dieudonne's theorem (which is an "iff" in char 0).

  - rank T = N => semisimple (det T != 0).
  - rank T < N => J(A) = ker T, dim J(A) = N - rank T. Genuinely non-semisimple.

Cost is higher than the F_p version because Python's fractions.Fraction is
slow (bignum rationals). Feasible at small N.
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from collections import deque
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def parse_basis_Q(path):
    """Parse f4ncgb char-0 basis output. Coefficients may be rationals like '1/5'."""
    polys = []
    with open(path) as f:
        for raw in f:
            line = raw.strip().rstrip(",")
            if not line:
                continue
            terms = []
            normalized = line.replace(" - ", " + -")
            for term in normalized.split(" + "):
                term = term.strip()
                if not term:
                    continue
                # Term forms (after sign-normalization):
                #   [-]<num>[/<denom>]*<mon1>*<mon2>...   (with coefficient)
                #   [-]<num>[/<denom>]                    (constant)
                #   [-]<mon1>*<mon2>...                   (implicit ±1)
                # Detect leading sign:
                sign = 1
                rest = term
                if rest.startswith("-"):
                    sign = -1
                    rest = rest[1:]
                # Try to consume a numeric coefficient (with optional fraction)
                m = re.match(r"^(\d+)(?:/(\d+))?(?:\*(.+))?$", rest)
                if m:
                    num = int(m.group(1))
                    den = int(m.group(2)) if m.group(2) else 1
                    coef = Fraction(sign * num, den)
                    mon = m.group(3)
                    if mon:
                        mon_tuple = tuple(v.strip() for v in mon.split("*"))
                    else:
                        mon_tuple = ()
                else:
                    # No leading digit: coefficient is sign*1, rest is the monomial
                    coef = Fraction(sign)
                    mon_tuple = tuple(v.strip() for v in rest.split("*"))
                terms.append((coef, mon_tuple))
            if not terms:
                continue
            polys.append((terms[0][1], terms))
    return polys


def collect_variables(polys):
    s = set()
    for _, terms in polys:
        for _, mon in terms:
            s.update(mon)
    return sorted(s)


def build_aho_corasick(lms, variables):
    children = [dict()]
    is_match = [False]
    match_pat = [None]
    for lm in lms:
        s = 0
        for v in lm:
            if v not in children[s]:
                children.append(dict())
                is_match.append(False)
                match_pat.append(None)
                children[s][v] = len(children) - 1
            s = children[s][v]
        is_match[s] = True
        match_pat[s] = lm
    nstates = len(children)
    fail = [0] * nstates
    delta = {}
    q = deque()
    for v in variables:
        if v in children[0]:
            c = children[0][v]
            fail[c] = 0
            delta[(0, v)] = c
            q.append(c)
        else:
            delta[(0, v)] = 0
    while q:
        s = q.popleft()
        for v in variables:
            if v in children[s]:
                c = children[s][v]
                fs = fail[s]
                while fs != 0 and v not in children[fs]:
                    fs = fail[fs]
                fail[c] = children[fs][v] if v in children[fs] and children[fs][v] != c else 0
                if is_match[fail[c]] and not is_match[c]:
                    is_match[c] = True
                    match_pat[c] = match_pat[fail[c]]
                delta[(s, v)] = c
                q.append(c)
            else:
                delta[(s, v)] = delta.get((fail[s], v), 0)
    return delta, is_match, match_pat


def enumerate_standard_monomials(delta, is_match, variables, max_total=None):
    result = [(0, ())]
    frontier = deque([(0, ())])
    while frontier:
        state, mon = frontier.popleft()
        for v in variables:
            t = delta[(state, v)]
            if is_match[t]:
                continue
            new_mon = mon + (v,)
            result.append((t, new_mon))
            frontier.append((t, new_mon))
            if max_total is not None and len(result) > max_total:
                return None
    result.sort(key=lambda sm: (len(sm[1]), sm[1]))
    return [mon for _, mon in result]


def find_leading_match(word, ac_delta, ac_is_match, ac_match_pat):
    state = 0
    for i, v in enumerate(word):
        if (state, v) not in ac_delta:
            return None
        state = ac_delta[(state, v)]
        if ac_is_match[state]:
            lm = ac_match_pat[state]
            end = i + 1
            start = end - len(lm)
            return (start, end, lm)
    return None


def reduce_to_basis_Q(word_terms, mon_idx, lm_to_tail, ac_delta, ac_is_match, ac_match_pat):
    """Reduce list of (Fraction, mon_tuple) to {basis_idx: Fraction}."""
    out = {}
    worklist = list(word_terms)
    while worklist:
        coef, w = worklist.pop()
        if coef == 0:
            continue
        if w in mon_idx:
            k = mon_idx[w]
            new = out.get(k, Fraction(0)) + coef
            if new == 0 and k in out:
                del out[k]
            else:
                out[k] = new
            continue
        m = find_leading_match(w, ac_delta, ac_is_match, ac_match_pat)
        if m is None:
            raise RuntimeError(f"Standard word {w} not in basis")
        start, end, lm = m
        u, v = w[:start], w[end:]
        for c, t_mon in lm_to_tail[lm]:
            worklist.append((coef * c, u + t_mon + v))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="path to f4ncgb basis file (char=0)")
    ap.add_argument("--max-dim", type=int, default=500)
    args = ap.parse_args()

    polys = parse_basis_Q(args.basis)
    variables = collect_variables(polys)
    print(f"[parse] {len(polys)} polys over Q, vars = {variables}", flush=True)

    lms = [lm for lm, _ in polys]
    ac_d, ac_m, ac_p = build_aho_corasick(lms, variables)
    standard = enumerate_standard_monomials(ac_d, ac_m, variables, max_total=args.max_dim)
    if standard is None:
        print(f"[abort] dim > {args.max_dim}")
        return 1
    N = len(standard)
    print(f"[dim] N = {N}", flush=True)

    lm_to_tail = {}
    for lm, terms in polys:
        lead = terms[0][0]
        inv = Fraction(1) / lead
        lm_to_tail[lm] = [(-c * inv, mon) for c, mon in terms[1:]]
    mon_idx = {mon: i for i, mon in enumerate(standard)}

    # Build generator left-mult matrices (sparse, over Q)
    print(f"[gen-L] building L_v for {len(variables)} generators...", flush=True)
    gen_L = {}
    for v in variables:
        col = {}
        for j, m_j in enumerate(standard):
            red = reduce_to_basis_Q([(Fraction(1), (v,) + m_j)],
                                    mon_idx, lm_to_tail, ac_d, ac_m, ac_p)
            if red:
                col[j] = red
        gen_L[v] = col
    print(f"[gen-L] done.", flush=True)

    def apply_gen(v, vec):
        out = {}
        col = gen_L[v]
        for j, c in vec.items():
            if j not in col:
                continue
            for k, ck in col[j].items():
                new = out.get(k, Fraction(0)) + c * ck
                if new == 0 and k in out:
                    del out[k]
                else:
                    out[k] = new
        return out

    print(f"[traces] computing tr(L_{{m_l}}) for l=0..{N-1}...", flush=True)
    traces = [Fraction(0)] * N
    for l, m_l in enumerate(standard):
        if not m_l:
            traces[l] = Fraction(N)
            continue
        s = Fraction(0)
        for j in range(N):
            vec = {j: Fraction(1)}
            for v in reversed(m_l):
                vec = apply_gen(v, vec)
                if not vec:
                    break
            s += vec.get(j, Fraction(0))
        traces[l] = s
    print(f"[traces] done.", flush=True)

    print(f"[T] filling N x N trace form...", flush=True)
    T = [[Fraction(0)] * N for _ in range(N)]
    for i in range(N):
        m_i = standard[i]
        for j in range(N):
            vec = {j: Fraction(1)}
            for v in reversed(m_i):
                vec = apply_gen(v, vec)
                if not vec:
                    break
            s = Fraction(0)
            for l, c in vec.items():
                s += c * traces[l]
            T[i][j] = s
    print(f"[T] done.", flush=True)

    # Gaussian elimination over Q to compute rank and det
    print(f"[rank] computing over Q...", flush=True)
    M = [row[:] for row in T]
    rank = 0
    sign = 1
    det_num = Fraction(1)
    for col in range(N):
        pivot = -1
        for r in range(rank, N):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            det_num = Fraction(0)
            continue
        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]
            sign = -sign
        det_num *= M[rank][col]
        inv = Fraction(1) / M[rank][col]
        for r in range(N):
            if r == rank or M[r][col] == 0:
                continue
            fac = M[r][col] * inv
            for c2 in range(col, N):
                M[r][c2] -= fac * M[rank][c2]
        rank += 1
    det = det_num * sign if det_num != 0 else Fraction(0)
    defect = N - rank
    print(f"\n[CHAR-0 RESULT]")
    print(f"  N         = {N}")
    print(f"  rank(T)   = {rank}    (over Q)")
    print(f"  deficit   = {defect}  (dim J(A) by Dieudonne in char 0)")
    print(f"  det T     = {det}")
    if defect == 0:
        print(f"\n[THEOREM] det T != 0 over Q  =>  Z_n^variant is SEMISIMPLE over Q.")
    else:
        print(f"\n[THEOREM] rank deficit {defect} over Q + char 0 Dieudonne")
        print(f"           =>  dim J(Z_n^variant) = {defect} over Q.")
        print(f"           =>  Z_n^variant is GENUINELY NON-SEMISIMPLE over Q.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
