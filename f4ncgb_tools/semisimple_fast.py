#!/usr/bin/env python3
"""
Fast semisimplicity check via trace form, using NumPy for the matrix work.

Strategy (same as semisimple_check.py but with linear algebra in NumPy):
  1. Parse basis polynomials and enumerate standard monomials via AC.
  2. Reduction routine: word -> linear combo of standard monomials (modular int).
  3. Build left-multiplication matrices L_v for each generator v (sparse).
  4. Iteratively compute L_{m_l} for each basis monomial m_l = v_1 ... v_r, as
     a NumPy dense matrix L_v_1 @ L_v_2 @ ... @ L_v_r.
  5. Compute t_l = trace(L_{m_l}) for each l.
  6. Compute structure constant tensor: for each (i, j), m_i * m_j = sum_l p_l^{ij} m_l.
  7. T[i, j] = sum_l p_l^{ij} * t_l. This is a single matrix-vector product if we
     view the structure constants as an N x N x N tensor flattened to an
     (N*N) x N matrix.
  8. Compute rank(T) over GF(p); semisimple iff rank == N (det T != 0).

This still scales as N^3 for the matrix multiplication chain, but uses NumPy
which is much faster than pure Python. Memory is dominated by the L_v matrices
(few of them, each N x N).
"""
from __future__ import annotations
import argparse
import re
import sys
import numpy as np
from collections import deque

P = 2147483647  # default Mersenne prime; matches f4ncgb default


def parse_basis_full(path):
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
                m = re.match(r"^(-?\d+)(?:\*(.+))?$", term)
                if m:
                    coef = int(m.group(1))
                    mon = m.group(2)
                    if mon:
                        mon_tuple = tuple(v.strip() for v in mon.split("*"))
                    else:
                        mon_tuple = ()
                else:
                    coef = 1
                    mon_tuple = tuple(v.strip() for v in term.split("*"))
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


def reduce_word_to_basis(word, mon_idx, lm_to_tail, ac_delta, ac_is_match, ac_match_pat, p):
    """Reduce a single word (with coef 1) to a linear combo of basis monomials.
    Returns a dict {basis_idx: coef}."""
    out = {}
    worklist = [(1, word)]
    while worklist:
        coef, w = worklist.pop()
        if coef == 0:
            continue
        if w in mon_idx:
            k = mon_idx[w]
            out[k] = (out.get(k, 0) + coef) % p
            if out[k] == 0:
                del out[k]
            continue
        m = find_leading_match(w, ac_delta, ac_is_match, ac_match_pat)
        if m is None:
            raise RuntimeError(f"Standard word {w} not enumerated. Incomplete GB?")
        start, end, lm = m
        u = w[:start]
        v = w[end:]
        tail = lm_to_tail[lm]
        for c, t_mon in tail:
            new_coef = (coef * c) % p
            new_word = u + t_mon + v
            worklist.append((new_coef, new_word))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="path to f4ncgb basis file")
    ap.add_argument("--char", type=int, default=2147483647)
    ap.add_argument("--max-dim", type=int, default=2500)
    ap.add_argument("--rank-only", action="store_true",
                    help="only print rank of trace form, not det")
    args = ap.parse_args()
    p = args.char

    polys = parse_basis_full(args.basis)
    print(f"[parse] {len(polys)} basis polynomials.", file=sys.stderr)
    variables = collect_variables(polys)
    print(f"[vars] {len(variables)}: {variables}", file=sys.stderr)

    lms = [lm for lm, _ in polys]
    ac_delta, ac_is_match, ac_match_pat = build_aho_corasick(lms, variables)

    standard = enumerate_standard_monomials(ac_delta, ac_is_match, variables,
                                            max_total=args.max_dim)
    if standard is None:
        print(f"[abort] dim > {args.max_dim}", file=sys.stderr)
        return 1
    N = len(standard)
    print(f"[dim] N = {N}", file=sys.stderr)

    lm_to_tail = {}
    for lm, terms in polys:
        lead_coef = terms[0][0]
        inv = pow(lead_coef, p - 2, p)
        tail = []
        for coef, mon in terms[1:]:
            tail.append(((-coef * inv) % p, mon))
        lm_to_tail[lm] = tail

    mon_idx = {mon: i for i, mon in enumerate(standard)}

    # Build left-mult matrices for each generator as N x N dense uint64 array.
    print(f"[gen-matrices] building L_v for each generator...", file=sys.stderr)
    gen_L = {}
    for v in variables:
        M = np.zeros((N, N), dtype=np.int64)
        for j, m_j in enumerate(standard):
            red = reduce_word_to_basis((v,) + m_j, mon_idx, lm_to_tail,
                                       ac_delta, ac_is_match, ac_match_pat, p)
            for i, c in red.items():
                M[i, j] = c
        gen_L[v] = M
    print(f"[gen-matrices] done. Memory: ~{len(variables) * N * N * 8 / 1e6:.1f} MB",
          file=sys.stderr)

    # For each standard monomial m_l = v_1 v_2 ... v_r, compute L_{m_l} = L_{v_1} @ L_{v_2} @ ... @ L_{v_r}.
    # Memoize via the trie of standard monomials: L_empty = I.
    # We process standard monomials in order of length so each one extends a shorter one
    # by appending one letter on the LEFT: m = v * m', and L_m = L_v @ L_{m'}.
    # Actually m_l corresponds to *left*-multiplication by m_l, and m_l acts as
    # L_{v_1 v_2 ... v_r} = L_{v_1} @ L_{v_2} @ ... @ L_{v_r} (compose).
    # We can build L_{m_l} from L_{m_l[:-1]} @ L_{m_l[-1]}, but each m_l[:-1] may
    # not be a standard monomial (might have reduced) -- so this doesn't work directly.
    # Simpler: m_l is a fixed tuple; L_{m_l} acts on basis by left-mult; for each j,
    # column j of L_{m_l} is reduce(m_l * m_j). We can build it incrementally:
    # start with L_acc = I, then for each letter v in m_l (left-to-right), L_acc = L_v @ L_acc.

    # We only need the traces t_l = tr(L_{m_l}); we don't need to store L_{m_l}.
    print(f"[traces] computing t_l = tr(L_{{m_l}}) for l = 0..{N-1}...", file=sys.stderr)
    traces = np.zeros(N, dtype=np.int64)
    for l, m_l in enumerate(standard):
        if not m_l:
            traces[l] = N % p  # tr(I) = N
            continue
        # L_{m_l} = product of L_v's for v in m_l (left-to-right):
        L = np.eye(N, dtype=np.int64)
        for v in m_l:
            L = (gen_L[v] @ L) % p
        traces[l] = int(np.trace(L) % p)
        if (l + 1) % 50 == 0:
            print(f"   {l+1}/{N}", file=sys.stderr)
    print(f"[traces] done.", file=sys.stderr)

    # Now T[i, j] = sum_l p_l^{ij} * t_l, where p_l^{ij} is the coefficient of m_l in m_i * m_j.
    # Approach: for each (i, j), reduce m_i * m_j and dot with traces.
    print(f"[trace form] computing T (NxN)...", file=sys.stderr)
    T = np.zeros((N, N), dtype=np.int64)
    for i in range(N):
        for j in range(N):
            red = reduce_word_to_basis(standard[i] + standard[j],
                                       mon_idx, lm_to_tail,
                                       ac_delta, ac_is_match, ac_match_pat, p)
            s = 0
            for l, c in red.items():
                s = (s + c * traces[l]) % p
            T[i, j] = s
        if (i + 1) % 25 == 0:
            print(f"   T row {i+1}/{N}", file=sys.stderr)
    print(f"[trace form] done.", file=sys.stderr)

    # Compute rank of T over GF(p) via Gaussian elimination.
    rank = 0
    M = T.copy()
    for col in range(N):
        pivot = -1
        for r in range(rank, N):
            if M[r, col] != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        if pivot != rank:
            M[[rank, pivot]] = M[[pivot, rank]]
        # Eliminate below
        inv = pow(int(M[rank, col]), p - 2, p)
        for r in range(N):
            if r != rank and M[r, col] != 0:
                fac = int(M[r, col] * inv % p)
                M[r] = (M[r] - fac * M[rank]) % p
        rank += 1
    print(f"[result] rank(T) = {rank} (out of N = {N})")
    if rank == N:
        print(f"[conclusion] det T != 0: semisimple at (q,t)=(5,7) in GF({p}).")
        print(f"             (And generically over Q(q,t), since semisimplicity is generic.)")
    else:
        defect = N - rank
        print(f"[conclusion] rank deficit {defect}: NOT semisimple at this specialization.")
        print(f"             Jacobson radical J(A) has dim {defect} here.")
        print(f"             (May be specialization-only or generic non-semisimplicity;")
        print(f"              re-run at independent (q,t) to distinguish.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
