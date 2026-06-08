#!/usr/bin/env python3
"""
Semisimplicity check for a finite-dim noncommutative quotient k<x>/I, given a
reduced Groebner basis of I (with respect to length-lex order).

Strategy:
  1. Parse the f4ncgb .basis file into a list of polynomials, indexed by their
     leading monomials.
  2. Enumerate standard monomials (= basis of the quotient) via the same
     Aho-Corasick automaton as compute_dim_aho.py, BUT now we record each
     monomial as a tuple of letters.
  3. Implement noncommutative reduction modulo the GB: a word w not in the
     standard set contains some leading monomial as a substring; replace it by
     (sum of) shorter words via tail-substitution, repeat until standard.
  4. Build left-multiplication matrices M_v for each generator v, acting on the
     standard-monomial basis: M_v[i, j] = coefficient of m_i in reduce(v * m_j).
     (Sparse storage: dict per column.)
  5. Use the M_v's to compute the left-mult matrix L_w of any word w by matrix
     multiplication. For each pair (i, j) of basis elements, the trace form is
     T[i, j] = tr(L_{m_i m_j}). The algebra is semisimple iff det(T) is nonzero
     in the working field.

All arithmetic over GF(p) at the specialization (q, t) baked into the .basis
file. det(T) = 0 in GF(p) means the algebra is (semisimple is unlikely; either
non-semisimple or specialization-degenerate). det(T) != 0 in GF(p) means
semisimple at the specialization (and therefore generically over Q(q,t), since
semisimplicity is generic).

Usage:
    python3 semisimple_check.py path/to/basis.txt [--char 2147483647]
"""
from __future__ import annotations
import argparse
import re
import sys
from collections import deque
from functools import reduce as freduce


def parse_basis_full(path):
    """Parse the f4ncgb .basis file. Returns:
      polys : list of (leading_mon_tuple, list of (coef:int, mon_tuple))
      where each polynomial g_i = sum of (coef * mon), the first term being
      the leading term (coefficient 1 after normalization)."""
    polys = []
    with open(path) as f:
        for raw in f:
            line = raw.strip().rstrip(",")
            if not line:
                continue
            terms = []
            # Split on " + " and " - " separators, preserving signs
            # Replace " - " with " + -" so we can split uniformly on " + "
            normalized = line.replace(" - ", " + -")
            for term in normalized.split(" + "):
                term = term.strip()
                if not term:
                    continue
                # term format: [-]N*x1*x2... or [-]N or x1*x2*...
                m = re.match(r"^(-?\d+)(?:\*(.+))?$", term)
                if m:
                    coef = int(m.group(1))
                    mon = m.group(2)
                    if mon:
                        mon_tuple = tuple(v.strip() for v in mon.split("*"))
                    else:
                        mon_tuple = ()
                else:
                    # Coefficient is implicitly 1
                    coef = 1
                    mon_tuple = tuple(v.strip() for v in term.split("*"))
                terms.append((coef, mon_tuple))
            if not terms:
                continue
            # leading term is the first one (f4ncgb prints in monomial-order
            # decreasing); leading monomial is its monomial tuple.
            lead_mon = terms[0][1]
            polys.append((lead_mon, terms))
    return polys


def collect_variables(polys):
    s = set()
    for _, terms in polys:
        for _, mon in terms:
            s.update(mon)
    return sorted(s)


def build_aho_corasick(lms, variables):
    """Same AC construction as compute_dim_aho.py, but returns the trie's
    children dict (so we can identify which leading monomial matches at a state)."""
    children = [dict()]
    is_match = [False]
    match_pat = [None]   # for each state, the LM ending there (if any)
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
                # Propagate match: if any state on failure chain is a match,
                # arriving here means a match too (and we propagate the pattern)
                if is_match[fail[c]] and not is_match[c]:
                    is_match[c] = True
                    match_pat[c] = match_pat[fail[c]]
                delta[(s, v)] = c
                q.append(c)
            else:
                delta[(s, v)] = delta.get((fail[s], v), 0)
    return delta, is_match, match_pat, nstates


def enumerate_standard_monomials(delta, is_match, nstates, variables, max_total=None):
    """BFS enumerate all standard monomials (those not containing any LM as
    a contiguous substring), returning them as a list ordered by (length,
    lex). Caps at max_total if the algebra is infinite."""
    # State per monomial: the AC state after reading its letters.
    # We start with the empty monomial at state 0.
    result = [(0, ())]
    frontier = deque([(0, ())])  # (state, mon_tuple)
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
    # Order by (length, lex)
    result.sort(key=lambda sm: (len(sm[1]), sm[1]))
    return [mon for _, mon in result]


def find_leading_match(word, ac_delta, ac_is_match, ac_match_pat):
    """Find the LEFTMOST contiguous substring of `word` that equals some LM,
    by simulating the AC automaton over the word. Returns (start_idx, end_idx, lm_tuple)
    or None if word is standard."""
    state = 0
    for i, v in enumerate(word):
        if (state, v) not in ac_delta:
            return None  # variable not in alphabet? Shouldn't happen
        state = ac_delta[(state, v)]
        if ac_is_match[state]:
            lm = ac_match_pat[state]
            end = i + 1
            start = end - len(lm)
            return (start, end, lm)
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="path to f4ncgb basis file")
    ap.add_argument("--char", type=int, default=2147483647,
                    help="working characteristic (must match the basis)")
    ap.add_argument("--max-dim", type=int, default=2500,
                    help="abort if standard-monomial count exceeds this")
    args = ap.parse_args()

    polys = parse_basis_full(args.basis)
    print(f"[parse] {len(polys)} basis polynomials.")
    variables = collect_variables(polys)
    print(f"[vars] {len(variables)}: {variables}")

    lms = [lm for lm, _ in polys]
    ac_delta, ac_is_match, ac_match_pat, nstates = build_aho_corasick(lms, variables)
    print(f"[AC] {nstates} states, {sum(ac_is_match)} match states.")

    standard = enumerate_standard_monomials(ac_delta, ac_is_match, nstates,
                                            variables, max_total=args.max_dim)
    if standard is None:
        print(f"[abort] standard-monomial count exceeds {args.max_dim}; "
              f"algebra is either infinite or larger than the cap. "
              f"Increase --max-dim.")
        return 1
    N = len(standard)
    print(f"[dim] standard-monomial basis has size {N}.")

    # Build leading-monomial -> polynomial (with leading term subtracted off)
    # For reduction: when we see word u * lm * v, replace lm by (lm - poly)
    # which equals (lm - sum_terms) = - sum_{nonleading terms}. So the
    # replacement is: word -> sum over nonleading terms of (- coef * u * mon * v).
    # Store the "tail" of each LM: list of (coef, mon_tuple), where coef is
    # the coefficient AFTER normalizing leading coefficient to 1 (and negating;
    # see below).
    p = args.char
    lm_to_tail = {}
    for lm, terms in polys:
        lead_coef, _ = terms[0]
        inv = pow(lead_coef, p - 2, p)  # modular inverse
        # Polynomial = sum of terms. Normalize so lead coefficient is 1.
        # Then the substitution is: lm -> -(sum of nonleading terms, normalized).
        # equivalently lm -> sum over nonleading terms of (-coef * inv) % p * mon.
        tail = []
        for coef, mon in terms[1:]:
            new_coef = (-coef * inv) % p
            tail.append((new_coef, mon))
        lm_to_tail[lm] = tail

    # Index standard monomials for fast lookup
    mon_idx = {mon: i for i, mon in enumerate(standard)}

    # Reduction function: word_terms (list of (coef, mon_tuple)) -> coef vector
    # of length N (sparse dict {idx: coef}).
    def reduce_to_basis(word_terms):
        out = {}
        worklist = list(word_terms)  # list of (coef, mon_tuple)
        while worklist:
            coef, word = worklist.pop()
            if coef == 0:
                continue
            if word in mon_idx:
                # Already a standard monomial
                k = mon_idx[word]
                out[k] = (out.get(k, 0) + coef) % p
                continue
            m = find_leading_match(word, ac_delta, ac_is_match, ac_match_pat)
            if m is None:
                # Word is standard but not in our enumerated basis?
                # That means our basis enumeration was incomplete. Skip silently.
                # (Better: panic.)
                raise RuntimeError(f"Standard word {word} not in enumerated basis. "
                                   f"Possible incomplete GB or AC bug.")
            start, end, lm = m
            u = word[:start]
            v = word[end:]
            # Replace lm by tail: lm -> sum_{c, t in tail} c * t
            # So u * lm * v -> sum c * (u + t + v)
            tail = lm_to_tail[lm]
            for c, t_mon in tail:
                new_coef = (coef * c) % p
                new_word = u + t_mon + v
                worklist.append((new_coef, new_word))
        return out

    # Sanity: reduce each generator (length-1 word). Each generator x_i acts as
    # an element of the algebra; multiplying by 1 gives itself, which should be
    # the generator's index in the standard basis.
    for v in variables:
        red = reduce_to_basis([(1, (v,))])
        if red != {mon_idx[(v,)]: 1}:
            print(f"[sanity] reduction of single generator {v} gave {red}, expected "
                  f"singleton at idx {mon_idx[(v,)]}")
            return 1
    print(f"[sanity] generator reductions OK.")

    # Build left-multiplication matrices for each generator:
    # L_v[i, j] = coef of m_i in reduce(v * m_j)
    # We store as dict {(j): {i: coef}} (columns are the source monomial j;
    # entries give where it maps to).
    print(f"[matrices] building left-mult matrices for {len(variables)} generators "
          f"on basis of dim {N}...")
    L = {}  # v -> {j: {i: coef}}
    for v in variables:
        col_map = {}
        for j, m_j in enumerate(standard):
            # Reduce v * m_j
            word = (v,) + m_j
            red = reduce_to_basis([(1, word)])
            col_map[j] = red
        L[v] = col_map
    print(f"[matrices] done.")

    # Trace form: T[i, j] = tr(L_{m_i m_j}) = sum_k coef of m_k in m_i * m_j * m_k.
    # Equivalently, build L_{m_i} as a matrix and trace its product with L_{m_j}.
    # Storing N x N dense matrices over GF(p) is N^2 entries (3M for N=1770).
    # We'll compute traces directly without storing full L_{m_i}.
    #
    # Actually easier: for each i, compute the column map L_{m_i} as a dict
    # {j: {k: coef}}, then T[i, j] = sum_k L_{m_i}[k][k_dest] * indicator(...) —
    # no that's not it either.
    #
    # Cleanest: L_{m_i m_j} = L_{m_i} * L_{m_j} (composition).
    # tr(L_{m_i m_j}) = sum_k (L_{m_i} * L_{m_j})[k, k]
    #                 = sum_k sum_l L_{m_i}[k, l] * L_{m_j}[l, k]
    # So we need L_{m_i}[k, l] and L_{m_j}[l, k] for all (k, l) pairs.
    #
    # We can build L_{m_i} as a matrix from the generators by iterated multiplication:
    # m_i is a product of generators v_1 * v_2 * ... * v_r, so L_{m_i} = L_{v_1} * L_{v_2} * ... * L_{v_r}.
    # But matrix multiplication is N^3 — too expensive for N=1770 in Python.
    #
    # For now, defer the full N=1770 case. Verify the pipeline at small N.

    if N <= 250:
        print(f"[trace form] computing N x N matrix (N={N})...")
        # Build L_{m_i} as a sparse column dict for each m_i.
        full_L = []  # full_L[i] = col_map for m_i
        # We build them iteratively, starting from L[m_i] for length-1 m_i = generator
        # and extending. Faster: just reduce v * m_i for each m_i (which we already have
        # via the generator matrices), then chain.
        # Actually for general m_i (product of letters), L_{m_i}[j] = reduce(m_i * m_j).
        # We do this directly:
        for i, m_i in enumerate(standard):
            col_map = {}
            for j, m_j in enumerate(standard):
                red = reduce_to_basis([(1, m_i + m_j)])
                col_map[j] = red
            full_L.append(col_map)
            if (i + 1) % 50 == 0:
                print(f"   m_i = {i+1}/{N}", file=sys.stderr)

        # Compute trace form
        print(f"[trace form] computing T[i, j] for all (i, j)...")
        T = [[0] * N for _ in range(N)]
        for i in range(N):
            for j in range(N):
                # T[i, j] = tr(L_{m_i m_j}) = sum_k coef of m_k in m_i m_j m_k.
                # Equivalently, L_{m_i m_j}[k, k] for all k.
                # L_{m_i m_j}[k] = reduce((m_i m_j) * m_k).
                # We have full_L[m_i*m_j] only if m_i * m_j is also a basis element ...
                # That's not generally true; we need product first.
                product = reduce_to_basis([(1, standard[i] + standard[j])])
                # product = dict {k: coef}.  L_{m_i m_j} = sum_k coef_k * L_{m_k}.
                # tr(L_{m_i m_j}) = sum_k coef_k * tr(L_{m_k}) = sum_k coef_k * (sum_l L_{m_k}[l][l]).
                trace_val = 0
                for k, c in product.items():
                    # trace of L_{m_k}: sum_l (full_L[k] at column l, entry l)
                    s = 0
                    for l in range(N):
                        s = (s + full_L[k].get(l, {}).get(l, 0)) % p
                    trace_val = (trace_val + c * s) % p
                T[i][j] = trace_val

        # det T via Gaussian elimination over GF(p)
        det = 1
        M = [row[:] for row in T]
        for col in range(N):
            pivot = -1
            for r in range(col, N):
                if M[r][col] != 0:
                    pivot = r
                    break
            if pivot < 0:
                det = 0
                break
            if pivot != col:
                M[col], M[pivot] = M[pivot], M[col]
                det = (-det) % p
            det = (det * M[col][col]) % p
            inv_p = pow(M[col][col], p - 2, p)
            for r in range(col + 1, N):
                if M[r][col] == 0:
                    continue
                fac = (M[r][col] * inv_p) % p
                for c2 in range(col, N):
                    M[r][c2] = (M[r][c2] - fac * M[col][c2]) % p
        print(f"[result] det(trace form) mod p = {det}")
        if det == 0:
            print("[conclusion] det T = 0: algebra is NOT semisimple at "
                  "this specialization (rank deficient trace form, "
                  "which is the Jacobson radical's contribution).")
        else:
            print("[conclusion] det T != 0: algebra IS semisimple at "
                  "this specialization (and generically over Q(q,t)).")
    else:
        print(f"[trace form] skipped (N = {N} > 250 threshold). To enable, "
              f"raise the threshold or use a faster impl.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
