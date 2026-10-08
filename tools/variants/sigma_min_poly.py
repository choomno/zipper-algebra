#!/usr/bin/env python3
"""Minimal polynomial of left-multiplication by sigma_1 = x_1 on a variant
zipper algebra Z_n^(a) (higher-framing variant).

Conjecture (universal at generic (q, t)):

    min_poly(sigma_1) = (sigma - q) * (sigma + q^{-1}) * (sigma^a - t)

so the number of distinct sigma_1-eigenvalues should be a + 2.

The .basis files we are given are already F_p-specialised at the numerical
point (q, t) = (5, 7) modulo p = 2^31 - 1 (a Mersenne prime). At this
specialisation, the predicted annihilator is

    (sigma - 5)(sigma + 1/5)(sigma^a - 7),

which has (a + 2) distinct roots over F_p as long as none of {5, -1/5, t^{1/a}}
coincide -- generically true.

Pipeline (mirrors compute_dim_aho.py / semisimple_v3.py):

1. parse_basis_Fp on the .basis file.
2. Build the AC automaton from leading monomials and enumerate the standard
   monomial basis of the quotient.
3. Build a sparse matrix M_{ij} for left-multiplication by x_1:
       x_1 * m_j  =  sum_i M_{ij} * m_i  (reduced via Groebner).
4. Find min poly of M over F_p by tracking the Krylov sequence
       v_k = M^k * 1   (1 = empty word = first standard monomial)
   and looking for the smallest linear dependence. For a generic cyclic
   vector this is the full min poly; we then verify on all other standard
   basis vectors that p(M) annihilates them. If not, we extend.
5. Factor the min poly over F_p (sympy.Poly factor_list).
6. Report dimension, min-poly degree, factorisation, # distinct eigenvalues,
   and the expected value a + 2 (a inferred from the file's "k" suffix).

The script also verifies the conjectured annihilating polynomial
   A(sigma) = (sigma - 5)(sigma + 1/5)(sigma^a - 7)
by evaluating A(M) on every standard basis vector and reporting whether the
result is identically zero.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Optional

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
if PARENT not in sys.path:
    sys.path.insert(0, PARENT)

import core

try:
    import sympy
except ImportError:  # pragma: no cover
    sympy = None  # graceful fallback; we factor via our own routines below


# ---------------------------------------------------------------------------
# Sparse vector / matrix helpers over F_p
# ---------------------------------------------------------------------------

SparseVec = dict[int, int]
SparseMat = dict[int, SparseVec]  # M[col j] = sparse column vector


def build_x1_matrix_Fp(
    standard, mon_idx, lm_to_tail, delta, is_match, match_pat, p: int
) -> SparseMat:
    """Return M as ``M[j] = {i: c}`` where x_1 * m_j  =  sum_i c * m_i.

    Internally we use the existing reduce_to_basis_Fp helper.
    """
    M: SparseMat = {}
    x1_word = ("x1",)
    for j, m_j in enumerate(standard):
        red = core.reduce_to_basis_Fp(
            [(1, x1_word + m_j)],
            mon_idx, lm_to_tail, delta, is_match, match_pat, p,
        )
        if red:
            M[j] = red
    return M


def matvec(M: SparseMat, v: SparseVec, p: int) -> SparseVec:
    """Compute (M @ v) mod p where M[j] is column j as sparse dict."""
    out: SparseVec = {}
    for j, vj in v.items():
        col = M.get(j)
        if not col:
            continue
        for i, mij in col.items():
            new = (out.get(i, 0) + vj * mij) % p
            if new == 0:
                if i in out:
                    del out[i]
            else:
                out[i] = new
    return out


def vec_add(a: SparseVec, b: SparseVec, scalar: int, p: int) -> SparseVec:
    """Return (a + scalar * b) mod p as a new sparse vector."""
    out = dict(a)
    for i, bi in b.items():
        new = (out.get(i, 0) + scalar * bi) % p
        if new == 0:
            if i in out:
                del out[i]
        else:
            out[i] = new
    return out


def vec_scale(v: SparseVec, scalar: int, p: int) -> SparseVec:
    if scalar == 0:
        return {}
    s = scalar % p
    if s == 0:
        return {}
    return {i: (vi * s) % p for i, vi in v.items()}


def vec_to_dense(v: SparseVec, N: int) -> list[int]:
    out = [0] * N
    for i, c in v.items():
        out[i] = c
    return out


# ---------------------------------------------------------------------------
# Minimal polynomial of M over F_p via Krylov + Gauss
# ---------------------------------------------------------------------------


def min_poly_on_vector(
    M: SparseMat, v0: SparseVec, N: int, p: int
) -> list[int]:
    """Smallest monic polynomial f(x) = x^d + c_{d-1} x^{d-1} + ... + c_0
    with f(M) * v0 = 0. Returned as ``[c_0, c_1, ..., c_{d-1}, 1]``.

    Algorithm: build Krylov sequence v_0, v_1 = M v_0, ..., v_k. Maintain a
    triangularised basis of {v_0, ..., v_{k-1}} (as dense rows). When v_k
    becomes a linear combination of the previous ones, read off the
    coefficients to get the dependence -> min poly on v_0.
    """
    # rows: list of (dense_vec, pivot_col, expansion). expansion[k] gives the
    # coefficients to recover this row from {v_0, ..., v_{rows_so_far}}.
    basis: list[list[int]] = []        # triangularised dense rows
    basis_pivots: list[int] = []       # pivot column for each basis row
    # We need to track how each Krylov vector v_k was reduced to a combination
    # of the basis. Let A_k denote the original v_k as a row; after reduction
    # we obtain coefficients lambda_{k,r} so that
    #    v_k - sum_r lambda_{k,r} * basis[r]_original_form = 0_after_reduction.
    # Easier: store both the original Krylov vector AND its reduction
    # coefficients against the triangular basis.

    # We store: krylov_orig[k] = dense v_k.
    # During reduction we keep coeffs telling us how to express v_k as a
    # sum of original Krylov vectors v_0, ..., v_{k-1} (modulo current
    # reduction).
    krylov_coeffs: list[list[int]] = []  # how to write reduced basis[r] in
                                          # terms of original v_0..v_r.

    v = vec_to_dense(v0, N)
    k = 0
    while True:
        # Take v_k and reduce against existing basis (triangular by pivot).
        # During reduction, we also reduce the "expression in original v_i"
        # vector ``expr`` (initially e_k).
        cur = list(v)  # copy
        expr = [0] * (k + 1)
        expr[k] = 1
        for r, row in enumerate(basis):
            pc = basis_pivots[r]
            coef = cur[pc]
            if coef == 0:
                continue
            # cur -= coef * row;  expr -= coef * krylov_coeffs[r] (zero-padded)
            for c in range(N):
                if row[c]:
                    cur[c] = (cur[c] - coef * row[c]) % p
            kr = krylov_coeffs[r]
            for i, kri in enumerate(kr):
                if kri:
                    expr[i] = (expr[i] - coef * kri) % p
        # Is cur the zero vector?
        if not any(cur):
            # Then expr[0..k] gives v_k = sum_{i<k} (-expr[i]) * v_i  (since
            # expr[k] = 1). I.e. v_k - sum_{i<k} (-expr[i] mod p ... )* v_i = 0
            # Wait: we started expr = e_k and subtracted to zero out cur.
            # If cur is zero, then v_k - (linear combo of basis rows) = 0.
            # Each basis row is itself a linear combo of earlier Krylov
            # vectors. Substituting via krylov_coeffs gives v_k as a linear
            # combo of v_0..v_{k-1}. The final ``expr`` vector encodes that
            # combo: we have
            #     0 = sum_{i<=k} expr[i] * v_i.
            # The polynomial is thus  f(x) = sum_i expr[i] * x^i; we want
            # it monic of degree k, so expr[k] should be 1 (it is, by
            # construction, since we never subtracted any basis row that had
            # nonzero contribution to coordinate k of expr). Coefficients of
            # f are expr[0..k].
            coeffs = expr[:k + 1]
            return coeffs
        # Otherwise pivot on the first nonzero entry of cur and add to basis.
        pivot = next(c for c, val in enumerate(cur) if val != 0)
        inv = pow(cur[pivot], p - 2, p)
        if inv != 1:
            for c in range(N):
                if cur[c]:
                    cur[c] = (cur[c] * inv) % p
            for i in range(k + 1):
                if expr[i]:
                    expr[i] = (expr[i] * inv) % p
        # Also reduce earlier basis rows (above row) so the basis stays
        # column-triangular and we save work later. (Optional but simplifies
        # reasoning; skip for speed; the back-substitution path above already
        # handles non-RREF.)
        basis.append(cur)
        basis_pivots.append(pivot)
        krylov_coeffs.append(expr)
        # Next Krylov vector
        v_sparse: SparseVec = {i: c for i, c in enumerate(v) if c}
        v_next_sparse = matvec(M, v_sparse, p)
        v = vec_to_dense(v_next_sparse, N)
        k += 1
        if k > N:
            raise RuntimeError("Krylov dimension exceeded N; should not happen.")


def poly_lcm(f: list[int], g: list[int], p: int) -> list[int]:
    """LCM of two monic polynomials over F_p, returned monic."""
    if not f:
        return list(g)
    if not g:
        return list(f)
    prod = poly_mul(f, g, p)
    gcd = poly_gcd(f, g, p)
    quot, _ = poly_divmod(prod, gcd, p)
    return poly_make_monic(quot, p)


def poly_mul(f: list[int], g: list[int], p: int) -> list[int]:
    out = [0] * (len(f) + len(g) - 1)
    for i, fi in enumerate(f):
        if fi == 0:
            continue
        for j, gj in enumerate(g):
            if gj == 0:
                continue
            out[i + j] = (out[i + j] + fi * gj) % p
    return out


def poly_strip(f: list[int]) -> list[int]:
    while f and f[-1] == 0:
        f.pop()
    return f


def poly_divmod(f: list[int], g: list[int], p: int) -> tuple[list[int], list[int]]:
    f = list(f)
    g = list(g)
    poly_strip(g)
    poly_strip(f)
    if not g:
        raise ZeroDivisionError
    inv_lead = pow(g[-1], p - 2, p)
    q = [0] * max(0, len(f) - len(g) + 1)
    while len(f) >= len(g) and any(f):
        deg_diff = len(f) - len(g)
        coef = (f[-1] * inv_lead) % p
        q[deg_diff] = coef
        for i, gi in enumerate(g):
            f[deg_diff + i] = (f[deg_diff + i] - coef * gi) % p
        poly_strip(f)
    return q, f


def poly_gcd(f: list[int], g: list[int], p: int) -> list[int]:
    f = list(f); g = list(g)
    poly_strip(f); poly_strip(g)
    while g:
        _, r = poly_divmod(f, g, p)
        poly_strip(r)
        f, g = g, r
    return poly_make_monic(f, p)


def poly_make_monic(f: list[int], p: int) -> list[int]:
    f = list(f)
    poly_strip(f)
    if not f:
        return f
    inv = pow(f[-1], p - 2, p)
    if inv == 1:
        return f
    return [(c * inv) % p for c in f]


def global_min_poly(M: SparseMat, N: int, p: int) -> list[int]:
    """Min poly of M on the whole space: LCM of min polys on each e_i.

    We start from e_0 (the empty-word standard monomial), get its min poly,
    then verify it annihilates every basis vector. If not, replace by LCM
    of its min poly and the offender's min poly, and retry.
    """
    # Start with e_0 (empty monomial = first standard monomial)
    mp = min_poly_on_vector(M, {0: 1}, N, p)
    for i in range(N):
        # Evaluate mp(M) * e_i; if nonzero, take min poly on e_i and LCM.
        v: SparseVec = {i: 1}
        # Horner: mp(M) v = sum_k mp[k] M^k v
        # easier: compute mp[d] v_d + mp[d-1] v_{d-1} + ... where v_k = M^k v
        d = len(mp) - 1
        # Horner from highest degree:
        # result = ((... (mp[d] * I + mp[d-1]) M + mp[d-2]) M + ... + mp[0])
        # i.e. iteratively: r = mp[d] * v; r = M*r + mp[d-1] * v; etc.
        r: SparseVec = vec_scale(v, mp[d], p)
        for k in range(d - 1, -1, -1):
            r = matvec(M, r, p)
            if mp[k]:
                r = vec_add(r, v, mp[k], p)
        if r:
            # Extend the min poly
            mp_i = min_poly_on_vector(M, {i: 1}, N, p)
            mp = poly_lcm(mp, mp_i, p)
    return mp


# ---------------------------------------------------------------------------
# Factorisation via sympy over F_p
# ---------------------------------------------------------------------------


def poly_eval(f: list[int], x: int, p: int) -> int:
    """Horner evaluation of polynomial f at x mod p."""
    r = 0
    for c in reversed(f):
        r = (r * x + c) % p
    return r


def poly_deriv(f: list[int], p: int) -> list[int]:
    if len(f) <= 1:
        return []
    return [(k * f[k]) % p for k in range(1, len(f))]


def poly_powmod(base: list[int], exponent: int, modulus: list[int], p: int) -> list[int]:
    """base^exponent mod modulus over F_p[x]."""
    result = [1]
    cur = list(base)
    _, cur = poly_divmod(cur, modulus, p)
    poly_strip(cur)
    e = exponent
    while e > 0:
        if e & 1:
            result = poly_mul(result, cur, p)
            _, result = poly_divmod(result, modulus, p)
            poly_strip(result)
        e >>= 1
        if e:
            cur = poly_mul(cur, cur, p)
            _, cur = poly_divmod(cur, modulus, p)
            poly_strip(cur)
    return result


def squarefree_part(f: list[int], p: int) -> list[int]:
    """Return f / gcd(f, f'), the squarefree part. Distinct eigenvalues =
    deg(squarefree_part) (since each distinct root contributes a degree-1
    factor of the squarefree part; over the algebraic closure, more generally
    each irreducible factor of degree d contributes d roots).
    """
    fp = poly_deriv(f, p)
    if not fp or not any(fp):
        # f' = 0: in char p, f is a polynomial in x^p; for small-degree min
        # polys at p ~ 2^31 this should never happen.
        return f
    g = poly_gcd(f, fp, p)
    q, _ = poly_divmod(f, g, p)
    return poly_make_monic(q, p)


def distinct_degree_factor(f: list[int], p: int, max_deg: int) -> dict[int, list[int]]:
    """Distinct-degree factorisation of a squarefree polynomial f over F_p.

    Returns ``{d: prod_of_irreducible_factors_of_degree_d}`` for d up to
    ``max_deg``.
    """
    result: dict[int, list[int]] = {}
    cur = list(f)
    # x mod f -- represented as [0, 1] (poly x)
    x_pow = [0, 1]  # x^{p^0} = x
    d = 0
    while len(cur) > 1:
        d += 1
        if d > max_deg:
            break
        # Raise x_pow to the p-th power mod cur
        x_pow = poly_powmod(x_pow, p, cur, p)
        # h(x) = x^{p^d} - x
        h = vec_sub(x_pow, [0, 1], p)
        g = poly_gcd(cur, h, p)
        if len(g) > 1:
            result[d] = g
            cur, _ = poly_divmod(cur, g, p)
            poly_strip(cur)
            cur = poly_make_monic(cur, p)
            # also reduce x_pow mod the new cur, so subsequent steps work
            _, x_pow = poly_divmod(x_pow, cur, p)
            poly_strip(x_pow)
    if len(cur) > 1:
        # remaining is irreducible of degree (len(cur) - 1)
        result[len(cur) - 1] = cur
    return result


def vec_sub(a: list[int], b: list[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        ai = a[i] if i < len(a) else 0
        bi = b[i] if i < len(b) else 0
        out[i] = (ai - bi) % p
    return out


def find_root_in_Fp(coeffs: list[int], p: int, q: int = 5, t: int = 7) -> Optional[int]:
    """Try a handful of obvious candidates ({q, -1/q, t^{1/a} for small a, ...})
    to identify an F_p-root of a linear or low-degree factor."""
    if len(coeffs) == 2:
        # x - r  ->  coeffs = [-r, 1]
        return (-coeffs[0]) % p
    # For higher-degree factors we don't bother finding individual roots;
    # we only need the count, which equals the degree (since the factor is
    # irreducible after squarefree-distinct-degree splitting).
    return None


def factor_over_Fp(coeffs: list[int], p: int) -> tuple[list[tuple[list[int], int]], list[int]]:
    """Factor a monic polynomial (low-to-high) over F_p.

    Returns (factors, roots_in_Fp) where each ``factor`` entry is
    ``(prod_of_irred_factors_of_some_degree_within_one_multiplicity, mult)``.

    For our purposes (counting distinct eigenvalues, verifying the conjecture),
    the squarefree factorisation + distinct-degree factorisation is enough:

      # distinct eigenvalues over alg. closure  =  deg(squarefree part of f).

    We still report the distinct-degree split of the squarefree part so the
    user can see how many F_p-roots there are vs higher-degree factors.
    """
    coeffs = poly_make_monic(list(coeffs), p)
    if len(coeffs) <= 1:
        return [], []
    # Squarefree factorisation by repeated extraction.
    factors: list[tuple[list[int], int]] = []
    cur = list(coeffs)
    mult = 1
    while len(cur) > 1:
        cur_d = poly_deriv(cur, p)
        if not any(cur_d):
            # cur is a polynomial in x^p; in our setting it shouldn't occur.
            factors.append((cur, mult))
            break
        g = poly_gcd(cur, cur_d, p)
        sf, _ = poly_divmod(cur, g, p)
        sf = poly_make_monic(sf, p)
        if len(sf) > 1:
            # break sf by distinct-degree factorisation
            dd = distinct_degree_factor(sf, p, max_deg=max(len(sf) - 1, 1))
            for d, gd in dd.items():
                # gd is product of irreducibles of degree d. If we want the
                # exact equal-degree split we'd run Cantor-Zassenhaus; for the
                # eigenvalue count we don't need to.
                factors.append((gd, mult))
        cur = poly_make_monic(g, p)
        mult += 1
    # Collect F_p roots: linear pieces of the squarefree distinct-degree split
    roots: list[int] = []
    for fac, m in factors:
        if len(fac) == 2:
            r = (-fac[0] * pow(fac[1], p - 2, p)) % p
            for _ in range(m):
                roots.append(r)
    return factors, roots


def fmt_poly(coeffs: list[int], p: int, var: str = "x") -> str:
    """Pretty-print a polynomial low-to-high, with coeffs in (-p/2, p/2]."""
    if not coeffs:
        return "0"
    parts = []
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        cc = c if c <= p // 2 else c - p
        if k == 0:
            parts.append(str(cc))
        elif k == 1:
            if cc == 1:
                parts.append(var)
            elif cc == -1:
                parts.append("-" + var)
            else:
                parts.append(f"{cc}*{var}")
        else:
            if cc == 1:
                parts.append(f"{var}^{k}")
            elif cc == -1:
                parts.append(f"-{var}^{k}")
            else:
                parts.append(f"{cc}*{var}^{k}")
    return " + ".join(parts).replace("+ -", "- ")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def infer_a_from_basename(path: str) -> Optional[int]:
    """Infer the framing parameter a from a filename like framing-n3-k2.basis."""
    base = os.path.basename(path)
    m = re.search(r"-k(\d+)\.basis", base)
    if m:
        return int(m.group(1))
    return None


def infer_n_from_basename(path: str) -> Optional[int]:
    base = os.path.basename(path)
    m = re.search(r"-n(\d+)-", base)
    return int(m.group(1)) if m else None


def evaluate_poly_on_M(
    M: SparseMat, N: int, p: int, poly: list[int]
) -> bool:
    """Return True iff ``poly(M) == 0`` on every standard basis vector."""
    d = len(poly) - 1
    for i in range(N):
        v: SparseVec = {i: 1}
        r = vec_scale(v, poly[d], p)
        for k in range(d - 1, -1, -1):
            r = matvec(M, r, p)
            if poly[k]:
                r = vec_add(r, v, poly[k], p)
        if r:
            return False
    return True


def annihilator_check(
    M: SparseMat, N: int, p: int, a: int, q: int = 5, t: int = 7
) -> dict:
    """Test several candidate annihilators of sigma_1 = x_1 on the algebra.

    Reports for each candidate whether it acts as zero on every standard
    basis vector.
    """
    q_mod = q % p
    q_inv = pow(q_mod, p - 2, p)
    out = {}
    # Original conjecture: (x - q)(x + 1/q)(x^a - t)
    hecke_orig = [(-1) % p, (q_inv - q_mod) % p, 1]
    framing = [0] * (a + 1)
    framing[0] = (-t) % p
    framing[a] = 1
    prod_orig = poly_mul(hecke_orig, framing, p)
    out["(x-q)(x+1/q)(x^a-t)"] = (evaluate_poly_on_M(M, N, p, prod_orig),
                                  len(prod_orig) - 1)
    # Refined empirical: (x - 1)(x + q)(x^a - t)
    hecke_refined = [(-q_mod) % p, (q_mod - 1) % p, 1]
    prod_refined = poly_mul(hecke_refined, framing, p)
    out["(x-1)(x+q)(x^a-t)"] = (evaluate_poly_on_M(M, N, p, prod_refined),
                                len(prod_refined) - 1)
    # Same with extra (x) factor
    prod_with_x = poly_mul(prod_refined, [0, 1], p)
    out["x * (x-1)(x+q)(x^a-t)"] = (evaluate_poly_on_M(M, N, p, prod_with_x),
                                    len(prod_with_x) - 1)
    return out


def process(basis_path: str) -> dict:
    print(f"\n=== {basis_path} ===")
    a = infer_a_from_basename(basis_path)
    n = infer_n_from_basename(basis_path)
    polys = core.parse_basis_Fp(basis_path)
    # Look up p from the .ms file if available; default to 2^31 - 1.
    p = 2**31 - 1
    ms_path = basis_path.replace(".basis", ".ms")
    if os.path.exists(ms_path):
        with open(ms_path) as f:
            lines = [ln.strip() for ln in f.readlines()]
            if len(lines) >= 2 and lines[1].isdigit():
                p = int(lines[1])
    print(f"  prime p = {p}, framing parameter a = {a}, n = {n}")
    lms = [lm for lm, _ in polys]
    variables = core.collect_variables(polys)
    print(f"  variables = {variables}, |basis| = {len(polys)}")

    delta, is_match, match_pat, nstates = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        delta, is_match, variables, max_total=20000
    )
    if standard is None:
        print("  dim > 20000 (truncated)")
        return {"path": basis_path, "a": a, "n": n, "dim": ">20000"}
    N = len(standard)
    print(f"  dim = {N}")
    mon_idx = {m: i for i, m in enumerate(standard)}
    lm_to_tail = core.lm_to_tail_Fp(polys, p)

    M = build_x1_matrix_Fp(
        standard, mon_idx, lm_to_tail, delta, is_match, match_pat, p
    )
    nnz = sum(len(col) for col in M.values())
    print(f"  M (left-mult by x_1): {len(M)} nonzero columns, {nnz} nonzero entries")

    mp = global_min_poly(M, N, p)
    deg = len(mp) - 1
    print(f"  min poly degree = {deg}")
    print(f"  min poly        = {fmt_poly(mp, p)}")

    factors, roots = factor_over_Fp(mp, p)
    distinct_roots = sorted(set(roots))
    sf = squarefree_part(mp, p)
    n_distinct_total = len(sf) - 1  # # distinct roots over alg. closure
    print(f"  squarefree-part degree (# distinct eigenvalues over alg. closure) = "
          f"{n_distinct_total}")
    print(f"  squarefree-distinct-degree split:")
    for fac, mult in factors:
        # Each (fac, mult) is the PRODUCT of all irreducibles of some shared
        # degree d, raised to multiplicity ``mult``. Recover d via the
        # distinct-degree pass on ``fac`` (a recompute on a small poly).
        if len(fac) <= 1:
            continue
        dd = distinct_degree_factor(fac, p, max_deg=max(len(fac) - 1, 1))
        for d_irr, fac_d in dd.items():
            count = (len(fac_d) - 1) // d_irr
            kind = "linear" if d_irr == 1 else f"degree-{d_irr}-irreducible"
            print(f"    {count} x ({kind}) at multiplicity {mult}: "
                  f"prod = ({fmt_poly(fac_d, p)})")
    # Check specific candidates for being F_p roots of the min poly
    q_val = 5 % p
    q_inv = pow(q_val, p - 2, p)
    neg_q_inv = (-q_inv) % p
    t_val = 7 % p
    candidates = {q_val: "q=5", neg_q_inv: "-q^{-1}=-1/5", t_val: "t=7"}
    if a is not None:
        # t^{1/a} candidates: any r with r^a = t. Try small integers.
        for r in range(-20, 21):
            rr = r % p
            if pow(rr, a, p) == t_val and rr not in candidates:
                candidates[rr] = f"t^{{1/{a}}}={r}"
    found = []
    for c, label in candidates.items():
        val = poly_eval(mp, c, p)
        if val == 0:
            found.append((c if c <= p // 2 else c - p, label))
    print(f"  F_p roots of the min poly among predicted candidates:")
    for c, label in found:
        print(f"    {label}: {c}  (root)")
    print(f"  # distinct F_p roots found among predicted candidates = {len(found)}")
    if a is not None:
        print(f"  total distinct eigenvalues (over alg. closure) = {n_distinct_total}")
        print(f"  expected per conjecture (a + 2)               = {a + 2}")
        # Annihilator check on multiple candidate polynomials
        ann_results = annihilator_check(M, N, p, a)
        print(f"  candidate annihilators (does poly(M) = 0 on the algebra?):")
        for name, (ok, deg_ann) in ann_results.items():
            marker = "ZERO " if ok else "NONZERO"
            print(f"    {marker} (deg {deg_ann}): {name}")
        orig_ok = ann_results["(x-q)(x+1/q)(x^a-t)"][0]
        refined_ok = ann_results["(x-1)(x+q)(x^a-t)"][0]
        with_x_ok = ann_results["x * (x-1)(x+q)(x^a-t)"][0]
        return {
            "path": basis_path, "a": a, "n": n, "dim": N,
            "min_poly_deg": deg,
            "n_distinct_eigenvalues": n_distinct_total,
            "expected": a + 2,
            "factors": [(fmt_poly(fac, p), mult) for fac, mult in factors],
            "annihilator_orig_zero": orig_ok,
            "annihilator_refined_zero": refined_ok,
            "annihilator_x_refined_zero": with_x_ok,
        }
    return {"path": basis_path, "a": a, "n": n, "dim": N,
            "min_poly_deg": deg,
            "factors": [(fmt_poly(fac, p), mult) for fac, mult in factors]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", nargs="+", help="one or more .basis files")
    ap.add_argument("--table", action="store_true",
                    help="print summary table at the end")
    args = ap.parse_args()

    results = []
    for path in args.basis:
        try:
            res = process(path)
            results.append(res)
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({"path": path, "error": str(e)})

    if args.table or len(args.basis) > 1:
        print("\n\n=== SUMMARY TABLE ===")
        hdr = (f"{'file':32s} {'n':>2s} {'a':>2s} {'dim':>5s} "
               f"{'minpoly_deg':>11s} {'n_eig':>5s} {'a+2':>4s} "
               f"{'(x-q)(x+1/q)(x^a-t)':>20s} {'(x-1)(x+q)(x^a-t)':>18s} "
               f"{'x*(x-1)(x+q)(x^a-t)':>20s}   factors")
        print(hdr)
        print("-" * len(hdr))
        for r in results:
            if "error" in r:
                print(f"{os.path.basename(r['path']):32s} ERROR: {r['error']}")
                continue
            facs = "; ".join(f"({f})^{m}" for f, m in r.get("factors", []))
            flag = ""
            if r.get("a") is not None and r.get("n_distinct_eigenvalues") is not None:
                if r["n_distinct_eigenvalues"] >= 5 and r["a"] >= 3:
                    flag = "  <-- 5+ eigenvalues (a >= 3)"
            print(
                f"{os.path.basename(r['path']):32s} "
                f"{r.get('n','?'):>2} "
                f"{r.get('a','?'):>2} "
                f"{r.get('dim','?'):>5} "
                f"{r.get('min_poly_deg','?'):>11} "
                f"{r.get('n_distinct_eigenvalues','?'):>5} "
                f"{r.get('expected','?'):>4} "
                f"{str(r.get('annihilator_orig_zero','?')):>20} "
                f"{str(r.get('annihilator_refined_zero','?')):>18} "
                f"{str(r.get('annihilator_x_refined_zero','?')):>20}   "
                f"{facs}{flag}"
            )


if __name__ == "__main__":
    main()
