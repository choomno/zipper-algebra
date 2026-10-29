#!/usr/bin/env python3
r"""Search for nilpotent two-sided ideal witnesses in supplied reduction models.

For the square-zero family, the historical candidate generator is
w_j=(sigma_j-t)X_{2,j}. This recipe calculates ideal products at an explicitly
chosen parameter pair. It uses exact rational arithmetic, but the supplied
solver basis is not certified complete here. Conclusions about the presented
algebra require that missing basis evidence; conclusions over Q(q,t) require
additional generic identities. A nonzero nilpotent element alone does not prove
nonsemisimplicity. A verified nonzero nilpotent two-sided ideal would do so.

The algorithms below retain the historical witness search. Run with
python -m experiments.variants.square_zero_ideal_nilpotence.
"""
from __future__ import annotations

import argparse
import os
import subprocess
from fractions import Fraction


from zipper_algebra import core, Presentation
from experiments.output import (capture_experiment_provenance, tag_solver_input,
                                tag_solver_basis)


# ---------------------------------------------------------------------------
# Polynomial primitives over Q in the free algebra (sparse {mon: coef})
# ---------------------------------------------------------------------------

Mon = tuple[str, ...]
Poly = dict[Mon, Fraction]  # sparse free-algebra polynomial


def poly_add(p: Poly, q: Poly) -> Poly:
    out = dict(p)
    for m, c in q.items():
        new = out.get(m, Fraction(0)) + c
        if new == 0 and m in out:
            del out[m]
        else:
            out[m] = new
    return out


def poly_mul(p: Poly, q: Poly) -> Poly:
    out: Poly = {}
    for m1, c1 in p.items():
        for m2, c2 in q.items():
            m = m1 + m2
            new = out.get(m, Fraction(0)) + c1 * c2
            if new == 0 and m in out:
                del out[m]
            else:
                out[m] = new
    return out


def poly_to_terms(p: Poly) -> list[tuple[Fraction, Mon]]:
    return [(c, m) for m, c in p.items()]


# ---------------------------------------------------------------------------
# Reduction wrapper
# ---------------------------------------------------------------------------


def reduce_poly_Q(
    p: Poly,
    mon_idx: dict[Mon, int],
    lm_to_tail: dict[Mon, list[tuple[Fraction, Mon]]],
    delta,
    is_match,
    match_pat,
) -> dict[int, Fraction]:
    """Reduce a free-algebra polynomial to its standard normal form
    expressed as a sparse vector ``{basis_idx: coef}`` over Q."""
    if not p:
        return {}
    return core.reduce_to_basis_Q(
        poly_to_terms(p), mon_idx, lm_to_tail, delta, is_match, match_pat,
    )


def vec_to_poly(vec: dict[int, Fraction], standard: list[Mon]) -> Poly:
    """Convert a sparse coordinate vector back into a free-algebra polynomial
    (each basis index maps to its standard monomial)."""
    return {standard[i]: c for i, c in vec.items() if c != 0}


# ---------------------------------------------------------------------------
# Witness w_j  =  (sigma_j - t) * X_{2,j}
#
# In the free algebra:
#   sigma_j  =  x_j                    (the positive braid generator)
#   S_2_j   =  q*y_j + (1 - q) - x_j   (the level-2 slider)
#   w_j     =  (x_j - t) * S_2_j
#           =   q*x_j*y_j  +  (1-q)*x_j  -  x_j*x_j
#               - t*q*y_j  - t*(1-q)    +  t*x_j
# ---------------------------------------------------------------------------


def witness_w(j: int, q: Fraction, t: Fraction) -> Poly:
    xj = f"x{j}"
    yj = f"y{j}"
    one_minus_q = Fraction(1) - q
    p: Poly = {}
    # (x_j - t) * (q*y_j + (1-q) - x_j)
    # = q*x_j*y_j + (1-q)*x_j - x_j*x_j
    #   - t*q*y_j - t*(1-q) + t*x_j
    def add(m: Mon, c: Fraction) -> None:
        new = p.get(m, Fraction(0)) + c
        if new == 0 and m in p:
            del p[m]
        elif new != 0:
            p[m] = new
    add((xj, yj), q)
    add((xj,), one_minus_q)
    add((xj, xj), Fraction(-1))
    add((yj,), -t * q)
    add((), -t * one_minus_q)
    add((xj,), t)
    return p


# ---------------------------------------------------------------------------
# Vector-space subspace tracking (row-reduced)
# ---------------------------------------------------------------------------


class Subspace:
    """A subspace of Q^N tracked as a list of pivot-reduced sparse vectors.
    Rows are dicts ``{column_index: Fraction}``. Pivot column for row i is
    ``pivots[i]`` (the smallest column index with a non-zero entry, after
    elimination)."""

    def __init__(self) -> None:
        self.rows: list[dict[int, Fraction]] = []
        self.pivots: list[int] = []

    def dim(self) -> int:
        return len(self.rows)

    def reduce(self, v: dict[int, Fraction]) -> dict[int, Fraction]:
        """Return v reduced modulo the current subspace (in place: v is copied)."""
        v = dict(v)
        for row, piv in zip(self.rows, self.pivots):
            if piv in v:
                c = v[piv] / row[piv]
                if c == 0:
                    continue
                for col, rc in row.items():
                    new = v.get(col, Fraction(0)) - c * rc
                    if new == 0 and col in v:
                        del v[col]
                    elif new != 0:
                        v[col] = new
        return v

    def add(self, v: dict[int, Fraction]) -> bool:
        """Try to add v; return True iff dim increased."""
        v = self.reduce(v)
        if not v:
            return False
        piv = min(v.keys())
        # Re-eliminate above-pivot in stored rows? Not necessary -- we only need
        # span. Keep upper-triangular structure by inserting in order.
        self.rows.append(v)
        self.pivots.append(piv)
        return True


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------


def fmt_poly(p: Poly, max_terms: int = 8) -> str:
    if not p:
        return "0"
    items = sorted(p.items(), key=lambda mc: (len(mc[0]), mc[0]))
    parts: list[str] = []
    for m, c in items[:max_terms]:
        ms = "*".join(m) if m else "1"
        parts.append(f"{c}*{ms}")
    if len(items) > max_terms:
        parts.append(f"... (+{len(items) - max_terms} more)")
    return " + ".join(parts)


def fmt_vec(vec: dict[int, Fraction], standard: list[Mon], max_terms: int = 8) -> str:
    if not vec:
        return "0"
    items = sorted(vec.items(), key=lambda ic: (len(standard[ic[0]]), standard[ic[0]]))
    parts: list[str] = []
    for i, c in items[:max_terms]:
        m = standard[i]
        ms = "*".join(m) if m else "1"
        parts.append(f"{c}*{ms}")
    if len(items) > max_terms:
        parts.append(f"... (+{len(items) - max_terms} more)")
    return " + ".join(parts)


def ensure_basis(n: int, q: int, t: int, basis_path: str | None) -> str:
    """If basis_path is given and exists, use it. Otherwise generate it from
    temporary inputs using Presentation.emit_ms and f4ncgb."""
    if basis_path and os.path.exists(basis_path):
        return basis_path
    if basis_path is None:
        basis_path = f"/tmp/square-zero-n{n}-q{q}-t{t}-c0.basis"
    if os.path.exists(basis_path):
        return basis_path
    producer = capture_experiment_provenance(__file__)
    ms_path = basis_path.replace(".basis", ".ms")
    print(f"[setup] generating {ms_path} (n={n}, q={q}, t={t}, square-zero)", flush=True)
    Presentation(n, untwist="square-zero").emit_ms(ms_path, q, t, 0, provenance=False)
    inputs = tag_solver_input(ms_path, producer=producer)
    print(f"[setup] running f4ncgb -> {basis_path}", flush=True)
    subprocess.run(
        [
            os.environ.get("F4NCGB", "f4ncgb"),
            "-m", "500",
            "-o", basis_path,
            "-i", ms_path,
        ],
        check=True,
    )
    tag_solver_basis(basis_path, producer=producer, inputs=inputs)
    return basis_path


# ---------------------------------------------------------------------------
# Two-sided ideal power computation
# ---------------------------------------------------------------------------


def build_LR_generators(
    variables: list[str],
    standard: list[Mon],
    mon_idx: dict[Mon, int],
    lm_to_tail,
    delta,
    is_match,
    match_pat,
) -> tuple[
    dict[str, dict[int, dict[int, Fraction]]],
    dict[str, dict[int, dict[int, Fraction]]],
]:
    """Build, for each generator v, the sparse matrices of left- and right-
    multiplication by v.

    ``gen_L[v][j]`` is the sparse column vector
        reduce( v * standard[j] )  in standard form, as {k: coef}.
    Then left-multiplication of a vector u by v is
        (Lv u)[k]  =  sum_j u[j] * gen_L[v][j][k].
    Symmetric for gen_R[v][j] = reduce(standard[j] * v).
    """
    gen_L: dict[str, dict[int, dict[int, Fraction]]] = {}
    gen_R: dict[str, dict[int, dict[int, Fraction]]] = {}
    for v in variables:
        cl: dict[int, dict[int, Fraction]] = {}
        cr: dict[int, dict[int, Fraction]] = {}
        for j, mj in enumerate(standard):
            red_l = core.reduce_to_basis_Q(
                [(Fraction(1), (v,) + mj)],
                mon_idx, lm_to_tail, delta, is_match, match_pat,
            )
            if red_l:
                cl[j] = red_l
            red_r = core.reduce_to_basis_Q(
                [(Fraction(1), mj + (v,))],
                mon_idx, lm_to_tail, delta, is_match, match_pat,
            )
            if red_r:
                cr[j] = red_r
        gen_L[v] = cl
        gen_R[v] = cr
    return gen_L, gen_R


def apply_L(gen_L, v: str, vec: dict[int, Fraction]) -> dict[int, Fraction]:
    """Compute v * vec."""
    out: dict[int, Fraction] = {}
    col = gen_L[v]
    for j, cj in vec.items():
        cj_vec = col.get(j)
        if cj_vec is None:
            continue
        for k, ck in cj_vec.items():
            new = out.get(k, Fraction(0)) + cj * ck
            if new == 0 and k in out:
                del out[k]
            elif new != 0:
                out[k] = new
    return out


def apply_R(gen_R, v: str, vec: dict[int, Fraction]) -> dict[int, Fraction]:
    """Compute vec * v."""
    out: dict[int, Fraction] = {}
    col = gen_R[v]
    for j, cj in vec.items():
        cj_vec = col.get(j)
        if cj_vec is None:
            continue
        for k, ck in cj_vec.items():
            new = out.get(k, Fraction(0)) + cj * ck
            if new == 0 and k in out:
                del out[k]
            elif new != 0:
                out[k] = new
    return out


def two_sided_ideal_basis(
    w_vec: dict[int, Fraction],
    variables: list[str],
    gen_L,
    gen_R,
) -> Subspace:
    """Compute V_1 = <w> as a two-sided ideal of A, given w_vec already in
    standard form.

    Algorithm: BFS closure. Start with the singleton {w_vec}. Maintain a
    Subspace V_1. For each new pivot vector enqueued, apply each L[v] and
    each R[v]; if the result reduces non-trivially modulo V_1, add it and
    enqueue. Terminates because V_1 is a subspace of the finite-dimensional
    algebra A.

    This is the fast path: instead of N^2 (a, b) pairs we do at most
    dim(V_1) * 2 * |variables| variable-applications, each costing O(supp).
    """
    V1 = Subspace()
    from collections import deque
    if not w_vec:
        return V1
    V1.add(dict(w_vec))
    queue: deque[dict[int, Fraction]] = deque([dict(w_vec)])
    while queue:
        u = queue.popleft()
        for v in variables:
            out_L = apply_L(gen_L, v, u)
            if out_L and V1.add(out_L):
                queue.append(out_L)
            out_R = apply_R(gen_R, v, u)
            if out_R and V1.add(out_R):
                queue.append(out_R)
    return V1


def build_monomial_left_op(
    a: int,
    standard: list[Mon],
    gen_L,
    cache: dict[int, dict[int, dict[int, Fraction]]],
) -> dict[int, dict[int, Fraction]]:
    """Return the cached column dict  M[j] = standard form of  standard[a] *
    standard[j]  for all j with non-zero result.

    Computed by walking the letters of standard[a] right-to-left and applying
    the appropriate L[v] to each basis column. We populate the cache lazily.
    """
    if a in cache:
        return cache[a]
    ma = standard[a]
    if not ma:
        # standard[a] = 1; left-mult is the identity on each basis vector.
        ident: dict[int, dict[int, Fraction]] = {
            j: {j: Fraction(1)} for j in range(len(standard))
        }
        cache[a] = ident
        return ident
    # Start from the identity (j -> {j: 1}) and apply L[letter] in reverse.
    cur: dict[int, dict[int, Fraction]] = {j: {j: Fraction(1)} for j in range(len(standard))}
    for letter in reversed(ma):
        nxt: dict[int, dict[int, Fraction]] = {}
        col = gen_L[letter]
        for j, jvec in cur.items():
            out: dict[int, Fraction] = {}
            for jj, cjj in jvec.items():
                cj_vec = col.get(jj)
                if cj_vec is None:
                    continue
                for k, ck in cj_vec.items():
                    new = out.get(k, Fraction(0)) + cjj * ck
                    if new == 0 and k in out:
                        del out[k]
                    elif new != 0:
                        out[k] = new
            if out:
                nxt[j] = out
        cur = nxt
    cache[a] = cur
    return cur


def algebra_product_cached(
    u_vec: dict[int, Fraction],
    v_vec: dict[int, Fraction],
    standard: list[Mon],
    gen_L,
    left_op_cache: dict[int, dict[int, dict[int, Fraction]]],
) -> dict[int, Fraction]:
    """Multiply coefficient vectors using cached left-multiplication operators
    for standard monomials. Each (a, j) lookup is an O(1) sparse vector,
    accumulated as u[a] * v[j] * M_a[j]."""
    out: dict[int, Fraction] = {}
    for a, ca in u_vec.items():
        M_a = build_monomial_left_op(a, standard, gen_L, left_op_cache)
        for j, cj in v_vec.items():
            col = M_a.get(j)
            if not col:
                continue
            scale = ca * cj
            if scale == 0:
                continue
            for k, ck in col.items():
                new = out.get(k, Fraction(0)) + scale * ck
                if new == 0 and k in out:
                    del out[k]
                elif new != 0:
                    out[k] = new
    return out


def compute_w_power_subspaces(
    w_poly: Poly,
    standard: list[Mon],
    mon_idx: dict[Mon, int],
    lm_to_tail,
    delta,
    is_match,
    match_pat,
    variables: list[str],
    gen_L,
    gen_R,
    left_op_cache: dict[int, dict[int, dict[int, Fraction]]],
    max_power: int = 4,
    log=lambda s: None,
    V1_precomputed: Subspace | None = None,
) -> list[Subspace]:
    """For a witness polynomial w (as a free-algebra Poly), compute
    successively the standard-form span of

        V_1 = span{ a * w * b : a, b in standard basis }   (= <w> as two-sided ideal)
        V_2 = span{ a * w * m * w * b : a, b in standard, m in standard }
            = span{ A * (w * m * w) * B : ... }, equivalently V_1 * V_1
        ...
        V_k = V_1^k = <w>^k

    A cleaner recurrence: <w>^k is the two-sided ideal generated by all
    products of k elements of <w>. Equivalently (since <w> = span_{l,r}
    l*w*r), <w>^k = span of  a_0 * w * m_1 * w * m_2 * w * ... * m_{k-1}
    * w * b  over standard a_0, m_i, b. Since each such product is
    expressible (working in the algebra A) as an algebra product
    iteratively, we compute:

        V_1 vectors = { reduce(a * w * b) for a, b in standard }.
        For k >= 2, V_k is spanned by
            { reduce(v * c) for v in V_{k-1}-spanning-set,
                                c in V_1-spanning-set }
        ... which is the same as  V_{k-1} * V_1 (as a left-module
        product in the algebra).

    But for the Wedderburn argument we need only: is <w>^k = 0?
    Equivalently, is (V_1)^k = 0 in A? Since <w>^k = (algebra ideal
    product) = sum of products v_1 v_2 ... v_k of elements of <w>, and
    since V_1 is itself a two-sided ideal in A, we have <w>^k = V_1^k as
    a subspace of A. So:

        V_k = span{ u * v : u in V_{k-1}, v in V_1 }.

    We track V_k as a Subspace and return the list [V_1, V_2, ..., V_K]
    where K is the smallest power with V_K = 0 (or max_power if no such
    power is found within the budget).
    """
    # --- step 1: reduce w to standard form ---
    w_vec = reduce_poly_Q(w_poly, mon_idx, lm_to_tail, delta, is_match, match_pat)
    log(f"  w (reduced) supp = {len(w_vec)}   {fmt_vec(w_vec, standard)}")
    if not w_vec:
        log("  w == 0; the witness is trivial, no Wedderburn argument possible.")
        return []

    # --- step 2: build V_1 = <w> as a two-sided ideal in A ---
    if V1_precomputed is not None:
        V1 = V1_precomputed
        log(f"  V_1 has dim {V1.dim()} (cached)")
    else:
        log('  building V_1 = <w> as two-sided ideal of A via BFS closure under L_v, R_v')
        V1 = two_sided_ideal_basis(w_vec, variables, gen_L, gen_R)
        log(f"  V_1 has dim {V1.dim()}")

    powers: list[Subspace] = [V1]
    if V1.dim() == 0:
        return powers

    # --- step 3: fast check whether V_1^2 = 0 via the w * m * w spanning set ---
    # V_1 = A * w * A  is spanned by  { a * w * b : a, b standard }, so
    #   (V_1)^2  =  span {  (a_1 * w * b_1) * (a_2 * w * b_2)  }
    #           =  span {  a_1 * w * (b_1 * a_2)_reduced * w * b_2  }
    #           =  span {  a_1 * (w * m * w) * b_2  : m standard, a_1, b_2 standard }
    #           =  sum_{m standard}  A * (w * m * w) * A .
    # So V_1^2 = 0  <=>  w * m * w = 0 in A  for every standard m.
    # This needs only N reductions (one per standard m), not dim(V_1)^2.
    if max_power >= 2:
        log(f"  checking V_2 = 0 via:  w * m * w  for all {len(standard)} standard m")
        V2 = Subspace()
        # w as a free-algebra polynomial in standard form
        w_poly_std = vec_to_poly(w_vec, standard)
        any_nonzero = 0
        for m_idx, m in enumerate(standard):
            # Build  w * m * w  as a free-algebra polynomial.
            # For m = () (the unit) this is w * w.
            if m:
                m_poly: Poly = {m: Fraction(1)}
                prod1 = poly_mul(w_poly_std, m_poly)
            else:
                prod1 = dict(w_poly_std)
            prod2 = poly_mul(prod1, w_poly_std)
            vec = reduce_poly_Q(
                prod2, mon_idx, lm_to_tail, delta, is_match, match_pat,
            )
            if vec:
                any_nonzero += 1
                V2.add(vec)
        log(f"  non-zero  w * m * w  for {any_nonzero} of {len(standard)} standard m;"
            f"  dim V_2 = {V2.dim()}")
        powers.append(V2)
        if V2.dim() == 0:
            return powers
        # fall through to higher powers using the cached operator path
    for k in range(3, max_power + 1):
        Vprev = powers[-1]
        if Vprev.dim() == 0:
            break
        log(f"  computing V_{k} = V_{k-1} * V_1 ({Vprev.dim()} x {V1.dim()} pair products)")
        Vk = Subspace()
        for u in Vprev.rows:
            for vv in V1.rows:
                vec = algebra_product_cached(
                    u, vv, standard, gen_L, left_op_cache,
                )
                if vec:
                    Vk.add(vec)
        log(f"  V_{k} has dim {Vk.dim()}")
        powers.append(Vk)
        if Vk.dim() == 0:
            break
    return powers


def compute_full_radical_candidate(
    V1_list: list[Subspace],
    log=lambda s: None,
) -> Subspace:
    """Given precomputed two-sided ideals  V_1[j] = <w_j>, return their sum
    as a Subspace (= the two-sided ideal generated by all w_j)."""
    I = Subspace()
    for j_idx, V1 in enumerate(V1_list):
        for row in V1.rows:
            I.add(dict(row))
        log(f"  after adding <w_{j_idx + 1}>  (dim {V1.dim()}):  running sum dim = {I.dim()}")
    return I


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, required=True, help="strand count")
    ap.add_argument("--q", type=int, default=5, help="q specialisation")
    ap.add_argument("--t", type=int, default=7, help="t specialisation")
    ap.add_argument("--basis", default=None,
                    help="path to f4ncgb .basis (auto-generates if absent)")
    ap.add_argument("--max-power", type=int, default=4,
                    help="largest k for which to test  <w>^k = 0")
    ap.add_argument("--max-dim", type=int, default=400,
                    help="abort if standard-monomial enumeration exceeds this")
    ap.add_argument("--also-radical", action="store_true",
                    help="also compute  dim <w_1, ..., w_{n-1}>  for comparison "
                         "with dim J(A)")
    ap.add_argument("--js", type=str, default=None,
                    help="comma-separated subset of j values to test (default: all)")
    args = ap.parse_args()

    def log(s: str) -> None:
        print(s, flush=True)

    basis_path = ensure_basis(args.n, args.q, args.t, args.basis)
    log(f"[parse] reading basis from {basis_path}")
    polys = core.parse_basis_Q(basis_path)
    variables = Presentation(args.n).var_names()
    log(f"[vars]  {variables}")
    log(f"[basis] {len(polys)} GB polynomials")

    lms = [lm for lm, _ in polys]
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        delta, is_match, variables, max_total=args.max_dim,
    )
    if standard is None:
        log(f"[abort] dim > {args.max_dim}")
        return 1
    N = len(standard)
    log(f"[dim]   N = dim_Q A = {N}")
    mon_idx = {m: i for i, m in enumerate(standard)}
    lm_to_tail = core.lm_to_tail_Q(polys)

    log('[ops]   building left/right generator matrices for each variable')
    gen_L, gen_R = build_LR_generators(
        variables, standard, mon_idx, lm_to_tail, delta, is_match, match_pat,
    )
    log('[ops]   L/R generator matrices ready')
    left_op_cache: dict[int, dict[int, dict[int, Fraction]]] = {}

    q = Fraction(args.q)
    t = Fraction(args.t)

    # Per-j analysis
    if args.js:
        js = [int(x) for x in args.js.split(",")]
    else:
        js = list(range(1, args.n))

    overall_proof = False
    smallest_index: dict[int, int] = {}  # j -> smallest k with <w_j>^k = 0
    w_polys_for_radical: list[Poly] = []
    V1_cache: list[Subspace] = []

    for j in js:
        log("")
        log('=' * 72)
        log(f"WITNESS j = {j}:   w_{j}  =  (sigma_{j} - t) * X_{{2,{j}}}")
        log('=' * 72)
        w_poly = witness_w(j, q, t)
        w_polys_for_radical.append(w_poly)
        log(f"  free-algebra form  w_{j}  =  {fmt_poly(w_poly)}")

        # Step 1: reduce w and check non-zero.
        w_vec = reduce_poly_Q(
            w_poly, mon_idx, lm_to_tail, delta, is_match, match_pat,
        )
        log(f"  reduced w_{j}      =  {fmt_vec(w_vec, standard)}")
        log(f"  supp(w_{j})        =  {len(w_vec)}")
        if not w_vec:
            log(f"  [warn] w_{j} reduces to 0 in A; no Wedderburn argument from this j")
            continue

        # Step 2: directly verify w^2 = 0 (the N1 identity).
        w_vec_poly = vec_to_poly(w_vec, standard)
        w_sq_prod = poly_mul(w_vec_poly, w_vec_poly)
        w_sq_vec = reduce_poly_Q(
            w_sq_prod, mon_idx, lm_to_tail, delta, is_match, match_pat,
        )
        log(f"  reduced w_{j}^2    =  {fmt_vec(w_sq_vec, standard)}    (supp = {len(w_sq_vec)})")
        if not w_sq_vec:
            log(f"  [ok] w_{j}^2 = 0  (N1 identity holds; w_{j} is an element-level nilpotent)")
        else:
            log(f"  [info] w_{j}^2 != 0; check the N1 writeup")

        # Step 3: compute the two-sided-ideal-power tower V_1, V_2, ...
        log("")
        log(f"  -- computing <w_{j}>^k as two-sided ideal of A --")
        powers = compute_w_power_subspaces(
            w_poly, standard, mon_idx, lm_to_tail, delta, is_match, match_pat,
            variables, gen_L, gen_R, left_op_cache,
            max_power=args.max_power, log=log,
        )
        if not powers:
            log('  [skip] no powers computed')
            continue
        V1_cache.append(powers[0])
        dims = [p.dim() for p in powers]
        log(f"  dim <w_{j}>^k for k = 1..{len(powers)}  =  {dims}")
        # smallest k with V_k = 0
        smallest = None
        for k, d in enumerate(dims, start=1):
            if d == 0:
                smallest = k
                break
        if smallest is not None:
            smallest_index[j] = smallest
            log("")
            log(f"  In this reduction model, <w_{j}> has dimension {dims[0]} and power {smallest} vanishes.")
            log("  A nonzero nilpotent two-sided ideal implies nonsemisimplicity only after")
            log("  the reduction model is justified as the intended specialized algebra.")
            overall_proof = True
        else:
            log(f"  [warn] <w_{j}>^k != 0 for all k <= {args.max_power}; bump --max-power")

    # Bonus radical check
    if args.also_radical and V1_cache:
        log("")
        log("=" * 72)
        log("BONUS  --  dim <w_1, ..., w_{n-1}> as two-sided ideal of A")
        log("=" * 72)
        I = compute_full_radical_candidate(V1_cache, log=log)
        log(f"  dim_Q <w_1, ..., w_{{n-1}}>  =  {I.dim()}")
        # Typical comparison points
        if args.n == 3:
            log('  prior trace-form dim J(A)    =  9')
            if I.dim() == 9:
                log("  [MATCH] the w_j generate the full radical; Wedderburn-Mal'cev")
                log('          splitting candidate identified.')
        elif args.n == 4:
            log('  prior trace-form dim J(A)    =  81 (per writeup)')
            if I.dim() == 81:
                log('  [MATCH] the w_j generate the full radical.')

    log("")
    log("=" * 72)
    log("REDUCTION-MODEL OBSERVATION")
    log("=" * 72)
    if overall_proof:
        log(f"  Nilpotent ideal witness found in the supplied reduction model for n={args.n}.")
        log(f"  Per-witness nilpotence indices:  {smallest_index}")
        log("  Interpretation for the presented specialized algebra requires a certified basis.")
        log("  This does not establish a generic conclusion over Q(q,t).")
        return 0
    else:
        log(f"  No witness with <w_j>^k = 0 (k <= {args.max_power}) was found.")
        log('  Increase --max-power or check the square-zero presentation.')
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
