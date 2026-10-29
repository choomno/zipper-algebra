"""Sampled reduction-map and ideal-product observations.

These computations use uncertified basis candidates. Kernel dimensions,
quotient dimensions and radical identifications need separate basis and
structural evidence before they can be inferred from the observed matrices.
"""
from __future__ import annotations
from .spec import Spec
from .pointalg import PointAlgebra
from zipper_algebra.linear import validate_field_characteristic


def _nullspace_cols_Fp(M, ncol, p):
    """Right null space of a matrix given as a dict {(row,col): val}; returns list of basis vectors
    (each a dict col->val) over F_p."""
    validate_field_characteristic(p)
    if p == 0:
        raise ValueError("a prime characteristic is required")
    # dense rows keyed by row index
    rows = {}
    for (r, c), v in M.items():
        rows.setdefault(r, {})[c] = v % p
    A = [rows.get(r, {}) for r in range(max((r for r, _ in M), default=-1) + 1)]
    # Gaussian elimination (row echelon), track pivot columns
    pivot_col = {}
    r = 0
    dense = [[row.get(c, 0) for c in range(ncol)] for row in A]
    nr = len(dense)
    for c in range(ncol):
        piv = next((i for i in range(r, nr) if dense[i][c] % p), None)
        if piv is None:
            continue
        dense[r], dense[piv] = dense[piv], dense[r]
        inv = pow(dense[r][c], p - 2, p)
        dense[r] = [(x * inv) % p for x in dense[r]]
        for i in range(nr):
            if i != r and dense[i][c] % p:
                f = dense[i][c]
                dense[i] = [(dense[i][j] - f * dense[r][j]) % p for j in range(ncol)]
        pivot_col[c] = r
        r += 1
        if r == nr:
            break
    free = [c for c in range(ncol) if c not in pivot_col]
    basis = []
    for fc in free:
        vec = {fc: 1}
        for c, pr in pivot_col.items():
            val = (-dense[pr][fc]) % p
            if val:
                vec[c] = val
        basis.append(vec)
    return basis


def pi_kernel_nilpotency(n: int, q: int, t: int, scratch: str, char: int, log=print):
    """Observe proposed quotient-map nullity and products at one specialization.

    The square-zero test is conditional on the supplied normal forms being
    correct for complete presentations; it supplies no radical verdict here.
    """
    validate_field_characteristic(char)
    if char == 0:
        raise ValueError("a prime characteristic is required")
    p = char
    Aj = PointAlgebra(Spec("square-zero", n), q, t, scratch, char)
    Al = PointAlgebra(Spec("linear", n), q, t, scratch, char)
    dj, dl = Aj.dim, Al.dim
    # Compute nullity of the proposed reduction map; the expected dimension difference is a consistency check.
    M = {}
    for j, b in enumerate(Aj.standard):
        for ii, c in Al.nf(b).items():
            M[(ii, j)] = c % p
    ker = _nullspace_cols_Fp(M, dj, p)
    log(f"[struct n={n}] observed standard-word counts: square-zero={dj}, linear={dl}; "
        f"reduction-map nullity={len(ker)} (expected count difference {dj-dl})")
    assert len(ker) == dj - dl, "observed reduction-map nullity differs from the word-count difference"

    def e_poly(i):  # e_i as free-algebra poly (dict word->coeff mod p)
        x, y = f"x{i}", f"y{i}"
        X2 = {(y,): q % p, (): (1 - q) % p, (x,): (-1) % p}
        pre = {(x,): 1, (): (-t) % p}  # (x_i − t)
        out = {}
        for wp, cp in pre.items():
            for w2, c2 in X2.items():
                w = wp + w2; out[w] = (out.get(w, 0) + cp * c2) % p
                if not out[w]: del out[w]
        return out

    es = {i: e_poly(i) for i in range(1, n)}
    # Consistency check: each displayed generator reduces to zero in the linear model.
    for i, e in es.items():
        vlin = {}
        for w, c in e.items():
            for k, ck in Al.nf(w).items(): vlin[k] = (vlin.get(k, 0) + c * ck) % p
        assert not any(v % p for v in vlin.values()), f"observed reduction of e_{i} is nonzero in the linear model"

    # Observe all displayed e_i · m · e_j reductions; an ideal-square conclusion requires complete bases.
    def nf_poly_word(poly, mid_word, poly2):
        out = {}
        for w1, c1 in poly.items():
            for w2, c2 in poly2.items():
                for k, ck in Aj.nf(w1 + mid_word + w2).items():
                    out[k] = (out.get(k, 0) + c1 * c2 * ck) % p
                    if not out[k]: del out[k]
        return out

    nonzero_pairs = 0
    checked = 0
    for i in es:
        for j in es:
            for m in Aj.standard:
                r = nf_poly_word(es[i], m, es[j]); checked += 1
                if any(v % p for v in r.values()):
                    nonzero_pairs += 1
    nilpotent2 = (nonzero_pairs == 0)
    log(f"[struct n={n}] checked {checked} displayed products e_i·m·e_j; "
        f"nonzero reductions: {nonzero_pairs}; observed reductions all vanish: {nilpotent2}")
    return {"n": n, "q": q, "t": t, "char": char,
            "status": "observed", "basis_complete": False,
            "standard_word_count_square_zero": dj, "standard_word_count_linear": dl,
            "reduction_map_nullity": len(ker),
            "observed_square_zero_products": nilpotent2,
            "products_checked": checked}
