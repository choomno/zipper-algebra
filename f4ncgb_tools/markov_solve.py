#!/usr/bin/env python3
"""
Solve the Markov-trace linear system on the tower Z_3 -> Z_4 for a given
variant (framing-power a or Jordan), over F_p.

The setup. A Markov trace 5-tuple (tau_3, tau_4, z, mu_+, mu_-) consists of
  - tau_n in Tr(Z_n^var) = null space of [Z_n, Z_n] as functionals,
  - scalars z, mu_+, mu_- in k.
The constraints are:
  (M2) compatibility: tau_4(m) = z*tau_3(m)  for each basis element m of Z_3
       (regarded inside Z_4 via the tower inclusion -- verified to be the
        literal inclusion of standard monomials);
  (M3+) Markov-sigma_3:    tau_4(m*sigma_3)   = mu_+*tau_3(m) for each m basis of Z_3;
  (M3-) Markov-sigma_3^-1: tau_4(m*sigma_3^-1) = mu_-*tau_3(m) for each m basis of Z_3.

In coordinates:
  Let V_3 = trace space of Z_3 (5-dim), V_4 = trace space of Z_4 (11-dim).
  Choose bases of V_n by null space of [Z_n, Z_n] matrix.
  Each tau_n is a vector in V_n. All the constraints above are linear in
  (tau_3, tau_4, z, mu_+, mu_-), so we get a single big matrix and we want
  its null space dimension.

  The unknowns: 5 + 11 + 3 = 19.
  The equations:
    (M2)  dim(Z_3) = 24 equations (one per basis monomial of Z_3).
    (M3+) dim(Z_3) = 24 equations.
    (M3-) dim(Z_3) = 24 equations.
  Total 72 equations.

  Many equations are linearly dependent (e.g. M2 contributes only dim(V_3)=5
  independent ones modulo the fact that tau_4 is a trace, etc.).

  The dim of the solution space is what we want.
"""
from __future__ import annotations
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from semisimple_check import (
    parse_basis_full, collect_variables, build_aho_corasick,
    enumerate_standard_monomials, find_leading_match,
)
from check_x3_nilpotent import reduce_to_basis


def build_commutator_matrix(N, standard, mon_idx, lm_to_tail,
                            ac_d, ac_m, ac_p, p):
    """[A,A] matrix: rows = m_i*m_j - m_j*m_i reduced, sparse dicts."""
    rows = []
    for i in range(N):
        for j in range(i + 1, N):
            v1 = reduce_to_basis([(1, standard[i] + standard[j])],
                                 mon_idx, lm_to_tail, ac_d, ac_m, ac_p, p)
            v2 = reduce_to_basis([(1, standard[j] + standard[i])],
                                 mon_idx, lm_to_tail, ac_d, ac_m, ac_p, p)
            diff = dict(v1)
            for k, c in v2.items():
                new = (diff.get(k, 0) - c) % p
                if new == 0:
                    if k in diff:
                        del diff[k]
                else:
                    diff[k] = new
            if diff:
                rows.append(diff)
    return rows


def reduced_echelon_kernel(rows, N, p):
    """Compute a kernel basis (over F_p) of a matrix given as sparse rows.

    Returns: list of N-length vectors (lists of ints mod p), one per kernel direction.
    Echelon reduce the rows in place; pivot columns track which columns are
    "constrained". Non-pivot columns parametrize the kernel.
    """
    # Convert sparse rows to dense and reduce
    M = []
    for r in rows:
        dense = [0] * N
        for k, c in r.items():
            dense[k] = c
        M.append(dense)
    rank = 0
    pivot_col = []  # pivot column for each pivot row
    col_of_row = [-1] * len(M)
    for col in range(N):
        pivot = -1
        for r in range(rank, len(M)):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]
        inv = pow(M[rank][col], p - 2, p)
        for r in range(len(M)):
            if r != rank and M[r][col] != 0:
                fac = (M[r][col] * inv) % p
                for c2 in range(col, N):
                    M[r][c2] = (M[r][c2] - fac * M[rank][c2]) % p
        # normalize pivot to 1
        for c2 in range(N):
            M[rank][c2] = (M[rank][c2] * inv) % p
        pivot_col.append(col)
        col_of_row[rank] = col
        rank += 1
    # Free columns
    pivot_set = set(pivot_col)
    free_cols = [c for c in range(N) if c not in pivot_set]
    # Each free col gives one kernel direction
    kernel = []
    for fc in free_cols:
        vec = [0] * N
        vec[fc] = 1
        # Solve for pivot positions: M_r * vec = 0 means
        #   M_r[pivot_col[r]] * vec[pivot_col[r]] + M_r[fc] * vec[fc] = 0
        # (since other positions are 0 in normalized rows for r-th pivot
        #  and fc isn't a pivot)
        # Actually after row reduction the row at index r has 1 at pivot_col[r]
        # and 0 at all other pivot cols, and various values at free cols.
        # So M[r][pivot_col[r]] * vec[pivot_col[r]] + sum_{fc in free} M[r][fc] * vec[fc] = 0.
        # We want vec[fc] = 1, vec[other_free] = 0, so M[r][pivot_col[r]] = 1
        # gives vec[pivot_col[r]] = - M[r][fc] mod p.
        for r in range(rank):
            vec[pivot_col[r]] = (-M[r][fc]) % p
        kernel.append(vec)
    return kernel, free_cols


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--basis3", required=True, help=".basis at n=3")
    ap.add_argument("--basis4", required=True, help=".basis at n=4")
    ap.add_argument("--variant", required=True, choices=['framing', 'jordan'])
    args = ap.parse_args()
    p = 2147483647

    # === Parse n=3 ===
    polys3 = parse_basis_full(args.basis3)
    vars3 = collect_variables(polys3)
    lms3 = [lm for lm, _ in polys3]
    ac3_d, ac3_m, ac3_p, _ = build_aho_corasick(lms3, vars3)
    std3 = enumerate_standard_monomials(ac3_d, ac3_m, None, vars3, max_total=200)
    N3 = len(std3)
    print(f"[Z_3] N = {N3}", flush=True)
    mi3 = {m: i for i, m in enumerate(std3)}
    lmt3 = {}
    for lm, terms in polys3:
        lead = terms[0][0]
        inv = pow(lead, p - 2, p)
        lmt3[lm] = [((-c * inv) % p, mon) for c, mon in terms[1:]]

    # === Parse n=4 ===
    polys4 = parse_basis_full(args.basis4)
    vars4 = collect_variables(polys4)
    lms4 = [lm for lm, _ in polys4]
    ac4_d, ac4_m, ac4_p, _ = build_aho_corasick(lms4, vars4)
    std4 = enumerate_standard_monomials(ac4_d, ac4_m, None, vars4, max_total=400)
    N4 = len(std4)
    print(f"[Z_4] N = {N4}", flush=True)
    mi4 = {m: i for i, m in enumerate(std4)}
    lmt4 = {}
    for lm, terms in polys4:
        lead = terms[0][0]
        inv = pow(lead, p - 2, p)
        lmt4[lm] = [((-c * inv) % p, mon) for c, mon in terms[1:]]

    # === Build commutator matrices and trace spaces ===
    print(f"[Z_3] building commutator matrix...", flush=True)
    comm3 = build_commutator_matrix(N3, std3, mi3, lmt3, ac3_d, ac3_m, ac3_p, p)
    print(f"[Z_3] {len(comm3)} commutator rows", flush=True)
    print(f"[Z_3] computing trace basis (null of commutator transpose)...", flush=True)
    trace3, free3 = reduced_echelon_kernel(comm3, N3, p)
    d3 = len(trace3)
    print(f"[Z_3] trace space dim = {d3}", flush=True)

    print(f"[Z_4] building commutator matrix...", flush=True)
    comm4 = build_commutator_matrix(N4, std4, mi4, lmt4, ac4_d, ac4_m, ac4_p, p)
    print(f"[Z_4] {len(comm4)} commutator rows", flush=True)
    print(f"[Z_4] computing trace basis (null of commutator transpose)...", flush=True)
    trace4, free4 = reduced_echelon_kernel(comm4, N4, p)
    d4 = len(trace4)
    print(f"[Z_4] trace space dim = {d4}", flush=True)

    # === Tower inclusion: Z_3 monomials -> Z_4 monomials ===
    # We've verified that std3 ⊆ std4. So for each m in std3, its index in std4 is mi4[m].
    incl = {}
    for m in std3:
        if m in mi4:
            incl[std3.index(m)] = mi4[m]
        else:
            print(f"WARNING: {m} in std3 not in std4. Tower inclusion broken?")
    print(f"[tower] {len(incl)} of {N3} std3 monomials embed into std4", flush=True)

    # === Set up the Markov linear system ===
    # Unknowns: tau_3 (d3 components, coefficients in some basis of V_3),
    #           tau_4 (d4 components), z (1), mu_+ (1), mu_- (1).
    # Total = d3 + d4 + 3.
    nvars = d3 + d4 + 3
    print(f"[Markov] total unknowns: {nvars} = {d3} (tau_3) + {d4} (tau_4) + 3 (z, mu_+, mu_-)")

    # tau_3 acts on Z_3 as: tau_3(m_i) = sum_l (tau_3 components)[l] * trace3[l][i]
    # Similarly tau_4(m_j) = sum_l (tau_4 components)[l] * trace4[l][j]

    # Now sigma_3 corresponds to the variable name 'x3' in the Z_4 alphabet
    # (the n-1=3-th generator).
    sigma_idx = 'x3'
    sigma_inv = 'y3'
    if sigma_idx not in vars4:
        print(f"ERROR: '{sigma_idx}' not in Z_4 variables {vars4}")
        return 1

    # For each basis element a of Z_3 (indexed i in std3), and for each
    # constraint type, build a row of the constraint matrix.
    # Constraint matrix columns indexed by [tau_3 (d3), tau_4 (d4), z, mu+, mu-].
    rows = []

    # (M2): tau_4(m_i^[Z_4]) - z * tau_3(m_i^[Z_3]) = 0 for each i in std3
    # tau_4(m_i^[Z_4]) = sum_l (tau_4)[l] * trace4[l][incl[i]]
    # tau_3(m_i^[Z_3]) = sum_l (tau_3)[l] * trace3[l][i]
    # But z * tau_3(m_i) is BILINEAR (z times tau_3), so this is nonlinear unless we
    # parametrize differently. Standard trick: introduce z as another variable but
    # constraint becomes (tau_4)[l] * trace4 - z * (tau_3)[l] * trace3, which is
    # quadratic in unknowns.
    #
    # We solve this by FIXING z = 1 (loss of generality: the trace can be rescaled).
    # Similarly fix mu+ and mu- as one branch of the curve (or treat them as scalars
    # we solve for AFTER imposing constraints).
    #
    # A cleaner approach: separate the linear part from the scaling.
    #
    # Approach: fix tau_3 first. Then (M2) becomes: tau_4 restricted to the image of
    # std3 equals z * tau_3, which is linear in tau_4 for fixed tau_3 and z.
    # (M3) similarly linear in tau_4 with fixed tau_3, mu+, mu-.
    #
    # We want the dim of pairs (tau_3, tau_4, scalars) up to overall scaling. To get
    # a TRUE linear system, we can ENUMERATE over tau_3 in a basis of V_3 (d3 cases),
    # each time solving for (tau_4, z, mu+, mu-) given fixed tau_3. This gives the
    # solution dimension for each tau_3; total Markov solutions = direct sum.
    #
    # BUT actually the right thing is: treat z, mu+, mu- as new unknowns and the
    # products z*tau_3, mu+*tau_3, mu-*tau_3 introduce quadratic terms. We
    # LINEARIZE by absorbing: let u_l := z * tau_3[l] (new vector u in V_3, also
    # d3 unknowns). Then (M2) becomes a LINEAR constraint tau_4 = u (on the
    # restricted indices). Similarly let v_l := mu+ * tau_3[l] and w_l := mu- * tau_3[l].
    # But u is constrained to be a SCALAR MULTIPLE of tau_3, similarly v, w. This
    # itself is a non-linear constraint.
    #
    # FINAL APPROACH: treat (tau_3, tau_4) as variables and impose:
    #   (M2): tau_4|_{Z_3} - z * tau_3 = 0
    # by fixing z = 1 first. Then the system is linear and we count solutions.
    # The freedom in z is one extra dimension on top of the linear solution.
    # Similar for mu+/mu-.
    #
    # For now: assume z is given (treat it implicitly).
    # The Markov system has, with z, mu+, mu- treated as parameters:
    #   24 equations from (M2) for each i in std3
    #   24 from (M3+)
    #   24 from (M3-)
    # All linear in (tau_3, tau_4) given (z, mu+, mu-).

    # CONCRETE SETUP: Set up the matrix with unknowns
    # x = [t_3[0], ..., t_3[d3-1], t_4[0], ..., t_4[d4-1], z, mu_+, mu_-]
    # And constraints are quadratic; we LINEARIZE the products as new variables:
    # Let u_l = z * t_3[l] for l in 0..d3-1: d3 new vars.
    # Let v_l = mu_+ * t_3[l]: d3 new vars.
    # Let w_l = mu_- * t_3[l]: d3 new vars.
    # Now (M2): tau_4(m_i) = u . trace3[*, i]  (linear in tau_4 and u)
    # (M3+): tau_4(m_i*sigma_3) = v . trace3[*, i]
    # (M3-): tau_4(m_i*sigma_3^{-1}) = w . trace3[*, i]
    # The "scaling consistency" u_l/t_3[l] all equal (= z), etc., are nonlinear
    # and we ignore them for the COUNT (gives upper bound on Markov dim).
    # An upper bound on Markov dim = solution_dim(linear system) (treating
    # u, v, w as independent of tau_3).
    # The actual Markov dim is at most this minus the nonlinear constraints,
    # but for the SOLUTION DIMENSION it's a useful first step.

    # New unknowns: (tau_3, tau_4, u, v, w) with dims (d3, d4, d3, d3, d3).
    nvars = d3 + d4 + 3 * d3  # = 5 + 11 + 15 = 31 for n=3->4

    print(f"[Markov-linearized] unknowns: {nvars}")
    print(f"  - tau_3 (V_3 components): {d3}")
    print(f"  - tau_4 (V_4 components): {d4}")
    print(f"  - u = z*tau_3:           {d3}")
    print(f"  - v = mu_+*tau_3:        {d3}")
    print(f"  - w = mu_-*tau_3:        {d3}")

    # Slices:
    s_tau3 = 0
    s_tau4 = d3
    s_u = d3 + d4
    s_v = d3 + d4 + d3
    s_w = d3 + d4 + 2 * d3
    assert s_w + d3 == nvars

    # Helper: compute tau_n vector dotted with a sparse v in std_n basis.
    # tau_n(vec) = sum_i vec[i] * tau_n(m_i)
    # tau_n(m_i) = sum_l tau_n[l] * trace_basis[l][i]
    # so tau_n(vec) = sum_{l, i} tau_n[l] * trace_basis[l][i] * vec[i]
    #              = sum_l tau_n[l] * (trace_basis[l] . vec)

    # ===== Build rows =====
    constraints = []

    def fresh_row():
        return [0] * nvars

    # (M2): for each i in std3, tau_4(m_i^{Z_4}) - sum_l u_l * trace3[l][i] = 0
    for i in range(N3):
        if i not in incl:
            continue
        j = incl[i]
        row = fresh_row()
        for l in range(d4):
            row[s_tau4 + l] = trace4[l][j] % p
        for l in range(d3):
            row[s_u + l] = (-trace3[l][i]) % p
        constraints.append(row)

    # (M3+): for each i in std3, tau_4(m_i * sigma_3) - sum_l v_l * trace3[l][i] = 0
    # Need to reduce m_i * sigma_3 in Z_4 to a sparse vector and dot with trace4
    print(f"[Markov] computing m_i * sigma_3 reductions in Z_4 for {N3} basis elements...", flush=True)
    sigma_reds = []
    for i, m in enumerate(std3):
        red = reduce_to_basis([(1, m + (sigma_idx,))],
                              mi4, lmt4, ac4_d, ac4_m, ac4_p, p)
        sigma_reds.append(red)

    for i in range(N3):
        red = sigma_reds[i]
        row = fresh_row()
        # tau_4(m_i * sigma_3) = sum_j red[j] * tau_4(m_j)
        #                     = sum_l tau_4[l] * (sum_j red[j] * trace4[l][j])
        for l in range(d4):
            dot = 0
            for j, c in red.items():
                dot = (dot + c * trace4[l][j]) % p
            row[s_tau4 + l] = dot
        for l in range(d3):
            row[s_v + l] = (-trace3[l][i]) % p
        constraints.append(row)

    # (M3-): for each i in std3, tau_4(m_i * sigma_3^{-1}) - sum_l w_l * trace3[l][i] = 0
    print(f"[Markov] computing m_i * sigma_3^{{-1}} reductions in Z_4...", flush=True)
    sigma_inv_reds = []
    for i, m in enumerate(std3):
        red = reduce_to_basis([(1, m + (sigma_inv,))],
                              mi4, lmt4, ac4_d, ac4_m, ac4_p, p)
        sigma_inv_reds.append(red)

    for i in range(N3):
        red = sigma_inv_reds[i]
        row = fresh_row()
        for l in range(d4):
            dot = 0
            for j, c in red.items():
                dot = (dot + c * trace4[l][j]) % p
            row[s_tau4 + l] = dot
        for l in range(d3):
            row[s_w + l] = (-trace3[l][i]) % p
        constraints.append(row)

    print(f"[Markov] {len(constraints)} total constraint rows on {nvars} unknowns")

    # === Rank of constraint matrix ===
    M = [r[:] for r in constraints]
    rank = 0
    for col in range(nvars):
        pivot = -1
        for r in range(rank, len(M)):
            if M[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]
        inv = pow(M[rank][col], p - 2, p)
        for r in range(len(M)):
            if r != rank and M[r][col] != 0:
                fac = (M[r][col] * inv) % p
                for c2 in range(col, nvars):
                    M[r][c2] = (M[r][c2] - fac * M[rank][c2]) % p
        rank += 1
    null_dim = nvars - rank
    print(f"[Markov] constraint matrix: {len(constraints)} x {nvars}")
    print(f"[Markov] rank = {rank}, null space dim = {null_dim}")
    print()
    print(f"[Markov-linearized upper bound] dim solution space = {null_dim}")
    print(f"  Interpretation: with u, v, w treated as independent of tau_3, this")
    print(f"  is the dimension of (tau_3, tau_4, u, v, w) tuples satisfying the")
    print(f"  linear part of (M2)(M3). The true Markov dim is at most this and")
    print(f"  generically equals this minus the constraint that u, v, w are scalar")
    print(f"  multiples of tau_3 (nonlinear). For the upper bound: null_dim.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
