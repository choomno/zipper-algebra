#!/usr/bin/env python3
r"""
Markov-trace linear system search for the framing and square-zero variants of the
zipper algebra towers (Z_3 -> Z_4).

Setup
-----
A Markov trace on the tower {Z_n} is data (tau_3, tau_4, z, mu_+, mu_-) with:

  (M1) tau_n is a linear functional on Z_n vanishing on [Z_n, Z_n].
       So tau_n is a vector of length N_n, in the kernel of the transpose of
       the commutator matrix.

  (M2) Compatibility: for every a in Z_3,
            tau_4(iota(a)) = z * tau_3(a).
       This is a system of N_3 linear equations in
       (tau_3, tau_4) once z is fixed (and a bilinear constraint otherwise).

  (M3) Markov property: for every a in Z_3,
            tau_4(iota(a) * x_3) = mu_+ * tau_3(a),
            tau_4(iota(a) * y_3) = mu_- * tau_3(a).
       Linear in tau_3, tau_4 once mu_+, mu_- are fixed.

We linearise the system by treating (tau_3, tau_4) as unknowns of total length
N_3 + N_4 and (z, mu_+, mu_-) as parameters. We will then sweep over a small
generic parameter set; for the linear-algebra "dimension" of the Markov-trace
space we use the standard trick:

  - Encode (tau_3, tau_4, z*tau_3, mu_+*tau_3, mu_-*tau_3) as ONE unknown
    vector w in R^{N_3 + N_4 + 3*N_3}, plus the trace-vanishing constraints
    on tau_3 and tau_4. We then write tau_3 = v_3 in the trace space T_3,
    tau_4 = v_4 in T_4 and z, mu_+, mu_- as scalars; this gives an unknown
    vector (v_3, v_4, z, mu_+, mu_-) of length
            dim T_3 + dim T_4 + 3
    The constraints (M2), (M3) become BILINEAR (z*v_3 etc.). To linearise we
    introduce w_z := z*v_3, w_p := mu_+*v_3, w_m := mu_-*v_3 as new unknowns
    of dimension dim T_3 each. Then (M2), (M3) become linear:
            tau_4(iota(a)) = w_z(a),
            tau_4(iota(a) * x_3) = w_p(a),
            tau_4(iota(a) * y_3) = w_m(a)
    But these new unknowns are constrained to be SCALAR MULTIPLES of v_3
    (rank-1 constraint), which is a nonlinear condition. The clean trick:

  - Fix tau_3 to be a generic basis element of T_3 (sweep over a basis), then
    z, mu_+, mu_- become honest unknowns (linear): w = (tau_4, z, mu_+, mu_-)
    in R^{N_4 + 3}, with constraints
        tau_4(iota(a)) - z * tau_3(a) = 0
        tau_4(iota(a) x_3) - mu_+ tau_3(a) = 0
        tau_4(iota(a) y_3) - mu_- tau_3(a) = 0
    plus tau_4 ∈ T_4 (commutator-vanishing).
    The space of (tau_4, z, mu_+, mu_-) is linear: a homogeneous system in
    N_4 + 3 unknowns with 3*N_3 + dim [Z_4,Z_4] constraints. We compute its
    nullity. Then we vary tau_3 over a basis of T_3 and sum/cone.

  - Equivalent simpler description: the full Markov-trace space (with both
    tau_3 and tau_4 as unknowns) is a linear subspace of T_3 x T_4 x R^3
    that is closed under (tau_3, tau_4, z, mu_+, mu_-) -> (lambda tau_3,
    lambda tau_4, z, mu_+, mu_-). So we fix tau_3 generically (or take a
    direct sum over a basis of T_3) and report the dimension.

Implementation chooses the second approach: set up the linear system
        Big_constraint @ (tau_4 || z || mu_+ || mu_- || tau_3) = 0
all together; tau_3 ∈ T_3, tau_4 ∈ T_4. The unknown vector has length
    N_3 + N_4 + 3
with constraints
    [trace-vanishing on Z_3]:      from commutator matrix of Z_3   (N_3 cols of tau_3)
    [trace-vanishing on Z_4]:      from commutator matrix of Z_4   (N_4 cols of tau_4)
    [compatibility]:               N_3 equations in (tau_4, tau_3, z)
    [Markov +]:                    N_3 equations in (tau_4, tau_3, mu_+)
    [Markov -]:                    N_3 equations in (tau_4, tau_3, mu_-)

These are BILINEAR in (z, tau_3), etc. So we LINEARISE by introducing the
"effective trace on Z_3 scaled by z (resp. mu_+, mu_-)" as auxiliary trace
functionals; equivalently, we sweep over (z, mu_+, mu_-) being scalars and
do the linear solve for fixed (z, mu_+, mu_-) is also fine if we treat
tau_3 itself as varying.

We use this concrete approach:
  - Treat tau_3 in T_3 and tau_4 in T_4 as the two trace unknowns.
  - For (M2), the constraint is: there exists z such that
        tau_4(iota(a)) - z * tau_3(a) = 0  for all a in Z_3.
    Equivalently, the vector u_2 := (tau_4 . iota) is a SCALAR multiple of
    tau_3 (as functionals on Z_3). I.e., u_2 and tau_3 are linearly dependent.
  - For (M3 +/-), similarly: u_+ := (tau_4 . (iota(.) * x_3)) is a scalar
    multiple of tau_3, and u_- := (tau_4 . (iota(.) * y_3)) is a scalar
    multiple of tau_3.

To compute the Markov-trace space DIMENSION (over Q) cleanly:
  - For each basis vector e_alpha of T_3, view the system
        u_2 = z_alpha * e_alpha, u_+ = p_alpha * e_alpha, u_- = m_alpha * e_alpha
    with (tau_4 in T_4, z_alpha, p_alpha, m_alpha) the unknowns (4*dim T_3 +
    dim T_4 + 3 unknowns total when we sum over alpha, but linearly).
  - The dim of the Markov-trace space equals the dim of solutions
    (tau_3, tau_4, z, mu_+, mu_-) modulo nothing (they form a vector space
    only when we LINEARISE z*tau_3 etc. as separate variables; see below).

We use the FINAL linearised formulation:
  Unknowns: x = (tau_4 || a || b || c) in R^{N_4 + 3*N_3}
    where a := z*tau_3 (length N_3), b := mu_+*tau_3, c := mu_-*tau_3.
  Constraints:
    (C1) tau_4 ∈ T_4: commutator matrix of Z_4 applied to tau_4 = 0
         (these are linear in tau_4).
    (C2) a, b, c lie in T_3 (i.e., they are trace functionals on Z_3).
    (C3) tau_4(iota(e_i)) = a_i             for all basis e_i of Z_3
    (C4) tau_4(iota(e_i) * x_3) = b_i        for all basis e_i of Z_3
    (C5) tau_4(iota(e_i) * y_3) = c_i        for all basis e_i of Z_3

Without rank-1 conditions on (a, b, c), the dimension of the linearised
solution space is an UPPER BOUND on the actual Markov-trace dimension. In
fact, every Markov trace tau_3 with (z, mu_+, mu_-) gives a solution
(tau_4, z*tau_3, mu_+*tau_3, mu_-*tau_3); conversely, if a, b, c are all
nonzero proportional to a common vector tau_3, we get a Markov trace. So
the gap is governed by the rank of (a, b, c) viewed as functionals on Z_3.

For our reporting: we compute dim of the linearised solution space, then
separately VERIFY whether the BMW pullback is a solution and whether
there's a solution outside its line.

Outputs:
  - dim of trace space T_3
  - dim of trace space T_4
  - dim of linearised solution space
  - whether the BMW Kauffman pullback is a solution
  - whether there's a solution beyond the BMW pullback (== dim of linearised > 1)


Experimental interpretation: dimension, minimality and quotient-structure
conclusions require a certified complete basis or a separate exact witness.
A successful solver exit and agreement at several specializations do not supply
that evidence. State the corresponding evidence when interpreting these observations.
"""
from __future__ import annotations
import argparse
import os
import sys
from fractions import Fraction


from zipper_algebra import core, Presentation


def reduce_word(word, mon_idx, lm_to_tail, ac_d, ac_m, ac_p):
    """Reduce a single word (tuple of variable names) to a sparse vector."""
    return core.reduce_to_basis_Q(
        [(Fraction(1), word)], mon_idx, lm_to_tail, ac_d, ac_m, ac_p,
    )


def build_algebra(basis_path, n, max_dim=500):
    """Parse a char-0 .basis file and assemble the AC automaton, standard-
    monomial basis, monomial index, and reduction-substitution tails into a
    single ``alg`` record passed around by the Markov-system builders below.
    """
    polys = core.parse_basis_Q(basis_path)
    variables = Presentation(n).var_names()
    lms = [lm for lm, _ in polys]
    ac_d, ac_m, ac_p, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        ac_d, ac_m, variables, max_total=max_dim,
    )
    if standard is None:
        raise RuntimeError(f"basis {basis_path} dim > {max_dim}")
    return {
        "variables": variables,
        "standard": standard,
        "mon_idx": {m: i for i, m in enumerate(standard)},
        "lm_to_tail": core.lm_to_tail_Q(polys),
        "ac_d": ac_d,
        "ac_m": ac_m,
        "ac_p": ac_p,
        "N": len(standard),
    }


def commutator_matrix(alg):
    """Return list of sparse rows: each row is m_i*m_j - m_j*m_i as a vector
    in the basis (a dict idx -> Fraction)."""
    N = alg["N"]
    standard = alg["standard"]
    rows = []
    for i in range(N):
        for j in range(i + 1, N):
            v1 = core.reduce_to_basis_Q(
                [(Fraction(1), standard[i] + standard[j])],
                alg["mon_idx"], alg["lm_to_tail"], alg["ac_d"], alg["ac_m"], alg["ac_p"],
            )
            v2 = core.reduce_to_basis_Q(
                [(Fraction(1), standard[j] + standard[i])],
                alg["mon_idx"], alg["lm_to_tail"], alg["ac_d"], alg["ac_m"], alg["ac_p"],
            )
            diff = dict(v1)
            for k, c in v2.items():
                new = diff.get(k, Fraction(0)) - c
                if new == 0 and k in diff:
                    del diff[k]
                else:
                    diff[k] = new
            if diff:
                rows.append(diff)
    return rows


def trace_space_basis(alg, comm_rows):
    """Return a list of length-N vectors (lists of Fractions), each a basis
    of the trace space T_n = (A/[A,A])^* viewed as a subspace of A^*.

    Concretely: T_n is the kernel of the transpose of the commutator matrix.
    Each commutator row r gives a linear constraint sum_k r_k tau_k = 0 on
    tau in A^*. So tau ∈ T_n iff M @ tau = 0 where M has rows r.
    """
    N = alg["N"]
    # Build dense matrix
    M = [[Fraction(0)] * N for _ in range(len(comm_rows))]
    for i, r in enumerate(comm_rows):
        for k, c in r.items():
            M[i][k] = c
    # Compute kernel of M (acting on N-dim column vectors)
    # Reduce M to RREF, identify free columns, build basis.
    R = [row[:] for row in M]
    nrows = len(R)
    pivot_cols = []
    r_idx = 0
    for col in range(N):
        pivot = -1
        for r in range(r_idx, nrows):
            if R[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        if pivot != r_idx:
            R[r_idx], R[pivot] = R[pivot], R[r_idx]
        inv = Fraction(1) / R[r_idx][col]
        for c2 in range(col, N):
            R[r_idx][c2] *= inv
        for r in range(nrows):
            if r != r_idx and R[r][col] != 0:
                fac = R[r][col]
                for c2 in range(col, N):
                    R[r][c2] -= fac * R[r_idx][c2]
        pivot_cols.append(col)
        r_idx += 1
    free_cols = [c for c in range(N) if c not in set(pivot_cols)]
    # Each free col gives a basis vector of the kernel
    basis = []
    pivot_set = set(pivot_cols)
    pivot_row_of = {c: i for i, c in enumerate(pivot_cols)}
    for fc in free_cols:
        v = [Fraction(0)] * N
        v[fc] = Fraction(1)
        for pc in pivot_cols:
            v[pc] = -R[pivot_row_of[pc]][fc]
        basis.append(v)
    return basis


def vec_dot(v_funct, v_alg):
    """Functional v_funct ∈ k^N applied to algebra vector v_alg (dict or list)."""
    if isinstance(v_alg, dict):
        s = Fraction(0)
        for k, c in v_alg.items():
            s += v_funct[k] * c
        return s
    else:
        s = Fraction(0)
        for k in range(len(v_alg)):
            if v_alg[k] != 0 and v_funct[k] != 0:
                s += v_funct[k] * v_alg[k]
        return s


def gaussian_rank_nullity(M, ncols):
    """Compute rank of M (list of rows, each list of Fractions of length ncols),
    and return rank and a basis of the right kernel."""
    R = [row[:] for row in M]
    nrows = len(R)
    pivot_cols = []
    r_idx = 0
    for col in range(ncols):
        pivot = -1
        for r in range(r_idx, nrows):
            if R[r][col] != 0:
                pivot = r
                break
        if pivot < 0:
            continue
        if pivot != r_idx:
            R[r_idx], R[pivot] = R[pivot], R[r_idx]
        inv = Fraction(1) / R[r_idx][col]
        for c2 in range(col, ncols):
            R[r_idx][c2] *= inv
        for r in range(nrows):
            if r != r_idx and R[r][col] != 0:
                fac = R[r][col]
                for c2 in range(col, ncols):
                    R[r][c2] -= fac * R[r_idx][c2]
        pivot_cols.append(col)
        r_idx += 1
    free_cols = [c for c in range(ncols) if c not in set(pivot_cols)]
    pivot_row_of = {c: i for i, c in enumerate(pivot_cols)}
    kernel_basis = []
    for fc in free_cols:
        v = [Fraction(0)] * ncols
        v[fc] = Fraction(1)
        for pc in pivot_cols:
            v[pc] = -R[pivot_row_of[pc]][fc]
        kernel_basis.append(v)
    return len(pivot_cols), kernel_basis


def build_markov_system(alg3, alg4, comm_rows3, comm_rows4, t3_basis, t4_basis,
                        sigma_big_x="x3", sigma_big_y="y3"):
    """Build the linearised Markov system. Unknown vector layout:

        x = (v_4_coords || v_3_coords || z || mu_+ || mu_-)
    where v_3 is parameterised by its T_3 coordinates (length dim T_3),
    v_4 by its T_4 coordinates (length dim T_4).

    The constraints are BILINEAR in (z, v_3), (mu_+, v_3), (mu_-, v_3); to
    have a TRULY linear system we substitute:
        a = z * v_3   (length N_3, an element of Z_3^*)
        b = mu_+ * v_3
        c = mu_- * v_3
    But then we should keep track that a, b, c are scalar multiples of the
    SAME v_3. We capture this by also keeping v_3 as an unknown, with the
    rank-1 condition being NONLINEAR; the linearised dim is an upper bound.

    To get the EXACT Markov-trace dim, we sweep over each fixed v_3 (=e_alpha)
    in T_3 and solve linearly for (v_4, z, mu_+, mu_-). This gives the dim
    of Markov traces FIXING v_3 = e_alpha; summing/conjoining the solutions
    over different e_alpha gives the full Markov-trace dim.

    But since "z, mu_+, mu_-" multiply v_3, the right object is:
    we solve for (v_4, w_z, w_p, w_m) where w_z, w_p, w_m are auxiliary
    trace functionals on Z_3, plus the rank-1 condition (w_z, w_p, w_m all
    proportional to one v_3) which is a Segre-variety condition.

    We solve the linearised system and report its dim.
    """
    # We parameterise:
    #   tau_4 in T_4: coord in R^{dim T4}; tau_4 = sum c4_i * t4_basis[i]
    #   v_3 in T_3:    tau_3 = sum c3_i * t3_basis[i]
    # Treat scalars z, mu_+, mu_- as unknowns; system is BILINEAR.
    # For each constraint vs basis-element a_k of Z_3:
    #       tau_4(iota(a_k)) = z * tau_3(a_k)
    # This is bilinear; for linear-system purposes we linearise as above by
    # introducing auxiliary unknowns:
    #       A_k := z * tau_3(a_k)     (scalar)  -- but A_k = sum c3_i e3i(a_k) * z
    # Instead use trick: substitute u = z * c3 (i.e. u_i = z * c3_i), and
    # similarly p_i = mu_+ * c3_i, q_i = mu_- * c3_i. These are honest linear
    # unknowns. Then the constraint becomes:
    #       sum_i c4_i * t4_basis[i](iota(a_k)) - sum_i u_i * t3_basis[i](a_k) = 0
    # for each k. Linear! Likewise for x_3 and y_3.
    #
    # The catch: u, p, q must all be SCALAR multiples of c3. That's the
    # rank-1 condition; we report linearised dim and then check rank-1.
    #
    # Unknown layout: x = (c4_1..c4_{d4}, u_1..u_{d3}, p_1..p_{d3}, q_1..q_{d3})
    # Total length: d4 + 3*d3.
    # (c3 itself is the "rank-1 projection" of (u, p, q); we report it separately.)
    #
    # In addition: NOTE the linearised system says nothing about c3 directly.
    # The full Markov-trace space dim is at most d3 + dim(linearised
    # solutions / ~), but since (u, p, q) = (z*c3, mu_+*c3, mu_-*c3) lies in
    # the rank-1 cone of (3*d3)-space, the actual Markov-trace dim equals the
    # number of independent (c3, z, mu_+, mu_-) tuples such that the rank-1
    # ansatz can be filled in linearly. We will compute upper bound and then
    # extract explicit solutions to check rank-1.

    d3 = len(t3_basis)
    d4 = len(t4_basis)
    N3 = alg3["N"]
    N4 = alg4["N"]
    standard3 = alg3["standard"]

    # Precompute: t4_basis[i](iota(a_k)) for all i, k.
    # iota(a_k) is a_k viewed as a word in the big algebra's variables (the
    # variables of Z_3 are a subset). We just reduce a_k in Z_4 (mon_idx_big,
    # lm_to_tail_big, ac_*_big).
    iota_red = []
    for k, a in enumerate(standard3):
        red = core.reduce_to_basis_Q(
            [(Fraction(1), a)],
            alg4["mon_idx"], alg4["lm_to_tail"], alg4["ac_d"], alg4["ac_m"], alg4["ac_p"],
        )
        iota_red.append(red)
    iota_x = []
    for k, a in enumerate(standard3):
        word = a + (sigma_big_x,)
        red = core.reduce_to_basis_Q(
            [(Fraction(1), word)],
            alg4["mon_idx"], alg4["lm_to_tail"], alg4["ac_d"], alg4["ac_m"], alg4["ac_p"],
        )
        iota_x.append(red)
    iota_y = []
    for k, a in enumerate(standard3):
        word = a + (sigma_big_y,)
        red = core.reduce_to_basis_Q(
            [(Fraction(1), word)],
            alg4["mon_idx"], alg4["lm_to_tail"], alg4["ac_d"], alg4["ac_m"], alg4["ac_p"],
        )
        iota_y.append(red)

    # Helper: t4_basis[i] is a vector in R^{N4}, viewing it as a functional on Z_4.
    # t4_basis[i](v) = sum_j t4_basis[i][j] * v[j].
    # t3_basis[i](a) = t3_basis[i][a_index] when a is a basis element.

    # Build constraints. Each constraint is a row of length d4 + 3*d3.
    rows = []  # list of lists of Fractions

    for sym_idx, iota_list in enumerate([iota_red, iota_x, iota_y]):
        # iota_list[k] is a dict idx_in_Z4 -> coeff, representing iota(a_k) [or *x_3, *y_3].
        for k in range(N3):
            row = [Fraction(0)] * (d4 + 3 * d3)
            # tau_4(iota(a_k)) - aux_sym(a_k) = 0
            # aux_sym is u, p, or q depending on sym_idx
            # tau_4 = sum c4_i * t4_basis[i]; tau_4(iota(a_k)) = sum_i c4_i * (t4_basis[i] . iota_list[k])
            for i in range(d4):
                val = Fraction(0)
                for j, cj in iota_list[k].items():
                    val += t4_basis[i][j] * cj
                if val != 0:
                    row[i] = val
            # - aux(a_k) where aux = sum u_i * t3_basis[i]
            # aux(a_k) = t3_basis_i[k] when a_k is the k-th basis element.
            offset = d4 + sym_idx * d3
            for i in range(d3):
                val = t3_basis[i][k]
                if val != 0:
                    row[offset + i] = -val
            # Append
            if any(x != 0 for x in row):
                rows.append(row)

    return rows, d3, d4


def kauffman_pullback_check(alg3, alg4, t3_basis, t4_basis,
                            sigma_big_x="x3", sigma_big_y="y3"):
    """Check the BMW Kauffman trace pulled back to Z_n.

    The Kauffman trace on BMW_n is the unique (up to scalar) Markov trace.
    We don't actually need to compute Kauffman in closed form to check our
    system. Instead we use a CONCRETE CRITERION:

    A linearised Markov solution (tau_4, u, p, q) where u, p, q are linearly
    independent vectors in T_3 (i.e., they span a 3-dim subspace) is NOT
    coming from a rank-1 ansatz, hence not from a single Markov trace.

    A linearised Markov solution where u, p, q are all proportional to one
    common v_3 in T_3 (in the rank-1 cone) corresponds to a genuine Markov
    trace (tau_3 = v_3, z = u/v_3, mu_+ = p/v_3, mu_- = q/v_3).

    The BMW Kauffman pullback would give a specific 1-parameter family of
    such genuine traces; if it exists then the rank-1 component of our
    solution space has dim >= 1. We test by:
      (1) computing the linearised solution space K (rank, nullity);
      (2) reporting whether any solution has u, p, q all in the same line.
    """
    return None  # placeholder; we report below


def analyse_variant(label, basis3_path, basis4_path):
    print(f"\n{'='*60}\n=== {label} ===\n{'='*60}", flush=True)
    print(f"[load] building Z_3 from {basis3_path}", flush=True)
    alg3 = build_algebra(basis3_path, 3)
    print(f"[load] N3 = {alg3['N']}, vars = {alg3['variables']}", flush=True)
    print(f"[load] building Z_4 from {basis4_path}", flush=True)
    alg4 = build_algebra(basis4_path, 4)
    print(f"[load] N4 = {alg4['N']}, vars = {alg4['variables']}", flush=True)

    # Trace spaces
    print(f"[T_3] computing commutator matrix...", flush=True)
    comm3 = commutator_matrix(alg3)
    print(f"[T_3] commutators: {len(comm3)}; computing kernel...", flush=True)
    t3_basis = trace_space_basis(alg3, comm3)
    d3 = len(t3_basis)
    print(f"[T_3] dim T_3 = {d3}", flush=True)

    print(f"[T_4] computing commutator matrix (slow)...", flush=True)
    comm4 = commutator_matrix(alg4)
    print(f"[T_4] commutators: {len(comm4)}; computing kernel...", flush=True)
    t4_basis = trace_space_basis(alg4, comm4)
    d4 = len(t4_basis)
    print(f"[T_4] dim T_4 = {d4}", flush=True)

    # Build Markov constraints
    print(f"[markov] building linearised Markov system...", flush=True)
    rows, d3_, d4_ = build_markov_system(alg3, alg4, comm3, comm4, t3_basis, t4_basis)
    ncols = d4 + 3 * d3
    print(f"[markov] matrix size: {len(rows)} rows x {ncols} cols", flush=True)
    rank, kernel = gaussian_rank_nullity(rows, ncols)
    print(f"[markov] rank = {rank}, nullity (linearised solution dim) = {len(kernel)}",
          flush=True)

    # Each kernel vector is (c4_1..c4_{d4}, u_1..u_{d3}, p_1..p_{d3}, q_1..q_{d3}).
    # Check rank-1 (Segre) condition: are u, p, q proportional in T_3?
    # For each kernel vector, extract (u, p, q) and check linear dependence.
    rank1_solutions = 0
    for v in kernel:
        u = v[d4:d4 + d3]
        p = v[d4 + d3:d4 + 2 * d3]
        q = v[d4 + 2 * d3:d4 + 3 * d3]
        # All three must be in the same line (over Q).
        # Stack as rows; compute rank.
        M = [u[:], p[:], q[:]]
        r, _ = gaussian_rank_nullity(M, d3)
        if r <= 1:
            rank1_solutions += 1
    print(f"[markov] kernel basis vectors that are rank-1 (genuine Markov-style) = "
          f"{rank1_solutions}/{len(kernel)}", flush=True)

    # GENERAL linearised solution Segre intersection: the variety of rank-1
    # tuples inside the kernel has dim = dim(rank-1 component). For our small
    # data we just report linearised dim, then separately check explicit
    # solutions.
    return {
        "label": label,
        "N3": alg3["N"], "N4": alg4["N"],
        "dim_T3": d3, "dim_T4": d4,
        "markov_rows": len(rows),
        "markov_cols": ncols,
        "markov_rank": rank,
        "linearised_nullity": len(kernel),
        "rank1_basis_count": rank1_solutions,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--framing-small", default="/tmp/bigelow-n3-a2-Q.basis")
    ap.add_argument("--framing-big", default="/tmp/bigelow-n4-a2-Q.basis")
    ap.add_argument("--square-zero-small", default="/tmp/n3-square-zero-Q.basis")
    ap.add_argument("--square-zero-big", default="/tmp/square-zero-n4.basis")
    ap.add_argument("--only", choices=["framing", "square-zero", "both"], default="both")
    args = ap.parse_args()

    results = []
    if args.only in ("framing", "both"):
        results.append(analyse_variant("FRAMING a=2", args.framing_small, args.framing_big))
    if args.only in ("square-zero", "both"):
        results.append(analyse_variant("SQUARE-ZERO", args.square_zero_small, args.square_zero_big))

    print(f"\n{'='*60}\n=== SUMMARY ===\n{'='*60}", flush=True)
    for r in results:
        print(r)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
