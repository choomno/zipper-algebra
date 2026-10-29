#!/usr/bin/env python3
r"""Linear relaxation of Markov-trace conditions on supplied rational algebras.

For the proposed inclusion iota: Z_3 -> Z_4, the Markov conditions compare
three restrictions of a trace tau_4 with scalar multiples of one trace tau_3:

    tau_4(iota(a))       = z    * tau_3(a)
    tau_4(iota(a) * x_3) = mu_+ * tau_3(a)
    tau_4(iota(a) * y_3) = mu_- * tau_3(a).

The scalars and tau_3 enter bilinearly. This program instead introduces three
independent auxiliary traces u, p, q on Z_3 and solves the resulting linear
system in their trace-space coordinates together with those of tau_4. It
reports this linear system's nullity. It then counts chosen kernel basis
vectors whose three auxiliary blocks have rank at most one.

That count depends on the chosen kernel basis. The code does not compute the
rank-one locus of all kernel combinations, a Kauffman trace pullback, or the
dimension of the nonlinear Markov-trace solution variety.

The supplied reduction models must come from certified complete presentations
to interpret these calculations as statements about the presented quotients.
"""
from __future__ import annotations
import argparse
from fractions import Fraction


from zipper_algebra import core, Presentation


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
    pivot_row_of = {c: i for i, c in enumerate(pivot_cols)}
    for fc in free_cols:
        v = [Fraction(0)] * N
        v[fc] = Fraction(1)
        for pc in pivot_cols:
            v[pc] = -R[pivot_row_of[pc]][fc]
        basis.append(v)
    return basis


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


def build_markov_system(alg3, alg4, t3_basis, t4_basis,
                        sigma_big_x="x3", sigma_big_y="y3"):
    """Build the linear relaxation in coordinates (c4, u, p, q).

    c4 parameterizes a trace on the larger reduction model. Each auxiliary
    block independently parameterizes a trace on the smaller model. The rows
    identify these blocks with the three restrictions of c4. No proportionality
    condition is imposed here; nullity is not a nonlinear trace-family dimension.
    """
    d3 = len(t3_basis)
    d4 = len(t4_basis)
    N3 = alg3["N"]
    standard3 = alg3["standard"]

    # Precompute: t4_basis[i](iota(a_k)) for all i, k.
    # iota(a_k) is a_k viewed as a word in the big algebra's variables (the
    # variables of Z_3 are a subset). We just reduce a_k in Z_4 (mon_idx_big,
    # lm_to_tail_big, ac_*_big).
    iota_red = []
    for a in standard3:
        red = core.reduce_to_basis_Q(
            [(Fraction(1), a)],
            alg4["mon_idx"], alg4["lm_to_tail"], alg4["ac_d"], alg4["ac_m"], alg4["ac_p"],
        )
        iota_red.append(red)
    iota_x = []
    for a in standard3:
        word = a + (sigma_big_x,)
        red = core.reduce_to_basis_Q(
            [(Fraction(1), word)],
            alg4["mon_idx"], alg4["lm_to_tail"], alg4["ac_d"], alg4["ac_m"], alg4["ac_p"],
        )
        iota_x.append(red)
    iota_y = []
    for a in standard3:
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


def analyse_variant(label, basis3_path, basis4_path):
    print(f"\n{'='*60}\n=== {label} ===\n{'='*60}", flush=True)
    print(f"[load] building Z_3 from {basis3_path}", flush=True)
    alg3 = build_algebra(basis3_path, 3)
    print(f"[load] N3 = {alg3['N']}, vars = {alg3['variables']}", flush=True)
    print(f"[load] building Z_4 from {basis4_path}", flush=True)
    alg4 = build_algebra(basis4_path, 4)
    print(f"[load] N4 = {alg4['N']}, vars = {alg4['variables']}", flush=True)

    # Trace spaces
    print('[T_3] computing commutator matrix...', flush=True)
    comm3 = commutator_matrix(alg3)
    print(f"[T_3] commutators: {len(comm3)}; computing kernel...", flush=True)
    t3_basis = trace_space_basis(alg3, comm3)
    d3 = len(t3_basis)
    print(f"[T_3] dim T_3 = {d3}", flush=True)

    print('[T_4] computing commutator matrix (slow)...', flush=True)
    comm4 = commutator_matrix(alg4)
    print(f"[T_4] commutators: {len(comm4)}; computing kernel...", flush=True)
    t4_basis = trace_space_basis(alg4, comm4)
    d4 = len(t4_basis)
    print(f"[T_4] dim T_4 = {d4}", flush=True)

    # Build Markov constraints
    print('[markov] building linearised Markov system...', flush=True)
    rows, d3_, d4_ = build_markov_system(alg3, alg4, t3_basis, t4_basis)
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
    print(f"[markov] kernel basis vectors with auxiliary rank at most one = "
          f"{rank1_solutions}/{len(kernel)}", flush=True)

    # This inspects the selected basis directions, not all kernel combinations.
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


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--framing-small", help="explicit basis file for this algebra")
    ap.add_argument("--framing-big", help="explicit basis file for this algebra")
    ap.add_argument("--square-zero-small", help="explicit basis file for this algebra")
    ap.add_argument("--square-zero-big", help="explicit basis file for this algebra")
    ap.add_argument("--only", choices=["framing", "square-zero", "both"], default="both")
    args = ap.parse_args(argv)
    families = ("framing", "square-zero") if args.only == "both" else (args.only,)
    missing = [f"--{family}-{size}" for family in families for size in ("small", "big")
               if getattr(args, f"{family.replace('-', '_')}_{size}") is None]
    if missing:
        ap.error("basis files required for the selected variants: " + ", ".join(missing))

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
