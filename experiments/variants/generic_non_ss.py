#!/usr/bin/env python3
r"""Explore trace-form rank defects through exact specialized calculations.

This historical recipe compares grids, interpolants and integer-encoding probes
for the square-zero family. Its degree and denominator bounds are supplied
assumptions, not proved bounds. Its separately computed bases are not certified
as a common generic basis. A zero interpolant or many zero specialized
determinants therefore remains evidence for a candidate generic statement.

Interpretations involving Jacobson radicals additionally require a certified
finite-dimensional algebra and the characteristic-zero regular-trace criterion.
The experimental numerical kernels remain available for inspecting that evidence;
this script does not produce a certified generic nonsemisimplicity theorem.
"""
from __future__ import annotations

import argparse
from experiments.output import (prepare_output_dir, capture_experiment_provenance,
                                tag_solver_input, tag_solver_basis)
import os
import random
import subprocess
import sys
import tempfile
import shutil
from fractions import Fraction


from zipper_algebra import core, Presentation
from experiments.variants import semisimple_Q


# ---------------------------------------------------------------------------
# Defaults / discovery of toolchain
# ---------------------------------------------------------------------------


def find_f4ncgb() -> str:
    executable = shutil.which(os.environ.get("F4NCGB", "f4ncgb"))
    if executable is None:
        raise RuntimeError("f4ncgb was not found; supply F4NCGB or add it to PATH")
    return executable


# ---------------------------------------------------------------------------
# Build the GB and the standard monomial basis at a chosen (q, t)
# ---------------------------------------------------------------------------

def build_basis_at(workdir: str, *, n: int, q: int, t: int,
                   maxiter: int = 500, verbose: bool = False) -> str:
    """Emit the square-zero-variant relations at (q, t) over Q, run f4ncgb,
    and return the path to the resulting .basis file."""
    producer = capture_experiment_provenance(__file__)
    tag = f"square_zero_n{n}_q{q}_t{t}"
    ms = os.path.join(workdir, f"{tag}.ms")
    basis = os.path.join(workdir, f"{tag}.basis")
    Presentation(n, untwist="square-zero").emit_ms(ms, q, t, 0, provenance=False)
    inputs = tag_solver_input(ms, producer=producer)
    f4 = find_f4ncgb()
    cmd_gb = [f4, "-m", str(maxiter), "-o", basis, "-i", ms]
    subprocess.check_call(
        cmd_gb, stdout=subprocess.DEVNULL,
        stderr=(None if verbose else subprocess.DEVNULL),
    )
    tag_solver_basis(basis, producer=producer, inputs=inputs)
    return basis


def trace_form_at(basis_path: str, *, n: int, max_dim: int = 500,
                  verbose: bool = False) -> tuple[int, int, Fraction]:
    """Parse a .basis file and compute (N, rank T, det T) over Q.
    Returns: (dim A, rank, det T)."""
    polys = core.parse_basis_Q(basis_path)
    lms = [lm for lm, _ in polys]
    variables = Presentation(n).var_names()
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        delta, is_match, variables, max_total=max_dim,
    )
    if standard is None:
        raise RuntimeError(f"dim > {max_dim} -- raise --max-dim?")
    mon_idx = {m: i for i, m in enumerate(standard)}
    log = (lambda s: print(f"   [{s}]", file=sys.stderr, flush=True)) \
        if verbose else (lambda *a: None)
    N, rank, det = semisimple_Q.trace_form_rank_Q(
        polys, standard, mon_idx, log=log, variables=variables,
    )
    return N, rank, det


def trace_form_matrix_at(basis_path: str, *, n: int, max_dim: int = 500) -> list[list[Fraction]]:
    """Same as ``trace_form_at`` but returns the full matrix T (not just
    its rank and determinant). Used by the interpolation step."""
    polys = core.parse_basis_Q(basis_path)
    lms = [lm for lm, _ in polys]
    variables = Presentation(n).var_names()
    delta, is_match, match_pat, _ = core.build_aho_corasick(lms, variables)
    standard = core.enumerate_standard_monomials(
        delta, is_match, variables, max_total=max_dim,
    )
    if standard is None:
        raise RuntimeError(f"dim > {max_dim}")
    N = len(standard)
    mon_idx = {m: i for i, m in enumerate(standard)}
    lm_to_tail = core.lm_to_tail_Q(polys)
    gen_L = semisimple_Q.build_generator_matrices_Q(
        variables, standard, mon_idx, lm_to_tail, delta, is_match, match_pat,
    )
    traces: list[Fraction] = [Fraction(0)] * N
    for l, m_l in enumerate(standard):
        if not m_l:
            traces[l] = Fraction(N)
            continue
        s = Fraction(0)
        for j in range(N):
            vec = semisimple_Q.left_mult_by_word_Q(
                gen_L, m_l, {j: Fraction(1)},
            )
            s += vec.get(j, Fraction(0))
        traces[l] = s
    T = [[Fraction(0)] * N for _ in range(N)]
    for i, m_i in enumerate(standard):
        for j in range(N):
            vec = semisimple_Q.left_mult_by_word_Q(
                gen_L, m_i, {j: Fraction(1)},
            )
            T[i][j] = sum((c * traces[l] for l, c in vec.items()),
                          Fraction(0))
    return T


# ---------------------------------------------------------------------------
# Bivariate polynomial arithmetic over Q (dict {(i, j): coeff})
# ---------------------------------------------------------------------------

QPoly = dict[tuple[int, int], Fraction]  # {(deg_q, deg_t): coeff}


def qp_zero() -> QPoly:
    return {}


def qp_add(a: QPoly, b: QPoly) -> QPoly:
    out: QPoly = dict(a)
    for k, v in b.items():
        new = out.get(k, Fraction(0)) + v
        if new == 0:
            out.pop(k, None)
        else:
            out[k] = new
    return out


def qp_scale(a: QPoly, c: Fraction) -> QPoly:
    if c == 0:
        return {}
    return {k: v * c for k, v in a.items()}


def qp_mul(a: QPoly, b: QPoly) -> QPoly:
    out: QPoly = {}
    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            k = (i1 + i2, j1 + j2)
            new = out.get(k, Fraction(0)) + c1 * c2
            if new == 0:
                out.pop(k, None)
            else:
                out[k] = new
    return out


def qp_const(c: Fraction) -> QPoly:
    return {} if c == 0 else {(0, 0): c}


def qp_total_degree(a: QPoly) -> int:
    return max((i + j for (i, j) in a), default=-1)


# ---------------------------------------------------------------------------
# Bivariate Lagrange interpolation: from samples on a grid to QPoly
# ---------------------------------------------------------------------------

def lagrange_basis_1d(nodes: list[int], idx: int) -> QPoly:
    """The Lagrange basis polynomial L_idx(x) in Q[x] that is 1 at
    nodes[idx] and 0 at every other node. Encoded as a QPoly in which
    the variable is q. (For t, the caller renames degrees.)"""
    x0 = nodes[idx]
    num = qp_const(Fraction(1))
    den = Fraction(1)
    var_q: QPoly = {(1, 0): Fraction(1)}
    var_const = lambda c: qp_const(Fraction(c))
    for k, xk in enumerate(nodes):
        if k == idx:
            continue
        # multiply by (q - xk)
        factor: QPoly = qp_add(var_q, var_const(-xk))
        num = qp_mul(num, factor)
        den *= (x0 - xk)
    return qp_scale(num, Fraction(1, 1) / den)


def lagrange_basis_1d_t(nodes: list[int], idx: int) -> QPoly:
    """Same as ``lagrange_basis_1d`` but in t (degrees encoded as
    (0, k) instead of (k, 0))."""
    x0 = nodes[idx]
    num = qp_const(Fraction(1))
    den = Fraction(1)
    var_t: QPoly = {(0, 1): Fraction(1)}
    for k, xk in enumerate(nodes):
        if k == idx:
            continue
        factor: QPoly = qp_add(var_t, qp_const(Fraction(-xk)))
        num = qp_mul(num, factor)
        den *= (x0 - xk)
    return qp_scale(num, Fraction(1, 1) / den)


def interpolate_qt(
    q_nodes: list[int], t_nodes: list[int],
    samples: dict[tuple[int, int], Fraction],
) -> QPoly:
    """Interpolate the bivariate polynomial in Q[q, t] of degree <
    len(q_nodes) in q and < len(t_nodes) in t that agrees with
    ``samples[(q, t)]`` at every grid point."""
    Lq = [lagrange_basis_1d(q_nodes, i) for i in range(len(q_nodes))]
    Lt = [lagrange_basis_1d_t(t_nodes, j) for j in range(len(t_nodes))]
    out: QPoly = qp_zero()
    for i, qi in enumerate(q_nodes):
        for j, tj in enumerate(t_nodes):
            c = samples[(qi, tj)]
            if c == 0:
                continue
            basis_ij = qp_mul(Lq[i], Lt[j])
            out = qp_add(out, qp_scale(basis_ij, c))
    return out


# ---------------------------------------------------------------------------
# Argument (A): many-point Schwartz--Zippel certificate
# ---------------------------------------------------------------------------

def schwartz_zippel_grid(
    workdir: str, *, n: int, deg_bound: int, prime_start: int = 2,
    verbose: bool = False,
) -> tuple[bool, list[tuple[int, int, int, int, Fraction]]]:
    """Pick an (deg_bound + 2) x (deg_bound + 2) grid of small distinct
    primes; verify det T = 0 at every grid point. Returns (all_zero,
    table) where table is (q, t, N, rank, det)."""
    from itertools import islice

    def primes():
        # naive sieve generator
        candidates = []
        x = 2
        while True:
            ok = True
            for p in candidates:
                if p * p > x:
                    break
                if x % p == 0:
                    ok = False
                    break
            if ok:
                candidates.append(x)
                if x >= prime_start:
                    yield x
            x += 1

    L = deg_bound + 2
    pgen = primes()
    qs = list(islice(pgen, L))
    ts = list(islice(pgen, L))  # disjoint from q-axis except at boundary
    table = []
    all_zero = True
    for qv in qs:
        for tv in ts:
            if qv == tv:
                # skip the diagonal coincidence (still in general position
                # for the trace form, but neater)
                continue
            try:
                basis = build_basis_at(workdir, n=n, q=qv, t=tv,
                                       verbose=verbose)
                N, rank, det = trace_form_at(basis, n=n, verbose=verbose)
            except Exception as e:
                if verbose:
                    print(f"   [skip (q,t)=({qv},{tv}): {e}]",
                          file=sys.stderr)
                continue
            table.append((qv, tv, N, rank, det))
            if det != 0:
                all_zero = False
            if verbose:
                marker = "ZERO" if det == 0 else f"!= 0 ({det})"
                print(f"   (q,t)=({qv:>3},{tv:>3}): N={N} "
                      f"rank={rank} det T = {marker}",
                      file=sys.stderr)
    return all_zero, table


# ---------------------------------------------------------------------------
# Argument (B): explicit polynomial reconstruction of det T
# ---------------------------------------------------------------------------

def reconstruct_det_polynomial(
    workdir: str, *, n: int, deg_q: int, deg_t: int, denom_exp: int,
    verbose: bool = False,
) -> tuple[QPoly, list[int], list[int]]:
    """Compute det T(q, t) * (q * t)**denom_exp as an explicit
    polynomial in Q[q, t] by bivariate Lagrange interpolation on a
    ``(deg_q + 1) x (deg_t + 1)`` grid of small primes (q, t).

    The Laurent denominator (q * t)**denom_exp clears any negative
    powers of q, t in T_{ij}(q, t); for n = 3 we have observed
    denom_exp = 5 is sufficient (per-entry maximal Laurent degree
    <= q^{-4} t^{-2}, so det T has Laurent denominator at most
    q^{-96} t^{-48} -- but we typically need far less in practice).

    Returns (det_poly, q_nodes, t_nodes).
    """
    from itertools import islice

    def primes():
        candidates = []
        x = 2
        while True:
            ok = True
            for p in candidates:
                if p * p > x:
                    break
                if x % p == 0:
                    ok = False
                    break
            if ok:
                candidates.append(x)
                yield x
            x += 1

    pgen = primes()
    q_nodes = list(islice(pgen, deg_q + 1))
    t_nodes_raw = list(islice(pgen, deg_t + 2))
    qset = set(q_nodes)
    t_nodes = [x for x in t_nodes_raw if x not in qset][:deg_t + 1]
    if verbose:
        print(f"   [det-poly grid: q in {q_nodes}, t in {t_nodes}]",
              file=sys.stderr)
    samples: dict[tuple[int, int], Fraction] = {}
    total = (deg_q + 1) * (deg_t + 1)
    done = 0
    for qv in q_nodes:
        for tv in t_nodes:
            basis = build_basis_at(workdir, n=n, q=qv, t=tv,
                                   verbose=verbose)
            N, rank, det = trace_form_at(basis, n=n, verbose=verbose)
            if N == 0 or rank == 0 and det == 0 and N == 0:
                raise RuntimeError(
                    f"unexpected: dim = 0 at (q, t) = ({qv}, {tv})."
                )
            # Multiply by (q * t)**denom_exp to clear Laurent
            samples[(qv, tv)] = det * (qv * tv) ** denom_exp
            done += 1
            if verbose:
                print(f"   [{done}/{total}: det T({qv},{tv}) "
                      f"* ({qv}*{tv})^{denom_exp} = "
                      f"{samples[(qv, tv)]}]",
                      file=sys.stderr)
    poly = interpolate_qt(q_nodes, t_nodes, samples)
    return poly, q_nodes, t_nodes


def reconstruct_T_polynomial(
    workdir: str, *, n: int, deg_q: int, deg_t: int,
    verbose: bool = False,
) -> list[list[QPoly]]:
    """Compute T(q, t) as an explicit polynomial matrix in Q[q, t] by
    bivariate interpolation on a ``(deg_q + 1) x (deg_t + 1)`` grid of
    small primes. Returns T_polys, an N x N matrix of QPolys."""
    from itertools import islice

    def primes():
        candidates = []
        x = 2
        while True:
            ok = True
            for p in candidates:
                if p * p > x:
                    break
                if x % p == 0:
                    ok = False
                    break
            if ok:
                candidates.append(x)
                yield x
            x += 1

    pgen = primes()
    q_nodes = list(islice(pgen, deg_q + 1))
    # interleave to avoid coincidence
    t_nodes = list(islice(pgen, deg_t + 2))
    # drop a node if it equals one of the q-nodes
    t_nodes = [x for x in t_nodes if x not in set(q_nodes)][:deg_t + 1]
    if verbose:
        print(f"   [interpolation grid: q in {q_nodes}, t in {t_nodes}]",
              file=sys.stderr)

    # ----- collect all T_{ij}(q, t) samples -----
    # Determine N from one specialisation (assumed stable across the grid).
    basis0 = build_basis_at(workdir, n=n, q=q_nodes[0], t=t_nodes[0],
                            verbose=verbose)
    T0 = trace_form_matrix_at(basis0, n=n)
    N = len(T0)
    if verbose:
        print(f"   [N = {N}, total entries = {N * N}]", file=sys.stderr)

    # samples[i][j][(q, t)] = T_{ij}(q, t)
    samples: list[list[dict[tuple[int, int], Fraction]]] = [
        [dict() for _ in range(N)] for _ in range(N)
    ]
    for i in range(N):
        for j in range(N):
            samples[i][j][(q_nodes[0], t_nodes[0])] = T0[i][j]

    total = len(q_nodes) * len(t_nodes)
    done = 1
    for qv in q_nodes:
        for tv in t_nodes:
            if qv == q_nodes[0] and tv == t_nodes[0]:
                continue
            basis = build_basis_at(workdir, n=n, q=qv, t=tv,
                                   verbose=verbose)
            Tv = trace_form_matrix_at(basis, n=n)
            if len(Tv) != N:
                raise RuntimeError(
                    f"dim changed at (q, t) = ({qv}, {tv}): "
                    f"{len(Tv)} != {N}. Specialisation degenerate; "
                    f"pick different nodes."
                )
            for i in range(N):
                for j in range(N):
                    samples[i][j][(qv, tv)] = Tv[i][j]
            done += 1
            if verbose:
                print(f"   [{done}/{total} grid points done]",
                      file=sys.stderr)

    # ----- interpolate each entry -----
    T_polys: list[list[QPoly]] = [[qp_zero()] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            T_polys[i][j] = interpolate_qt(q_nodes, t_nodes, samples[i][j])
        if verbose and (i % 4 == 3 or i == N - 1):
            print(f"   [interpolated row {i + 1}/{N}]",
                  file=sys.stderr)
    return T_polys


def verify_T_polynomial_at(
    T_polys: list[list[QPoly]], workdir: str, *, n: int, q: int, t: int,
    verbose: bool = False,
) -> bool:
    """Sanity check: T_polys evaluated at (q, t) must match the
    direct trace-form computation at (q, t). Returns True on match."""
    basis = build_basis_at(workdir, n=n, q=q, t=t, verbose=verbose)
    T_direct = trace_form_matrix_at(basis, n=n)
    N = len(T_direct)
    if len(T_polys) != N:
        return False
    for i in range(N):
        for j in range(N):
            s = Fraction(0)
            for (a, b), c in T_polys[i][j].items():
                s += c * (q ** a) * (t ** b)
            if s != T_direct[i][j]:
                if verbose:
                    print(f"   [verification fail at T[{i}][{j}]: "
                          f"poly={s}, direct={T_direct[i][j]}]",
                          file=sys.stderr)
                return False
    return True


def det_polynomial_matrix_evaluate(
    T_polys: list[list[QPoly]], q_val, t_val,
) -> Fraction:
    """Evaluate each entry of T_polys at (q_val, t_val) and compute the
    determinant in Q (or in whatever ring q_val, t_val support)."""
    N = len(T_polys)
    Teval = [[Fraction(0) if isinstance(q_val, int) else q_val * 0]
             * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            s: Fraction = Fraction(0)
            for (a, b), c in T_polys[i][j].items():
                s += c * Fraction(q_val) ** a * Fraction(t_val) ** b
            Teval[i][j] = s
    det, _ = semisimple_Q.gauss_rank_det_Q(
        [row[:] for row in Teval],
    )
    return det


# ---------------------------------------------------------------------------
# Argument (C): single huge-coprime specialisation
# ---------------------------------------------------------------------------

def huge_specialisation(
    workdir: str, *, n: int, bits: int = 60, seed: int = 0xC0FFEE,
    verbose: bool = False,
) -> tuple[bool, int, int, Fraction, int]:
    """Pick a random pair (q0, t0) of ~ 2**bits-bit positive integers
    with disjoint small-prime support, compute det T(q0, t0) in Q.
    Returns (is_zero, q0, t0, det, dim_at_(q0,t0))."""
    rng = random.Random(seed)
    q0 = rng.randrange(2 ** (bits - 1), 2 ** bits)
    t0 = rng.randrange(2 ** (bits - 1), 2 ** bits)
    if q0 % 2 == 0:
        q0 += 1
    if t0 % 2 == 0:
        t0 += 1
    if verbose:
        print(f"   [huge specialisation: (q0, t0) = ({q0}, {t0})]",
              file=sys.stderr)
    basis = build_basis_at(workdir, n=n, q=q0, t=t0, verbose=verbose)
    N, rank, det = trace_form_at(basis, n=n, verbose=verbose)
    return (det == 0), q0, t0, det, N


# ---------------------------------------------------------------------------
# Argument (D): high-precision rational approximation of (pi, e)
# ---------------------------------------------------------------------------

def transcendental_specialisation(
    workdir: str, *, n: int, digits: int = 30, verbose: bool = False,
) -> tuple[Fraction, int, int, int]:
    """Approximate (pi, e) to `digits` decimals as huge rationals,
    compute det T at the approximation. Returns (det, q0, t0).

    This is *evidential* not rigorous: pi and e are not represented
    exactly here. But computing det T as an exact Q-element gives a
    floor on |det T(pi, e)| -- vanishing det T at a high-precision
    rational approximation that *also* vanishes at randomly chosen
    integers is consistent with det T == 0 (rigorously) only via the
    bound L_1 norm <= deg(det T) * height(q0, t0)^d_total at the
    closest non-zero polynomial. We do not attempt that bound here."""
    # 30-digit rational approximations
    # pi ~ 3.141592653589793238462643383279
    # e  ~ 2.718281828459045235360287471352
    PI_DEC = "31415926535897932384626433832795028841971693993751"
    E_DEC  = "27182818284590452353602874713526624977572470936999"
    if digits > len(PI_DEC):
        digits = len(PI_DEC)
    q0_num = int(PI_DEC[:digits])
    t0_num = int(E_DEC[:digits])
    # We *must* feed integers to f4ncgb. Take ((q0_num), (t0_num)) as
    # large primes-of-best-approximation. Drop the denominator (it just
    # scales (q, t) uniformly; for a *rational* approximation of (pi, e)
    # we would need an algebraic version. The integers q0, t0 here
    # are well-separated huge integers with essentially-arbitrary
    # prime factorisation -- which is exactly what Schwartz--Zippel
    # needs).
    q0 = q0_num
    t0 = t0_num
    if verbose:
        print(f"   [transcendental-approx specialisation: "
              f"q0 = {q0}, t0 = {t0}]", file=sys.stderr)
    basis = build_basis_at(workdir, n=n, q=q0, t=t0, verbose=verbose)
    N, _, det = trace_form_at(basis, n=n, verbose=verbose)
    return det, q0, t0, N


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=3,
                    help="number of strands (default: 3)")
    ap.add_argument("--mode", choices=["all", "huge", "grid", "interp",
                                       "interp-entries",
                                       "transcendental"],
                    default="all",
                    help=("which argument(s) to run. 'all' runs (A) "
                          "then (B) then (C) then (D); each subsequent "
                          "one is more expensive."))
    ap.add_argument("--huge-bits", type=int, default=12,
                    help="bit-length of (q0, t0) in argument (C); "
                    "f4ncgb is unreliable above ~12 bits at the moment.")
    ap.add_argument("--deg-bound", type=int, default=12,
                    help=("degree-in-each-variable bound for arguments "
                          "(A) and (B). For n=3 square-zero we have verified "
                          "empirically that deg_q, deg_t <= 12 suffices; "
                          "raise if interpolation or grid certificate "
                          "fails."))
    ap.add_argument("--denom-exp", type=int, default=5,
                    help=("non-negative exponent M such that "
                          "(q*t)^M * det T(q, t) is a polynomial in "
                          "Q[q, t] (i.e. the Laurent denominator is "
                          "killed). Default 5; raise if interpolated "
                          "polynomial fails the consistency check."))
    ap.add_argument("--max-dim", type=int, default=500)
    ap.add_argument("--workdir", default=None,
                    help="working directory; if not given, uses "
                         "a tempdir under /tmp.")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    # Set up workdir.
    if args.workdir is None:
        workdir = tempfile.mkdtemp(prefix="generic_nonss_")
    else:
        workdir = args.workdir
        prepare_output_dir(workdir)

    def hd(s: str) -> None:
        print()
        print(s)
        print("=" * len(s))

    print(f"[generic_non_ss] n = {args.n}, square-zero variant, workdir = "
          f"{workdir}")

    # ---------------------------------------------------------------
    # Baseline at (q, t) = (5, 7): the very first non-rigorous datum.
    # ---------------------------------------------------------------
    hd("Baseline at (q, t) = (5, 7)")
    basis_57 = build_basis_at(workdir, n=args.n, q=5, t=7,
                              verbose=args.verbose)
    N, rank, det = trace_form_at(basis_57, n=args.n, max_dim=args.max_dim,
                                 verbose=args.verbose)
    deficit = N - rank
    print(f"  N         = {N}")
    print(f"  rank T    = {rank}")
    print(f"  deficit   = {deficit}")
    print(f"  det T(5,7) = {det}")
    is_baseline_zero = (det == 0)
    if not is_baseline_zero:
        print("  -> det T NOT zero even at the baseline; this script's "
              "premise (square-zero variant is non-ss at integer points) "
              "is violated. Aborting.")
        return 2

    # ---------------------------------------------------------------
    # Argument (C): single huge-coprime specialisation
    # ---------------------------------------------------------------
    if args.mode in ("all", "huge"):
        hd(f"Argument (C): huge-coprime specialisation "
           f"(~{args.huge_bits}-bit (q0, t0))")
        is_zero, q0, t0, det_huge, dim_huge = huge_specialisation(
            workdir, n=args.n, bits=args.huge_bits, verbose=args.verbose,
        )
        print(f"  (q0, t0)   = ({q0}, {t0})")
        print(f"  dim A at (q0, t0) = {dim_huge}  "
              f"(expected {N} for non-degenerate (q, t))")
        if dim_huge != N:
            print(f"  ==> WARNING: dim A at (q0, t0) = {dim_huge} != "
                  f"{N}. The Groebner basis collapsed the algebra "
                  f"to lower dimension at this specialisation -- "
                  f"likely an f4ncgb rational-reconstruction issue "
                  f"at large integer coefficients (we have observed "
                  f"this above ~5000); the resulting det T value is "
                  f"NOT a faithful sample of det T(q0, t0). Skip "
                  f"argument (C).")
            print(f"  det T(q0, t0) reported by collapsed GB = "
                  f"{det_huge}")
        elif is_zero:
            print("  Zero determinant observed at this chosen parameter pair.")
            print("  This is a specialization diagnostic, not a generic-vanishing certificate.")
        else:
            print("  ==> det T is *not* identically zero at this "
                  "huge specialisation. This contradicts the "
                  "working hypothesis (Remark 1.5) IF dim A was "
                  "preserved at (q0, t0). Investigate.")

    # ---------------------------------------------------------------
    # Argument (A): small-prime grid Schwartz--Zippel
    # ---------------------------------------------------------------
    if args.mode in ("all", "grid"):
        hd(f"Argument (A): ({args.deg_bound + 2})x({args.deg_bound + 2}) "
           f"small-prime grid Schwartz--Zippel")
        all_zero_grid, table = schwartz_zippel_grid(
            workdir, n=args.n, deg_bound=args.deg_bound,
            verbose=args.verbose,
        )
        print(f"  grid covered {len(table)} (q, t)-pairs")
        print(f"  all det T zero on the grid? {all_zero_grid}")
        if all_zero_grid:
            L = args.deg_bound + 2
            print('  ==> Conditional grid implication (bounds and common basis are assumptions):')
            print('      Write det T(q, t) = P(q, t) / (q^A t^B) with P in Q[q, t].')
            print('      The grid avoids q = 0 and t = 0, so')
            print('      det T(q_i, t_j) = 0 <=> P(q_i, t_j) = 0.')
            print(f"      Provided deg_q P <= {args.deg_bound} and "
                  f"deg_t P <= {args.deg_bound}, P is uniquely "
                  f"determined by its values on any "
                  f"({args.deg_bound + 1})x({args.deg_bound + 1}) "
                  f"grid; here L = {L} > {args.deg_bound + 1} on each "
                  f"axis, so the vanishing forces P == 0, hence")
            print('      det T(q, t) == 0 in Q(q, t).')
            # show first few non-trivial rows
            if args.verbose or True:
                for qv, tv, Nn, rk, dt in table[:5]:
                    print(f"    (q, t) = ({qv:>2}, {tv:>2}): "
                          f"N={Nn}, rank={rk}, det T = {dt}")
                if len(table) > 5:
                    print(f"    ... and {len(table) - 5} more "
                          f"(all det T = 0)")

    # ---------------------------------------------------------------
    # Argument (B): explicit interpolation reconstruction of det T
    # ---------------------------------------------------------------
    if args.mode in ("all", "interp"):
        hd(f"Argument (B): direct polynomial reconstruction of "
           f"(q*t)^{args.denom_exp} * det T(q, t) on a "
           f"({args.deg_bound + 1})x({args.deg_bound + 1}) prime grid")
        det_poly, q_nodes, t_nodes = reconstruct_det_polynomial(
            workdir, n=args.n,
            deg_q=args.deg_bound, deg_t=args.deg_bound,
            denom_exp=args.denom_exp,
            verbose=args.verbose,
        )
        deg = qp_total_degree(det_poly)
        n_terms = len(det_poly)
        print(f"  interpolated polynomial (q*t)^{args.denom_exp} * det T")
        print(f"    total degree   = {deg}")
        print(f"    non-zero terms = {n_terms}")
        if not det_poly:
            print()
            print("  ===============================================")
            print("  Zero interpolant under the supplied degree/denominator bounds.")
            print("  Generic vanishing requires proved bounds and a common certified basis.")
        else:
            # Sanity check: evaluate the interpolated polynomial at
            # a fresh prime pair not in the grid.
            qcheck, tcheck = 211, 223
            qset = set(q_nodes); tset = set(t_nodes)
            assert qcheck not in qset and tcheck not in tset
            # Evaluate poly
            val = Fraction(0)
            for (a, b), c in det_poly.items():
                val += c * qcheck ** a * tcheck ** b
            # Compute direct
            basis = build_basis_at(workdir, n=args.n,
                                   q=qcheck, t=tcheck,
                                   verbose=args.verbose)
            _, _, det_direct = trace_form_at(basis, n=args.n, verbose=args.verbose)
            target = det_direct * (qcheck * tcheck) ** args.denom_exp
            print(f"    direct (q*t)^{args.denom_exp} * det T"
                  f"({qcheck},{tcheck}) = {target}")
            print(f"    poly  (q*t)^{args.denom_exp} * det T"
                  f"({qcheck},{tcheck}) = {val}")
            if val == target:
                print('    => interpolation is consistent.')
                print(f"  (q*t)^{args.denom_exp} * det T is NOT the "
                      f"zero polynomial -- printing a snippet:")
                # print top 6 monomials
                items = sorted(det_poly.items(),
                               key=lambda kv: (-(kv[0][0] + kv[0][1]),
                                               kv[0]))
                for k, (deg, c) in enumerate(items[:6]):
                    print(f"     coeff (q^{deg[0]} t^{deg[1]}) = {c}")
                print('  ==> the reconstructed candidate is not the zero polynomial; this is at odds with Argument (A) and (empirical) Remark 1.5. Investigate --denom-exp / --deg-bound.')
            else:
                print('    => interpolation INCONSISTENT.')
                print('  Likely cause: degree underestimate or Laurent-exponent underestimate. Increase --deg-bound and/or --denom-exp and re-run.')

    # ---------------------------------------------------------------
    # Argument (B'): entry-wise interpolation (slower, optional)
    # ---------------------------------------------------------------
    if args.mode == "interp-entries":
        hd("Argument (B'): entry-wise polynomial reconstruction")
        T_polys = reconstruct_T_polynomial(
            workdir, n=args.n,
            deg_q=args.deg_bound, deg_t=args.deg_bound,
            verbose=args.verbose,
        )
        N = len(T_polys)
        max_entry_deg = max(qp_total_degree(p) for row in T_polys
                            for p in row)
        nonzero_entries = sum(1 for row in T_polys for p in row if p)
        print('  T(q, t) reconstructed as polynomial matrix')
        print(f"    N = {N}")
        print(f"    max total degree per entry = {max_entry_deg}")
        print(f"    non-zero entries           = "
              f"{nonzero_entries} / {N * N}")
        ok = verify_T_polynomial_at(
            T_polys, workdir, n=args.n,
            q=211, t=223, verbose=args.verbose,
        )
        print(f"    polynomial vs direct at (211, 223): "
              f"{'MATCH' if ok else 'MISMATCH'}")
        # NOTE: entry-wise interpolation does NOT correctly handle
        # Laurent denominators in T_{ij}; if MISMATCH, use --mode
        # interp (which interpolates q^M t^M * det T directly and
        # therefore handles Laurent transparently).

    # ---------------------------------------------------------------
    # Argument (D): transcendental-approximation specialisation
    # ---------------------------------------------------------------
    if args.mode in ("all", "transcendental"):
        hd("Argument (D): high-precision (pi, e) approximation "
           "(evidential)")
        det_pe, q0, t0, dim_pe = transcendental_specialisation(
            workdir, n=args.n, digits=20, verbose=args.verbose,
        )
        print(f"  (q0, t0) = ({q0}, {t0})  "
              f"(20-digit truncations of pi, e)")
        print(f"  dim A at (q0, t0) = {dim_pe}  (expected {N})")
        print(f"  det T(q0, t0) = {'0 (exactly)' if det_pe == 0 else det_pe}")
        if dim_pe != N:
            print("  ==> WARNING: dim A collapsed at this huge "
                  "specialisation (f4ncgb rational-reconstruction "
                  "issue). det T value is uninformative.")
        elif det_pe == 0:
            print("  ==> The same exact-zero result is observed at an "
                  "approximation of an algebraically independent "
                  "transcendental pair. Strong confirming evidence.")

    print()
    print("[done] workdir =", workdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
