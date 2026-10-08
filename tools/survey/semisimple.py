"""Semisimplicity via the trace form (char-0 Dieudonne criterion): the algebra is semisimple iff
the left-regular trace form is nondegenerate. Generic (over Q(q,t)) invariants are read at random
(q,t) points (Schwartz-Zippel: generic rank = rank at a random point).

The trace-form RANK is an exact integer invariant. Two computation modes:
  * char=p (default, SCALE): build the trace form over GF(p) at (q,t); rank over GF(p). Fast (no
    Fraction blowup). Rigorous for the rank of the exactly-defined form at a generic point/good
    prime; cross-checked against exact-Q at small dim.
  * char=0 (VERIFY): build the trace form over Q (exact Fractions); rank via reduction mod two
    large primes (must agree). The load-bearing reference; slow past dim ~200.
"""
from __future__ import annotations
from fractions import Fraction as Fr
from .pointalg import PointAlgebra
from .spec import Spec

_PRIMES = (2147483647, 2147483629)
DEFAULT_CHAR = 2147483647


def _rank_Fp(rows, p):
    """Rank over F_p of an integer matrix (entries already reduced mod p)."""
    n = len(rows); m = len(rows[0]) if n else 0
    A = [[v % p for v in r] for r in rows]
    rank = 0
    for col in range(m):
        piv = next((r for r in range(rank, n) if A[r][col] % p), None)
        if piv is None:
            continue
        A[rank], A[piv] = A[piv], A[rank]
        inv = pow(A[rank][col], p - 2, p)
        A[rank] = [(v * inv) % p for v in A[rank]]
        for r in range(n):
            if r != rank and A[r][col] % p:
                f = A[r][col]
                A[r] = [(A[r][c] - f * A[rank][c]) % p for c in range(m)]
        rank += 1
        if rank == n:
            break
    return rank


def rank_of(G, char):
    """Exact rank of the trace form. char=0: reduce mod two primes (must agree). char=p: rank/GF(p)."""
    if char == 0:
        ranks = set()
        for p in _PRIMES:
            # Reducing an exact rational mod p is undefined when p divides its denominator.
            # pow(0, p-2, p) == 0 would silently turn such an entry into 0 and the rank could
            # only come out too small -- the safe direction for a semisimplicity certificate,
            # but a wrong number with no error. Refuse instead (notes/reptheory_survey/
            # qqt-strategy.md H1).
            for row in G:
                for e in row:
                    if e.denominator % p == 0:
                        raise ValueError(f"trace-form entry {e} has denominator divisible by "
                                         f"p={p}; choose another prime")
            rows = [[(e.numerator % p) * pow(e.denominator % p, p - 2, p) % p for e in row] for row in G]
            ranks.add(_rank_Fp(rows, p))
        if len(ranks) != 1:
            raise RuntimeError(f"prime disagreement on rank: {ranks}")
        return ranks.pop()
    return _rank_Fp(G, char)


def semisimplicity(spec: Spec, scratch: str, points=((5, 7), (13, 17)), char: int = DEFAULT_CHAR):
    """dim + generic rank of the trace form at the given generic points.
    Returns dict: dim, rank, radical_dim, semisimple, symmetric, points, char."""
    results = []
    dim = None
    for (q, t) in points:
        A = PointAlgebra(spec, q, t, scratch, char=char)
        G = A.trace_form()
        symm = all(G[i][j] == G[j][i] for i in range(A.dim) for j in range(A.dim))
        r = rank_of(G, char)
        results.append((q, t, A.dim, r, symm))
        dim = A.dim
    rank = max(r for *_, r, _ in results)  # generic rank = max over sampled generic points
    return {
        "spec": spec.name, "dim": dim, "rank": rank, "radical_dim": dim - rank,
        "semisimple": rank == dim, "symmetric": all(s for *_, s in results),
        "points": [(q, t, d, r) for (q, t, d, r, _) in results], "char": char,
    }
