"""Per-point algebra: emit + f4ncgb -> Groebner basis (exact Q with char=0, or a large prime
field with char=p); standard-monomial basis + memoized normal-form reducer; and the left-regular
TRACE FORM at the point.

G_{ij} = Tr(L_{m_i m_j}) on the standard basis {m_k}: needs the products NF(m_i m_j) (structure
constants). tau(m_u) = Tr(L_{m_u}) = sum_j [coeff of m_j in NF(m_u m_j)]; G_{ij} = <NF(m_i m_j),tau>.

The trace-form RANK is an exact integer invariant. Computing it over a large prime field GF(p) at a
generic (q,t) equals the rank over Q for a good prime (cross-checked against exact-Q at small dim);
this is what makes the survey scale (Q Fraction arithmetic blows up past dim ~200).
"""
from __future__ import annotations
import os, subprocess
from fractions import Fraction as Fr
from f4ncgb_tools import core
from f4ncgb_tools.research import theta_min_poly as _T
from .spec import Spec

F4 = os.environ.get("F4NCGB", "/home/choom/src/f4ncgb/build-release/f4ncgb")
_HERE = os.path.dirname(os.path.abspath(__file__))
EMIT = os.path.normpath(f"{_HERE}/../zipper_in_ms.py")


def point_gb(spec: Spec, q: int, t: int, scratch: str, char: int = 0, maxiter: int = 200000):
    """Emit at (q,t) and return the reduced GB. char=0 -> exact Q (parse_basis_Q); char=p -> GF(p)."""
    tag = f"{spec.name}-{q}-{t}-c{char}"
    ms = os.path.join(scratch, f"pa-{tag}.ms"); bs = ms[:-3] + ".basis"
    subprocess.run(["python3", EMIT, "--n", str(spec.n), "--untwist", spec.emit_untwist,
                    "--framing-power", str(spec.framing_power), "--q", str(q), "--t", str(t),
                    "--char", str(char), "--max-untwist-level", str(spec.max_untwist_level),
                    "--out", ms], check=True, capture_output=True)
    subprocess.run([F4, "-m", str(maxiter), "-i", ms, "-o", bs], check=True, capture_output=True)
    return core.parse_basis_Q(bs) if char == 0 else core.parse_basis_Fp(bs)


class PointAlgebra:
    """Algebra at (q,t) over Q (char=0) or GF(char): standard basis, reducer, leading monomials."""
    def __init__(self, spec: Spec, q: int, t: int, scratch: str, char: int = 0):
        self.spec, self.q, self.t, self.char = spec, q, t, char
        polys = point_gb(spec, q, t, scratch, char)
        self.lms = [lm for lm, _ in polys]
        ring = "Q" if char == 0 else "Fp"
        self.dim, self.nf, self.standard, self._zero, self._norm = _T.make_reducer(polys, ring)

    def product_table(self):
        """prod[i][j] = NF(m_i m_j) as {basis_idx: coef in the ring}."""
        std, nf, d = self.standard, self.nf, self.dim
        return [[nf(std[i] + std[j]) for j in range(d)] for i in range(d)]

    def trace_form(self):
        """Return the trace-form Gram matrix G (dim x dim) over the ring."""
        d, zero, norm = self.dim, self._zero, self._norm
        prod = self.product_table()
        tau = [norm(sum((prod[u][j].get(j, zero) for j in range(d)), zero)) for u in range(d)]
        G = [[norm(sum((c * tau[idx] for idx, c in prod[i][j].items()), zero))
              for j in range(d)] for i in range(d)]
        return G
