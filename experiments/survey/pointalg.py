"""Specialized reduction and trace-matrix observations from solver output.

The supplied polynomials have no completeness certificate here. Standard-word
counts and normal forms are therefore observations of this reduction model,
not certified invariants of the presented quotient or its generic fibre.
"""
from __future__ import annotations
import os, subprocess
from zipper_algebra import core
from zipper_algebra.linear import validate_field_characteristic
from experiments.research import theta_min_poly as _T
from .spec import Spec

F4 = os.environ.get("F4NCGB", "f4ncgb")


def point_gb(spec: Spec, q: int, t: int, scratch: str, char: int = 0, maxiter: int = 200000):
    """Emit the presentation and parse solver output; completion is uncertified."""
    validate_field_characteristic(char)
    tag = f"{spec.name}-{q}-{t}-c{char}"
    ms = os.path.join(scratch, f"pa-{tag}.ms"); bs = ms[:-3] + ".basis"
    spec.presentation().emit_ms(ms, q, t, char)
    subprocess.run([F4, "-m", str(maxiter), "-i", ms, "-o", bs], check=True, capture_output=True)
    return core.parse_basis_Q(bs) if char == 0 else core.parse_basis_Fp(bs)


class PointAlgebra:
    """Specialized reduction model from an uncertified solver basis candidate."""
    def __init__(self, spec: Spec, q: int, t: int, scratch: str, char: int = 0):
        validate_field_characteristic(char)
        self.spec, self.q, self.t, self.char = spec, q, t, char
        self.basis_complete = False
        polys = point_gb(spec, q, t, scratch, char)
        self.lms = [lm for lm, _ in polys]
        ring = "Q" if char == 0 else "Fp"
        self.dim, self.nf, self.standard, self._zero, self._norm = _T.make_reducer(
            polys, ring, char=char, variables=spec.presentation().var_names())

    def product_table(self):
        """prod[i][j] = NF(m_i m_j) as {basis_idx: coef in the ring}."""
        std, nf, d = self.standard, self.nf, self.dim
        return [[nf(std[i] + std[j]) for j in range(d)] for i in range(d)]

    def trace_form(self):
        """Compute the trace matrix of the supplied standard-word reduction model."""
        d, zero, norm = self.dim, self._zero, self._norm
        prod = self.product_table()
        tau = [norm(sum((prod[u][j].get(j, zero) for j in range(d)), zero)) for u in range(d)]
        G = [[norm(sum((c * tau[idx] for idx, c in prod[i][j].items()), zero))
              for j in range(d)] for i in range(d)]
        return G
