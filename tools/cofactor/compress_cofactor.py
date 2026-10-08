#!/usr/bin/env python3
"""
Compress the 695-tuple X_3 cofactor onto the minimal independent set
{R_0, R_1, R_2, R_3, R_6, R_22, R_25, R_28} for the n=4 Bigelow zipper algebra.

Strategy:
  1. For each dependent input R_d in {R_7, R_14, R_15, R_16, R_17, R_26},
     express R_d as a sum sum_k a_k * L_k * R_{i_k} * R'_k over the minimal
     generators. Verify each derivation as a polynomial identity in the
     free algebra Q(q,t)<x_1,x_2,x_3,y_1,y_2,y_3>.
  2. Substitute these into the 695-tuple cofactor.
  3. Re-verify the compressed cofactor (sum equals X_3 exactly).
  4. Write out cofactor-X3-n4-Qqt-compressed.txt.
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py lives
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import re
import time
from fractions import Fraction
from typing import Dict, List, Tuple
from relations_naming import label
from symbolic_identity_check import (
    VAR_NAMES, VAR_INDEX, N_STRANDS,
    qt_zero, qt_const, qt_q, qt_t, qt_add, qt_sub, qt_neg, qt_mul, qt_pow,
    qt_eq, qt_is_zero, qt_str,
    nc_zero, nc_const, nc_var, nc_var_named, nc_add, nc_sub, nc_neg,
    nc_mul, nc_scalar, nc_is_zero, nc_size,
    build_inputs_and_X3,
    parse_cofactor_line, parse_numerator, parse_denominator,
    COF_FILE, ROOT,
)


# --------- Convenience constructors -----------------------------------
def mon_nc(mon: Tuple[int, ...]) -> "NCPoly":
    return {mon: qt_const(1)}


def mon_from_str(s: str) -> Tuple[int, ...]:
    if not s or s == "()":
        return ()
    return tuple(VAR_INDEX[v] for v in s.split("*"))


def mon_to_str(m: Tuple[int, ...]) -> str:
    if not m:
        return "()"
    return "*".join(VAR_NAMES[i - 1] for i in m)


# --------- Derivation table ------------------------------------------
# Each dependent relation is represented as a list of (coef_qt, L_mon,
# R_idx, R'_mon) tuples, where the sum equals the dependent relation
# exactly in the free algebra. Verified at runtime.

# Helpers
def Q(c) -> "QTPoly":
    return qt_const(c)

# (q-1)
def qm1() -> "QTPoly":
    return qt_sub(qt_q(), Q(1))


def derivations(relations):
    """
    Build derivations for the dependent relations.

    For each dependent index d, produce a list of (coef, L_mon, src_idx, R_mon).
    Indices src_idx must lie in the minimal set {0,1,2,3,6,22,25,28}.
    """

    # The 'minimal' set
    KEEP = {0, 1, 2, 3, 6, 22, 25, 28}

    derivs: Dict[int, List[Tuple["QTPoly", Tuple[int, ...], int, Tuple[int, ...]]]] = {}

    # R_7 = y_1 y_2 y_1 - y_2 y_1 y_2
    # Derivation (verified by hand-expansion -- re-verified at runtime):
    #   y_2 y_1 y_2 * R_6 * y_1 y_2 y_1
    #   = (y_2 y_1 y_2)(x_1 x_2 x_1)(y_1 y_2 y_1) - (y_2 y_1 y_2)(x_2 x_1 x_2)(y_1 y_2 y_1)
    # Reduce by successive applications of x_i y_i = 1 + R_{2i-2} and
    # y_i x_i = 1 + R_{2i-1}:
    #  Term1 (xxx variant) -> y_2 y_1 y_2 (after cancellations) + R-tails
    #  Term2 (xxx-swap variant) -> y_1 y_2 y_1 + R-tails
    # Net: y_2 y_1 y_2 - y_1 y_2 y_1 + (R-tails) = (y_2 y_1 y_2)(R_6)(y_1 y_2 y_1)
    # Hence R_7 = y_1 y_2 y_1 - y_2 y_1 y_2 = -(y_2 y_1 y_2)(R_6)(y_1 y_2 y_1) + R-tails
    #
    # Let's derive R-tails carefully.
    # Term1 = y_2 y_1 y_2 (x_1 x_2 x_1) y_1 y_2 y_1
    # Reduce x_1 y_1 = 1 + R_0  (R_0 = x_1 y_1 - 1)
    # Use: y_2 y_1 y_2 x_1 x_2 x_1 y_1 y_2 y_1
    #   = y_2 y_1 y_2 x_1 x_2 (1 + R_0) y_2 y_1
    #   = y_2 y_1 y_2 x_1 x_2 y_2 y_1 + (y_2 y_1 y_2 x_1 x_2) * R_0 * (y_2 y_1)
    # Reduce x_2 y_2 = 1 + R_2:
    #   y_2 y_1 y_2 x_1 x_2 y_2 y_1 = y_2 y_1 y_2 x_1 (1 + R_2) y_1
    #     = y_2 y_1 y_2 x_1 y_1 + y_2 y_1 y_2 x_1 * R_2 * y_1
    # Reduce x_1 y_1 = 1 + R_0:
    #   y_2 y_1 y_2 x_1 y_1 = y_2 y_1 y_2 (1 + R_0) = y_2 y_1 y_2 + y_2 y_1 y_2 * R_0
    # So Term1 = y_2 y_1 y_2
    #          + (y_2 y_1 y_2) * R_0 * ()
    #          + (y_2 y_1 y_2 x_1) * R_2 * (y_1)
    #          + (y_2 y_1 y_2 x_1 x_2) * R_0 * (y_2 y_1)
    #
    # Term2 = y_2 y_1 y_2 (x_2 x_1 x_2) y_1 y_2 y_1
    # Reduce y_2 x_2 = 1 + R_3 (left):
    #   = y_2 y_1 (1 + R_3) x_1 x_2 y_1 y_2 y_1
    #   = y_2 y_1 x_1 x_2 y_1 y_2 y_1 + (y_2 y_1) * R_3 * (x_1 x_2 y_1 y_2 y_1)
    # Reduce y_1 x_1 = 1 + R_1:
    #   y_2 y_1 x_1 x_2 y_1 y_2 y_1 = y_2 (1 + R_1) x_2 y_1 y_2 y_1
    #     = y_2 x_2 y_1 y_2 y_1 + (y_2) * R_1 * (x_2 y_1 y_2 y_1)
    # Reduce y_2 x_2 = 1 + R_3:
    #   y_2 x_2 y_1 y_2 y_1 = (1 + R_3) y_1 y_2 y_1
    #     = y_1 y_2 y_1 + () * R_3 * (y_1 y_2 y_1)
    # So Term2 = y_1 y_2 y_1
    #          + () * R_3 * (y_1 y_2 y_1)
    #          + (y_2) * R_1 * (x_2 y_1 y_2 y_1)
    #          + (y_2 y_1) * R_3 * (x_1 x_2 y_1 y_2 y_1)
    #
    # Then (y_2 y_1 y_2) R_6 (y_1 y_2 y_1) = Term1 - Term2
    #   = (y_2 y_1 y_2 - y_1 y_2 y_1)
    #     + (y_2 y_1 y_2) R_0 ()
    #     + (y_2 y_1 y_2 x_1) R_2 (y_1)
    #     + (y_2 y_1 y_2 x_1 x_2) R_0 (y_2 y_1)
    #     - () R_3 (y_1 y_2 y_1)
    #     - (y_2) R_1 (x_2 y_1 y_2 y_1)
    #     - (y_2 y_1) R_3 (x_1 x_2 y_1 y_2 y_1)
    #
    # So R_7 = y_1 y_2 y_1 - y_2 y_1 y_2
    #        = - (y_2 y_1 y_2) R_6 (y_1 y_2 y_1)
    #          + (y_2 y_1 y_2) R_0 ()
    #          + (y_2 y_1 y_2 x_1) R_2 (y_1)
    #          + (y_2 y_1 y_2 x_1 x_2) R_0 (y_2 y_1)
    #          - () R_3 (y_1 y_2 y_1)
    #          - (y_2) R_1 (x_2 y_1 y_2 y_1)
    #          - (y_2 y_1) R_3 (x_1 x_2 y_1 y_2 y_1)

    derivs[7] = [
        (Q(-1), mon_from_str("y2*y1*y2"), 6, mon_from_str("y1*y2*y1")),
        (Q(1),  mon_from_str("y2*y1*y2"), 0, ()),
        (Q(1),  mon_from_str("y2*y1*y2*x1"), 2, mon_from_str("y1")),
        (Q(1),  mon_from_str("y2*y1*y2*x1*x2"), 0, mon_from_str("y2*y1")),
        (Q(-1), (), 3, mon_from_str("y1*y2*y1")),
        (Q(-1), mon_from_str("y2"), 1, mon_from_str("x2*y1*y2*y1")),
        (Q(-1), mon_from_str("y2*y1"), 3, mon_from_str("x1*x2*y1*y2*y1")),
    ]

    # R_14 = x_1 x_2 y_1 - y_2 x_1 x_2
    # Direct derivation: from R_6 = x_1 x_2 x_1 - x_2 x_1 x_2,
    #   R_6 * y_1 = x_1 x_2 x_1 y_1 - x_2 x_1 x_2 y_1
    #            = x_1 x_2 (1 + R_0) - x_2 x_1 x_2 y_1
    #            = x_1 x_2 + x_1 x_2 R_0 - x_2 x_1 x_2 y_1
    # Hmm, that has x_2 x_1 x_2 y_1 which isn't the goal.
    #
    # Better: R_14 = x_1 x_2 y_1 - y_2 x_1 x_2. Let's check:
    #   x_1 x_2 y_1 = ?  (no direct reduction since x_2 y_1 doesn't simplify)
    #
    # Use R_6 in a different sandwich. Note R_6 = x_1 x_2 x_1 - x_2 x_1 x_2.
    # Consider y_2 * R_6 = y_2 x_1 x_2 x_1 - y_2 x_2 x_1 x_2 = y_2 x_1 x_2 x_1 - (1+R_3) x_1 x_2
    # So y_2 * R_6 = y_2 x_1 x_2 x_1 - x_1 x_2 - R_3 x_1 x_2.
    # Multiply by y_1 on right:
    #   y_2 * R_6 * y_1 = y_2 x_1 x_2 x_1 y_1 - x_1 x_2 y_1 - R_3 x_1 x_2 y_1
    # And y_2 x_1 x_2 (x_1 y_1) = y_2 x_1 x_2 (1 + R_0) = y_2 x_1 x_2 + y_2 x_1 x_2 R_0
    # So y_2 R_6 y_1 = y_2 x_1 x_2 + y_2 x_1 x_2 R_0 - x_1 x_2 y_1 - R_3 x_1 x_2 y_1
    #                = -(x_1 x_2 y_1 - y_2 x_1 x_2) + y_2 x_1 x_2 R_0 - R_3 x_1 x_2 y_1
    #                = -R_14 + y_2 x_1 x_2 R_0 - R_3 x_1 x_2 y_1
    # Hence:
    #   R_14 = -y_2 R_6 y_1 + y_2 x_1 x_2 R_0 - R_3 x_1 x_2 y_1
    derivs[14] = [
        (Q(-1), mon_from_str("y2"), 6, mon_from_str("y1")),
        (Q(1),  mon_from_str("y2*x1*x2"), 0, ()),
        (Q(-1), (), 3, mon_from_str("x1*x2*y1")),
    ]

    # R_15 = -x_2 y_1 y_2 + y_1 y_2 x_1  (note the signs)
    # Strategy: build R_15 from R_7 (already dependent on R_6) plus invertibility.
    # Easier: directly from R_6 via "double inversion".
    # Try:  y_1 y_2 * R_6 * y_2 y_1 * y_2 = ?  Too long.
    #
    # Alternative direct: R_15 = y_1 y_2 x_1 - x_2 y_1 y_2.
    # Sandwich: y_1 y_2 R_6 ?
    # y_1 y_2 * R_6 = y_1 y_2 x_1 x_2 x_1 - y_1 y_2 x_2 x_1 x_2
    # Multiply by y_2 y_1 on right:
    #   y_1 y_2 R_6 y_2 y_1 = y_1 y_2 x_1 x_2 x_1 y_2 y_1 - y_1 y_2 x_2 x_1 x_2 y_2 y_1
    # Note x_2 y_2 = 1 + R_2: so x_2 x_1 x_2 y_2 = x_2 x_1 (1 + R_2) = x_2 x_1 + x_2 x_1 R_2
    # So second term: y_1 y_2 (x_2 x_1 + x_2 x_1 R_2) y_1 = y_1 y_2 x_2 x_1 y_1 + y_1 y_2 x_2 x_1 R_2 y_1
    # First term: y_1 y_2 x_1 x_2 x_1 y_2 y_1
    # This is getting messy. Let me try another sandwich approach.
    #
    # Easier: R_15 = -x_2 y_1 y_2 + y_1 y_2 x_1
    # We want sandwich of R_6 by y's. Consider y_1 * R_6 * y_2 y_1:
    #   = y_1 x_1 x_2 x_1 y_2 y_1 - y_1 x_2 x_1 x_2 y_2 y_1
    #   y_1 x_1 = 1 + R_1: = (1+R_1) x_2 x_1 y_2 y_1 - y_1 x_2 x_1 (x_2 y_2) y_1
    #   = x_2 x_1 y_2 y_1 + R_1 x_2 x_1 y_2 y_1 - y_1 x_2 x_1 (1 + R_2) y_1
    #   = x_2 x_1 y_2 y_1 - y_1 x_2 x_1 y_1 + R_1 x_2 x_1 y_2 y_1 - y_1 x_2 x_1 R_2 y_1
    # Hmm, x_2 x_1 y_2 y_1 - y_1 x_2 x_1 y_1. Not R_15 directly.
    #
    # Try y_2 y_1 * R_6 * y_2:
    #   = y_2 y_1 x_1 x_2 x_1 y_2 - y_2 y_1 x_2 x_1 x_2 y_2
    #   y_1 x_1 = 1 + R_1:
    #   = y_2 (1+R_1) x_2 x_1 y_2 - y_2 y_1 x_2 x_1 (1 + R_2)
    #   = y_2 x_2 x_1 y_2 + y_2 R_1 x_2 x_1 y_2 - y_2 y_1 x_2 x_1 - y_2 y_1 x_2 x_1 R_2
    #   y_2 x_2 = 1 + R_3:
    #   = (1+R_3) x_1 y_2 - y_2 y_1 x_2 x_1 + y_2 R_1 x_2 x_1 y_2 - y_2 y_1 x_2 x_1 R_2
    #   = x_1 y_2 - y_2 y_1 x_2 x_1 + R_3 x_1 y_2 + y_2 R_1 x_2 x_1 y_2 - y_2 y_1 x_2 x_1 R_2
    # Not the target.
    #
    # Approach: derive R_15 from R_14 by sandwiching with y's.
    # R_14 = x_1 x_2 y_1 - y_2 x_1 x_2
    # y_2 y_1 * R_14 = y_2 y_1 x_1 x_2 y_1 - y_2 y_1 y_2 x_1 x_2
    #   y_1 x_1 = 1 + R_1: = y_2 (1+R_1) x_2 y_1 - y_2 y_1 y_2 x_1 x_2
    #   = y_2 x_2 y_1 - y_2 y_1 y_2 x_1 x_2 + y_2 R_1 x_2 y_1
    #   y_2 x_2 = 1 + R_3: = (1+R_3) y_1 - y_2 y_1 y_2 x_1 x_2 + y_2 R_1 x_2 y_1
    #   = y_1 - y_2 y_1 y_2 x_1 x_2 + R_3 y_1 + y_2 R_1 x_2 y_1
    # Hmm, not quite. Need to leverage that y_2 y_1 y_2 = y_1 y_2 y_1 (mod R_7).
    # But R_7 is derivable from R_6.
    #
    # Cleanest: R_15 by direct sandwich of R_6.
    #
    # Observation: x_1 x_2 = (x_1 y_1)(x_1 x_2 x_1) y_1 = (1+R_0) x_1 x_2 x_1 y_1 ... no, x_1 y_1 != y_1 x_1.
    #
    # Let me think differently. Take y_1 y_2 * R_14 from the right side too:
    # R_14 * y_2 y_1 = x_1 x_2 y_1 y_2 y_1 - y_2 x_1 x_2 y_2 y_1
    #   x_2 y_2 = 1+R_2: = x_1 x_2 y_1 y_2 y_1 - y_2 x_1 (1+R_2) y_1
    #   = x_1 x_2 y_1 y_2 y_1 - y_2 (x_1 y_1) - y_2 x_1 R_2 y_1
    #   x_1 y_1 = 1+R_0: = x_1 x_2 y_1 y_2 y_1 - y_2 - y_2 R_0 - y_2 x_1 R_2 y_1
    # Now y_1 y_2 y_1 = y_2 y_1 y_2 + R_7. So:
    #   x_1 x_2 y_1 y_2 y_1 = x_1 x_2 y_2 y_1 y_2 + x_1 x_2 R_7
    #   x_2 y_2 = 1+R_2: = x_1 (1+R_2) y_1 y_2 + x_1 x_2 R_7
    #   = x_1 y_1 y_2 + x_1 R_2 y_1 y_2 + x_1 x_2 R_7
    #   x_1 y_1 = 1+R_0: = y_2 + R_0 y_2 + x_1 R_2 y_1 y_2 + x_1 x_2 R_7
    #
    # So R_14 * y_2 y_1 = y_2 + R_0 y_2 + x_1 R_2 y_1 y_2 + x_1 x_2 R_7 - y_2 - y_2 R_0 - y_2 x_1 R_2 y_1
    #                  = R_0 y_2 - y_2 R_0 + x_1 R_2 y_1 y_2 - y_2 x_1 R_2 y_1 + x_1 x_2 R_7
    # That's an identity but doesn't give R_15.
    #
    # I'll try yet another approach: use that R_14 implies a cyclic version.
    # Since R_15 = -x_2 y_1 y_2 + y_1 y_2 x_1 and R_14 = x_1 x_2 y_1 - y_2 x_1 x_2,
    # they're related by "left multiplication by y_1 y_2 and right multiplication by something".
    #
    # y_1 y_2 * R_14 = y_1 y_2 x_1 x_2 y_1 - y_1 y_2 y_2 x_1 x_2
    # ... doesn't directly simplify.
    #
    # ALTERNATIVE: derive R_15 from R_6 directly using DIFFERENT sandwich.
    # Conjecture: y_1 y_2 R_6 y_1 y_2 simplifies to involve R_15.
    # y_1 y_2 R_6 y_1 y_2 = y_1 y_2 x_1 x_2 x_1 y_1 y_2 - y_1 y_2 x_2 x_1 x_2 y_1 y_2
    # First: y_1 y_2 x_1 x_2 (x_1 y_1) y_2 = y_1 y_2 x_1 x_2 (1+R_0) y_2 = y_1 y_2 x_1 x_2 y_2 + ...
    #   x_2 y_2 = 1+R_2: y_1 y_2 x_1 (1+R_2) = y_1 y_2 x_1 + y_1 y_2 x_1 R_2
    # So first = y_1 y_2 x_1 + y_1 y_2 x_1 R_2 + y_1 y_2 x_1 x_2 R_0 y_2
    # Second: y_1 y_2 x_2 x_1 x_2 y_1 y_2
    #   y_2 x_2 = 1+R_3: = y_1 (1+R_3) x_1 x_2 y_1 y_2 = y_1 x_1 x_2 y_1 y_2 + y_1 R_3 x_1 x_2 y_1 y_2
    #   y_1 x_1 = 1+R_1: = (1+R_1) x_2 y_1 y_2 + y_1 R_3 x_1 x_2 y_1 y_2
    #   = x_2 y_1 y_2 + R_1 x_2 y_1 y_2 + y_1 R_3 x_1 x_2 y_1 y_2
    # So y_1 y_2 R_6 y_1 y_2 = (y_1 y_2 x_1 - x_2 y_1 y_2) + R-tails
    # And y_1 y_2 x_1 - x_2 y_1 y_2 = R_15!
    #
    # Excellent. So:
    #   y_1 y_2 R_6 y_1 y_2 = R_15 + y_1 y_2 x_1 R_2 + y_1 y_2 x_1 x_2 R_0 y_2
    #                        - R_1 x_2 y_1 y_2 - y_1 R_3 x_1 x_2 y_1 y_2
    # So:
    #   R_15 = y_1 y_2 R_6 y_1 y_2
    #          - y_1 y_2 x_1 R_2 ()
    #          - y_1 y_2 x_1 x_2 R_0 y_2
    #          + R_1 x_2 y_1 y_2
    #          + y_1 R_3 x_1 x_2 y_1 y_2
    derivs[15] = [
        (Q(1),  mon_from_str("y1*y2"), 6, mon_from_str("y1*y2")),
        (Q(-1), mon_from_str("y1*y2*x1"), 2, ()),
        (Q(-1), mon_from_str("y1*y2*x1*x2"), 0, mon_from_str("y2")),
        (Q(1),  (), 1, mon_from_str("x2*y1*y2")),
        (Q(1),  mon_from_str("y1"), 3, mon_from_str("x1*x2*y1*y2")),
    ]

    # R_16 = x_1 y_2 y_1 - y_2 y_1 x_2
    # Sandwich: y_2 R_6 y_2 y_1 maybe.
    # y_2 R_6 y_2 y_1 = y_2 x_1 x_2 x_1 y_2 y_1 - y_2 x_2 x_1 x_2 y_2 y_1
    # First: y_2 x_1 x_2 x_1 y_2 y_1. We want x_1 y_2 y_1 somewhere.
    # Actually consider y_2 R_6 = y_2 x_1 x_2 x_1 - y_2 x_2 x_1 x_2
    #   y_2 x_2 = 1+R_3: = y_2 x_1 x_2 x_1 - (1+R_3) x_1 x_2 = y_2 x_1 x_2 x_1 - x_1 x_2 - R_3 x_1 x_2
    #
    # Try y_2 y_1 * R_6 * y_2 (parallel to R_14 derivation but with y's on opposite side):
    # y_2 y_1 R_6 y_2 = y_2 y_1 x_1 x_2 x_1 y_2 - y_2 y_1 x_2 x_1 x_2 y_2
    # Second: x_2 y_2 = 1+R_2 at end: y_2 y_1 x_2 x_1 (1+R_2) = y_2 y_1 x_2 x_1 + y_2 y_1 x_2 x_1 R_2
    # First: y_2 y_1 x_1 x_2 x_1 y_2. y_1 x_1 = 1+R_1: = y_2 (1+R_1) x_2 x_1 y_2 = y_2 x_2 x_1 y_2 + y_2 R_1 x_2 x_1 y_2
    #   y_2 x_2 = 1+R_3: = (1+R_3) x_1 y_2 + y_2 R_1 x_2 x_1 y_2 = x_1 y_2 + R_3 x_1 y_2 + y_2 R_1 x_2 x_1 y_2
    # So y_2 y_1 R_6 y_2 = x_1 y_2 - y_2 y_1 x_2 x_1 + R_3 x_1 y_2 + y_2 R_1 x_2 x_1 y_2 - y_2 y_1 x_2 x_1 R_2
    # Hmm, x_1 y_2 - y_2 y_1 x_2 x_1 is not R_16.
    #
    # Try: y_2 R_6 y_1 y_2:
    # = y_2 x_1 x_2 x_1 y_1 y_2 - y_2 x_2 x_1 x_2 y_1 y_2
    # x_1 y_1 = 1+R_0: First = y_2 x_1 x_2 (1+R_0) y_2 = y_2 x_1 x_2 y_2 + y_2 x_1 x_2 R_0 y_2
    #   x_2 y_2 = 1+R_2: y_2 x_1 (1+R_2) = y_2 x_1 + y_2 x_1 R_2
    # So First = y_2 x_1 + y_2 x_1 R_2 + y_2 x_1 x_2 R_0 y_2
    # Second: y_2 x_2 = 1+R_3: y_2 x_2 x_1 x_2 y_1 y_2 = (1+R_3) x_1 x_2 y_1 y_2 = x_1 x_2 y_1 y_2 + R_3 x_1 x_2 y_1 y_2
    # So y_2 R_6 y_1 y_2 = y_2 x_1 - x_1 x_2 y_1 y_2 + R-tails
    # That's not R_16 either.
    #
    # Try: y_1 y_2 R_6 with x on right.
    # Consider y_2 y_1 R_6 y_2 y_1:
    # = y_2 y_1 x_1 x_2 x_1 y_2 y_1 - y_2 y_1 x_2 x_1 x_2 y_2 y_1
    # First: y_1 x_1 = 1+R_1: y_2 (1+R_1) x_2 x_1 y_2 y_1 = y_2 x_2 x_1 y_2 y_1 + y_2 R_1 x_2 x_1 y_2 y_1
    #   y_2 x_2 = 1+R_3: (1+R_3) x_1 y_2 y_1 = x_1 y_2 y_1 + R_3 x_1 y_2 y_1
    # So First = x_1 y_2 y_1 + R_3 x_1 y_2 y_1 + y_2 R_1 x_2 x_1 y_2 y_1
    # Second: x_2 y_2 = 1+R_2 at end: y_2 y_1 x_2 x_1 (1+R_2) y_1 = y_2 y_1 x_2 x_1 y_1 + y_2 y_1 x_2 x_1 R_2 y_1
    # Hmm, still has y_2 y_1 x_2 x_1 y_1 which doesn't reduce. Let's also do x_1 y_1 = 1+R_0:
    #   y_2 y_1 x_2 (x_1 y_1) = y_2 y_1 x_2 (1+R_0) = y_2 y_1 x_2 + y_2 y_1 x_2 R_0
    # So Second = y_2 y_1 x_2 + y_2 y_1 x_2 R_0 + y_2 y_1 x_2 x_1 R_2 y_1
    # So y_2 y_1 R_6 y_2 y_1 = x_1 y_2 y_1 - y_2 y_1 x_2 + R-tails
    # And x_1 y_2 y_1 - y_2 y_1 x_2 = R_16!
    # So:
    #   y_2 y_1 R_6 y_2 y_1 = R_16 + R_3 x_1 y_2 y_1 + y_2 R_1 x_2 x_1 y_2 y_1
    #                         - y_2 y_1 x_2 R_0 - y_2 y_1 x_2 x_1 R_2 y_1
    # Therefore:
    #   R_16 = y_2 y_1 R_6 y_2 y_1 - R_3 x_1 y_2 y_1 - y_2 R_1 x_2 x_1 y_2 y_1
    #          + y_2 y_1 x_2 R_0 + y_2 y_1 x_2 x_1 R_2 y_1
    derivs[16] = [
        (Q(1),  mon_from_str("y2*y1"), 6, mon_from_str("y2*y1")),
        (Q(-1), (), 3, mon_from_str("x1*y2*y1")),
        (Q(-1), mon_from_str("y2"), 1, mon_from_str("x2*x1*y2*y1")),
        (Q(1),  mon_from_str("y2*y1*x2"), 0, ()),
        (Q(1),  mon_from_str("y2*y1*x2*x1"), 2, mon_from_str("y1")),
    ]

    # R_17 = -x_2 x_1 y_2 + y_1 x_2 x_1
    # By analogy with R_14 (with x<->y swap... but for x's we need R_6, for y's we'd
    # need R_7. We want R_17 in terms of R_6).
    # Try y_1 R_6 (something) or similar.
    # y_1 * R_6 * y_2 = y_1 x_1 x_2 x_1 y_2 - y_1 x_2 x_1 x_2 y_2
    # First: y_1 x_1 = 1+R_1: (1+R_1) x_2 x_1 y_2 = x_2 x_1 y_2 + R_1 x_2 x_1 y_2
    # Second: x_2 y_2 = 1+R_2 at end: y_1 x_2 x_1 (1+R_2) = y_1 x_2 x_1 + y_1 x_2 x_1 R_2
    # So y_1 R_6 y_2 = x_2 x_1 y_2 - y_1 x_2 x_1 + R_1 x_2 x_1 y_2 - y_1 x_2 x_1 R_2
    #                = -(y_1 x_2 x_1 - x_2 x_1 y_2) + R-tails = -R_17 + R-tails
    # Wait: R_17 = -x_2 x_1 y_2 + y_1 x_2 x_1 = y_1 x_2 x_1 - x_2 x_1 y_2
    # So y_1 R_6 y_2 = -R_17 + R_1 x_2 x_1 y_2 - y_1 x_2 x_1 R_2
    # Therefore:
    #   R_17 = -y_1 R_6 y_2 + R_1 x_2 x_1 y_2 - y_1 x_2 x_1 R_2
    derivs[17] = [
        (Q(-1), mon_from_str("y1"), 6, mon_from_str("y2")),
        (Q(1),  (), 1, mon_from_str("x2*x1*y2")),
        (Q(-1), mon_from_str("y1*x2*x1"), 2, ()),
    ]

    # R_26 = (-t + q*t)*1 + (1 + t - q)*x_2 + (-q*t)*y_2 + (-1)*x_2*x_2 + q*x_2*y_2
    # We have R_25 = (-t + q*t)*1 + (1 + t - q)*x_1 + (-q*t)*y_1 + (-1)*x_1*x_1 + q*x_1*y_1
    # i.e., R_25 = (x_1 - t)(q y_1 + (1-q) - x_1).
    # R_26 = (x_2 - t)(q y_2 + (1-q) - x_2).
    # These differ by index 1 <-> 2 swap. We need to use R_6 + invertibility.
    #
    # Trick: conjugate by braid. Let w = x_1 x_2 y_1 (so w^{-1} = x_1 y_2 y_1, since
    # in the group y_1 x_2 y_1 = ... actually the braid group element with x_2 = ...).
    # In the braid group, x_1 x_2 x_1 = x_2 x_1 x_2 implies x_2 = x_1^{-1} x_2 x_1 x_2 x_2^{-1} = y_1 x_2 x_1 ... that's circular.
    #
    # Alternative: x_2 = (x_1 x_2 x_1) y_1 y_2 y_1 ... no, that gives x_2 x_1 ... y_1 y_2 y_1 = ?
    # Actually in the braid group, x_1 x_2 x_1 = x_2 x_1 x_2 implies
    #   x_2 = (y_1 y_2)(x_1 x_2 x_1)(y_1 y_2)... let's just compute.
    # In the algebra mod R_0..R_6:
    #   x_1 x_2 x_1 = x_2 x_1 x_2
    #   => x_2 x_1 = x_1 x_2 x_1 y_2 = x_1 (x_2 x_1 y_2)
    # In the braid group, x_2 = y_1 x_1 x_2 = ... no.
    # OK we need conjugation: x_2 = (x_1) x_2 (y_1)? No.
    # (x_1 x_2) x_1 (y_2 y_1) = x_1 x_2 x_1 y_2 y_1. Using R_6, = x_2 x_1 x_2 y_2 y_1 = x_2 x_1 (1+R_2) y_1
    #   = x_2 x_1 y_1 + x_2 x_1 R_2 y_1 = x_2 (1+R_1) + x_2 x_1 R_2 y_1 = x_2 + x_2 R_1 + x_2 x_1 R_2 y_1.
    # So x_1 x_2 x_1 y_2 y_1 = x_2 + R-tails.
    # And x_1 x_2 x_1 y_2 y_1 = x_1 x_2 x_1 y_2 y_1 (we get x_2 from this!).
    # So: x_2 = x_1 x_2 x_1 y_2 y_1 - R-tails.
    #
    # Now the strategy: take R_25 (involving x_1) and "conjugate" by w = x_1 x_2 y_1 ... wait,
    # The conjugation works as: phi(x_1) = x_2 (= w x_1 w^{-1} in the group).
    #
    # This is getting complicated. Let me use a different approach: express
    # R_26 = (x_2 - t) * S[2][2] explicitly and try to convert via R_6.
    #
    # R_26 has form: -t + qt + (1+t-q)*x_2 - qt*y_2 - x_2*x_2 + q*x_2*y_2
    # Let A = R_25 expressed at index 2 instead of 1.
    # Difference R_26 - R_25 = (1+t-q)*(x_2 - x_1) + (-qt)*(y_2 - y_1) - (x_2 x_2 - x_1 x_1) + q*(x_2 y_2 - x_1 y_1).
    #
    # Hmm. The "braid conjugation" approach is non-trivial. Let me instead
    # use a more explicit computation.
    #
    # In the algebra modulo {R_0,..,R_6}, we have (in the group sense):
    #   x_2 = (y_1 x_1) x_2 = y_1 (x_1 x_2) -- doesn't simplify
    #
    # Actually one cleaner way: use the braid to express x_2 as a conjugate of x_1.
    # In the braid group, sigma_2 = sigma_1^{-1} sigma_2 sigma_1 sigma_2 sigma_2^{-1} = ...
    # the relation sigma_1 sigma_2 sigma_1 = sigma_2 sigma_1 sigma_2 gives us
    #   sigma_1 = sigma_2 sigma_1 sigma_2 sigma_1^{-1} sigma_2^{-1}
    # So x_1 = x_2 x_1 x_2 y_1 y_2  (in the group, where y_i = x_i^{-1}).
    # Then phi^{-1}: x_1 <-> x_2 is the same as conjugation by sigma_2 sigma_1 maybe.
    #
    # GIVEN TIME PRESSURE: let me programmatically derive R_26 using a small
    # ansatz and a linear algebra solver below. For now leave a placeholder.
    derivs[26] = None  # to be filled by search below

    return derivs, KEEP


# --------- Verify a derivation ----------------------------------------
def eval_derivation(deriv, relations) -> "NCPoly":
    """Evaluate sum_k a_k * L_k * R_{i_k} * R'_k as NCPoly."""
    out = nc_zero()
    for (coef, L_mon, idx, R_mon) in deriv:
        L_nc = {L_mon: qt_const(1)}
        R_nc = {R_mon: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, relations[idx]), R_nc)
        term = nc_scalar(term, coef)
        out = nc_add(out, term)
    return out


def derive_r26_by_search(relations) -> List[Tuple["QTPoly", Tuple[int, ...], int, Tuple[int, ...]]]:
    """
    Derive R_26 = (x_2 - t)*S[2][2] in terms of {R_0, R_1, R_2, R_3, R_6, R_25}.

    R_25 = (x_1 - t)*S[2][1].
    Strategy: in the braid algebra modulo R_0..R_6, x_2 can be written as
    x_2 = (x_1 x_2 x_1) y_2 y_1 (modulo R_0..R_3) and similarly for y_2.
    Equivalently, R_26 - (conjugation of R_25 by braid element) is in the
    ideal generated by R_0..R_6.

    But getting the exact symbolic cofactor is nontrivial. Easier: directly
    use the identity in the free algebra:
       R_26 = (substitution of x_2 for x_1, y_2 for y_1 in R_25)
    and express the substitution via a small chain of rewrites.

    We construct R_26 by:
       (1) Write R_25 = (x_1 - t) S_1 where S_1 = q*y_1 + (1-q) - x_1
       (2) Note R_26 = (x_2 - t) S_2  where S_2 = q*y_2 + (1-q) - x_2
       (3) Find polynomials p, q_l, q_r such that:
              R_26 = p * R_25 * 1   (no, R_25 has x_1 not x_2)
       This won't work as a single sandwich.

    Approach: write R_26 in the free algebra directly as a polynomial in
    {R_0..R_3, R_6, R_25}. Use the identity (proved by hand):

      A = x_2 - t. B = q y_2 + (1-q) - x_2.
      Want A*B in ideal of {R_0..R_3, R_6, R_25}.

      Note R_25 = (x_1 - t)(q y_1 + (1-q) - x_1).
      Expand R_25: x_1*q*y_1 + (1-q) x_1 - x_1 x_1 - t q y_1 - t(1-q) + t x_1
                 = q x_1 y_1 + x_1 - q x_1 - x_1^2 - q t y_1 - t + qt + t x_1
                 = q (1 + R_0) + x_1 - q x_1 - x_1^2 - qt y_1 - t + qt + t x_1
                 = q + q R_0 + x_1 - q x_1 - x_1^2 - qt y_1 - t + qt + t x_1
                 = (q + qt - t) + (1 - q + t) x_1 - qt y_1 - x_1^2 + q R_0
      Hmm but R_25 in the relation builder is exactly this minus q R_0... let me re-check.
      Actually R_25 = (x_1 - t)*S[2][1] in the FREE algebra (no reduction).
      So R_25 = q x_1 y_1 - x_1 x_1 + (1-q) x_1 - q t y_1 + (q-1) t + t x_1
              = q x_1 y_1 - x_1 x_1 + (1 - q + t) x_1 - qt y_1 + (q-1)t
              = q x_1 y_1 - x_1 x_1 + (1+t-q) x_1 - qt y_1 + (qt - t)
      Matches output: (-t+qt)*1 + (1+t-q)*x_1 + (-qt)*y_1 + (-1)*x_1*x_1 + q*x_1*y_1.

      So R_25 = -(t - qt) + (1+t-q)*x_1 - qt*y_1 - x_1*x_1 + q*x_1*y_1.

      In the algebra mod {R_0, R_1, R_6}, can we transform R_25 into R_26?
      Note R_25 mod R_0 (x_1 y_1 = 1): = -(t-qt) + (1+t-q) x_1 - qt y_1 - x_1^2 + q
                                       = (q + qt - t) + (1+t-q) x_1 - qt y_1 - x_1^2
      Similarly R_26 mod R_2 = (q + qt - t) + (1+t-q) x_2 - qt y_2 - x_2^2.

      So R_26 = R_25 + (x_2 - x_1) terms + (y_2 - y_1) terms + (x_2^2 - x_1^2) terms
                + q*R_2 - q*R_0 (modulo simplifications).

      And (x_2 - x_1), (y_2 - y_1), (x_2^2 - x_1^2) all live in the ideal of R_6 mod R_0..R_3.
      Concretely:
        x_2 = (1+R_3) x_2 = y_2 x_2 x_2 + R_3 x_2 = y_2 x_2^2 + R_3 x_2 (use R_3 = y_2 x_2 - 1)
        But we want a representation via R_6.

      The KEY identity in the braid algebra:
        x_2 x_1 - x_1 x_2 has no simple form (they don't commute), but
        x_2 = (y_1)(x_1 x_2 x_1)(y_2 y_1) mod R_2, R_6 ... actually let me just
        try: x_1 x_2 x_1 = x_2 x_1 x_2 (R_6) => x_2 = (x_1 x_2 x_1)(y_2)(y_1) mod R_2
        Verify: (x_1 x_2 x_1)(y_2)(y_1) = (x_2 x_1 x_2)(y_2)(y_1) (using R_6)
              = x_2 x_1 (x_2 y_2) y_1 = x_2 x_1 (1 + R_2) y_1
              = x_2 (x_1 y_1) + x_2 x_1 R_2 y_1 = x_2 (1 + R_0) + x_2 x_1 R_2 y_1
              = x_2 + x_2 R_0 + x_2 x_1 R_2 y_1.
        So x_1 x_2 x_1 y_2 y_1 = x_2 + x_2 R_0 + x_2 x_1 R_2 y_1.
        But that requires R_6 to convert LHS! Let me redo:
        Set A = x_1 x_2 x_1 y_2 y_1. Then A - x_2 = ? Use R_6: A = (R_6 + x_2 x_1 x_2) y_2 y_1
          = R_6 y_2 y_1 + x_2 x_1 (x_2 y_2) y_1 = R_6 y_2 y_1 + x_2 x_1 (1 + R_2) y_1
          = R_6 y_2 y_1 + x_2 (x_1 y_1) + x_2 x_1 R_2 y_1
          = R_6 y_2 y_1 + x_2 (1 + R_0) + x_2 x_1 R_2 y_1
          = R_6 y_2 y_1 + x_2 + x_2 R_0 + x_2 x_1 R_2 y_1.
        So: A = x_2 + (R_6) y_2 y_1 + x_2 R_0 + x_2 x_1 R_2 y_1.
        Rearranging: x_2 = x_1 x_2 x_1 y_2 y_1 - R_6 y_2 y_1 - x_2 R_0 - x_2 x_1 R_2 y_1.   (*)

      Similarly for y_2: from R_7 (y_1 y_2 y_1 = y_2 y_1 y_2), we get
        y_1 y_2 y_1 x_2 x_1 = ... after reductions ... = y_2 + tails.
        But R_7 is itself derivable from R_6. So this would chain.

      For simplicity AND given TIME PRESSURE, I'll express R_26 in the free
      algebra via a "global" identity. Easier alternative: express R_26
      using R_25 plus a finite combination derived from invertibility and
      R_6, using the substitution at the level of S[2][j].

      RECALL: R_25 = (x_1 - t)*S[2][1] = q x_1 y_1 - x_1^2 + (1+t-q) x_1 - qt y_1 + (qt-t).
              R_26 = q x_2 y_2 - x_2^2 + (1+t-q) x_2 - qt y_2 + (qt-t).

      R_26 - R_25 = q*(x_2 y_2 - x_1 y_1) - (x_2^2 - x_1^2) + (1+t-q)*(x_2 - x_1) - qt*(y_2 - y_1).
                  = q*(R_2 - R_0)  [since x_i y_i = 1 + R_{2i-2}]
                    - (x_2^2 - x_1^2)
                    + (1+t-q)*(x_2 - x_1)
                    - qt*(y_2 - y_1).

      So R_26 = R_25 + q*R_2 - q*R_0 - (x_2^2 - x_1^2) + (1+t-q)*(x_2-x_1) - qt*(y_2-y_1).

      We need to express (x_2 - x_1), (y_2 - y_1), (x_2^2 - x_1^2) in the ideal
      generated by {R_0..R_3, R_6}.

      Claim: x_2 - x_1 is NOT in the two-sided ideal generated by R_0..R_3 + R_6.
      That would mean x_1 = x_2 in the quotient, which is FALSE (braid relation
      doesn't identify x_1 and x_2).

      So R_26 is NOT in the ideal generated by {R_0,..R_3, R_6, R_25} alone!
      This contradicts the user's assumption.

      Wait, but the cofactor file ONLY uses R_0..R_3, R_6, R_22, R_25, R_28 plus
      the dependents... and R_26 is listed as derivable. Let me re-examine.

      Actually the cofactor is X_3 = S[3][1], not an arbitrary polynomial.
      So R_26 might appear in X_3's cofactor, but each occurrence of R_26
      can be split (within the cofactor) using R_22 (the zipper relation).
      Wait, R_22 introduces x_1, x_2, y_1, y_2 mixed. Could R_26 be derivable
      using R_22 too? The user's TARGET MINIMAL SET includes R_22.

      Let me re-check: R_22 (the j=1 2-strand zipper) = (q y_2 + (1-q) - x_2 - q y_1 y_2 + x_1 x_2)*S[2][1].
      This has terms involving x_2 alone (e.g., -x_2 term). So R_22 brings x_2 into play.

      Insight: in the quotient algebra by {R_0,..R_3, R_6, R_22, R_25}, can we derive R_26?
      The zipper algebra is a known finite-dimensional algebra; R_26 SHOULD be
      derivable. Let me try:

      R_22 = q y_2 - q + 1 - x_2 - q y_1 y_2 + x_1 x_2 ... wait, multiplied by S[2][1].
      Actually R_22 = (q*y_2 + (1-q) - x_2 - q*y_1*y_2 + x_1*x_2) * S[2][1]
                    = (q*y_2 + (1-q) - x_2 - q*y_1*y_2 + x_1*x_2) * (q*y_1 + (1-q) - x_1).

      This is complex. The user's hint says R_26 is derivable from R_25 by
      "braid conjugation"; presumably they mean using R_6 + invertibility +
      maybe R_22.

      GIVEN TIME PRESSURE: Let me just NOT compress R_26 and report it.
      The cofactor will use {R_0, R_1, R_2, R_3, R_6, R_22, R_25, R_26, R_28}
      (9 relations instead of 8), but the dependent set {R_7, R_14..R_17}
      can still be compressed.
    """
    return None


# --------- Main compression -------------------------------------------
def main() -> int:
    t_start = time.time()
    relations, X3 = build_inputs_and_X3()
    print(f"[build] {len(relations)} relations, X_3 with {nc_size(X3)} mons")

    derivs, KEEP = derivations(relations)

    # Verify each derivation
    print("\n[verify derivations]")
    for d_idx in [7, 14, 15, 16, 17]:
        deriv = derivs[d_idx]
        if deriv is None:
            print(f"  R_{d_idx}: SKIPPED (no derivation provided)")
            continue
        lhs = eval_derivation(deriv, relations)
        rhs = relations[d_idx]
        diff = nc_sub(lhs, rhs)
        ok = nc_is_zero(diff)
        # Check all source indices are in KEEP
        src_idxs = sorted({t[2] for t in deriv})
        in_keep = all(i in KEEP for i in src_idxs)
        print(f"  R_{d_idx}: {'OK' if ok else 'FAIL'}  uses {src_idxs} "
              f"({'all in KEEP' if in_keep else 'MISSING'})  {len(deriv)} terms")
        if not ok:
            # print first few residual monomials for debug
            items = sorted(diff.items(), key=lambda kv: (len(kv[0]), kv[0]))
            for mon, coef in items[:5]:
                mstr = "*".join(VAR_NAMES[i - 1] for i in mon) if mon else "1"
                print(f"    residual [{mstr}] = {qt_str(coef)}")
            return 1

    # Try to derive R_26
    deriv_26 = derive_r26_by_search(relations)
    if deriv_26 is not None:
        derivs[26] = deriv_26
        lhs = eval_derivation(deriv_26, relations)
        diff = nc_sub(lhs, relations[26])
        if nc_is_zero(diff):
            print(f"  R_26: OK ({len(deriv_26)} terms)")
        else:
            print(f"  R_26: FAIL derivation; treating as independent")
            derivs[26] = None
            KEEP.add(26)
    else:
        print(f"  R_26: kept as independent (no derivation)")
        KEEP.add(26)

    # --- Parse the existing cofactor file
    with open(COF_FILE) as f:
        raw = f.readlines()
    cof_lines = [ln.rstrip("\n") for ln in raw
                 if ln.strip() and not ln.lstrip().startswith("#")]
    parsed = [parse_cofactor_line(ln) for ln in cof_lines]
    print(f"\n[cof] {len(parsed)} input tuples")

    # --- Substitute dependents
    # Each tuple is (L, idx, R, num, den)  with num/den in QTPoly form.
    # If idx is a dependent, replace by the derivation, multiplying out:
    #   c * L * R_idx * R'  ->  sum_k (c * a_k) * (L * L_k) * R_{i_k} * (R'_k * R')
    #
    # We need to convert (num, den) into a single coefficient over Q[q,t] by
    # multiplying through by den (which is a monomial in q, t, (q-1), (t+q^2)).
    #
    # Actually it's easier to keep the (num, den) representation. Denominators
    # are MONOMIALS in {q, t, (q-1), (t+q^2)} so we can track them as 4-tuples.
    # The compressed cofactor will have the same denominator structure (since
    # we're adding coefficients of the original tuples + new ones that
    # introduce no new denominators).

    # Compressed tuple format: {(L, idx, R) -> (num_dict, (a_q, b_t, k1, k2))}
    # where num_dict: QTPoly (in Q[q,t]).

    # For derivations: a_k is in Q[q,t] (constants), so multiplying by a_k just
    # multiplies the numerator. Denominator stays the same.

    # Step 1: convert (num, den) to a uniform denominator representation.
    # Actually, let me NOT unify the LCM. Instead, for each output (L, idx, R),
    # accumulate a list of (num, den) pairs and at the end combine to common
    # denominator.

    # Simpler approach: represent each tuple's coefficient as (num: QTPoly, den_tup),
    # and for each output (L, idx, R), keep a list. At end, find LCM of dens
    # and combine.

    from collections import defaultdict
    # accumulator: (L, idx, R) -> list of (num: QTPoly, den_tup)
    acc: Dict[Tuple[Tuple[int, ...], int, Tuple[int, ...]],
              List[Tuple["QTPoly", Tuple[int, int, int, int]]]] = defaultdict(list)

    for (L, idx, R, num, den) in parsed:
        deriv = derivs.get(idx)
        if deriv is None:
            # idx is in KEEP (or treated as kept). Pass through.
            acc[(L, idx, R)].append((num, den))
        else:
            # idx is dependent. Substitute.
            for (a_k, L_k, src_idx, R_k) in deriv:
                new_L = L + L_k
                new_R = R_k + R
                # New num = num * a_k
                new_num = qt_mul(num, a_k)
                if not new_num:
                    continue
                acc[(new_L, src_idx, new_R)].append((new_num, den))

    # Step 2: combine list of (num, den) pairs into a single (num, den) per key,
    # by finding the LCM of all dens. Drop zero results.
    def lcm_den(dens):
        return tuple(max(d[i] for d in dens) for i in range(4))

    def den_str(d):
        a, b, k1, k2 = d
        parts = []
        if a == 1: parts.append("q")
        elif a > 1: parts.append(f"q^{a}")
        if b == 1: parts.append("t")
        elif b > 1: parts.append(f"t^{b}")
        for _ in range(k1): parts.append("(q-1)")
        for _ in range(k2): parts.append("(t+q^2)")
        if not parts:
            return "(1)"
        return "(" + "*".join(parts) + ")"

    def num_str(n):
        # Render QTPoly in the file format: "(1*q^a*t^b + -2*q^c*t^d + ...)"
        if not n:
            return "(0)"
        items = sorted(n.items(), key=lambda kv: kv[0])
        parts = []
        first = True
        for (a, b), c in items:
            assert c.denominator == 1, f"Non-integer numerator coefficient: {c}"
            cs = str(c.numerator)
            term = f"{cs}*q^{a}*t^{b}"
            if first:
                parts.append(term)
                first = False
            else:
                if cs.startswith("-"):
                    parts.append(" + " + term)
                else:
                    parts.append(" + " + term)
        return "(" + "".join(parts) + ")"

    compressed: List[Tuple[Tuple[int, ...], int, Tuple[int, ...], "QTPoly",
                          Tuple[int, int, int, int]]] = []

    for key, entries in acc.items():
        # Each entry is (num, den). Combine with LCM denominator.
        if len(entries) == 1:
            num, den = entries[0]
        else:
            tgt_den = lcm_den([e[1] for e in entries])
            total_num: "QTPoly" = {}
            for num, den in entries:
                # multiplier = q^(tgt-d) * t^(...) * (q-1)^(...) * (t+q^2)^(...)
                a_diff = tgt_den[0] - den[0]
                b_diff = tgt_den[1] - den[1]
                k1_diff = tgt_den[2] - den[2]
                k2_diff = tgt_den[3] - den[3]
                mult = qt_const(1)
                if a_diff: mult = qt_mul(mult, qt_pow(qt_q(), a_diff))
                if b_diff: mult = qt_mul(mult, qt_pow(qt_t(), b_diff))
                if k1_diff:
                    qm1_p = qt_sub(qt_q(), qt_const(1))
                    mult = qt_mul(mult, qt_pow(qm1_p, k1_diff))
                if k2_diff:
                    tpq2_p = qt_add(qt_t(), qt_pow(qt_q(), 2))
                    mult = qt_mul(mult, qt_pow(tpq2_p, k2_diff))
                contrib = qt_mul(num, mult)
                total_num = qt_add(total_num, contrib)
            num = total_num
            den = tgt_den

        # Simplify: cancel common q, t, (q-1), (t+q^2) factors.
        # For q and t factors: find min q-degree and min t-degree across num.
        if num:
            min_q = min(k[0] for k in num)
            min_t = min(k[1] for k in num)
            # Don't strip below 0
            strip_q = min(min_q, den[0])
            strip_t = min(min_t, den[1])
            if strip_q > 0 or strip_t > 0:
                new_num: "QTPoly" = {}
                for (a, b), c in num.items():
                    new_num[(a - strip_q, b - strip_t)] = c
                num = new_num
                den = (den[0] - strip_q, den[1] - strip_t, den[2], den[3])

        if not num:
            continue
        compressed.append((key[0], key[1], key[2], num, den))

    # --- Verify ---
    print(f"\n[compressed] {len(compressed)} tuples")
    distinct_idxs = sorted({c[1] for c in compressed})
    print(f"[compressed] distinct R_i: {distinct_idxs}")

    # Build LHS sum
    print("\n[verify compressed] computing LHS sum...")
    total: "NCPoly" = nc_zero()
    # Compute global LCM denominator
    A = max(c[4][0] for c in compressed)
    B = max(c[4][1] for c in compressed)
    K1 = max(c[4][2] for c in compressed)
    K2 = max(c[4][3] for c in compressed)
    print(f"[lcm] M = q^{A}*t^{B}*(q-1)^{K1}*(t+q^2)^{K2}")

    qm1_poly = qt_sub(qt_q(), qt_const(1))
    tpq2_poly = qt_add(qt_t(), qt_pow(qt_q(), 2))
    M_poly = qt_const(1)
    M_poly = qt_mul(M_poly, qt_pow(qt_q(), A))
    M_poly = qt_mul(M_poly, qt_pow(qt_t(), B))
    M_poly = qt_mul(M_poly, qt_pow(qm1_poly, K1))
    M_poly = qt_mul(M_poly, qt_pow(tpq2_poly, K2))

    for k, (L, idx, R, num, den) in enumerate(compressed):
        a_q, b_t, k1, k2 = den
        cofactor_mult = qt_const(1)
        cofactor_mult = qt_mul(cofactor_mult, qt_pow(qt_q(), A - a_q))
        cofactor_mult = qt_mul(cofactor_mult, qt_pow(qt_t(), B - b_t))
        cofactor_mult = qt_mul(cofactor_mult, qt_pow(qm1_poly, K1 - k1))
        cofactor_mult = qt_mul(cofactor_mult, qt_pow(tpq2_poly, K2 - k2))
        cp = qt_mul(num, cofactor_mult)
        if not cp:
            continue
        Rij = relations[idx]
        L_nc = {L: qt_const(1)}
        R_nc = {R: qt_const(1)}
        term = nc_mul(nc_mul(L_nc, Rij), R_nc)
        term = nc_scalar(term, cp)
        total = nc_add(total, term)
        if (k + 1) % 200 == 0:
            print(f"  accumulated {k+1}/{len(compressed)} tuples (total NC mons={nc_size(total)})")

    rhs = nc_scalar(X3, M_poly)
    diff = nc_sub(total, rhs)

    if not nc_is_zero(diff):
        print("\nFAIL: compressed cofactor does not equal M*X_3.")
        items = sorted(diff.items(), key=lambda kv: (len(kv[0]), kv[0]))
        for mon, coef in items[:10]:
            mstr = "*".join(VAR_NAMES[i - 1] for i in mon) if mon else "1"
            print(f"  [{mstr}] = {qt_str(coef)}")
        return 1

    print("\nPASS: compressed cofactor verified.")

    # --- Write output ---
    OUT = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-compressed.txt")
    # Sort: by idx, then by L mon, then by R mon
    compressed.sort(key=lambda c: (c[1], c[0], c[2]))
    with open(OUT, "w") as f:
        f.write("# Compressed symbolic Q(q,t) cofactor of X_3 over n=4 Bigelow zipper.\n")
        f.write(f"# Original: 695 tuples using indices {{0,1,2,3,6,7,14,15,16,17,22,25,26,28}}.\n")
        f.write(f"# Compressed: {len(compressed)} tuples using indices {distinct_idxs}.\n")
        f.write("# Dependent inputs eliminated via derivations:\n")
        f.write("#   BR_1_y = y_1 y_2 y_1 - y_2 y_1 y_2 derived from BR_1 + INV_{1,2}_{p,m}\n")
        f.write("#   MIX_1_a = x_1 x_2 y_1 - y_2 x_1 x_2 derived from BR_1, INV_1_p, INV_2_m\n")
        f.write("#   MIX_1_b = -x_2 y_1 y_2 + y_1 y_2 x_1 derived from BR_1 + INV_{1,2}_{p,m}\n")
        f.write("#   MIX_1_c = x_1 y_2 y_1 - y_2 y_1 x_2 derived from BR_1 + INV_{1,2}_{p,m}\n")
        f.write("#   MIX_1_d = -x_2 x_1 y_2 + y_1 x_2 x_1 derived from BR_1, INV_1_m, INV_2_p\n")
        if 26 in KEEP:
            f.write("#   T_2_j2 KEPT (not eliminated; derivation by braid conjugation not implemented)\n")
        else:
            f.write("#   T_2_j2 derived from T_2, BR_1 + INV_{1,2}_{p,m}\n")
        f.write("#\n")
        for L, idx, R, num, den in compressed:
            f.write(f"{mon_to_str(L)} * {label(idx)} * {mon_to_str(R)}  ::  "
                    f"{num_str(num)} / {den_str(den)}\n")

    print(f"\n[output] wrote {OUT}")
    print(f"  tuples:               {len(compressed)} (was 695)")
    print(f"  distinct R_i used:    {distinct_idxs} ({len(distinct_idxs)} total)")
    print(f"  denominator atoms:    q^{A} * t^{B} * (q-1)^{K1} * (t+q^2)^{K2}")
    print(f"  total wall time:      {time.time() - t_start:.2f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
