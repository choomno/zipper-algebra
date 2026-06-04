"""Canonical Z/T naming for the 31 defining relations of Z_4.

The integer indices R_0..R_30 used internally by f4ncgb / extract_cofactor.py
are mapped to readable names of the form
    Z_k        level-k zipper at anchor j=1
    Z_k_j<n>   level-k zipper at anchor j=n  (n >= 2)
    T_k        level-k untwisting at anchor j=1
    T_k_j<n>   level-k untwisting at anchor j=n  (n >= 2)
    INV_i_p    invertibility x_i y_i - 1
    INV_i_m    invertibility y_i x_i - 1
    BR_i       x-braid for adjacent generators sigma_i, sigma_{i+1}
    BR_i_y     y-braid for adjacent generators sigma_i^{-1}, sigma_{i+1}^{-1}
    FC_<sgn>   far-commutativity, sigma_1 vs sigma_3 (sgn in {xx,yy,yx,xy})
    MIX_i_<x>  length-3 mixed-sign braid identity, family i (x in {a,b,c,d})

This is the convention used in the cofactor .txt reports.
"""
from __future__ import annotations
from typing import Dict

RELATION_LABEL: Dict[int, str] = {
    # Invertibility, R_0..R_5
    0: "INV_1_p",  # x1*y1 - 1
    1: "INV_1_m",  # y1*x1 - 1
    2: "INV_2_p",  # x2*y2 - 1
    3: "INV_2_m",  # y2*x2 - 1
    4: "INV_3_p",  # x3*y3 - 1
    5: "INV_3_m",  # y3*x3 - 1
    # Braid + far-commutativity, R_6..R_13
    6: "BR_1",       # x1*x2*x1 - x2*x1*x2
    7: "BR_1_y",     # y1*y2*y1 - y2*y1*y2
    8: "FC_xx",      # x1*x3 - x3*x1
    9: "FC_yy",      # y1*y3 - y3*y1
    10: "FC_yx",     # y1*x3 - x3*y1
    11: "FC_xy",     # x1*y3 - y3*x1
    12: "BR_2",      # x2*x3*x2 - x3*x2*x3
    13: "BR_2_y",    # y2*y3*y2 - y3*y2*y3
    # Length-3 mixed sign, R_14..R_21
    14: "MIX_1_a",   # x1*x2*y1 - y2*x1*x2
    15: "MIX_1_b",   # y1*y2*x1 - x2*y1*y2
    16: "MIX_1_c",   # x1*y2*y1 - y2*y1*x2
    17: "MIX_1_d",   # y1*x2*x1 - x2*x1*y2
    18: "MIX_2_a",   # x2*x3*y2 - y3*x2*x3
    19: "MIX_2_b",   # y2*y3*x2 - x3*y2*y3
    20: "MIX_2_c",   # x2*y3*y2 - y3*y2*x3
    21: "MIX_2_d",   # y2*x3*x2 - x3*x2*y3
    # Zipper relations
    22: "Z_2",       # level-2 zipper at j=1
    23: "Z_2_j2",    # level-2 zipper at j=2
    24: "Z_3",       # level-3 zipper at j=1
    # Untwisting relations
    25: "T_2",       # level-2 untwisting at j=1: (x_1 - t) * S[2][1]
    26: "T_2_j2",    # level-2 untwisting at j=2
    27: "T_2_j3",    # level-2 untwisting at j=3
    28: "T_3",       # level-3 untwisting at j=1
    29: "T_3_j2",    # level-3 untwisting at j=2
    30: "T_4",       # level-4 untwisting at j=1
}


def label(idx: int) -> str:
    """Return the readable label for relation index idx."""
    return RELATION_LABEL[idx]


# Reverse mapping
INDEX_BY_LABEL: Dict[str, int] = {v: k for k, v in RELATION_LABEL.items()}


def parse_label(s: str) -> int:
    """Parse a label back to an integer index.

    Accepts either a Z/T-style label (e.g. 'Z_2', 'T_3_j2', 'INV_1_p')
    or the legacy 'R_<idx>' form (so old cofactor files still parse).
    """
    s = s.strip()
    if s.startswith("R_") and s[2:].isdigit():
        return int(s[2:])
    if s in INDEX_BY_LABEL:
        return INDEX_BY_LABEL[s]
    raise ValueError(f"Unknown relation label: {s!r}")


def header_legend() -> str:
    """Returns a multi-line header explaining each label."""
    lines = [
        "# Relation labels (Z = zipper, T = untwisting; subscript = level, j = anchor):",
        "#   Z_2          level-2 zipper      (q*y_{j+1} + (1-q) - x_{j+1}",
        "#                                     - q*y_j*y_{j+1} + x_j*x_{j+1}) * S[2][j] | j=1",
        "#   Z_2_j2       same with j=2",
        "#   Z_3          level-3 zipper (j=1, n>=4): contains x_3, y_3",
        "#   T_2          level-2 untwist (x_1 - t) * S[2][1]",
        "#   T_2_j2,_j3   level-2 untwist at j=2, j=3",
        "#   T_3          level-3 untwist (x_2*x_1 - q^2*y_2*y_1) * S[3][1]",
        "#   T_3_j2       level-3 untwist at j=2",
        "#   T_4          level-4 untwist at j=1, n=4",
        "# Bookkeeping (algebra of B_n itself):",
        "#   INV_i_{p,m}  x_i*y_i - 1, y_i*x_i - 1",
        "#   BR_i         x_i*x_{i+1}*x_i - x_{i+1}*x_i*x_{i+1}",
        "#   BR_i_y       same with y in place of x",
        "#   FC_{xx,yy,yx,xy}  far-commutativity between sigma_1 and sigma_3",
        "#   MIX_i_{a,b,c,d}   length-3 mixed-sign identities",
    ]
    return "\n".join(lines) + "\n"
