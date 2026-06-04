#!/usr/bin/env python3
"""
Rewrite the 161-tuple three-relation cofactor with relation labels in the
Z_i (zipper) / T_j (untwisting) convention:
    R_22 (level-2 zipper)      ->  Z_2
    R_25 (level-2 untwisting)  ->  T_2
    R_28 (level-3 untwisting)  ->  T_3

Input:  logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-three-rels.txt
Output: logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-ZT.txt
"""
from __future__ import annotations
import os

ROOT = "/home/choom/gap_scripts/zipper-algebra"
SRC_FILE = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-three-rels.txt")
OUT_FILE = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-ZT.txt")

RELABEL = {"R_22": "Z_2", "R_25": "T_2", "R_28": "T_3"}


def main():
    with open(SRC_FILE) as f:
        lines = f.readlines()

    out_lines = []
    out_lines.append("# Three-relation Q(q,t) X_3 cofactor at n=4 (Bigelow zipper algebra).\n")
    out_lines.append("# Relation names use the Z_i (zipper) / T_j (untwisting) convention:\n")
    out_lines.append("#   Z_2 = level-2 zipper      (original index R_22):\n")
    out_lines.append("#         (q*y_2 + (1-q) - x_2 - q*y_1*y_2 + x_1*x_2) * S[2][1]\n")
    out_lines.append("#   T_2 = level-2 untwisting  (original index R_25):\n")
    out_lines.append("#         (x_1 - t) * S[2][1]\n")
    out_lines.append("#   T_3 = level-3 untwisting  (original index R_28):\n")
    out_lines.append("#         (x_2*x_1 - q^2*y_2*y_1) * S[3][1]\n")
    out_lines.append("#\n")
    out_lines.append("# Identity holds in the group algebra Q(q,t)[B_3]. The free-algebra\n")
    out_lines.append("# multipliers L_j and R_j are presented as canonical B_3 words in\n")
    out_lines.append("# {x_1, x_2, y_1, y_2} where x_i = sigma_i and y_i = sigma_i^{-1}.\n")
    out_lines.append("#\n")

    # Pass through header counts and global denominator from source
    src_idx_for_hist = None
    src_idx_for_denom = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if not s.startswith("#"):
            break
        if "Tuple count:" in s:
            out_lines.append(s.replace("#", "#", 1) + "\n")
        elif "R_i histogram:" in s:
            # rewrite histogram with relabeled keys
            #   {22: 37, 25: 122, 28: 2}  ->  {Z_2: 37, T_2: 122, T_3: 2}
            mapped = (
                s.replace("R_i histogram:", "Relation histogram:")
                 .replace("22:", "Z_2:")
                 .replace("25:", "T_2:")
                 .replace("28:", "T_3:")
            )
            out_lines.append(mapped + "\n")
        elif "Global denominator" in s:
            out_lines.append(s + "\n")
        elif "IDENTITY" in s or "sum_j" in s or "where L_j" in s:
            out_lines.append(s + "\n")
        elif "Verified" in s:
            out_lines.append(s + "\n")

    out_lines.append("#\n")

    # Pass through cofactor lines, replacing R_22, R_25, R_28 with Z_2, T_2, T_3.
    n_relabeled = 0
    for ln in lines:
        s = ln.rstrip("\n")
        if not s or s.startswith("#"):
            continue
        # The format is: L * R_i * R  ::  N / D
        relabeled = s
        for old, new in RELABEL.items():
            # Match " R_xx " (with adjacent spaces) to avoid false hits.
            relabeled = relabeled.replace(f" {old} ", f" {new} ")
        out_lines.append(relabeled + "\n")
        n_relabeled += 1

    with open(OUT_FILE, "w") as f:
        f.writelines(out_lines)

    print(f"Wrote {OUT_FILE}")
    print(f"Relabeled {n_relabeled} cofactor tuples.")


if __name__ == "__main__":
    main()
