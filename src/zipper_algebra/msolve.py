"""
Generate an f4ncgb input file (msolve format) for the zipper presentation
specialized at chosen integer (q, t) over either Q (characteristic 0) or GF(p).

Usage:
    python3 zipper_to_ms.py --n 3 --q 5 --t 7 --char 2147483647         --untwist bigelow --out zipper-3-2147483647-5-7.ms

    # Or over Q with char=0:
    python3 zipper_to_ms.py --n 3 --q 5 --t 7 --char 0 --out zipper-3-Q-5-7.ms

The resulting .ms file feeds directly into f4ncgb:
    f4ncgb -p -o basis.txt zipper-3-...ms      # -p enables cofactor proof
"""
from __future__ import annotations

import argparse

from .relations import build_relations, append_zero_sliders, append_target

__all__ = ["main"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, required=True, help="number of strands")
    ap.add_argument("--q", type=int, default=5, help="q specialization (integer)")
    ap.add_argument("--t", type=int, default=7, help="t specialization (integer)")
    ap.add_argument("--char", type=int, default=0,
                    help="coefficient characteristic (0 for Q, prime for GF(p))")
    ap.add_argument("--untwist", choices=["bigelow", "silly_question", "square-zero", "none"],
                    default="bigelow",
                    help="'bigelow': (sigma^a - t)X_2 = 0 (the higher-framing variant). "
                         "'silly_question': (sigma^2 - t^2)X_2 = 0. "
                         "'square-zero': (sigma - t)^2 X_2 = 0. "
                         "'none': no level-2 framing relation at all.")
    ap.add_argument("--framing-power", type=int, default=1, metavar="K",
                    help="In 'bigelow' mode only: use the level-2 relation "
                         "(sigma_j^K - t) * S[2][j] = 0 instead of the "
                         "canonical K=1 case. This specifies an annihilating relation; "
                         "quotient comparisons require separate verification.")
    ap.add_argument("--theta-power", type=int, default=1, metavar="M",
                    help="At level three, M > 1 uses (Theta_3^M - q^(2M)) X_3 = 0. "
                         "Higher levels retain the canonical relation; M=1 is canonical.")
    ap.add_argument("--max-zipper-level", type=int, default=None, metavar="L",
                    help="If set, include only zipper relations at level i <= L "
                         "(i.e., omit zipper relations at levels > L). Default: all.")
    ap.add_argument("--max-untwist-level", type=int, default=None, metavar="L",
                    help="If set, include only untwisting relations at level i <= L "
                         "(i.e., omit untwisting relations at levels > L). Default: all.")
    ap.add_argument("--zero-from", type=int, default=None, metavar="K",
                    help="Force X_i = 0 for all i >= K by appending S[i][j] = 0 "
                         "as extra relations for every valid (i,j). Use to test "
                         "whether forcing a level-K vanishing makes the algebra "
                         "finite-dimensional (e.g., --zero-from 4 forces X_4 = 0).")
    ap.add_argument("--out", required=True, help="output .ms path")
    ap.add_argument("--append-target", type=int, default=None, metavar="K",
                    help="add the K-strand slider as a defining relation of an augmented "
                         "presentation. This cannot prove its membership in the original ideal.")
    ap.add_argument("--target", help=argparse.SUPPRESS)
    args = ap.parse_args()
    if args.target is not None:
        ap.error("--target was removed: appending a target does not prove membership in "
                 "the original ideal. Use --append-target K only to define an augmented presentation.")

    var_names, relations, S = build_relations(args.n, args.q, args.t, args.char,
                                              args.untwist,
                                              framing_power=args.framing_power,
                                              theta_power=args.theta_power,
                                              max_zipper_level=args.max_zipper_level,
                                              max_untwist_level=args.max_untwist_level)
    if args.zero_from is not None:
        added = append_zero_sliders(relations, S, args.zero_from, args.n)
        print(f"Appended {added} zero-slider relations for i >= {args.zero_from}")
    if args.append_target is not None:
        target_idx = append_target(relations, S, args.append_target)
        print(f"Augmented presentation: added S[{args.append_target}][1] "
              f"as defining relation #{target_idx} (0-based).")

    with open(args.out, "w") as f:
        f.write(",".join(var_names) + "\n")
        f.write(str(args.char) + "\n")
        for i, rel in enumerate(relations):
            sep = "," if i < len(relations) - 1 else ""
            f.write(rel.to_msolve_string(var_names) + sep + "\n")

    print(f"Wrote {len(relations)} relations in {len(var_names)} variables to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
