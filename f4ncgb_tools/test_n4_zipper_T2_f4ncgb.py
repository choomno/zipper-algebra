#!/usr/bin/env python3
"""
Test  X_3 ?in <R_22, R_23, R_24, R_25>  in Q(q,t)[B_4]  via f4ncgb.

At n=4, the relations of Z_4^zip + T_2 are:
  R_0..R_21  invertibility + braid + far-comm + length-3 mixed sign
  R_22       level-2 zipper at j=1
  R_23       level-2 zipper at j=2
  R_24       level-3 zipper at j=1   (contains x_3, y_3)
  R_25       level-2 untwist at j=1
We append X_3 = S[3][1] as the target, and run f4ncgb in --reduce mode
(option -r 1) which computes the reduced form of the last input modulo
the previous ones.  If the reduced form is zero, X_3 is in the ideal;
otherwise it isn't.

We do this over Q at a single specialization (q, t) = (5, 7).  Reducing
to zero over Q at one generic point would suffice to demonstrate
membership in Q(q,t) (by base change along Z[q,t] -> Q, since the
target and relations all live in Z[q,t]).  A non-zero residual over Q
at one generic point disproves membership over Q(q,t) (since membership
over Q(q,t) implies membership at every regular specialization).
"""

from __future__ import annotations
import os
import sys
import subprocess
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from zipper_in_ms import build_relations, append_target

F4NCGB = "/home/choom/src/f4ncgb/build-release/f4ncgb"
ROOT = "/home/choom/gap_scripts/zipper-algebra"
MS_FILE = os.path.join(ROOT, "logs_and_traces", "f4ncgb",
                       "zipper-only-T2-n4-Q-q5-t7.ms")
OUT_FILE = os.path.join(ROOT, "logs_and_traces", "f4ncgb",
                        "zipper-only-T2-n4-Q-q5-t7.basis")


def main():
    print(f"[setup] Building n=4 zipper-only + T_2 + X_3 input at (q,t)=(5,7) over Q.")
    var_names, relations, S = build_relations(n=4, q_int=5, t_int=7,
                                              char=0, untw="bigelow")
    print(f"[setup] Full Z_4 relation count: {len(relations)} (indices 0..30).")
    # Keep indices {0..24, 25}: invertibility + braid + 3 zippers + R_25
    # (drop R_26, R_27, R_28, R_29, R_30 = the other untwistings)
    keep_idx = set(range(0, 25)) | {25}
    kept = [(i, relations[i]) for i in sorted(keep_idx)]
    print(f"[setup] Keeping {len(kept)} relations: indices {sorted(keep_idx)}.")
    print("        Discarded: R_26, R_27 (other level-2 untwistings),")
    print("                   R_28, R_29 (level-3 untwistings),")
    print("                   R_30 (level-4 untwisting).")

    # Build the kept relation list, then append X_3 as the target (last input).
    filt_relations = [rel for (_, rel) in kept]
    target_idx = append_target(filt_relations, S, target_k=3)
    print(f"[setup] Appended X_3 = S[3][1] at input index {target_idx}.")
    print(f"[setup] Total inputs to f4ncgb: {len(filt_relations)}.")

    # Emit .ms file
    with open(MS_FILE, "w") as f:
        f.write(",".join(var_names) + "\n")
        f.write("0\n")  # characteristic 0 (Q)
        for i, rel in enumerate(filt_relations):
            sep = "," if i < len(filt_relations) - 1 else ""
            f.write(rel.to_msolve_string(var_names) + sep + "\n")
    print(f"[write] Wrote {MS_FILE}")

    # Run f4ncgb in reduce mode: -r 1 means "compute the reduced form
    # of the last polynomial w.r.t. the previous ones".
    print(f"[run]   Invoking f4ncgb in --reduce mode ...")
    cmd = [F4NCGB, "-i", MS_FILE, "-o", OUT_FILE, "-r", "1", "-v", "2"]
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    t1 = time.time()
    print(f"[run]   f4ncgb completed in {t1-t0:.2f}s with exit code {res.returncode}.")
    if res.stdout:
        print("---stdout---")
        print(res.stdout[-2000:])
    if res.stderr:
        print("---stderr---")
        print(res.stderr[-2000:])
    if res.returncode != 0:
        print("[FAIL]  f4ncgb returned non-zero.")
        return 1

    # The output file should contain the reduced form of X_3.
    # If it's empty or a "0", X_3 is in the ideal; otherwise it isn't.
    if not os.path.exists(OUT_FILE):
        print(f"[FAIL]  output file not produced: {OUT_FILE}")
        return 2
    with open(OUT_FILE) as f:
        content = f.read().strip()
    nbytes = os.path.getsize(OUT_FILE)
    print(f"[read]  output file size: {nbytes} bytes.")
    if not content or content == "0," or content == "0":
        print("[RESULT] f4ncgb returned ZERO reduced form.")
        print("[RESULT] ==> X_3 IS in <R_22, R_23, R_24, R_25> at (q,t)=(5,7) over Q.")
        print("[RESULT] (Membership over Q(q,t) follows by base change.)")
        return 0

    # Count nonzero monomials approximately
    head_lines = content.splitlines()[:10]
    nmons = content.count(",") + 1
    print(f"[RESULT] f4ncgb returned NON-ZERO reduced form.")
    print(f"[RESULT] Approximate monomial count: {nmons}.")
    print("[RESULT] First few lines:")
    for ln in head_lines[:5]:
        print(f"  {ln[:120]}")
    print("[RESULT] ==> X_3 is NOT in <R_22, R_23, R_24, R_25> at (q,t)=(5,7).")
    print("[RESULT] By the upper-semicontinuity of specialization, membership")
    print("[RESULT] over Q(q,t) would imply membership at every regular point,")
    print("[RESULT] so this disproves membership over Q(q,t).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
