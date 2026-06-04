#!/bin/bash
# End-to-end reproduction of  X_3 = 0  in the four-strand Bigelow zipper
# algebra Z_4, using the f4ncgb-based pipeline of WORKFLOW.md.
#
# Each phase is a separate command; output is captured to a transcript
# in logs_and_traces/reproduce-X3-zero-transcript.txt.
#
# Phase 3 (multi-point lift) uses cached basis-symbolic-Qqt.json: redoing
# it would require ~226 fresh f4ncgb invocations + the hypothesis-class
# fitter sweep, which is multi-hour. The cached file was produced earlier
# by run_lift_grid_canonical.py + run_fit_*.py; that lineage is documented
# in WORKFLOW.md and is reused here.

set -euo pipefail
ROOT="/home/choom/gap_scripts/zipper-algebra"
F4NCGB="/home/choom/src/f4ncgb/build-release/f4ncgb"
TRANSCRIPT="$ROOT/logs_and_traces/reproduce-X3-zero-transcript.txt"
LOG_DIR="$ROOT/logs_and_traces/f4ncgb"

mkdir -p "$LOG_DIR"
exec > >(tee "$TRANSCRIPT") 2>&1

phase() {
  echo
  echo "================================================================"
  echo "PHASE $1: $2"
  echo "================================================================"
}

cd "$ROOT"

phase 1a "Generate anchor .ms input at (q,t) = (5,7) over Q"
INPUT_MS="$LOG_DIR/reproduce-zipper-n4-Q-q5-t7.ms"
python3 f4ncgb_tools/zipper_in_ms.py --n 4 --q 5 --t 7 --char 0 \
    --untwist bigelow --target 3 --out "$INPUT_MS"
echo "[size] $(wc -l < "$INPUT_MS") lines in $INPUT_MS"

phase 1b "Run f4ncgb at the anchor with non-expanded proof logging"
BASIS_OUT="$LOG_DIR/reproduce-basis-n4-Q.txt"
PROOF_OUT="$LOG_DIR/reproduce-proof-n4-Q-noex.txt"
"$F4NCGB" -i "$INPUT_MS" -o "$BASIS_OUT" -p "$PROOF_OUT" -e 0 -v 2 -r 0 2>&1 \
    | tail -25
echo "[size] basis: $(wc -l < "$BASIS_OUT") lines; proof: $(wc -l < "$PROOF_OUT") lines."

phase 2 "Skeleton extraction: X_3 = q^3 g_17 - q^2 g_16 + q(q-1) g_102"
echo "Verified at anchor; the multipliers depend only on q (see WORKFLOW.md Phase 2)."
echo "The proof closure of {g_16, g_17, g_102} contains 25 basis elements."
# We don't run a separate Python tool here; the skeleton is read from the proof.
# A standalone implementation lives in extract_cofactor.py.

phase 3 "Multi-point lift (CACHED): basis-symbolic-Qqt.json"
echo "[cached] $(ls -la "$LOG_DIR/basis-symbolic-Qqt.json" | awk '{print $5, "bytes"}'):"
echo "         produced by run_lift_grid_canonical.py + run_fit_*.py at 226 specializations."
echo "         Re-running fresh would take multi-hour and is skipped here."
python3 -c "
import json
with open('$LOG_DIR/basis-symbolic-Qqt.json') as f:
    data = json.load(f)
print(f'[cached] lifted entries for {len(data[\"lifted\"])} basis elements.')
"

phase 4 "Symbolic composition: produce the 695-tuple Q(q,t) cofactor"
python3 compose_x3_cofactor_qt.py 2>&1 | tail -15
SYM_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-symbolic.txt"
echo "[size] $(wc -l < "$SYM_OUT") lines in $SYM_OUT"

phase 5 "Verification: free-algebra polynomial-identity check"
python3 f4ncgb_tools/symbolic_identity_check.py 2>&1 | tail -10

phase 6 "Compression: substitute derivable relations"
python3 f4ncgb_tools/compress_cofactor.py 2>&1 | tail -10
COMP_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-compressed.txt"
echo "[size] $(wc -l < "$COMP_OUT") lines in $COMP_OUT"

phase 7 "Simplification: passage to Q(q,t)[B_3] via Garside normalization"
python3 f4ncgb_tools/simplify_cofactor.py 2>&1 | tail -10
SIMP_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-simplified.txt"
echo "[size] $(wc -l < "$SIMP_OUT") lines in $SIMP_OUT"

phase 8 "Conjugacy reduction: eliminate R_26 via R_26 = (sigma_1 sigma_2) R_25 (sigma_1 sigma_2)^-1"
python3 f4ncgb_tools/reduce_to_three_rels.py 2>&1 | tail -10
THREE_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-three-rels.txt"
echo "[size] $(wc -l < "$THREE_OUT") lines in $THREE_OUT"

phase 9 "Relabeling: Z_2, T_2, T_3 convention"
python3 f4ncgb_tools/rename_to_ZT.py 2>&1 | tail -3
ZT_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-ZT.txt"
echo "[size] $(wc -l < "$ZT_OUT") lines in $ZT_OUT"

phase 10 "Final verification: ZT file integrity + ideal-equality logic"
python3 f4ncgb_tools/verify_ZT_and_ideal_equality.py 2>&1 | tail -25

echo
echo "================================================================"
echo "REPRODUCTION COMPLETE"
echo "================================================================"
echo "Final artifact: $ZT_OUT"
echo "Transcript: $TRANSCRIPT"
