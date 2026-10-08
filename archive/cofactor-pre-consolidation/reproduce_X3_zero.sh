#!/bin/bash
# End-to-end reproduction of  X_3 = 0  in the four-strand Bigelow zipper
# algebra Z_4, using the f4ncgb-based pipeline of WORKFLOW.md.
#
# Each phase is a separate command; output is captured to a transcript
# in logs_and_traces/reproduce-X3-zero-transcript.txt.
#
# Phase 3 (multi-point lift) uses cached basis-symbolic-Qqt.json. Redoing it
# means ~226 f4ncgb invocations: MEASURED 12.4 s per point, so ~47 min serial
# and ~8 min at the driver's 6-way parallelism -- not multi-hour. (The
# "multi-hour" figure in earlier versions of this comment belonged to the
# abandoned hypothesis-class fitter sweeps, which carried a 600 s budget PER
# SLOT across 73 slots; those are not the production path.) Produced earlier
# by tools/cofactor/run_lift_grid_canonical.py + run_fit_*.py; that lineage is documented
# in WORKFLOW.md and is reused here.

set -euo pipefail
# Repo-relative, and machine values from settings.toml — see ../DECISIONS.md.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
eval "$(python3 "$ROOT/../repo_settings.py" --sh)"
F4NCGB="$BRAID_F4NCGB"
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
# NO --target. Two distinct routes exist and must not be mixed:
#   --target + "-e on"  : f4ncgb emits the target's expanded cofactor itself.
#                         REJECTED -- extract_cofactor.py records that f4ncgb 1.0.1's
#                         expanded proofs "agree with the basis only modulo the ideal
#                         in many cases".
#   no target + "-e 0"  : the production route. 31 inputs -> 137-element basis, whose
#                         non-expanded proofs we re-expand ourselves. This is what
#                         extract_cofactor.py and symbolic_lift_x3.py consume.
# Until 2026-09-16 this line passed --target while running -e 0, producing a
# 131-element basis that no phase read.
INPUT_MS="$LOG_DIR/zipper-n4-Q-q5-t7.ms"
python3 tools/zipper_in_ms.py --n 4 --q 5 --t 7 --char 0 \
    --untwist bigelow --out "$INPUT_MS"
echo "[size] $(wc -l < "$INPUT_MS") lines in $INPUT_MS"

phase 1b "Run f4ncgb at the anchor with non-expanded proof logging"
# canonical names: these ARE the inputs phases 2 and 3 read
BASIS_OUT="$LOG_DIR/basis-n4-Q.txt"
PROOF_OUT="$LOG_DIR/proof-n4-Q-noex.txt"
# Cached on the .ms content + the exact f4ncgb flags, so phases 4-10 can be
# re-run without redoing the head of the pipeline. REFRESH=1 forces recompute.
F4_ARGS="-e 0 -v 2 -r 0"
CACHE_KEY=$(python3 tools/cache.py key --file "$INPUT_MS" --tool "$F4NCGB" --tag f4ncgb -- $F4_ARGS)
if [ "${REFRESH:-0}" = "0" ] \
   && python3 tools/cache.py get --key "$CACHE_KEY" --into "$BASIS_OUT" "$PROOF_OUT"; then
  echo "[cache] HIT  $CACHE_KEY  -- f4ncgb skipped (set REFRESH=1 to recompute)"
else
  "$F4NCGB" -i "$INPUT_MS" -o "$BASIS_OUT" -p "$PROOF_OUT" $F4_ARGS 2>&1 | tail -25
  python3 tools/cache.py put --key "$CACHE_KEY" -- "$BASIS_OUT" "$PROOF_OUT"
  echo "[cache] STORED $CACHE_KEY"
fi
echo "[size] basis: $(wc -l < "$BASIS_OUT") lines; proof: $(wc -l < "$PROOF_OUT") lines."

phase 2 "Skeleton extraction: reduce the target against the anchor basis"
# This used to be four lines of prose asserting X_3 = q^3 g_17 - q^2 g_16 +
# q(q-1) g_102 with no tool behind it, so a change in phase 1 could not disagree
# with it. The extractor derives the skeleton from the phase-1 artifacts and
# writes it as data (skeleton-n4-Q.json) for phase 3 to consume, instead of
# phase 3 hard-coding the basis indices.
TARGET="${TARGET:-3,1}"
python3 tools/cofactor/extract_cofactor.py --n "${N:-4}" --target "$TARGET" \
        --q "${ANCHOR_Q:-5}" --t "${ANCHOR_T:-7}" 2>&1 | tail -12
SKEL="$LOG_DIR/skeleton-n4-Q.json"
python3 -c "
import json, sys
s = json.load(open('$SKEL'))
idx = s['basis_indices']
print(f\"[skeleton] {s['target']} = {s['n_steps']} step(s) on basis elements {idx}, \"
      f\"remainder {s['remainder_terms']} terms\")
sys.exit(0 if s['remainder_terms'] == 0 else 1)
" || { echo "FAIL: target does not reduce to zero against the anchor basis"; exit 1; }

phase 3 "Multi-point lift (CACHED): basis-symbolic-Qqt.json"
echo "[cached] $(ls -la "$LOG_DIR/basis-symbolic-Qqt.json" | awk '{print $5, "bytes"}'):"
echo "         produced by tools/cofactor/lift_until_verified.py (sample -> lift -> compose ->"
echo "         symbolic check, growing the grid until the check passes) at 226 specializations."
echo "         Re-running fresh: ~226 points at 12.4 s each = ~8 min at 6-way parallelism."
echo "         REFRESH_LIFT=1 regenerates it instead of reading the cache."
if [ "${REFRESH_LIFT:-0}" = "1" ]; then
  python3 tools/cofactor/lift_until_verified.py 2>&1 | tail -20
fi
# The cache is only sound if it covers the skeleton phase 2 just derived. Without
# this check a change in phase 1 or 2 is silently ignored here, which is exactly
# how the pipeline came to ship a cofactor fit on a different skeleton.
python3 -c "
import json
lift = json.load(open('$LOG_DIR/basis-symbolic-Qqt.json'))
skel = json.load(open('$SKEL'))
have = {int(k) for k in lift['lifted']}
want = set(skel['basis_indices'])
print(f'[cached] lifted entries for {len(have)} basis elements; '
      f'skeleton needs {sorted(want)}')
missing = want - have
assert not missing, f'cached lift does not cover skeleton elements {sorted(missing)}'
" || { echo "FAIL: cached lift is stale with respect to phase 2"; exit 1; }

phase 4 "Symbolic composition: produce the 695-tuple Q(q,t) cofactor"
python3 tools/cofactor/compose_x3_cofactor_qt.py 2>&1 | tail -15
SYM_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-symbolic.txt"
echo "[size] $(wc -l < "$SYM_OUT") lines in $SYM_OUT"

phase 5 "Verification: free-algebra polynomial-identity check"
python3 tools/cofactor/symbolic_identity_check.py 2>&1 | tail -10

phase 6 "Compression: substitute derivable relations"
python3 tools/cofactor/compress_cofactor.py 2>&1 | tail -10
COMP_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-compressed.txt"
echo "[size] $(wc -l < "$COMP_OUT") lines in $COMP_OUT"

phase 7 "Simplification: passage to Q(q,t)[B_3] via Garside normalization"
python3 tools/cofactor/simplify_cofactor.py 2>&1 | tail -10
SIMP_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-simplified.txt"
echo "[size] $(wc -l < "$SIMP_OUT") lines in $SIMP_OUT"

phase 8 "Conjugacy reduction: eliminate R_26 via R_26 = (sigma_1 sigma_2) R_25 (sigma_1 sigma_2)^-1"
python3 tools/cofactor/reduce_to_three_rels.py 2>&1 | tail -10
THREE_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-three-rels.txt"
echo "[size] $(wc -l < "$THREE_OUT") lines in $THREE_OUT"

phase 9 "Relabeling: Z_2, T_2, T_3 convention"
python3 tools/cofactor/rename_to_ZT.py 2>&1 | tail -3
ZT_OUT="$LOG_DIR/cofactor-X3-n4-Qqt-ZT.txt"
echo "[size] $(wc -l < "$ZT_OUT") lines in $ZT_OUT"

phase 10 "Final verification: ZT file integrity + ideal-equality logic"
python3 tools/cofactor/verify_ZT_and_ideal_equality.py 2>&1 | tail -25

echo
echo "================================================================"
echo "REPRODUCTION COMPLETE"
echo "================================================================"
echo "Final artifact: $ZT_OUT"
echo "Transcript: $TRANSCRIPT"
