#!/usr/bin/env bash
# Reproduce the X_3 cofactor certificate. Kept as the historical entry point; the pipeline is now
# tools/cofactor/pipeline.py, which does the same ten phases for ANY presentation and target.
#
#   bash zipper-algebra/reproduce_X3_zero.sh          # the X_3 run, in the historical file names
#   python3 -m tools.cofactor.pipeline --help         # everything else
#
# The 500-line bash driver this replaced is archive/cofactor-pre-consolidation/reproduce_X3_zero.sh;
# notes/REMOVED.md says what each of the scripts it orchestrated became. Acceptance test for any
# change here: notes/pipeline_test/acceptance-20260917/README.md.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python3 -m tools.cofactor.pipeline --legacy-x3 "$@"
