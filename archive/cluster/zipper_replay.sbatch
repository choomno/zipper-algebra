#!/bin/bash
# Replay an existing anchor traced basis at the sample points listed inside
# trace_replay.gap. Cheap compared to the anchor SGrobnerTrace run, since each
# sample is a linear solve over GF(p) reusing the saved trace structure.
#
# Usage:
#   sbatch --export=ALL,ANCHOR_FILE=logs_and_traces/grobner_zipper_bigelow-trace-4-\(2^31-1\)-5-7.gap \
#          cluster/zipper_replay.sbatch
#
# To change the sample points themselves, edit the `samples` list in
# trace_replay.gap and commit; this script does not parameterize them.
#
#SBATCH --job-name=zipper_replay
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=YOUR_EMAIL@ucsb.edu

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs_and_traces

GAP_HEAP=${GAP_HEAP:-140G}

# ANCHOR_FILE is read by trace_replay.gap (with a sensible default if unset).
echo "Host:        $(hostname)"
echo "Date:        $(date)"
echo "PWD:         $(pwd)"
echo "ANCHOR_FILE: ${ANCHOR_FILE:-<default from trace_replay.gap>}"
echo "GAP:         $(which gap)"
echo "GAP_HEAP:    $GAP_HEAP"
echo "---"

gap -o $GAP_HEAP -K $GAP_HEAP trace_replay.gap
