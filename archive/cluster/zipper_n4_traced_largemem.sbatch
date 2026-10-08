#!/bin/bash
# Fallback submission for the n=4 traced run on UCSB's pod cluster if the
# batch-partition version (150G) OOMs. Uses largemem (4 nodes; ~710 GB
# user-addressable per node despite being marketed as 1 TB).
# Usage:
#   sbatch --export=ALL,TRACED=true cluster/zipper_n4_traced_largemem.sbatch
#
#SBATCH --job-name=zipper_n4_traced
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=largemem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=700000M
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_50,TIME_LIMIT_80
#SBATCH --mail-user=YOUR_EMAIL@ucsb.edu

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs_and_traces

export N=${N:-4}
export Q=${Q:-5}
export T=${T:-7}
export PSTR=${PSTR:-2^31-1}
export UNTWIST=${UNTWIST:-bigelow}
export TRACED=${TRACED:-false}

GAP_HEAP=${GAP_HEAP:-650G}

echo "Host:     $(hostname)"
echo "Date:     $(date)"
echo "PWD:      $(pwd)"
echo "Params:   N=$N Q=$Q T=$T PSTR=$PSTR UNTWIST=$UNTWIST TRACED=$TRACED"
echo "GAP:      $(which gap)"
echo "GAP_HEAP: $GAP_HEAP"
echo "---"

gap -o $GAP_HEAP -K $GAP_HEAP zipper.gap
