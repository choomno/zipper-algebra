#!/bin/bash
# Submit one traced bigelow n=4 SGrobnerTrace run on UCSB's pod cluster.
# Pod's batch nodes have 40 cores and 192 GB RAM each, so we can ask for
# generous memory without needing largemem.
#
# Usage:
#   sbatch --export=ALL,TRACED=true cluster/zipper_n4_traced.sbatch
# Override params:
#   sbatch --export=ALL,N=4,Q=7,T=11,TRACED=true cluster/zipper_n4_traced.sbatch
# Override memory:
#   sbatch --mem=180G --export=ALL,TRACED=true,GAP_HEAP=170G \
#          cluster/zipper_n4_traced.sbatch
#
#SBATCH --job-name=zipper_n4_traced
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_50,TIME_LIMIT_80
#SBATCH --mail-user=YOUR_EMAIL@ucsb.edu

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs_and_traces

# Defaults; override via --export=N=...,Q=...,etc at submit time.
export N=${N:-4}
export Q=${Q:-5}
export T=${T:-7}
export PSTR=${PSTR:-2^31-1}
export UNTWIST=${UNTWIST:-bigelow}
export TRACED=${TRACED:-false}

# GAP workspace size. Default is 140G to stay safely under the 150G --mem
# allocation while leaving room for OS/IO buffers. Bump up if --mem is bumped.
GAP_HEAP=${GAP_HEAP:-140G}

echo "Host:     $(hostname)"
echo "Date:     $(date)"
echo "PWD:      $(pwd)"
echo "Params:   N=$N Q=$Q T=$T PSTR=$PSTR UNTWIST=$UNTWIST TRACED=$TRACED"
echo "GAP:      $(which gap)"
echo "GAP_HEAP: $GAP_HEAP"
echo "---"

gap -o $GAP_HEAP -K $GAP_HEAP zipper.gap
