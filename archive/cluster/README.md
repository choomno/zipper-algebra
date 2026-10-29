# Cluster submission scripts

Submission scripts for running the zipper algebra GAP script on UCSB's Pod cluster (SLURM). Cluster use was ultimately made unnecessary by use of f4ncgb, but in principle these scripts could be ported to run those calculations as well if it were needed.

Pod's regular batch nodes have 40 cores and 192 GB RAM each; largemem has 4 fat nodes with 1 TB RAM each.

## Workflow

1. **Traced run** — one `SGrobnerTrace` at a chosen (Q, T, PSTR) point. Produces a saved traced basis file.
2. **Trace replay** —  runs `trace_replay.gap`, which iterates internally over its `samples` list (linear-solve replays using the trace structure). The `samples` list inside `trace_replay.gap` is the one place to edit when adding evaluation points.

## Files

| File | Partition | Memory | Purpose |
|---|---|---:|---|
| `zipper_n4_traced.sbatch` | batch | 150 G | Anchor SGrobnerTrace. Default first attempt. |
| `zipper_n4_traced_largemem.sbatch` | largemem | 700000 M | Same, on largemem (4 nodes; ~710 GB user-addressable, marketed as 1 TB). Fallback if batch OOMs. |
| `zipper_replay.sbatch` | batch | 150 G | Runs `trace_replay.gap` against an anchor file. |

## Traced run

```bash
# Traced bigelow n=4 anchor at q=5, t=7, p=2^31-1 — the heavy job
sbatch --export=ALL,TRACED=true cluster/zipper_n4_traced.sbatch

# Override params at submit time
sbatch --export=ALL,N=4,Q=7,T=11,TRACED=true cluster/zipper_n4_traced.sbatch

# Override memory (script reads GAP_HEAP for the gap workspace flag)
sbatch --mem=180G --export=ALL,TRACED=true,GAP_HEAP=170G \
       cluster/zipper_n4_traced.sbatch
```

Default env values (when unset): `N=4, Q=5, T=7, PSTR=2^31-1, UNTWIST=bigelow, TRACED=false`. Explicitly setting `TRACED=true` is required for an anchor run.

## Trace replay

Edit `samples` inside `../trace_replay.gap` to add/remove evaluation points, commit, then:

```bash
sbatch --export=ALL,ANCHOR_FILE=logs_and_traces/grobner_zipper_bigelow-trace-4-\(2^31-1\)-5-7.gap \
       cluster/zipper_replay.sbatch
```

If `ANCHOR_FILE` is omitted, `trace_replay.gap` falls back to its hard-coded default path.

## Monitoring

```bash
squeue -u $USER
sstat -j JOBID --format=JobID,MaxRSS,AveRSS --noheader
sacct -j JOBID --format=JobID,State,Elapsed,MaxRSS,ReqMem,ExitCode
tail -f logs_and_traces/logfile_zipper_bigelow-traced-4-\(2^31-1\)-5-7.txt
scancel JOBID
```

Pod-specific partition state:

```bash
sinfo -p batch    -h -o "%t %D" | sort | uniq -c
sinfo -p largemem -h -o "%t %D" | sort | uniq -c
```

## Before submitting

1. Edit `--mail-user=YOUR_EMAIL@ucsb.edu` in each `.sbatch` (or delete the `mail-*` lines).
2. Confirm `gap` is on `PATH` (`~/.bashrc` should have `export PATH="$HOME/gap-4.15.1:$PATH"`).
3. Parameters consumed from environment:
   - `zipper.gap`: `N`, `Q`, `T`, `PSTR`, `UNTWIST`, `S`, `TRACED`
   - `trace_replay.gap`: `ANCHOR_FILE`
   - SLURM scripts: `GAP_HEAP` for the gap workspace size
