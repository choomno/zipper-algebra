# Census reuse with independent verification

See [the computation graph](zipper-computation-graph.md) for the node and edge
definitions, exact-output invariant, and independent verification argument.
The [performance measurements](zipper-reuse-performance.md) include distinct
inputs, repeated requests, cold-input regressions, and the Bloom filter trial.

`experiments.zipper_exhaustive_reuse` evaluates the current procedure through its
proved input reconstruction interface and native completed-computation cache.
The existing `NativeSession.reduce` and literal census runner retain full
standalone execution, exact replay, optional traces, and their existing storage.

```sh
python -m experiments.zipper_exhaustive_reuse /path/to/reuse-census \
  --targets '{"2":10,"3":10,"4":10,"5":10,"6":10,"7":10,"8":10}' \
  --workers 8 --wall-seconds 3600 \
  --limits '{"max_steps":1000000,"max_terms":100000,"expansion_budget":1000000}'
```

The runner completes its startup integrity audit and worker initialization,
then admits jobs for the specified time and drains admitted work. Startup time
does not consume the admission window. Its summary reports `audit_seconds`,
`setup_seconds` (including the audit), `admission_and_drain_seconds`, and the
total `wall_seconds` separately.
Resource limits remain effective inside each producer, independent verifier,
and ambient reconstruction. A conservative bound on Cartesian products,
polynomial support, coefficient bits, and route lengths is checked before
reconstruction allocates those products. A refusal records an incomplete
attempt; it does not count as completed coverage. Increase limits and rerun
under the resulting new policy namespace when needed.

Use `--words requests.json` for explicit `[rank, word]` pairs. Without this flag,
the runner enumerates full-support nonpeelable words directly by colexicographic
ordinal, fairly across ranks. Rank two uses the 19 positive/negative exponent
representatives through length ten that still require primitive calculations.
The identity is included by requesting rank one with target zero. No general
inverse cancellation, connected-support commutation, or braid relation quotient
is used. `--mode standalone` changes only the **leaf backend** to the original
producer and checker; ambient inputs still use the reconstruction interface.

Each worker owns a persistent native producer, independent native verifier, and
SQLite shard. Cache reuse is local to that worker. The coordinator's
`census.sqlite3` stores completed/incomplete/failed requested jobs, shard
references, append-only attempts, and enumeration cursors. Cursor advancement
and job admission commit together. Interrupted jobs are retried on restart;
`--retry-incomplete` additionally retries recorded incomplete and failed jobs.
Persisted rank rotation prevents short resumed runs from favoring the first
rank. `--audit-only` checks stored identities, root/dependency bindings, payload
digests, and referenced shard existence. This integrity audit does not claim
fresh mathematical replay of historical data. Read transactions hold the
coordinator snapshot and each subsequently opened shard snapshot. Completed
coordinator rows refer to previously committed, immutable shard records, so
concurrent census appends do not change the records being audited. Source
identities are checked across startup before any new jobs are admitted.

`--cache-mib` applies separately to each transition and completion cache in
both native processes. At 256 MiB this permits up to 1 GiB of estimated cache
payload per worker, in addition to active calculations, materialized outputs,
and allocator overhead. Each native process also has a rule-selection cache;
the verifier additionally caches validated bound semilinear records. Each of
these caches is capped at the smaller of 4 MiB and its configured cache budget,
with at most 4,096 entries. Shared immutable completion payloads are charged
once per cache, while every retained key has its own charge. Worker count and
cache size are independent controls; the cache estimate is not a process RSS
limit.

The result store keys completed leaves by exact request, producer source,
catalog, execution policy, limits, independent verifier identity, and Python
reconstruction source hashes. Cached content is checked on every load and
independently replayed once per fresh process before reuse. The bounded cache
of checked digests does not retain unbounded certificate objects. A verifier
resource refusal is retryable and preserves an existing leaf; a mathematical
mismatch fails the job.

Input witnesses have an independent checker. Ambient output assembly is also
checked independently by following that witness tree, multiplying exact
Laurent polynomials, restoring endpoint labels and routes, and comparing full
canonical terminal sums. Root rows store witnesses, output digests, and leaf
references rather than duplicating certificates. Reconstructed results retain
their own schema and provenance; they are not full ambient execution traces.

Counts distinguish completed requested roots, distinct verified leaves per
shard, and failed/incomplete attempts. Shard counts may duplicate a leaf across
workers. Neither reconstructed output counts nor native DAG evaluation counts
are standalone rewrite-step distributions. Full outputs remain materializable
from the stored dependencies with `ReuseCensus.evaluate`; fresh processes replay
needed leaves before reconstructing them.

Every reuse database is a new sidecar. The runner refuses a database containing
legacy census tables and does not open historical literal-census databases.
Different policies coexist without overwriting historical results.
