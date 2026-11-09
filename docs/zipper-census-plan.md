# A census focused on depth in ranks three and four

The next census targets **the rank-three length-nine primitive layer**, with
**252,497 inputs**. The historical literal census already completed every word
in B3 through length eight and B4 through length six. Repeating those lower
layers with a newer implementation would be recertification, not new input
coverage. There is no separate rank-two campaign. Higher ranks remain comparison
samples rather than an exhaustive rectangle through eight strands and length ten.

The historical store contains 16,863 completed B3 length-nine literal words and
67,611 B4 length-seven words, as well as those complete lower layers. Its source
epochs, configurations, and independent replay records remain intact. Those
results are historical evidence, not automatically certified rows in the new
queue policy. The new length-nine campaign establishes one complete primitive
layer under its own recorded policy; its overlap with old results must not be
counted as newly covered words.

Of the B3 length-nine primitive requests, 16,351 already appear in the legacy
literal prefix and **236,146 do not**. The corresponding B4 length-seven
figures are 45,083 previously recorded and 153,254 absent, out of 198,337
primitive inputs. All four historical source epochs bind the same catalog
digest, but their implementation, scheduler, limits, and verifier identities
differ. The pipeline therefore retains their provenance rather than silently
importing them as current-policy certifications.

Here a primitive input means a full-support literal word whose largest bottom
strand has an undercrossing. Support separation, translation, and the proved
initial singleton peeling reconstruction explain the omitted ambient words.
These are sufficient computation inputs, not braid classes, trace classes,
distinct algebra elements, or minimal output-equivalence classes. Full ambient
coverage also requires the lower-rank reconstruction data; the primitive
counts alone are not a certificate for every ambient word. No unproved
inverse-cancellation or distant-commutation identification is used.

## Measurements used for planning

A pilot at revision `357b55e` selected one pseudorandom ordinal from each of
32 equal population intervals in each of five primitive layers, using seed
20261109. All 160 inputs completed with exact independent verification and
committed SQLite storage. The pilot used one persistent producer/verifier pair,
256 MiB per native cache, tcmalloc in the verifier, and limits of 250,000
transition evaluations, 20,000 terms, and expansion budget 1,000,000. It ran
while other project work was active, without an exclusive CPU reservation.

| Layer | Exact primitive inputs | Mean seconds | Maximum sampled seconds | Linear projection, worker-hours | Compressed leaf payload projection |
|---|---:|---:|---:|---:|---:|
| B3 length 8 | 61,937 | 0.513 | 3.61 | 8.83 | 0.216 GB |
| B3 length 9 | 252,497 | 1.182 | 5.73 | 82.92 | 0.918 GB |
| B3 length 10 | 1,022,593 | 1.494 | 4.50 | 424.30 | 3.465 GB |
| B4 length 6 | 28,379 | 0.281 | 3.19 | 2.22 | 0.114 GB |
| B4 length 7 | 198,337 | 1.113 | 20.25 | 61.31 | 0.948 GB |

The projections multiply the sample mean by the exact layer population. They
are planning estimates, not upper bounds or population confidence intervals.
The B4 length-seven median is only 0.155 seconds: one 20.25-second input
contributes more than half that sample's elapsed work. A larger sample is
needed before treating its mean as stable. Storage figures count compressed
leaf payloads only; queue records, indices, attempts, roots, journals, and
independent source epochs need additional space. The new asynchronous queue
must have its own storage growth measured rather than inheriting these figures.

The completed census prefixes are substantially easier than these deeper
stratified inputs. For example, the observed colex prefix at B3 length seven
averaged 0.084 seconds. Extrapolating that rate to long inputs would be misleading.
The pilot's verifier accounts for approximately 79–83% of total elapsed work.
Separating computation from certification improves scheduling but does not
remove that cost if every result is eventually independently checked.

## Checkpoints and resource budgets

1. **Complete B3 length nine.** The 252,497-input primitive layer projects to
   83 worker-hours
   for the combined synchronous route. Its measured producer and verifier
   components project to 13.49 and 68.87 hours respectively. One compute plus
   one certification worker therefore needs approximately 69 hours for the
   bottleneck stage, before tail and coordination allowances; it does not
   provide the ideal 42-hour result of two fully utilized synchronous workers.
   With one compute and two certification workers, the ideal bottleneck is
   approximately 34.5 hours. Allocate a 72-hour window for that split
   (96 hours if limited to one worker of each role) and 10 GB of additional
   database allowance. Treat a broader pilot of at least 512 stratified
   length-nine inputs and initial sustained production as a budget checkpoint:
   revise the worker split or window if certification throughput or unresolved
   tails are worse than projected. This is a budget, not a completion guarantee.
2. **Reassess B3 length ten and B4 length seven.** B3 length ten alone projects
   to 424 combined worker-hours. The asynchronous duration depends on the
   worker split and sustained service rates. It is a possible later campaign,
   not part of the first completion promise. Rebenchmark after verifier changes and after
   measuring sustained worker scaling. Do not infer linear speedup from the
   hardware's logical CPU count.
3. **Explore longer words separately.** For B3 lengths 12, 16, and 24, use
   finite stratified and structured cohorts with explicit seeds, per-input
   limits, and an aggregate exploration budget. Keep these cohorts in a
   separate database and label them as samples. They cannot fill missing exhaustive
   layers or support claims of exhaustive coverage.

Reserve at most 10% of the available worker budget for longer exploratory or
higher-rank comparison cases until the first complete checkpoint is certified.
Choose the compute/certify worker split from measured service rates. If a layer
requires total producer work C and verifier work V, the ideal stage time with
c compute and v certification workers is max(C/c, V/v), with storage,
initialization, contention, and tails charged additionally. Adding producers
while certification falls behind only creates a larger provisional backlog.
The backlog allowance pauses claims by compute workers, not the separate
enumeration/enqueue operation. Drain and certify existing candidates before
expanding the target. The disk allowances above are operational budgets, not
hard quotas enforced by the runner; monitor the whole database and journals.

The one-compute/two-certifier split is justified by stage balance and durable
recovery, not by a claimed speedup over three synchronous workers. At ideal
utilization, three synchronous workers would divide the 83-hour estimate by
three, while the proposed split is limited by its two certifiers. Measured
throughput, memory, and backlog determine which deployment is preferable.

Every checkpoint reports computed, certified, incomplete, and mismatching inputs
separately. Use a dedicated target-only database; do not enqueue additional
`--words` requests there. With that restriction, completion requires exactly
252,497 certified requests and zero queued, running, provisional, unresolved,
or rejected requests. The current `status` command reports aggregate counts,
not an independent audit of an arbitrary mixed request set. Halt expansion and
investigate any mismatch. Retain attempts when a limit is increased; larger
limits create a new policy. Do not sum policy counts to assert completion:
cross-policy coverage requires a separate explicit request-by-request audit.
An unresolved hard input remains a coverage hole; neither sampling nor finishing
the rest of its layer makes the layer complete. Certification is exact and
independent, with no sampling-based promotion of provisional results.

For checkpoint one, use a frozen checkout and a new database:

```sh
python -m experiments.zipper_async_census init /path/to/depth.sqlite3 \
  --limits '{"max_steps":250000,"max_terms":20000,"expansion_budget":1000000}'
python -m experiments.zipper_async_census enqueue /path/to/depth.sqlite3 \
  --targets '{"3":9}' --min-length 9 --admit-max 252497
python -m experiments.zipper_async_census run /path/to/depth.sqlite3 \
  --compute-workers 1 --certify-workers 2 --cache-mib 256 \
  --max-pending 128 --max-pending-mib 1024 --wall-seconds 259200 --drain \
  --allocator /usr/lib64/libtcmalloc_minimal.so.4.6.5
python -m experiments.zipper_async_census status /path/to/depth.sqlite3
```

The allocator path is specific to the measured host; authenticate the available
library on another machine. See [the asynchronous runner](zipper-async-census.md)
for leases, retries, source bindings, and the exact scope of candidate and
certification counts. This initial pipeline evaluates primitive requests
directly; it does not yet certify their ambient reconstruction witnesses.

## Why depth remains useful

Complete primitive layers show no saturation of the recorded outputs. B3
lengths two through six have respectively 5, 28, 91, 281, and 876 distinct
output digests. Their maximum term counts are 12, 22, 41, 67, and 106; recorded
twist magnitudes reach 1, 2, 3, 4, and 5. At length six, 785 output digests were
absent from the observed shorter primitive layers. In B4 the complete layers
three through five have 21, 190, and 1,089 distinct outputs, with maximum term
counts 58, 116, and 284.

These counts concern the full recorded formal outputs. They do not establish
linear independence, new algebra elements, or interesting identities. They do
show that the observed output data have not become constant merely because
the rank is small. Longer ordered prefixes must not be mistaken for complete
layers when comparing such growth.

The exact counts are reproducible with `NonPeelableWords(n, length).count` in
`peeling_enumeration.py`. The execution and coverage contracts are described in
[the census documentation](zipper-census-reuse.md), and the mathematical scope
of the reductions in [the input/output interface](zipper-current-normal-forms.md).
