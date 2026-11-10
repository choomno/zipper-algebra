# Frontier and completed-continuation DAG scheduling

The retained-session comparison at base revision `9ead70b` completed the same
512 B3 length-nine literal requests in both schedules. Independent checking
accepted every result, and all 512 complete canonical formal outputs matched
as bytes, including their exact coefficients and diagram records.

For this observed traversal, producer plus independent-checker calls totalled
**650.129 seconds for the frontier and 411.129 seconds for the DAG**, an observed
frontier/DAG ratio of **1.581**. This supports retaining the current DAG option
for this measured workload. It supplies no population-wide dispatch rule.

The [machine-readable evidence](benchmarks/schedule-2026-11-09.json) records the
exact 512-input manifest, parameters, source/build hashes and measured totals.
Full certificates and canonical outputs remain in the benchmark SQLite evidence.

## Retained B3 length-nine results

| Cost or outcome | Frontier | DAG |
|---|---:|---:|
| Complete and independently checked inputs | 512/512 | 512/512 |
| Equal complete canonical outputs | 512 | 512 |
| Producer calls, seconds | 130.910 | 82.765 |
| Independent checker calls, seconds | 519.219 | 328.365 |
| Producer plus checker calls, seconds | 650.129 | 411.129 |
| Result-storage stage, seconds | 2.563 | 2.391 |
| Timed evaluation through committed storage, seconds | 654.655 | 415.500 |
| Native pair setup, seconds | 0.269 | 0.270 |
| Whole benchmark process elapsed, seconds | 664.834 | 424.237 |
| Whole process and children CPU, seconds | 658.199 | 418.794 |
| Sampled peak process-tree RSS, decimal GB | 1.957 | 2.032 |
| Session generations | 1 | 1 |
| Timeouts or resource refusals | 0 | 0 |

Call timers include Python encoding/decoding, IPC and native work. Checker
timers also include source authentication and candidate validation. The storage
stage includes result JSON encoding, compression and committed SQLite writes;
the enclosing timer also includes canonical output encoding. Whole-process
elapsed time additionally includes source guards, setup and diagnostic work.
RSS is sampled process-tree resident memory, includes retained caches, and
counts shared resident pages in each process that maps them. It is an observed
arm maximum, rather than an incremental memory cost per request.

The checker artefacts have different scopes. The frontier checker independently
replays the full standalone certificate, including execution fields. The DAG
checker independently verifies the complete formal output; it does not certify
a standalone schedule, trace or producer resource counters. Both establish the
full-output equality target used here. Call totals measure the actual supported
production and checking paths, including that difference in checking scope.

The observed ratio also varies across the frozen ascending-ordinal traversal.
For successive 32-input blocks, the first block has frontier/DAG call ratio
0.982, cases 225–256 have 1.566, and the last block has 1.716. Across individual
requests the median ratio is 1.603, with empirical 10th/90th percentiles
0.985/2.574 and range 0.190–23.334. These are descriptive ratios, with no sampling
uncertainty interpretation. Word suffix, retained-cache history and wall time
all vary with traversal position. The evidence includes every block, per-input
costs/output digests and cache snapshots at block boundaries.

## Conditions and interpretation

The exact input list is the earlier frozen 512-input B3 length-nine pilot: one
pseudorandom ordinal per equal population interval, seed `20261109`, from the
252,497 full-support nonpeelable requests. It is evaluated in frozen manifest
order. The benchmark keeps each arm's producer and independently implemented
verifier caches across distinct requests; checking time is therefore amortised
within that arm. The retained regime is specific to this ordering and sample
spacing, which differ from consecutive census enumeration.

Both arms use one sequential producer/checker pair pinned to CPU 6, limits
`max_steps=250000`, `max_terms=20000`, `expansion_budget=1000000`, and budgets
of 256 MiB and 262,144 entries per configured native cache. Transition and
completion caches have separate budgets. The DAG step limit bounds newly
expanded intrinsic transitions; the frontier uses its standalone step guard.
Equal numeric limits do not give equal logical step counts or identical
resource outcomes. Cache estimates exclude container
and allocator overhead. Each process has a 4 GiB address-space ceiling. The
frontier uses its normal batch size 64; the DAG uses its dependency evaluator.
Traces are disabled.

The producer has no explicit allocator override in either arm; its allocator
identity and full inherited environment were not separately recorded. The verifier uses
`/usr/lib64/libtcmalloc_minimal.so.4.6.5`; its metadata authenticates the file's
SHA-256 and confirms that the library is mapped. Producer and verifier binary,
source, wrapper and catalog identities match across the two arms. The active
census continued on its separate frozen source and database.

There is one fixed-order traversal per arm, DAG followed by frontier, under
concurrent host load. Before/after load and pressure snapshots are retained.
CPU affinity supplies a common allocation; it does not isolate shared caches,
I/O, memory bandwidth or frequency. The table describes measured finite costs.
It gives no replicated causal estimate, uncertainty interval, census throughput,
termination guarantee or asymptotic bound. Native verifier changes made after
this baseline can change the relative costs and require their own comparison.

## Longer words with fresh native sessions

A separate matched native benchmark supplies ten cold inputs: the previously
selected stratum-zero word at each B3/B4 length 10, 12, 16 and 24, plus historical
B5 and B8 probes from the earlier cold-optimization cohort. Only eight of the
declared sixteen long probes ran: one word from the lower ordinal half of
each layer; the upper-half word in each layer was left unrun. Each input gets a fresh producer/verifier pair and a 60-second
whole-call deadline. Deterministic limits, CPU 6, cache budgets, allocators and
4 GiB per-process address-space ceilings match the preceding policy. No OS or
build-cache flush is performed. Source `9ead70b` and native candidate `1d7cf7d`
are reported as separate within-source schedule comparisons.

| Source and scope | Frontier call seconds | DAG call seconds | Complete matches |
|---|---:|---:|---:|
| Baseline, five completed B3/B4 probes | 25.731 | 37.047 | 5 |
| Candidate, five completed B3/B4 probes | 22.691 | 32.637 | 5 |
| Baseline, including completed B5/B8 probes | 86.136 | 78.530 | 7 |
| Candidate, including completed B5/B8 probes | 78.604 | 70.580 | 7 |

These are producer plus independent-checker calls for common complete inputs.
For the candidate, their enclosing evaluation-through-storage totals are
84.472 seconds frontier and 76.219 seconds DAG; for the baseline they are
91.958 and 84.566 seconds. The larger B8 gain changes the sign of the aggregate
relative to the five completed low-rank probes. Combining these rows into a
single rank-independent recommendation would obscure that behavior.

The per-input call costs retain that distinction:

| Literal input | Baseline frontier | Baseline DAG | Candidate frontier | Candidate DAG |
|---|---:|---:|---:|---:|
| b3-l10-stratum0 | 0.542 | 1.192 | 0.529 | 1.038 |
| b3-l12-stratum0 | 1.194 | 3.443 | 1.046 | 3.094 |
| b3-l16-stratum0 | 10.753 | 19.302 | 9.479 | 17.402 |
| b4-l10-stratum0 | 12.343 | 11.319 | 10.809 | 9.630 |
| b4-l12-stratum0 | 0.898 | 1.791 | 0.829 | 1.473 |
| hard-b5 | 24.384 | 30.100 | 23.039 | 26.848 |
| hard-b8 | 36.021 | 11.383 | 32.873 | 11.095 |

B5 favors the frontier in both traversals; B8 favors the DAG. The B8 probe was
already present in the historical optimization cohort, so it is an intentional
comparison case, rather than a newly sampled estimate of higher-rank behavior.
The byte-identical producer's matched DAG call total differs by 8.9% between
the baseline and candidate traversals. This observed run variation gives no
precision estimate for the roughly 10% seven-input aggregate difference.

All four arms complete the same seven inputs. B3 length 24 and B4 lengths 16
and 24 remain incomplete in every arm. The frontier reaches its collected-term
guard at `max_terms=20000`; the DAG reaches its configured `max_steps=250000`
guard on newly expanded intrinsic transitions. These are
different resource refusals, and their partial outputs are excluded from exact
completion comparisons. No input reaches the 60-second deadline. The frontier
also independently checks incomplete results; the DAG incomplete path does
not perform that replay. Attempt costs for these unresolved requests therefore
have different checking scopes.

The three refused inputs incur the following stage costs:

| Refused-input stage, seconds | Baseline frontier | Baseline DAG | Candidate frontier | Candidate DAG |
|---|---:|---:|---:|---:|
| Producer | 10.722 | 35.309 | 10.578 | 33.991 |
| Independent partial-result checking | 35.815 | not run | 33.837 | not run |
| Storage | 2.010 | 0.010 | 1.971 | 0.010 |

The following ten-attempt totals describe workflows with different refusal
checking scopes. They cannot rank the schedules. Omitting frontier partial
checking reverses the sign of their elapsed difference, without completing
another input.

| Cost over all ten attempted inputs | Baseline frontier | Baseline DAG | Candidate frontier | Candidate DAG |
|---|---:|---:|---:|---:|
| Evaluation through committed storage, seconds | 140.506 | 119.885 | 130.860 | 110.221 |
| Cumulative native pair setup, seconds | 2.842 | 2.945 | 2.789 | 2.823 |
| Whole process elapsed, seconds | 148.597 | 127.759 | 138.857 | 117.934 |
| Sampled peak process-tree RSS, decimal GB | 1.791 | 3.169 | 1.790 | 3.167 |

These totals retain the cost of all three refusals; they are not times to solve
ten requests. Setup is cumulative over ten fresh pairs and cannot be allocated
to just the seven commonly completed inputs from the existing timer. Whole-
process time also includes teardown and diagnostic bookkeeping; their separate
costs were not timed. The frontier RSS maximum occurs at the refused B4 length-
16 probe, and the DAG maximum at the completed historical B8 probe. The exact
ten inputs, status pairs, guard diagnostics, timings, session generations and
both source manifests are included in the machine-readable evidence. Guard
thresholds and available counters are retained. For a rolled-back frontier
batch, only the largest committed frontier size survives; the attempted failing
size was not retained. All seven common outputs also match as full bytes across
all six native benchmark arms, including the repeated DAG arms.

Within each source, schedule arms have a fixed order and one traversal. Native
DAG baseline/candidate repetition, recorded in the native optimization report,
is a separate comparison. The cold longer cohort and retained B3 length-nine
cohort differ in words, lengths and cache regime; their differences do not
isolate the effect of cache retention. These measurements justify keeping both
schedules available, without an automatic rule based only on rank.

The earlier declared 16-probe long manifest was not run as a separate cohort;
using this matched ten-input native benchmark avoided duplicate measurements.

## Cold evidence and the declared subset

An earlier cold frontier traversal was interrupted after 84 completed inputs;
one additional row remains recorded as running. Its 84 stored-result calls sum
to 71.436 seconds, including 12.492 seconds of production and 58.317 seconds of
checking. The cumulative setup timer is 25.360 seconds and includes initialising
the interrupted row's native pair. There is no finished whole-process timer.
This is unpaired exploratory evidence and is excluded from the schedule ratio.

Before any cross-schedule outcome comparison, the cold plan was amended to a
fixed subset of frozen manifest indices `0,4,...,508`. That 128-input pair was
left unrun to prioritise the full retained pair and bounded longer-word work.
Its cancellation was recorded separately after the primary result became
available. The systematic subset is distinct from a newly sampled uniform
128-stratum cohort. Cold here means fresh in-memory native sessions with the usual OS and
build caches, whose startup is recorded separately.

## Strict comparison and reproduction

`experiments.zipper_schedule_report` requires both arms to have completed,
identical one-pass resource/cache configurations, full epoch source manifests,
driver hashes and native/interface identities. Each measurement must retain
its run identity. Every completed output must have independent acceptance,
a valid canonical encoding and a matching stored digest. The report records per-arm verification scopes, run timestamps, input/output
digests, the consolidated SQLite file hash and report-generator/test hashes.
It opens the source database read-only and refuses a changed database or
pending WAL pages. It then compares full output bytes and retains all status pairs, elapsed costs,
diagnostics and timeout-reset session generations. A timeout's phase is
unknown when the older worker did not preserve partial stage timers.

For a benchmark database created by the existing frozen-epoch driver:

```sh
PYTHONPATH=src:. python -m experiments.zipper_schedule_report \
  benchmark.sqlite3 FRONTIER_RUN DAG_RUN --cases-file cases.json > schedule-report.json
```

The report orients by backend even when run identifiers are supplied in reverse.
Its focused tests reject unfinished or failed arms, source/driver/configuration
and identity mismatches, missing independent acceptance, corrupt payloads and
unequal complete outputs; they retain timeout costs without claiming equal
completion. Twenty focused tests and Ruff passed for the reporting addition.
