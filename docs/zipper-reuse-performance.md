# Measured census reuse performance

The completed-continuation cache and checked input recipes reduce repeated
computation, but their benefit depends on the inputs. Cold branching inputs can
cost more than the standalone scheduler because they retain more intermediate
terminal sums. These measurements do not establish a uniform speedup or make
the remaining 42,454,828,932 primitive inputs through eight strands and length
ten tractable by themselves.

The benchmark harness is `experiments.zipper_census_reuse_benchmark`. It freezes
input words, source identities, execution settings and exact outputs in SQLite.
Timing includes producer and verifier process transfers, exact independent
checking, output canonicalization, reconstruction, compression and committed
storage. Process setup and diagnostic reporting are recorded separately.
Comparisons use the same completed inputs; unresolved cases are reported as
such. The standalone baseline already includes the verifier's normalization of
external contexts in its one-step cache.

## Distinct inputs

Measurements on an Intel Core Ultra 7 265H, with persistent native sessions,
bounded caches, and the same allocator configuration for both routes:

| Input sample | Standalone | Census reuse | Interpretation |
|---|---:|---:|---|
| 28 primitive colexicographic prefix inputs | 6.099 s | 5.705 s | 1.069 times as fast |
| Six hard or related inputs | 52.092 s | 38.847 s | 1.341 times as fast |
| Eleven commonly completed full-support inputs | 7.063 s | 9.596 s | Reuse takes 36% longer |
| Cold five-strand, length-six hard input | 24.042 s | 34.199 s | Reuse takes 42% longer |

The primitive prefix sample contains sixteen inputs in rank three at length
ten, eight in rank four at length eight, and four in rank eight at length eight.
It follows actual primitive enumeration order; it is not a random estimate of
the entire census. Peak sampled process-tree RSS was 360 MB for standalone and
373 MB for reuse on this sample. Summed RSS can count shared pages more than once.

The full-support cohort contains eighteen inputs, one at each rank three
through eight and length eight through ten. At a 10,000-evaluation guard,
standalone completed eleven and reuse completed fourteen. Evaluation counts
measure actual work in each scheduler and are not identical logical step
budgets. The three additional reuse completions were checked against standalone
runs with a larger budget. Two further positive/mixed rank-eight, length-ten
probes remained unresolved in every route at the tested limits.

For the six hard or related inputs, a second pass repeating exactly the same
roots took 39.089 s standalone and 1.653 s with census reuse. This measures
repeated requests and is not an estimate for an enumeration of distinct words.
Every compared completed output agreed exactly.

These timings precede the subsequent verifier payload-size optimization
described below. The recorded source manifests distinguish the measurements.

## Evaluated rules

Stored producer telemetry gives the following leading rule frequencies. Cache
hits do not evaluate rules, so these are actual DAG work counts, not the rule
distribution of a full standalone trace.

| Cohort | Producer evaluations | Leading rules and shares |
|---|---:|---|
| Hard/context, three unique primitive leaf calls | 146,467 | `control.absorb_internal_tape` 26.38%; `control.turnaround` 16.73%; `control.layer_absorb_internal` 10.28% |
| Primitive colex, 28 leaf calls | 23,553 | `control.turnaround` 17.32%; `control.absorb_internal_tape` 11.02%; `zip.internal.i2.negative.I` 9.06%; `control.extract` 8.71% |

These distributions already differ between the two input samples. The
benchmark retained aggregate independent-verifier cache counters but no
verifier rule histogram; its per-rule distribution cannot be recovered from
the producer's counts. The census runner records producer and verifier work
separately for subsequent campaigns.

## Bloom filter experiment

A Bloom filter can reject a cache lookup that is certainly absent. A positive
answer still requires an exact key comparison before a completed continuation
can be reused. The filter cannot identify equivalent geometric states or
replace the exact completion cache.

An isolated prototype replayed actual completed-cache query, insertion and
eviction streams from rank-three length-ten, rank-five length-six, and
rank-eight length-eight computations. It compared the existing exact index,
a blocked Bloom prefilter, and exact fingerprint/index alternatives with full
collision checks. The Bloom prefilter made cache-operation replay about 10–23%
slower on these streams. Even eliminating completed-cache lookup time entirely
would save only about 0.5% of combined producer/verifier time on the hard
five-strand case. Key construction and output payload costs were excluded from
the index-only replay and measured separately in the actual engines.

No Bloom filter is enabled. The larger measured opportunity was repeated
serialization of the same immutable verifier output solely to determine its
cache charge: 1,276 of 1,596 size computations repeated a still-live output
pointer in the eight-strand case. Those repeats consumed 2.437 s of 6.658 s
native verifier time. Each immutable result now carries its exact serialized
byte count and shares that count along identity continuations. This preserves
exact cache accounting and avoids the repeated work.

An adjacent before/after comparison used direct DAG execution with independent
verification, output encoding and committed storage, with 256 MiB caches and
the same settings. Both outputs and all recorded cache-accounting counters
agreed exactly:

| Input | Verifier before / after | End-to-end before / after |
|---|---:|---:|
| Rank five, `[-4,-3,1,2,3,4]` | 26.667 / 26.124 s | 32.138 / 31.783 s |
| Rank eight, `[4,1,6,7,2,5,6,3]` | 9.416 / 7.319 s | 17.619 / 15.633 s |

The eight-strand verifier is 1.287 times as fast and its end-to-end calculation
is 1.127 times as fast in this comparison. The combined two-input elapsed time
falls from 49.757 to 47.416 s. The five-strand transition cost still dominates.
Peak sampled RSS did not decrease; this optimization avoids repeated work,
without claiming a memory saving. These are individual paired measurements,
not population confidence intervals.

See [the computation graph](zipper-computation-graph.md) for the output-equality
argument and [the census runner](zipper-census-reuse.md) for storage, limits,
independent replay, and restart behavior.

## Cold-workload profiling and second optimization round

This round compares the preceding implementation at revision `0758be0` with
unique shared-payload accounting, bounded catalog-selection caches, cached
independently validated semilinear records, fewer output copies, and coefficient
conversion during fresh JSON parsing. The rule catalog, exact arithmetic,
transition order, input equivalences, and completed-output contract are unchanged.

The main cold regression had two distinct causes. On the five-strand input
`[-4,-3,1,2,3,4]`, the old completion cache evaluated 133,601 continuations
although only 109,380 exact intrinsic keys occurred. Charging the same immutable
output separately for each referring key caused excessive eviction. Charging
each shared payload once eliminates all 24,221 repeated evaluations in this
measurement and reduces evictions from 92,756 to 5,541. More outputs remain in
memory. This fixes avoidable recomputation but does not eliminate the difference
between the DAG and standalone schedulers: the latter evaluates 64,476 distinct
one-step transitions, with collection and cancellation between batches. Its
116,863 logical steps are a different work measure.

Independent-verifier profiling identified geometric controller replay and rule
selection as the main costs. Coefficient multiplication occupied less than 1%
of the instrumented profile. A subsequent instrumented pass counted 105,793
bound semilinear-rule derivations, taking 1.96 seconds including validation and
record construction. The final verifier computes and validates each cached
record itself: on this input it performs 298 derivations and reuses the records
105,495 times, with 1.63 MB of estimated record residency. Complete typed
environments are part of the key, and shared ownership keeps a borrowed record
alive across nested eviction. These profiling times are attribution evidence,
not uninstrumented speed measurements.

The eight-strand input `[4,1,6,7,2,5,6,3]` produces 8,994 terminal terms and a
21.65 MB JSON certificate. Python processing is material there: a fresh parse
followed by a copying coefficient conversion took 1.937 seconds, versus 0.674
seconds when the conversion occurs during parsing. These are three-run medians
for the parsing stage alone. Each parsed result remains independently mutable;
the general decoder for existing Python objects retains its nonmutating API.

### Matched end-to-end measurements

The maintained `experiments/zipper_cold_epoch_benchmark.py` driver freezes two
source epochs and applies the same benchmark helper to each. It records source
and native-binary identities and compares the full independently checked
canonical output bytes. The cold cohort uses fresh sessions and three distinct,
unchanged primitive requests. The locality cohort uses retained sessions and 28
distinct primitive requests in colexicographic order. Neither cohort counts
repeating the same root as census progress.

Two runs used opposite epoch order, baseline/candidate then candidate/baseline.
The table reports medians of those two measurements. Each benchmark process was
pinned to CPU 6 of an Intel Core Ultra 7 265H, with 256 MiB per native cache and
authenticated tcmalloc in both verifiers. The two-worker live census continued
throughout; the benchmark had no exclusive CPU reservation. Timing includes
transfers, authentication, independent exact checking, reconstruction,
canonicalization, compression, and committed SQLite writes. Initialization and
diagnostic reporting are excluded. Two observations do not supply population
confidence intervals or a census-wide forecast.

| Workload | Previous | Optimized | Speed ratio |
|---|---:|---:|---:|
| Cold B5, `[-4,-3,1,2,3,4]` | 41.258 s | 29.938 s | 1.378 |
| Cold B8, `[4,1,6,7,2,5,6,3]` | 21.632 s | 18.491 s | 1.170 |
| Cold B3, `[-1,-2]` repeated five times | 1.318 s | 1.015 s | 1.299 |
| Sum of the three cold requests | 64.208 s | 49.444 s | 1.299 |
| 28 distinct primitive colex requests | 7.142 s | 5.555 s | 1.286 |

Both cold-cohort comparisons agree closely: 63.815 to 49.647 seconds and
64.601 to 49.242 seconds. All completed cross-epoch outputs agree exactly.
Peak sampled process-tree RSS increases from 2.34 to 3.15 GB for the cold cohort;
on B5 alone it increases from 0.95 to 1.44 GB. On the locality cohort the change
is 0.372 to 0.382 GB. These are decimal GB, and summed RSS can double-count shared
pages. The configured cache budget estimates retained data; it is not an RSS
limit.

A separate bounded four-input probe completes the same two inputs in both
epochs: B7 length eight takes 1.557 to 1.370 seconds, and B8 length nine takes
6.052 to 5.143 seconds. The positive and mixed B8 length-ten probes remain
incomplete at the 10,000-evaluation guard. Faster refusal is not completed-word
throughput, and no length-ten completion improvement is established here.

### Remaining bottlenecks and scheduler choice

The optimized census with standalone leaves takes 24.842 seconds on the cold B5
case and 0.856 seconds on the negative B3 case, still faster than DAG reuse for
those individual requests. On the output-heavy B8 case it takes 41.951 seconds,
versus approximately 18.49 seconds for DAG reuse. These standalone references
are single runs using the same input reconstruction, independent verification,
storage, and limits. The default remains DAG reuse; `--mode standalone` remains
available. The measurements do not justify a uniform scheduler replacement or
an input-based dispatch rule.

On optimized cold B5, independent verification still takes approximately 25
seconds of the 30-second total. On B8, producer and verifier calls together take
about 10.8 seconds of the 18.5-second total; output conversion, canonicalization,
and storage remain material. Future work should target controller-state
construction and exact output materialization, and investigate when early
collection can reduce distinct DAG work. None of these measurements changes the
algorithm's proved word identifications or its asymptotic complexity argument.
