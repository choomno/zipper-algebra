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
