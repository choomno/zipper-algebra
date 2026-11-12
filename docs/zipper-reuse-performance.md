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

## Restart integrity audit

A checkpoint restart exposed a separate cost: the full integrity audit of
39,578 stored jobs ran beyond a ten-minute admission window before any new word
was admitted. The runner had started that clock before its startup audit and
worker initialization. The admission clock now starts after setup; total wall
time and the audit, setup, and admission/drain durations remain separately
reported. This prevents a successful but long audit from consuming the entire
word-processing window.

The audit also normalized each stored terminal sum four times through nested
validation helpers. It now performs one full structural normalization, validates
the same source, request, limits, and dependency bindings, and hashes the
already canonical output with the existing canonical digest helper. Fresh JSON
parsing uses the coefficient object hook. Tests check digest equality with the
former public path, including arbitrary-width coefficients, Unicode extension
data, collection, and cancellation, and reject corrupted payloads and bindings.
The full audit remains mandatory at startup. It is an integrity check;
independent mathematical replay still occurs when a stored leaf is reused.

A read-only sample of the newest 256 stored jobs was measured in both orders:
6.349 to 1.945 seconds for baseline/candidate, then 1.949 seconds candidate and
6.367 seconds baseline. This is approximately 3.27 times as fast for that audit
sample. Separate instrumented profiles counted 1,024 versus 256 normalization
calls and attributed most of the baseline cost to normalization rather than
enumeration ordinals. These audit measurements are separate from the native
word-processing benchmarks above and do not predict a whole-store restart time.

## Typed internal verifier records and fewer temporary JSON objects

The next verifier revision replaces internal frame-transport, band-transport,
operation-step, and processed-tail-event objects with typed records. Word arrays
remain in the existing exact representation. Independent local identity,
interval, coverage, and endpoint checks still run; JSON witnesses are constructed
when a trace needs them. The public JSON witness constructors and checkers remain
available and provide a differential oracle for the typed paths. This is a
conversion of those internal records, not a conversion of the entire verifier
state engine.

Two adjacent changes remove avoidable copying. Continuation identities are
serialized directly while omitting exactly the same coefficient and depth
fields, rather than copying the whole state and deleting fields. Concatenation
borrows its input words for the duration of the call, avoiding the deep copies
previously introduced by an initializer list of owning JSON objects. The
resulting word remains independently owned. No producer implementation is used
by these changes, and no mathematical checks or exact coefficient operations
are removed.

Matched runs compared frozen revision `357b55e` against a frozen candidate source
snapshot, recording complete producer/verifier build identities and output
digests. Each run used fresh native sessions, CPU 6, 256 MiB per cache, and the
same authenticated tcmalloc allocator. The hard pair was run baseline/candidate
and then candidate/baseline. Other work continued on the machine; the core was
not exclusively reserved. The hard rows below are medians of two observations
per version; the rank-three rows are one observation per version.

| Input | Previous verifier | Typed verifier | Speed ratio |
|---|---:|---:|---:|
| B5, `[-4,-3,1,2,3,4]` | 28.339 s | 21.466 s | 1.320 |
| B8, `[4,1,6,7,2,5,6,3]` | 9.167 s | 8.207 s | 1.117 |
| B3, `[-1,-2]` repeated five times | 0.944 s | 0.638 s | 1.479 |
| B3, `[-1,-2]` repeated six times | 3.697 s | 3.178 s | 1.163 |

These are independent verification-call timings including request/response
conversion and process transfers. They exclude producer execution, startup,
output canonicalization, and SQLite storage; they must not be presented as
end-to-end census speedups. All calls completed exact independent checks, and
canonical output digests agreed across the source versions. The two hard-case
verifier medians sum to 37.506 versus 29.673 seconds, approximately 1.264 times
as fast. These small cohorts do not establish population confidence intervals
or predict the hard tail of a complete census.

A sampled process-tree RSS measurement in the second order found essentially
unchanged residency: B5 1.435 versus 1.435 GB, B8 3.000 versus 2.998 GB. The
rank-three measurements were 0.282 versus 0.283 GB at length ten and 0.354 versus
0.352 GB at length twelve. These sums include the Python driver and producer as
well as verifier; shared pages can be counted more than once. The low-rank
length-twelve example establishes completion for that word, not feasibility of
all rank-three words through length twelve.

The focused transition/replay/adversarial suite passed 77 tests, including exact
standalone trace comparisons and completed-cache modes. New native tests compare
typed witnesses against legacy JSON through rank five with both orientations,
frame directions, and several powers, and reject corrupted typed witnesses.
A supplemental corpus checks accepted four-hop whole-band and whole-frame
crossings of width at least three through rank seven, with both signs, frame
directions, seven powers, and mirrored bands. It compares canonical witness
bytes and exact exception types/messages, and mutates corresponding typed and
JSON witnesses at each step. Continuation-key tests compare thousands of states
and malformed boundary cases,
including exact exception type and message, large coefficients, Unicode,
overflowing exponents, and retained versus omitted fields. An attempted optional
AddressSanitizer/UndefinedBehaviorSanitizer supplement could not link because the
installed GCC sanitizer runtime was missing; it is not counted as validation.

## Count-only checks and compact restoration frames

The next independent verifier revision replaces expression construction used
only for a token count with a checked observer. It applies every reflection and
context shift to scalar coordinates without materializing concatenated or
mapped word arrays. Arithmetic operations retain their original order for each coordinate;
a composed map could conceal an intermediate overflow. Context routes receive
exactly the shifts from their outer contexts. On refusal, the observer calls
the original materializer so malformed states retain the original exception
and diagnostic precedence. Coefficient support, exact integer bit limits and
all restored-output checks remain in place.

DAG frames now own a typed restoration record containing rank, scalar, ordered
contexts and endpoint labels. These fields move from each owned child before
the remaining complete control state is normalized. After a frame expands,
its unit JSON is released. Its complete continuation key, restoration data,
children, collected output, cycle membership and cache bookkeeping remain.
Standalone replay borrows stable frontier entries until atomic commit;
standalone traces retain their original predecessor records and microsteps.

The implementation continues to use the independent JSON controller and rule
dispatch. A complete typed verifier would also require independent token,
controller, boundary-factor, context and formal-output representations with
canonical encoders. The JSON implementation supplies a differential oracle for
such a conversion. The current measurements do not isolate the fraction of
runtime attributable to JSON representation and do not predict a numerical
benefit for a complete conversion.

Unary transition fusion also remains a separate experiment. Each original
microstep has a continuation key and cycle check, consumes a transition on a
cache miss, accounts staged children against the active-term bound, validates
restoration and scalar multiplication, and admits completed continuations in a
specific bounded-cache order. Fusing transitions or composing nonidentity maps
must retain these events and their intermediate refusals, or define a distinct
schedule/resource contract. This revision reduces retained frame data without
changing transition counts or omitting intermediate output checks.

### Matched longer-word measurements

The frozen baseline is `9ead70b`; the production native change is `1d7cf7d`.
The cohort has one predeclared sample each of B3 and B4 lengths 10, 12, 16 and
24, plus the previously documented hard B5 and B8 words. All six arms used
fresh producer and verifier sessions per input, CPU 6, authenticated tcmalloc,
256 MiB / 262,144-entry cache budgets, 250,000 transition evaluations, 20,000
terms, expansion budget 1,000,000, a 60-second per-input deadline, and a 4 GiB
virtual-address-space guard per process. The census continued on the shared
host. Session/build setup is excluded from per-input timings.

The DAG order was baseline/candidate, then candidate/baseline. DAG rows below
are medians of the two observations for each build and input. Standalone ran
candidate/baseline once. Checker time includes authentication, encoding,
transfer, independent replay and exact comparison. Elapsed time additionally
includes normalization, result/output serialization, compression and durable
payload insertion. It is not just native execution time.

| Input | DAG checker, baseline / candidate | Standalone checker, baseline / candidate | DAG elapsed, baseline / candidate |
|---|---:|---:|---:|
| B3 length 10 | 0.923 / 0.832 s | 0.441 / 0.425 s | 1.128 / 1.046 s |
| B3 length 12 | 2.809 / 2.546 s | 0.986 / 0.839 s | 3.413 / 3.160 s |
| B3 length 16 | 15.832 / 14.315 s | 8.779 / 7.698 s | 19.050 / 17.505 s |
| B4 length 10 | 8.788 / 7.635 s | 10.158 / 8.669 s | 10.847 / 9.558 s |
| B4 length 12 | 1.390 / 1.200 s | 0.728 / 0.664 s | 1.699 / 1.484 s |
| B5 `[-4,-3,1,2,3,4]` | 23.713 / 21.507 s | 19.683 / 18.007 s | 29.359 / 26.756 s |
| B8 `[4,1,6,7,2,5,6,3]` | 8.138 / 7.976 s | 23.944 / 20.215 s | 16.681 / 16.155 s |

All seven completed inputs have byte-identical full canonical output payloads
across all six arms, with independent replay checking in every completed run.
The remaining B3 length-24, B4 length-16 and B4 length-24 inputs are incomplete
in every arm: DAG production reaches the intrinsic transition-evaluation
limit; standalone reaches the collected-term limit. None reaches the wall
clock deadline or reports an allocation/address-space refusal. These are
explicitly excluded from completed-case speed ratios and establish no longer
exhaustive coverage.

For the seven common completions, sums of the per-input DAG medians are
61.593 / 56.011 seconds for checking (1.100 ratio), and 82.177 / 75.665 seconds
elapsed (1.086). The unchanged producer also varies, 14.641 / 13.933 seconds
(1.051); host/order effects account for part of the observed difference.
Standalone checking totals 64.719 / 56.516 seconds (1.145), and elapsed totals
91.958 / 84.472 seconds (1.089), while its unchanged producer is slower in the
candidate arm, 21.417 / 22.088 seconds. These are small finite matched
observations, with two DAG observations and one standalone observation per
build/input, rather than a census-wide throughput estimate.

Sampled tree RSS gives no material memory reduction. The per-input DAG median
peaks, baseline / candidate, are 0.113 / 0.112 GiB (B3 length 10),
0.180 / 0.179 GiB (B3 length 12), 0.593 / 0.593 GiB (B3 length 16),
0.465 / 0.464 GiB (B4 length 10), and 0.141 / 0.140 GiB (B4 length 12).
Hard B5 is 1.345 / 1.345 GiB, and B8 is 2.942 / 2.950 GiB. The maximum
completed standalone samples are 1.558 / 1.575 GiB. Incomplete-case samples
are approximately 0.66–0.71 GiB for DAG and 1.07–1.67 GiB for standalone.
Tree RSS can count shared pages more than once and is distinct from the
per-process virtual-address-space guard.

The direct C++ observer oracle compares 12,000 deterministic randomized states
through ranks 1–16, signed coordinate extrema, malformed fields and cancelling
reflections/shifts against the original materializer. Full C++ replay,
completion reuse, transition, wrapper and resource-bound tests cover standalone
traces, atomic rollback, cache restoration/eviction and tampering. Independent
final-source reviews found no confirmed correctness defect. Remaining measured
work is still dominated by checking; no profile here assigns that cost to a
particular JSON operation. Future performance work includes avoiding quadratic
context-map setup, exposing valid-state fallback frequency, and auditing
redundant checks before attempting their removal.
