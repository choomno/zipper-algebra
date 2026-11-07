# Exact zipper algorithm reference

The `zipper_algebra.zipping_algorithm` package implements the variant that
continues layering after a same-size correction and restarts zipping after
slider growth. It scans bottom labels, performs upward and downward zipping,
transports framed sliders through processed words, layers components, and
recursively processes the complement. The coefficient ring is
$\mathbb{Z}[q,q^{-1}]$. Signed generator integers describe literal braid words
in algebraic product order; sliders remain explicit tokens during computation.

The mathematical analysis establishes uniform literal dispatch, whole-band
transport, preservation of the packetwise ribbon hierarchy, and termination
in the augmented nested terminal family for this variant. These statements
hold over $\mathbb{Z}[q,q^{-1}]$ without polynomial untwisting relations. The
thesis's printed restart schedule restarts after both a same-size correction
and slider growth. The termination and complexity results here supply no bound
for that printed schedule. Terminal twist collection and finite boundary
routing are included in the analyzed procedure. Termination or
confluence under arbitrary choices of algebraic rewrites also requires a
separate argument.

The executable implementation retains exact stopped states when computation
limits or defensive control checks refuse a transition. Finite executions,
independent identity checks, and installed-package tests provide separate
evidence about the implementation.

## Run a literal braid word

The default `run` command uses the packed C++ engine. It retains braid words,
exact Laurent coefficients, continuation states, collection, and scheduling in
native memory for the complete calculation. All eligible slider sizes remain
enabled. The Python reference independently verifies saved results.

```sh
python -m zipper_algebra.zipping_algorithm run --n 3 --word '[1,2,1,2]' --batch-size 64 > run.json
python -m zipper_algebra.zipping_algorithm check run.json
```

`run-native` is an explicit alias for `run`. Both accept `--trace`, `--limits`,
`--batch-size`, `--cache-bytes`, and `--cache-entries`. Their cache budgets count
estimated stored payload; process memory also includes frontiers and temporary
values. C++17, Boost headers, and OpenSSL development files are required. See
[the native engine guide](ZIPPING_NATIVE.md) for build and Python usage.

In Python, `NativeSession` in `zipper_algebra.zipping_algorithm.native` provides
`reduce`, `check`, and persistent transition caches. Its producer uses intrinsic
transition keys and restores caller contexts, endpoint labels, coefficients,
and depth exactly. Each literal input has its own formal sum.

The independent verifier executes the Python reference and retains a separate
cache of transitions it computes itself. For cache lookup it sets the coefficient
to one and depth to zero, then restores their exact values after replay. Keys
retain the full caller geometry, outer contexts, endpoint labels, resource limits,
and trace setting. The verifier does not factor out contexts or labels, and no
native transition or producer cache entry populates its cache. Its catalog clones
rule trees from trusted bytes serialized inside the verifier process; external
certificates are never deserialized with that mechanism. Compact replay omits
unused turnaround certificate construction while retaining the geometric checks.
Exact comparisons include outputs, pending states, diagnostics, deterministic
metrics, and optional traces. `session.check(result, cold=True)` uses uncached
Python collected replay for an additional audit.
Independent word processes can each retain a native session and verifier; each
native session currently runs one thread.

Explicit Python reference commands remain available for historical certificates
and differential testing: `run-reference` executes the uncollected procedure,
`run-collected` collects after each batch, and `run-memoized` uses the Python
intrinsic transition cache. Their module APIs are unchanged. The old CLI `run`
behavior is now named `run-reference`; the default output schema is
`zipper-native-run-v1`. The generic `check` command dispatches by certificate
schema and requires the matching source and build identities.

The uncollected reference example is:

```python
import json
from zipper_algebra.zipping_algorithm.algorithm import AlgorithmSession

session = AlgorithmSession()
certificate = session.reduce([1, 2], n=3)
assert certificate["status"] == "complete"
assert certificate["input_word"] == [1, 2]
assert session.check(json.loads(json.dumps(certificate)))
```

`AlgorithmSession` reuses the frozen catalog and isolates each input queue.
`algorithm.reduce_word` provides the same operation without retaining a session.
Every input letter is retained initially, including adjacent inverse pairs.
No polynomial untwisting relation is applied.

The three run statuses have distinct meanings:

| Status | Contract |
| --- | --- |
| `complete` | The pending queue is empty; `output` contains the resulting linear combination and its component reconstruction data. |
| `incomplete` | A computation limit stopped execution; `output` and `pending` retain the entire represented linear combination. |
| `control_obligation` | A defensive control check refused the next transition or geometric recognition; the diagnostic includes the exact stopped state. |

Each output term retains its coefficient, full braid-and-slider word, and component
records. These records include original endpoint labels, component words,
external routing braids, orientation, nested twist powers, and finite boundary
collars with their routing and parity data. Reconstruction uses these data
together. Terminal normalization converts mirrored coupons to the normal
coupon, transports both boundary hierarchies, and collects the input twist
powers into the output chain. Its value can be reconstructed from two bounded
exterior words and one integer twist vector per component. Ambient endpoint
subsets are separate data determined during extraction. Component records
also retain the inherited collar words and residual routes separately. Full
traces additionally retain their factorization and normalization witnesses;
these are not recoverable from the compressed value record alone.

For a component of size $k\geq2$, each exterior word has length at most
$\binom{k}{3}+\binom{k}{2}-1$: at most $\binom{k}{3}$ letters from inherited
fans and $\binom{k}{2}-1$ from residual root crossings. This bound concerns
components generated by the direct-continuation scheduler, whose fans occur
at most once per growth size on each side. Extraction checks the projected
per-fan bounds, their sum, and the final exterior-word bounds before accepting
the component. The public normalizer can also
accept supplied permutation-braid boundary factors satisfying its broader
validation contract; that validation alone does not certify the sharper
scheduler fan bound.

The standard reconstruction from bare unordered
endpoint sets $(T_{\mathrm{bottom}},T_{\mathrm{top}})$ and full nested-lasso
twist data $\tau$, using a fixed canonical coupon, cannot generally replace
these enriched records. Routing data are required for that reconstruction.
This obstruction concerns that reconstruction; alternative encodings and
dimension bounds require their own arguments.

The checker regenerates the deterministic run and compares the complete saved
result, including source and catalog hashes. It detects altered inputs, control
states, coefficients, outputs, metrics, and trace records. This replay check
relies on the local algebraic identities and geometric transitions implemented
by the scheduler. Independent local original-relator checking is described
below.

## Limits, traces, and continuation

```python
from zipper_algebra.zipping_algorithm.algorithm import reduce_word, resume
from zipper_algebra.zipping_algorithm.interpreter import Limits

stopped = reduce_word([1, 2], n=3, limits=Limits(max_steps=2))
assert stopped["status"] == "incomplete"
assert stopped["pending"]
continued = resume(stopped, limits=Limits(max_steps=1000))
assert continued["status"] == "complete"
```

`Limits` bounds transitions, live terms, word size, coefficient support, integer
bit size, and local expansion. A refused transition preserves its predecessor.
`resume` checks and regenerates the saved deterministic prefix before extending
it under the new limits. A limit hit supplies no vanishing or nonmembership
conclusion.

Full traces are optional. `trace=False` is the default for general, upward, and
two-strand runs. Compact results retain the input, exact pending work, output,
source identities, and streaming metrics. `trace=True` additionally retains
predecessors, successors, rule bindings, coefficients, and control evidence.
Both modes can be checked after JSON reload, and continuation preserves the
selected mode.

```sh
python -m zipper_algebra.zipping_algorithm run --n 3 --word '[1,2]' --trace > run.json
python -m zipper_algebra.zipping_algorithm check run.json
```

For a bounded compact run, supply limits as JSON:

```sh
python -m zipper_algebra.zipping_algorithm run --n 3 --word '[1,2]' --limits '{"max_steps":2}' > stopped.json
python -m zipper_algebra.zipping_algorithm check stopped.json
```

The first command exits `3` because its run is incomplete. Checking that exact
checkpoint exits `0`. Producer commands exit `0` on completion of their stated
scope, `2` on invalid input, and `3` on incomplete computation or an unsupported
control case. A successful `check` certifies its recorded scope, including
valid checkpoints with pending work.

## Reference complexity

For rank $n$ and input length $L$, set $A=8(n+1)^4$, $P=(n+1)^2$, and
$D=A^{3P+4}(L+n^2+1)$. In the uncollected reference model, the rewrite tree has
at most $(4^{D+1}-1)/3$ nodes. The time bound, including exact Artin checks, is
$12^{D+1}\operatorname{poly}(n,D)$. Thus the proved upper bound is
$\exp(O_n(L))$ at fixed rank. This bound controls the variant that continues
layering after a same-size correction, its enriched terminal output, and its
intermediate growth. It supplies no bound for the printed restart schedule or
for additional polynomial untwisting. The size estimate includes terminal
normalization, retained fan collars, finite routes, and conversion of mirrored
coupons to the default orientation. Default execution limits can stop a run
before completion; the bound describes the procedure with sufficient budgets.
Exact replay has the same asymptotic bound. Collection and shared structures
require their own analysis; timing fits are separate from the bound.

The collected executor uses the same local transitions. Unfolding a merged
state into its contributing histories recovers the corresponding reference
branches; scalar collection creates no new geometric transition. Its
conservative fixed-rank upper bound remains exponential in input length,
including the cost of keyed collection. The specific numerical expression
$12^{D+1}\operatorname{poly}(n,D)$ above describes the uncollected reference
cost model. No polynomial bound on intermediate collected states is established.

The linear retained-word bound and terminal compression give
$O_n((L+1)^{n-1})$ possible compressed terminal value records for words of
length at most $L$. This counts boundary words, endpoint subsets, component
order, and twist vectors; it does not count saved histories or uncollected
occurrences. It gives no polynomial bound on intermediate states or runtime.
Exhaustive enumeration additionally visits every word in the signed alphabet:
there are $[2(n-1)]^\ell$ inputs of length $\ell$, including unreduced words.

### Runtime and parallel execution

The collected executor keys pending terms by their complete continuation state,
including phase, word, frames, endpoint labels, boundary factors, and extraction
contexts. It excludes the scalar coefficient and keeps maximum path depth as
diagnostic data. Identical completed records are collected too. Coefficients
are exact sparse Laurent polynomials; addition retains negative coefficients
and deletes zeros. Collection uses an incremental dictionary, updating affected
states after each bounded batch.

Each worker applies one transition to a unit-coefficient state. Its result is
scaled by the incoming coefficient before collection. This uses linearity over
$\mathbb Z[q,q^{-1}]$; the semilinear transformation used to derive a local BM
rule does not act on the incoming coefficient. A fixed batch size determines
the schedule independently of worker count. The coordinator commits a whole
batch only after all transitions and collected limits pass; refusal retains
the preceding formal sum exactly. `max_steps` counts committed term
transitions. Term-count and coefficient-support limits apply after collection;
primitive transitions retain word, integer-bit, and expansion limits and use
local term/support floors of four/two. Batch traces retain unit rewrite
evidence, weighted successors, and collection updates. Serial replay verifies
the same certificate produced by multiple workers.

This follows the term-processing architecture described for
[FORM's parallel implementations](https://arxiv.org/html/1006.2099v1): independent
term chunks followed by collection. The
[FORM manual](https://tueda.github.io/form-doc/latest/manual/) discusses early
cancellation, worker buffers, and the cost of dispatch and final merging.
The implementation here reuses the zipper transition functions and adds no
external rewriting-engine dependency.

Runtime depends strongly on the word as well as its length and strand count.
Record completion status, transition count, output and pending term counts,
and producer and checker times together. A stopped prefix is not a timing for
the complete calculation. CPU time measures computation in a worker; wall time
also reflects contention, startup, and transfer. Comparisons should retain the
source hashes, limits, trace mode, and number of workers.

Independent input words can also run in separate processes without changing the
arithmetic. Coefficients remain sparse Laurent polynomials with arbitrary-size
integer coefficients. The uncollected executor processes one word sequentially;
the collected executor also supports workers within one word. More workers can
increase campaign throughput, subject to memory, certificate transfer, and
the single database writer. This does not improve the per-word complexity bound.

Rule selection filters the validated catalog before copying matching rules.
It preserves rule order and returns independent copies, avoiding repeated
copies of all 128 catalog records at each dispatch. This reduces interpreter
overhead without changing the rewrite schedule, coefficients, or trace.

## Frozen rules and independent local verification

`src/zipper_algebra/zipping_algorithm/data/manifest.json` records the snapshot
version, hashes, coefficient domain, and partial author-review status. There
are 32 base schemata and 128 oriented records. The original catalog digest and
the packaged digest differ because source metadata uses filenames in place of
workstation paths. Executable expressions and mathematical provenance are
retained.

```python
from zipper_algebra.zipping_algorithm import Catalog

catalog = Catalog.bundled()
instance = catalog.instantiate(
    "zip.nearest.negative.i_ge_3.I",
    {"n": 5, "i": 3, "j": 1, "k": 2, "l": 1},
)
```

Symmetry images retain the original upward parameter binding. Mirroring
preserves factor order; rotation reverses factor order and reflects indices.
The explicit descending friend-block order is part of the format. Instantiation
checks parameter guards and formal syntax; the scheduler supplies geometric
control conditions. Slider expansion is optional and bounded.

The `verification` module separately expands a bounded local instance and
checks an explicit cofactor identity against the original defining zipper
relators in the braid-group algebra. It uses the faithful Artin action to
compare braid words exactly, with $q$ indeterminate. It accounts for reflected
or rotated relators explicitly. Its certificate domain differs from the
suite's free-associative cofactor certificates over an ordered input
presentation.

```sh
python -m zipper_algebra.zipping_algorithm verify-rule single.foe.right.positive.I --parameters '{"n":2,"i":1,"j":1,"l":1}' > local-rule.json
python -m zipper_algebra.zipping_algorithm check local-rule.json
```

`verification.certify_rule` and `check_local_certificate` expose the same
operations in Python. The CLI accepts `--words` for named braid contexts and
`--limits` for verification budgets. A bounded local check establishes that
instance's identity; it does not establish dispatch coverage or all-input
termination.

## Control operations

The `upward` module starts with explicit unprocessed and processed words, a
slider, its friends, and an exterior tail. It selects a canonical rule from the
incoming crossing and requires the complete tape crossing. Every coefficient
branch retains its own control data, including the exceptional two-strand
branches. `upper_pass_complete` concerns one upward traversal.

The `framed` controller handles higher normal and mirrored sliders, both
half-twist frames, and protected braid decoration inside an inherited child
tape. `turnaround` transports a completed framed pass through typed processed
tails. It follows the slider's changing anchor through whole-band crossings
and unwraps recorded flype packets with their half twists retained. The
scheduler preserves the compact structural records required by these
operations when full tracing is disabled.

The `extraction` module recognizes supplied hereditary ribbon structure,
normalizes nested crossings, and retains reconstruction data. Its component
splitter checks terminal port sets, exterior overcrossings, the reindexed
complement, and both boundary routing braids. Coupon normalization combines
upper and lower twist data, returns the normal coupon, and retains finite
boundary routing. Mirror conversion inserts both full-component half twists
before recognition; their signs and algebraic positions are retained.
The general scheduler builds and checks these witnesses before recursive
continuation; an unsupported context returns `control_obligation`.

The narrower `transport` module checks whole-band crossings, disjoint crossings,
singleton transport, and internal two-strand twists:

```python
from zipper_algebra.zipping_algorithm.transport import transport_band, check_transport

witness = transport_band((-2, -1), n=3, size=2, start=1)
assert witness.end == 2
check_transport(witness)
```

Its checker verifies every local naturality identity and full word coverage.
Partial band crossings and higher internal crossings requiring frames are
refused with their offset and remaining word.

The `two_strand` module provides a complete procedure in ranks one and two.
Its output records singleton diagrams and two-strand sliders with their integer
twist powers. Its checker additionally substitutes the defining slider and
checks equality in the Laurent group algebra of the two-strand braid group.
That independent equality check includes pending terms in incomplete runs.
Neither this procedure nor the general scheduler asserts a canonical quotient
normal form or an equality decision procedure.

## Exact contextual replay

An explicit schedule can also be executed independently of the scheduler:

```python
from zipper_algebra.zipping_algorithm.interpreter import (
    Branch, Request, RunState, execute, check_certificate,
)

slider = ("slider", 2, 1, False, False, False)
initial = RunState(3, (Branch((2, slider, -1)),))
request = Request(
    branch=0, offset=0, rule_id="zip.outer.positive.I",
    parameters={"n": 3, "i": 2, "j": 1, "k": 0, "l": 2},
)
certificate = execute(initial, [request])
assert check_certificate(certificate)
```

The interpreter checks literal contexts, guards, complete exterior words, and
coefficient propagation. It isolates sibling control records and retains the
explicit transition history. `schedule_exhausted` means that the supplied
schedule ended. Supplied control evidence is recorded data; replay alone does
not prove scheduler reachability or completion.

## Exhaustive literal words

`enumeration.iter_words(n, length)` yields every word in the signed braid
generators at that length. Inverse pairs and braid-equivalent words remain
separate inputs. The general scheduler has its own exhaustive driver:

```python
from zipper_algebra.zipping_algorithm.general_enumeration import exhaust_general
from zipper_algebra.zipping_algorithm.interpreter import Limits

summary = exhaust_general(
    "runs/general-literal-words", ranks=(3,), wall_seconds=10,
    max_length=1, limits=Limits(max_steps=5000, max_terms=1024), trace=False,
)
```

```sh
python -m zipper_algebra.zipping_algorithm exhaust-general runs/general-literal-words-cli --ranks '[3]' --seconds 10 --max-length 1 --limits '{"max_steps":5000,"max_terms":1024}'
```

The SQLite database stores exact words, ordinals, compact or traced run
certificates, hashes, metrics, and producer/checker timings. Every run is
reloaded and checked before storage. The summary distinguishes completed
length layers from a partial layer. The wall-time budget is soft and checked
between bounded individual word computations. The first incomplete or
`control_obligation` word remains at the cursor. Larger limits can extend a
bounded checkpoint; resolving a control obligation requires an implementation
change and a fresh source identity.

Resume audits the saved prefix and binds its rank order, trace mode, sources,
catalog, configurations, and certificates. Use a fresh output directory after
source or catalog changes. The separate `enumeration.exhaust` and CLI `exhaust`
retain their ranks-one-and-two procedure and database format; they report
unsupported requested ranks explicitly. Higher-rank enumeration uses
`exhaust_general` or `exhaust-general`.

### Parallel experiments with one SQLite writer

On Unix, the repository runner computes independent words in worker processes
and commits them through one writer in the same length/rank/ordinal order:

```sh
python experiments/zipper_exhaustive_parallel.py runs/literal-words --ranks 2 3 4 5 6 7 8 --workers 4 --wall-seconds 1800 --max-words 8463
```

`--max-words` caps the total stored rows, including earlier invocations; omit
it to continue until the time budget or first unresolved word. `--max-length`
is optional. `--limits` accepts the same JSON `Limits` fields as the sequential
driver, and `--trace` retains full traces. The default is compact certificates.

The runner can extend an existing `exhaust_general` database with matching
rank order, trace mode, source, and catalog identity. Every saved row is audited
before new work begins. Both producer and checker timings are stored, and each
worker JSON-reloads and replays its result before it reaches the writer. A
worker exception or incomplete result leaves the exact checked prefix; later
speculative results are discarded. At most twice the worker count is in flight.
The soft time limit excludes the initial audit and may be exceeded while
bounded in-flight calculations finish.

The runner records its own source hash, checks source identities before commits,
and uses an exclusive file lock to reject a second parallel writer. Do not run
the sequential driver concurrently against that database: it does not use the
parallel runner's lock. Source changes require a new campaign. SQLite WAL permits
read-only progress queries while computation continues:

```sql
SELECT n, length, count(*) AS visited,
       sum(status = 'complete') AS complete
FROM words GROUP BY n, length ORDER BY length, n;
```

Selected longer-word calculations can share the same SQLite file without being
counted as exhaustive coverage. `experiments/zipper_probe_store.py` imports
reports into the separate `long_word_probes` and `probe_imports` tables. It
replays certificates by default, checks their source and summary identities,
and imports the batch atomically. The explicit `--trust-recorded-replay` option
is for reports already JSON-reloaded and replay-checked in a trusted run; the
database records that distinction. Attempts that ended without a checkpoint
remain metadata and are never counted as checked probe rows.

The process runner is a repository experiment, rather than an installed-package
entry point. Its tests cover out-of-order worker completion, interruption at an
incomplete word, exceptions, resume, source binding, exclusive writer ownership,
and a real process-pool run. The Unix-lock tests are skipped on Windows.

### Exhaustive native campaigns

```sh
python experiments/zipper_exhaustive_collected.py runs/native --execution native --targets '{"2":10,"3":10,"4":10,"5":10,"6":10,"7":10,"8":10}' --workers 8 --wall-seconds 86400
```

This runner enumerates every literal signed word in increasing length and
ordinal within each rank, including unreduced inputs. Each independent-word
worker retains a C++ session and a separate exact Python verifier. Each rank has
its own consecutive cursor; its first incomplete word pauses that rank while
other ranks continue. Increasing computation limits retries that exact word.
Targets can be extended on resume without changing the source identity.

One parent writes `words.sqlite3`, with compressed JSON certificates, hashes,
inputs, timings, limits, replay provenance, and per-rank cursors. Every new
certificate undergoes JSON reload and independent exact replay before storage.
Verifier transitions are reused only from its own prior computations; producer
and verifier cache counters and memory budgets are recorded separately.
Resume checks all saved hashes and input/configuration/cursor bindings;
`--replay-audit` additionally repeats the saved arithmetic. `summary.json` is
updated during the run and distinguishes completed layers from partial ones.
The wall budget stops admission between words; bounded in-flight work may
finish afterward. An exclusive lock prevents a second writer. Unregistered
source changes are refused; explicit version migration is described below.

The versioned runner defaults to `--execution native`. Explicit
`--execution collected` and `--execution memoized` retain historical reproduction
paths. Production and verification cache limits are exposed separately; consult
`--help` for the current options. It limits
speculative occupancy per rank as well as globally. At the wall deadline it
stops submissions and drains admitted consecutive results before exiting.
Snapshots distinguish running, paused, and completed campaigns. `--status`
and `--audit-only` are read-only.

To upgrade an existing database, first stop its writer. Explicit `--migrate`
audits its saved bindings, creates a backup, and registers a new source version
atomically. Historical identities, configurations, and completed certificates
are retained. A retried unresolved certificate is archived before replacement.
Every row is checked against its own source version; an old certificate is not
silently attributed to the new executor. An arithmetic replay of historical
certificates requires their matching source implementation.

`experiments/zipper_rule_statistics.py` produces JSON and CSV distributions
from a consistent read-only SQLite snapshot. Its primary strata are exact
strand count and word length, with layer coverage and source-version counts.
It separates control actions, transport moves, and catalog rewrites. Logical
rule counts describe committed term transitions; actual evaluations and cache
reuse are separate runtime measurements. Frequencies alone do not measure
time spent in a rule.

### Selected runtime sweeps across source versions

```sh
python experiments/zipper_runtime_sweep.py runs/literal-words/words.sqlite3 --ranks 2 3 4 5 6 7 8 --lengths 4 10 12 --workers 4 --per-word-wall 30 --wall-seconds 1800
```

For collected calculations with parallel terms within each word:

```sh
python experiments/zipper_runtime_sweep.py runs/literal-words/words.sqlite3 --execution collected --ranks 3 --lengths 10 --workers 1 --term-workers 4 --batch-size 64 --per-word-wall 300
```

Outer workers handle independent words; term workers belong to one collected
calculation. Collected mode defaults to one outer worker. Its recorded CPU
time includes the coordinator and reaped term-worker processes. A per-word
timeout terminates that word's process group, including its term workers.

This deterministic sample includes positive, negative, alternating, and seeded
random signed words. It visits length 10 first and records whether every
generator position occurs. It does not extend the exhaustive prefix.
The `runtime_campaigns` and `runtime_cases` tables bind each campaign to its
algorithm, catalog, runner, configuration, and exact input words, so new source
versions can share a database with historical results. Resume verifies that
identity and the saved certificate hashes, then visits pending cases only.

Each input runs in a fresh spawned process. Workers retain exact arithmetic
and JSON-reload and replay their certificates before storage. The parent
alone writes SQLite, with at most `--workers` inputs in flight. Producer and
checker CPU and wall times are recorded separately. The per-word wall cap
includes startup, setup, replay, and transfer. A timeout records its phase and
does not supply an accepted certificate or completed-word timing. A producer
that completed before a replay timeout is recorded as such in the diagnostic.
The campaign wall budget stops new submissions; running cases can continue to
their individual caps. Use a new campaign for a different budget or worker
count. The runner's JSON summary reports completion fractions and censored
attempts separately.

## Validation

```sh
python -m pytest tests/test_zipping_*.py -q
python -m ruff check src/zipper_algebra/zipping_algorithm tests/test_zipping_*.py
python -m mypy
```

Regressions cover literal whole-word runs, both directions, framed transport,
finite routing during extraction, recursive continuation, optional traces,
checkpoint limits, exact local original-relator checks, exhaustive-prefix
audits, and tampering after JSON reload. Finite tests and deterministic replay
provide implementation evidence for their recorded inputs. The mathematical
guarantees concern the implemented layering variant and its augmented
terminal family. The standard reconstruction from unordered endpoint sets
and full nested lassos cannot generally replace that output; other encodings
and dimension counts require separate arguments. The printed restart schedule
has no complexity bound established here.
