# Exact zipper algorithm reference

The `zipper_algebra.zipping_algorithm` package implements a deterministic
reconstruction of the thesis's zipper procedure on literal braid words. It
scans bottom labels, performs upward and downward zipping, transports framed
sliders through processed words, layers components, and recursively processes
the complement. The coefficient ring is $\mathbb{Z}[q,q^{-1}]$. Signed generator
integers describe braid words in algebraic product order; sliders remain
explicit tokens during computation.

The mathematical analysis of the corrected prescribed schedule establishes
uniform literal dispatch, whole-band transport, preservation of the packetwise
ribbon hierarchy, and termination in the augmented nested terminal family.
These statements hold over $\mathbb{Z}[q,q^{-1}]$ without polynomial untwisting
relations. They concern this schedule; termination or confluence under
arbitrary choices of algebraic rewrites requires a separate argument.

The executable implementation retains exact stopped states when computation
limits or defensive control checks refuse a transition. Finite executions,
independent identity checks, and installed-package tests provide separate
evidence about the implementation.

## Run a literal braid word

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
together. The terminal-family theorem retains this finite boundary data.
Compression to bare $(T_{\mathrm{bottom}},T_{\mathrm{top}},\tau)$ records, or to
component sizes and twist powers alone, remains unproved. A dimension count
based on such compression also remains unproved.

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
$12^D\operatorname{poly}(n,D)$. Thus the proved upper bound is
$\exp(O_n(L))$ at fixed rank. This bound controls the prescribed schedule and
its intermediate growth. Collection, shared structures, optional traces, and
certificate checking have their own measured costs; timing fits are separate
from the bound.

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
upper and lower twist data and retains finite boundary routing and orientation.
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
guarantees concern the corrected prescribed schedule and its augmented
terminal family; bare-tuple compression and the associated dimension count
remain separate open questions.
