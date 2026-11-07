# Smaller census baselines

Counts below cover strand counts 2 through 8 and word lengths at most 10.
They exclude one globally shared identity case. They count sufficient distinct
computations, not exact equivalence classes of all possible formal outputs.
Different ambient inputs still need their endpoint embeddings and reconstruction
contexts. None of these counts is a wall-clock speedup estimate.

| Baseline | Nonempty computations | Procedure |
|---|---:|---|
| Full-support literal inputs | 56,133,831,430 | Existing controller, after support separation and translation |
| Full-support literal inputs with no initial singleton peel | 42,454,830,949 | Existing controller, full-output reconstruction |
| Full-support distant-commutation normal forms | 323,420,326 | Revised procedure required |
| Full-support reduced-trace normal forms | 105,841,628 | Revised procedure required |
| Reduced-trace normal forms with no initial singleton peel | 87,500,985 | Revised recursive procedure required |
| Also exclude decreasing signed Artin triples | 47,216,310 | Revised terminating policy; not full braid classes |
| Also use Artin moves exposing inverse cancellation | **46,047,034** | Smallest sufficient domain obtained here |

The reduced-trace/peeling row is 485.19 times smaller than the proved
existing-controller residual baseline; the additional local Artin filter
increases that ratio to 899.16. The final cancellation-exposing filter gives
921.99 times fewer jobs than the current-controller baseline, or 1,219.05
times fewer than the full-support literal baseline. These are not minimality
claims. The filters impose selected directed Artin rewrites, not full braid equality.

## Current-controller coverage

Recursively split generator support into maximal consecutive intervals,
retain each block's literal subsequence, and translate its first index to 1.
Retain the original embeddings. If a resulting full-support word on n strands
has original bottom strand n over at every incident crossing, delete that
strand, retain its top port, and recurse on the literal complement.

Every nontrivial support split lowers rank on each child, as does every peel.
Neither operation increases word length. The recursion therefore terminates
at empty inputs or at full-support words whose largest bottom strand has an
undercrossing. These are exactly the inputs counted in row two.

Compose the [block reconstruction](zipper-block-separation.md) and
[peeling reconstruction](zipper-descending-peeling.md) from the leaves upward.
Under those theorems' completion and normalization hypotheses, this recovers
the existing procedure's full terminal output. In particular, the general
block argument retains its own collection/normalization requirements; the
unary peeling prefix needs only a completed collected continuation.
Resource-limited statuses, traces, and transition counts are not identified.

`NonPeelableWords(n,length)` in `peeling_enumeration.py` implements exact
counting, rank/unrank, and ordinal interval enumeration. It uses **colex**
order: signed integer lexicographic order on reversed literal words. Its
ordinals are incompatible with the existing left-to-right census checkpoints.
No runner configuration or old SQLite results have been migrated.

The counting state is remaining length, missing generator support, the current
position p of original bottom strand n, and whether an undercrossing occurred.
Before any undercrossing, an adjacent generator permits exactly one over sign,
and a far generator permits both signs. Thus the removable count is a finite
weighted position/support automaton. Subtract it from the full-support literal
count. After an undercrossing, only missing-support cardinality matters.

At B8 length 10 this gives 30,346,444,800 full-support literal inputs minus
8,073,424,591 peelable inputs, leaving 22,273,020,209 residual computations.

## Revised procedure using reduced traces

Identify inverse cancellation and distant-generator commutation, without
imposing the adjacent three-letter Artin relation. This is the right-angled
Artin group whose commuting pairs satisfy |i-j|>=2. Its reduced traces admit
unique normal representatives. General shortlex-normal-form algorithms are
studied by [Diekert, Kausch, and Lohrey](https://arxiv.org/abs/1201.3174).
The count construction here is independent of their implementation. The
literature check for this comparison was limited to that paper's abstract and
primary-source search results; no survey of braid-group implementations was
performed.

For completeness, reduction on traces deletes inverse pairs with only
independent occurrences between them. It decreases length. Disjoint deletions
commute; overlapping deletions sharing an occurrence involve the same absolute
generator, and the ambiguity a,a^-1,a leaves the same a. The independent
intervening occurrences can be moved aside. Hence trace reduction is locally
confluent and terminating. The reduced trace is unique. Choose its least
linear extension, applying this convention to the reversed input to obtain a
colex normal word NF(w). Reversing is just a choice of normal-form order;
the final word still represents w.

Every operation is valid in the braid group and length never increases.
Consequently `w -> NF(w) -> existing zipper procedure` gives an algebraically
valid representative with the same allowed terminal conditions whenever the
existing procedure does. Inputs with equal NF have identical subsequent
calculations by construction. This does **not** assert equality with the old
controller's formal output on each original input.

To obtain the last table row, apply this policy recursively:

1. Normalize the current word to NF.
2. Split disconnected support and translate the blocks; recurse on proper blocks.
3. If original bottom n is always over, peel it and recurse on the core.
4. Otherwise invoke the existing zipper procedure on this residual input.

Reconstruction at steps 2 and 3 uses the recorded routes and embeddings.
Normalizing a peeled core again may change its old formal output. Thus this
is a revised recursive procedure, not merely a proof about one-shot NF followed
by the old controller. The algebraic equalities and the geometric extraction
reconstruction still yield valid terminal expressions. Each recursive edge
lowers rank; each normalization preserves or lowers length.

Any existing finite-dimensional spanning conclusion based on these terminal
conditions transfers to this revised procedure under the same hypotheses.
It proves none of the remaining compression, untwisting, or coefficient-ring
hypotheses by itself. Nor does this preprocessing establish arbitrary geometric
layer-partition persistence.

## Joint normal-form and peeling count

Read the prospective normal word right to left. For each signed letter a,
let I(a) be its independent letters. Maintain forbidden-letter masks F and C:

    F(pa) = (F(p) intersect I(a)) union {b in I(a): b<a}
    C(pa) = (C(p) intersect I(a)) union {-a}.

Reject a if it is in F or C. F excludes movable lexicographic inversions;
C excludes inverse pairs cancellable across a commuting suffix. Combine these
masks with used-support, the position of bottom strand n, and its undercrossing
flag. The accepted full-support words count reduced traces; requiring an
undercrossing gives residual jobs for the revised recursive procedure.
Literal peeling counts cannot be subtracted from trace counts: the joint
automaton is necessary. At B8 length 10 it gives 23,245,824 reduced traces
and 17,887,002 nonpeelable representatives.

An independent check uses the clique growth series. If
C_G(x)=sum_cliques(-x)^|clique|, signed graph-group growth is
1/C_G(2t/(1+t)). Apply inclusion-exclusion over missing generators for full
support. The two counting methods agree through all ranks and lengths in the
table. Small domains were additionally checked by enumerating raw words,
explicitly canceling trace inverse pairs, and greedily selecting least linear
extensions, with a full original-label scan for the peeling predicate.

## Decreasing local Artin rewrites

For adjacent generator indices i,j, the six signed identities

    (e*i, f*j, g*i) = (g*j, f*i, e*j)

hold for e,f,g in {-1,+1}, except the two alternating triples (+,-,+)
and (-,+,-). They follow by multiplying or inverting the Artin relation;
tests also check them with the exact Artin action on a free group. The two
excluded triples fail that same action check.

Work on the reversed word, so ordinary lexicographic comparisons implement
colex order of the literal input. Starting from reduced-trace normal form,
replace the first lexicographically decreasing valid triple and normalize
again. If normalization cancels letters, length decreases. Otherwise the
lexicographic order decreases. Thus (length,colex) strictly decreases at each
iteration and the procedure terminates. Each step preserves the braid.

Its fixed words are exactly the reduced-trace normal words having no such
forbidden contiguous triple. An additional suffix automaton counts them;
its state is the longest suffix which is a proper prefix of a forbidden
triple. Combining it with the existing support and peeling states gives
57,461,962 full-support fixed words and 47,216,310 residual jobs through the
stated limits. At B8 length10 the residual count is 11,449,648.

Use this normalization in the recursive revised procedure above. It preserves
algebraic validity and the allowed terminal family under the same hypotheses.
This local policy is deterministic and terminating, but is not asserted
confluent for the braid presentation. Equivalent braids may still produce
different fixed words. It is a sufficient workload domain, not a count of
braid elements or a proof of geometric-layer invariance. The experiment's
reference normalizer is for small-input checks, not a production implementation.

## Artin moves exposing cancellation

For each valid triple equality L=R, in both directions, also use

    L (-last(R)) -> R without its last letter
    (-first(R)) L -> R without its first letter.

These are length-four to length-two rules obtained by multiplying the same
Artin identity by an inverse letter and canceling. They strictly decrease
length even when the intermediate triple move increases lexicographic order.
Resolve competing rules deterministically and normalize the reduced trace
after each replacement. The same length/colex termination proof applies.

Extend the forbidden-substring automaton to these length-four left sides.
The resulting exact count is 56,064,932 full-support fixed words and
**46,047,034 residual jobs**, plus one identity. The last refinement removes
2.48 percent of the triple-only jobs. B8 length10 alone leaves **11,309,325**
residual jobs. These are still sufficient computation inputs for a revised
policy, not full braid classes or exact classes of the existing outputs.

## Reproduction and storage

Run from the repository root with the development environment:

```sh
PYTHONPATH=src:. python experiments/zipper_reduced_trace_counts.py \
  --database /path/to/baseline-counts.sqlite3 \
  --json-output /path/to/counts.json
```

The research database stores source-bound run definitions, exact arbitrary-
precision decimal counts, and the current-controller comparison. Repeated
runs with the same definition must agree; different source versions have
different run IDs. Select one run ID when totaling rows. The script checks
source drift before committing, and checks SQLite integrity. This database
contains aggregate counts, not millions of newly verified zipper executions.
The production census remains separate.

For the smallest domain above, run:

```sh
PYTHONPATH=src:. python experiments/zipper_artin_baseline_counts.py \
  --length-four --database /path/to/artin-counts.sqlite3 \
  --json-output /path/to/artin-counts.json
```

Omit `--length-four` to reproduce the triple-only comparison. The independent
small-domain tests verify each local identity using the faithful Artin action
on a free group, check representative idempotence and braid equality, and
compare fixed-word counts against exhaustive enumeration. These checks support
the structural proofs; no parameter specialization is used.
