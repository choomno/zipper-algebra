# Essential signed-word enumeration

The census defaults to literal words whose generator support is the full
interval `1,...,n-1`, retaining every order and every inverse pair. The
[separation theorem](zipper-block-separation.md) reconstructs outputs for idle
inclusions, translations, and shuffles of disconnected blocks from these
computations. It preserves each block's literal subsequence.

An optional mode enumerates lexicographic distant-commutation classes. That
combinatorial normal form is correct, but it does not identify the existing
procedure's outputs or costs: see the [exact counterexample](TRACE_EXECUTION_COUNTEREXAMPLE.md).
The two enumeration domains are recorded separately in SQLite identities.

## Input coverage theorem

Let `w` be a literal word in the signed generators of `B_n`. Define its support
as the set of absolute generator indices appearing in `w`. Partition this set
into maximal consecutive intervals `I_1,...,I_r`, ordered from left to right.
For each interval `I_a=[s_a,t_a]`, retain the subsequence of letters with indices
in `I_a`, in its original order, and translate every absolute index down by
`s_a-1`, preserving signs. Call the resulting word `v_a`. It belongs to
`B_(t_a-s_a+2)` and contains every generator index of that group.

Then

\[
 w = \iota_1(v_1)\cdots\iota_r(v_r),
 \qquad \iota_a(\sigma_j^{\pm1})=\sigma_{j+s_a-1}^{\pm1},
\]

using only distant-generator commutations. Indeed, distinct support intervals
are separated by a missing index. Every letter in one interval therefore
commutes with every letter in another. Stable sorting occurrences by interval
produces exactly the displayed product. The inclusions act on disjoint blocks
of strands. Strands outside those blocks remain idle. The empty word gives an
empty product and retains its ambient strand count.

Within each `v_a`, one may further choose its lexicographically least word under
these same commutations. Treat the positive and negative generators as separate
letters, ordered by their signed integer values. Put an edge from occurrence
`i` to occurrence `j`, for `i<j`, whenever their absolute indices differ by at
most one. The commutation class consists precisely of the labeled linear
extensions of this dependency order. Repeatedly selecting the least available
letter gives a canonical representative. Equal-letter occurrences are always
dependent, so two occurrences of the same letter cannot be available together;
the least available labeled occurrence is unique.
Neither adjacent inverse cancellation nor the adjacent Artin braid relation is
used. In particular, `sigma_1 sigma_1^-1` remains a word of length two.

Every factor has length at most the input length and rank at most the ambient
rank. Therefore a complete representative census for ranks through `N` and
lengths through `L`, together with this construction, covers every literal input
in those bounds as an algebraic input. It retains all connected support words
that happen to simplify by other braid relations; no minimal-braid-rank claim
is made.

## Explicit witnesses and checking

`word_coverage.decompose(word, n=n, commute=True)` returns immutable blocks with
normalized words, embedding offsets, and original occurrence positions. The
concatenated positions specify a permutation of every input occurrence. Blocks
are concatenated in support order; their concatenation need not be the global
lexicographic minimum (for example, the blocks of `(3,-5)` remain in that order).
Canonicality is required separately within each block.
`check_coverage` independently checks the permutation, each embedding, full
support, separation of blocks, and canonicality. Its commutation check requires
every inverted pair of occurrences to have distant indices. This condition is
both necessary and sufficient. Necessity follows because the relative order of
a dependent pair cannot change under adjacent independent swaps. For
sufficiency, bubble sorting the permutation exchanges only those inverted
pairs. `adjacent_swaps` expands the certificate into this
sequence of elementary commutations.

The witness retains the exact input word and every occurrence, so the
representative computation does not erase the literal-word provenance. It
also permits independent verification without importing the zipper producer
or its rewrite catalog.

## Direct enumeration theorem

Fix rank `n`, put `m=n-1`, and let

\[
 A_m=\{-m,\ldots,-1,1,\ldots,m\},\qquad
 I(a)=\{b\in A_m:||a|-|b||\ge 2\}.
\]

Order letters by signed integer value. For a prefix `p`, define its forbidden
set `F(p)` to contain `a` exactly when `p` has an occurrence of a letter `b>a`
such that every letter from that occurrence through the end of `p` belongs to
`I(a)`. Appending `a` in this situation creates a movable inversion: it can
commute left past that entire suffix, giving a lexicographically smaller word.
The initial set is empty. Directly from this definition,

\[
 F(pa)=(F(p)\cap I(a))\cup\{b\in I(a):b<a\}.
\]

An old witness for a forbidden `b` survives precisely when the new letter `a`
commutes with `b`; the new final occurrence of `a` is itself a witness precisely
when `b<a` and the two letters commute.

**Canonicality.** A word is lexicographically least in its commutation class
if and only if no prefix extension appends a forbidden letter. Necessity
follows by moving that letter left. For sufficiency, use the dependency order
of occurrences above. If a word is not the lexicographically least labeled
linear extension, consider its first deviation from greedy choice of the least
available letter. A smaller available occurrence lies farther to the right.
Every occurrence between the chosen position and that occurrence must commute
with it: a preceding dependent occurrence would be its predecessor in the
poset, contradicting availability. It therefore gives precisely a movable
inversion in the word. Equal-letter occurrences cannot create an ambiguity,
since they are dependent and their order is fixed. Thus the forbidden-set
condition accepts exactly one word per commutation class. Every prefix of an
accepted word satisfies the same condition.

**Completion counts.** Let `M` be the set of generator indices still missing
from the prefix, and let `r` be the number of remaining letters. Define

\[
 \begin{aligned}
 C(0,M,F)&=\mathbf 1_{M=\varnothing},\\
 C(r,M,F)&=\sum_{a\in A_m\setminus F}
 C\bigl(r-1,M\setminus\{|a|\},
       (F\cap I(a))\cup\{b\in I(a):b<a\}\bigr),\quad r>0.
 \end{aligned}
\]

If `|M|>r`, the count is zero. Induction on `r` proves that this recurrence
counts exactly the accepted completions having full support: their first
letters partition them into disjoint summands and the displayed update gives
the remaining conditions. The desired layer size is
`C(L,{1,...,m},empty)`. The implementation stores sets as bit masks and caches
these integer counts per layer.

At each prefix, the completions partition into consecutive lexicographic
intervals, one for each allowed next letter in increasing order, with sizes
given by the recurrence. Unranking subtracts the sizes of earlier intervals
and enters the unique interval containing the ordinal. Ranking adds those
earlier sizes while following the supplied word, rejecting forbidden letters
or incomplete support. Induction on remaining length proves that these maps
are mutually inverse between the interval of integers from zero to one less
than the layer size and the canonical full-support words. Interval iteration
uses the same disjoint partition and prunes subtrees outside the requested
ordinal interval. It neither enumerates earlier ordinals nor filters the full
literal-word universe.

When commutations are disabled, only the missing-support condition remains.
Writing `D(r,k)` for its completion count with `k` missing indices gives

\[
 D(r,k)=2(m-k)D(r-1,k)+2kD(r-1,k-1),
 \qquad D(0,k)=\mathbf 1_{k=0}.
\]

This counts literal full-support words with both signs distinct. Rank one has
one empty word and no positive-length words in either mode.

## Independent count check

For a subset `S` of generator indices, form the signed alphabet on `S`, retaining
the original index distances. Its clique polynomial is

\[
 P_S(t)=\sum_{J\subseteq S\;:\;J\text{ has no adjacent indices}}
                 (-2t)^{|J|}.
\]

A clique contains at most one signed letter of each index, and choosing signs
accounts for the factor `2^{|J|}`. Its trace growth series is `1/P_S(t)`.
For completeness, in the completed trace-monoid algebra multiply the formal
sum of all traces by the signed sum of clique traces. For a nonempty trace,
its clique right divisors are precisely the subsets of the maximal occurrences
in its dependency poset. Distinct maximal occurrences commute pairwise, so
the coefficient is the alternating subset sum `(1-1)^d=0`, where `d>0` is their
number. The empty trace has coefficient one. Specializing each letter to `t`
gives the claimed inverse. Each coefficient is finite, so this formal argument
needs no analytic convergence assumption.

Inclusion-exclusion over omitted generator indices therefore gives the
full-support count independently of the prefix automaton:

\[
 E_{m,L}=\sum_{S\subseteq\{1,\ldots,m\}}
 (-1)^{m-|S|}[t^L]\frac1{P_S(t)}.
\]

For the complete interval,
`P_{1,...,m}(t)=sum_k (-2t)^k binom(m-k+1,k)`; subsets with gaps must retain
those gaps in their independence relation. The tests compute the subset
polynomials and their inverse coefficients directly, then compare against
prefix-automaton counts through rank eight at length ten. Small layers are
also checked independently by exhaustively exploring adjacent commuting swaps.

The resulting counts summed over lengths zero through ten are:

| Strand count | Full-support commutation classes |
| --- | ---: |
| 2 | 2,046 |
| 3 | 1,394,008 |
| 4 | 19,606,864 |
| 5 | 62,662,528 |
| 6 | 95,367,296 |
| 7 | 89,418,240 |
| 8 | 54,969,344 |
| **Total** | **323,420,326** |

Adding one global empty-word case gives **323,420,327** representative inputs.
These are enumeration counts. They establish neither that the original zipper
controller factors over disconnected support nor that all omitted literal
executions terminate or have equivalent traces, costs, or layered outputs.

## Computations and coverage claims

`experiments/zipper_exhaustive_essential.py` maintains a separate experimental
`essential.sqlite3`. Its default domain is the global identity plus full-support
literal words in ranks two through eight, through length ten. This is
56,133,831,430 nonempty inputs, plus one global identity. The optional
`--commutation-classes` mode instead has the trace counts tabulated above;
it is an experimental output-coverage domain for the existing procedure. Each
saved representative is checked by the independent C++ verifier. The database
binds the enumeration source, engine builds, limits, and trace setting; it stores
resumable per-rank prefixes and archives incomplete attempts on retry.

For a bounded validation run:

```sh
PYTHONPATH=src:. python experiments/zipper_exhaustive_essential.py runs/essential-pilot \
  --targets '{"1":0,"2":3,"3":3}' --workers 2 --wall-seconds 60
```

Use `--trace` to retain calculation traces, and `--audit-only` to check saved
hashes and bindings through a read-only connection. This audit does not repeat
mathematical replay; add `--replay-audit` for independent C++ replay in the
same read transaction. Verification failures retain their inputs, available
certificates, and diagnostics; `--retry-failures` explicitly retries them.
Resume requires the same source identity; trace settings
belong in separate campaigns. The runner leaves the legacy literal census
untouched and does not import its certificates or declare omitted inputs checked.

If a representative has a verified algebraic derivation, its translation is a
valid derivation in the ambient algebra whenever the defining identities are
stable under the standard strand embedding. This uses the induced algebra
homomorphism and requires no injectivity assertion for the quotient algebras.
Multiplying the factor derivations and
prepending the checked commutations then gives a derivation for the original
input. Its output is a product of translated factor outputs. Declaring this
product to be the output of the existing procedure additionally requires
its component order and routing. The separation theorem proves precisely
that reconstruction when each block retains its original literal order.

An input witness alone is not a certificate of zipper reduction, termination,
normal form, complexity, or deterministic scheduler equivalence. In particular,
two commutation-equivalent inputs can take different paths and have different
costs under the existing scheduler. Representative counts must be reported
separately from literal inputs directly executed, and a completed representative
must have its own independently checked computation certificate. Existing
literal-word SQLite rows retain their original meanings.

The [reference separation proof](zipper-block-separation.md) covers catalog
rewrites, phase changes, packet contiguity, layering, extraction, normalization,
and coordinate routing. Its strong-normalization argument justifies exact
collection. It identifies routed terminal sums, not resource counters or
bounded stopping statuses. Native implementation validation remains distinct
from the reference-transition theorem.

## Background

For the standard trace-monoid interpretation and inverse clique-polynomial
counting identity, see Bouillard and Mairesse, *Möbius inversion formula for the
trace group* (2004), §§2–3,
[primary text](https://www.numdam.org/item/10.1016/j.crma.2004.10.017.pdf).
Here the alphabet has two letters for each generator index and independence is
exactly `abs(abs(a)-abs(b)) >= 2`; it is a trace **monoid** on signed letters,
without cancellation. The elementary input factorization and witness argument
above provide the coverage proof directly.

## Why the stronger controller invariant matters

The [separation theorem](zipper-block-separation.md) projects every control
field and preserves the literal subsequence on every other block. In addition
to monochromatic crossings and coupons, it needs the [reachable packet
grammar](zipper-packet-contiguity.md): required crossing packets are emitted
intact and remain accessible throughout growth, both turns, and layering.

Monochromaticity alone is insufficient. A synthetic mirrored size-three
framed state in rank four, at position two, with left word `(1,2,3)`, tail
`(1,)`, output power `-1`, input power `1`, and policy `bar_mirror` admits
`zip.flype.negative.BM`. Translating by two indices to rank six and inserting
spectator `1` inside the packet gives left word `(3,4,1,5)`, coupon position
four, and tail `(3,)`. The ambient controller reports a missing full-tape rule,
although deleting the spectator recovers the first state. Moving the spectator
to give `(1,3,4,5)` restores the rule. Python and C++ transitions agree.

The obstructed state splits a required packet and violates the stronger
reachable-state invariant. It is not a counterexample from a literal input.
The packet and controller proofs establish separation on reachable states;
they make no deletion-equivariance claim for arbitrary well-formed framed
states. This distinction is separate from the actual literal-input trace
counterexample, which has connected support.
