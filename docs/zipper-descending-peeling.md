# Descending initial components and arbitrary retained cores

This is a source-level theorem for the existing reference controller. It gives a sufficient condition for equality of full formal outputs. The unrestricted single-use far-commutation conjecture remains open.

Fix 1 <= m <= n. Label the bottom endpoints 1,...,n and follow these labels through the literal input w from right to left. At a letter e*i, write a,b for the current labels at positions i,i+1. The overstrand is a when e=+1 and b when e=-1. Assume:

(P_m) At every crossing for which max(a,b)>m, the larger label is over.

Let v be the literal B_m word obtained by deleting all strands with bottom labels >m, retaining the signs and order of the remaining crossings and renumbering their positions. Let p_k be the original top position reached by bottom label k, for k=n,...,m+1. No condition is imposed on v.

## Theorem

Provided resources suffice, the reference computation on w begins by extracting the singleton components n,n-1,...,m+1, each with coefficient one and empty internal word. It then executes precisely the computation on v, with a fixed prefix of extraction contexts, an inherited depth offset, and the order-preserving relabeling of its top endpoints onto the unoccupied original top positions. Its bottom endpoint labels remain 1,...,m.

Consequently the full collected terminal output (coefficients, words, and ordered component records) is determined by the smaller computation Z_m(v) and the finite data

    (n,m,p_n,p_(n-1),...,p_(m+1)).

In particular, two words satisfying P_m and having identical v and identical ordered peeled top positions have equal full collected output. This theorem identifies terminal outputs, not event traces or input provenance.

## Proof of the initial extraction sequence

The controller starts its single scan at the rightmost bottom port, hence at label n. Its scan reads w right to left. `_single_step` declares a crossing offensive exactly when the active strand is adjacent and under. Every crossing incident to label n has n over by P_m. Nonincident crossings are far steps. Thus every transition is unary with unchanged coefficient; the processed literal word is unchanged, and the active position follows label n. The empty-left control transition creates the size-one slider at its top position.

The extraction routine applies `split_layered`. Every mixed crossing has the selected singleton over, so its explicit recognition condition holds. The selected internal word becomes empty. Its complement is precisely the literal deletion of strand n: `split_layered` retains every crossing between complement strands with its sign and order and renumbers positions by the number of deleted ports to its left. Its routes are the canonical stable-partition routes determined by its top and bottom ports.

Extraction restarts `_single` on this complement, with boundaries empty and world_phi false, appending the singleton context. Deletion preserves the original-label identities, over/under choices, and order of all crossings between surviving strands. Therefore P_m continues to hold. The largest remaining bottom label is n-1. Repeat down to m+1. This proves the stated prefix and exact literal remaining word v.

## Continuation and complete output fields

After the prefix, bottom labels are 1,...,m. The top-label list is the increasing list Q of original top positions other than p_n,...,p_(m+1). It need not equal 1,...,m. Relative to a standalone computation on v, the larger computation has only these differences: its frozen prefix contexts, its depth offset, and the replacement of each original top port r by Q[r].

The reference geometric transitions do not inspect the values of the top/bottom label lists. These lists are read when extraction records original ports, and updated by order-preserving deletion. Contexts are frozen and appended; they are read to construct the terminal ambient word. Hence induction on transitions gives exactly the same relative geometric states and coefficients. The relabeling and appending a common prefix are injective on state and terminal keys, so they preserve equality classes used in collection. This is a conjugacy of the collected continuations, not an appeal to confluence or to linear independence of terminal diagrams.

For completeness, the prefix contexts can be reconstructed from the stated finite data. Before extracting k, its local bottom port is k; its local top port is the position of p_k in the increasing list of top positions remaining after extracting n,...,k+1. Its rank is k, size is one, component_word is empty, and its intrinsic record has empty twist powers and routes and normal orientation. The context boundary routes are `_over_route(k,(k,))` below and inverse `_over_route(k,(local_top,))` above. Append the smaller computation's contexts after applying the Q-label relabeling to original top ports. Applying `_lift` to this combined list constructs exactly the full output word. Thus the theorem covers route words and all component fields, not only twists or endpoints. Q is determined by the peeled positions.

For this unary prefix, successful finite collected continuation on v suffices. A separate assumption that an uncollected expansion tree is finite is unnecessary. Caps must allow the larger prefix, lifted contexts, and total steps; a limited run is not equated with a completed run.

## Exact logical transition count

Let L(v) be the length of the literal retained word. For each discarded crossing c, let b(c)=max(a,b) be its larger original bottom label. Let T count controller transitions, including scan-completion and extraction, using the same continuation scheduler and sufficient resources. Then

    T(w) = T(v) + 2(n-m) + (n-m)L(v)
           + sum_{discarded c} (n-b(c)+1).

Indeed a retained crossing survives every peeled scan and is read n-m times before the continuation. A discarded crossing survives and is read exactly during scans n,n-1,...,b(c), and is then deleted. Each peeled component has one scan-completion and one extraction transition. The continuation counts agree by the injective state correspondence. This counts logical work, not runtime, allocation costs, or verifier work.

## Recognizing the maximal prefix

During one backwards scan, record every crossing at which the larger original label is under. Set m equal to the maximum larger label of these crossings, with default m=1 if there are none. Then P_m holds, and no smaller m>=1 has P_m. Thus it finds the maximal initial block of singleton extractions certified by this criterion.

The retained core and top positions can be constructed in O(n+L(w)) time: maintain full labels at current positions, and relative positions of labels1,...,m within their own subsequence. At a both-retained crossing emit its sign times the smaller of their relative positions and swap those positions. At any other crossing emit nothing. Always swap the full labels at the crossing. Reverse the emitted chronological letters at the end. The final full-label permutation gives p_k. This is a literal deletion operation; no braid equivalence test is needed.

## Corollary for far-commuting letters

Suppose w and w' differ by swapping two adjacent distant letters. These two crossing events involve disjoint pairs of labels. Swapping their order preserves the attached label pairs, signs and overstrands, all later attached label pairs, and the final permutation. Therefore P_m and the ordered peeled top positions agree for w and w'.

If at least one of these two crossing events has an endpoint label >m, it disappears in the deletion. The retained words are then literally identical: projecting the two disjoint crossing events onto retained strands either removes both or leaves a single identical event. The theorem gives equal full output. The crossing-label multiset and retained length also agree, so the exact transition-count formula gives equal logical step counts. This conclusion does not require the swapped generator indices to occur only once elsewhere.

For example, in B4 the pair `(3,1,-3)` and `(1,3,-3)` satisfies P_2 and has retained core `(1)`. Both crossings with generator index3 involve the peeled label4 passing over; their signs differ because its side changes. The swapped index3 occurs twice, so this example lies outside the single-use hypothesis. The theorem applies to arbitrarily complicated retained cores in the same way.

If both events involve only retained labels, their projections are still disjoint adjacent crossings in the retained strand order; thus they are distant generators there. The theorem reduces the equality question to their smaller-rank projected pair. It does not settle that pair. In particular, uniqueness of the original positional generator indices does not itself establish P_m or prove the general conjecture.

## Source audit

- `src/zipper_algebra/zipping_algorithm/algorithm.py`: `_single`, `_single_step`, `_extract_step`, `_lift`.
- `src/zipper_algebra/zipping_algorithm/extraction.py`: `_over_route`, `split_layered`.
- `docs/zipper-block-separation.md`: reconstruction of frozen extraction contexts from ordered original ports.

Independent reviewer checked the coefficient-one prefix, top-label relabeling, injectivity on collection keys, and the exact transition formula. No source code changes are proposed here.

## Exact checks with a branching core

In B4 consider `(1,-3,1,3,1,-3,2,1)` and
`(-3,1,1,3,1,-3,2,1)`. Both swapped indices occur repeatedly.
Bottom strand 4 always passes over and ends at top position 3. Deleting it
leaves the same B3 core `(1,1,1,2,1)`, whose calculation takes 637 logical
transitions and has 35 terminal terms. Both ambient calculations take 647
transitions and have identical full 35-term outputs. The surviving top-label
map is `(1,2,4)`, so reconstructing the routes is necessary.

Six pairs with branching cores at ranks 4, 5, and 6, together with their six
cores, gave 18 exact certificates. All passed independent C++ replay; full
output equality and the transition formula were checked. These finite checks
support the source audit; the general conclusion follows from the proof above.

## Size of the first-peel family

Exactly 8,073,424,591 of the 30,346,444,800 full-support signed B8 words of
length 10 have bottom strand 8 always over. This is about 26.6 percent.
Their core has lower rank before any nontrivial zipper is created. This count
is not a runtime saving estimate: these words may be relatively cheap.

A finite dynamic program proves the count. Its state is `(p,S)`, where p is
the current position of bottom strand n and S is the set of generator indices
used. Start with `(n,empty)` and weight 1. At each step choose i in
`{1,...,n-1}` and add i to S. If i is not adjacent to position p, both signs
are allowed and contribute weight 2 without changing p. If i is adjacent,
exactly one sign puts the distinguished strand over, contributing weight 1
and moving p across that crossing. Sum the weights with full S after L steps.
Other strand labels are unnecessary for this single-strand condition.
The census runner does not yet use this theorem as an input reduction.
