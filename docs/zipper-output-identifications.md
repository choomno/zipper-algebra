# Identifications of zipper outputs and calculations

The statements below distinguish equality of formal output, correspondence of
calculations, and equality of cost. General trace invariance is false; see the
[exact counterexample](TRACE_EXECUTION_COUNTEREXAMPLE.md).

## Status of the identifications

| Identification | Established conclusion |
|---|---|
| Already layered individual strands, same endpoint permutation | Equal full formal output; exact transition formula below. |
| Far commutations within that family | Equal output and logical transition count. |
| Idle-strand inclusions and translations | See the source-level [separation argument](zipper-block-separation.md) and its stated proof dependencies. |
| Disconnected support blocks | See the ordered product and reconstruction statement in the [separation argument](zipper-block-separation.md). |
| Initial descending strands with the same literal core and peeled top ports | Equal full output; [peeling theorem and exact cost](zipper-descending-peeling.md). |
| Far swap disappearing under that peeling | Equal output and logical cost, even with repeated indices. |
| Arbitrary geometrically separated internal braids | Original partition need not persist: [exact three-strand counterexample](zipper-geometric-layer-persistence.md). |
| Already zipped topmost active slider | Frozen by extraction; [precise reached-state statement](zipper-geometric-layer-persistence.md#preservation-of-an-already-completed-slider). |
| Far swap of indices each used once | Unproved output-invariance candidate; not calculation correspondence. |
| Adjacent inverse cancellation | Unproved output-invariance candidate; transition counts need not agree. |
| Reflection, reversal, crossing-sign inversion | No general termwise output or calculation symmetry; examples below. |

The theorem in the next sections is a source-level proof for the existing
procedure. The Python and C++ singleton/extraction paths were independently
inspected. It needs neither the higher-slider packet lemma nor a general
normal-form uniqueness result.

## Definition and recognition

Let n >= 1 and L >= 0. Label bottom endpoints 1,...,n from left to right. Read the literal word from
right to left, transporting these labels through each crossing. Immediately
before reading signed generator e*i, let a,b be the labels in positions i,i+1.
The overstrand is a if e=+1 and b if e=-1, matching algorithm._single_step's
`active_over = adjacent and ((position == index) == (letter > 0))`.

Require that the overstrand's original bottom label be max(a,b) at every
crossing. Thus the sign must be positive exactly when a>b. Swap a,b and
continue. This is an O(n+L) recognition pass with O(n) stored labels, counting index
operations at unit cost.
Call the resulting top-position list of bottom labels p; it determines the
endpoint permutation. This condition is about the given literal diagram,
not an arbitrary representative of its braid element.

## Theorem

For every word satisfying this condition, the existing algorithm:

1. uses only singleton far/close steps, singleton scan-completion, and
   singleton extraction; no foe branch or size>=2 slider is created;
2. emits one coefficient-one terminal term, with all n components singletons;
3. gives exactly the same full terminal context/word encoding for any two
   such words with the same n and endpoint permutation p;
4. takes exactly

       T(w) = 2*n + sum_crossings (n - max(a,b) + 1)

   logical term transitions, where a,b are the original bottom labels at
   each crossing in the backwards recognition scan.

These claims concern successful execution with resource caps large enough
for the finite calculation. Runtime, trace contents and labelscan metrics
are not identified. Different lengths may have the same output but different
T. The algorithm still retains and scans inverse pairs; this is not an
assertion of upfront free reduction.

## Proof

The initial singleton follows original bottom strand n. Every crossing not
incident with this strand is far. At every incident crossing it is the
larger original bottom label and therefore passes over. Hence `_single_step`
never takes its offensive/foe branch: it merely moves the last `left` letter
to `processed`, updating singleton position at an incident crossing.
Each original crossing is scanned once. When the scan is empty,
`control.single_scan_complete` creates the size-one coupon on this strand.

`split_layered` accepts because the selected strand passes over the whole
complement. Its mixed crossings are deleted from the complement, while all
crossings between unselected strands remain in their original order. The
size-one branch of `_extract_step` emits the empty intrinsic component word
and empty twist vector. It deletes original bottom label n and its top
endpoint, then starts the singleton scan on the complement. Every surviving
crossing still has its larger original bottom label over its smaller label.
This proves the induction step. The active bottom labels are n,n-1,...,1.

At stage k, a crossing between original labels a,b is still present exactly
when max(a,b)<=k. It is scanned once at that stage and is deleted during the
extraction at stage max(a,b). Therefore it contributes exactly
n-max(a,b)+1 singleton transitions. Every stage additionally has one
scan-completion transition and one extraction transition, giving 2*n.
There is always a single pending branch of coefficient one, so collection
and batch scheduling neither remove nor multiply any of these steps.

For output identity, the component extracted at stage k has original bottom
port {k} and original top port {p^{-1}(k)}. These ordered original endpoint
sets determine its local ports among the surviving labels. The intrinsic
component word and every size-one normalization field are fixed. The
extraction overroutes are deterministic stable partitions of those local
ports, and `_lift` deterministically reconstructs the output word from the
ordered contexts. Thus n and p determine every terminal field. No uniqueness
claim for general zipper normal forms is used.

## Restricted trace invariance

A distant commutation exchanges two adjacent crossings on disjoint current
positions. In the backward label scan, the two associated transpositions
commute and act on disjoint label pairs. Thus each crossing keeps its sign
and original bottom-label pair, and the final endpoint permutation is unchanged.
Consequently the recognition condition is invariant on a signed trace class.
Within this family, distant commutation preserves BOTH the full terminal output
and the exact logical transition count above. Calculation traces can still
list independent scans in different orders. This is a genuine restricted
version of trace invariance, rather than a claim about arbitrary words.

## Enumeration

Choose an arbitrary unsigned sequence i_1,...,i_L in {1,...,n-1}. Starting
from bottom labels and scanning right to left, each required sign is uniquely
determined. Conversely, forgetting the signs of a recognized word returns
that sequence. Hence there are exactly (n-1)^L such literal words of length
L. Restricting to full generator support gives

       sum_{j=0}^{n-1} (-1)^j binom(n-1,j) (n-1-j)^L.

Their outputs have at most n! distinct endpoint permutations. Counts by
output permutation are the numbers of length-L walks from the identity in
the adjacent-transposition Cayley graph of S_n. A support bitmask enforces
full support. Attach the edge weight n-max(a,b)+1 to recover exact transition
counts, or a generating polynomial to obtain the full cost distribution.

This family is a fraction 2^-L of the signed words of any fixed unsigned
support condition. It identifies a substantial set of easy outputs but is
not a universal substitute for the remaining words. For n=8,L=10 there are
29,635,200 full-support words in this family, all with directly determined
outputs and step counts.

The proof uses `algorithm.py::_single`, `_single_step`, `_extract_step` at
size one, and `_lift`; `extraction.py::split_layered` and `_over_route`; and
`collected.py::_transition` and `collect_terms`. The corresponding native
paths are `algorithm.cpp::single_step` and the size-one extraction branch,
plus `extraction.cpp::split_layered` and `over_route`.

## Recorded exact checks

All 5,705 unsigned sequences in ranks 1 through 6 and lengths 0 through 5 were
signed by the recognizer, reduced natively, and independently replayed by the
C++ verifier. Every run had one coefficient-one all-singleton output; every
transition count matched the formula. Full serialized outputs were identical
within each of the 272 observed `(n, permutation)` groups. These calculations
check the proof rather than justify its general conclusion.

The adjacent-transposition dynamic program gives 29,635,200 full-support B8
words of length 10 in this family, with 2,780 distinct outputs and transition
counts between 26 and 65. Its total agrees with support inclusion-exclusion.
The dynamic program is defined as follows. A state is `(p,S)`, consisting of
the current bottom-label permutation and the set of generator indices used.
Initially it is `(identity, empty)`, with count one and cost `2n`. Choosing
index i swaps positions i and i+1, adds i to S, and adds edge weight
`n-max(p_i,p_(i+1))+1`. At each length add counts of identical states. Retain
states with full S for the full-support count; p determines the output.
Keeping a polynomial in a cost variable gives the exact cost distribution,
or keeping minima and maxima gives the reported cost range.

This is a family of easy computations. Its word count is not a prediction of
a proportional reduction in total census runtime.

## A larger, unproved far-commutation candidate

A candidate suggested by exact computations is:

> Swapping adjacent distant letters preserves the collected formal output
> when each swapped absolute generator index occurs exactly once in the
> entire input word. Other indices may repeat.

The tested complete terminal sums agree on all 196 such B4 pairs through
length five, all 288 signed squarefree full-support B5 pairs, 240 structured
B4 pairs at lengths six through eight, and 288 additional B5 pairs of length
five. The complete B4 scan reduced all 9,324 words of lengths two through five
and compared 3,940 far-swap pairs; 98 pairs outside the candidate restrictions
had unequal outputs. The equality comparison used the full output words and
ordered component records with exact Laurent coefficients. The surviving
family scans used the native producer; they are bounded evidence, not a
general theorem or a second independently checked census.

This candidate concerns output only. In 102 of 484 tested B4/B5 candidate
pairs, even the histograms of non-far rules differ. For example `(2,-3,-1)`
and `(2,-1,-3)` give equal terminal sums, but one uses a whole-band layering
step, a farthest-negative zipper rule, and an outer-positive rule where the
other uses an additional left-negative zipper rule. Both calculations passed
independent Python and C++ replay. Thus the difference is more than a
reordering of harmless far scans.

The globally squarefree-index family is a natural first proof target and
includes nontrivial branching computations. For example `(1,3,2)` and
`(3,1,2)` produce the same 27-term output, with sliders up to size three, in
211 and 210 transitions. If output invariance were proved, its minimal
full-support length layer would have `2^(2n-3)` signed trace classes for
n >= 2. This count follows from signs and orientations of adjacent pairs in
the generator path; it does not prove the output-invariance candidate.

Simpler restrictions already fail: even all-positive words `(1,3,2,3)` and
`(3,1,2,3)` have different outputs (686 versus 1,446 transitions), and limiting
a swap to a prefix or suffix does not suffice. All-negative B4 words survived
142 tested far swaps through length five as output equalities, but 42 of
these pairs had different transition counts; this is also only a candidate.

## Further proof attempt and broader candidate

The single-use conjecture remains open. A further sampled stress campaign
completed 1,355 distinct pairs with identical exact full outputs and left 16
pairs unresolved under resource limits. It includes B5 and B6 words of lengths
6 through 10. Six selected heavy certificates passed both independent C++ and
Python replay. These counts are for this campaign, deduplicated across its
five databases; they are not claimed disjoint from all earlier searches.

For `i<j-1`, write the conjectured pair as
`U sigma_i^e sigma_j^f V` and `U sigma_j^f sigma_i^e V`.
Both U and V omit i,j and lie in the interval parabolic product
`B_i x B_(j-i) x B_(n-j)`. The missing statement is invariance of the
controller's terminal formal sum under interchanging these two boundary
attachments, with arbitrary permitted internal histories of all three blocks.
The attachments exchange original strands; subsequent crossings need not
preserve the original color blocks. Thus the disconnected-support theorem
does not supply this relative statement. Growing sliders can interact with
both attachments, and rewrites can recreate the supposedly unique indices.
Input uniqueness alone is not an inductive invariant. An independent Opus
proof review also left this relative interchange statement unresolved.

A stronger candidate allows `U C V`, where U,V omit i,j and C contains only
signed i,j letters, possibly repeated. Reorder C while preserving each
index's signed subsequence. This includes the single-use conjecture.
A sampled campaign completed 429 pairs with equal full outputs; one pair
remained resource-limited. Of the equal pairs, 343 had different transition
counts. No general proof has been obtained.

There is a separate proved sufficient condition: the
[descending-prefix peeling theorem](zipper-descending-peeling.md). If at least
one of the swapped crossings disappears when the certified descending bottom
strands are deleted, the retained words agree literally, so full outputs and
logical transition counts agree. This allows repeated indices and arbitrary
branching cores, but does not contain every single-use case. When both swapped
crossings survive, it reduces the question to the retained core.

## Free inverse cancellation: evidence, not a theorem

The engine does not freely reduce raw inputs before scanning. In B2,
`(1,-1)`, `(-1,1)`, and the empty word have identical formal outputs but take
6, 29, and 4 transitions respectively.

A native search inserted every signed adjacent inverse pair at every position
of every B3 and B4 word of lengths zero through three. All 7,162 comparisons
completed and had equal full collected formal outputs, ignoring only output
list order. Component records, routing words, and Laurent coefficients were
included. A larger 12,282-comparison search including B3 base length four
found no discrepancy in the coarser endpoint/twist-vector projection. These
searches were not independently replayed in full; selected basic pairs were.

This suggests output invariance under free reduction. A proof would require
compatibility of inverse-pair elimination with zipper transitions, turnaround,
layering, extraction, and exact branch cancellation. It would not imply equal
calculation costs. The current enumeration retains inverse pairs.

## Ordinary geometric symmetries are not execution symmetries

Let R reflect indices `i -> n-i`, V reverse word order, and B invert every
crossing sign. Each nonidentity combination already has a connected B3
length-two counterexample to transition and terminal-term count invariance:

| Operation | Input | Transformed input | Transitions | Terminal terms |
|---|---|---|---:|---:|
| B | `(-2,-1)` | `(2,1)` | 9 / 45 | 1 / 7 |
| V | `(-1,2)` | `(2,-1)` | 25 / 49 | 3 / 6 |
| R | `(1,-2)` | `(2,-1)` | 23 / 49 | 3 / 6 |
| VB | `(-2,-1)` | `(1,2)` | 9 / 90 | 1 / 12 |
| RB | `(-2,-1)` | `(1,2)` | 9 / 90 | 1 / 12 |
| RV | `(-2,1)` | `(2,-1)` | 19 / 49 | 3 / 6 |
| RVB | `(-2,-1)` | `(2,1)` | 9 / 45 | 1 / 7 |

All these complete computations passed independent C++ replay. Substituting
`q -> q^-1` or multiplying slider terms by units cannot change a transition
count. Different term counts rule out a simple bijective relabeling of
terminal records with nonzero scalar multipliers. They do not disprove the
known algebraic semilinear symmetry, or an identity obtained through further
normalization. The procedure privileges the largest remaining bottom endpoint
and treats its undercrossings offensively; these transformations generally
change that computational choice.

## Exact continuation identifications

The collected scheduler already identifies equal complete continuation states
apart from their scalar coefficient and diagnostic maximum depth. Its unit
continuation is linear in that scalar. The native transition cache also
removes frozen extraction contexts and renames original endpoint labels
locally, restoring them on return. These are exact state-level reuse rules.
They do not identify arbitrary initial words before enough computation has
been performed to reach the matching continuations.
