# Separation of disconnected strand blocks

**Reference-transition theorem.** The existing procedure separates disconnected
strand blocks: shuffling their letters while preserving each block's literal
subsequence preserves the full routed terminal formal sum. The proof uses
[reachable packet contiguity P](zipper-packet-contiguity.md) and
[matched catalog equivariance C](zipper-catalog-equivariance.md), together with
the controller, normalization, and routing arguments below. It supplies no
quotient by far commutations inside a connected block.

The conclusion concerns the existing procedure on the original
literal input. It preserves the order of letters inside every input block.
It does not depend on first commuting the input into consecutive factors.

The source-level case analysis refers to the Python reference controller in
`algorithm.py`, its upward/framed controllers, and `extraction.py`. Applying
the same statement to a native implementation additionally requires that
implementation to realize these transitions; finite parity tests alone do not
prove correspondence for all reachable states.

The geometric projection is:

```text
original endpoint colors       A A A | B B B B | C C
                                  \   B-only   /
D_B: delete other colors          B B B B
intrinsic endpoint numbers        1 2 3 4
original endpoint labels          4 5 6 7
```

The vertical bars are cuts with no crossing between colors. The intrinsic
numbers change on deletion; original endpoint labels supply the reconstruction
map. Frozen overroutes are handled separately from the active diagram.

## Reachable packet interface

For a width-t child at anchor s, read the incoming word from the right. Its
literal tokens are internal single letters I, far single letters F, and the
two complete band crossings

```text
I:  ±l, s <= l < s+t-1          F: ±l, l < s-1 or l >= s+t
L:  ±(s+t-2,...,s-1)           R: ±(s,...,s+t-1)
        anchor s -> s-1                 anchor s -> s+1

... | F | I | R | F | L |  -> active child
            whole tokens, with a common sign inside L or R
```

Width one permits arbitrary signed input words because all its tokens are
single letters. The needed packet lemma says the existing X2 traversal emits
complete width-two parent packets, turnaround retains their cuts, and growth
uses those parent packets as the next child's packets. Reflection transfers
the scanning convention. Consequently a foreign-color letter can lie between
required packets but cannot be interleaved inside a required multi-letter
packet in a reachable higher-slider state. The symbolic all-family proof and
endpoint-anchor induction are in [P](zipper-packet-contiguity.md).
Monochromatic crossings alone would not supply this property.

## Coordinate convention and simulation invariants

Partition the original strands into consecutive color intervals, with no input
crossing between colors. Idle strands can be singleton colors. At any stage,
ignore frozen extraction contexts and decode the active physical word as
`work_expression` does (including `zip_down` and `world_phi`). Give its top
endpoints the original colors of `top_labels`.

For a color C, deletion D_C removes the other strands, preserving order. On a
same-color generator in a block with current offset a it sends signed index j
to the same signed index j-a. Original top and bottom labels are retained as
an ordered label map. Coordinate reflection is compatible with deletion:

    D_C(phi_N(w)) = phi_m(D_C(w)),

where the coloring on the reflected diagram is reflected too and m is the
number of surviving C-strands. For an index j in [a+1,a+m-1], its reflected
index N-j minus the reflected offset N-a-m is m-(j-a), proving this identity
generator by generator; reversal handles concatenations. The identical
calculation sends a slider anchor to the smaller-rank reflected anchor.

The induction uses three assertions, established by the controller cases below:

1. Every decoded active crossing joins equal colors, and the active coupon is
   monochromatic. The top and bottom color lists agree with their original
   label maps. Colors occupy consecutive intervals, because a same-color
   transposition leaves the color sequence fixed.
2. Projection gives the standalone active-color pass up to identity steps;
   corresponding children have identical Laurent edge coefficients, intrinsic
   frames/payloads, and boundary factors. Every foreign color's literal word
   projection is unchanged.
3. Every marked tape and layer side has the reachable packet grammar P, with
   cuts retained on deleting foreign singleton packets.

Assertion 3 is lemma P. The following source case analysis proves assertions
1 and 2 inductively from P and catalog covariance C. Merely
preserving monochromatic crossings would not imply packet contiguity.

## Lemma C: catalog covariance and foreign crossings

The full structural proof and its executable exact polynomial audit are
linked in [catalog equivariance C](zipper-catalog-equivariance.md).

Let `T_a` translate a word on m strands into positions `[a+1,a+m]` of
`N=a+m+b` strands. It preserves generator signs, coupon sizes, and orientation
flags. Reflected coordinates obey

    R_N T_b = T_a R_m,

so normal I/Mphi rule parameters `j,l` increase by a, whereas the base
parameters of M/phi/BM rules increase by b. Using the same offset for both
orientations would be incorrect when a and b differ.

The structural audit of the frozen catalog
`27c5dd71d09718a5561d18c802f714e07abac87890dc87169215369541462470`
checks all coordinate, dimension, sign, and coefficient expressions as exact
integer polynomials. Coordinate expressions have derivatives `(0,1)` in
`(partial_n, partial_j+partial_l)` for I/Mphi and `(1,-1)` for M/phi; dimension and sign expressions
have zero derivatives, and coefficients depend only on size i. Induction on
the expression constructors gives literal expansion covariance. The 20 BM
identities inherit it under their semilinear transformation, covering all 148
executable rule IDs. Named word bindings are translated as well.

Guards need a separate qualification. Padding can make an ambient-boundary
guard true where its smaller-rank counterpart is false. Such an extra candidate
requires a left or right boundary crossing in its literal left-hand side.
A separated input contains no such crossing, so that candidate does not match.
Thus the invariant object is the ordered list of guard-and-literal-matched
candidates, not the list satisfying guards alone. Far crossings in a different
color have only the far rule available: its single coefficient is one and it
moves that crossing from the incoming to the processed word without changing
local packet geometry. Deletion makes it an identity step. The branch-order
claim uses the controllers' raw RHS expansion order. Generic catalog display
serialization is not translation invariant; the singleton branching families'
comparisons are decided by size/token form before variable coordinates and do
retain their controller order.

## Projection of complete controller state

In raw branch coordinates the reflection parity is `world_phi XOR
(phase == zip_down)`; boundary words use `world_phi` alone. Delete foreign
letters from `left`, `processed`, `tail`, `output_payload`, boundary words,
and boundary-factor words. Preserve each factor's growth-size tag. Translate
the coupon anchor and `initial_slider.j`. Preserve size, friend count, mirror
and policy, frame powers, `completed_sides`, and `world_phi`. Recompute tape
length markers from the projected segments, and recompute the internal branch
phase from whether its projected incoming prefix is empty. The outer algorithm
phase is unchanged.

Delete foreign far events from `processed_events`. Retained events have their
`n,j,l` parameters translated in their own normal/reflected rule coordinates;
their rule, RHS branch index, size, friend count, sign, and mirror policy are
unchanged. Every reached `_new_zip` and turn initializer starts with friend
count zero, empty events, and `initial_slider.processed_length=0`, so this last
field has an unambiguous projected value. Restrict the endpoint label maps to
the color using its fixed original increasing correspondence. Depth and trace
hashes are outside the state relation.

Completed contexts are frozen. Lemma 3 below determines their route map from
the ordered intrinsic components and original endpoint sets.

## Physical color closure

The union of a size-i coupon and its k friends is the connected interval
`[j,j+i+k-1]` in normal coordinates, or `[j-k,j+i-1]` when mirrored. Each friend
block meets the preceding one and the root. Thus its crossings and coupon
force this entire interval to have one color. The child tapes and their frame
powers are internal to it. A foreign letter is strictly far from the entire
interval, including in reflected coordinates.

By C, an active matched rule has exactly the translated standalone RHS;
its output is contained in that color interval. A foreign step has one unit
far branch. Tape and band transport can move an anchor only across an existing
complete packet, whose crossing strands have the active color. Layer growth
annexes precisely such a strand and places its new frames internally.
Reflections preserve physical colors. Lemma 2 gives order-preserving deletion
at extraction. This simultaneous induction establishes color closure for
all reached states, using P to ensure that the multi-letter windows consumed
by the induction are whole packets.

## Upward and framed scans

Internal child absorption uses a relative interval test, so it projects to
the same absorption and cannot consume a foreign letter. The word stays
literal while its marker moves into `output_payload`.

For an active-color token, P supplies the same atomic packet after deletion.
Frame transport compares indices relative to the marked interval and retains
signs and frame exponents. It commutes with translation. C gives the same
matched rule and RHS branches, with the same Laurent coefficients.
`_branch_shape` and `_successor_control` use those translated expressions, so
anchors, friend counts, frame powers, residual tails, and markers agree. The
payload displacement `before.end-out` is unchanged because both anchors shift
by the same amount. The frame checks therefore agree as literal equalities.

A foreign singleton has only the unit far rule. It moves from the incoming
suffix to the processed prefix and changes no active packet geometry. Its
framed event is deleted by the projection. The projected internal phase may
already be complete when only foreign incoming letters remain; recomputing that
phase as specified above makes each such step an identity. Each step removes
one foreign singleton from the incoming prefix, so there are finitely many
before the next active packet or completion. P excludes a foreign singleton
from the interior of a required active packet.

## Lemma T: transport and actual candidate windows

A foreign letter F at index l is strictly outside the marked size-t interval
at s: `l<s-1` or `l>=s+t`. `frames._local` tests this far case first and
returns the same endpoint s and the one-letter word F on both sides,
independently of frame power, sign, and transport side. `transport_frame`
therefore consumes exactly F and emits exactly F; `check_frame` reconstructs
these same data and accepts. `transport_band` and its checker use the same far
condition and unchanged anchor. Thus an extra foreign step introduces no
frame or band recognition failure.

For a framed foreign step, `before.end=out=out2`, `after.end=inp=inp2`, and the
transformed residual is F. Its core and retained frame words are unchanged,
so the subsequent `_parts` validation reconstructs exactly its stored word.
C provides the admissible far guard.

The framed controller matches candidate windows **after** frame transport.
Translation commutes with this parser: all far/internal interval tests
translate; an internally reflected index `2*s+t-2-l` translates by the same
offset; and complete left/right tape blocks and their endpoints translate.
P supplies a contiguous ambient packet for the local candidate. Conversely,
an ambient successful window cannot hide a foreign letter in its interior:
T retains that letter, while active frame transport cannot turn it into an
active-color crossing or create a crossing across the block boundary. C
excludes added boundary candidates. Every ambient successful window therefore
projects to a local successful window, and every local one has its ambient
copy. A different shorter or longer ambient match would contradict the
standalone candidate's uniqueness. This establishes equality of the actual
guard-plus-transported-literal candidate sets, including the residual-frame
and `_parts` checks. Any active structural failure corresponds under the
projection; added foreign steps have no non-resource failure.

## Turnaround and phase changes

For X2, the processed word is, by P, a stream of complete width-two packets,
permitted internal singleton twists, and foreign far letters. Band transport
is the identity on anchor and orientation for the foreign letters. Removing
them leaves precisely the standalone transport stream; all other operations
translate. The turn retains the processed word literally.

For higher sliders, turnaround folds the recorded events in order. A retained
event has the same catalog result, pre/post anchors, frame increments,
internal reconstruction W, and residual R after projection. A foreign far
event has zero frame increment, unchanged anchor/mirror, empty W, and a
single-letter residual F. It contributes only F to the exterior and changes
no fold variable. Removing it therefore commutes with the fold, including the
final initial-anchor, mirror, frame, policy, and processed-length checks.

The wrapper removes the entire recorded friend fan and shifts the active
payload by the same signed friend count in both ranks. Up-to-down change uses
phi on the untouched tail and new body, swaps input/output frames, and reflects
the anchor. Since phi is an anti-involution,

    phi(phi(tail) X phi(body)) = body X tail

with the corresponding coupon and frame coordinates. The projection identity
therefore proves compatibility with the decoded physical expression. The
down turn uses the same identity and saves the reflected fan at the physical
bottom boundary, with the same growth-size tag. This establishes the claimed
state relation through both phase changes.

## Foreign words are preserved literally

Delete all active-color strands from the decoded physical local expression,
including its boundaries and excluding frozen extraction contexts. Single and
zip far steps transfer the same F from an incoming suffix to a processed
prefix, preserving its foreign subsequence. Active rewrite words, frames,
payloads, and friend fans contain only the active color.

For a framed turn, an event tail `T=C R` is replaced by `W R`; C and W are
active-colored. Its foreign subsequence is unchanged. Event tails and their
events are prepended in the same order, so concatenation preserves the full
foreign word. X2 retains its processed word, and the anti-involution identity
above handles the wrappers without reversing the physical foreign order.
Layer far steps likewise preserve the concatenated foreign word; side changes
leave the decoded expression unchanged. Extraction retains every foreign
crossing with order-preserving reindexing. Thus when a block is exhausted,
the next block receives its original literal subsequence, rather than a new
commutation representative.

## Singleton scan and idle inclusions

The singleton scan selects the largest remaining bottom position. A foreign
crossing is far, so it moves from `left` to the front of `processed` with unit
coefficient. Both the active-color projection and that foreign color's full
word remain unchanged. An active-color crossing has the same relative index,
active-over predicate, sign, and side after deletion. Catalog covariance then
gives the same close/foe children and Laurent coefficients, including the new
size-two zipper when present. A completed scan goes to extraction of the same
single component after projection.

As a special case, embed an m-strand word with a idle strands to its left and
b to its right. Each right idle strand undergoes only a finite singleton far
scan and extraction before any strand of the word is selected. Its complement
is the same translated input. The corresponding left idle scans occur after
the word's intrinsic components have been extracted. Padding introduces no
letters into the active tapes, so it cannot disrupt a packet suffix. All
active position and anchor formulas translate, and reflection exchanges the
two padding offsets. Subject to the common controller and extraction lemmas,
this gives the contiguous-interval inclusion correspondence.

Ambient output routes are reconstructed rather than translated tokenwise.
For example a terminal `X_(2,1)` term in B2 becomes, for a left-padded copy in
B3, the routed word `sigma_1 sigma_2 X_(2,1) sigma_2^-1 sigma_1^-1`.
Its intrinsic component has original endpoint set `{2,3}`. Lemma 3 states the
exact reconstruction relation used in the theorem.

## Lemma 1: layering commutes with color deletion

Consider a reached `layer` state, with the simulation invariants above. Delete all colors except
that of its coupon, and project the coordinate-dependent fields. Then each
layer transition is either an identity on this projected state or exactly the
standalone layer transition with the same scalar multiplier and ordered
children. The foreign-color projections remain unchanged.

Proof. Write the controller's word as `unprocessed payload X right`.
`payload` consists of letters internal to X, so every payload letter belongs
to its color; its stored length is unchanged by deletion. Consider first the
guard `_outside_over(unprocessed)`.

The function scans from the right with membership in X. A crossing of a
foreign color exchanges two false membership entries and contributes no
outside-over event. A crossing of the active color updates the restricted
membership list exactly as in the smaller rank, with the same sign and
predicate. Therefore the guard has the same truth value after deletion.

If the guard is false, the side-change transition applies phi, swaps the two
boundary words, reverses and reflects their factor lists, resets payload
length, toggles `world_phi`, and sets `completed_sides=1`. Each operation
commutes with deletion by the reflection identity above. On the second side,
both states instead enter extraction. Thus premature stopping because of
spectators cannot occur.

If the guard is true, the last unprocessed letter belongs to one of these
cases:

* A foreign-color letter is geometrically far from X, because color blocks
  are consecutive and adjacent crossings inside a different block cannot
  meet X. `control.layer_far` transfers that one letter to the other side of
  X, across the active internal payload. Deleting foreign colors makes the
  before and after states identical. Within that foreign color it does not
  cross another foreign letter, so its own word projection is unchanged.
* An active-color internal letter is absorbed into the payload. Its local
  index and the payload increment coincide after deletion.
* An active-color far letter gives the identical far transfer at the reduced
  index.
* A letter adjacent to the coupon belongs to the boundary of a full packet.
  By P the whole crossing packet is a contiguous suffix, and
  every one of its crossings is monochromatic. In particular its outside
  strand has the coupon's color. Deletion translates all packet indices by
  the same offset and preserves its size, sign, direction, and contiguity.
  The payload shift is the same signed index shift. The nonoffensive branch
  moves the anchor by the same amount. In the offensive case the continuation
  and growth coefficients depend only on size and direction, so they agree;
  the new mirror bit, frame powers, and size+1 coupon agree as well. The
  growth packet lies in the same color, and the untouched prefix and opposite
  word have unchanged foreign projections. The separate packet-growth lemma
  supplies the zipper grammar required by the recursive growth branch.

These are all branches of `_layer_step`. Between active transitions only
finitely many foreign far steps occur, since each removes the final foreign
letter from the finite unprocessed word. No transition creates a foreign
letter. This proves the claimed stuttering correspondence for layering. ∎

## Lemma 2: extraction and recursive restart commute with deletion

Suppose a layer-complete state satisfies the color invariant. Its selected
component is monochromatic. The extracted intrinsic component word, terminal
normalization data, original endpoint sets, and remaining word of that color
are exactly those obtained by standalone extraction. Each other color retains
its previous literal word. The complement again consists of consecutive
color blocks, with their surviving labels in order.

Proof. `_bottom_ports` starts with membership in the coupon and transports it
through the right word. Each transport crosses equal colors, so all selected
bottom endpoints have the coupon's color. Transport through the full word
similarly gives selected top endpoints of that color. In `split_layered`, the
scan has three cases at a crossing: two selected endpoints contribute a
component generator; two unselected endpoints contribute a complement
generator; one selected endpoint contributes no generator to either projection
and is checked to pass over the complement.

Deleting foreign strands first cannot change any of these cases for active
color crossings. The component index is the count of earlier selected strands
plus one, so it is identical whether foreign strands were removed first or
not. For a retained complement crossing, ordered deletion of selected strands
and ordered deletion of foreign strands commute. Foreign crossings all have
two unselected endpoints and retain their order and signed indices within their
own color. The sole coupon becomes the same size-k coupon at index one.
This proves equality of `component_word`, of the selected local ports after
projection, and of all remaining color projections.

The physical `world_phi` correction does not alter the argument by the
coordinate identity. For size one, normalization emits the same empty
component and empty twist vector. For size two, it sums the same component
signed generators and emits the same power and default coupon. For size at
least three, the body-only `inner_split` obeys the same deletion argument.
Thus the whole/inner boundary differences `route_left` and `route_right` are
identical intrinsic words. `_project_factors` retains precisely crossings
with both endpoints selected and numbers them by the number of earlier
selected strands. Ordered deletion therefore gives identical projected
boundary factors, including unchanged growth-size tags and omission of the
same empty factors. All arguments to `normalize_canonical_slider` are now
identical: component word, upper/lower factors, and the component-local
normalization limits. It consequently emits identical twist powers,
orientations, intrinsic top/bottom routes, hierarchy lists, and parities.
No claim that arbitrary equal braids have identical normalization is needed;
the inputs to this deterministic function are literally identical.

The original endpoint lists are obtained by indexing the same selected
endpoints in the ordered original label lists. Restart deletes exactly those
entries, clears both boundary words and factor lists, resets `world_phi`, and
calls `_single` on the complement. Ordered deletion cannot interleave color
blocks. Frozen outer routes are stored only in the new extraction context;
they are not inserted into the next active word. This proves the restart
claim. ∎

## Lemma 3: intrinsic ordered records determine all frozen routes

A terminal context list is determined uniquely by (a) the original ambient
rank, and (b) the ordered list of intrinsic component words and complete
terminal records after removing only `top_ports` and `bottom_ports`. In
particular this retained information includes the original top/bottom port
sets. Component order must be retained; a sorted multiset alone does not
establish this lemma.

Proof. Start with the two lists [1,...,N] of surviving original endpoint
labels. For the next record, its original top port set determines local top
ports uniquely by looking up positions in the current top list; similarly at
the bottom. The current rank is the length of either list and size is the
cardinality of the selected set. The extraction code sets

    top_route = inverse(_over_route(current_rank, local_top_ports)),
    bottom_route = _over_route(current_rank, local_bottom_ports).

`_over_route` is a deterministic stable partition, so these two words are
uniquely determined. The intrinsic component word is supplied. Hence every
field of this context is recovered. Delete its original labels from the two
lists and repeat. Finally `_lift((), contexts)` is deterministic, recovering
the entire serialized terminal word. No assertion that routing words contain
only one color is made: they can cross color blocks, but their role is frozen
and fully recoverable from the endpoint data. ∎

## Strong normalization of the uncollected branch process

Fix the original rank N. At a continuation let m be its remaining rank, i the
active coupon size (one in singleton scanning), and d its phase priority:

| Phase | d |
|---|---:|
| singleton scan or upward scan | 6 |
| upward pass complete | 5 |
| downward scan | 4 |
| downward pass complete | 3 |
| first layering side | 2 |
| second layering side | 1 |
| extraction | 0 |

Let U be the actual unprocessed prefix and k the friend count. Set

\[
\Phi=\begin{cases}
2|U|+\mathbf1_{\text{next crossing offensive}}, & \text{singleton scanning},\\
N|U|+(N-1)k, & \text{an X2 pass},\\
|U|, & \text{a framed pass of size at least three},\\
\text{unprocessed side length}, & \text{layering},\\
0, & \text{extraction}.
\end{cases}
\]

Then `(m,N-i,d,Phi)` strictly decreases lexicographically on every successor
branch. A safe singleton step consumes a letter. Its same-size offensive
correction changes the offensive indicator from one to zero; the scalar branch
consumes a letter. A new zipper increases i. For i at least three, each successor sets its incoming marker to
`len(urest)`, strictly shorter than U; no prefix is returned to the unprocessed
word. Thus Phi decreases without requiring a bound on new frame powers or a
claim about the change in k. Internal-payload absorption also consumes a letter.

The X2 exceptions require the actual successor parser. The fourth branches
of nearest-negative and left-positive return k prefix letters to U and reset
k to zero, after consuming one letter. Thus

    Delta|U| = k-1,   Delta k = -k,
    Delta Phi = N(k-1)-(N-1)k = k-N <= -2.

The valid interval gives `k <= N-2`. In nearest-negative, the returned prefix
is a signed singleton followed by a width-one block starting at j+2, while a
recognized fan would start at j+1. In left-positive, the width-one block starts
at j+1 and the new slider is at j-1, whose fan would start at j. These are
therefore literal returned letters, not recognized friend fans. All other X2 branches consume one letter and increase k by at most one,
giving `Delta Phi <= -N+(N-1)=-1`. Independently, `(|U|+k,|U|)` decreases
lexicographically on every X2 branch; no assertion that |U| alone always
decreases is needed.

Turnaround and pass completion decrease d. Same-size layer steps consume a
letter/packet or absorb an internal letter, decreasing its unprocessed side
length. Side change decreases d. Layer growth increases i and hence decreases
N-i even when it introduces a longer prefix. Extraction decreases m and permits
the other entries to reset. This uses the actual continuation of same-size
layering branches; restarting zipper passes after each correction would require
a different argument.

The tuple has nonnegative integer entries. Every call has finitely many
children and finite auxiliary exact computations. There is consequently no
infinite branch; finite branching then implies a finite uncollected tree.
Collection and cancellation are unnecessary for this termination argument.
This is a well-foundedness proof, not an aggregate running-time bound.

## Theorem: exact ordered product of terminal sums

Let each original color block's standalone input be its original subsequence,
translated to start at generator one. Preserve this subsequence's literal
order; this statement does not identify far-commutation-equivalent words
inside a block. Use the analyzed reference schedule that continues layering
after same-size corrections and restarts zipping after growth, with resource
caps large enough for the finite computations. The existing reachable
packetwise ribbon and terminal-normalization contracts are those of the
[reference procedure](ZIPPING_REFERENCE.md); this is not an assertion that
normalization accepts arbitrary enriched words. Here its active-color inputs
are identical before and after projection, so those contracts transfer
without any additional normalization assumption on spectators.

The preceding lemmas establish state simulation, exact scalars, and finite
uncollected expansion trees. The conclusion concerns full terminal sums;
finite-budget prefixes and resource-limit statuses are not identified.

Then the ambient terminal formal sum is obtained as follows. Choose one
terminal term from each block's standalone terminal sum. Multiply their
Laurent coefficients. Translate every original endpoint label in a block term
by the block's original offset. Concatenate its ordered intrinsic records in
right-to-left block order. Reconstruct all ambient contexts and the output
word using Lemma 3. Sum these terms and collect identical reconstructed
outputs. Idle singleton colors use their unique empty-input term.

Proof. `_single` begins at the largest current bottom position. Since colors
are consecutive intervals in the ordered surviving labels, its color is the
rightmost nonempty color. Every subsequent coupon remains in this color,
by color closure and Lemma 1; extraction removes only this color
by Lemma 2. Restart therefore continues in that color until it is empty,
then proceeds to the next color to its left. This proves block order without
assuming that each grown coupon retains the particular original maximal
endpoint. Within the current color, projection simulation identifies the
ordered transitions, their scalar edges, and the ordered extracted records
with the standalone calculation, up to foreign identity steps. All other
color words remain unchanged. Induct on the number of remaining colors and
apply Lemma 3 to lift the intrinsic output records to the exact ambient
encoding.

Collection does not invalidate this statement. In `interpreter.apply_rule`,
RHS grouping is independent of the predecessor coefficient and each retained
successor coefficient is its product with the grouped local coefficient.
Nonzero Laurent polynomials have nonzero products; the zero case is the formal
linear extension, omitted by collected execution rather than passed as an
invalid branch. Singleton rules, layer scalars, and extraction use the same
homogeneous formulas. The collected wrapper explicitly executes a unit
continuation and scales its successor edges by the incoming coefficient.

Induction on the strict continuation rank therefore assigns each unit state
its finite terminal formal sum. Replacing a weighted state by its successors
preserves this valuation, and collecting equal states preserves it by
linearity. The distributive law gives the product of the block sums, including
negative coefficients and cancellations. The final operation is the pushforward

\[
 \sum_{H_r,\ldots,H_1}
 c_{H_r}\cdots c_{H_1}
 [\operatorname{Route}(H_r\mathbin{\Vert}\cdots\mathbin{\Vert}H_1)],
\]

where each H is an ordered intrinsic history and Route is Lemma 3's
deterministic ambient reconstruction. Route need not be injective: exact
collection is performed after routing, using the complete output word and
contexts. This proves equality of the resulting formal sums. It does not
identify step counts, collected frontier sizes, batch boundaries, caches,
traces, or bounded stopping statuses. ∎

## Consequences and scope

An arbitrary shuffle of disconnected block words has the same reconstructed
ambient terminal sum as any other shuffle preserving each block's literal
subsequence. Hence the connected-support standalone computations determine
all such disconnected inputs. The right-to-left component order and original
endpoint maps are part of this reconstruction. Adding idle strands or
translating a contiguous block is the special case with singleton empty-input
factors, giving the inclusion corollary.

The claim excludes far commutations inside a connected block. Those can alter
the deterministic terminal formal sum and require a different coverage claim.
The reference controller's traces and work counters also differ under foreign
identity steps; the theorem concerns the exact routed formal sum after full
computation. Native implementation correspondence remains a separate
source-validation requirement.

The proof dependencies are the source-level [packet lemma P](zipper-packet-contiguity.md),
[catalog lemma C](zipper-catalog-equivariance.md), the transport lemma T above,
and the existing reachable normalization contracts of the reference procedure.
P includes the X2 exceptional branches, event-tail cuts, endpoint fans, and
growth anchors. The strong-normalization rank makes collection a finite
linear calculation. Native source validation and finite campaigns are separate
evidence about executable implementations, rather than premises of this
reference-transition theorem.
