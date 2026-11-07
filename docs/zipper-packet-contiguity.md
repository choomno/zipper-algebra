# Literal packet contiguity in the zipper procedure

This proof concerns the fixed rewrite catalog and the literal scheduler in
[`algorithm.py`](../src/zipper_algebra/zipping_algorithm/algorithm.py).
It establishes the packet invariant needed when comparing a run with one in
which disconnected spectator crossings have been deleted. It does **not** by
itself establish that comparison, a representative-counting formula, equality
of resource-limit outcomes, or equality of traces.

The statement is about literal words, including unreduced words. No commuting
or braid simplification is used to obtain the packet boundaries below.

## Statement and notation

For a marked band of width $t\geq1$, beginning at strand $s$, define an
**incoming packet stream** by reading a literal word from right to left. At the
current anchor $s$, its next token is one of:

| Token | Literal word, all displayed signs equal to $e\in\{\pm1\}$ | New anchor |
|---|---|---|
| Internal $I$ | $\sigma_l^e$, $s\leq l<s+t-1$ | $s$ |
| Far $F$ | $\sigma_l^e$, $l<s-1$ or $l\geq s+t$ | $s$ |
| Left $L$ | $\sigma_{s+t-2}^e\cdots\sigma_{s-1}^e$ | $s-1$ |
| Right $R$ | $\sigma_s^e\cdots\sigma_{s+t-1}^e$ | $s+1$ |

Each complete $L/R$ token has exactly $t$ letters. For $t=1$, every valid
literal braid word is a stream of singleton tokens. An **outgoing stream** is
read from left to right: an ascending block beginning at $s-1$ moves left,
and a descending block beginning at $s+t-1$ moves right. Reading the same
stream in the opposite direction reverses its sequence of anchor changes.
Write these two grammars as $\mathcal C_{\leftarrow}(t,s)$ and
$\mathcal C_{\rightarrow}(t,s)$. Internal tokens are allowed in both.

**Packet-contiguity lemma.** In every state reached from literal input by the
fixed scheduler, an active $X_i$ pass with $i\geq3$ has an incoming exterior
word in $\mathcal C_{\leftarrow}(i-1,s)$, where $s$ is its recorded output
child anchor. Its untouched opposite word has the outgoing child grammar at
the pass's *initial* input anchor. Every framed rewrite consumes precisely one
incoming token. Completed traversals and layering preserve the corresponding
whole-parent streams, and growth installs them as the child streams of the
next slider. Consequently a required multi-letter tape block is never obtained
by joining pieces across two token boundaries.

The full simultaneous invariant, including the two layer sides, is given
below. The lemma is conditional only on reachability by the specified
procedure, not on a chosen endpoint partition. It does not apply to arbitrary
syntactically valid framed states.

For compact formulas set

$$
 A_h(a,e)=\sigma_a^e\cdots\sigma_{a+h-1}^e,
 \qquad
 D_h(a,k)=A_h(a+k-1,+1)\cdots A_h(a,+1).
$$

A superscript $-1$ on a word denotes its braid inverse, including reversal of
letter order. A zero length or zero count denotes the empty word.

## Incoming dispatch for $i\geq3$

In normal orientation put $t=i-1$ and $s=j+k+1$. The input child anchor is
$j$. The following table exhausts all incoming tokens and both signs:

| Token/support | Rule family |
|---|---|
| $I$ | `control.absorb_internal_tape` |
| $F$, $l<j-1$ or $l\geq s+t$ | `zip.far` |
| $F$, $l=j-1$ | `zip.left` |
| $F$, $l=j$, necessarily $k\geq1$ | `zip.nearest` |
| $F$, $j<l<j+k$, necessarily $k\geq2$ | `zip.close` |
| $L$, $k=0$ | `zip.flype` |
| $L$, $k\geq1$ | `zip.farthest` |
| $R$ | `zip.outer` |

Internal absorption occurs before catalog matching. A one-letter candidate
cannot consume just the boundary letter of an $L/R$ token: the marked-tape
transport requires all $t\geq2$ letters. Conversely each full-block catalog
prefix ends at one of the two boundary indices $s-1,s+t-1$, so it cannot
match a far singleton. Thus the table gives exactly one matching family.

All RHS branches have the anchor change of their consumed token. The far,
close, left, and nearest families leave $s$ fixed. For example the moving
nearest branch has $j'=j+1,k'=k-1$, hence $j'+k'+1=s$. Farthest moves the
child to $s-1$: its mirrored branches have $j'=j+k,k'=0$, whose mirrored
output anchor is $j'=s-1$. Outer moves it to $s+1$: its moving normal
branches have $j'=j+k+1,k'=0$, with output anchor $j'+1=s+1$. Flype has
the same left move with $k=0$. The other branches change $k$ by the
corresponding unit.

Only a far step or one complete tape step occurs in accepted incoming frame
transport. These steps leave their literal words unchanged for every integer
frame power. In particular transporting an arbitrary protected child payload
requires only the net anchor shift; no internal conjugation is omitted.

The $X_2$ exception is treated separately: at $k=0$, `internal.i2` has
priority and `upward.step` excludes flype.

## Complete emitted parent tails

Here $i\geq3$, and the table is for normal `.I` rules. Strip the retained
output frame, friend/slider core, and retained input frame. For a branch that
changes frame powers, further strip its initial outgoing flype block $C$.
The remaining word is $R$. Every listed $R$, read from the left as a
width-$i$ parent stream, takes the post-slider anchor back to the pre-slider
anchor $j$. Branch numbers are zero-based catalog order.

| Family/branch | Post-anchor | Remaining $R$ |
|---|---:|---|
| far $e$ | $j$ | $\sigma_l^e$ |
| close $e$ | $j$ | $\sigma_{l+i-1}^e$ |
| nearest positive | $j+1$ | $A_i(j,+1)$ |
| left negative | $j-1$ | $A_i(j-1,+1)^{-1}$ |
| nearest negative 0 | $j$ | $A_i(j,+1)^{-1}A_i(j,-1)$ |
| nearest negative 1 | $j$ | empty |
| nearest negative 2 | $j+1$ | $A_i(j,-1)$ |
| left positive 0 | $j$ | empty |
| left positive 1 | $j$ | $A_i(j-1,+1)A_i(j-1,-1)^{-1}$ |
| left positive 2 | $j-1$ | $A_i(j-1,-1)^{-1}$ |
| flype $e$ | $j$ | empty; $C=A_{i-1}(j,e)$ |
| farthest negative | $j$ | empty |
| outer positive | $j$ | empty |
| farthest positive 0 | $j+k$ | $D_i(j+1,k-1)A_{i+k-1}(j,-1)$ |
| farthest positive 1 | $j+k$ | $A_{k-1}(j,-1)D_i(j,k)$ |
| farthest positive 2 | $j$ | empty |
| outer negative 0 | $j+k+1$ | $D_i(j+1,k)A_{i+k}(j,-1)$ |
| outer negative 1 | $j+k+1$ | $A_k(j,-1)D_i(j,k+1)$ |
| outer negative 2 | $j$ | empty |

For farthest-positive branches 0 and 1, the stripped block is
$C=A_{i-1}(j+k,+1)$. A flype packet is restored at turnaround as
$W=A_{i-1}(j',-e)^{-1}$ at its relevant parent anchor $j'$. This is an
internal word of the width-$i$ parent, followed by the displayed $R$.

These are literal decompositions for arbitrary admissible parameters:

- Far and close tails are far parent singletons. The nearest/left tails are
  one full parent crossing or a pair returning to the initial anchor.
- In farthest-positive branch 0, $D_i(j+1,k-1)$ moves left $k-1$ times
  from $j+k$ to $j+1$. The following long word is $A_i(j,-1)$, moving
  to $j$, followed by $k-1$ far letters of indices $j+i,\ldots,j+i+k-2$.
- In farthest-positive branch 1, the first $k-1$ letters are far left of
  anchor $j+k$; $D_i(j,k)$ then moves left $k$ times.
- The outer-negative formulas are the same two decompositions with one extra
  whole-parent crossing and one extra far letter. They work also at $k=0$.

The definitions of `block` and `fat` in the catalog give exactly these words;
the assertions are not inferred from bounded-rank tests.

## The separate input-frame residual lemma

The previous parent-tail calculation is not sufficient for
[`framed.step`](../src/zipper_algebra/zipping_algorithm/framed.py). It also
transports the *unstripped* residual $T$, including any outgoing flype block
$C$, at child width $t=i-1$. It reads $T$ **from right to left**, starting
at the old input anchor $j$, and must end at the new input anchor with the
literal word unchanged for every stored integer input power.

Five elementary parsing identities suffice:

1. $A_{t+1}(a,e)$, read at $a$, has one far letter of index $a+t$, then
   one right packet; it ends at $a+1$.
2. $A_{t+1}(a,e)^{-1}$, read at $a+1$, has one left packet, then one far
   letter of index $a+t$; it ends at $a$.
3. $A_{t+r}(a,e)$, read at $a$, has $r$ far letters followed by one
   right packet, for every $r\geq0$.
4. $D_{t+1}(a,k)$, read at $a$, repeats identity 1 at anchors
   $a,a+1,\ldots,a+k-1$, ending at $a+k$.
5. $A_r(a,e)$, read at $a+r+1$, consists entirely of far letters.

None has an internal step. Thus frame parity never changes the emitted word.
For completeness the endpoint calculations for every family are:

| Family/branch | Right-to-left traversal of $T$, starting at $j$ |
|---|---|
| far, close | one far letter, end $j$ |
| nearest positive | identity 1, end $j+1$ |
| left negative | identity 2, end $j-1$ |
| nearest negative 0 | identity 1 then 2, $j\to j+1\to j$ |
| nearest negative 1 | empty, end $j$ |
| nearest negative 2 | identity 1, end $j+1$ |
| left positive 0 | empty, end $j$ |
| left positive 1 | identity 2 then 1, $j\to j-1\to j$ |
| left positive 2 | identity 2, end $j-1$ |
| flype $e$ | $T=A_t(j,e)$, one right packet, end $j+1$ |
| farthest negative, outer positive | empty, end $j$ |
| farthest positive 0 | $kF,R$ from the long $A$, then $k-1$ uses of identity 1, then $C=R$; end $j+k+1$ |
| farthest positive 1 | $k$ uses of identity 1, then $k-1$ far letters, then $C=R$; end $j+k+1$ |
| farthest positive 2 | empty, end $j$ |
| outer negative 0 | $(k+1)F,R$, then $k$ uses of identity 1; end $j+k+1$ |
| outer negative 1 | $k+1$ uses of identity 1, then $k$ far letters; end $j+k+1$ |
| outer negative 2 | empty, end $j$ |

The endpoint is precisely $j'$ for a normal successor and $j'+1$ for a
mirrored successor, independently of its friend count. Therefore this runtime
frame check follows from the catalog formulas for every admissible parameter
and every integer input power; its success is not assumed from reachability.

## Mirror and bar-mirror variants

The exact symbolic audit and coefficient-equivariance proof are in
[Catalog equivariance](zipper-catalog-equivariance.md). Reflection maps
$\sigma_r^e\mapsto\sigma_{n-r}^e$ without reversing product order, and a
width-$t$ interval anchor $s$ to $s^*=n-s-t+2$. It sends each $I/F$
token to the same type and interchanges entire $L/R$ tokens, preserving
literal cuts. For a reflected slider at $J=n-i-j+2$, its anchors
$(J-k,J+1)$ are precisely the reflected normal child anchors.

The actual `.M` records agree term-by-term with this transformation, retaining
half twists as canonical frame atoms. Reflection of a chosen expanded
half-twist word need not be literally the chosen word at the reflected
interval: the equality is a braid identity. This does not change any packet
cut, since these frames remain separate internal words.

The normalized semilinear map in
[`semilinear.py`](../src/zipper_algebra/zipping_algorithm/semilinear.py) reverses
crossing signs and Laurent exponents, sends a canonical half twist to its
inverse, and leaves the slider AST unchanged after removing its common
scalar factor from both sides. It introduces neither a barred nor a reversed
slider flag. Outside the retained frame atoms, `.BM` is literal sign reversal
of `.M`, so all packet decompositions and input-frame calculations transfer.

The only nonzero base frame increments are:

| Family/branch | Output increment | Input increment |
|---|---:|---:|
| flype positive 0 | $-1$ | $+1$ |
| flype negative 0 | $+1$ | $-1$ |
| farthest positive 0, 1 | $-1$ | $+1$ |

`.M` keeps these powers; `.BM` negates both. All base mirrored outputs have
zero friends. Consequently a BM branch that becomes normal never needs a
negative normal friend fan. Base normal outputs become mirrored, with the
negative mirrored fan constructed by the controller. These observations
verify the actual retained-core shapes, not just the support of their words.
They discharge the reflection hypotheses of the input-frame lemma, which
therefore holds for `.M` and `.BM` as well as `.I`.

## The $X_2$ base case and strong normalization

The size-two controller uses singleton child tokens. It excludes flype and
uses `internal.i2` at $k=0$; normalizing mirrored $X_2$ changes no literal
braid letters, since the two defining polynomials coincide.

All ordinary outgoing tails have the preceding parent-tail decomposition at
$i=2$. A mirrored branch has empty width-one half twists, and its initial
$C=A_1(j',+1)$ is an internal singleton of the width-two parent. The
`internal_i2` transport case consumes it before the remaining complete parent
packets. The fourth nearest-negative, farthest-positive, and outer-negative
branches have empty outgoing tails. The fourth left-positive branch, at
post-anchor $j-1$, emits $A_2(j-1,-1)^{-1}$, one full rightward parent
packet returning to $j$. This table proves restoration of the pass's initial
anchor in the X2 turnaround; that path has no separate runtime initial-anchor
check.

The controller recognizes a pre-slider block as a friend fan only if its
length is one and its start is the new slider start plus one. All unmatched
prefix letters are returned to the unprocessed word. The controller computes
`processed_length` only from the nodes after the slider, so none of these
returned prefix letters enters the processed-tail stream. Since the child has
width one, these returned words remain valid streams without any hypothesis
on their internal order.

Here is also an explicit local termination measure. Let $u=|U|$, with
ambient strand count $n\geq2$, and let $k$ be the recorded friend count.
Each step consumes one old $U$ letter. Let $e$ be the number of RHS letters
returned to $U$, and $k'$ the recognized new friend count:

| Family | $(e,k')$ in catalog branch order |
|---|---|
| far, close, either sign | $(0,k)$ |
| nearest positive | $(0,k-1)$ |
| left negative | $(0,k+1)$ |
| internal.i2, either sign | $(0,0)$, with $k=0$ |
| nearest negative.i2 | $(0,k),(0,k),(0,k-1),(k,0)$ |
| left positive.i2 | $(0,k),(0,k),(0,k+1),(k,0)$ |
| farthest negative | $(0,k-1)$ |
| farthest positive.i2 | $(0,0),(0,0),(0,k-1),(0,k)$ |
| outer positive | $(0,k+1)$ |
| outer negative.i2 | $(0,0),(0,0),(0,k+1),(0,k)$ |

For the fourth nearest-negative branch, the leading singleton plus its
nonmatching count-$k-1$ block return $k$ letters. For the fourth
left-positive branch, the nonmatching count-$k$ block returns $k$ letters.
These are the only exceptions to $e=0$. Thus

$$
 M=n u+(n-1)k
$$

strictly decreases: ordinary branches change it by at most $-1$; exceptional
branches change it by $k-n\leq-2$, since $k\leq n-2$. Every uncollected
X2 pass path has at most its initial $M$ rewrites. The RHS branching is
finite, so the whole uncollected pass tree is finite. This is a path bound,
not a bound on the total number of terms of the full algorithm.

## Simultaneous invariant through turnaround and layering

The four word clauses needed for induction are:

1. **Zip incoming:** $U\in\mathcal C_{\leftarrow}(i-1,\operatorname{out}(c))$.
   Protected output payload is internal at that anchor and is separate from
   $U$.
2. **Zip opposite:** $T\in\mathcal C_{\rightarrow}(i-1,\operatorname{inp}(c_0))$, where $c_0$ is the pass's initial slider data.
   $T$ is unchanged during the pass. The processed event stream lies between
   the current input frame and $T$; the current input anchor need not equal
   the initial one.
3. **Layer left:** $L=UP\in\mathcal C_{\leftarrow}(m,a)$, where $P$ is
   internal at the current parent anchor $a$, $U$ has the same incoming
   grammar, and `payload_length` is exactly $|P|$.
4. **Layer right:** $R\in\mathcal C_{\rightarrow}(m,a)$ at that same current
   parent anchor. This holds throughout layering, not merely at entry.

Processed tails and their events are prepended in the same order during a
pass. Their cuts therefore match exactly. Turnaround reads the events in
reverse chronological order, removes each complete outgoing flype block
$C$, and appends the entire $WR$ to the exterior. Here $W$ is internal
to the current parent and $R$ is its full packet stream from the tail table.
Their anchors compose. The retained frames and protected payload are internal
parent words. The complete recorded friend fan is removed as a prefix before
the next terminator is installed; moving the payload through that fan is a
uniform index shift. The fan is not assumed to be a whole-parent packet.

At the first turnaround, the untouched opposite word becomes the incoming
word for the down pass under $\phi$, which reverses product order and
reflects indices. The produced up-body becomes the untouched opposite tail of
the down pass. The second turnaround produces the other layer side. Reversal
and reflection preserve token cuts and transfer incoming to outgoing
conventions at the reflected anchors.

For a layer step:

- An internal singleton is removed from $U$ and prepended to $P$.
- A far singleton $f$ is removed from $U$, and the right word becomes
  $fR$; this is a far first token at the same anchor.
- A complete crossing $B$ moves the parent from $a$ to $a+d$. The shifted
  payload is internal at $a+d$. Reading $B$ from the left at $a+d$ moves
  back to $a$, so $BR\in\mathcal C_{\rightarrow}(m,a+d)$. Uniform sign
  reversal of $B$ for the nongrowth branch preserves this statement.
- A side change sends $UPXR$ to
  $\phi(R)\phi(X)\phi(P)\phi(U)$. At the reflected anchor, $\phi(R)$ is an
  incoming stream, and $\phi(P)\phi(U)$ is an outgoing stream. Resetting
  `payload_length` to zero does not discard its literal internal letters.
  The argument does not require the old unscanned $U$ to be empty.

A growth branch removes the entire offensive width-$m$ packet. The remaining
prefix is an incoming old-child stream at $a+d$, and the unchanged opposite
word is an outgoing old-child stream at $a$. The new $X_{m+1}$ has exactly
these anchors:

| Growth | New slider | Output child anchor | Input child anchor |
|---|---|---:|---:|
| right | normal, start $a$ | $a+1$ | $a$ |
| left | mirrored, start $a-1$ | $a-1$ | $a$ |

The shifted internal payload is installed separately, and any added half
twist is internal to that child. Thus all four word clauses are preserved. Internal tokens permitted in the
opposite/right grammar are consumed by framed or layer absorption when that
side becomes incoming. They do not enter `transport_band` above width two;
its inputs there are the complete parent tails from the emitted-tail table.

## Entry, extraction, and conclusion of the induction

The actual `single.*.I` records contain one size-one or size-two normal slider
per RHS term, with pure braid letters elsewhere. A foe term creates $X_2$
or retains $X_1$, never a higher slider. All initial X2 input words therefore
have the automatic width-one grammar, even for arbitrary interleavings of
spectator letters in the original input.

The contextual interpreter constructs successors literally as
`prefix + instantiated_rhs + suffix`. Collection adds coefficients of
identical words; it never reorders letters within a word. The X2 table produces
whole-parent streams, the two turnarounds give the two layer sides, and growth
supplies the next higher child's streams. This proves the induction on slider
growth.

After extraction, `split_layered` builds the complement by appending only
signed integer generators whose two endpoints are outside the selected
component. Slider tokens enter only the component word. Its complementary
word is therefore a pure braid word, and `_extract_step` restarts the singleton
scan on it. This re-enters the width-one base case. Only this syntactic fact
about extraction is needed here; its reconstruction identity and canonical
component form are separate obligations.

A disconnected spectator crossing can consequently appear between tokens in
a reachable higher-slider input, but cannot split a multi-letter token that
the controller next needs. The initial raw word need not already be packetized
at higher width: its packets are emitted intact during the earlier traversals.
This proves packet contiguity for the existing procedure. The relationship
between two complete runs under spectator deletion requires the additional
catalog, scheduling, and extraction arguments in the separation proof.
