# Catalog equivariance under interval embeddings

Scope: frozen catalog SHA256
27c5dd71d09718a5561d18c802f714e07abac87890dc87169215369541462470,
Python catalog.py/algorithm.py/upward.py/framed.py/semilinear.py. This is a
local catalog/controller-selection lemma, conditional on packet contiguity
and the absence of crossings across component boundaries. It does not by
itself prove global separation or its preservation through extraction.

## Coordinates

Put m strands into the consecutive interval [a+1,a+m] of N=a+m+b strands,
a,b>=0. Write T_a(sigma_r^eps)=sigma_(r+a)^eps, and translate slider starts
by a, retaining slider size and all orientation flags. Pure named words
must also be translated. Parameters i,k,eps do not change.

For ordinary I or Mphi catalog records use j_N=j_m+a, l_N=l_m+a.
For M or phi records, whose parameters are the ORIGINAL, unreflected
coordinates, use j_N=j_m+b, l_N=l_m+b. Their actual output coordinates
then translate by a. Equivalently R_N T_b=T_a R_m, with
R_n(sigma_r)=sigma_(n-r), and
R_n(X_i at j)=reflected X_i at n-i-j+2.
BM has the same coordinate convention as M.

For example framed.py computes j_base=N-i-j_actual+2. If actual j is
increased by a, then j_base increases by b, not a. It computes l_base
similarly. This distinction is essential when left and right padding differ.

## Uniform structural audit

[The symbolic audit](../tests/test_zipping_catalog_equivariance.py) evaluates every integer AST as an exact
multivariate polynomial over Z[n,i,j,k,l,eps], then differentiates. The audit checks:
  128 frozen records; 1718 word nodes;
  952 coordinate expressions (start/index);
  1532 size/length/count/sign expressions;
  816 coefficient expressions (integer/q_power).
Each coordinate has (d/dn, d/dj+d/dl) equal to (0,1) for I/Mphi,
and (1,-1) for M/phi. Each dimension/sign expression has both derivatives
zero. Every coefficient depends on i alone (or is constant).
These polynomial identities prove the requisite substitutions for arbitrary
integer parameters, not only tested ranks. In particular coefficients and
loop counts do not depend on ambient padding.

Induction on the word AST now proves literal expansion equivariance:
identity, sigma, slider, fat, block and half_twist use these expressions;
concat preserves order; inverse/reverse/bar commute with translation;
reflections commute using the identity above. Named words require translated
bindings. Expansion preserves insertion order. Translation is injective on
words, so exact collection/cancellation and first-occurrence order commute.
The 20 BM rules follow from their source M records: semilinear transformation
only changes signs, q exponents, and inverse wrappers, and common units
are functions of i. Thus the result covers all 148 executable rule IDs.

## Matched candidates

Raw guard-selected lists are NOT invariant. For example j>=2 may fail at
local j=1 and hold after padding. The correct theorem is equality of
GUARD-AND-LITERAL-LHS-MATCHED candidates for a separated local tape.

The complete structural audit of noninvariant guard atoms found exactly:
  j>=1, j<=n, l>=1, l<n,
  j-1>=1, j-1<n, j<n,
  i+j+k-1<=n, j>=2, i+j<=n, i+j+k<=n.
All other guard atoms are translation invariant; Boolean combinations are
unchanged once these boundary guards are dealt with. Basic interval guards
hold in both environments for an actual local packet and local crossing.
The remaining single/zip guards express availability of an adjacent strand
outside a packet: single left/right and zip left/outer. Padding can only
add such availability. When local availability fails, the added LHS contains
a generator across the active block boundary, absent from a separated tape.
Thus these added records cannot literally match. Conversely every locally
matched record remains admissible and matches after embedding. The exclusion
of flype at i=2 is unchanged. The literal occurrence lengths and offsets in
the local tape are unchanged. Reflection transfers the same argument to M,
phi, Mphi, and BM. Layer boundary guards require the separate layering lemma.

For avoidance of doubt, the zip outer LHS consumes a fat crossing ending at
j+k+i-1; at the right local boundary this is precisely the forbidden boundary
generator. The zip left LHS consumes sigma_(j-1), the corresponding left one.

## Foreign crossings

For an ordinary packet F_(i,j,k) X_i at interval [j,j+i+k-1], a crossing
strictly within a different separated block has l<j-1 or l>=j+i+k.
The far LHS is sigma_l^eps F X; its RHS is F X sigma_l^eps, a single term
with coefficient exactly 1 and no q power. Both crossing signs have records.
All other normal zip LHS prefixes have their last unprocessed crossing at:
  close: j<l<j+k;
  internal/flype/nearest: j;
  left: j-1;
  farthest: j+k;
  outer: j+k+i-1.
Consequently none can match a foreign final unprocessed crossing. The far
rule of the actual sign is unique, preserves i,j,k and all local packet
geometry, and merely moves that foreign crossing to the processed tape.
Reflection and semilinear bar give exactly the same conclusion in the other
orientations/policies. For a singleton the same proof uses its support [j,j]
and the guard l<j-1 or l>j; algorithm._single_step chooses far directly.
Deleting the foreign crossing turns this step into an identity/stutter on
the projected local calculation. This requires the packet-contiguity lemma;
no claim is made for a packet crossing multiple component boundaries.

## Controller ordering

Word expansion and raw RHS term collection preserve order under translation.
This is the ordering used by the zip/upward/framed rule application.
Generic catalog.serialize sorts by repr((q,word)); that DISPLAY ORDER is NOT
translation invariant across decimal digit boundaries. Concrete example:
zip.left.positive.i2.I with (n,i,j,k,l)=(9,2,2,0,3), left padding 8 and
right padding 8 changes serialized RHS order. This does not change the formal
polynomial, and a theorem about generic serialized term order would be false.

The single controller does call instantiate/serialize. Its only branching
families are the four single.foe.*.I records. Direct inspection shows their
serialized comparisons are decided before any varying coordinate: first by
q-power, then by slider versus crossing first-token syntax, or by slider size
1 versus 2. Thus their serialized term order, and the first-occurrence order
of grouped words used by _single_step, ARE invariant under all padding.
Single close/far have one unit branch. Hence this serialization caveat does
not invalidate the local single/zip controller order lemma.

## Scope

A global theorem still needs packet contiguity, preservation of separated
blocks through frame transports/turnaround/layering/extraction, and a
scheduler/collection argument. Global trace parameters, labels, hashes,
iteration counts, and resource-budget outcomes are not asserted identical:
foreign far steps are projected away, and global interleavings may differ.
The result supplies exact local coefficients, literal rewrite words,
matched-rule choices, and local branch order after projection.

## Mirrored packets and retained canonical frames

The second symbolic test compares all 20 `.I` zipper records with their `.M`
records, including all 60 LHS/RHS term pairs, guards, coefficient lists, and
term order. It verifies the following transformation after exact polynomial
normalization of coordinate expressions:

- `sigma(r,e)` becomes `sigma(n-r,e)`.
- A slider of size `h` at `s` moves to `n-s-h+2`, toggling normal/mirrored
  orientation and retaining false reversed/barred flags.
- `fat(L,s,e)` becomes `inverse(fat(L,n-s-L+1,-e))`.
- A descending `block(L,s,c,e)` becomes
  `inverse(block(L,n-s-L-c+2,c,-e))`.
- A half twist of size `h` at `s` becomes the canonical half twist at
  `n-s-h+2`. Concatenation preserves order; inverse transforms its child.

The fat/block transformations are literal word equalities, including empty
blocks. Half twists are retained **canonical frame atoms**: reflection of the
expanded chosen half-twist word is generally only braid-equal to the canonical
reflected expansion. No literal packet claim relies on identifying these
expanded words. The actual catalog stores the canonical frame, which the
controller expands explicitly, internally to the marked child and separately
from the packet cuts.

For a child of width `t` at `s`, the reflected start is `s*=n-s-t+2`.
A descending incoming left packet becomes the ascending right packet at `s*`,
and conversely. Internal/far singletons remain internal/far, and lengths and
cuts remain unchanged. The normal output/input anchors `(j+k+1,j)` become
`(J-k,J+1)`, with `J=n-i-j+2`, exactly the mirrored controller anchors.

For `BM`, semilinear transport negates crossing signs and Laurent exponents.
Every term has one noninverted rank-`i` slider. The common scalar
`(-1)^(i-1) q^(-i(i-1)/2)` cancels from both sides, so the slider AST remains
unchanged: no barred or reversed flag is introduced. The test explicitly
checks slider-node equality, sign inversion, and preservation of the remaining
crossing/block fields for every mirrored zipper term. The integer-expression
parser rejects Boolean values, unknown variables/operators, and malformed
arities rather than silently interpreting them.

The identity `B(Delta)=inverse(Delta)` again uses braid equality; the inverse
canonical frame is retained explicitly. Outside frames, `BM` is literal sign
inversion of `M`, preserving packet supports and cuts. The only nonzero base
frame increments are:

| Base family | RHS indices, starting at zero | Output/input increments |
| --- | --- | --- |
| flype positive | 0 | `(-1,+1)` |
| flype negative | 0 | `(+1,-1)` |
| farthest positive, either size case | 0, 1 | `(-1,+1)` |

`M` preserves these increments at reflected anchors; `BM` negates both.
Existing frame powers are increased by the increments. All base branches
producing a mirrored coupon have zero friends. Consequently their `BM` images
have normal coupons with no friends; no case requires negative normal friends.
Other `BM` outputs have the negative mirrored friends stored by the controller.
For size two, width-one frames expand to empty words. These observations
transfer the packet and anchor calculations to `M/BM`; they do not replace
the separate reachability and global scheduling arguments.
