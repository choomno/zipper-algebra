# Current inverse-cancellation invariance: proof audit

Status: OPEN. No general proof or counterexample obtained in this audit.
The desired theorem is equality of fully collected terminal sums, including
all words and ordered component records, for U s_i^eps s_i^-eps V and UV.
It is stronger than their algebraic equality. Current Artin-relation output
noninvariance already rules out deriving it merely from algebraic soundness.

## Singleton critical pair: what really cancels

Write x=sigma_i and E=X_{2,i}. Depending on the side, either sign may
be offensive to the active singleton. The two
crossings of an inverse pair have the SAME overstrand, so the active strand
can be under at both; they are not necessarily two harmless close steps.
The identities are

 x=(1-q)-E+q x^-1,
 x^-1=(1-q^-1)+q^-1 E+q^-1 x.

Expanding the first offensive crossing and then its reversed-crossing child
at the other crossing gives opposite E contributions, and the coefficients
of the surviving one-crossing singleton branch sum to zero. However the E
contributions are encoded as DIFFERENT zipper startups: one has the internal
crossing adjacent on the unprocessed side, and the other has that crossing
in its exterior tail. Algebraically these are xE and Ex. For X2 they commute,
but equality of completed algorithm outputs from these two control states
has not been proved for arbitrary surrounding U,V.

Even the remaining unit-coefficient singleton branch carries the now-over
inverse pair in its processed tape. Removing that pair requires a second
compatibility statement when future scanning creates/grows a zipper.

Thus the first nontrivial proof obligations already at singleton birth are:

1. Internal-crossing placement covariance for X2 startup: transferring an
   internal crossing between unprocessed tape and exterior tail, across E,
   preserves the final full output in every admissible surrounding context.
2. Processed over-bigon deletion: erasing a pair carried behind an active
   singleton remains valid through all its possible future births/growths.

The existing rank-two theorem settles these without spectators. General
ambient rank introduces the friends, frames, turnaround and layering data.
Neither ordinary algebraic commutation xE=Ex nor matching decoded physical
expressions alone identifies deterministic completions: the controller sees
processed/tail boundaries and uses them in turnaround.

## Nonlocal layering obstruction to a stuttering proof

algorithm.py::_outside_over scans the WHOLE unprocessed side before deciding
whether to change sides/extract or execute another layering step. Its Boolean
is not invariant under free reduction. An exact direct check is

 _outside_over((-1,-2,2,1),n=3,size=2,start=2) == True,
 _outside_over((),n=3,size=2,start=2) == False.

The first word freely reduces to empty. Its two-strand selected band {2,3}
is crossed over by the outside strand, then crossed back; the completed band
has no net braid interaction. The lookahead nevertheless sees an offending
outside crossing. The full word is a pair of inverse whole-band crossings.

For a layer state containing this tape, the original controller can take
far/internal scan steps and then a band-growth rewrite, whereas erasing the
tape permits immediate side-change/extraction. A proof must show cancellation
of the extra higher-slider descendants after zipping, turnaround, and internal
normalization. The raw state relations cannot be a simple stuttering equality.

This is a controller-state obstruction, NOT an initial-word counterexample.
I have not proved that this particular displayed layer state is reachable from
an ordinary initial word with all the algorithm's inherited boundary grammar.
Any theorem stated for arbitrary geometrically admissible layer states is
stronger than the original initial-word conjecture. A valid proof may instead
restrict its relation to genuinely reachable typed states and use their extra
invariants. The source-level Boolean calculation itself is exact.

## Sufficient route to a general proof

A relative diamond/bisimulation theorem would suffice: inverse-pair erasure
in a reachable typed state must commute with completion of each catalog and
control phase, possibly through a finite linear combination of intermediate
states whose contributions cancel. It must include:

- singleton birth and the processed-versus-tail placements above;
- internal, far, friend and foe/farthest rules, including pairs split by a
  whole-band packet boundary;
- both reflected passes and frame transport at turnaround;
- the nonlocal layering lookahead and slider-growth case;
- extraction/canonical internal normalization and ordered original endpoints.

Termination plus algebraic soundness alone is insufficient. The missing
statement is equality in the free module of terminal records, before mapping
those records to the zipper algebra. No all-family relative-diamond proof was
completed in this audit, so free reduction must not be enabled as a proved
current-output preprocessing rule.

## Targeted finite sanity checks

Native completed outputs in B3 agree for base (1,-2) and the words
(2,-2,1,-2), (2,1,-1,-2,1,-2), (-1,-2,2,1,1,-2), and
(2,1,-1,-2,2,1,-1,-2,1,-2). They take respectively 23, 67, 236, 124, 685
logical steps and each has 3 full terminal terms. This was a producer-only
check using the full output normalizer, not an independent replay campaign.
It illustrates that cancellation can require substantial branching and
phase changes. It does not prove any general invariant.

## Additional bounded search

An additional campaign found 247 distinct completed inverse-insertion pairs
with equal full outputs and 26 distinct resource-limited pairs; no mismatch
was found. Tests include full-support B5/B6 inputs through insertion length 10
and every inverse insertion into the positive B3 Artin pair. Four selected
higher-rank certificates passed both independent C++ and cold Python replay.
The remaining equality observations used the exact current native producer.
Certificates and source identities were retained in separate research SQLite
databases. This evidence does not enable an arbitrary-rank census quotient.

## Conditional savings

If the conjecture is proved, free reduction followed recursively by support
separation and peeling would leave 27,052,259,653 requests through ranks 2–8
and length 10, plus identity. This is 36.28 percent fewer than the current proved
42,454,828,932 requests, a factor of only 1.569. No connected far commutations
or Artin relations are included. It remains a billions-of-inputs domain.

The exact count uses a right-to-left automaton with state(last signed letter,
used-support mask,position of original bottom n). Reject a new letter if it
is the inverse of the previous letter. Position0 is an absorbing marker after
an undercrossing; otherwise an adjacent overcrossing changes position, and a
far crossing leaves it unchanged. Sum full-support states with position0.
Free reduction may disconnect support and peeling may expose new cancellations,
so repeat the process on the lower-rank cores. Every recursive edge lowers
rank and no operation increases length. This coverage argument remains
conditional on the missing current-output theorem.

All 77 rank/length rows were independently checked against full-support
inclusion-exclusion and small raw-word enumeration with full original-label
tracking. Rank-two pure-sign representatives already give the 19 proved jobs.
At B8 length 10 alone the conditional residual count is 15,646,753,584.
Near the minimum full-support length n-1 there are few repeated occurrences,
which limits the benefit of inverse cancellation at high rank.

The active current-output interface is unchanged. This investigation proves
neither a new arbitrary-rank identification nor a modified-algorithm result.
