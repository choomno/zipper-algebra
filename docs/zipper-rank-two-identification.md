# Current-algorithm rank-two output identification

Status: proved for the existing reference controller, with sufficient resources.
No input normalization is inserted into the controller. Full terminal output
is compared after exact coefficient collection and ignoring only term-list order.

Let R=Z[q,q^-1], s=sigma_1 and E=X_2=q s^-1+1-q-s. For every literal word w
on {s,s^-1}, let e(w) be its exponent sum. The completed rank-two output is
uniquely determined by e(w). Therefore arbitrary adjacent inverse cancellation,
and arbitrary permutation of the signed letters, preserve its full output.
They do not preserve traces, transition counts or resource-limited status.

## Algebraic uniqueness

There is a direct sum of R-modules

  R[s,s^-1] = R{1,s^-1} direct-sum E R[s,s^-1].

Existence is Laurent division. Cancel highest exponents >=1 using the leading
coefficient -1 of E; these steps introduce no exponent below -1. Then cancel
lowest exponents <=-2 using its lowest coefficient q, a unit; these steps
introduce no exponent above0. Both processes terminate.

For uniqueness, if a nonzero finite Laurent polynomial f has lowest/highest
exponents l,h, then Ef has nonzero extreme exponents l-1,h+1. Their distance
is h-l+2>=2. It therefore cannot equal a polynomial supported on {-1,0}.
Moreover multiplication by E is injective, by the same extreme-coefficient
argument. Consequently the family {1,s^-1} union {s^k E : k in Z} is an
R-basis. This independence is special to the unrestricted rank-two Laurent
algebra; it is not an independence assertion for general terminal sliders.

As an optional closed formula, the remainder of s^e is A_e+B_e s^-1 with

  B_e = q (1-(-q)^e)/(1+q), A_e = (1+q(-q)^e)/(1+q).

Both quotients belong to R for every integer e: their numerators are divisible
by 1+q as Laurent polynomials. This is not a parameter specialization or a
requirement to invert 1+q. The finite Laurent quotient then gives the
coefficients of s^k E. Division by the two unit extreme coefficients is a
simpler exact implementation than rational expressions.

## Why this is the current full output, not just an algebraic equality

Starting a rank-two singleton scan, every all-singleton terminal branch has
its active original bottom label2 always over. Its successive signs alternate
-1,+1,-1,+1,..., so its terminal exponent is0 or-1. The selected top port is
respectively2 or1. Both remaining singletons have empty component words;
the ordered top/bottom records and stable-partition routes are fixed by those
ports. In the actual serialization, identity is the fixed routed word(1,-1),
and the negative crossing is(-1). The identity routing is not casually dropped
when comparing full records.

A foe branch births X2. There are no outside strands, friends, layer-growth
or higher sliders in rank2. The X2 internal rule precedes the flype rule;
all remaining crossings commute internally. On extraction algorithm._extract_step
uses the sum of internal signed letters, emits s^k E with normal orientation,
twist_powers=[k], empty routes, full original portsets{1,2}, and no further
component. Thus k determines every field of its terminal record and word.

The singleton process terminates: a close step consumes a letter; a foe
identity branch consumes a letter; a foe reversed-crossing branch immediately
has that reversed crossing close and consumes it. Its X2 child has only finite
internal scans and extraction. Every finite input therefore has a finite tree.
Exact collection cannot change its evaluated Laurent element. The defining
single-foe identities and internal commutations are exact in R[s,s^-1], so
soundness gives the element s^e(w). All terminal terms lie in the basis above.
Uniqueness fixes every coefficient, and the fixed record construction fixes
all full output fields. This proves the claim independently of finite tests.

Source points: algorithm.py::_single_step, _extract_step (size1/size2), _lift;
upward.py::step/_step (i2 internal priority); extraction.py::split_layered and
_over_route. two_strand.py gives the same explicit complete schedule/family.

## Recursive reuse and scope

When extraction leaves exactly two strands, algorithm._extract_step initializes
an ordinary rank-two singleton scan with empty boundaries, world_phi=False,
fixed frozen contexts, and inherited top/bottom label lists. The local
calculation sees only relative positions; label lists are used to restore
original endpoints at extraction. For any two such remainders with equal
exponent sum and the SAME frozen contexts and label maps, the continuation
outputs agree in all fields. Fixed context append plus injective endpoint
renaming preserves the local output equality; deterministic _lift restores
the ambient word. This yields output-only memoization by exponent sum at those
recursive interfaces. One must retain the contexts/maps rather than equating
all ambient outputs sharing an exponent.

The support-block theorem likewise transports this identification to an
isolated two-strand support block and to its translated inclusions. It does
not establish inverse-cancellation invariance for a two-strand substring
interacting with other strands in a connected block. General-rank inverse
cancellation and the unique-index far-swap conjecture remain unproved here.

## Baseline count

At cutoffL10, NonPeelableWords(2,ell) totals2036 over ell<=10. Replace these
with exponent representatives1,...,10 and-2,...,-10:19 jobs. Exponent0 is
identity, exponent-1 is the already determined single-peel output. The current
safe residual baseline therefore decreases by2017 from42,454,830,949 to
42,454,828,932, plus one identity. This small change should not be presented
as a major standalone census saving. Reuse of rank-two remainders inside
larger computations may be more valuable; that runtime effect is unmeasured.

## Exact checks

All2,047 literal signed rank-two words of lengths0 through10 completed. Full
serialized terminal sums (including all routing words and component fields,
ignoring only term-list order) matched the pure-sign representative for their
exponent sum. Every certificate passed independent C++ replay;24 selected
certificates also passed cold Python replay. The2,047 certificates, output
hashes, step counts, verifierflags and source/build metadata are consolidated
in a source-bound research SQLite database; integrity_check passed. The report
records script hashes, verifier identities, and exact scope. These
checks support the source proof, and do not justify its general conclusion.
