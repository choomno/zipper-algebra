# Complete Gröbner bases and certified normal forms

An exact cofactor identity proves ideal membership. A complete Gröbner basis
additionally decides membership by reduction and identifies a vector-space
basis of the quotient. The explicit verifier in
[`verification/groebner.py`](../src/zipper_algebra/verification/groebner.py)
checks these additional assertions over the declared coefficient field.
It accepts the exact fields described in the
[parameter-field guide](PARAMETER_FIELDS.md): the rationals, checked prime
fields, and rational-function fields over either. Generic parameter fields are
the principal research setting; specialization checks concern the specified
fibres and their denominator conditions.

## Data and mathematical checks

Let `presentation` describe the free algebra `F` and the ordered input relations
`r_i`, and let `I` be their two-sided ideal. A `GroebnerCertificate` contains the
full `presentation.content_hash()`, an ordered tuple of nonzero candidate basis
polynomials, and one tuple of `ProofTerm` values for each polynomial. A proof
term is a coefficient times a left word, a referenced polynomial and a right
word. Its `kind` is `"input"` or `"basis"`; its `index` is zero-based. Context
words use one-based generator indices and `()` for the empty word.

Every local derivation must be a literal polynomial equality in `F`. References
to other basis polynomials can point forward, but their directed dependency
graph must be acyclic. References are validated even when their coefficients
are zero. Induction along the graph then proves that every candidate belongs
to `I`. Checking these local equalities avoids expanding shared derivations
into a potentially enormous list of cofactors of the original relations.

The verifier uses degree-lexicographic order: compare word lengths first, then
the tuples of indices in the declared alphabet. It computes leading words and
normalizes leading coefficients over the declared field. It then checks:

1. Every original input relation reduces to zero. Together with the derivations,
   this proves equality of the two-sided input and candidate ideals.
2. Every inclusion ambiguity reduces to zero, including distinct rules with
   equal leading words.
3. Every proper nonempty suffix-prefix overlap reduces to zero, for ordered
   pairs of rules including self-overlaps.

The degree-lex order decreases at each rewrite. Resolution of these ambiguities
therefore gives unique normal forms by Bergman’s Diamond Lemma [1]. Words
containing no leading word as a subword form a basis of `F/I`. A candidate
polynomial proved to be a nonzero constant gives the unit ideal and dimension
zero.

`prune=True` selects reduction rules by leading-word divisibility. It still
checks the complete derivation graph, all original input relations and all
critical ambiguities for the selected rules. Leading-word divisibility alone
does not establish that a discarded rule is redundant.

## Python interface

```python
from zipper_algebra.verification.groebner import (
    GroebnerCertificate, ProofTerm, VerificationLimit, verify_groebner,
)

# presentation and certificate must describe the same field, alphabet and input.
verified = verify_groebner(
    presentation, certificate,
    max_reductions=2_000_000, max_words=100_000, prune=False,
)
remainder = verified.normal_form(target)
dimension = verified.dimension
standard_words = verified.standard_words
report = verified.report
```

This is an API sketch: the caller supplies the presentation, certificate and
target. The verifier constructs no candidate basis. `ProofTerm` has positional
fields `(coefficient, left, kind, index, right)`, with an exact `FieldElement`
coefficient. `GroebnerCertificate` has positional fields
`(presentation_hash, basis, derivations)`. Its `to_json` and `from_json` methods
use schema `zipper-groebner-certificate-v1`, recording the complete domain-bound
polynomials and proof terms. Arithmetic equality of fractions does not erase
the literal denominator data used in the presentation hash.

Successful verification returns a `VerifiedReduction`. A zero normal form
proves membership in the original defining ideal, and a nonzero normal form
proves nonmembership over the declared field. `normal_form` requires exactly
the verified free-algebra parent, including for a zero target. The defining data
are immutable. `report` returns a copy of the verification counts; later normal
forms share the memoized reductions and remaining work budget, but do not
change those recorded counts.

The two limits have different meanings:

| Outcome | Meaning |
| --- | --- |
| `VerificationLimit` | The rewrite/ambiguity work budget ended; the requested operation did not finish. |
| `dimension` is an integer | Standard-word enumeration finished after completion was proved; its length is the exact dimension. |
| `dimension is None`, `standard_words is None` | Completion was proved, but enumeration exceeded `max_words`; dimension is undetermined. |
| `TypeError` or `ValueError` | Data are malformed, incompatible, or fail a required identity or ambiguity check. |

The work budget counts uncached word rewrites and tested ambiguities. Exhaustion
during a later normal-form call does not revoke the previously established
completion. An enumeration limit gives no assertion of infinite dimension.

## Solver text and other verification routes

[`certificate_from_solver`](../src/zipper_algebra/verification/f4ncgb.py) accepts
`(presentation, basis_text, proof_text)` and decodes numeric f4ncgb output into
a candidate `GroebnerCertificate`. Rational numerical coefficients are embedded
exactly in the declared field. Parsing supplies neither a membership proof nor
a completion verdict; call `verify_groebner` to establish those conclusions.
This decoder is not a producer for bases over arbitrary parameter fields. The
[four-strand square-zero example](FOUR_STRAND_SQUARE_ZERO.md) supplies a
separate local Singular producer over `Q(q,t)` with independent GAP checks.

The established discovery pipeline continues to record `basis_complete=False`.
Its proof checks establish membership of returned polynomials by local acyclic
derivations; a successful solver exit is not a completion test. The explicit
completed-basis interface described here is a separate operation. The existing
[X₃ reproduction route](X3_REPRODUCTION.md) produces an exact generic cofactor
identity and does not require a complete generic Gröbner basis.

GAP’s GBNP package provides an independent implementation: its
[`IsGrobnerBasis`](https://gap-packages.github.io/gbnp/doc/chap3.html) tests a
supplied set over its supported fields. Its quotient-dimension and basis
operations require an actual Gröbner basis [2]. Such a check establishes
completion for the supplied ideal; comparison with the original presentation
also requires evidence for both ideal inclusions. This is distinct from the
existing GAP exporter, which checks a generic cofactor identity by expansion.

## References

1. George M. Bergman, *The diamond lemma for ring theory*, Advances in
   Mathematics **29** (1978), 178–218, Theorem 1.2.
   [DOI](https://doi.org/10.1016/0001-8708(78)90010-5);
   [author’s corrections and updates](https://math.berkeley.edu/~gbergman/papers/updates/diamond.html).
2. A. M. Cohen and J. W. Knopper, *GBNP*, version 1.1.0, GAP package.
   [Official function manual](https://gap-packages.github.io/gbnp/doc/chap3.html).
