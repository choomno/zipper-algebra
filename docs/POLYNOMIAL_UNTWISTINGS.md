# Explicit polynomial untwistings

`Presentation.from_untwistings` accepts independent rational-coefficient
polynomials at selected strand levels. This extends the preset constructors with
an exact Python interface. For parameter-dependent coefficients over arbitrary
ordered rational-function fields, use the domain-bound `AlgebraPresentation`
interface in the [parameter-field guide](PARAMETER_FIELDS.md). The established
`Presentation` interface retains its symbolic `q,t` format.

```python
from fractions import Fraction
from zipper_algebra import Presentation, UntwistingPolynomial

presentation = Presentation.from_untwistings(3, {
    2: UntwistingPolynomial([Fraction(1, 2), -3, 1]),
    3: UntwistingPolynomial([2, 0, 1]),
})
presentation.emit_ms("input.ms", q=5, t=7, char=0)
```

Coefficients are listed in ascending order: `[a, b, 1]` means
`a + b*lambda + lambda**2`. Integers and `Fraction` values are exact; floats are
rejected. Trailing zeros are removed, but a nonmonic polynomial is not divided
by its leading coefficient. A zero polynomial is rejected; omit the level to
impose no relation there. Constant nonzero polynomials are permitted.

At level 2 the relation is `p_2(sigma_j) X_{2,j} = 0`. At higher levels it is
`p_i(Theta_{i-1,j}) X_{i,j} = 0`, where `Theta_{i-1,j}` is the shifted
lasso twist on the same `i` strands. All valid shifts are included. No missing higher polynomial
is inserted, and no additional relation is inferred by this constructor.

The definitions and finiteness hypotheses are those of Moos' *The Zipper
Algebra*, Definition 2.9 and Theorem 3.41; see the
[mathematical reference](mathematical-specification/mathematical-specification.pdf)
and [bibliography](../README.md#references).

## Raw presentations and finiteness hypotheses

The constructor accepts nonmonic polynomials and zero constant terms as raw
relations. Accepting such input makes no finiteness claim. The untwisting
theorem requires the appropriate relations at every level, monicity and unit
constant terms, together with an invertible `q` and the zipper relations.
Missing relations may instead be proved as consequences in a separate computation.

Numeric coefficient specialization follows the selected ring. Reducing rational
coefficients modulo a modulus requires every denominator to be invertible.
A leading or constant coefficient can become zero after specialization; the
polynomial's generic degree, monicity and constant term do not certify the
specialized theorem hypotheses. The numeric backend also permits residue rings;
a field-dependent algorithm must check its field hypotheses separately.

## Identity and verification

Explicit families use the direct polynomial relations above. Existing presets
retain their literal defining polynomials, relation order and presentation
hashes. The shortened higher relators of the presets are equivalent in the braid
quotient to full-twist relations, but are different free polynomials. An existing
certificate cannot be transferred between these presentations merely because
they define isomorphic algebras.

The existing cofactor APIs accept the resulting `Presentation` and bind
certificates to its ordered defining relations. The command-line presentation
flags still select presets; use the Python interface for an explicit family.

`UntwistingFamily` stores an immutable copy of a level-to-polynomial mapping.
Its `to_json()` and `from_json()` methods encode canonical decimal level names
and coefficient lists, using integers and rational strings such as `"1/2"`.
The family contains only polynomial data; retain the strand count, zipper-level
choices and additional relations when reconstructing the full presentation.

The tests include independent free-word expansions, shifted lasso words, and
companion-matrix examples for two strands. For monic `p` of degree `d`, with
`q` and `p(0)` units, the two-strand quotient is
`R[sigma, sigma^-1] / (p(sigma)*(sigma-1)*(sigma+q))`, free of rank `d+2`.
This supplies a known-answer check of the constructor without assuming that a
solver's basis is complete.
