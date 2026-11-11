# Exact parameter fields and worked quotients

The Python API supports rational-function fields with an arbitrary ordered list
of independent parameters, over the rationals or a checked prime field. It also
supports free algebras, literal presentations and exact cofactor identities over
these fields. No external algebra system is required for this interface.

```python
from zipper_algebra import AlgebraPresentation, RationalFunctionField

field = RationalFunctionField(("q", "a", "b"))
q, a, b = (field.parameter(name) for name in field.parameters)
presentation = AlgebraPresentation.zipper(
    2, field=field, q=q, untwistings={2: [b, a, 1]},
)
```

The coefficient list specifies the quadratic `b + a*lambda + lambda**2`.
The constructor includes every valid shift, and omitted levels remain omitted.
It accepts raw nonmonic polynomials and zero constant terms. The
[compiled reference](mathematical-specification/mathematical-specification.pdf)
states the finite-generation hypotheses and proves the two-strand calculations
used below. The family and theorem are from Moos' thesis, Definition 2.9 and
Theorem 3.41; see the [references](../README.md#references).

Trailing zero coefficients are removed from the input polynomial. A zero
polynomial is rejected; omit the level to impose no untwisting there.

## Small computations

After installing the package, run from the checkout:

```sh
python examples/polynomial_quotients.py
```

These computations use exact cofactor identities and companion representations;
they do not run a Gröbner-basis search.

| Computation | Coefficient field | Verified dimension |
| --- | --- | --- |
| General quadratic with coefficients `b,a,1` | `Q(q,a,b)` | 4 |
| Its specialization `q=2, a=1, b=3` | `GF(5)` | 4 |
| Repeated-root polynomial `(lambda-u)**2`, with `q=2` | `Q(u)` | 4 |
| Polynomial `s*lambda**2 + lambda + 1`, with `q=2` | `Q(s)` | 4 |
| The same polynomial at `s=0` | `Q` | 3 |
| Polynomial `lambda**2-s`, with `q=2` | `Q(s)` | 4 |
| The same polynomial at `s=0` | `Q` | 2 |

Every computation checks the literal defining relations in a companion
representation, checks that powers of the generator applied to a cyclic vector
give a coordinate basis, and verifies a cofactor identity giving the matching
upper bound. The repeated-root example also checks a nonzero square-zero
element, including an exact cofactor for its square.

The last example uses the invertibility of the braid generator to remove a
square factor. Its dimension change has a different cause from the preceding
degree-drop example. Neither generic dimension nor a generic basis is assumed
to survive specialization.

To save the defining data and witnesses in one JSON file, with a separate
source/output provenance record:

```sh
python examples/polynomial_quotients.py --output runs/polynomial-quotients/result.json
```

The `runs/` directory is ignored. The saved cofactor identities can be loaded
with `CofactorIdentity.from_json` and checked against the saved
`AlgebraPresentation.from_json` data. These examples are small enough to run in
the ordinary test suite.

## Domains and specialization

`RationalFunctionField(())` is the rational field; setting `characteristic=5`
gives the prime field. Positive characteristics are verified by exact trial
division. Parameter order belongs to the domain: `Q(a,b)` and the separately
declared `Q(b,a)` require an explicit change of coordinates. Arithmetic rejects
mixed domains, including operations involving zero. Floats and booleans are
rejected as coefficients.

`field.polynomial` accepts a sparse mapping from exponent tuples to exact
integers or `Fraction` values. Ordinary arithmetic then constructs rational
functions. Equality uses cross multiplication. Fractions retain their
polynomial denominators without a polynomial gcd calculation; they are
unhashable because equal values can have different stored representatives.

`element.specialize(assignments, characteristic=...)` evaluates the stored
expression at a complete exact parameter point. It rejects undefined rational
constants and denominator zeros. For example, `(s*s-1)/(s-1)` equals `s+1` in
the function field, but its stored denominator excludes `s=1`. Specialization
applies to these expressions on their common domain of definition; it is not
a homomorphism from the entire rational-function field to the specialized field.
Arithmetic may change stored representatives, including cancellation to zero.

The interface implements fields. It does not certify ideal membership over a
polynomial coefficient ring, a smaller localization, or an algebraic extension.
Such conclusions require their own declared domains and unit checks.

## Presentations and certificates

`FreeAlgebra` records the coefficient field and ordered generators. Its
immutable elements retain that parent even when zero. Words use one-based
generator indices; multiplication concatenates words without braid rewriting.

`AlgebraPresentation` records the ordered literal relations and distinct labels.
Its full SHA-256 includes the versioned domain, alphabet and coefficient data.
Equivalent rational-function expressions can conservatively produce different
literal presentation hashes. No old certificate hash is reused.
The hash uses the validated, normalized JSON representation rather than the
bytes of a particular JSON file; spacing and equivalent scalar spellings do
not change it.

`CofactorIdentity` stores a target and triples `(left, relation_index, right)`,
where both factors are free polynomials and relation indices start at zero.
`verify(presentation)` checks an exact identity over the declared field. Invalid
references or mismatched presentation identities raise errors; a well-formed
false identity returns `False`. `save(path)` writes JSON and automatic provenance.
JSON readers validate exact scalars, domains, word indices and references.

The established `Presentation`, `UntwistingPolynomial` and text `Cofactor`
interfaces retain their existing rational-input and symbolic `q,t` behavior.
Their solver, reconstruction and independent GAP/Singular routes continue to
use that format. The new domain-bound interface currently provides exact
construction and verification through Python; solver-backed production and
generalized reconstruction have further integration work ahead.
