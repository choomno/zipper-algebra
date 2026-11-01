# Certified element and orbit minimal polynomials

Given an element `a` and a seed `v` in a verified finite algebra, the operation
computes the monic polynomial of least degree satisfying $m(a)v=0$. The default
seed is the unit, giving the minimal polynomial of `a` itself. An explicit seed
gives the minimal polynomial on its cyclic orbit. These can have different
degrees: a relation on one slider need not annihilate the full algebra.

The input must be a `VerifiedFiniteAlgebra` returned by `verify_word_basis`.
An arbitrary representation can be a proper quotient and does not establish
annihilation in the original presented algebra. See the
[finite-algebra guide](FINITE_ALGEBRAS.md) for the required spanning and
independence proofs.

## Production and independent checking

```python
from zipper_algebra.invariants import minimal_polynomial
from zipper_algebra.verification.krylov import (
    MinimalPolynomialCertificate, verify_minimal_polynomial,
)

# A is a verified finite algebra; a and v belong to A.
full = minimal_polynomial(A, a, threads=2)
relative = minimal_polynomial(A, a, seed=v, threads=2)

# Coefficients are exact field elements, in ascending power order.
assert len(relative.coefficients) == relative.degree + 1
reloaded = MinimalPolynomialCertificate.from_json(relative.to_json())
verify_minimal_polynomial(A, reloaded)
```

The producer uses [FiniteFlow](RECONSTRUCTION.md) over the rationals and rational
function fields. It forms the Krylov columns $v,av,a^2v,...$ on demand, grows
prefixes geometrically, and reconstructs exact kernels with the existing
certified matrix interface. Once a prefix is dependent, its exact cyclic rank
identifies the first dependency degree. A shorter prefix is reconstructed when
needed. No full regular multiplication matrix or cubic multiplication table
is formed.

For degree `d`, the independent verifier recomputes the first `d+1` columns from
the certified algebra. It checks the complete kernel identity and both rank
bounds, requires rank `d` and a one-dimensional kernel, and checks the monic
relation against its kernel vector. The last column is a combination of the
preceding columns, so rank `d` proves those preceding columns independent.
This proves minimality as well as annihilation. The
[compiled specification](mathematical-specification/mathematical-specification.pdf)
gives the argument precisely.

The checker does not import FiniteFlow. It supports the fields and modular
minor witnesses of the independent kernel verifier, including supplied
positive-characteristic certificates. Nonzero-orbit production currently
requires characteristic zero. A zero seed, including the zero algebra, has
minimal polynomial `1` and uses a direct certificate without FiniteFlow.
Parameter-field verification requires the symbolic extra.

## Bounds and data

`max_degree` bounds the orbit search and is capped by the verified dimension.
The default is that dimension; an exhausted smaller degree budget raises
`VerificationLimit`. It returns no minimal polynomial from an unfinished
search. `max_operations` applies to each multiplication; `max_products` applies
to each independent kernel check. Verification recomputes the required orbit
with fresh budgets. The native prime, reconstruction-degree, thread, starting
prime and shift options have the meanings in the reconstruction guide.
These settings do not impose process-level wall-time or memory limits.

`MinimalPolynomialCertificate` retains the full presentation hash, ordered word
basis, element and seed coordinates, monic coefficients, and kernel witness.
Its strict `zipper-minimal-polynomial-certificate-v1` JSON encoding stores the
element and seed sparsely. `to_json` and `from_json` do not write files. Saved
workflow outputs should also retain [source and input provenance](PROVENANCE.md).

## A generic four-strand computation

Use the generic 186-dimensional square-zero algebra from the
[four-strand finite-basis reproduction](FINITE_ALGEBRAS.md#reproducing-the-four-strand-finite-algebra).
After producing its presentation and word-basis certificate:

```python
import json
from pathlib import Path
from zipper_algebra.algebra import AlgebraPresentation
from zipper_algebra.verification.finite import WordBasisCertificate, verify_word_basis
from zipper_algebra.invariants import minimal_polynomial

run = Path("runs/four-strand-finite-basis")
p = AlgebraPresentation.from_json(json.loads((run / "presentation.json").read_text()))
c = WordBasisCertificate.from_json(
    json.loads((run / "word-basis-certificate.json").read_text()))
A = verify_word_basis(p, c)
q, t = (A.field.parameter(name) for name in ("q", "t"))
x, y = (p.algebra.generator(name) for name in ("x1", "y1"))
X2 = q*y + (1-q) - x
a = A.element(A.coordinates(x))
v = A.element(A.coordinates(X2))
w = A.element(A.coordinates((x-t)*X2))
full = minimal_polynomial(A, a)
relative = minimal_polynomial(A, a, seed=v)
nilpotent = minimal_polynomial(A, w)
```

A local run computed and independently checked, over $K=\mathbb{Q}(q,t)$:

| Calculation | Certified minimal polynomial |
| --- | --- |
| `x1` in the full algebra | $(z-t)^2(z-1)(z+q)$ |
| Left multiplication by `x1` on the orbit of `X2` | $(z-t)^2$ |
| The element `(x1-t)*X2` | $z^2$ |

The degree-two polynomial in the last row also proves that this element is
nonzero. These results were computed from the verified actions; the displayed
polynomials were not supplied as reconstruction inputs. All three serialized
certificates were reloaded and checked independently. This local example is
outside CI; bounded noncommutative controls exercise the routine tests.
