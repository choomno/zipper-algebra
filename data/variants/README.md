# Exact historical matrix data

`square-zero-n3-regrep-qt.json` preserves a historical matrix export labelled
“left regular representation of the square-zero three-strand quotient over
Q(q,t).” Its status is `historical-matrix-export`. The recorded dimension 24
is the matrix size and the number of listed monomials. Verifying a basis or
representation claim requires the intended presentation and its mathematical
hypotheses separately.

## Format

The schema is `zipper-regular-representation-v1`:

- `coefficient_field` declares the rational-function field `Q(q,t)`, parameter
  order `q,t`, and characteristic zero.
- `standard_monomials` is the ordered list of 24 words. A word is a list of
  generator names; the empty list is the unit.
- `generators` fixes the order `x1,x2,y1,y2`.
- `rational_functions` is a shared pool of exact coefficients.
- `matrices[g].entries` is a list of triples `[row, column, function_index]`.
  Coordinates and pool indices start at zero. Unlisted entries are zero.
- `L_g[row,column]` is the coefficient of `standard_monomials[row]` in the
  product `g * standard_monomials[column]`. Thus multiplication by a generator
  acts on columns.

Each pooled rational function has `numerator` and `denominator` polynomials.
A polynomial is a list of `[[degree_q, degree_t], coefficient]` pairs. Degrees
are nonnegative integers, terms are ordered lexicographically by their degree
pairs, and coefficients are exact integer or reduced rational strings. An
empty polynomial list means zero. For example:

```json
{
  "numerator": [[[0, 0], "-1"], [[1, 0], "1"]],
  "denominator": [[[0, 1], "1"]]
}
```

This represents `(q-1)/t`. Numerator and denominator are coprime over Q[q,t];
the denominator is nonzero and has leading coefficient one in lexicographic
order with `q` before `t`. Specialization requires the denominator of each
used coefficient to remain nonzero. Pool indices have no mathematical meaning.

The file contains four 24-by-24 matrices, 560 stored nonzero entries and 104
pooled coefficients. It is plain JSON data, with no executable expressions or
pickle objects.

## Provenance and conversion

The source is `archive/variants-data/square-zero-n3-regrep-qt.py` at commit
`0db1dd9e2a33fb0e5d07e216c8a5354911446e97`, a 67,477-byte generated export
with SHA-256
`ef01808166a1ee5656f00e41e3665b440f067c37d68a3e0b764c1e080dc84787`.
The JSON retains those values and the original multiplication convention.

Conversion read only the source's expected symbolic constructors through a
strict AST interpreter, normalized the coefficients exactly over Q[q,t], and
wrote the sparse matrix records. A separately implemented fraction-field reader
checked all 24 basis words and all 2,304 matrix entries against the original
source. This establishes fidelity of the format conversion. The data remain a
historical export rather than a new representation or dimension proof.

The repository's MIT license applies to this data and its documentation.
