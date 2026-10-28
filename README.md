# zipper-algebra

Research software for zipper algebras and their untwisting quotients: exact
noncommutative computations motivated by quantum topology and braid group
representations.

Bigelow introduced the zipper algebra as a quotient of the braid group algebra,
with the Hecke and BMW algebras as familiar quotients. He asked whether the
three-strand slider vanishes and which additional relations produce
finite-dimensional algebras ([Bigelow][1], §8, Questions 6–7). Moos's thesis develops a
family of such quotients, defined by polynomial untwisting relations on the
sliders. Its main finiteness theorem follows from the zipper algorithm and
provides the mathematical setting for this project (Moos, Theorem 3.41).

The software supports the study of that family: which defining relations are
independent, when sliders vanish, and how the algebra changes with the
untwisting polynomials. It constructs presentations, expands recorded proofs,
composes explicit two-sided ideal-membership certificates, reconstructs
candidate rational coefficients from exact samples, and verifies the resulting
identities.

Read the [compiled mathematical reference](docs/mathematical-specification/mathematical-specification.pdf)
([LaTeX source](docs/mathematical-specification/main.tex)) for the definitions,
finite-generation theorem, connection with BMW algebras, and precise software
contracts. It distinguishes the general polynomial family from the selected
presentations currently implemented.

## Mathematical setting

At each slider level one chooses a monic polynomial with invertible constant
term and imposes the corresponding untwisting relation. The theorem states
that the resulting algebra is finitely generated over the coefficient ring;
over a field, it is finite-dimensional. The layered-zipper spanning result and
reduction of twist exponents explain this conclusion (Moos, §3.5). It leaves
the dimensions, bases and structure of the individual quotients to further study.

The linear untwisting case leads back to the BMW algebra through an explicit
slider-vanishing identity under the stated parameter hypotheses (Moos,
Theorem 4.1 and Corollary 4.2). The package's
certificate tools make statements of this kind reproducible, with the
presentation and coefficient domain specified. General untwisting polynomials,
the geometric zipper algorithm, dimensions and trace calculations extend beyond
the maintained package's present scope.

## Install

Python 3.11 or newer is required. From a checkout:

```sh
python -m pip install .
zipper-algebra describe --n 4
```

The maintained algebra and certificate code uses Python's standard library.
Optional denominator-factor inference uses Sympy:

```sh
python -m pip install '.[inference]'
```

For development, use `python -m pip install -e .`. The installed package is named
`zipper_algebra`; `python -m zipper_algebra` provides the same commands as
`zipper-algebra`.

## Verify a certificate

From the repository root, check the retained historical free-algebra fixture:

```sh
zipper-algebra verify \
  logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt \
  --n 4 --target X_3
```

Add `--pit` for the second exact check using a degree-bounded product grid, or
`--json` for a machine-readable report. The explicit presentation options are
checked against the certificate's presentation hash. Exit codes distinguish a
valid identity, a false identity, and invalid input; see [WORKFLOW](WORKFLOW.md).
This retained large certificate is a regression fixture, rather than a claim
that it is the latest mathematical proof.

The same operation is available in Python:

```python
from zipper_algebra import Cofactor, Presentation
from zipper_algebra.cofactor.verify import check

presentation = Presentation(n=4)
certificate = Cofactor.load(
    "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt",
    pres=presentation,
)
report = check(certificate, presentation, target="X_3", verbose=False)
assert report["pass"]
```

## Mathematical interpretation

An accepted certificate gives an explicit identity over the specified rational
function field. Its presentation and denominator conditions are part of the
statement. Verification compares polynomials in the free associative algebra;
words representing the same braid are not automatically identified. A
certificate simplified in a braid group algebra therefore needs the appropriate
checker or explicit braid and inverse relation terms.

Reconstruction supplies candidates. A fit at sampled parameter values becomes
a membership result only after exact identity verification. A nonzero remainder
from an uncertified basis establishes no nonmembership conclusion, even at the
sampled parameter value.
See the specification for these distinctions and for the relation formulas.

## Tests and development

After an editable install:

```sh
python -m unittest discover -s tests -v
```

[TESTING](docs/TESTING.md) describes the CI suite, historical compatibility
checks, package checks and optional integrations.
[CONTRIBUTING](CONTRIBUTING.md) records review and commit practice, including
independent adversarial review for complex changes.

## Repository map

| Location | Purpose |
| --- | --- |
| `src/zipper_algebra/` | Maintained Python implementation and command line interface. |
| `tests/` | Regression tests, installation checks and immutable fixtures. |
| `docs/` | Mathematical specification, public scope and validation records. |
| `tools/` | Compatibility entry points and experimental solver/research scripts. |
| `logs_and_traces/`, `ms_inputs/` | Historical scientific artifacts and inputs, retained with their provenance. |
| `gap_scripts/`, `singular_run/`, `archive/` | Independent or older computational paths and saved experiments. |

The installed wheel contains the maintained package, not the private recovery
workspace or the research archive. Solver execution, variant surveys, dimension
and trace calculations remain experimental while their completion and exactness
requirements are audited. [INDEX](INDEX.md) maps those sources;
[RESEARCH_AGENDA](RESEARCH_AGENDA.md) records the broader publication programme.

The current release work and remaining decisions are tracked in
[PUBLIC_SCOPE](docs/PUBLIC_SCOPE.md). History reconstruction and attribution
history cleanup are recorded separately from the software changes.

## References

[1]: https://web.math.ucsb.edu/~bigelow/publications/10.pdf
[2]: docs/mathematical-specification/mathematical-specification.pdf

1. Stephen Bigelow, *Braid groups and Iwahori–Hecke algebras*, in *Problems on
   Mapping Class Groups and Related Topics*, Proceedings of Symposia in Pure
   Mathematics **74** (2006), 285–299.
   [Author's copy][1]; [arXiv:math/0505064](https://arxiv.org/abs/math/0505064).
2. Choomno Moos, *The Zipper Algebra*, M.A. thesis, University of California,
   Santa Barbara, September 2026. Cited by section and theorem numbers in the
   recovered manuscript. The [mathematical reference][2] supplies the relevant
   definitions, theorem statement and version information.
3. Joan S. Birman and Hans Wenzl, *Braids, link polynomials and a new algebra*,
   Transactions of the American Mathematical Society **313** (1989), 249–273.
   [DOI](https://doi.org/10.1090/S0002-9947-1989-0992598-X).
4. Jun Murakami, *The Kauffman polynomial of links and representation theory*,
   Osaka Journal of Mathematics **24** (1987), 745–758.

## License

The software is available under the [MIT license](LICENSE).
