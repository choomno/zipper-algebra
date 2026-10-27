# zipper-algebra

Exact algebraic computations for quantum topology research on zipper relations,
sliders, and untwisting quotients of braid group algebras. The questions include
which slider relations follow from a chosen presentation and how changes to the
untwisting relations affect quotients related to BMW algebras.

The Python package constructs presentations, expands recorded proofs, composes
explicit two-sided ideal-membership certificates, reconstructs candidate
coefficients from exact samples, and verifies the resulting identities. The
[mathematical specification](docs/mathematical-specification/main.tex) defines
the coefficient domains, conventions, hypotheses and exact conclusions.

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

## License

The software is available under the [MIT license](LICENSE).
