# zipper-algebra

Research software for zipper algebras and their untwisted quotients. This software conducts exact
noncommutative computations motivated by quantum topology and braid group
representations.

Bigelow introduced the zipper algebra as a quotient of the braid group algebra,
with the Hecke and BMW algebras as familiar quotients. He asked whether the
three-strand slider vanishes and which additional relations produce
finite-dimensional algebras ([Bigelow][1], §8, Questions 6–7). Moos' thesis develops a
family of such quotients, defined by polynomial untwisting relations on the
sliders. Its main finiteness theorem follows from the zipper algorithm and
provides the mathematical setting for this project (Moos, Theorem 3.41).

This repository supports the study of that family. It's components compute: which defining relations are
independent, when sliders vanish, and how the algebra changes with a choice of
untwisting polynomials. It constructs presentations, expands recorded proofs,
composes explicit two-sided ideal-membership certificates via cofactor representations, reconstructs
candidate rational coefficients from exact samples, and verifies the resulting
identities.

Read the [mathematical reference](docs/mathematical-specification/mathematical-specification.pdf) for the definitions,
finite-generation theorem, connection with BMW algebras, and precise input requirements and verification guarantees. It distinguishes the general polynomial family from the selected presentations currently implemented.

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

## Algebra solver and citation

The original Gröbner-basis and cofactor-discovery computations used
[f4ncgb](https://gitlab.sai.jku.at/f4ncgb/f4ncgb), developed by Maximilian
Heisinger and Clemens Hofstadler. The original project is the authoritative source
for the solver and its documentation.

The certificate audit below runs from saved data with the Python package.
For fresh solver computations, the [setup guide](patches/f4ncgb/README.md)
records the upstream revision, parser fix, dependencies and maintained GitHub copy.

## X_3 cofactor computation

Instructions to reproduce the $X_3$ certificate claimed in Moos' thesis is
provided. The full process is as follows: input the defining relations,
compute a basis with nonexpanded proofs using f4ncgb, reduce, repeat at
parameter samples, reconstruct the relevant proof coefficients, compose the
cofactor expression, and verify the resulting identity exactly. Note that 
cofactor representations are not unique in their term count or coefficients.

The [reproduction guide](docs/X3_REPRODUCTION.md) identifies the
input, program and output at every stage, including the original presentation
input and the files needed for reconstruction. It distinguishes those stages
from replaying a saved symbolic lift and states the current driver's workspace
dependencies and missing historical proof logs.

### Audit the result

Run from the repository root after installation. This checks the saved 695-term
cofactor identity exactly, using both symbolic accumulation and a degree-bounded
product grid:

```sh
python -m zipper_algebra describe --n 4 --json
python -m zipper_algebra verify \
  logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt \
  --n 4 --target X_3 --pit --json > x3-verification.json
```

Success means exit code `0` and `"pass": true` in the overall, `exact`, and `pit`
reports in `x3-verification.json`; the exact report records 695 terms. The
presentation options select the four-strand canonical untwisting presentation.
The command checks a generic rational-function identity, with specialization
conditions explained in the mathematical reference. It needs no algebra solver.

To reproduce the certificate body from the saved symbolic coefficients and
reduction skeleton, then replay the retained numeric proof expansion:

```sh
python -m unittest \
  tools.cofactor.tests_compose_io_verify.TestCompose.test_composed_body_is_byte_identical_to_the_shipped_certificate \
  -v
python -m unittest tools.cofactor.tests_proof_skeleton -v
```

The first command must report one passing test and byte-for-byte agreement;
the second must report 16 passing tests. These replay saved proof data; they do
not rerun the original sampling and interpolation.

The [reproduction guide](docs/X3_REPRODUCTION.md) maps every substantial source,
input and output to its role, including the later compact certificate and the
independent specialized calculations. Its [file inventory](docs/x3-reproduction-files.json)
records exact paths, sizes and SHA-256 hashes. Check those against the checkout:

```sh
python ci/check_x3_manifest.py
```

Fresh discovery uses the experimental solver driver and needs additional
configuration. The [f4ncgb parser patch and build instructions](patches/f4ncgb/README.md)
are included for that route. The guide states the remaining portability and
archival gaps. The retained large certificate is auditable computational
evidence, rather than the latest mathematical proof.

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

## Tests and development

After an editable install:

```sh
python -m unittest discover -s tests -v
```

[TESTING](docs/TESTING.md) describes the CI suite, historical compatibility
checks, package checks and optional integrations.

## Repository map

| Location | Purpose |
| --- | --- |
| `src/zipper_algebra/` | Maintained Python implementation and command line interface. |
| `tests/` | Regression tests, installation checks and immutable fixtures. |
| `docs/` | Mathematical reference, reproduction guide and testing instructions. |
| `tools/` | Compatibility entry points and experimental solver/research scripts. |
| `logs_and_traces/`, `ms_inputs/` | Historical scientific artifacts and inputs, retained with their provenance. |
| `gap_scripts/`, `singular_run/`, `archive/` | Independent or older computational paths and saved experiments. |

The installed wheel contains the maintained package, not the private recovery
workspace or the research archive. Solver execution, variant surveys, dimension
and trace calculations remain experimental while their completion and exactness
requirements are audited. [INDEX](INDEX.md) maps those sources.

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

5. Maximilian Heisinger and Clemens Hofstadler, *f4ncgb: High Performance
   Gröbner Basis Computations in Free Algebras*, in *Computer Algebra in
   Scientific Computing (CASC 2025)*, Lecture Notes in Computer Science
   **16235**, Springer Nature Switzerland (2025), 79–97.
   [DOI](https://doi.org/10.1007/978-3-032-09645-6_5).

## License

The software is available under the [MIT license](LICENSE).
