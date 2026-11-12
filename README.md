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

This repository supports the study of that family. Its components investigate
dependencies among defining relations, slider vanishing, and how the algebra changes with a choice of
untwisting polynomials. It constructs presentations, expands recorded proofs,
composes explicit two-sided ideal-membership certificates via cofactor representations, reconstructs
candidate rational coefficients with FiniteFlow, and verifies the resulting
identities.

Read the [mathematical reference](docs/mathematical-specification/mathematical-specification.pdf) for the definitions,
finite-generation theorem, connection with BMW algebras, and precise input requirements and verification guarantees. It distinguishes the general polynomial family from the selected presentations currently implemented.

The [parameter-field interface](docs/PARAMETER_FIELDS.md) supports exact
free-algebra presentations and cofactor identities over arbitrary rational-function
fields. The [completed-basis verifier](docs/GROEBNER_VERIFICATION.md) checks
derivations and all critical ambiguities before providing certified normal forms
and, when standard-word enumeration finishes, the exact quotient dimension.
These interfaces run with the Python package alone; producing a candidate basis
is a separate computation. [Finite word-basis certificates](docs/FINITE_ALGEBRAS.md)
also give certified coordinates and products from sparse generator actions and
explicit spanning and independence proofs. Parameter-field verification on this
route uses the optional symbolic arithmetic dependency.

The [FiniteFlow reconstruction interface](docs/RECONSTRUCTION.md) constructs
right kernels over rational-function fields and independently verifies both
rank bounds and the complete generic kernel identity. Its native dependency
and the maintained fork's reliability controls have separate installation
instructions; checking supplied parameter-field kernel certificates uses the symbolic
extra and does not load FiniteFlow.

The [four-strand square-zero computation](docs/FOUR_STRAND_SQUARE_ZERO.md)
works over $\mathbb{Q}(q,t)$: the quotient has dimension 186 and contains a
nonzero square-zero two-sided ideal. Its reproduction checks connect the
original presentation to the completed basis and verify the ideal calculation.
The [finite word-basis route](docs/FINITE_ALGEBRAS.md#reproducing-the-four-strand-finite-algebra)
additionally verifies its 186-word basis and sparse regular actions with the
maintained Python API, providing exact coordinates and products.

## Install

Use Linux with Python 3.11 or newer. From a checkout:

```sh
python -m pip install .
zipper-algebra describe --n 4
```

Presentation construction and the existing cofactor and Gröbner checks use
Python's standard library. Finite word-basis verification over parameter fields
uses SymPy's exact fraction fields:

```sh
python -m pip install '.[symbolic]'
```

Optional denominator-factor inference also uses SymPy, through
`python -m pip install '.[inference]'`.

For development, install the tools under [Tests and development](#tests-and-development).
The installed package is named `zipper_algebra`; `python -m zipper_algebra`
provides the same commands as `zipper-algebra`.

For **fresh Gröbner-basis and cofactor discovery**, also build the external
[f4ncgb solver](https://gitlab.sai.jku.at/f4ncgb/f4ncgb), by Maximilian Heisinger
and Clemens Hofstadler. Follow the [solver installation instructions](patches/f4ncgb/README.md)
for the pinned upstream revision, parser patch, C++ dependencies and build commands.
Python package installation does not install this executable. Sample data are
stored in SQLite. The default solver interface uses temporary files; an
[optional SQLite-enabled solver build](patches/f4ncgb/README.md#optional-direct-sqlite-io)
accepts `--solver-io sqlite` for direct database input and output. The discovery
command takes the solver, output directory, cache and limits explicitly;
the [reproduction guide](docs/X3_REPRODUCTION.md) explains its stages and outputs.
Fresh cofactor discovery additionally requires the symbolic extra and the
[maintained FiniteFlow build](docs/RECONSTRUCTION.md#installation). It obtains
word contexts from f4ncgb proofs and reconstructs their coefficients with FiniteFlow.

Presentation construction, checking saved certificates, recomposition from saved
proof data and the maintained tests run without the solver. Optional factor
inference needs the Sympy extra above.

## Algebra solver and citation

The original Gröbner-basis and cofactor-discovery computations used
[f4ncgb](https://gitlab.sai.jku.at/f4ncgb/f4ncgb), developed by Maximilian
Heisinger and Clemens Hofstadler. Solver-backed discovery uses a parser-patched
build; the patch corrects line endings and identifier validation. Their project
is the authoritative source for the solver and its documentation; their
[CASC paper](https://doi.org/10.1007/978-3-032-09645-6_5) is listed below.

The certificate audit below runs from saved data with the Python package.
For fresh solver computations, the [setup guide](patches/f4ncgb/README.md)
records the upstream revision, parser fix, dependencies and maintained GitHub copy.

## X_3 cofactor computation

This repository documents the production and verification of a computational
cofactor certificate for the $X_3$ identity discussed in Moos' thesis.
The full process is as follows: input the defining relations,
compute a basis with nonexpanded proofs using f4ncgb at one anchor, reduce the
target, expand the proof graph into word contexts, form the generic linear
system, reconstruct its kernel with FiniteFlow, and verify the resulting
cofactor identity exactly. Note that
cofactor representations are not unique in their term count or coefficients.

The [reproduction guide](docs/X3_REPRODUCTION.md) identifies the
input, program and output at every stage, including the original presentation
input and the files needed for reconstruction. Fresh production generates its own
anchor proof and word contexts; checking the curated certificate is a separate route.

### Produce a new certificate

Install the symbolic extra, build f4ncgb, and install the
[FiniteFlow dependency](docs/RECONSTRUCTION.md#installation) using the instructions above.
Then run from any working directory:

```sh
zipper-algebra discover --n 4 --target X_3 \
  --solver /absolute/path/to/f4ncgb \
  --solver-revision 789504d7c7da6af6946faecb73421b037766a5b8 \
  --output runs/x3 --cache runs/solver-cache \
  --anchor 41,43 --max-primes 5 --reconstruction-degree 24 --threads 1 --json
```

The command emits the defining presentation and validates every numeric solver
proof. It expands the relevant proof graph into contextual occurrences of the
original relations, then forms a matrix over $\mathbb{Q}(q,t)$. FiniteFlow
reconstructs a complete kernel of that matrix. Independent exact checks prove
the kernel identity and rank, then verify the full free-algebra cofactor identity.
Saved cofactors, lifts and denominator guesses are not inputs.

The result is recorded in `runs/x3/result.json`. Success has `"status": "verified"`;
the report names an immutable certificate bundle and its hashes. Exit code `3`
means discovery was inconclusive, for example because contextual support was
insufficient or a computation limit was exhausted. The full four-strand
search can be expensive; a small three-strand example and the exact source map
are in the guide.

A successful run stores its anchor proof in one database and retains a compact
bundle of support, matrix, kernel and certificate records. The guide explains the
[stored products and sample export](docs/X3_REPRODUCTION.md#read-the-result).

Generated computation files also identify their source commit, local changes
and source-file hashes. See [output provenance](docs/PROVENANCE.md) for bundle
records, standalone sidecars and installed-package build origins.

### Audit the result

Run from the repository root after installation. The curated example contains a
129-term cofactor identity. Check it exactly using symbolic accumulation and a
degree-bounded product grid:

```sh
python -m zipper_algebra describe --n 4 --json
python -m zipper_algebra verify \
  examples/x3/certificate.txt \
  --n 4 --target X_3 --pit --json > x3-verification.json
```

Success means exit code `0` and `"pass": true` in the overall, `exact`, and `pit`
reports in `x3-verification.json`; the exact report records 129 terms. The
presentation options select the four-strand canonical untwisting presentation.
The command checks a generic rational-function identity, with specialization
conditions explained in the mathematical reference. It needs no algebra solver.

For an independent generic identity check in GAP, export a standalone script
with `python -m zipper_algebra.verification.gap`; the
[reproduction guide](docs/X3_REPRODUCTION.md#independent-algebra-system-checks)
gives the complete commands and expected result.

The [reproduction guide](docs/X3_REPRODUCTION.md) maps every substantial source,
input and output to its role. Its [file inventory](docs/x3-reproduction-files.json)
records exact paths, sizes and SHA-256 hashes. Check those against the checkout:

```sh
python ci/check_x3_manifest.py
```

The [f4ncgb parser patch and build instructions](patches/f4ncgb/README.md)
provide the external Gröbner dependency for fresh discovery. FiniteFlow supplies
the scalar reconstruction stage. The curated 129-term certificate was produced
by the earlier sampled fitting route. Fresh production can yield a different
valid cofactor representation; see the reproduction guide for the current run.
Full four-strand computations remain outside CI.

The same operation is available in Python:

```python
from zipper_algebra import Cofactor, Presentation
from zipper_algebra.cofactor.verify import check

presentation = Presentation(n=4)
certificate = Cofactor.load(
    "examples/x3/certificate.txt",
    pres=presentation,
)
report = check(certificate, presentation, target="X_3", verbose=False)
assert report["pass"]
```

## Explicit polynomial untwistings

Use `Presentation.from_untwistings` to supply independent exact polynomials at
selected strand levels. The [polynomial-input guide](docs/POLYNOMIAL_UNTWISTINGS.md)
explains rational coefficients, omitted relations, specialization and verification.
The existing preset presentations and certificate route remain available.

## Tests and development

See the [pull request guidelines](https://github.com/choomno/zipper-algebra/blob/main/docs/PULL_REQUESTS.md)
for branch targets, review scope, and how to report evidence.

Install the development tools and optional inference dependency, then run the
same static checks and test suite used in CI:

```sh
python -m pip install -e '.[dev,inference,symbolic]'
python -m ruff check src tests experiments examples ci setup.py
python -m mypy
python -m pytest tests -q --durations=15
```

Ruff checks the public Python code; Mypy checks all maintained package modules,
including unannotated function bodies. CI requires at least **99.5% line coverage
and 99.5% branch coverage**, measured separately with no exclusions. Each test
has a 60-second hang guard.

CI also checks wheel and source-package installation and compiles the mathematical
reference. Its combined native integration job is configured to exercise both
solver I/O modes on a small case once the pinned FiniteFlow revision is published. Full X₃ discovery and larger research computations run locally.

[TESTING](docs/TESTING.md) explains the exact algebraic and failure-path
regressions, coverage commands, gradual type-checking limits and optional
GAP/Singular checks.

## Repository map

| Location | Purpose |
| --- | --- |
| `src/zipper_algebra/` | Maintained Python implementation and command line interface. |
| `tests/` | Algebraic regressions, constructed test inputs and installation checks. |
| `docs/` | Mathematical reference, reproduction guide and testing instructions. |
| `experiments/` | Quotient surveys, representations and structural research recipes. |
| `examples/` | Executable exact computations and the curated X_3 certificate. |
| `singular_run/` | Independent specialized membership check in Singular. |

Retired research material remains available in Git history.

The installed wheel contains the maintained package. The execution layer
validates solver proofs and exact certificate identities. Variant surveys, dimension and trace recipes
retain separate completion and exactness requirements. [INDEX](INDEX.md) maps those sources.

## References

[1]: https://web.math.ucsb.edu/~bigelow/publications/10.pdf
[2]: docs/mathematical-specification/mathematical-specification.pdf

1. Stephen Bigelow, *Braid groups and Iwahori–Hecke algebras*, in *Problems on
   Mapping Class Groups and Related Topics*, Proceedings of Symposia in Pure
   Mathematics **74** (2006), 285–299.
   [Author's copy][1]; [arXiv:math/0505064](https://arxiv.org/abs/math/0505064).
2. Choomno Moos, *The Zipper Algebra*, M.A. thesis, University of California,
   Santa Barbara, September 2026. The [mathematical reference][2] supplies
   the relevant definitions, theorem statement and section references.
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

6. George M. Bergman, *The diamond lemma for ring theory*, Advances in
   Mathematics **29** (1978), 178–218.
   [DOI](https://doi.org/10.1016/0001-8708(78)90010-5);
   [author’s corrections](https://math.berkeley.edu/~gbergman/papers/updates/diamond.html).
7. A. M. Cohen and J. W. Knopper, *GBNP*, version 1.1.0, GAP package.
   [Function manual](https://gap-packages.github.io/gbnp/doc/chap3.html).

8. Wolfram Decker, Gert-Martin Greuel, Gerhard Pfister and Hans Schönemann,
   *Singular*, version 4.4.1, computer algebra system.
   [Official project site](https://www.singular.uni-kl.de/);
   [citation information](https://www.singular.uni-kl.de/index.php/how-to-cite-singular).

9. Aaron Meurer et al., *SymPy: symbolic computing in Python*, PeerJ Computer
   Science **3** (2017), e103.
   [DOI](https://doi.org/10.7717/peerj-cs.103);
   [official citation metadata](https://docs.sympy.org/latest/citing.html).

10. Tiziano Peraro, *FiniteFlow: multivariate functional reconstruction using
    finite fields and dataflow graphs*, Journal of High Energy Physics
    **2019** (2019), no. 7, article 031.
    [DOI](https://doi.org/10.1007/JHEP07(2019)031);
    [arXiv:1905.08019](https://arxiv.org/abs/1905.08019).

## License

The software is available under the [MIT license](LICENSE).
