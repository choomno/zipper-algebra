# Tests and continuous integration

Use Python 3.11 or newer. Install the package before running the source tests:

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
python -m unittest tools.cofactor.tests_proof_skeleton -v
```

The maintained tests need no solver, personal configuration, sibling research
folders, populated cache or network access after installation. They cover exact
and modular arithmetic, basis parsing, presentation validation, proof expansion,
certificate validation, contextual composition, repeated reduction steps, exact
interpolation, held-out validation, inconclusive outcomes, package imports and
CLI exit codes. Tests use temporary directories for output.

The historical comparisons use the archived implementation and immutable
fixtures under `tests/fixtures/legacy_x3/`. Their hashes and original records are
listed in the adjacent manifest. The retained large free-algebra certificate is
a compatibility fixture; tests do not identify it as the latest mathematical proof.

## Maintained package coverage

The target is **at least 99.5% line coverage and 99.5% branch coverage**, checked
separately before rounding. The measured source is every Python file in
`src/zipper_algebra/`. Experimental scripts and archived implementations have
separate test requirements and are outside this package measurement.

The complete local suite measured **99.90% lines (3088/3091)** and
**99.72% branches (1424/1428)** across 231 maintained tests and 16 historical
proof comparisons. No source lines or branches are excluded. The tests cover
exact identities, coefficient domains, parsing and validation, reduction,
reconstruction, noncommuting contexts, CLI behaviour and optional inference.
These percentages describe exercised code; the mathematical assertions and
independent expected results in the tests remain essential.

To measure the complete default suite and historical proof comparisons:

```sh
python -m pip install -e '.[inference]' coverage==7.17.0
python -m coverage erase
python -m coverage run -m unittest discover -s tests -v
python -m coverage run -m unittest tools.cofactor.tests_proof_skeleton -v
python -m coverage combine
python -m coverage json -o coverage.json
python ci/check_coverage.py coverage.json
```

The checked-in configuration enables branch measurement and Python subprocess
coverage, including isolated CLI invocations. The checker requires the complete
maintained file list, rejects excluded source lines, and enforces both thresholds.
The configuration clears both line-exclusion and partial-branch-exclusion patterns.
The CI coverage job installs the real optional inference dependency and runs
these checks. For a local annotated report, run `python -m coverage html` after
combining the data, then open `htmlcov/index.html`.

## Distribution checks

Build both formats, then inspect their contents:

```sh
python -m pip install build
python -m build
python ci/check_distribution.py dist
```

In a fresh environment, install the built wheel and run
`python ci/check_installed.py` from the checkout. That script imports the package
and runs commands from a separate temporary working directory in isolated Python
processes. It rejects an editable installation for this check. It distinguishes
a valid certificate, a false identity and malformed input.

`python ci/check_sdist.py dist` extracts the source distribution, rebuilds its
wheel and runs the source regressions. Run it in a disposable environment: it
replaces that environment's installed development package with the rebuilt wheel.

The wheel contains only the maintained package and distribution metadata. The
source archive includes public documentation, tests and selected regression
fixtures. Private histories, recovered manuscripts, machine configuration,
caches and experimental result collections are excluded.

## GitHub Actions

The workflow runs on pushes, pull requests and manual dispatch. Its matrix
checks Linux on Python 3.11, 3.12, 3.13 and 3.14, plus Windows on Python 3.13.
It builds the wheel and source archive, tests the installed wheel, runs the
regressions and rebuilds the source archive. A separate job measures maintained-package coverage and checks the optional
Sympy factorization extra. Documentation checks validate local links and compile
the mathematical specification.

The compiled [mathematical reference](mathematical-specification/mathematical-specification.pdf)
is committed beside its source and included in the source distribution. With
`pdflatex` installed, rebuild it from the repository root using:

```sh
python ci/build_docs.py
```

The helper runs two passes, keeps auxiliary files in a temporary directory and
rejects unresolved references or overfull boxes before replacing the PDF. It
uses the source revision's timestamp when Git history is available; PDF bytes
can still differ across TeX and font versions. CI rebuilds the reference and
retains it as a downloadable artifact on the workflow run.

Jobs have time limits, read-only repository permissions and no publication step.
The workflow does not require project secrets or an algebra solver. A local test
pass and a GitHub-hosted run are separate evidence; workflow run pages record the
hosted results.

## Experimental integrations

The older `python -m tools.tests`, solver sampling, dimension and trace scripts
retain separate dependency and completion requirements. They are outside the
maintained CI boundary. Large discovery computations need explicit tool versions
and completion evidence before their mathematical conclusions can be accepted.

See the [mathematical specification](mathematical-specification/main.tex) for
the coefficient domains, assumptions and conclusions of the maintained operations.
