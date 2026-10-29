# Tests and continuous integration

Use Python 3.11 or newer and install the package before running source tests:

```sh
python -m pip install -e '.[inference]'
python ci/run_tests.py --suite all
```

Ordinary unittest discovery remains supported:

```sh
python -m unittest discover -s tests -v
```

The complete suite includes unit, integration and reproduction tests, including
all 16 historical proof comparisons. Separate selections are:

```sh
python ci/run_tests.py --suite fast
python ci/run_tests.py --suite reproduction
```

These selections partition the full suite without overlap or omissions. The fast
selection includes all non-reproduction cases, including execution-boundary,
CLI and installation regressions. The reproduction selection checks the 695-term
certificate fixtures and frozen historical implementation. The curated 129-term certificate lives under
`examples/x3/`. No default test selection requires the optional bulk archive.

## What the tests establish

Small tests check independently specified identities, coefficient domains,
parsing, invalid data, noncommuting contexts, repeated occurrences, interpolation
rank and held-out checks. They are separate from compatibility with the large
retained witness.

Process-boundary tests use an explicitly identified deterministic producer to
exercise timeout, unsuccessful exit, malformed/truncated proof output, cache
identity, interrupted publication and resumption. They do not pretend to be
Gröbner-basis computations. Native f4ncgb integration is a separate CI job,
described below.

Historical comparisons use immutable fixtures under `tests/fixtures/legacy_x3/`
and a frozen implementation under `tests/reference/`. That implementation shares
the coefficient kernel with the package. Its purpose is historical equivalence;
the separate symbolic and degree-grid checks have their own stated independence.
The retained 695-term witness and its numeric bases are computational regression
cases. The current example is checked as a free-algebra identity using both
symbolic accumulation and a degree-bounded product grid.

## Maintained package coverage

The requirement is **at least 99.5% line coverage and 99.5% branch coverage**,
checked separately before rounding. Every Python file under `src/zipper_algebra/`
is measured, including the execution layer. No source lines or partial branches
are excluded. Research recipes have their own regression checks and are not
part of the maintained package measurement.

```sh
python -m pip install -e '.[inference]' coverage==7.17.0
python -m coverage erase
python -m coverage run ci/run_tests.py --suite all
python -m coverage combine
python -m coverage json -o coverage.json
python ci/check_coverage.py coverage.json
```

The configuration measures Python subprocesses as well as the parent process.
The gate requires the complete maintained file list and rejects exclusions.
It compares raw line and branch ratios with the threshold. The CI log records
the current counts; a percentage displayed after rounding is not enough.

For an annotated local report, run `python -m coverage html` after combining
the data. Coverage describes execution; the mathematical assertions and
independent expected results remain essential.

## Distribution checks

```sh
python -m pip install build
python -m build
python ci/check_distribution.py dist
```

Install the built wheel in a fresh environment and run
`python ci/check_installed.py` from the checkout. It imports the package and
executes commands outside the source tree in isolated Python processes.

`python ci/check_sdist.py dist` extracts the source archive, rebuilds and installs
its wheel, then runs the fast suite and a full retained-certificate recomposition.
Use a disposable environment: this check replaces its installed package.

The wheel contains only package code and distribution metadata. The source
archive includes reader documentation, tests, examples, selected fixtures and
experimental recipes needed by those tests. It omits the full research archive,
private working records, caches, binary matrix pickles and machine configuration.

## GitHub Actions

The platform matrix checks Linux on Python 3.11, 3.12, 3.13 and 3.14, and Windows
on Python 3.13. It runs the fast suite and isolated wheel checks; the Linux 3.13
job also rebuilds the source archive. The full suite, including historical
comparisons, runs in the coverage job with the real optional Sympy dependency.

The native solver job builds the pinned parser-fixed f4ncgb source on Ubuntu,
installs the zipper-algebra wheel, and exercises a small two-generator test
presentation with `X_2=0` explicitly imposed. It checks native input/output,
proof identities, coefficient reconstruction, the result bundle and both
certificate checks outside this workspace. Full four-strand discovery and
other research computations run locally. A local full four-strand run produced
the curated 129-term certificate from 100 fresh parameter points, with both exact
checks passing; see the [reproduction guide](X3_REPRODUCTION.md). The workflow
and `ci/check_solver.py` record the pinned source and commands.

The documentation job validates local links, deterministic inventory
regeneration and recorded hashes, then compiles the mathematical reference.
Workflow run pages record actual hosted results; configured jobs and local
checks alone do not establish a hosted pass.

## Mathematical documentation and audit inventory

```sh
python ci/check_docs.py
python ci/build_x3_manifest.py --check
python ci/check_x3_manifest.py
python ci/build_docs.py
```

The compiler runs two passes, checks references and overfull boxes, and replaces
the output PDF only after success. The committed
[mathematical reference](mathematical-specification/mathematical-specification.pdf)
and its source are included in the source archive. CI also retains its newly
compiled PDF as an artifact. PDF bytes can vary across compiler/font versions.

The inventory covers the curated current source, example, fixtures and test/build
support. The [reproduction guide](X3_REPRODUCTION.md) identifies the inputs and
implementation used by each current route. Retired material remains in Git history.

