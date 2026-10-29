# Tests and continuous integration

Use Linux with Python 3.11 or newer and install the package before running source tests:

```sh
python -m pip install -e '.[dev,inference]'
python -m ruff check src tests experiments examples ci setup.py
python -m pytest tests
```

The suite uses pytest 9, native pytest fixtures and parametrization.
Standard-library `unittest.mock` remains available for mocking process boundaries.

The pinned Ruff version checks package code, tests, research recipes, examples, CI
scripts and the build script. The `F` and `E9` rule families detect unused imports
and local variables, undefined names, repeated definitions and fatal syntax errors.
The same command runs in the CI coverage job before the tests.

The suite covers unit and integration tests, including installed interfaces,
execution boundaries and the 129-term certificate under `examples/x3/`. It
constructs its remaining inputs locally and does not replay archived research
collections.

## What the tests establish

Small tests check independently specified identities, coefficient domains,
parsing, invalid data, noncommuting contexts, repeated occurrences, interpolation
rank and held-out checks.

Process-boundary tests use an explicitly identified deterministic producer to
exercise timeout, unsuccessful exit, malformed/truncated proof output, cache
identity, interrupted publication and resumption. They do not pretend to be
Gröbner-basis computations. Native f4ncgb integration is a separate CI job,
described below.

Source-provenance tests use temporary Git repositories and installed-package
copies to check clean, staged, unstaged and unversioned states, build origins
and later source changes. SQLite tests check concurrent readers and writers,
rollback, crash recovery, corruption rejection, historical generations and
byte-exact exports. Native-adapter tests check database job isolation, process
failure and cleanup, and rejection of changed or incomplete solver output.
Cached samples are independently rechecked as polynomial proof identities.

The current example is checked as a free-algebra identity using symbolic
accumulation and a degree-bounded product grid. Proof-graph regressions use
independently specified polynomial identities, including forward references.

## Mathematical and failure-path evidence

Tests use literal noncommuting polynomial identities, exact rational and prime-field
examples, deliberate false certificates, and storage/publication failures. Expected
values are specified independently rather than obtained by calling the operation
under test a second time.

Degree-bound regressions include false identities with target, numerator and common
denominator factors vanishing on an undersized grid. A mixed-denominator identity
and its perturbed version test the clearing complements. These guard the reason
that the finite grid check establishes an exact polynomial identity.

Controlled local fault checks have exercised wrong inverses, commuted words, ignored
held-out values, bypassed proof checks, altered metadata and stale storage selections.
They found three missing degree-bound adversaries, now covered by explicit regressions.
This is evidence for those selected faults, not an exhaustive correctness argument.

## Runtime bounds

Use `python -m pytest tests/unit -q --durations=10` for quick algebraic checks.
The full software suite also exercises CLI, database and package boundaries.
Each test, including its fixtures, has a 60-second Linux signal timeout; this is
a hang guard, not a performance target. Review `--durations=15` output to find
slow cases. Solver processes and native GAP checks also have explicit subprocess
timeouts.

CI uses small fixed inputs and a bounded two-strand native solver integration.
Larger presentation searches and X₃ discovery remain optional local computations;
termination of a general noncommutative basis computation is not assumed.

## Maintained package coverage

The requirement is **at least 99.5% line coverage and 99.5% branch coverage**,
checked separately before rounding. Every Python file under `src/zipper_algebra/`
is measured, including the execution layer. No source lines or partial branches
are excluded. Research recipes have their own regression checks and are not
part of the maintained package measurement.

```sh
python -m pip install -e '.[dev,inference]' coverage==7.17.0
python -m coverage erase
python -m coverage run -m pytest tests
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
its wheel, then runs the current tests and checks the installed certificate interface.
Use a disposable environment: this check replaces its installed package.

The wheel contains only package code and distribution metadata. The source
archive includes reader documentation, tests, examples and
experimental recipes needed by those tests. It omits generated computation data,
private working records, caches, binary matrix pickles and machine configuration.

Both distributions embed the source build origin. The source-archive check
rebuilds below an unrelated Git checkout and requires that origin to survive;
package-source changes are rejected when they would inherit stale metadata.

## GitHub Actions

The platform matrix checks Linux on Python 3.11, 3.12, 3.13 and 3.14. It runs the unit and integration tests and isolated wheel checks; the Linux 3.13
job also rebuilds the source archive. The full suite, including integration
comparisons, runs in the coverage job with the real optional Sympy dependency.

The native solver job builds the pinned f4ncgb source with SQLite enabled on Ubuntu,
installs the zipper-algebra wheel, and exercises a small two-generator test
presentation with `X_2=0` explicitly imposed through both file and native SQLite
interfaces. It checks native input/output,
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

The inventory covers the current source, certificate example and test/build
support. The [reproduction guide](X3_REPRODUCTION.md) identifies the inputs and
implementation used by each current route. Retired material remains in Git history.


## Independent algebra-system checks

The GAP exporter has small tests for exact coefficient emission, generator and
context ordering, atom polynomials, invalid inputs and provenance. The
[reproduction guide](X3_REPRODUCTION.md#independent-algebra-system-checks) gives
the local GAP check of the generic certificate. This optional external check
performs identity expansion, with no solver search.

The [Singular calculation](../singular_run/README.md) is an optional local
cross-check of specialized X₃ membership over the rationals. Its degree bound
and parameter specialization are explicit; generic cofactor verification uses
the symbolic and exact product-grid checks described above.
