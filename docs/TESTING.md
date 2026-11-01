# Tests and continuous integration

Use Linux with Python 3.11 or newer and install the package before running source tests:

```sh
python -m pip install -e '.[dev,inference,symbolic]'
python -m ruff check src tests experiments examples ci setup.py
python -m mypy
python -m pytest tests
```

The suite uses pytest 9, native pytest fixtures and parametrization.
Standard-library `unittest.mock` remains available for mocking process boundaries.

The pinned Ruff version checks package code, tests, research recipes, examples, CI
scripts and the build script. The `F` and `E9` rule families detect unused imports
and local variables, undefined names, repeated definitions and fatal syntax errors.
The same command runs in the CI coverage job before the tests.

Mypy checks every module under `src/zipper_algebra`, including the bodies of
unannotated functions, and is pinned in the development extra. CI runs the same
`python -m mypy` command. Collection and scalar annotations describe exact
coefficients, parameter points, proof terms and stored artifacts. The optional
Sympy dependency lacks type information and is treated as an external untyped
boundary. Other unannotated interfaces can still carry `Any`, so this is gradual
type checking; it complements the exact arithmetic and failure-path tests.

The suite covers unit and integration tests, including installed interfaces,
execution boundaries and the 129-term certificate under `examples/x3/`. It
constructs its remaining inputs locally and does not replay archived research
collections.

## What the tests establish

Discovery tests check fresh proof expansion, generic matrix formation,
normalization of kernel vectors, rejection of insufficient supports and false
coefficients, and immutable output publication. The production route uses one
anchor proof rather than sample-by-sample coefficient fitting.

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

Solver-reuse tests run distinct targets against one presentation and check that
only one external calculation occurs, including through a shared cache and
separate output roots. Invalid targets fail before cache lookup; rehashed false
proofs remain unusable. Changes to input, solver configuration, executable bytes
or source provenance prevent reuse. Earlier target-bound records remain readable
without being treated as current cache entries.

The current example is checked as a free-algebra identity using symbolic
accumulation and a degree-bounded product grid. Proof-graph regressions use
independently specified polynomial identities, including forward references.

## Mathematical and failure-path evidence

Tests use literal noncommuting polynomial identities, exact rational and prime-field
examples, deliberate false certificates, and storage/publication failures. Expected
values are specified independently rather than obtained by calling the operation
under test a second time.

The parameter-field tests compare independently expanded three-parameter
polynomials, reject mixed domains even for zero, and test denominator zeros and
characteristic changes. Small [acceptance controls](PARAMETER_FIELDS.md) combine
free-algebra cofactor identities with cyclic companion representations to prove
matching upper and lower dimension bounds. They include repeated roots, a
nonzero square-zero element, and two different exceptional specializations.

Completed-basis tests check independent finite-algebra answers, both ideal
inclusions, inclusion ambiguities, equal leading words and proper overlaps,
including self-overlaps. Deliberately incomplete rewrite systems must fail even
when every displayed basis polynomial has an honest input derivation. Proof
graph tests reject cycles, missing references and invalid references with zero
coefficients; a heavily shared acyclic graph checks that validation uses local
identities without expanding the graph into a cofactor list. Domain mismatches,
the unit ideal, optional rule pruning, work-budget exhaustion and incomplete
standard-word enumeration have distinct expected outcomes. Numeric solver-text
decoding is tested separately from the completion verdict. See the
[completed-basis guide](GROEBNER_VERIFICATION.md).

Finite word-basis tests specify matrix-unit and parameter-field controls
independently of a Gröbner producer. They reject a cyclic representation of a
proper quotient, dependent proposed words with honest spanning proofs, corrupted
action columns, and defining relations that fail away from the cyclic vector.
Zero-algebra tests require an explicit unit-ideal proof. Further checks cover
noncommuting products, sparse action orientation, strict graph and JSON data,
mixed parents including zero, immutable coordinates, and independent operation
limits. Working-field tests compare exact fraction arithmetic, cancellation,
parameter order and prime-field conversions while checking that literal
coefficient data are preserved. Presentation identity remains bound to those
original literal records. See
[finite word-basis verification](FINITE_ALGEBRAS.md).

Degree-bound regressions include false identities with target, numerator and common
denominator factors vanishing on an undersized grid. A mixed-denominator identity
and its perturbed version test the clearing complements. These guard the reason
that the finite grid check establishes an exact polynomial identity.

Controlled local fault checks have exercised wrong inverses, commuted words, ignored
held-out values, bypassed proof checks, altered metadata and stale storage selections.
They found three missing degree-bound adversaries, now covered by explicit regressions.
This is evidence for those selected faults, not an exhaustive correctness argument.

The local four-strand finite-basis reproduction also exercises the maintained
verifier on freshly generated generic data: 1,116 action identities, 1,563
proof nodes, 5,766 original matrix-relation columns and 186 cyclic words.
It remains a local research check, outside CI; the
[finite-algebra guide](FINITE_ALGEBRAS.md) gives the full commands and recorded
counts. Small independent controls supply the routine regressions.

Minimal-polynomial tests use independently specified Krylov columns and modular
minor inverses for matrix algebras, upper triangular matrices, nilpotent and
quadratic algebras, and generic and prime coefficient fields. They reject a
nonminimal annihilator even when its supplied kernel is valid, distinguish a
relative orbit from the full element, and check scaled alternative kernel layouts,
zero algebras, sparse JSON, field/basis binding and budgets. Separate producer tests
exercise exact prefix growth, truncation, native options and stale-answer controls.
Small actual native tests check generic parameter fields, relative nilpotence and
bad-prime recovery. The 186-dimensional four-strand example remains local; see the
[minimal-polynomial guide](MINIMAL_POLYNOMIALS.md).

Generic-production tests use literal noncommutative algebra presentations,
independent small exact kernel oracles, shared and forward derivations, and
separate native reconstruction controls. They check rational parameter points,
compacted zero-input mappings, canceled generic tails, lost leading words,
anchor-compatible kernel combinations, poles, incomplete bases, false proofs,
and independent action certificates. Publication tests corrupt saved data and
exercise interrupted attempts and concurrent immutable generations. Tiny native
controls check both file and SQLite transport through an installed package;
three-strand research examples stay outside CI. See the
[generic production guide](GENERIC_PRODUCTION.md).

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
python -m pip install -e '.[dev,inference,symbolic]' coverage==7.17.0
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

The native solver job is configured to build the pinned f4ncgb source with SQLite
enabled on Ubuntu, install the zipper-algebra wheel, and exercise a small two-generator test
presentation with `X_2=0` explicitly imposed through both file and native SQLite
interfaces. Its checks cover native input/output,
proof identities, FiniteFlow reconstruction, saved kernel witnesses, immutable
result bundles and both certificate checks outside this workspace. Full four-strand discovery and
other research computations run locally. The curated 129-term certificate was produced by the earlier sampled fitting
route, with both exact checks passing; see the [reproduction guide](X3_REPRODUCTION.md). The workflow
and `ci/check_solver.py` record the pinned source and commands.

The combined native solver/reconstruction job is also configured to build the exact
maintained FiniteFlow revision pinned in the workflow and run the small installed
matrix reconstruction tests and both discovery I/O modes outside the checkout. These tests check parameterized and rational
systems, prime and shift recovery, repeated calls and serialized native sessions,
then independently check the complete kernel identities. The workflow fetches the published native revision identified in the
[reconstruction guide](RECONSTRUCTION.md).
The fork also tests modular failures, strict original-system validation and
relearning through its C API, with and without FLINT.

Unit tests separately reject malformed matrix/witness data, false rank bounds,
nonzero generic residuals, undefined specializations, failed native statuses,
incomplete outputs, and incorrect learned metadata. Controlled native faults
exercise these rejection paths; the installed-binding tests exercise the actual
solver. The larger X3 and 186-dimensional matrix checks remain local computations.

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

The [four-strand square-zero example](FOUR_STRAND_SQUARE_ZERO.md) produces a
generic derivation certificate with Singular, then uses GAP/GBNP to check the
original presentation, every local proof identity, completion, dimension and
the square-zero ideal calculation over `Q(q,t)`. This research run stays local.
Its bounded parser, proof-mutation and process tests run with the software
suite; the small GAP tests skip when GAP/GBNP is unavailable.

The GAP exporter has small tests for exact coefficient emission, generator and
context ordering, atom polynomials, invalid inputs and provenance. The
[reproduction guide](X3_REPRODUCTION.md#independent-algebra-system-checks) gives
the local GAP check of the generic certificate. This optional external check
performs identity expansion, with no solver search.

The [Singular calculation](../singular_run/README.md) is an optional local
cross-check of specialized X₃ membership over the rationals. Its degree bound
and parameter specialization are explicit; generic cofactor verification uses
the symbolic and exact product-grid checks described above.
