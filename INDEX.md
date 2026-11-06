# Source map

The maintained implementation is installed from `src/zipper_algebra/`.

| Area | Purpose |
| --- | --- |
| `presentation.py`, `relations.py`, `untwisting.py` | Algebra presentations, ordered generators, relations and targets. |
| `coefficients.py`, `polynomial.py`, `msolve.py` | Coefficient domains, free polynomials and solver-input emission. |
| `fields.py`, `algebra.py` | Arbitrary parameter fields, domain-bound free algebras, literal presentations and exact cofactor identities. |
| `matrices.py`, `reconstruction.py`, `verification/linear.py` | Exact sparse matrices, FiniteFlow kernel reconstruction, and independent generic rank/kernel certificates. |
| `invariants.py`, `verification/krylov.py` | Certified element and relative-orbit minimal polynomials, with independent Krylov/rank verification. |
| `core.py`, `linear.py` | Basis parsing, standard-word automata, reduction and exact linear algebra. |
| `cofactor/proof.py`, `cofactor/skeleton.py` | Proof graphs, contextual reductions and dependency closure. |
| `cofactor/discovery.py` | Proof-graph word contexts, generic cofactor matrices and exactly checked FiniteFlow solutions. |
| `cofactor/inference.py`, `cofactor/univariate.py`, `cofactor/lift.py` | Standalone sampled denominator inference and coefficient fitting; separate from current discovery. |
| `cofactor/compose.py`, `cofactor/cofactor_io.py`, `cofactor/artifacts.py` | Certificate composition, text interchange and validated mathematical records. |
| `cofactor/ncpoly.py`, `cofactor/verify.py` | Symbolic and degree-bounded grid checks of free-algebra identities. |
| `verification/finite.py` | Finite word-basis certificates, verified sparse regular actions, coordinates and products. |
| `verification/derivations.py` | Shared local acyclic proofs of exact two-sided ideal identities. |
| `verification/working_field.py` | Internal normalized arithmetic for repeated exact parameter-field calculations. |
| `verification/groebner.py` | Exact derivation checks, both ideal inclusions, complete critical-ambiguity checks and certified normal forms. |
| `verification/f4ncgb.py` | Decode numeric solver basis/proof text into a candidate domain-bound certificate. |
| `verification/gap.py` | Export standalone generic cofactor identities for independent expansion in GAP. |
| `zipping_algorithm/`, `native/zipper_engine/` | Literal-word zipper engines, exact collection, packed C++ execution, and independent Python replay. See the [native engine guide](docs/ZIPPING_NATIVE.md). |
| `execution/` | Solver configuration, execution, SQLite storage, native dependency identities and verified publication. |
| `provenance.py` | Source revisions, local changes and input/source hashes for outputs. |
| `cli.py` | Installed commands and their exit statuses. |

The coefficient and polynomial modules provide the canonical Python imports.
`zipper-algebra` is the command line interface. Execution records use version-2
envelopes; older JSON envelope versions are rejected. The domain-bound algebra
interface has separate, versioned JSON formats described in the
[parameter-field guide](docs/PARAMETER_FIELDS.md). The
[completed-basis guide](docs/GROEBNER_VERIFICATION.md) describes proof data,
completion criteria, normal forms and resource limits.

## Tests and examples

`tests/unit/` and the top-level test modules exercise algebraic identities,
reconstruction, proof graphs, invalid inputs and execution boundaries.
`tests/integration/` checks installed interfaces, distributions, example output
and publication outcomes. Test data are constructed locally by those tests.
`examples/x3/certificate.txt` is the retained 129-term certificate, with commands
and expected results in its adjacent README.
`examples/polynomial_quotients.py` supplies small acceptance controls with
independently known answers: companion representations, cofactor identities and
exceptional fibers. These exercise the parameter-field interface.

The [four-strand square-zero example](docs/FOUR_STRAND_SQUARE_ZERO.md) produces
a candidate generic basis and local derivations using Singular; the adjacent GAP
checker independently verifies the original presentation, completion and the
nonzero square-zero ideal. Its generated outputs remain under ignored `runs/`.

The [Singular check](singular_run/README.md) independently tests specialized
X₃ membership over the rationals. It supplies a separate computation, with its
degree bound and parameter specialization stated explicitly.

## Exploratory algorithms

`experiments/research/` contains lasso-twist polynomial searches.
`experiments/survey/` contains quotient sampling and structural calculations.
`experiments/variants/` contains dimension, trace, representation and ideal
recipes. These computations have their own completion and exactness hypotheses;
a sampled observation does not establish a generic theorem.

Generated computation outputs live in ignored run directories. Retired scripts,
datasets and interfaces remain available in Git history.

- [README](README.md): mathematical purpose, installation and principal commands.
- [Mathematical reference](docs/mathematical-specification/mathematical-specification.pdf): definitions, hypotheses and verification conclusions.
- [Four-strand square-zero computation](docs/FOUR_STRAND_SQUARE_ZERO.md): generic dimension, certified basis and a nonzero square-zero ideal.
- [Minimal polynomials](docs/MINIMAL_POLYNOMIALS.md): exact element and orbit relations with minimality witnesses.
- [Generic production](docs/GENERIC_PRODUCTION.md): arbitrary parameter fields, fresh f4ncgb acquisition, shared proof reconstruction and certified actions.
- [Certified reconstruction](docs/RECONSTRUCTION.md): native setup, exact matrix/kernel schemas, sampling controls and verification.
- [Finite algebras](docs/FINITE_ALGEBRAS.md): word-basis witnesses, sparse actions, coordinates and exact products.
- [Completed-basis verification](docs/GROEBNER_VERIFICATION.md): derivations, ideal equality, critical ambiguities and normal forms.
- [Reproduction guide](docs/X3_REPRODUCTION.md): fresh production and certificate checks.
- [Workflow](WORKFLOW.md): supported operations and result statuses.
- [Testing](docs/TESTING.md): regression tests, coverage and CI.
- [Output provenance](docs/PROVENANCE.md): source revisions and byte-bound output records.
