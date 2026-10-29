# Source map

The maintained implementation is installed from `src/zipper_algebra/`.

| Area | Purpose |
| --- | --- |
| `presentation.py`, `relations.py`, `untwisting.py` | Algebra presentations, ordered generators, relations and targets. |
| `coefficients.py`, `polynomial.py`, `msolve.py` | Coefficient domains, free polynomials and solver-input emission. |
| `core.py`, `linear.py` | Basis parsing, standard-word automata, reduction and exact linear algebra. |
| `cofactor/proof.py`, `cofactor/skeleton.py` | Proof graphs, contextual reductions and dependency closure. |
| `cofactor/inference.py`, `cofactor/univariate.py`, `cofactor/lift.py` | Denominator inference and exact coefficient reconstruction. |
| `cofactor/compose.py`, `cofactor/cofactor_io.py`, `cofactor/artifacts.py` | Certificate composition, text interchange and validated mathematical records. |
| `cofactor/ncpoly.py`, `cofactor/verify.py` | Symbolic and degree-bounded grid checks of free-algebra identities. |
| `verification/gap.py` | Export standalone generic cofactor identities for independent expansion in GAP. |
| `execution/` | Solver configuration, execution, sampling, SQLite storage and verified publication. |
| `provenance.py` | Source revisions, local changes and input/source hashes for outputs. |
| `cli.py` | Installed commands and their exit statuses. |

The coefficient and polynomial modules provide the canonical Python imports.
`zipper-algebra` is the command line interface. New bound mathematical records
use version-2 envelopes; older JSON envelope versions are rejected.

## Tests and examples

`tests/unit/` and the top-level test modules exercise algebraic identities,
reconstruction, proof graphs, invalid inputs and execution boundaries.
`tests/integration/` checks installed interfaces, distributions, example output
and publication outcomes. Test data are constructed locally by those tests.
`examples/x3/certificate.txt` is the retained 129-term certificate, with commands
and expected results in its adjacent README.

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
- [Reproduction guide](docs/X3_REPRODUCTION.md): fresh production and certificate checks.
- [Workflow](WORKFLOW.md): supported operations and result statuses.
- [Testing](docs/TESTING.md): regression tests, coverage and CI.
- [Output provenance](docs/PROVENANCE.md): source revisions and byte-bound output records.
