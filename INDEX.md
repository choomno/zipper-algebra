# Repository map

The installed package provides current presentation and certificate operations.
Research recipes and historical data remain separately identifiable.

## Maintained implementation

All source below is under `src/zipper_algebra/`.

| Location | Purpose |
| --- | --- |
| `coefficients.py`, `polynomial.py` | Exact coefficient backends and free polynomial arithmetic. |
| `presentation.py`, `relations.py` | One authority for the implemented zipper presentations and relation order. |
| `core.py` | Basis parsing, leading-word matching and reduction utilities. |
| `msolve.py` | Solver-format emission and the explicit augmented-presentation operation. |
| `cofactor/atoms.py` | Denominator data and validation. |
| `cofactor/univariate.py`, `cofactor/lift.py`, `cofactor/inference.py` | Exact reconstruction, held-out validation and optional factor inference. |
| `cofactor/proof.py`, `skeleton.py`, `compose.py` | Proof graphs, contextual reductions and symbolic composition. |
| `cofactor/cofactor_io.py`, `cofactor/artifacts.py` | Historical readers and new presentation-bound artifact envelopes. |
| `cofactor/verify.py`, `cofactor/ncpoly.py` | Exact symbolic and degree-bounded product-grid checks. |
| `execution/` | Explicit configuration, validated solver runs, sampling and bounded discovery. |
| `execution/store.py` | Transactional SQLite sample storage and byte-exact inspection exports. |
| `linear.py` | Exact rational rank and validated prime-field rank. |
| `provenance.py` | Source revisions, local changes, exact file hashes and output sidecars. |
| `cli.py` | Installed describe, emit, verify, compose and discover commands. |

`zipper_in_ms.py` retains the documented Python imports and module command as a
small compatibility facade. Historical `tools.*` imports and launchers have
been retired; use the installed package and commands.

## Tests, examples and research recipes

| Location | Purpose |
| --- | --- |
| `tests/unit/` | Mathematical arithmetic, parsing, validation and reconstruction regressions. |
| `tests/integration/` | CLI, installation, coverage and computation-outcome checks. |
| `tests/reproduction/` | Saved certificate and historical proof comparisons. |
| `tests/reference/` | Frozen historical implementation used by those comparisons. |
| `tests/fixtures/legacy_x3/` | Immutable 695-term witness, lift, skeleton, canonical anchor proof and numeric bases. |
| `examples/` | Executable small examples and the curated 129-term X_3 certificate. |
| `experiments/variants/` | Specialized quotient, trace, radical and representation calculations. |
| `experiments/research/` | Slider/lasso minimal-polynomial investigations. |
| `experiments/survey/` | Per-point observations with explicit coefficient domains and completion limits. |
| [experiments/x3](experiments/x3/README.md) | Compact certificate transformations and braid-algebra checks. |

Experimental observations are not generic theorems or automatically certified
finite-dimensional-algebra data. Their unique mathematical computations remain
available while shared arithmetic and presentation construction use the package.

## Scientific data and documentation

`data/research-archive.json` records the pinned historical artifact paths and
SHA-256 hashes. Restore the archive on demand with
`python -m experiments.restore_data --output .research-data`, or select the
`x3-history`, `variants` or `gap` group with `--group`.
The manifest’s `previous_aliases` field records historical duplicate names. Default tests and fresh
production need no restored bulk data.
The [Singular producer](singular_run/README.md) and `archive/gap/` preserve
source programs for independent specialized calculations; their saved logs and
bases are in the optional archive. The reconstruction manuscript source remains
at `reconstruction_theory.tex`.

- [README](README.md): mathematical purpose, installation and principal commands.
- [Mathematical reference](docs/mathematical-specification/main.tex): definitions, hypotheses and exact conclusions.
- [Reproduction guide](docs/X3_REPRODUCTION.md): fresh production, saved-data replay and historical provenance.
- [Workflow](WORKFLOW.md): operations and result statuses.
- [Testing](docs/TESTING.md): mathematical regressions, coverage and CI.
- [Output provenance](docs/PROVENANCE.md): source revisions, dirty builds and byte-bound output records.
