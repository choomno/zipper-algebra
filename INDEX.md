# Repository map

Updated 27 October 2026 for the public package. Historical audits and experiments
remain available, with the status recorded at the time they were written.
Current installation and test instructions are in [README](README.md) and
[TESTING](docs/TESTING.md).

## Maintained package

| Path | Purpose |
| --- | --- |
| `src/zipper_algebra/presentation.py` | Presentation parameters, relation labels and symbolic content hashes. |
| `src/zipper_algebra/zipper_in_ms.py` | Exact coefficient and noncommutative polynomial arithmetic; presentation builder and input emitter. |
| `src/zipper_algebra/core.py` | Basis parsing, leading-word automata and reduction utilities. |
| `src/zipper_algebra/cofactor/atoms.py` | Structured denominator factors and optional inference. |
| `src/zipper_algebra/cofactor/cofactor_io.py` | Certificate data, validation, targets and legacy readers. |
| `src/zipper_algebra/cofactor/proof.py` | Basis keys, proof references and expansion. |
| `src/zipper_algebra/cofactor/skeleton.py` | Recorded contextual reduction steps. |
| `src/zipper_algebra/cofactor/lift.py` | Exact interpolation and candidate reconstruction. |
| `src/zipper_algebra/cofactor/compose.py` | Expansion and composition preserving both word contexts and repeated occurrences. |
| `src/zipper_algebra/cofactor/verify.py` | Symbolic and degree-bounded grid checks of free-algebra identities. |
| `tests/` | Default regression and package tests. |
| `tests/fixtures/legacy_x3/` | Immutable inputs for historical compatibility comparisons, with hashes. |

The modules under `tools/` that correspond to these package modules are
compatibility entry points. They share the maintained implementation rather than
maintaining a second copy of the arithmetic.

## Experimental execution and research

| Path | Purpose and status |
| --- | --- |
| `tools/cofactor/pipeline.py` | Solver sampling, reconstruction and verification driver; execution remains experimental. |
| `tools/cofactor/paths.py`, `samples.py`, `tools/cache.py` | Historical workspace configuration, sample collection and caching. |
| `tools/cofactor/{compress_cofactor,simplify_cofactor,reduce_to_three_rels,rename_to_ZT,verify_ZT_and_ideal_equality}.py` | Specialized historical transformations and checks. Their ambient algebra and dependencies require individual review. |
| `tools/auxiliary/` | Diagnostics for earlier certificates. |
| `tools/variants/` | Dimension, trace, radical, inclusion and representation experiments. |
| `tools/research/`, `tools/survey/` | Parameter surveys and candidate higher relations. |
| `tools/tests.py` | Older solver integration suite, separate from the default CI tests. |
| `gap_scripts/`, `singular_run/` | Independent computational routes and saved outputs. |
| `reproduce_X3_zero.sh` | Historical solver workflow entry point. |

Experimental script names, saved process output and agreement at several
specializations do not establish a completed generic basis or a generic theorem.
Use the evidence requirements in the mathematical specification when interpreting
those records.

## Documentation

| Document | Purpose |
| --- | --- |
| [Mathematical specification](docs/mathematical-specification/main.tex) | Objects, conventions, relation formulas, certificate domains and verification conclusions. |
| [WORKFLOW](WORKFLOW.md) | Current package and certificate workflow. |
| [TESTING](docs/TESTING.md) | Default tests, package checks and optional integrations. |
| [PUBLIC_SCOPE](docs/PUBLIC_SCOPE.md) | Public contents, completed cleanup and remaining work. |
| [RESEARCH_AGENDA](RESEARCH_AGENDA.md) | Broader mathematical questions and publication objectives. |
| [AUDIT_spec_vs_generic](AUDIT_spec_vs_generic.md), [LEDGER_claims](LEDGER_claims.md) | Historical claim-by-claim evidence and gaps. |
| [GENERALITY_CENSUS](GENERALITY_CENSUS.md) | Historical implementation audit, with dated paths and line numbers. |
| [REPTHEORY_STRATEGY](REPTHEORY_STRATEGY.md) | Historical representation-theoretic programme and dependency assessment. |
| [SURVEY_WORKFLOW](SURVEY_WORKFLOW.md), [NIGHT_LOG](NIGHT_LOG.md) | Historical plans, experiments and candidate patterns. |
| [SOURCE_IMPORT](SOURCE_IMPORT.md) | Provenance and extent of the private source recovery. |

Older document versions and original paths remain in Git history and the
preserved private baseline. They are evidence of earlier development, not current
instructions to run long computations or modify the public repository.
