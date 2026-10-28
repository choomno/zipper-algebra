# Reproducing the three-strand slider cofactor result

The maintained route checks an explicit cofactor identity for the three-strand
slider, named `X_3` by the command line interface. It uses the four-strand
canonical presentation and the retained 695-term free-algebra certificate.
The commands below require a checkout of this repository, Python 3.11 or newer,
and no external algebra solver.

The [mathematical reference](mathematical-specification/mathematical-specification.pdf)
defines the sliders, presentation, coefficient field, denominator hypotheses
and exact conclusion. In particular, the identity is over the rational function
field in the two parameters; specialization requires the displayed denominators
to remain nonzero. The reference also explains the relationship with the
thesis's result and its BMW interpretation.

## 1. Install and verify the saved identity

Run these commands from the repository root:

```sh
python -m pip install .
python ci/check_x3_manifest.py
python -m zipper_algebra describe --n 4 --json
python -m zipper_algebra verify \
  logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt \
  --n 4 --target X_3 --pit --json
```

The verification command must exit with code `0`. Its JSON has top-level
`"pass": true`, an `exact` report with `"pass": true` and `"tuples": 695`,
and a `pit` report with `"pass": true`. The first check compares the cleared
identity symbolically; the second checks a product grid whose degree bounds
make it an exact identity test. The two checks share the presentation and input
validation, while using separate arithmetic for accumulating the identity.
The equivalent installed entry point is `zipper-algebra verify`.

This verifies the saved cofactor result. It requires the certificate and the
installed presentation and verifier modules. The sampled bases, old solver
executables, cache, thesis source and compression scripts are unnecessary for
this step. Exit codes `1` and `2` distinguish a false identity from invalid
input; see [WORKFLOW](../WORKFLOW.md).

## 2. Recompose the certificate from the saved symbolic proof

The following test starts from the immutable symbolic lift and reduction
skeleton, expands their basis references, composes the contextual terms, and
compares the resulting 695-term certificate body byte for byte with the saved
certificate:

```sh
python -m unittest \
  tools.cofactor.tests_compose_io_verify.TestCompose.test_composed_body_is_byte_identical_to_the_shipped_certificate \
  -v
```

The expected result is `Ran 1 test` followed by `OK`. Test outputs are temporary;
the saved research files are not overwritten. This establishes reproducibility
of composition from the supplied symbolic coefficients. It does not rerun the
sampling or interpolation that originally produced them.

To check the preceding proof-expansion and reduction stages against the archived
implementation, run the 16 comparisons:

```sh
python -m unittest tools.cofactor.tests_proof_skeleton -v
```

These load the retained anchor basis and immutable nonexpanded proof, check
proof expansion and reduction, and compare with the historical implementation.
The expected result is `Ran 16 tests` followed by `OK`. The broader suite is
`python -m unittest discover -s tests -v`; its commands and scope are described
in [TESTING](TESTING.md).

## 3. Follow the compact certificate's provenance

The historical transformation sequence is retained under
`logs_and_traces/f4ncgb/`:

| File | Contribution |
| --- | --- |
| `cofactor-X3-n4-Qqt-symbolic.txt` | The 695-term identity checked above, with explicit braid and inverse relation terms in the free algebra. |
| `cofactor-X3-n4-Qqt-compressed.txt` | Substitution of derivations for dependent input relations; this intermediate expression has 1077 terms. |
| `cofactor-X3-n4-Qqt-simplified.txt` | Braid normalization and collection yield 166 terms over four relations. |
| `cofactor-X3-n4-Qqt-three-rels.txt` | Conjugating the shifted two-strand relation gives 161 terms over three relations. |
| `cofactor-X3-n4-Qqt-ZT.txt` | Compatibility copy using the canonical `Z_2`, `T_2`, `T_3` labels. |

The last three files describe equality in the braid group algebra. They need
braid normalization and cannot be passed unchanged to the free-algebra CLI
checker. Their recorded common denominator is `t^4*(q-1)*(t+q^2)`.
The thesis uses a further simplified denominator, so its statement and these
retained files should be read with their respective denominator conditions.
The 695-term fixture is a reproducible computational witness; it is not the
latest mathematical proof.

The transformation code is
`tools/cofactor/compress_cofactor.py`, `simplify_cofactor.py`,
`reduce_to_three_rels.py`, and `rename_to_ZT.py`, supported by
`symbolic_identity_check.py` and `relations_naming.py` in the same directory.
`verify_ZT_and_ideal_equality.py` checks the relabelled result and the relation
comparison. `tools/auxiliary/independent_verify_three_rels.py` supplies a
separately written parser and summation for the compact result, sharing the
relation builder, polynomial arithmetic and Garside helpers.
These are historical scripts with fixed workspace paths. Their role here is
recorded provenance; they require path and dependency checks before rerunning.

## 4. Rerunning discovery with an external solver

The historical entry point is `reproduce_X3_zero.sh`, which invokes
`tools/cofactor/pipeline.py --legacy-x3`. The driver uses `paths.py`, `samples.py`
and `tools/cache.py`, with the maintained proof, reconstruction, composition
and verification modules. It can generate new inputs and sampled proofs,
reconstruct candidate coefficients, and invoke the specialized transformation
steps.

A fresh solver rerun is not currently a portable supported reproduction route.
The execution layer imports a `repo_settings` module from the sibling research
workspace, needs a configured `f4ncgb` executable, and retains fixed paths in the
transformation scripts. The public repository contains saved bases and
specialized cofactor samples, but not the full collection of per-point
nonexpanded proof logs expected by this driver. Several archived scripts also
retain paths from before their move into `archive/`. Installing the package
alone does not provide that execution environment.

The original successful transcript is
`logs_and_traces/reproduce-X3-zero-transcript.txt`. It records a cached anchor
computation and reuse of saved lifting data; it is not evidence of an entirely
fresh discovery run. The archived predecessor driver and its scripts are listed
below for tracing that computation. A future portable discovery procedure
should specify solver revision and arguments, configuration, sample inputs and
proof logs, and verify the final identity exactly.

### Retained f4ncgb parser patch

The [f4ncgb patch notes](../patches/f4ncgb/README.md) describe the retained
`patches/f4ncgb/msolve-parser-fix.patch` and its upstream baseline. The patch
corrects parser handling of CRLF input and unknown variable names. It is
unnecessary for the solver-free certificate checks above; correctly formed
LF input does not necessarily exercise the affected cases. The patch's
applicability is recorded separately from compilation and regression execution.

## Files contributing to reproduction

The [machine-readable inventory](x3-reproduction-files.json) expands all file
families below into exact repository-relative paths, with byte sizes, SHA-256
hashes and a role for each file. It covers source code, scientific inputs and
outputs, fixtures and test/build support. This guide and documentation being
edited alongside it are excluded from the hash list. The old fixture manifest
records the separate provenance of its three immutable files.

| Files or family | Role |
| --- | --- |
| `src/zipper_algebra/{__init__,__main__,cli,presentation,zipper_in_ms,core}.py` | Package entry points, presentation construction, exact arithmetic and basis parsing. |
| `src/zipper_algebra/cofactor/*.py` | Denominator representation, certificate input, proof expansion, recorded reductions, interpolation, composition and exact checking. |
| Matching `tools/` and `tools/cofactor/` modules | Compatibility imports used by the retained tests and execution scripts. |
| `tools/cofactor/{pipeline,paths,samples}.py`, `tools/cache.py` | Solver orchestration, sample and artifact naming, and cached solver outputs. |
| The transformation and compact-checker scripts named in section 3 | Intermediate certificate construction and braid-algebra verification. |
| `tests/fixtures/legacy_x3/{basis-symbolic-Qqt.json,skeleton-n4-Q.json,proof-n4-Q-noex.txt,manifest.json}` | Immutable symbolic proof, reduction skeleton and anchor proof, with provenance hashes. |
| `logs_and_traces/f4ncgb/basis-n4-Q.txt`, `cofactor-X3-n4-Q.txt`, `skeleton-n4-Q.json` | Original specialized anchor basis, numeric cofactor and reduction record. |
| `logs_and_traces/f4ncgb/reproduce-{zipper-n4-Q-q5-t7.ms,basis-n4-Q.txt,proof-n4-Q-noex.txt}` | Saved inputs and outputs of the recorded reproduction run. |
| `logs_and_traces/f4ncgb/basis-n4-Q-q*-t*.txt` and `cofactor-X3-n4-Q-q*-t*.txt` | Retained specialized bases and numeric cofactor expressions at parameter samples. |
| `logs_and_traces/f4ncgb/basis-symbolic-Qqt.{json,txt}` | Other retained symbolic lifting records. The JSON and text rendering come from different revisions and are not interchangeable with the immutable fixture. |
| `logs_and_traces/f4ncgb/lift_samples/coefficients.json` and `lift_samples_canonical/coefficients.json` | Two saved collections, each containing 212 numeric samples with 695 coefficient slots. |
| `ms_inputs/qt_sweep/*.ms` | Specialized solver inputs and reference inputs, including `zipper-n4-Q.ms` and `zipper-n4-X3-Q.ms`. |
| The five `cofactor-X3-n4-Qqt-*.txt` files in section 3 | Saved symbolic identity and each stage leading to the compact result. |
| `logs_and_traces/reproduce-X3-zero-transcript.txt` | Recorded execution, term counts and checks. |
| `archive/cofactor-pre-consolidation/*.py` and `reproduce_X3_zero.sh` in that directory | Historical extraction, sample collection, symbolic lifting, composition and predecessor driver. |
| `archive/gap/{zipper,traced_X3_witness_n3,reduce_X3_against_B,trace_replay,test_n4_presentation}.gap` | Independent GBNP presentation, traced witness, saved-basis reduction and historical replay/comparison routes. |
| `archive/gap/traces/grobner_zipper_bigelow-untraced-{3,4}-*.gap` | Saved finite-field bases for the independent GAP route. |
| `logs_and_traces/f4ncgb/X3-membership-singular-independent.txt` and `singular_X3_{Q_d12,Q_d14,GFp_d10,GFp_d14}.log` | Singular/Letterplace membership evidence at the specialized parameter pair `(5,7)`. |
| `tests/*.py` and `tools/cofactor/tests_*.py` | Recomposition and historical comparisons, together with input-validation, arithmetic, context, reconstruction, installation and CLI regressions. |
| `patches/f4ncgb/msolve-parser-fix.patch` and `patches/f4ncgb/LICENSE` | Retained upstream parser correction and its license, with baseline and application instructions in the adjacent README. |
| `pyproject.toml`, `MANIFEST.in`, `.github/workflows/ci.yml`, `ci/check_{distribution,installed,sdist}.py` | Installation, inclusion of reproduction fixtures, distribution inspection and CI execution. |

The GAP scripts need GAP and GBNP. Saved-basis reduction and trace replay refer
to old `logs_and_traces/gap_traces/` paths; the retained untraced bases now live
under `archive/gap/traces/`, and the replay's traced anchor is not among those
saved bases. The Singular files are specialized evidence, not a generic
rational-function proof; their successful-run scripts are not stored alongside
the logs. Failed Singular attempts and unrelated variant experiments are not
required by the maintained reproduction steps and are excluded from this
inventory.
