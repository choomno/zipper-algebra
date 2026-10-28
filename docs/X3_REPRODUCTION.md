# Producing and auditing the X_3 cofactor certificate

The production route begins with the defining presentation and ends with a new
verified cofactor expression. The stages below identify the source, initial
data, intermediate files and verification at each step. Separate commands
then audit the saved computation without rerunning the solver or interpolation.

The [mathematical reference](mathematical-specification/mathematical-specification.pdf)
defines the sliders, presentation, coefficient field and specialization
conditions. The original solver computations used
[f4ncgb](https://gitlab.sai.jku.at/f4ncgb/f4ncgb), by Maximilian Heisinger and
Clemens Hofstadler; its requested
[CASC 2025 citation](https://doi.org/10.1007/978-3-032-09645-6_5) is included in
the reference and repository README. The
[solver setup](../patches/f4ncgb/README.md) records its source, parser patch and
build instructions.

## The computation and its conclusion

The initial data are the four-strand canonical presentation, the slider target
`X_3`, an anchor parameter pair, and a solver executable. The final data are a
finite cofactor expression over the defining relations, together with an exact
verification report. The retained witness has 695 collected terms. A newly
computed valid witness can have different coefficients, contexts or term count;
the mathematical success condition is its exact identity, not agreement with
those historical choices.

The presentation and coefficient conventions are defined in the compiled
mathematical reference. Its denominator and specialization conditions apply to
the final identity. The following route uses characteristic zero throughout
coefficient reconstruction. The relation-emission input does **not** append the
target as another defining relation.

## Source and execution dependencies

- Package source: `src/zipper_algebra/presentation.py`,
  `zipper_in_ms.py`, `core.py`, and the `cofactor/` modules named below.
- Experimental orchestration: `tools/cofactor/{pipeline,samples,paths}.py`,
  `tools/cache.py`, plus compatibility imports under `tools/`.
- External solver: f4ncgb, whose original upstream is JKU's GitLab. The recorded
  public setup pins upstream `f152bb0d040d7947d4296f9d68ca849ee135629c`
  (v1.1.3), applies `patches/f4ncgb/msolve-parser-fix.patch`, and builds the
  executable with a C++20 compiler, CMake, Boost, FLINT and GMP. See the solver
  setup for exact build commands and patch version. The parser regression build
  is separate from an execution of this cofactor-production route.
- Historical solver identity: the recovered settings describe the research
  executable as v1.0.1. In the retained solver Git checkout that tag resolves to
  `236e7a362ae5584525952ba735f9147037203ff7`. This identifies the tag, not the
  exact original executable bytes. Valid-input comparisons between versions do
  not by themselves establish identical basis numbering or proof traces.
- Workspace-only dependencies in the current driver:
  the surrounding research workspace’s `repo_settings.py` and `settings.toml`.
  These are outside this software repository. The settings loader is
  required even if lower-level functions receive explicit binary and timeout
  arguments. Its configured recovered solver executable is absent on the current
  machine.
- Optional factor inference uses Sympy; an editable installation with
  `python -m pip install -e '.[inference]'` supplies that extra. Verification
  itself uses the Python standard library.

The current driver exposes presentation, target, anchor, characteristic,
`--start`, `--step`, `--max`, `--jobs`, `--atoms`, `--refresh`, and
`--legacy-x3` arguments. It has **no** `--root`, `--output`, `--binary`,
`--f4ncgb`, solver-flags or solver-timeout command-line argument.

When the workspace module is available, its actual environment overrides are
`BRAID_BINARIES_F4NCGB`, `BRAID_PATHS_OUTPUT`, `BRAID_PATHS_CACHE`, and
`BRAID_LIMITS_TIMEOUT_SECONDS`. The settings file itself is optional because
that module contains defaults. These variables do not remove the module-import
dependency.

## 1. Construct and emit the defining relations

Maintained functions: `Presentation(4)` and `Presentation.emit_ms`, which call
`build_relations_ring` through the exact numeric coefficient backend. The
canonical presentation contains 31 defining relations in six free generators;
its currently emitted content hash is `091bbb92271ab870`.

A standalone emission command, runnable after package installation, is:

```sh
mkdir -p reproduction/x3
python -m zipper_algebra.zipper_in_ms \
  --n 4 --q 5 --t 7 --char 0 --untwist bigelow \
  --out reproduction/x3/input-5-7-c0.ms
```

Input: the presentation formulae and the integer anchor values. Output: an
msolve-format input whose first line specifies generator order, second line
specifies characteristic zero, and remaining lines contain only the 31
relations. The file is generated from the presentation; the symbolic-lift JSON
and cofactor certificate are not inputs.

Historical output: `logs_and_traces/f4ncgb/zipper-n4-Q-q5-t7.ms`.
A saved rerun copy also exists as `reproduce-zipper-n4-Q-q5-t7.ms`.
Current driver: `pipeline.run`, phase 1, calls `pres.emit_ms`.

## 2. Compute the anchor basis and nonexpanded proof

The recorded anchor command is equivalent to:

```sh
"$F4NCGB" \
  -i reproduction/x3/input-5-7-c0.ms \
  -o reproduction/x3/basis-5-7.txt \
  -p reproduction/x3/proof-5-7-noex.txt \
  -e 0 -v 2 -r 0
```

Here `F4NCGB` must name the built executable. The `-p` output records the
nonexpanded proof for each computed basis element. The production route expands
these references in Python; it does not request an expanded target cofactor
from the solver.

Current functions: `samples.SampleSet.run_anchor` -> `run_f4ncgb` -> subprocess.
The wrapper validates input grammar, records the basis and proof, and uses
`tools/cache.py` to cache both outputs by input bytes, basename, solver binary
identity and flags. `--refresh` forces the anchor computation; it does not
force regeneration of every sampled proof already on disk.

Recorded anchor artifacts: `basis-n4-Q.txt` and the originally named
`proof-n4-Q-noex.txt`, each with 137 entries in the historical transcript.
The original proof filename is now absent; its retained contents are available
in `tests/fixtures/legacy_x3/proof-n4-Q-noex.txt`, with a separate reproduction
copy at `logs_and_traces/f4ncgb/reproduce-proof-n4-Q-noex.txt`.

The wrapper's default flags are exactly `-e 0 -v 2 -r 0`; the settings loader's
`f4ncgb_maxiter` and `f4ncgb_maxdeg` values are not added to that command.
Upstream's default iteration limit is 10. A successful solver exit is not a
basis-completion certificate; a verified positive cofactor identity does not
require a completeness claim.

## 3. Parse proofs, expand references, and derive the slider reduction

Maintained functions:

- `cofactor.proof.Basis.load` and `Proof.load` parse the ordered basis and
  nonexpanded proof using the same generator order as the presentation.
- `Proof.expand_all` exposes numeric substitution of proof references, preserving
  left and right word contexts. The archived extractor explicitly uses this
  expansion and checks every basis identity. The current driver parses the
  proof for its dependency closure at this point; its symbolic expansion occurs
  later in `compose.expand_over_inputs`.
- `cofactor.skeleton.reduce` reduces the slider constructed by the presentation.
- `Skeleton.from_reduction` records every occurrence, multiplier, word context,
  leading-word key and transitive proof dependency; `Skeleton.save` writes it.

Input: the presentation and the actual anchor basis/proof. Output:
`skeleton-S3_1.json` in the current directory layout, or
`logs_and_traces/f4ncgb/skeleton-n4-Q.json` in the legacy layout.
The historical reduction has three steps on basis indices 16, 17 and 102,
zero remainder, and a proof closure with 25 basis elements. These indices are
observed anchor data, not required constants for a different solver run.

Historical producer: `archive/cofactor-pre-consolidation/extract_cofactor.py`
formerly lived under `tools/cofactor/`; its `main` also verifies that every
expanded basis proof reconstructs its stated basis polynomial and writes the
specialized numeric cofactor `cofactor-X3-n4-Q.txt`.
The 16 retained proof-comparison tests exercise parsing, expansion and reduction
against the immutable anchor proof. They are a replay of this stage, not a run
of the solver.

## 4. Generate parameter samples and collect compatible proof coefficients

Current functions: `SampleSet.next_points`, `emit`, `run_point`, `grow`,
`shape_census`, `modal_shapes`, and `stable` in `tools/cofactor/samples.py`.
For each selected pair the presentation is re-emitted and f4ncgb produces its
own basis and nonexpanded proof. Points are selected by the spiral procedure,
filtered by denominator conditions and the retained sampling exclusions.
Proof elements are matched by leading-word keys; the ordered contextual proof
shapes must agree on the closure needed by the target. A nonmodal anchor may
be replaced and its slider reduction recomputed.

Current per-point outputs:

- `input-Q-T-c0.ms`
- `basis-Q-T.txt`
- `proof-Q-T-noex.txt`

Legacy equivalents:

- `ms_inputs/qt_sweep/zipper-n4-Q-qQ-tT.ms`
- `logs_and_traces/f4ncgb/basis-n4-Q-qQ-tT.txt`
- `logs_and_traces/f4ncgb/proof-n4-Q-qQ-tT-noex.txt`

Historical producer: `lift_until_verified.py` calls
`extract_cofactor_qt.py --q Q --t T`. The latter emits the input, runs the
solver, parses and expands its proof, checks expansion, reduces the target,
verifies the specialized cofactor and writes
`cofactor-X3-n4-Q-qQ-tT.txt`. Its stored subprocess command uses `-e off`
without the current driver's explicit verbosity and reduction flags; that is
historical source behaviour, not the recommended portable command above.

The public tree retains many specialized bases and numeric cofactors, but
**none of the per-point nonexpanded proof files are present**. These proof
files are the crucial inputs for replaying coefficient collection. Their
absence does not prevent a fresh discovery run from regenerating them once the
execution layer and solver are available.

## 5. Reconstruct the relevant basis-proof coefficients

Current function: `cofactor.lift.lift(skeleton, samples, ...)`.
It collects rational values from every compatible slot in the needed proof
closure, reserves held-out points, infers or loads denominator factors, and
reconstructs numerator coefficients with exact arithmetic. Rank and held-out
checks determine whether a candidate can proceed. A failed fit or exhausted
sampling limit is inconclusive.

Input: the anchor skeleton, anchor and sampled bases/proofs, and optionally a
factor file from `--atoms`. Output: `lift-S3_1.json`, including its factors,
anchor, held-out count, reconstructed entries and factor-inference report.
The current driver embeds inferred factors and the report in this lift JSON;
it does not itself write separate `atoms.json` or `atoms-report.json` outputs.
The legacy factor prior/input is retained as
`logs_and_traces/f4ncgb/atoms-n4-Q.json`.

Historical producer: `symbolic_lift_x3.py` reads
`proof-n4-Q-noex.txt`, `skeleton-n4-Q.json`, and the per-point proof files.
Its `transitive_deps`/`load_stable_samples` functions select the closure and
compatible samples; `fit_rational` reconstructs each coefficient. It writes
`basis-symbolic-Qqt.json` and `.txt`, with old JSON keys `needed`/`lifted`.

The current saved `logs_and_traces/f4ncgb/basis-symbolic-Qqt.json` has the newer
`Lift` schema and records 99 stable points and 9 held-out points for 25 basis
elements. The immutable fixture
`tests/fixtures/legacy_x3/basis-symbolic-Qqt.json` retains the older
`needed`/`lifted` schema used to regenerate the 695-term fixture. The historical
shell driver's cached-lift check expects the old schema and cannot consume the
current JSON unchanged.

`run_lift_grid.py` and `run_lift_grid_canonical.py` instead collect whole
numeric cofactors into `lift_samples/coefficients.json` and
`lift_samples_canonical/coefficients.json`; each saved collection has 212
samples and 695 coefficient slots. These are related historical experiments,
**not the files read by the final basis-proof coefficient lifter**.

## 6. Compose the rational cofactor expression

Current functions: `pipeline._multiplier_samples`,
`cofactor.compose.lift_multipliers`, `expand_over_inputs`, and `compose`.
The driver re-reduces the target at compatible samples and reconstructs one
multiplier per contextual reduction occurrence. Composition recursively expands
the symbolic proof closure over the input relations and substitutes the result
into the slider reduction. Identical full contexts may then be collected.

Input: the reconstructed `Lift`, the slider `Skeleton`, the chosen presentation,
and compatible multiplier samples. Output:
`cofactor-S3_1-symbolic.txt` in the current layout, or
`cofactor-X3-n4-Qqt-symbolic.txt` in the legacy layout.

Historical producer: `compose_x3_cofactor_qt.py` reads the old
`basis-symbolic-Qqt.json`, expands the 25 required elements, and uses the
historical three symbolic multipliers to produce the 695-term expression.
The maintained recomposition test recreates that certificate body byte for byte
from the immutable old lift and skeleton, with no solver or interpolation.
That useful replay begins at this stage; it is not full production from the
presentation.

## 7. Verify and finish, or enlarge the sample set

Current producer: `pipeline.run`, phase 5, invokes `cofactor.verify.check` on
the newly composed certificate and actual presentation. A passing exact identity
finishes the production loop. Otherwise the driver grows the sample target and
retries within `--max`; lack of a verified certificate is reported as
inconclusive, with exit code 3. Historical producer:
`symbolic_identity_check.py` follows composition; the old
`lift_until_verified.py` repeats sampling, lifting and composition until this
exact check passes. Its old interpretations of failed fits have since been
corrected in the maintained implementation.

Run the maintained CLI again on the produced file with `--n 4 --target X_3
--pit --json` to save both the symbolic report and the independent degree-bound
grid report. This final report, rather than a sample fit, solver exit or term
count, establishes the cofactor identity. Retain the input, solver revision and
flags, all bases and proof logs, skeleton, lift, certificate and reports together.

## Current driver invocation in a configured research workspace

This is the complete existing driver command, conditional on the external
workspace configuration module being importable and the solver being built:

```sh
BRAID_BINARIES_F4NCGB=/absolute/path/to/f4ncgb \
BRAID_PATHS_OUTPUT=/absolute/path/to/x3-production \
BRAID_PATHS_CACHE=/absolute/path/to/x3-cache \
BRAID_LIMITS_TIMEOUT_SECONDS=1800 \
python -m tools.cofactor.pipeline \
  --n 4 --untwist bigelow --target 3,1 --anchor 5,7 --char 0 \
  --start 100 --step 40 --max 260 --jobs 1
```

The generic layout places files under
`/absolute/path/to/x3-production/cofactor/n4-bigelow/`. The final symbolic
certificate is `cofactor-S3_1-symbolic.txt`; the skeleton is
`skeleton-S3_1.json` and reconstructed proof is `lift-S3_1.json`.
`--legacy-x3` instead uses the historical tracked artifact names and invokes
the compression stages in place. Prefer a separate output directory for fresh
production. The conditional command therefore describes the existing research execution
environment; it is not a supported procedure for a standalone public checkout.

## Audit the retained symbolic identity

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

## Recompose the retained certificate from its symbolic proof

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

## Follow the compact certificate’s provenance

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
| The transformation and compact-checker scripts named in the compact-certificate section | Intermediate certificate construction and braid-algebra verification. |
| `tests/fixtures/legacy_x3/{basis-symbolic-Qqt.json,skeleton-n4-Q.json,proof-n4-Q-noex.txt,manifest.json}` | Immutable symbolic proof, reduction skeleton and anchor proof, with provenance hashes. |
| `logs_and_traces/f4ncgb/basis-n4-Q.txt`, `cofactor-X3-n4-Q.txt`, `skeleton-n4-Q.json` | Original specialized anchor basis, numeric cofactor and reduction record. |
| `logs_and_traces/f4ncgb/reproduce-{zipper-n4-Q-q5-t7.ms,basis-n4-Q.txt,proof-n4-Q-noex.txt}` | Saved inputs and outputs of the recorded reproduction run. |
| `logs_and_traces/f4ncgb/basis-n4-Q-q*-t*.txt` and `cofactor-X3-n4-Q-q*-t*.txt` | Retained specialized bases and numeric cofactor expressions at parameter samples. |
| `logs_and_traces/f4ncgb/atoms-n4-Q.json` | Historical denominator-family prior for reconstruction. |
| `logs_and_traces/f4ncgb/zipper-n4-Q-q5-t7.ms` | Original emitted anchor presentation, distinct from the separately retained reproduction copy. |
| `logs_and_traces/f4ncgb/basis-symbolic-Qqt.{json,txt}` | Other retained symbolic lifting records. The JSON and text rendering come from different revisions and are not interchangeable with the immutable fixture. |
| `logs_and_traces/f4ncgb/lift_samples/coefficients.json` and `lift_samples_canonical/coefficients.json` | Two saved collections, each containing 212 numeric samples with 695 coefficient slots. |
| `ms_inputs/qt_sweep/*.ms` | Specialized solver inputs and reference inputs, including `zipper-n4-Q.ms` and `zipper-n4-X3-Q.ms`. |
| The five `cofactor-X3-n4-Qqt-*.txt` files in the compact-certificate section | Saved symbolic identity and each stage leading to the compact result. |
| `logs_and_traces/reproduce-X3-zero-transcript.txt` | Recorded execution, term counts and checks. |
| `archive/cofactor-pre-consolidation/*.py` and `reproduce_X3_zero.sh` in that directory | Historical extraction, sample collection, symbolic lifting, composition and predecessor driver. |
| `archive/gap/{zipper,traced_X3_witness_n3,reduce_X3_against_B,trace_replay,test_n4_presentation}.gap` | Independent GBNP presentation, traced witness, saved-basis reduction and historical replay/comparison routes. |
| `archive/gap/traces/grobner_zipper_bigelow-untraced-{3,4}-*.gap` | Saved finite-field bases for the independent GAP route. |
| `logs_and_traces/f4ncgb/X3-membership-singular-independent.txt` and `singular_X3_{Q_d12,Q_d14,GFp_d10,GFp_d14}.log` | Singular/Letterplace membership evidence at the specialized parameter pair `(5,7)`. |
| `tests/*.py` and `tools/cofactor/tests_*.py` | Recomposition and historical comparisons, together with input-validation, arithmetic, context, reconstruction, installation and CLI regressions. |
| `patches/f4ncgb/msolve-parser-fix.patch` and `patches/f4ncgb/LICENSE` | Retained upstream parser correction and its license, with baseline and application instructions in the adjacent README. |
| `pyproject.toml`, `MANIFEST.in`, `.coveragerc`, `.github/workflows/ci.yml`, `ci/check_{coverage,distribution,installed,sdist,x3_manifest}.py` | Installation, fixture inclusion, complete package coverage, inventory and distribution checks, and CI execution. |

The GAP scripts need GAP and GBNP. Saved-basis reduction and trace replay refer
to old `logs_and_traces/gap_traces/` paths; the retained untraced bases now live
under `archive/gap/traces/`, and the replay's traced anchor is not among those
saved bases. The Singular files are specialized evidence, not a generic
rational-function proof; their successful-run scripts are not stored alongside
the logs. Failed Singular attempts and unrelated variant experiments are not
required by the maintained reproduction steps and are excluded from this
inventory.
