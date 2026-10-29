# Producing and checking the X_3 cofactor identity

There are three useful routes:

| Route | Starting data | What it establishes |
| --- | --- | --- |
| Fresh production | Presentation, target, solver executable and computation limits | A newly computed cofactor identity, if exact verification succeeds. |
| Recomposition | Saved symbolic lift and contextual reduction skeleton | The retained certificate body can be reconstructed from those coefficients. |
| Certificate check | Installed package and one certificate | The stated free-algebra identity holds exactly. |

The [mathematical reference](mathematical-specification/mathematical-specification.pdf)
defines the algebras, coefficients, localization and verification conclusions.
The [source inventory](x3-reproduction-files.json) records the curated current
files. Bulk historical experiments and samples are available through the pinned
archive described below.

## Fresh production

Install the inference extra and build
[f4ncgb](https://gitlab.sai.jku.at/f4ncgb/f4ncgb), by Maximilian Heisinger and
Clemens Hofstadler. Their requested
[CASC citation](https://doi.org/10.1007/978-3-032-09645-6_5), original source,
pinned revision and parser patch are documented in the
[solver setup](../patches/f4ncgb/README.md).

```sh
python -m pip install '.[inference]'
zipper-algebra discover --n 4 --target X_3 \
  --solver /absolute/path/to/f4ncgb \
  --solver-revision 789504d7c7da6af6946faecb73421b037766a5b8 \
  --output runs/x3 --cache runs/solver-cache \
  --anchor 5,7 --start 100 --step 40 --max 260 --jobs 1 --json
```

After installation, the command works outside the source checkout. It uses
explicit configuration and generates all samples and proofs itself. The only
initial mathematical data are the selected defining presentation and target.
Saved lifts, cofactors, a populated cache and the old research workspace are
not inputs.

A local fresh run at source revision
`d30b36aa2da0819dfa28954bd9223ef72163edf5` used 100 fresh parameter points,
99 compatible samples and 107 reconstructed coefficient slots. It produced a
129-term certificate in 172.205 seconds including validation; reported F4 time
was about 5.36 seconds and coefficient reconstruction about 2.521 seconds. Exact
symbolic verification and a 25 by 7 product grid both passed. The curated
[example](../examples/x3/certificate.txt) is this free-algebra certificate.
These are local results; CI exercises a smaller native solver integration.
Timing and coefficient choices can vary across runs.

A small existing consequence
exercises the same route without adding its target as a defining relation:

```sh
zipper-algebra discover --n 3 --target '(Theta_3-q^2)X_3' \
  --solver /absolute/path/to/f4ncgb \
  --solver-revision 789504d7c7da6af6946faecb73421b037766a5b8 \
  --output runs/three-strand --start 12 --step 12 --max 96 --json
```

Here the full-twist relation follows from the canonical third-level
untwisting relation and inverse relations. This example is a software
integration case, not a new structural theorem.

### Read the result

A verified result has exit code `0` and `"status": "verified"` in
`result.json`. The report names the certificate and the immutable product
files, with their hashes. They are published together under
`verified/<generation>/`:

- `certificate.txt`: readable cofactor expression.
- `certificate.json`: versioned certificate bound to its presentation and target.
- `skeleton.json`: contextual reduction occurrences and proof closure.
- `lift.json`: reconstructed basis-proof coefficients.
- `multipliers.json`: reconstructed multipliers for those ordered occurrences.
- `provenance.json`: source revision, local changes and hashes of the five products.

`samples.sqlite3` stores every sample's input, basis, nonexpanded proof, process
output and validation record. Sample keys and generation identifiers in
`result.json` locate these records. The cache uses the same database format and
reuses matching, validated samples. The default `--solver-io files` uses
temporary files only while a solver call is active. The
[optional SQLite-enabled solver](../patches/f4ncgb/README.md#optional-direct-sqlite-io)
supports `--solver-io sqlite`: input, basis and proof bytes pass directly through
the database. Set `--solver-revision` to that executable's actual source revision.
Both routes use the same proof validation and sample storage. SQLite may create
transient journal files during transactions.

To inspect one sample, take its key and optional generation from the result:

```sh
zipper-algebra export-sample --store runs/x3/samples.sqlite3 \
  --key SAMPLE_KEY --generation GENERATION --output runs/sample-inspection --json
```

This restores the five original solver files and their record, retaining their
original producer provenance. The [provenance guide](PROVENANCE.md) explains the
revision, local-change and source-hash fields.

With a single successful product generation, discovery keeps ten files:
one database, six product files, the candidate skeleton, the result and the
directory's ignore file. The sample count changes the database's contents,
not the number of persistent sample files. A separate shared cache adds its own
database and ignore file. Refreshed product generations remain available.
The earlier directory layout used 607 pipeline files for a 100-sample run:
600 sample files and seven products or reports. The previously reported total
of 609 also included two reports added after that run; those were not pipeline
outputs. The new database layout avoids that per-sample file growth.

A result with exit code `3` is inconclusive. A nonzero remainder from an
uncertified basis, an exhausted fit budget or unavailable compatible proof
shapes does not establish nonmembership. Invalid input has exit code `2`.
Process completion and proof validation do not assert that a basis is complete.

### The mathematical stages and their implementation

All paths below are relative to `src/zipper_algebra/`.

| Stage | Initial data and output | Principal source |
| --- | --- | --- |
| 1. Construct the presentation | Parameters become an ordered alphabet, relations, labels and symbolic identity. The target is not automatically adjoined. | `presentation.py`, `relations.py` |
| 2. Emit and solve at an anchor | The defining polynomials become msolve input. f4ncgb emits basis polynomials and nonexpanded proofs. Every basis identity is reconstructed against the actual inputs before publication. | `execution/runner.py`, `execution/config.py`, `execution/artifacts.py`, `execution/native_io.py`, `execution/store.py` |
| 3. Record the target reduction | Ordered left/right contexts, basis keys and multipliers become a reduction skeleton with a dependency closure. | `cofactor/proof.py`, `cofactor/skeleton.py` |
| 4. Obtain compatible samples | Fresh parameter points produce validated proofs; leading-word keys and contextual proof shapes align their coefficients. A nonmodal anchor may be replaced. | `execution/sampling.py` |
| 5. Reconstruct and compose | Exact fitting proposes rational coefficients. Held-out values, rank and denominator checks constrain the candidate. Composition retains every contextual occurrence. | `cofactor/lift.py`, `univariate.py`, `inference.py`, `compose.py` under `cofactor/` |
| 6. Check and publish | The complete free-algebra identity is verified exactly before an immutable product generation is selected. Failed attempts grow sampling within the requested cap. | `cofactor/verify.py`, `cofactor/artifacts.py`, `execution/pipeline.py` |

`coefficients.py` and `polynomial.py` supply production arithmetic.
`cofactor/ncpoly.py` supplies the separate symbolic accumulator used by the
checker. `cli.py` connects the installed commands. This source map identifies
responsibilities; the inventory also lists the underlying arithmetic and
test/build support.

To inspect just the input, without invoking a solver:

```sh
zipper-algebra emit --n 4 --q 5 --t 7 --char 0 \
  --out input-5-7.ms --json
```

The canonical four-strand input has 31 relations in six free generators and
presentation hash `091bbb92271ab870`. The historical anchor input is available
in the `inputs` archive group. Discovery never automatically appends its target.
Explicit zeroing and extra relations are part of the selected presentation and
can make membership immediate; the resulting certificate concerns that presentation.

## Check the curated example

From the source checkout after installation:

```sh
zipper-algebra verify \
  examples/x3/certificate.txt \
  --n 4 --target X_3 --pit --json
```

Success is exit code `0`, overall `"pass": true`, and passing `exact` and
`pit` reports. The curated example has 129 collected terms. New valid
certificates can differ in coefficients, contexts and term count.

The first checker compares the cleared identity symbolically. The second uses
a product grid large enough for its derived degree bounds, so it is also exact.
They share input validation and presentation construction while accumulating
the identity with separate arithmetic. Neither route needs the old solver,
sample collection or thesis source. A false identity has code `1`; invalid data
has code `2`.

## Recompose the historical regression witness

The older 695-term witness is retained as
`tests/fixtures/legacy_x3/cofactor-X3-n4-Qqt-symbolic.txt`, alongside its frozen
symbolic lift, skeleton, anchor proof and numeric bases. This comparison starts
from those lift and skeleton fixtures and recreates the historical certificate
body byte for byte:

```sh
python -m unittest \
  tests.reproduction.test_x3.TestCompose.test_composed_body_is_byte_identical_to_the_shipped_certificate \
  -v
python -m unittest tests.reproduction.test_legacy_proofs -v
```

Expect one passing recomposition test and 16 passing historical proof comparisons.
The latter compare the retained 137-entry anchor with a frozen historical
implementation under `tests/reference/`. That implementation shares coefficient
machinery with the package; it is an implementation comparison rather than wholly
independent arithmetic. Outputs are temporary.

The installed `zipper-algebra compose` command accepts a lift, skeleton,
symbolic multiplier file and explicit presentation/target. It checks the result
exactly before replacing its output. New bound artifacts include the coefficient
domain, ordered generators and relations, full presentation digest, target
polynomial and applicable order. Historical readers remain explicit.

## Historical calculations

The `x3-history` archive group includes the historical compression sequence:

| File under the restored `logs_and_traces/f4ncgb/` | Role |
| --- | --- |
| `cofactor-X3-n4-Qqt-symbolic.txt` | 695-term free-algebra witness, also retained as a regression fixture. |
| `cofactor-X3-n4-Qqt-compressed.txt` | 1077-term substitution of dependent relations. |
| `cofactor-X3-n4-Qqt-simplified.txt` | 166-term braid-normalized expression. |
| `cofactor-X3-n4-Qqt-three-rels.txt` | 161-term expression over three relations. |

The last two use equality in a braid group algebra. Their common denominator
is `t^4*(q-1)*(t+q^2)`. They require the corresponding checker and cannot be
passed unchanged to a literal free-word checker.
[Compact-certificate commands](../experiments/x3/README.md) reproduce and check
those transformations in separate output files. They retain a separately
written parser/summation with its shared arithmetic and normalization dependencies
identified.

The [Singular producer](../singular_run/README.md) preserves the 31-relation
specialized reduction at `(q,t)=(5,7)`. GAP/GBNP source programs remain under
`archive/gap/`; their saved bases are in the `gap` archive group. Their specialized
conclusions, degree bounds and old path/completion assumptions remain distinct
from a generic rational-function certificate.

## Historical provenance and file inventory

The [inventory](x3-reproduction-files.json) records exact paths, sizes, SHA-256
hashes and roles for retained source, data and test/build support. Regenerate or
check it with:

```sh
python ci/build_x3_manifest.py --check
python ci/check_x3_manifest.py
```

The current checkout keeps one certificate example and the immutable data used
by the regression suite. Bulk inputs, logs, saved bases and variant survey results
are recorded in [data/research-archive.json](../data/research-archive.json), pinned
to Git revision `d30b36aa2da0819dfa28954bd9223ef72163edf5` with original paths and
SHA-256 hashes. Restore all groups or a selected group into a separate directory:

```sh
python -m experiments.restore_data --output .research-data
python -m experiments.restore_data --group x3-history --output .research-data
python -m experiments.restore_data --group variants --output .research-data
python -m experiments.restore_data --group gap --output .research-data
```

The restore command verifies the recorded hashes. Paths below that output
directory follow the historical repository layout. The default tests, checking
the curated example and fresh production need no restored archive.
The archive manifest’s `previous_aliases` field records the historical names of
byte-identical duplicates. The two old coefficient-collection names referred to
the same 212 samples and 695 coefficient slots.

The `reproduce-` input/basis/proof triplet is a different experiment:
it has **32 relations including the target itself**, and 131 basis/proof entries.
Its proof uses the appended target. It is not the canonical 31-relation,
137-entry computation and cannot establish membership in that original ideal.
The canonical 137-entry proof remains under `tests/fixtures/legacy_x3/`.

Old per-point nonexpanded proof logs were not recovered. Consequently the
original numerical fitting cannot be replayed solely from retained samples;
fresh production generates new proofs. The current symbolic JSON and historical
text rendering also come from different revisions and are not interchangeable
with the frozen old-schema lift.

Retired extraction/grid/lifting drivers remain accessible at the
[phase 1 source baseline](https://github.com/Choomno/zipper-algebra/tree/0db1dd9e2a33fb0e5d07e216c8a5354911446e97/archive/cofactor-pre-consolidation).
The inventory records these pinned historical sources separately. Current
production uses the installed implementation described above.
