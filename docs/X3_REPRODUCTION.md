# Producing and checking the X_3 cofactor identity

The main routes are:

| Route | Starting data | What it establishes |
| --- | --- | --- |
| Fresh production | Presentation, target, solver executable and computation limits | A newly computed cofactor identity, if exact verification succeeds. |
| Python certificate checks | Installed package and one certificate | The stated free-algebra identity holds, by symbolic accumulation and an exact product grid. |
| Independent GAP check | Exported standalone script | GAP verifies the generic identity over `Q(q,t)` using its own arithmetic and noncommutative expansion. |

The [mathematical reference](mathematical-specification/mathematical-specification.pdf)
defines the algebras, coefficients, localization and verification conclusions.
The [source inventory](x3-reproduction-files.json) records the curated current
files needed to inspect the current implementation and its tests.

## Fresh production

Install the symbolic extra and the
[FiniteFlow dependency](RECONSTRUCTION.md#installation), then build
[f4ncgb](https://gitlab.sai.jku.at/f4ncgb/f4ncgb), by Maximilian Heisinger and
Clemens Hofstadler. Their requested
[CASC citation](https://doi.org/10.1007/978-3-032-09645-6_5), original source,
pinned revision and parser patch are documented in the
[solver setup](../patches/f4ncgb/README.md).

```sh
python -m pip install '.[symbolic]'
zipper-algebra discover --n 4 --target X_3 \
  --solver /absolute/path/to/f4ncgb \
  --solver-revision 789504d7c7da6af6946faecb73421b037766a5b8 \
  --output runs/x3 --cache runs/solver-cache \
  --anchor 41,43 --max-primes 5 --reconstruction-degree 24 --threads 1 --json
```

After installation, the command works outside the source checkout. It uses
explicit configuration and generates its anchor proof and contextual support itself. The only
initial mathematical data are the selected defining presentation and target.
Saved lifts, cofactors, a populated cache and the old research workspace are
not inputs.

The curated [example](../examples/x3/certificate.txt) is the 129-term
free-algebra certificate produced by the earlier sampled coefficient-fitting
route. Current discovery instead solves the generic contextual system directly
with FiniteFlow. Cofactor representations are not unique; a fresh run can have
different coefficients and term counts. A fresh local run with anchor `(41,43)` generated 149 contextual columns,
an augmented matrix of shape `178 × 150`, generic rank 146 and a four-column
kernel. It produced a 120-term cofactor, passing symbolic verification and the
exact 17-by-6 product grid. Full X₃ computations remain local research checks,
outside CI.

A small existing consequence
exercises the same route without adding its target as a defining relation:

```sh
zipper-algebra discover --n 3 --target '(Theta_3-q^2)X_3' \
  --solver /absolute/path/to/f4ncgb \
  --solver-revision 789504d7c7da6af6946faecb73421b037766a5b8 \
  --output runs/three-strand --threads 1 --json
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
- `support.json`: contextual input-relation triples and ordered free-word rows.
- `matrix.json`: the generic matrix `[A|-h]` over `Q(q,t)`.
- `kernel.json`: the complete independently checked kernel and rank witness.
- `reconstruction.json`: FiniteFlow controls and hashes of its loaded implementation.
- `solver-runs.json`: presentation/target binding and the exact sample keys,
  immutable generations and artifact hashes recorded by this discovery. Its
  database paths are relative to the discovery output directory.
- `provenance.json`: source revision, local changes and hashes of the other bundle files.

`samples.sqlite3` stores every sample's input, basis, nonexpanded proof, process
output and validation record. Sample keys and generation identifiers in
`result.json` and the bundle's `solver-runs.json` locate these records. The cache uses the same database format and
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

One successful product generation keeps thirteen persistent files: one database,
nine bundle files, the candidate skeleton, `result.json`, and the root ignore
file. Temporary SQLite journals are excluded from this count. The report lists
the nine bundle paths and their hashes; these include support, augmented matrix,
kernel, cofactor, source provenance and dependency identities. A separate shared cache
adds its own database and ignore file. Refreshed generations remain available,
including their links to the original solver generations.

A result with exit code `3` is inconclusive. A nonzero remainder from an
uncertified basis, insufficient contextual support or an exhausted reconstruction budget does not establish nonmembership. Invalid input has exit code `2`.
Process completion and proof validation do not assert that a basis is complete.

### The mathematical stages and their implementation

All paths below are relative to `src/zipper_algebra/`.

| Stage | Initial data and output | Principal source |
| --- | --- | --- |
| 1. Construct the presentation | Parameters become an ordered alphabet, relations, labels and symbolic identity. The target is not automatically adjoined. | `presentation.py`, `relations.py` |
| 2. Emit and solve at an anchor | The defining polynomials become msolve input. f4ncgb emits basis polynomials and nonexpanded proofs. Every basis identity is reconstructed against the actual inputs before publication. | `execution/runner.py`, `execution/config.py`, `execution/artifacts.py`, `execution/native_io.py`, `execution/store.py` |
| 3. Record the target reduction | Ordered left/right contexts, basis keys and multipliers become a reduction skeleton with a dependency closure. | `cofactor/proof.py`, `cofactor/skeleton.py` |
| 4. Generate word contexts | Expand the relevant proof DAG, propagating left/right contexts to the original input relations. Deduplicate triples without using sampled coefficients as symbolic answers. | `cofactor/discovery.py` |
| 5. Form and reconstruct the generic matrix | Expand each contextual defining relation, append the negative target, reconstruct a full kernel with FiniteFlow, and check its full generic identity and both rank bounds. | `cofactor/discovery.py`, `matrices.py`, `reconstruction.py`, `verification/linear.py` |
| 6. Check and publish | Normalize a kernel column with nonzero target coordinate into a cofactor, check the complete free-algebra identity, and atomically select an immutable product generation. | `cofactor/verify.py`, `cofactor/artifacts.py`, `execution/pipeline.py`, `execution/dependencies.py` |

`fields.py`, `verification/working_field.py`, `coefficients.py` and
`polynomial.py` supply production arithmetic.
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
presentation hash `091bbb92271ab870`. Discovery never automatically appends its target.
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

## Independent algebra-system checks

[GAP](https://www.gap-system.org/) can independently expand the supplied generic identity over `Q(q,t)` in a
free associative algebra. With GAP installed, run from the checkout:

```sh
mkdir -p runs/x3 &&
python -m zipper_algebra.verification.gap \
  examples/x3/certificate.txt --n 4 --target X_3 > runs/x3/check.g &&
gap -A -q -r -T --quitonbreak runs/x3/check.g
```

Expect `PASS: exact cofactor identity over Q(q,t)` and exit code `0`. A nonzero
residual prints `FAIL` and exits `1`; GAP errors also cause a nonzero exit.
This checks the full generic identity, without sampling or a Gröbner-basis
computation. The generated script is standalone: running it requires GAP but
no Python package or additional data files.

The exporter shares zipper-algebra's input parser and defining presentation.
GAP supplies separate rational-function arithmetic and noncommutative expansion.
The script's comment header records the producing source revision, local changes,
source hashes and input-certificate hash. Mathematical denominator hypotheses
remain those stated in the specification.

The [Singular calculation](../singular_run/README.md) independently tests X₃
membership at `(q,t)=(5,7)` over the rationals with degree bound 12. It uses the
canonical 31 defining relations. This is a specialized membership computation;
the generic certificate is checked by the symbolic and product-grid routes above.

## Standalone coefficient-fitting records

The older fitting and composition APIs remain available for explicitly supplied
data. They are separate from the FiniteFlow discovery path above.


The installed `zipper-algebra compose` command accepts a lift, skeleton,
symbolic multiplier file and explicit presentation/target. It checks the result
exactly before replacing its output. New bound artifacts include the coefficient
domain, ordered generators and relations, full presentation digest, target
polynomial and applicable order. Bound JSON records use version-2 envelopes;
older envelope versions are rejected.

## Source inventory

The [inventory](x3-reproduction-files.json) records exact paths, sizes, SHA-256
hashes and roles for retained source, data and test/build support. Regenerate or
check it with:

```sh
python ci/build_x3_manifest.py --check
python ci/check_x3_manifest.py
```

The maintained reproduction path uses the curated example and newly generated
solver inputs and outputs. Retired datasets and drivers remain in Git history.
Current production uses the installed implementation described above.
