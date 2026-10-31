# Certificate workflow

The [mathematical reference](docs/mathematical-specification/mathematical-specification.pdf)
states the objects, coefficient domains, hypotheses and conclusions. The installed
commands handle presentation descriptions, relation emission, discovery, composition
and exact checking.

## Describe or emit a presentation

```sh
zipper-algebra describe --n 4 --json
zipper-algebra emit --n 4 --q 5 --t 7 --char 0 --out input.ms --json
```

Both use one relation builder. The presentation hash identifies the ordered
symbolic defining data, independently of specialization. Labels alone do not
identify an algebra. Discovery never automatically adjoins its target. Explicit
`--extra` and `--zero-from` options are part of the selected defining presentation;
they can make membership tautological. A certificate then concerns that augmented
presentation.

## Produce a certificate

Install the symbolic extra, build f4ncgb using the
[solver setup](patches/f4ncgb/README.md), and install the
[FiniteFlow dependency](docs/RECONSTRUCTION.md#installation). Then:

```sh
zipper-algebra discover --n 4 --target X_3 \
  --solver /absolute/path/to/f4ncgb \
  --solver-revision 789504d7c7da6af6946faecb73421b037766a5b8 \
  --output runs/x3 --cache runs/solver-cache \
  --anchor 41,43 --max-primes 5 --reconstruction-degree 24 --threads 1 --json
```

The [reproduction guide](docs/X3_REPRODUCTION.md) explains the stages and
provides a smaller example. Initial inputs are the presentation, target and
explicit execution configuration. Samples and proofs are generated and managed
automatically. Sample inputs, solver outputs and validation records are stored
in one `samples.sqlite3` database. The default `--solver-io files` uses temporary
files for active solver calls. With the
[optional SQLite-enabled solver](patches/f4ncgb/README.md#optional-direct-sqlite-io),
`--solver-io sqlite` reads and writes directly in that database. Use the source
revision of the selected executable in `--solver-revision`.

Each solver run has a versioned identity containing the actual input, coefficient
domain, ordered generators, order, executable digest, flags and limits.
The target belongs to the subsequent reduction and certificate, so different
targets can reuse the same presentation calculation. The Python function
`zipper_algebra.execution.run_basis_point` runs the presentation without selecting
a target; `run_point` additionally validates a named target before using the same
solver cache. Each discovery result retains its own target identity and links to
the shared solver key and immutable generation. Context supports, augmented
matrices, kernels and cofactors are specific to each target.
Its outputs are parsed and every numeric proof identity is checked before the
run is published. Atomic records distinguish usable proof data from a claim
of basis completeness.

New solver requests use `zipper-solver-request-v3`. Earlier target-bound v2
records remain readable and exportable. They are separate cache entries and
are not relabelled as new requests.

Reconstruction expands one fresh proof graph into contextual occurrences of the
original relations. Their generic coefficients form an augmented linear system
whose kernel is produced by FiniteFlow. The complete exact kernel identity,
independence and opposite rank bound are checked before extracting a cofactor.
The separate free-algebra check of the complete cofactor identity is required
before publication. No fitted denominator family or repeated solver samples
are required.

A successful result atomically selects an immutable product bundle with its hashes.
`result.json` reports the presentation/target, anchor, sampled runs and certificate
path. The bundle's `solver-runs.json` preserves its sample keys and immutable
generations after subsequent runs replace `result.json`. An interrupted
publication preserves the preceding valid result.

The result and solver records retain the producer's commit, recorded local
changes and exact source-file hashes. The final bundle includes a shared
`provenance.json`; standalone saves receive adjacent sidecars. See
[output provenance](docs/PROVENANCE.md) for field meanings and build origins.

## Verify or recompose

```sh
zipper-algebra verify \
  examples/x3/certificate.txt \
  --n 4 --target X_3 --pit --json
```

The checker validates the presentation, target, words, labels, coefficient domain
and denominator data before arithmetic. It compares the cleared identity exactly.
The optional product-grid check derives sufficient degree bounds. A rational
identity specializes wherever its displayed denominators remain nonzero.

New discovery bundles retain their augmented matrix and complete kernel witness,
so they can also be checked independently with `ExactMatrix.from_json`,
`KernelCertificate.from_json` and `verify_kernel`; see the
[reconstruction interface](docs/RECONSTRUCTION.md). Run `verify` on the bundle's
`certificate.txt` to check the cofactor independently of FiniteFlow.

The separate `compose` command remains available for explicitly supplied
skeleton/lift/multiplier records from the coefficient-fitting API. Those records
are no longer products of `discover`. Composition checks the full identity
before replacing its output.

| Exit code | Meaning |
| --- | --- |
| `0` | Requested exact checks passed, or discovery produced a verified identity. |
| `1` | A well-formed certificate or composed candidate is false. |
| `2` | Input, presentation, target, configuration or file is invalid. |
| `3` | Discovery is inconclusive within its support, reconstruction or execution limits. |

A nonzero remainder from an uncertified basis establishes no nonmembership
conclusion. Process success does not prove a completed basis. Sample agreement
does not prove a rational-function identity.

## Development and verification

[TESTING](docs/TESTING.md) describes the unit and integration tests,
coverage gate, isolated distributions and native solver integration.
[INDEX](INDEX.md) maps the source and experimental recipes.
The curated example is a 129-term free-algebra certificate. Unit and integration
tests use explicit algebraic identities, constructed proof graphs and temporary
solver records. Retired experiments and datasets remain in Git history.
