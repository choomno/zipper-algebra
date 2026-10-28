# Certificate workflow

The public package supports presentation construction, proof expansion,
contextual composition, exact reconstruction of candidate coefficients, and
verification of free-algebra membership identities. The
[mathematical specification](docs/mathematical-specification/main.tex) states the
objects, hypotheses and conclusions; this page describes how to use the software.

## Installation and presentation

Install with Python 3.11 or newer:

```sh
python -m pip install .
zipper-algebra describe --n 4
zipper-algebra describe --n 4 --json
```

A presentation includes the strand count, untwisting choice and its parameters,
level limits, optional zero-slider relations and named additional relations.
The symbolic relation content determines its hash. Identical readable labels
can belong to different presentations, so labels alone do not identify an algebra.
Use `zipper-algebra describe --help` for the available options.

## Checking an existing certificate

The retained historical free-algebra fixture can be checked without running a
Gröbner-basis solver:

```sh
zipper-algebra verify \
  logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt \
  --n 4 --target X_3 --pit --json
```

The checker validates the coefficient domain, denominator factors and exponents,
words, relation labels, presentation identity and target before arithmetic.
It compares the cleared polynomial identity exactly. `--pit` additionally checks
it on a product grid with degrees computed from the certificate and target.
The grid checker evaluates complementary denominator powers directly, so zeros
of denominator factors do not cause undefined division in the cleared identity.

Exit codes are:

| Code | Meaning |
| --- | --- |
| `0` | The requested identity checks passed. |
| `1` | The certificate is well formed but its identity is false. |
| `2` | The input, presentation, target or file is invalid. |

A generic membership identity specializes wherever its displayed denominators
remain nonzero. Agreement at a few parameter values is a different kind of
evidence. Certificates simplified using braid equivalence cannot be passed to a
literal free-word checker without accounting for the corresponding relations.
The large fixture above is retained for regression testing; it is not presented
as the latest mathematical proof.

## Producing candidates

The reusable stages have the following requirements and conclusions:

1. Construct the chosen presentation and record its content hash.
2. Parse a basis and its proof references, with explicit coefficient domain and
   variable order. Proof expansion permits forward references and rejects cycles.
3. Record a reduction skeleton with complete left and right word contexts.
4. Reconstruct candidate coefficients from exact samples. Fit rank, denominator
   validity and sample agreement are checked. Failure to fit is inconclusive.
5. Compose each contextual occurrence, retaining repeated uses of the same basis
   element and rejecting mismatched multiplier counts.
6. Verify the complete resulting identity exactly.

The installed package provides these algebraic stages. The historical solver
orchestration remains under `tools/cofactor/`; it also needs external binaries
and workspace configuration. It is not part of the installed command line
interface. A successful process exit or an output file alone does not certify a
completed basis. A nonzero remainder from an uncertified basis establishes no
nonmembership conclusion, even at the specialization being tested. The
experimental driver reports an inconclusive result with exit code `3` when it
has no verified certificate within its sampling cap.

## Historical transformations and reproduction

The [X3 reproduction guide](docs/X3_REPRODUCTION.md) gives the tested audit
commands, the complete file inventory and the status of external dependencies.

The original computational path used sampled solver proofs, lifted coefficients,
a large free-algebra certificate, and subsequent specialized compression and
braid-algebra transformations. Its scripts and scientific outputs are retained
in `tools/`, `logs_and_traces/` and `archive/`. Their dependencies and equality
conventions differ from those of the maintained checker.

Earlier versions of the ten-phase narrative are available in Git:

```sh
git log --oneline -- WORKFLOW.md
```

Those records include obsolete paths, unsupported generic conclusions and checks
that were not rerun during the present cleanup. Use it to trace provenance;
use the current test instructions to validate the maintained package.

## Development checks

```sh
python -m pip install -e .
python -m unittest discover -s tests -v
```

[TESTING](docs/TESTING.md) describes the CI suite and separate historical
compatibility checks. [CONTRIBUTING](CONTRIBUTING.md) records the requirements for
mathematical changes and independent adversarial review.
