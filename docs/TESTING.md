# Testing the certificate boundary

From the repository root, run:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

This suite uses Python's standard library. It requires no solver, sibling
research folders, personal settings, pre-existing caches or network access.
It tests certificate validation, symbolic and grid verification, contextual
composition, repeated multiplier occurrences, sample alignment, package imports
and compatibility with the historical free-algebra certificate. Test output is
written only to temporary directories. The longer grid check typically takes
several seconds.

The independent historical proof-expansion comparison can also be run:

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest tools.cofactor.tests_proof_skeleton -v
```

It uses the archived implementation as a reference and the immutable fixtures
under `tests/fixtures/legacy_x3/`. Their hashes and original source records are
listed in the adjacent manifest. The existing large certificate remains in its
historical artifact directory. Compatibility checks preserve its body bytes;
they do not identify it as the latest mathematical proof.

The older `python3 -m tools.tests` suite runs solver integrations and variant
calculations. It still depends on the former workspace configuration and an
installed solver. It is excluded from the default suite pending the execution
and configuration audit. Reconstruction tests also need their own fixture
repair and correctness work. Passing the certificate suite does not imply that
these optional integrations have run or that all research scripts are supported.

See the [mathematical specification](mathematical-specification/main.tex) for
the domains, assumptions and exact conclusions of the maintained checks.
