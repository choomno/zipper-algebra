# A small X_3 output example

`certificate.txt` is a 129-term two-sided ideal-membership certificate for
`X_3`, computed from the canonical four-strand presentation over `Q(q,t)`.
It is the single current output example retained for inspection without
installing or running the software. The historical data required for regression
tests are kept separately under `tests/fixtures/legacy_x3/`.

Each line records an exact rational coefficient, a left word, a labelled
defining relation and a right word. Their sum equals the target in the free
algebra. The displayed common denominator is
`q^5*t^2*(q-1)*(t+q^2)`; specialization of this expression requires it to be
nonzero. The mathematical reference describes this conclusion precisely.

After installation, check the example with:

```sh
zipper-algebra verify examples/x3/certificate.txt --n 4 --target X_3 --pit --json
```

Expected results are exit code `0`, overall `"pass": true`, 129 terms in the
exact report, and a passing `25` by `7` degree-bounded product grid. Timing
fields depend on the machine. Both checks use exact arithmetic.

The example came from a fresh local run of revision
`d30b36aa2da0819dfa28954bd9223ef72163edf5`, using parser-patched f4ncgb revision
`789504d7c7da6af6946faecb73421b037766a5b8`, 100 fresh parameter samples and
the original 31 defining relations. The presentation hash is
`091bbb92271ab870`. Its SHA-256 is
`f110a1aee32f8dde234bdc65fa0f57221d287655aaa930e50cc24ad6f20af071`.
It was independently rechecked symbolically and on the exact product grid.
The certificate is a witness, with no claim of minimality or uniqueness.

See the [reproduction guide](../../docs/X3_REPRODUCTION.md) for fresh production
and the [mathematical reference](../../docs/mathematical-specification/mathematical-specification.pdf)
for definitions and verification guarantees.
