# Compact X₃ certificate experiments

These commands reproduce certificate compression and check its
braid-group-algebra identity. They use the installed `zipper_algebra` package
and run from the source checkout. They do not run a Gröbner-basis solver.

The input is the 695-term free-algebra certificate retained as
`tests/fixtures/legacy_x3/cofactor-X3-n4-Qqt-symbolic.txt`. The stages produce
1077, 166 and 161 terms. The last two expressions use equality in
`Q(q,t)[B_3]`; they are not literal free-word identities. Their recorded common
denominator is `t^4*(q-1)*(t+q^2)`.

```sh
python -m pip install -e .
python -m experiments.x3.compress_cofactor --output runs/x3-compact/compressed.txt
python -m experiments.x3.simplify_cofactor \
  --input runs/x3-compact/compressed.txt --output runs/x3-compact/simplified.txt
python -m experiments.x3.reduce_to_three_rels \
  --input runs/x3-compact/simplified.txt --output runs/x3-compact/three-rels.txt
python -m experiments.x3.verify_three_relations \
  --input runs/x3-compact/three-rels.txt
python -m experiments.x3.independent_verify_three_rels \
  --input runs/x3-compact/three-rels.txt
```

Each transformation checks its resulting identity. The two final commands
must exit with code zero. The last command supplies a separately written parser
and summation; both routes share presentation construction, polynomial
arithmetic and Garside normalization. The checks do not establish a claim
about replacing `T_3` with `Z_3` based on old specialized dimension calculations.

Compression defaults to the retained 695-term fixture. Every later stage and
checker requires an explicit input path; writers also require an explicit output
path. The generated 1077-, 166- and 161-term expressions are not tracked. The
commands above place them under the ignored `runs/x3-compact/` directory.

| Module | Purpose |
| --- | --- |
| `compress_cofactor.py` | Substitute checked derivations of selected defining relations. |
| `simplify_cofactor.py` | Normalize braid words and collect contextual terms; contains the retained B₃ Garside helpers. |
| `reduce_to_three_rels.py` | Eliminate the shifted two-strand relation by checked conjugation. |
| `verify_three_relations.py` | Check the compact identity and the computed ideal containments. |
| `independent_verify_three_rels.py` | Separately parse and sum the compact expression. |
| `legacy_format.py` | Read the historical four-factor text format and use maintained polynomial arithmetic. |
| `relations_naming.py` | Derive labels from the actual presentation. |
| `cancel_qm1.py` | Attempt exact cancellation in a supplied fixed expression; requires an output path. |
| `check_denom_minimal.py` | Print selected-factor and sampled-gcd diagnostics for an explicit input. These probes do not certify polynomial minimality. |

The old numeric extraction, grid collection, symbolic lifting and shell drivers
have been replaced by the maintained proof/reconstruction code and portable
execution layer. Their original implementations remain in the
[phase 1 source baseline](https://github.com/choomno/zipper-algebra/tree/0db1dd9e2a33fb0e5d07e216c8a5354911446e97/archive/cofactor-pre-consolidation).
In particular, the
[historical driver](https://github.com/choomno/zipper-algebra/blob/0db1dd9e2a33fb0e5d07e216c8a5354911446e97/archive/cofactor-pre-consolidation/reproduce_X3_zero.sh)
records the earlier target-appending mistake and corrected 31-relation emission.
The frozen legacy extraction comparison is retained under `tests/reference/`.

The two historical grid collectors assembled whole-cofactor coefficients.
Those data have a different role from the basis-proof coefficient reconstruction
used by the maintained producer.
