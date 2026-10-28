# Independent Singular experiments

These scripts use Singular's `freegb.lib` for separate noncommutative calculations.
The repository's [MIT license](../LICENSE) applies to the scripts; Singular and
its libraries retain their own licenses.

## Specialized X3 membership

`x3-membership-Q.sing` uses the original 31 defining relations at `q=5`, `t=7`
over the rational numbers, with Letterplace degree bound 12. The six generators
are `x1,x2,x3,y1,y2,y3`. It computes `J = twostd(I)` and prints
`NF(X3, J)`. Run from this directory with Singular installed:

```sh
Singular -q x3-membership-Q.sing
```

The expected normal form is `0`. The recorded degree-12 run in the optional
`x3-history` archive group reports that reduction, followed by an unsuccessful
attempt to extract cofactors. Restore its log and numeric inputs with:

```sh
python -m experiments.restore_data --group x3-history --output .research-data
python -m experiments.restore_data --group inputs --output .research-data
```

Run these commands from the repository root. The log is then at
`.research-data/logs_and_traces/f4ncgb/singular_X3_Q_d12.log`. This
producer preserves the reduction stage and ends before that later attempt.
Its conclusion is specialized membership at `(5,7)`; the generic cofactor
identity has a separate [verification route](../docs/X3_REPRODUCTION.md).

The source is `archive/singular-failed/zipper-n4-X3-Q-db12.sing` at commit
`0db1dd9e2a33fb0e5d07e216c8a5354911446e97`. The presentation, target,
`twostd` call and reduction are unchanged. The defining relations match the
canonical specialized input, restored at
`.research-data/logs_and_traces/f4ncgb/zipper-n4-Q-q5-t7.ms`.
This extracted producer has been checked against those sources but has not
been rerun during the present cleanup.

The old directory also contained degree-bound 10 and 18 attempts and two
variations of degree-14 cofactor extraction. Their original `.sing`, `.out`
and `.err` files remain available in that commit's `archive/singular-failed/`
directory. Successful reduction observations and failed cofactor extraction
are distinct stages of those records.
