# Experimental variant-quotient calculations

These scripts explore dimensions, trace forms, radical candidates, tower maps and
Markov-trace conditions for specified untwisting quotients. They are research
programs outside the maintained installed package. Their historical inputs
and outputs are optional external research data listed in
[`data/research-archive.json`](../../data/research-archive.json).

| Script | Computation |
| --- | --- |
| `compute_dim_aho.py` | Leading-word avoidance counts with a declared alphabet. |
| `check_x3_nilpotent.py` | Specialized nilpotence search. |
| `semisimple_Q.py`, `semisimple_v3.py` | Trace-form calculations over the rationals or a finite field. |
| `trace_space.py`, `trace_space_fast.py` | Trace-space linear algebra in their stated coefficient domain. |
| `tower_check.py` | Tests of proposed maps at chosen parameters. |
| `markov_search.py` | Linear systems for candidate Markov traces. |
| `square_zero_*`, `explicit_rep_n3.py`, `generic_non_ss.py` | Specialized representation and ideal calculations. |

A dimension or normal-form conclusion requires a completed basis of the intended
ideal, with the coefficient field and specialization stated. Raising an engine
iteration limit does not prove completion. A calculation at several parameter
values does not establish the corresponding generic statement. In particular,
rank agreement modulo primes and generic non-semisimplicity require separate
justification.

Run recipes as `python -m experiments.variants.<module>` from this checkout
after installing the package. Solver recipes use `F4NCGB` or the `f4ncgb` on PATH.
See [TESTING](../../docs/TESTING.md) for the maintained test suite.

For file-based standard-word calculations, supply the declared strand count:

```sh
python -m experiments.variants.compute_dim_aho path/to/basis.txt --n 4 --max-deg 20
```

An explicit `--vars x1,x2,y1,y2` can supply another complete ordered alphabet
in the avoidance-count recipe. Missing generators cannot be recovered from the
terms that happen to occur in a basis file. Trace/rank recipes likewise take
`--n`; their quotient-algebra interpretation still requires completed-basis
evidence. The tower and Markov recipes explicitly use their stated three- and
four-strand alphabets.

The symbolic square-zero representation and trace-form recipes share the ignored
`runs/variants/` directory, separate from restored archives. They retain their legacy
Sympy-pickle cache format. Generate those local caches with the representation
producer before running the trace-form recipe:

```sh
python -m experiments.variants.square_zero_regular_rep_qt
python -m experiments.variants.square_zero_traceform_qt
```

Set `ZIPPER_EXPERIMENT_OUTPUT` to use another input and output directory for both
recipes. Historical results can separately be restored for provenance and
inspection with `python -m experiments.restore_data --group variants`. Restoration
copies verified bytes without loading Python objects. An exact JSON matrix export
is included in that archive under `data/variants/` for inspection without using
pickle. The recorded calculations retain the interpretation limits above.
