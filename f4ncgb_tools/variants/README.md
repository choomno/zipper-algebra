# Variant zipper algebras

Computational tooling for the variant-untwisting quotients of the
Bigelow zipper algebra studied in
`braid-research/variant_dimensions.tex`. Each script reads a Gröbner
basis emitted by `f4ncgb` and answers one question about the resulting
quotient algebra (dimension, Hilbert series, semisimplicity, trace
space, tower inclusion, Markov-trace space).

All scripts import shared utilities from `f4ncgb_tools/core.py` (the
parsers, the Aho–Corasick automaton, the noncommutative reducer).

## Pipeline

```
                                                                    char-0 path
                                                                    ↗
inputs --> zipper_in_ms --> f4ncgb --> .basis --> compute_dim_aho   trace_space
.ms (n, untwist, a)                              check_x3_nilpotent  semisimple_Q
                                                  semisimple_v3     trace_space_fast
                                                  trace_space_fast  tower_check
                                                                    markov_search
```

1. **`zipper_in_ms.py`** (at the parent `f4ncgb_tools/` level, shared)
   emits the `f4ncgb` `.ms` input for the variant of choice.

2. **`f4ncgb`** computes the noncommutative Gröbner basis. Run with
   `--maxiter 200` or higher — the default of 10 silently truncates the
   F4 saturation loop and produces incomplete bases (see the
   `--maxiter` discussion in `variant_dimensions.tex`).

3. Downstream analysis scripts (this folder):

   | Script | Question answered |
   |---|---|
   | `compute_dim_aho.py` | Hilbert series and total dim |
   | `check_x3_nilpotent.py` | Is $X_3$ nilpotent in the quotient? |
   | `semisimple_v3.py` | $F_p$ trace-form rank & det |
   | `semisimple_Q.py` | $\mathbb{Q}$ trace-form rank & det (char-0 theorem) |
   | `trace_space.py` | $\dim_\mathbb{Q} A/[A,A]$ |
   | `trace_space_fast.py` | $\dim_{F_p} A/[A,A]$, online row reduction |
   | `tower_check.py` | Does $Z_{n-1}^{\mathrm{var}} \hookrightarrow Z_n^{\mathrm{var}}$? |
   | `markov_search.py` | Markov-trace linear system on $Z_3 \to Z_4$ |

## Tests

```
python3 -m f4ncgb_tools.tests
```

The suite at `f4ncgb_tools/tests.py` pins the following numerical facts:

- `BMW_3` has dim 15 with 4 trace functionals.
- `BMW_3` is semisimple over $\mathbb{Q}$ ($\det T \ne 0$, rank 15).
- Jordan-3 has dim 24 with $\mathrm{rank}\,T = 15$ over $\mathbb{Q}$
  (so $\dim J = 9$).
- Round-trip parsers and Aho–Corasick automaton match expectations.

Tests must pass after any change to `core.py` or to a variant script.

## Reference

See `braid-research/variant_dimensions.tex` for the full mathematical
context (algebraic background, results, semisimplicity proofs,
Markov-trace setup). The §"Reproducible computational protocols"
section there gives the exact command sequence for each result.
