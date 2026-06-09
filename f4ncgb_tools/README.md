# `f4ncgb_tools/`

Python pipeline around the `f4ncgb` noncommutative Gröbner-basis engine
(Hofstadler & Heisinger, 2025). Two distinct workflows, each in its own
subfolder, share a single core utility module.

| File / folder | Purpose |
|---|---|
| `core.py`            | Shared utilities: `.basis` parser (F_p, Q), Aho–Corasick automaton, standard-monomial enumerator, polynomial reducer. The single source of truth for Gröbner-machinery routines. |
| `zipper_in_ms.py`    | Presentation emitter for the zipper-algebra family. Writes an `.ms`-format input for `f4ncgb` from a tuple `(n, untwist, framing-power, theta-power, char, …)`. Used by both pipelines. |
| `tests.py`           | Regression test suite. Pins numerical facts (BMW_3 dim 15, Jordan-3 rank 15 over Q, …). Run with `python3 -m f4ncgb_tools.tests`. |
| `cofactor/`          | **Project A — explicit $X_3 = 0$ cofactor pipeline.** Phase 2 (skeleton extraction) through Phase 8 (passage to ZT-form), plus the Phase 3 grid drivers and fitters, the Phase 4 composers, and the Phase 5 verifiers. See `cofactor/README.md`. |
| `variants/`          | **Project B — variant zipper algebras.** Hilbert-series, semisimplicity, trace-space, tower-inclusion, and Markov-trace tools for the higher-framing and Jordan-block variants of the canonical untwisting tower. See `variants/README.md`. |

## Quick start

```bash
# Generate a presentation, build the GB, compute its dimension.
python3 f4ncgb_tools/zipper_in_ms.py --n 4 --q 5 --t 7 --char 2147483647 \
    --untwist bigelow --framing-power 2 --out /tmp/example.ms
/path/to/f4ncgb --maxiter 500 -o /tmp/example.basis /tmp/example.ms
python3 f4ncgb_tools/variants/compute_dim_aho.py /tmp/example.basis
```

## Running the tests

```bash
python3 -m f4ncgb_tools.tests
```

The fixture-setup step generates four small bases (BMW_3 and Jordan-3,
each in F_p and Q) — about half a second in total. Tests must pass after
any change to `core.py` or any moved/edited script; they're the regression
guard against future refactors.

## Dependencies

- Python 3.11+ (uses `int | None`-style annotations and `Fraction`).
- The `f4ncgb` binary (Hofstadler & Heisinger 2025) somewhere on disk.
  Scripts default to `/home/choom/src/f4ncgb/build-release/f4ncgb` and
  accept a `--f4ncgb` override on the relevant CLIs.
- No external Python libraries; only the standard library
  (`fractions`, `collections`, `re`, `subprocess`, `unittest`).

## Conventions

- Coefficient rings are `F_p` with `p = 2^31 − 1` by default, or `Q` via
  `--char 0`. The variant pipeline uses both; the cofactor pipeline uses
  `Q` end-to-end after the Phase-1 anchor.
- The presentation emitter (`zipper_in_ms.py`) generates `.ms` files with a
  fixed monomial order (length–lex). The Aho–Corasick automaton built on
  the leading-monomial set is the central data structure in `core.py`.
- All sys-path bookkeeping in `cofactor/` scripts is uniform and handled by
  a 6-line block at the top of each script (see `cofactor/README.md`); to
  add a new cofactor script that imports `zipper_in_ms` or a sibling,
  copy that block verbatim.
