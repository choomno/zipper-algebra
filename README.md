# `zipper-algebra` — computational pipelines for the Bigelow zipper algebra

Phase-one work is recorded in [the public scope inventory](docs/PUBLIC_SCOPE.md),
[the mathematical specification](docs/mathematical-specification/main.tex), and
[the current test instructions](docs/TESTING.md). The maintained certificate
checks run without an installed solver or the surrounding private workspace.
The historical descriptions below retain superseded commands and claims and are
being reconciled with the supported software boundary.

The current objective is complete papers and reusable public code. See
[`RESEARCH_AGENDA.md`](RESEARCH_AGENDA.md) for the research scope, unresolved questions,
source recovery, and release criteria. The historical project descriptions follow.

> **Trust tiers (see `INDEX.md` for the full map).**
> 1. **Shared/trusted** — `tools/{zipper_in_ms,core,tests}.py`.
> 2. **Project A** (frozen main tool) — `tools/cofactor/`, the exact ℚ(q,t) `X_3=0`
>    cofactor pipeline; `reproduce_X3_zero.sh` reruns it.
> 3. **Project B** (results **distrusted**, code partly reusable) — `tools/variants/`;
>    see `AUDIT_spec_vs_generic.md` / `LEDGER_claims.md`.
> 4. **Research phase** (2026-07+) — `tools/research/`, `gap_scripts/theta3_*.gap`;
>    see `NIGHT_LOG.md`.
>
> Superseded material lives under `archive/`. **Exact char-0 (ℚ(q,t) symbolic, or the
> f4ncgb grid + rational reconstruction) is load-bearing; a lone mod-p spec is a scout only.**

Two distinct computational projects share this repository. They overlap in
infrastructure (the Gröbner-basis machinery built on `f4ncgb`) but address
distinct mathematical questions. After the most recent reorganisation each
project lives in a single subfolder of `tools/`, with shared
utilities at `tools/`'s root.

## Project A — explicit $X_3 = 0$ cofactor

**Goal.** Produce a fully-symbolic two-sided cofactor expressing $X_3 = 0$
over $\mathbb{Q}(q,t)\bigl[D^{-1}\bigr]$ in $\widetilde{\mathcal{Z}}_4$
(the canonical untwisting quotient of the zipper algebra). Operationally
a phase-by-phase pipeline going from a single-anchor `f4ncgb` Gröbner
basis to a symbolic $\mathbb{Q}(q,t)$ tuple list whose summation reduces
to $X_3$ in the free algebra.

**Documentation.**
- `WORKFLOW.md` — phase-by-phase operational recipe.
- `braid-research/X3_cofactor_n4.tex` — formal writeup.

**Code (all in `tools/cofactor/`).** See its README for the
phase-to-file mapping; the headlines are:

| Phase | Tools |
|---|---|
| 1 (input emission)              | `tools/zipper_in_ms.py` (shared) + `f4ncgb` binary |
| 2 (skeleton extraction)         | `cofactor/extract_cofactor.py` |
| 3 (multi-point lift, fitters)   | `cofactor/run_lift_grid*.py`, `cofactor/run_fit_*.py`, `cofactor/run_slice_fit.py`, `cofactor/extend_lift_grid.py`, `cofactor/run_holdout_test.py` |
| 4 (symbolic composition)        | `cofactor/compose_x3_cofactor_qt.py`, `cofactor/write_*_cofactor.py` |
| 5 (verification)                | `cofactor/symbolic_identity_check.py`, `cofactor/verify_*.py`, `cofactor/independent_verify_three_rels.py`, `gap_scripts/` (independent GAP verification path), `singular_run/` (Singular letterplace path) |
| 6 (compression)                 | `cofactor/compress_cofactor.py`, `cofactor/simplify_cofactor.py` |
| 7 (passage to ZT-form)          | `cofactor/rename_to_ZT.py`, `cofactor/reduce_to_three_rels.py`, `cofactor/relations_naming.py`, `cofactor/verify_ZT_and_ideal_equality.py` |
| auxiliary                       | `cofactor/cancel_qm1.py`, `cofactor/check_denom_minimal.py`, `cofactor/spot_check_q31_t37.py`, `cofactor/compare_skeletons.py` |
| shell driver                    | `reproduce_X3_zero.sh` |

**Data.**
- `ms_inputs/qt_sweep/` — ~230 `.ms` inputs for Phase 3's multi-point grid.
- `gap_scripts/` — GAP/GBNP scripts for Phase 5's independent verification path.
- `singular_run/` — letterplace Singular outputs (Phase 5).
- `cluster/` — cluster job descriptions used to parallelise Phase 3.
- `logs_and_traces/f4ncgb/` — the trace logs and the final cofactor files.
- `logs_and_traces/gap_traces/` — GAP traces.

## Project B — variant zipper algebras

**Goal.** Compute dimensions, semisimplicity, trace spaces and Markov-trace
linear systems for the variant-untwisting quotients of $\mathcal{Z}_n$ at
small $n$ (higher-framing $\sigma_j^a X_2 = t X_2$ and square-zero
$(\sigma_j - t)^2 X_2 = 0$). The variants strictly enlarge BMW while
remaining finite-dimensional; the framing case is generically semisimple,
the square-zero case is generically non-semisimple (the first such example
found, with $\dim J = \dim - \dim\mathrm{BMW}_n$ at $n=3,4$).

**Documentation.**
- `braid-research/variant_dimensions.tex` — formal writeup. Contains
  algebraic background for a topology audience, the trace-form
  criterion, computed dimensions, conjectures, and the §"Reproducible
  computational protocols" with exact command sequences for each result.

**Code.**
- `tools/zipper_in_ms.py` — presentation emitter (shared).
- `tools/core.py` — shared utilities (parsers, Aho–Corasick, reducer).
- `tools/variants/` — per-question scripts; see its README.

**Data.**
- `logs_and_traces/variants/` — `.basis` artifacts from the verified runs.

## Shared infrastructure (`tools/`)

- `core.py` — single source of truth for `.basis` parsing (F_p and Q),
  Aho–Corasick construction on the leading-monomial set, the standard-
  monomial enumerator, and the noncommutative reducer.
- `zipper_in_ms.py` — presentation emitter; CLI flags select the variant
  (canonical Bigelow, higher framing, silly-question, square-zero), the
  characteristic, and which higher-level relations to include.
- `tests.py` — regression suite. Run with `python3 -m tools.tests`;
  12 tests covering parser correctness, automaton invariants, reduction
  correctness in both rings, and end-to-end verdicts on BMW_3 and square-zero-3.

## Running the tests

```bash
python3 -m tools.tests
```

All 12 regression tests should pass in well under a second on a laptop.

## Quick start: rebuild one variant Hilbert series

```bash
python3 tools/zipper_in_ms.py --n 4 --q 5 --t 7 --char 2147483647 \
    --untwist bigelow --framing-power 2 --out /tmp/example.ms
/path/to/f4ncgb --maxiter 500 -o /tmp/example.basis /tmp/example.ms
python3 tools/variants/compute_dim_aho.py /tmp/example.basis
```

## Coordinates

Both projects work at $(q, t) = (5, 7)$ in $\mathbb{F}_p$ with
$p = 2^{31} - 1$ as their default specialisation; both also accept
`--char 0` to lift to $\mathbb{Q}$ for the rigorous statements (Dieudonné's
trace-form criterion, semisimplicity verdicts, etc.).

## Repository layout (current)

```
zipper-algebra/
├── README.md                        ← this file
├── WORKFLOW.md                      ← Project A complete record (math + pipeline + verification)
├── reproduce_X3_zero.sh             ← Project A reproducibility driver
├── tools/
│   ├── README.md                    ← pipeline overview
│   ├── core.py                      ← shared utilities (parsers, AC, reducer)
│   ├── zipper_in_ms.py              ← presentation emitter (shared)
│   ├── tests.py                     ← regression suite (12 tests)
│   ├── cofactor/                    ← Project A code
│   │   ├── README.md
│   │   └── 40 .py files
│   └── variants/                    ← Project B code
│       ├── README.md
│       └── 8 .py files
├── ms_inputs/qt_sweep/              ← Project A multi-point grid inputs (~230 .ms)
│   └── README.md
├── gap_scripts/                     ← GAP/GBNP scripts (Project A, untouched .gap)
├── singular_run/                    ← Singular letterplace outputs (Project A)
├── cluster/                         ← cluster job specs
└── logs_and_traces/
    ├── f4ncgb/                      ← Project A traces and final cofactor files
    ├── gap_traces/                  ← Project A GAP traces
    └── variants/                    ← Project B .basis files
```
