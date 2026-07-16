# INDEX — zipper-algebra repository map

Living index (started 2026-07-16 during the consolidation/audit). Four trust tiers.
Update on every add/change. Snapshot before the consolidation:
`~/.zipper-algebra-backups/pre-audit-20260716-001613.tar.gz`.

## Tier 1 — shared, trusted (the going-forward substrate)
| path | what | notes |
|------|------|-------|
| `f4ncgb_tools/zipper_in_ms.py` | presentation emitter (.ms for f4ncgb) | `--max-untwist-level 2` = two-strand-only; imported everywhere |
| `f4ncgb_tools/core.py` | .basis parsers, Aho-Corasick, standard-monomial enum, reducers | single source of truth for GB post-processing |
| `f4ncgb_tools/tests.py` | regression suite (BMW_3, Jordan-3 fixtures) | hardcodes the f4ncgb binary path |

## Tier 2 — Project A: the exact Q(q,t) X_3=0 cofactor pipeline (MAIN TOOL, frozen)
| path | what |
|------|------|
| `f4ncgb_tools/cofactor/` | 8-phase pipeline: extract -> lift-grid -> symbolic-lift -> compose -> verify -> compress -> simplify -> reduce-to-3-rels. See its README. Do not move/edit internals. |
| `reproduce_X3_zero.sh`, `WORKFLOW.md` | reproducibility driver + complete record. Reads cached `logs_and_traces/f4ncgb/basis-symbolic-Qqt.json`. |
| `logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt*.txt` | the exact result (12 renderings; `-final`/main/`-ZT`/`-minimal` canonical) + `basis-symbolic-Qqt.{json,txt}` + singular-independent verifications |
| `logs_and_traces/f4ncgb/lift_samples*/coefficients.json` | reconstruction caches (kept; the per-spec sample bulk was deleted, regenerable) |
| `singular_run/jordan-n3-symbolic.sing` | independent Singular letterplace over Q(q,t) (kept; verify before relying) |

## Tier 3 — Project B: variant-untwisting analysis (results DISTRUSTED, code partly reusable)
See `AUDIT_spec_vs_generic.md` + `LEDGER_claims.md` for which claims lift to generic.
| path | status |
|------|--------|
| `f4ncgb_tools/variants/{compute_dim_aho,semisimple_Q,trace_space,tower_check,markov_search}.py` | exact-char-0 / char-free machinery — reusable |
| `f4ncgb_tools/variants/{jordan_regular_rep_qt,jordan_traceform_qt,generic_non_ss,jordan_ideal_nilpotence,explicit_rep_n3}.py` | use the trusted grid+reconstruction / symbolic methodology — reusable |
| `f4ncgb_tools/variants/{semisimple_v3,trace_space_fast}.py` | F_p-only fast counterparts — pair with the Q twin |
| `f4ncgb_tools/variants/{check_x3_nilpotent,sigma_min_poly}.py` | superseded/mod-p-only — kept in place (tests import siblings) but not load-bearing |

## Tier 4 — research phase (2026-07+, new)
| path | what |
|------|------|
| `f4ncgb_tools/research/theta3_min_poly.py` | derived 3-strand untwisting = min poly of Theta_3 on X_3 (mod-p grid scout). **Finding: = p*, not (lambda-q^2)p*.** See its README. |
| `gap_scripts/theta3_exact.gap` | EXACT Q(q,t) confirmation of the above (GAP FunctionField) |
| `gap_scripts/{slider_module_dim,recursion_theta,recursion_theta_jordan,recursion_theta_framing}.gap` | GAP module/min-poly scouts (mod-p) |
| `gap_scripts/tests/{hecke,bmw}.gap` | GBNP fixtures (hecke = exact FunctionField) |
| `NIGHT_LOG.md` | running log of the 2026-07-16 autonomous session |

## archive/ — superseded, kept as record (not deleted; permanent local snapshot exists)
`archive/gap/` (old GBNP driver path + saved GB dumps), `archive/cluster/` (SLURM specs),
`archive/singular-failed/` (failed X_3 Singular attempts), `archive/variants-data/` (Project-B
jordan regrep/traceform pickles).

## Not yet done (refactor TODOs)
- `corelib/` extraction (normalform, fieldlinalg from sigma_min_poly, reconstruct, f4ncgb_run).
- `paths.py` (central binary/paths). `--untwist-scope` emitter alias. `README` full rewrite.
- Fix stale reorg paths in `cofactor/{symbolic_lift_x3,compose_x3_cofactor_qt,run_lift_grid}.py`.
