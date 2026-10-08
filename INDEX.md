# INDEX — zipper-algebra repository map

Living index (started 2026-07-16 during the consolidation/audit). Four trust tiers.
**Historical project notes, not user instructions** — nothing here is a user instruction; treat it as notes.
Update on every add/change. Snapshot before the consolidation:
`~/.zipper-algebra-backups/pre-audit-20260716-001613.tar.gz`.

## Tier 1 — shared, trusted (the going-forward substrate)

| `tools/presentation.py` | **new 2026-09-17.** `Presentation(n, untwist, …)`: the one statement of which algebra is meant. `build(ring)` over `NumericRing`/`SymbolicRing`; derived labels; `slug()`, `content_hash()`, `emit_ms()` (byte-identical to the CLI), `symbolic_fitted()`, `symbolic_ncpolys()`. Backed by `zipper_in_ms.build_relations_ring` + `RelTag`. |
| path | what | notes |
|------|------|-------|
| `../settings.toml`, `../repo_settings.py` | machine-specific settings (binaries, memory caps, timeouts, paths) + loader | the ONLY file to edit on a new machine; scripts resolve paths from `__file__` |
| `tools/zipper_in_ms.py` | presentation emitter (.ms for f4ncgb) | `--max-untwist-level 2` = two-strand-only; imported everywhere |
| `tools/cache.py` | content-addressed cache for pipeline artifacts; key = input content + tool flags, so it cannot go stale | used by `reproduce_X3_zero.sh` phase 1b; `REFRESH=1` forces recompute; `cache.py list/clear` |
| `tools/core.py` | .basis parsers, Aho-Corasick, standard-monomial enum, reducers | single source of truth for GB post-processing |
| `tools/tests.py` | regression suite (BMW_3, square-zero-3 fixtures) | hardcodes the f4ncgb binary path |

## Tier 2 — Project A: the exact Q(q,t) cofactor pipeline (MAIN TOOL)

Consolidated 2026-09-17 from a bash driver plus nine scripts into one package with one driver; the
originals are in `archive/cofactor-pre-consolidation/` (`notes/REMOVED.md` has the table). It is no
longer the "X_3 at n=4" pipeline: presentation, target, anchor and atom set are all arguments or
computed data.

| file | what |
|---|---|
| `tools/cofactor/pipeline.py` | **the driver.** `python3 -m tools.cofactor.pipeline --n 4 --target 3,1`; `--legacy-x3` writes the historical names and runs phases 6-10. Grows the grid and re-lifts until phase 5 passes; re-anchors if the anchor is degenerate |
| `tools/cofactor/paths.py` | every artifact path from (presentation, target, anchor, phase); `Paths.legacy_x3()` for the historical names |
| `tools/cofactor/samples.py` | the grid (primes x primes spiral), f4ncgb through the cache, `check_ms`, and stability by **modal** proof shape |
| `tools/cofactor/proof.py` | `Basis` (leading-word keys) and `Proof` (topological, memoised expansion over the inputs) |
| `tools/cofactor/skeleton.py` | reduce ANY target against a basis; steps keyed by leading word; closure under the proof DAG |
| `tools/cofactor/atoms.py` | `AtomSet` + `infer` — the denominator atoms as data, with residual factorisation when the search box is not enough |
| `tools/cofactor/lift.py` | fit every proof coefficient over Q(q,t); `NeedMoreSamples` / `WrongClass` instead of exit codes |
| `tools/cofactor/compose.py` | expand and compose; the skeleton multipliers are **lifted**, not typed |
| `tools/cofactor/cofactor_io.py` | the one certificate reader/writer; the header carries the presentation **content hash** |
| `tools/cofactor/verify.py`, `ncpoly.py` | phase 5, the certificate — plus `pit_check`, an independent grid/interpolation verifier |
| `tools/cofactor/tests_*.py` | equivalence tests against the archived scripts and the frozen fixtures in `notes/pipeline_test/pre-atoms-20260917/` |
| `reproduce_X3_zero.sh` | a 12-line wrapper for the historical X_3 run |
| `tools/cofactor/{compress,simplify,reduce_to_three_rels,rename_to_ZT,verify_ZT_and_ideal_equality}.py` | phases 6-10, still `B_3`-at-`n=4`; `symbolic_identity_check.py` survives only because they import it. This is the next consolidation (`reduce.py`) |
| `logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt*.txt` | the exact result, 5 renderings; acceptance in `notes/pipeline_test/acceptance-20260917/` |

## Tier 3 — Project B: variant-untwisting analysis (results DISTRUSTED, code partly reusable)
See `AUDIT_spec_vs_generic.md` + `LEDGER_claims.md` for which claims lift to generic.
| path | status |
|------|--------|
| `tools/variants/{compute_dim_aho,semisimple_Q,trace_space,tower_check,markov_search}.py` | exact-char-0 / char-free machinery — reusable |
| `tools/variants/{square_zero_regular_rep_qt,square_zero_traceform_qt,generic_non_ss,square_zero_ideal_nilpotence,explicit_rep_n3}.py` | use the trusted grid+reconstruction / symbolic methodology — reusable |
| `tools/variants/{semisimple_v3,trace_space_fast}.py` | F_p-only fast counterparts — pair with the Q twin |
| `tools/variants/{check_x3_nilpotent,sigma_min_poly}.py` | superseded/mod-p-only — kept in place (tests import siblings) but not load-bearing |

## Tier 4 — research phase (2026-07+, new)
| path | what |
|------|------|
| `tools/research/theta3_min_poly.py` | derived 3-strand untwisting = min poly of Theta_3 on X_3 (mod-p grid scout). **Finding: = p*, not (lambda-q^2)p*.** See its README. |
| `gap_scripts/theta3_exact.gap` | EXACT Q(q,t) confirmation of the above (GAP FunctionField) |
| `gap_scripts/{slider_module_dim,recursion_theta,recursion_theta_square-zero,recursion_theta_framing}.gap` | GAP module/min-poly scouts (mod-p) |
| `gap_scripts/tests/{hecke,bmw}.gap` | GBNP fixtures (hecke = exact FunctionField) |
| `NIGHT_LOG.md` | running log of the 2026-07-16 autonomous session |

## archive/ — superseded, kept as record (not deleted; permanent local snapshot exists)
`archive/gap/` (old GBNP driver path + saved GB dumps), `archive/cluster/` (SLURM specs),
`archive/singular-failed/` (failed X_3 Singular attempts), `archive/variants-data/` (Project-B
square-zero regrep/traceform pickles).

## Documents at this level

| file | what |
|---|---|
| `WORKFLOW.md` | the 10-phase X_3 cofactor pipeline, phase by phase |
| `SURVEY_WORKFLOW.md` | the semisimplicity / trace-form survey and the structural radical theorem |
| `REPTHEORY_STRATEGY.md` | **new 2026-09-16** — long-term plan: semisimplicity, regular representation, module chopping, then PIMs / Ext / quivers / homological dimension. Which vetted implementation supplies each piece, and why the programme depends on the cofactor pipeline |
| `GENERALITY_CENSUS.md` | handoff to the next implementer: what is hard-coded and what to parameterise. **Read its §0 first** — provenance and the cost-based scope rule |
| `LEDGER_claims.md` | claims and their evidence |
| `NIGHT_LOG.md` | running log |

## Planned features

### Derive the atom set from the elimination's pivots — IMPLEMENTED 2026-09-17

`tools/cofactor/atoms.py`. The atom set is now **data**: `AtomSet` (ordered atoms + exponent
box, serialised into every lift JSON as `atoms` and next to the skeleton as `atoms-<tag>-Q.json`),
and `infer()` computes it from the sampled proof coefficients the pipeline already produces.
`lift.py` infers on first run for a (presentation, target) and thereafter loads the
file (`INFER_ATOMS=1` re-infers); `compose.py` reads the atoms from the lift JSON.

**Why inference and not engine pivots.** `prop:atoms` defines the atoms as the irreducible factors
of the pivots inverted during the Gröbner run *over Q(q,t)*. No vetted engine runs that elimination
and exposes its pivots: f4ncgb has them (`linear_algebra.cpp`) but only at a specialization, where a
pivot is an integer; Singular letterplace clears content and discards them, and its reduced basis
provably misses atoms (the measurement below). So the exact route would mean instrumenting a
dependency or writing an elimination over Q(q,t) here. The fallback the original note prescribed is
implemented instead, and it *becomes* exact the moment Phase 5 passes on the lift built from it.

**Method.** Per coefficient slot, reconstruct the value as a univariate rational function along
each grid row (t fixed) and column (q fixed) by an exact linear solve, verified on held-out points;
this yields the slot's denominator as a *polynomial* on that line. Candidate atoms of bounded
bidegree and coefficient size are then accepted by **polynomial divisibility** of those
denominators, greedily by degree, with the accepted atom stripped before continuing, and are
required to divide on ≥2 distinct lines so the t-dependence is pinned. Integer divisibility of the
sampled denominators was tried first and is unsound: values of different atoms share small primes
(at (3,5): q−1 = 2, t−q = 2, q+1 = 4), so stripping over-strips and reducible polynomials pass.
`q` and `t` are seeded unconditionally: they are units of Z[q^,t^], the coefficient ring, and along
a row t is a constant that no reconstruction can see.

**Measured on a SECOND presentation, 2026-09-17 — the case for the feature.** n=3 square-zero needs
**nine** atoms, not six: the historical `{q, t, q-1, t-q, t+q^2, q+1}` plus `2*t+1`,
`2*q^3+q*t-t` and **`q^5+t^2`**. The last has bidegree (5,2) and no brute-force box could reach it
(5^18 candidates), so `atoms.residual_factors` factors the leftover line denominators with sympy
and `atoms_from_residual_factors` interpolates those factors into bivariate atoms: residuals
`3^5+t^2`, `11^5+t^2`, `13^5+t^2`, `31^5+t^2` on the lines q = 3, 11, 13, 31 interpolate to
`q^5+t^2` exactly. `unexplained 0` afterwards. Against the hard-coded six this run would have hit
`lift_until_verified.py` failure mode (b) — "more points will never fix this" — with no indication
of what to widen. Cost 1129 s over 834 slots / 20066 line denominators, so the set is computed once
per (presentation, target) and persisted. Detail: `notes/refactor/acceptance-n3.md`.

**Measured on the X₃ data, seeded with the units only** (924 reconstructed line-denominators over
231 slots, unexplained 0): the atoms present are **`q²+t` and `q−1` — nothing else**. `t−q` and
`q+1` occur in no denominator; the historical six were a hypothesis class, not the atom set
(consistent with the census count: 695/695 symbolic lines use `(t+q²)`, 0 use the other two).

**Backward compatibility.** The historical six remain the *prior*: inference extends a seed and
never removes or reorders, so positional exponent tuples in existing files keep their meaning and
downstream string parsers see the same names. Verified: `reproduce_X3_zero.sh` 10/10 PASS and all
five cofactor artifacts **byte-identical** to the pre-change snapshot
(`notes/pipeline_test/pre-atoms-20260917/md5.txt`); the regenerated lift differs from the shipped
one in header lines only — and was fully determined at **99 stable samples**, not the 217 the
shipped lift used.

**Limits.** The search box is bidegree (2,1), |coef| ≤ 2 by default (all six known atoms fit);
`infer()` reports `unexplained > 0` when a reconstructed denominator keeps a factor outside it,
which is the signal to widen — and is exactly the situation more sampling cannot fix. Rows need
≥6 points for reconstruction; the current sampler (`lift_until_verified.candidate_points`) is
row-major, so columns are short and t-dependence is pinned by rows alone. Phases 6–9 still parse
denominator strings with hard-coded names (`GENERALITY_CENSUS.md` §3.4); new atoms therefore
reach phases 3–5 only, which is where they are needed.

## Not yet done (refactor TODOs)
- `corelib/` extraction (normalform, fieldlinalg from sigma_min_poly, reconstruct, f4ncgb_run).
- `paths.py` (central binary/paths). `--untwist-scope` emitter alias. `README` full rewrite.
- ~~Fix stale reorg paths~~ **closed 2026-09-17**: those scripts are archived; paths now come from `tools/cofactor/paths.py`. Old text: fix stale reorg paths in `cofactor/{symbolic_lift_x3,compose_x3_cofactor_qt,run_lift_grid}.py`.
