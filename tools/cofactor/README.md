# `cofactor/` — Project A: explicit $X_3 = 0$ cofactor pipeline

Computational tooling for the explicit two-sided $\mathbb{Q}(q,t)$ cofactor
of $X_3 = 0$ in $\mathcal{Z}_4 / \langle$canonical untwisting tower$\rangle$,
documented in `braid-research/X3_cofactor_n4.tex` and operationalised by
the top-level `WORKFLOW.md`. The pipeline writes
$$X_3 \;=\; \sum_{j=1}^{695} c_j(q,t)\; L_j\, R_{i_j}\, R'_j$$
in the free algebra, where each $c_j \in \mathbb{Q}(q,t)$ has denominator
in the atom set $\{q,\, t,\, q-1,\, t-q,\, t+q^2,\, q+1\}$ and each
$R_{i_j}$ is one of the defining relations of $\mathcal{Z}_4$.

After the most recent reorganisation, the entire pipeline lives in this
folder. The shared presentation emitter `zipper_in_ms.py` (Phase 1) lives
one directory up, at `tools/zipper_in_ms.py`.

## Phase-to-file map

| Phase | Tools | Output |
|---|---|---|
| **1.** Single-anchor traced GB | (uses `../zipper_in_ms.py` + `f4ncgb` binary externally) | `.ms` input, `.trace` log |
| **2.** Skeleton extraction | `extract_cofactor.py` | three-step skeleton $X_3 = q^3 B_{17} - q^2 B_{16} + q(q-1)B_{102}$ at the $(5,7)$ anchor |
| **3.** Multi-point lift (grid + fitters) | `run_lift_grid.py`, `run_lift_grid_canonical.py`, `extend_lift_grid.py`, `run_holdout_test.py`, `run_fit_polynomial.py`, `run_fit_classes_sweep.py`, `run_fit_coord_transforms.py`, `run_slice_fit.py`, `run_deep_sweep.py`, `run_stubborn_sweep.py`, `run_ultra_deep_sweep.py` | per-coefficient $\mathbb{Q}(q,t)$ fits |
| **4.** Symbolic composition | `compose_x3_cofactor_qt.py`, `write_honest_cofactor.py`, `write_unified_cofactor.py`, `extract_cofactor_qt.py` | `cofactor-X3-n4-Qqt-symbolic.txt` |
| **5.** Verification | `symbolic_identity_check.py` (decisive free-algebra identity), `verify_honest_at_qt.py`, `verify_unified_at_qt.py`, `verify_three_rels_at_23_29.py`, `independent_verify_three_rels.py`, `run_independent_verify.py`, `run_verify_symbolic_largeqt.py`, `merge_and_verify.py`; plus Singular path via `ms_to_singular.py`, `ms_to_singular_lpdiv.py`; plus GAP path in `../../gap_scripts/` | verification log |
| **6.** Compression onto independent generators | `compress_cofactor.py`, `simplify_cofactor.py` | `cofactor-X3-n4-Qqt-compressed.txt` |
| **7.** Passage to $\{Z_2, T_2, T_3\}$-form | `reduce_to_three_rels.py`, `rename_to_ZT.py`, `relations_naming.py`, `verify_ZT_and_ideal_equality.py`, `verify_cofactor.py` | `cofactor-X3-n4-Qqt-three-rels.txt`, `cofactor-X3-n4-Qqt-ZT.txt` |
| **8.** Garside normalisation auxiliaries | (in `simplify_cofactor.py`) | (in-place rewriting) |
| Auxiliary / diagnostics | `cancel_qm1.py`, `check_denom_minimal.py`, `spot_check_q31_t37.py`, `compare_skeletons.py`, `symbolic_lift_x3.py`, `test_n4_zipper_T2_f4ncgb.py` | denominator-atom optimisations, sanity probes |

## How the imports work

Every script in this folder runs as a standalone Python program (it has a
`main()` and `if __name__ == "__main__"`). They commonly import each other
(`from extract_cofactor import VAR_NAMES`, `from symbolic_identity_check
import ...`) and also `zipper_in_ms` from one directory up. A uniform
`sys.path`-setup block near the top of each script adds both `HERE` (this
folder) and `PARENT` (`tools/`, where `zipper_in_ms.py` and `core.py`
live) so that all these imports resolve regardless of the invocation cwd.

## Dependency cheat-sheet (most-depended-on first)

- `relations_naming.py` — labels for the 31 input relations of $\mathcal{Z}_4$.
- `symbolic_identity_check.py` — bivariate-polynomial NC engine + parser
  for the cofactor file format; imported by `cancel_qm1.py`,
  `compress_cofactor.py`, `reduce_to_three_rels.py`, `simplify_cofactor.py`,
  `spot_check_q31_t37.py`, `verify_three_rels_at_23_29.py`,
  `verify_ZT_and_ideal_equality.py`.
- `simplify_cofactor.py` — Garside normal form in $B_3$; imported by
  `cancel_qm1.py`, `verify_ZT_and_ideal_equality.py`,
  `verify_three_rels_at_23_29.py`.
- `extract_cofactor.py` — the Phase-2 entry point; also provides
  `VAR_NAMES` and `parse_noex_proof` used by `merge_and_verify.py`,
  `symbolic_lift_x3.py`, `verify_honest_at_qt.py`, `verify_unified_at_qt.py`.

## Reference

`braid-research/X3_cofactor_n4.tex` is the writeup; the parent-directory
`WORKFLOW.md` gives the operational recipe phase by phase.
