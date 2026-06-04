# Workflow: extracting a symbolic two-sided cofactor over $\mathbb{Q}(q,t)$

This document records the computational pipeline used to produce the $X_3 = 0$ cofactor at $n=4$, focused on **what to run, in what order, and why**. The mathematical setting is in `REPORT.md` and the writeup is in `X3_cofactor_n4.tex`; this file is the operational recipe.

## Why the obvious approach fails

The straightforward idea — "run GAP/GBNP's `SGrobnerTrace` over $\mathbb{Q}(q,t)$, reduce $X_3$, extract the cofactor" — does not work, for two reasons:

1. **No general-purpose noncommutative Gröbner-basis engine supports working over a function field $\mathbb{Q}(q,t)$ directly.** GBNP and Singular's letterplace can be told to use $\mathbb{Q}$ as the coefficient ring but not $\mathbb{Q}(q,t)$ as a transcendental extension; f4ncgb supports prime fields and $\mathbb{Q}$ but not the function field either. Intermediate polynomial coefficients in $q, t$ grow rapidly during S-pair reductions, and exact arithmetic over $\mathbb{Q}(q,t)$ is prohibitively expensive even if it were implemented.

2. **Per-specialization runs don't naively interpolate.** Running over $\mathrm{GF}(p)$ with $(q,t)$ specialized to integers gives a cofactor over $\mathrm{GF}(p)$ at that one point. Running at multiple points and trying to read off a single $\mathbb{Q}(q,t)$-rational cofactor coefficient-by-coefficient generally **fails**, because different runs may select different representatives within the same syzygy fiber. The cofactor itself is not a well-defined function of $(q,t)$ — only its image in the ideal is.

The workflow below sidesteps both obstacles by lifting a **reduction skeleton** rather than the cofactor itself.

## The pipeline

```
            (Phase 1)       (Phase 2)       (Phase 3)       (Phase 4)
            single-anchor   skeleton        multi-point     symbolic
              GB+trace      extraction      proof lift      composition
                v               v               v               v
   inputs --> f4ncgb --> three-step --> rational lifts --> 695-tuple
   .ms file   trace        sequence      of 25 basis        Q(q,t) cofactor
                                         elements
                                                                v
                                                          (Phase 5)
                                                          verification
                                                          (free-algebra
                                                          poly identity)
                                                                v
                                                          (Phases 6-8)
                                                          algebraic
                                                          reductions
```

### Phase 1: single-anchor traced Gröbner basis

**Tool.** `f4ncgb` (Hofstadler & Heisinger 2025) with `--trace`.

**Input.** `.ms`-format file listing the relations of $\mathcal{Z}_n$ at a single specialization $(q_0, t_0)$ over $\mathrm{GF}(p)$ for a large prime $p$ (we used $p = 2^{31}-1$, $(q_0, t_0) = (5, 7)$).

**Generator.** `f4ncgb_tools/zipper_in_ms.py` — produces the `.ms` file from a Python implementation of the slider/zipper/untwisting construction. CLI:
```
python3 zipper_in_ms.py --n 4 --q 5 --t 7 --char 2147483647 \
    --untwist bigelow --target 3,1 --out inputs.ms
```
The `--target i,j` flag appends `S[i][j]` as the last "input relation" so the basis indices remain stable for the cofactor extractor.

**Output.**
- A Gröbner basis $B_0, \ldots, B_{136}$ (137 elements at $n=4$).
- A **non-expanded proof log** — for each $B_k$, a list of four-tuples $(c_{k,\ell}, L_{k,\ell}, i_{k,\ell}, R_{k,\ell})$ where each entry's index points either at another (earlier) basis element or at an original input relation. This is the key feature that distinguishes the f4ncgb pipeline from GBNP: GBNP's traces are *expanded all the way to inputs*, which loses the intermediate structure we exploit in Phase 3.

**Wall time.** Approximately 0.16 seconds at this anchor. Memory peak under 20 MB.

### Phase 2: reduction skeleton extraction

**Goal.** Express $X_3$ as a finite sequence of one-step reductions against the basis $B$, producing a tuple list $X_3 = \sum_\ell c_\ell L_\ell B_{k_\ell} R_\ell$ with each $c_\ell$ a **scalar in the prime field** (not yet in $\mathbb{Q}(q,t)$).

**Tool.** `f4ncgb_tools/extract_cofactor.py`. Reads the trace file from Phase 1, identifies the target, runs the one-step reductions, accumulates the four-tuples.

**Observed at $n=4$.** Three reduction steps:
$$
X_3 \;=\; q^3\, B_{17} \;-\; q^2\, B_{16} \;+\; q(q-1)\, B_{102}.
$$

**Crucial property.** The three multipliers $q^3, -q^2, q(q-1)$ are functions of $q$ alone (no $t$). This holds because the leading monomial coefficients of $X_3, B_{16}, B_{17}, B_{102}$ in the degree-lex order are themselves polynomials in $q$ — a consequence of the slider recursion being independent of $t$ (see §6.1 of `X3_cofactor_n4.tex`).

This property is what makes Phase 3 work. It is **algebra-dependent**; if you take a different target where the multipliers genuinely depend on $t$, this phase still produces a skeleton, but Phase 3 has to lift the multipliers themselves, which is significantly harder.

### Phase 3: multi-point lift of the proof closure

**Goal.** For each basis element $B_k$ appearing in the skeleton (transitively closed: 25 basis elements at $n=4$), recover a $\mathbb{Q}(q,t)$-rational expression for each coefficient $c_{k,\ell}$ in $B_k$'s non-expanded proof.

**Procedure.**
1. Run Phase 1 at $\sim 226$ different specializations $(q_0, t_0)$ on an integer grid chosen to avoid the candidate denominator vanishing loci $\{q=0, t=0, q-1=0, t-q=0, t+q^2=0, q+1=0\}$.
2. For each of the 25 basis elements in the proof closure, and each coefficient $c_{k,\ell}$ in its non-expanded proof, collect the 226 values.
3. Fit each value sequence to a rational function in $(q, t)$ by hypothesis-class search:
   - Search over candidate denominator factorizations on the atom set $\mathcal{A} = \{q, t, q-1, t-q, t+q^2, q+1\}$.
   - Use exact-rational three-split cross-validation: split the 226 points into train/validate/test, fit a candidate model on train, score on validate, and accept only if validate agrees exactly with test.
4. Output: a symbolic $\mathbb{Q}(q,t)$ expression for every coefficient in every proof of the 25 basis elements.

**Tools.** Several supporting scripts in the repo:
- `run_lift_grid.py`, `run_lift_grid_canonical.py` — Phase 1 driver across many specializations.
- `run_fit_polynomial.py`, `run_fit_classes_sweep.py`, `run_slice_fit_v*.py` — rational-function fitters with cross-validation.
- `extend_lift_grid.py`, `run_holdout_test.py` — grid expansion and holdout testing.

**Why this works where direct interpolation fails.** Phase 3 interpolates the *proof log entries* of the basis elements, not the cofactor of $X_3$ itself. The proof log entries are well-defined functions of $(q, t)$ in the appropriate region because f4ncgb's reduction is deterministic given a fixed monomial order and input ordering. By contrast, the cofactor of $X_3$ at different specializations can differ by syzygies (the "syzygy fiber" problem of §2.7 of the writeup).

**Required atom.** The fitting search will discover that the atom $t+q^2$ is required to fit the proof of $B_{68}$; the other five atoms have known algebraic meaning, $t+q^2$ does not (see §6.2 of `X3_cofactor_n4.tex`).

### Phase 4: symbolic composition

**Goal.** Substitute each lifted basis-element proof into the reduction skeleton from Phase 2 and multiply out.

**Tool.** `compose_x3_cofactor_qt.py`. Operates symbolically over $\mathbb{Q}(q,t)$ using exact bivariate-polynomial arithmetic (Python `fractions.Fraction` + dictionary representation; no floating point).

**Output.** `logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt`, a file with 695 lines, each declaring one tuple $(L_j, R_{i_j}, R'_j, c_j)$. Coefficients $c_j \in \mathbb{Q}(q,t)$ have denominators on the atom set $\mathcal{A}$ and integer numerators.

### Phase 5: verification

**Three independent verifications, no shared code:**

1. **Free-algebra polynomial identity** (the decisive one). Script `f4ncgb_tools/symbolic_identity_check.py`. Multiply out $\sum_j c_j(q,t) L_j R_{i_j} R'_j - X_3$ in $\mathbb{Q}(q,t)\langle x_i, y_i\rangle$ symbolically and check that the result is identically zero. No Gröbner basis machinery involved. Time: 0.13 s.

2. **Singular letterplace.** Translate to letterplace representation (`ms_to_singular.py`, `ms_to_singular_lpdiv.py`), run `twostd` to compute a two-sided Gröbner basis, then `reduce(X_3, GB)`. Done at letterplace degree bounds $\{10, 12, 14\}$ over both $\mathbb{Q}$ and $\mathrm{GF}(2^{31}-1)$.

3. **GAP/GBNP.** `SGrobner` on the inputs, then `StrongNormalFormNP(X_3, B)`. Done at $(q,t) \in \{(5,7), (11,13)\}$ over $\mathrm{GF}(2^{31}-1)$. Driven by `gap_scripts/zipper.gap`.

Independent of Phase 5, optional sanity checks include pointwise consistency at additional fresh specializations and a multi-point rational-interpolation cross-check (`run_independent_verify*.py`).

### Phase 6: compression onto an independent generating set

**Goal.** Eliminate derivable input relations from the cofactor.

**Tool.** `f4ncgb_tools/compress_cofactor.py`. Substitutes derivations of $R_7, R_{14}, R_{15}, R_{16}, R_{17}$ (the $y$-braid and the four mixed-sign length-3 relations) in terms of $\{R_0, R_1, R_2, R_3, R_6\}$.

**Output.** `cofactor-X3-n4-Qqt-compressed.txt`, 1077 tuples in 9 input indices $\{0, 1, 2, 3, 6, 22, 25, 26, 28\}$. Common denominator simplifies from $M = q^5 t^4 (q-1)(t+q^2)$ to $M' = t^4 (q-1)(t+q^2)$.

### Phase 7: passage to the group algebra of $B_n$ via Garside normalization

**Goal.** Replace each two-sided multiplier by its Garside normal form in the braid group, dropping tuples that lie in the invertibility + braid ideal.

**Tool.** `f4ncgb_tools/simplify_cofactor.py`. Implements left-greedy Garside normal form in $B_3$ ($\Delta = \sigma_1\sigma_2\sigma_1$, positive atoms $\{\sigma_1, \sigma_2, \sigma_1\sigma_2, \sigma_2\sigma_1\}$). Re-buckets tuples by their normalized multipliers and target index, merging duplicates and dropping cancellations.

**Output.** `cofactor-X3-n4-Qqt-simplified.txt`, 166 tuples in 4 slider relations $\{R_{22}, R_{25}, R_{26}, R_{28}\}$. The denominator $D = t^4(q-1)(t+q^2)$ is unchanged.

### Phase 8: conjugacy reduction

**Goal.** Eliminate $R_{26}$ using the $B_3$-conjugacy $R_{26} = (\sigma_1\sigma_2)\, R_{25}\, (\sigma_1\sigma_2)^{-1}$.

**Tool.** `f4ncgb_tools/reduce_to_three_rels.py`. Substitutes $(L, R_{26}, R) \mapsto (L \sigma_1 \sigma_2, R_{25}, \sigma_2^{-1}\sigma_1^{-1} R)$ on each $R_{26}$-tuple, re-normalizes in $B_3$, merges duplicates with existing $R_{25}$-tuples.

**Output.** `cofactor-X3-n4-Qqt-three-rels.txt`, 161 tuples in 3 input relations $\{R_{22}, R_{25}, R_{28}\}$.

### Optional: relabeling

`f4ncgb_tools/rename_to_ZT.py` renames the three relations to $\{Z_2, T_2, T_3\}$ (Z = zipper, T = untwisting, subscript = level). Output: `cofactor-X3-n4-Qqt-ZT.txt`.

## Reusability for a different target

The pipeline is parametric on **the target polynomial** ($X_3 = S[3][1]$ in our case) and **the relation set** ($\mathcal{Z}_n$ at fixed $n$). Re-running it for a different target or a different $n$ requires:

- **Phase 1**: regenerate the `.ms` file via `zipper_in_ms.py --target i,j --n N`. The same f4ncgb invocation works.
- **Phase 2**: rerun `extract_cofactor.py`. The reduction may take more or fewer steps depending on the target's depth.
- **Phase 3**: the multi-point grid and atom set may need expansion if the new target's proof closure introduces new denominator atoms. The fitter handles this if you broaden the atom search; in the worst case, you discover that the cofactor coefficients are not rational on a finite atom set, and the whole approach breaks (this did not happen at $n = 4$ for $X_3$).
- **Phase 4**: rerun `compose_x3_cofactor_qt.py`. No changes if the skeleton extraction succeeded.
- **Phases 5–8**: apply as-is. Phase 7's Garside normalization is currently specialized to $B_3$; for larger braid groups it generalizes (the Garside structure exists for every $B_n$) but the script would need to be extended.

The bottleneck-by-far is Phase 3, both in wall time (running f4ncgb at 226 points) and in fragility (if the atom set is wrong or the fitting fails on a coefficient, the symbolic lift cannot be assembled). When extending to new targets, **start by checking that the reduction skeleton is independent of $t$ (or of whichever parameter is hard to interpolate)**; if it is, the rest of the pipeline lifts cleanly.

## What the alternative would have looked like, and why it isn't better

A pure GBNP traced computation at $(q_0, t_0) = (5, 7)$ over $\mathrm{GF}(2^{31}-1)$ produces a basis of cardinality 89 (different from f4ncgb's 137 because of a different monomial order) together with a fully-expanded trace expressing each basis element directly over the original inputs. Reducing $X_3$ against this basis and unrolling the trace produces a cofactor over $\mathrm{GF}(2^{31}-1)$ at $(5,7)$ — a per-point certificate, not a $\mathbb{Q}(q,t)$ cofactor.

To get from this to a $\mathbb{Q}(q,t)$ cofactor, you would have to:
(a) Run GBNP at many $(q_0, t_0)$ values (slow: GBNP's Mora–Buchberger is ~10–100× slower than f4ncgb at our $n=4$ size).
(b) Lift each cofactor coefficient to $\mathbb{Q}(q,t)$ by multi-point fitting — **but** the cofactor coefficients are *not* well-defined functions of $(q,t)$ across runs because GBNP's expanded trace bakes in syzygy choices. You would have to find a stable structure (an analogue of our reduction skeleton + non-expanded proof closure) inside the trace; GBNP does not expose this.

Step (b) is the obstacle: it is the same obstacle that f4ncgb's non-expanded proof log was designed to solve. f4ncgb is the tool of record for this pipeline; GAP/GBNP is suitable only for verification at single $(q_0, t_0)$ points (Phase 5).

## File map

| Phase | Driver | Outputs |
|---|---|---|
| 1 | `f4ncgb_tools/zipper_in_ms.py` + `f4ncgb` binary | `*.ms` input, `*.trace` log |
| 2 | `f4ncgb_tools/extract_cofactor.py` | reduction skeleton |
| 3 | `run_lift_grid*.py`, `run_fit_*.py`, `run_slice_fit_v*.py` | per-coefficient $\mathbb{Q}(q,t)$ fits |
| 4 | `compose_x3_cofactor_qt.py` | `cofactor-X3-n4-Qqt-symbolic.txt` |
| 5 | `f4ncgb_tools/symbolic_identity_check.py`, `gap_scripts/zipper.gap`, Singular | verification log |
| 6 | `f4ncgb_tools/compress_cofactor.py` | `cofactor-X3-n4-Qqt-compressed.txt` |
| 7 | `f4ncgb_tools/simplify_cofactor.py` | `cofactor-X3-n4-Qqt-simplified.txt` |
| 8 | `f4ncgb_tools/reduce_to_three_rels.py` | `cofactor-X3-n4-Qqt-three-rels.txt` |
| --- | `f4ncgb_tools/rename_to_ZT.py` | `cofactor-X3-n4-Qqt-ZT.txt` |
