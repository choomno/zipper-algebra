# research/ — quarantined research-phase computations (2026-07+)

New, exact-char-0-oriented work. Imports only `core`, `zipper_in_ms`, and (currently)
`variants.sigma_min_poly` for the F_p linear-algebra helpers (TODO: hoist those into
`corelib/`). It does **not** depend on `cofactor/`. Trust the results here only at the
stated level (mod-p scout vs. exact Q(q,t)).

## theta3_min_poly.py — the derived three-strand untwisting

**Question (finitization program).** Impose only the two-strand untwisting `p_2(sigma_1) X_2 = 0`
plus all zipper relations — no three-strand relation. Is `Z_3` finite-dimensional, and what
minimal polynomial does the lasso twist `Theta_3 = sigma_1 sigma_2^2 sigma_1` satisfy on the
slider `X_3`? That min poly *is* the derived (least-imposed) three-strand untwisting.

**Result.**
- Finite-dimensional for all three `p_2` (dim stable across a (q,t) grid, so generic):
  linear `sigma-t` -> 19; `(sigma-t)^2` -> 32; `sigma^2-t` -> 32.
- The min poly of `Theta_3` on `X_3` is **exactly `p*(lambda)`** (the reflection of `p_2`),
  not the module ceiling `(lambda-q^2) p*`:

  | p_2 | p* | min poly of Theta_3 on X_3 |
  |-----|----|----------------------------|
  | sigma - t    | lambda + t     | lambda + t |
  | (sigma - t)^2 | (lambda + t)^2 | (lambda + t)^2 |
  | sigma^2 - t  | lambda^2 - t   | lambda^2 - t |

  `X_3 != 0` in all three. The algebra drops the rigid `(lambda - q^2)` factor that the
  *module* retains (`thm:modelp` proves only that every module relation is divisible by
  `(lambda - q^2) p*`). This confirms the conjectural sharpening (`conj:recursion`) and
  explains Conjecture 1: imposing the canonical `Theta_3 = q^2` kills `X_3` exactly when
  `q^2` is not a root of `p*`, i.e. `p(-q^2) != 0`.

**Reproduce.**
- Mod-p scout / grid (fast, speculative):
  `python3 -m f4ncgb_tools.research.theta3_min_poly --case all --scratch /tmp`
- Exact over Q(q,t) (load-bearing): `gap gap_scripts/theta3_exact.gap`
  (linear confirmed: dim 19, X_3 != 0, `(Theta_3 + t) X_3 = 0`; jordan/framing dim 32).

**Trust.** The grid runner is mod-p at 2^31-1 — a scout. The exact statement is the GAP
FunctionField(Rationals,[q,t]) computation. A self-contained exact path (char-0 grid +
rational reconstruction with a degree bound) is a TODO once `corelib/reconstruct.py` lands.
