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

## theta_min_poly.py — general n (supersedes theta3_min_poly.py)

Same question at arbitrary `n`: min poly of `Theta_k = x1..x_{k-1} x_{k-1}..x1` on
`X_k = S[(k,1)]` in the two-strand-only quotient `Z_n`. Per-point f4ncgb basis (mod p, or
`--char 0` for exact Q), then an incremental Krylov min poly with a **memoized** normal-form
reducer — the stock `core.reduce_to_basis_*` has no memoization and blows up past n=3
(a single `Theta_4 * <len-6 monomial>` reduction takes ~40 s; memoized + short-word
incremental it is milliseconds). Over Q it also checks the min poly against the tower formula.

**n=4 result (the interesting one; n=3 collapses to p*, n=4 does not).** Two-strand-only `Z_4`
is finite (dim 214 / 404 / 404). `Theta_3` on `X_3` is unchanged (= p*), but `Theta_4` on `X_4`
is strictly larger:

  g_4(L) = **p*(L) · q^{2d} p*(L / q^2)**,  d = deg p*  (spectrum = roots(p*) ∪ q^2·roots(p*))

  | p_2 | min poly of Theta_4 on X_4 |
  |-----|----------------------------|
  | sigma - t    | (L + t)(L + q^2 t) |
  | (sigma - t)^2 | (L + t)^2 (L + q^2 t)^2 |
  | sigma^2 - t  | (L^2 - t)(L^2 - q^4 t) |

So the sharpening-to-p* is special to level 3; at level 4 one q^2-dilated copy of p* survives
(first concrete appearance of `conj:recursion`'s q-power side). **Tower conjecture:**
`g_k(L) = prod_{j=0}^{k-3} q^{2 j d} p*(L / q^{2 j})`. Verified: mod-p at 3 points, exact Q at 2
points (all `== tower? YES`), for n=4. n=5 test pending (GB ~1648 elements — large).

**Reproduce.**
- Mod-p, n=4:  `python3 -m f4ncgb_tools.research.theta_min_poly --n 4 --scratch /tmp`
- Exact Q + tower check:  `... --n 4 --char 0 --grid 5,7 3,10 --scratch /tmp`

See the repo `NIGHT_LOG.md` "n=4" section for the full write-up and verification status.
