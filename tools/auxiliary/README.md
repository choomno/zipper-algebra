# reduction — post-hoc work on a finished cofactor

**Not part of the pipeline.** `cofactor/` produces and verifies the ℚ(q,t) cofactor of X₃; these operate
on that finished result. They are outside the transitive closure of `reproduce_X3_zero.sh`, so removing
them cannot affect reproducibility.

| script | what it does |
|---|---|
| `cancel_qm1.py` | Cancels the spurious (q−1) factor from the 161-tuple three-relation cofactor — every numerator carries it. |
| `check_denom_minimal.py` | Exhaustive proof that `t⁴(t+q²)` is the minimal common denominator of the 161-tuple cofactor. |
| `independent_verify_three_rels.py` | Re-parses `cofactor-X3-n4-Qqt-three-rels.txt` from scratch and verifies it, deliberately **not** using the symbolic engine that produced it. |

They import from `../cofactor/` (`relations_naming`, `symbolic_identity_check`, `simplify_cofactor`).
