# Historical certificate diagnostics

These scripts inspect earlier compressed certificates. Their inputs and
assumptions are specific to those artifacts; they are outside the maintained
installed package.

| Script | Historical purpose |
| --- | --- |
| `cancel_qm1.py` | Cancel a common numerator/denominator factor in an earlier certificate. |
| `check_denom_minimal.py` | Test the displayed common denominator for that fixed expression. |
| `independent_verify_three_rels.py` | Reparse and check a compressed three-relation expression by an independent route. |

Denominator minimality for a fixed expression is not minimality among every
possible certificate. Compressed braid-algebra expressions and literal
free-algebra certificates use different equality checks. The current package
interface and validation instructions are in [WORKFLOW](../../WORKFLOW.md) and
[TESTING](../../docs/TESTING.md).
