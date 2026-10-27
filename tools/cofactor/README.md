# Cofactor compatibility and execution tools

Maintained certificate code lives in `src/zipper_algebra/cofactor/`. The matching
modules in this directory are compatibility entry points. New code should import
from `zipper_algebra.cofactor` and use the installed `zipper-algebra verify`
command. See [WORKFLOW](../../WORKFLOW.md) for a complete example.

## Maintained algebraic stages

| Module | Purpose |
| --- | --- |
| `atoms` | Structured denominator factors and optional inference. |
| `cofactor_io` | Validated certificate data, targets, readers and writers. |
| `proof` | Basis keys and expansion of acyclic proof references. |
| `skeleton` | Recorded reductions with full word contexts. |
| `lift` | Exact interpolation and candidate coefficient reconstruction. |
| `compose` | Context-preserving substitution and occurrence-based multipliers. |
| `verify`, `ncpoly` | Exact symbolic and degree-bounded grid checks in the free algebra. |

The [mathematical specification](../../docs/mathematical-specification/main.tex)
defines their input domains and conclusions. A candidate reconstruction is
accepted as a membership result only after exact identity verification.

## Experimental and historical tools

`pipeline.py`, `samples.py` and `paths.py` retain the solver sampling workflow and
its workspace configuration. Specialized compression and normalization scripts
operate on historical certificate formats, including braid-group-algebra
expressions. They are separate from the installed verification command and need
individual dependency and completion checks before use.

`tests_*.py` includes regressions and comparisons with archived implementations.
The default suite and immutable fixture locations are described in
[TESTING](../../docs/TESTING.md). Historical free-algebra certificates remain
regression inputs rather than the current statement of the mathematical proof.
