# Compatibility tools and experiments

The maintained implementation lives in `src/zipper_algebra/`. Install the package
before using the compatibility entry points in this directory:

```sh
python -m pip install -e .
```

Prefer the public `zipper_algebra` imports and `zipper-algebra` commands for new
work. Corresponding modules here forward to the same maintained implementation,
so they do not create independent copies of coefficient classes or arithmetic.

| Location | Purpose |
| --- | --- |
| `presentation.py`, `zipper_in_ms.py`, `core.py` | Compatibility imports and historical entry points. |
| `cofactor/` | Compatibility modules, regression comparisons and experimental solver orchestration. |
| `variants/`, `research/`, `survey/` | Experimental dimension, representation and relation searches. |
| `auxiliary/` | Diagnostics for older certificates. |
| `tests.py` | Historical integration suite requiring external solver configuration. |

The package's algebraic and certificate operations need no solver. Experimental
execution scripts have additional dependencies and may still use the former
workspace settings. Their results need the stated presentation, coefficient
domain and completion evidence; saved outputs are not automatically generic
proofs.

Use [TESTING](../docs/TESTING.md) for the default CI tests and
[WORKFLOW](../WORKFLOW.md) for the current certificate workflow.
