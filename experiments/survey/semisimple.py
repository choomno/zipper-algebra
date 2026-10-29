"""Exact trace-matrix observations at specified parameter points.

The solver output used here carries no basis-completeness proof. Consequently
these observations do not certify dimensions, radical dimensions, generic ranks
or semisimplicity of the presented quotient. In positive characteristic a
regular-trace form may be degenerate even for a semisimple algebra.
"""
from __future__ import annotations

from zipper_algebra.linear import matrix_rank, validate_field_characteristic
from .pointalg import PointAlgebra
from .spec import Spec

DEFAULT_CHAR = 2147483647
TRACE_SCHEMA = "zipper-survey-trace-observations-v1"


def trace_observations(spec: Spec, scratch: str, points=((5, 7), (13, 17)),
                       char: int = DEFAULT_CHAR):
    """Collect specialized trace observations; no structural verdict is inferred.

    Every point retains its own standard-word count and matrix rank; ranks
    from different points
    are never combined with an unrelated dimension.
    """
    validate_field_characteristic(char)
    observations = []
    for q, t in points:
        model = PointAlgebra(spec, q, t, scratch, char=char)
        gram = model.trace_form()
        observations.append({
            "q": q, "t": t, "standard_word_count": model.dim,
            "trace_rank": matrix_rank(gram, characteristic=char),
            "trace_form_symmetric": all(gram[i][j] == gram[j][i]
                                        for i in range(model.dim) for j in range(model.dim)),
            "basis_complete": False,
        })
    if not observations:
        raise ValueError("at least one parameter point is required")
    return {
        "schema": TRACE_SCHEMA, "spec": spec.name,
        "status": "observed", "basis_complete": False,
        "coefficient_domain": "Q" if char == 0 else f"GF({char})", "char": char,
        "standard_word_counts_agree": len({p["standard_word_count"] for p in observations}) == 1,
        "points": observations,
    }
