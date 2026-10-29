#!/usr/bin/env python3
"""Compatibility imports and command for the separated presentation emitter.

Coefficient arithmetic, free polynomials, relation construction and msolve
emission live in ``coefficients``, ``polynomial``, ``relations`` and ``msolve``.
This module retains the established imports and ``python -m`` command.
"""
from zipper_algebra.coefficients import Coef, QTCoef, NumericRing, SymbolicRing
from zipper_algebra.polynomial import Poly
from zipper_algebra.relations import (
    RelTag, build_relations, build_relations_ring, append_zero_sliders, append_target,
)
from zipper_algebra.msolve import main

__all__ = [
    "Coef", "QTCoef", "NumericRing", "SymbolicRing", "Poly", "RelTag",
    "build_relations", "build_relations_ring", "append_zero_sliders", "append_target",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
