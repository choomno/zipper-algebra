"""Decode numeric f4ncgb proof data for independent completion verification.

This parser returns a candidate certificate. Only ``verify_groebner`` establishes
ideal equality and completeness. Rational output coefficients are interpreted
exactly in the declared field; no parameter-specialization conclusion is inferred.
"""
from __future__ import annotations

from typing import cast

from ..algebra import AlgebraPresentation
from ..coefficients import Coef
from ..cofactor.proof import Basis, Proof
from ..execution.artifacts import _check_syntax
from .groebner import GroebnerCertificate, ProofTerm

__all__ = ["certificate_from_solver"]


def certificate_from_solver(presentation: AlgebraPresentation, basis_text: str,
                            proof_text: str) -> GroebnerCertificate:
    """Read numeric basis/proof text as a domain-bound derivation graph.

    Shared and forward references remain unexpanded. The syntax reader checks
    complete terms before the historical tokenizer runs. The returned data must
    be verified against the supplied literal presentation before using normal
    forms or reporting dimension.
    """
    if not isinstance(presentation, AlgebraPresentation):
        raise TypeError("solver certificate decoding requires an AlgebraPresentation")
    if not isinstance(basis_text, str) or not isinstance(proof_text, str):
        raise TypeError("basis and proof data must be decoded text")
    names = list(presentation.algebra.generators)
    _check_syntax(basis_text, names, proof=False)
    _check_syntax(proof_text, names, proof=True)
    basis = Basis.from_text(basis_text, 0, names)
    proof = Proof.from_text(proof_text, names)
    algebra, field = presentation.algebra, presentation.algebra.field
    # Basis.from_text constructs numeric Coef values; the wider Poly annotation
    # also serves the separate symbolic q,t relation builder.
    polynomials = tuple(algebra.element({word: field.scalar(cast(Coef, coefficient).v)
                                         for word, coefficient in polynomial.terms.items()})
                        for polynomial in basis.polys)
    derivations = tuple(tuple(ProofTerm(field.scalar(coefficient), left, kind, index, right)
                              for coefficient, left, kind, index, right in row)
                        for row in proof.terms)
    return GroebnerCertificate(presentation.content_hash(), polynomials, derivations)
