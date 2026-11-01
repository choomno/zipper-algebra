"""Acquire shared numerical proofs for arbitrary literal parameter presentations.

Numerical acquisition supplies candidate supports for generic reconstruction.
Published samples prove only input-ideal identities at their recorded rational
point; neither completion nor a generic conclusion follows from acquisition.
"""
from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction

from ..algebra import AlgebraPresentation
from ..provenance import capture_provenance, validate_provenance
from ..verification.derivations import DerivationGraph, DerivationTerm, verify_derivations
from ..verification.f4ncgb import certificate_from_solver
from ..verification.groebner import GroebnerCertificate, ProofTerm
from .artifacts import ArtifactError, digest_bytes, file_identity, validate_source_snapshot
from .config import RunConfig
from .runner import SolverRun, _run_request

__all__ = ["run_algebra_point", "solver_certificate"]


def _specialize(presentation, values):
    if not isinstance(presentation, AlgebraPresentation):
        raise TypeError("acquisition requires an AlgebraPresentation")
    if presentation.algebra.field.characteristic:
        raise ValueError("the acquisition runner supports characteristic zero only")
    if (not isinstance(values, Mapping)
            or set(values) != set(presentation.algebra.field.parameters)
            or any(type(value) not in (int, Fraction) for value in values.values())):
        raise ValueError("the point must give an exact rational value for every declared parameter")
    point = [[name, str(Fraction(values[name]))] for name in presentation.algebra.field.parameters]
    return point, presentation.specialize(values)


def _input(presentation):
    algebra = presentation.algebra
    if not algebra.generators:
        raise ValueError("numerical solver acquisition requires at least one generator")
    indices = [index for index, relation in enumerate(presentation.relations) if not relation.is_zero()]
    if not indices:
        raise ValueError("no nonzero specialized relations: no numerical proof acquisition is needed")
    lines = []
    for index in indices:
        terms = []
        for word, coefficient in presentation.relations[index].terms:
            # Specialization gives constants. Division retains exact rational
            # input even if the stored field expression is not normalized.
            scalar = Fraction(coefficient.numerator[0][1]) / Fraction(coefficient.denominator[0][1])
            terms.append(str(scalar) + ("*" + "*".join(algebra.generators[i - 1] for i in word)
                                        if word else ""))
        lines.append("+".join(terms).replace("+-", "-"))
    text = ",".join(algebra.generators) + "\n0\n" + ",\n".join(lines) + "\n"
    return indices, text.encode()


def _decode(presentation, files, indices, limits):
    """Map native compact relation indices back to the literal ordered inputs."""
    try:
        parsed = certificate_from_solver(presentation, files["basis.txt"].decode("utf-8"),
                                         files["proof.txt"].decode("utf-8"))
        if not parsed.basis:
            raise ValueError("nonempty input requires a nonempty candidate basis")
        derivations = []
        for row in parsed.derivations:
            mapped = []
            for term in row:
                if term.kind == "input":
                    if term.index >= len(indices):
                        raise ValueError("proof input index exceeds the emitted relation list")
                    term = ProofTerm(term.coefficient, term.left, term.kind, indices[term.index], term.right)
                mapped.append(term)
            derivations.append(tuple(mapped))
        result = GroebnerCertificate(presentation.content_hash(), parsed.basis, tuple(derivations))
        graph = DerivationGraph(result.presentation_hash, result.basis,
                                tuple(tuple(DerivationTerm(term.coefficient, term.left,
                                                          "node" if term.kind == "basis" else "input",
                                                          term.index, term.right) for term in row)
                                      for row in result.derivations))
        verify_derivations(presentation, graph, **limits)
        return result
    except (KeyError, ValueError, TypeError, ZeroDivisionError) as error:
        raise ArtifactError(f"invalid algebra solver proof: {error}") from error


def run_algebra_point(presentation: AlgebraPresentation, values: Mapping[str, int | Fraction],
                      config: RunConfig, *, refresh: bool = False, producer=None,
                      max_operations: int = 2_000_000, max_nodes: int = 100_000) -> SolverRun:
    """Run or reuse f4ncgb at a complete exact rational parameter point.

    Zero specialized input rows are omitted because the native solver cannot
    reliably handle them. Their original positions are recorded explicitly and
    proof references are remapped before independent bounded validation. Raw
    native artifacts are retained unchanged. Every cache hit is reverified.
    """
    for limit in (max_operations, max_nodes):
        if type(limit) is not int or limit < 1:
            raise ValueError("proof validation limits must be positive integers")
    point, specialized = _specialize(presentation, values)
    indices, input_bytes = _input(specialized)
    producer = capture_provenance() if producer is None else validate_provenance(producer)
    validate_source_snapshot(producer)
    limits = {"max_operations": max_operations, "max_nodes": max_nodes}
    request = {
        "schema": "zipper-solver-request-v4",
        "presentation": {"sha256": presentation.content_hash(), "field": presentation.algebra.field.to_json()},
        "point": point, "specialized_domain": "Q",
        "specialized_presentation_hash": specialized.content_hash(),
        "input_relation_indices": indices, "proof_limits": limits,
        "monomial_order": config.monomial_order, "input_sha256": digest_bytes(input_bytes),
        "producer": producer,
        "solver": {"binary": str(config.binary), "identity": file_identity(config.binary),
                   "revision": config.solver_revision, "io_mode": config.solver_io,
                   "prefix": [{"path": str(path), **file_identity(path)} for path in config.command_prefix],
                   "flags": list(config.arguments()), "timeout_seconds": config.timeout_seconds},
    }
    def validate(files):
        return len(_decode(specialized, files, indices, limits).basis)
    return _run_request(config, input_bytes, request, validate, refresh=refresh)


def solver_certificate(run: SolverRun, presentation: AlgebraPresentation) -> GroebnerCertificate:
    """Recheck a general acquisition and return its remapped numerical proof.

    The result is bound to the full specialized input presentation, including
    zero relations. It is a candidate for a completion check, or a source of
    proof supports for generic reconstruction; it makes neither assertion.
    """
    request = run.record["request"]
    if (request.get("schema") != "zipper-solver-request-v4"
            or request["presentation"] != {"sha256": presentation.content_hash(),
                                           "field": presentation.algebra.field.to_json()}):
        raise ArtifactError("sample does not belong to this literal algebra presentation")
    try:
        values = {name: Fraction(value) for name, value in request["point"]}
        point, specialized = _specialize(presentation, values)
        indices, input_bytes = _input(specialized)
        if (point != request["point"] or indices != request["input_relation_indices"]
                or specialized.content_hash() != request["specialized_presentation_hash"]
                or digest_bytes(input_bytes) != request["input_sha256"]
                or run.files["input.ms"] != input_bytes):
            raise ValueError("sample's specialized input differs from its declared presentation")
        return _decode(specialized, run.files, indices, request["proof_limits"])
    except (KeyError, ValueError, TypeError, ZeroDivisionError) as error:
        raise ArtifactError(f"invalid algebra sample: {error}") from error
