"""Fresh generic finite-algebra production with retained, byte-bound evidence."""
from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
import json
import os
from pathlib import Path
import tempfile
import time

from ..algebra import AlgebraPresentation
from ..action_production import word_basis_from_rewrites
from ..cofactor.artifacts import write_json_atomic
from ..production import reconstruct_basis
from ..provenance import capture_provenance, write_bundle_provenance
from ..verification.finite import WordBasisCertificate
from ..verification.groebner import GroebnerCertificate, verify_groebner
from ..verification.derivations import VerificationLimit
from .artifacts import (ArtifactError, canonical_json, digest_bytes, file_identity,
                        prepare_output_root, validate_source_snapshot)
from .config import RunConfig
from .dependencies import reconstruction_provenance
from .general import run_algebra_point, solver_certificate
from .pipeline import _generation_matches
from .runner import ExecutionError

__all__ = ["produce_algebra"]


def _publish(root, documents, producer):
    destination = root / "algebras"
    destination.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".publish-", dir=destination) as temporary:
        stage = Path(temporary) / "generation"
        stage.mkdir()
        for name, document in documents.items():
            write_json_atomic(stage / name, document)
        loaded = {name: json.loads((stage / name).read_text()) for name in documents}
        if canonical_json(loaded) != canonical_json(documents):
            raise ArtifactError("serialized algebra production differs from verified data")
        # Exercise the public decoders on the exact saved certificate encodings.
        presentation = AlgebraPresentation.from_json(loaded["presentation.json"])
        generic = GroebnerCertificate.from_json(loaded["rewrite-candidate.json"])
        finite = WordBasisCertificate.from_json(loaded["word-basis-certificate.json"])
        if {generic.presentation_hash, finite.presentation_hash} != {presentation.content_hash()}:
            raise ArtifactError("serialized certificates differ from the saved presentation")
        write_bundle_provenance(stage, documents, producer=producer)
        validate_source_snapshot(producer)
        files = {path.name: file_identity(path) for path in stage.iterdir()}
        generation = digest_bytes(canonical_json(files))
        final = destination / generation
        if final.exists():
            if not _generation_matches(final, files):
                raise ArtifactError("existing algebra generation has different content")
        else:
            try:
                os.replace(stage, final)
            except OSError:
                if not _generation_matches(final, files):
                    raise
    return {"generation": generation, "files": files,
            "artifacts": {name: str(final / name) for name in files}}


def produce_algebra(presentation: AlgebraPresentation, point: Mapping[str, int | Fraction],
                    config: RunConfig, *, refresh=False, max_nodes=100000,
                    max_operations=10000000,
                    max_reductions=2000000, max_words=100000,
                    max_primes=5, max_degree=24, threads=1, start_prime=0,
                    shift=None, verify_completion: bool = False) -> dict:
    """Run f4ncgb, reconstruct a generic basis and certify finite generator actions.

    Results and numerical solver data stay in the configured ignored run folder.
    A verified result requires an independently checked finite word-basis
    certificate. ``verify_completion`` additionally requests complete Gröbner
    reduction; it is not required to certify this finite algebra.
    Numerical sampling never establishes a generic dimension. Solver timeout is
    process-level; algebraic reconstruction limits are operation/search bounds.
    """
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(config, RunConfig):
        raise TypeError("production requires an AlgebraPresentation and RunConfig")
    if type(verify_completion) is not bool:
        raise ValueError("verify_completion must be Boolean")
    root = config.output_root
    prepare_output_root(root)
    started = time.monotonic()
    stage = "solver"
    report = {"schema": "zipper-algebra-production-v1",
              "presentation_hash": presentation.content_hash(), "status": "running",
              "verified": False, "dimension": None, "word_basis_verified": False,
              "groebner_complete": None}
    timings = {}
    # Invalidate a previous selected result before any expensive work. A killed
    # process then leaves an unfinished attempt, never a stale verified verdict.
    write_json_atomic(root / "result.json", report)
    producer = capture_provenance()
    report["producer"] = producer
    try:
        before = time.monotonic()
        sample = run_algebra_point(presentation, point, config, refresh=refresh, producer=producer,
                                   max_operations=max_operations, max_nodes=max_nodes)
        timings["acquisition"] = time.monotonic() - before
        report["solver"] = {"key": sample.key, "store": "samples.sqlite3",
                            "generation": sample.record["generation"], "hit": sample.hit,
                            "seconds": sample.seconds, "request": sample.record["request"]}
        stage = "reconstruction"
        dependencies = reconstruction_provenance()
        options = {"max_nodes": max_nodes, "max_operations": max_operations,
                   "max_primes": max_primes,
                   "max_degree": max_degree, "threads": threads, "start_prime": start_prime,
                   "shift": None if shift is None else list(shift)}
        before = time.monotonic()
        reconstructed = reconstruct_basis(presentation, solver_certificate(sample, presentation),
                                           point=point, **options)
        timings["reconstruction"] = time.monotonic() - before
        stage = "actions"
        before = time.monotonic()
        word_basis = word_basis_from_rewrites(presentation, reconstructed.certificate,
                         max_nodes=max_nodes, max_operations=max_operations,
                         max_words=max_words)
        timings["actions_and_certification"] = time.monotonic() - before
        completion = None
        if verify_completion:
            stage = "completion"
            before = time.monotonic()
            try:
                completion = verify_groebner(presentation, reconstructed.certificate,
                               max_reductions=max_reductions, max_words=max_words, prune=True).report
            except VerificationLimit as error:
                # The finite algebra is already certified. An exhausted
                # optional cross-check does not revoke that mathematical result.
                completion = {"status": "inconclusive", "reason": str(error)}
            timings["completion"] = time.monotonic() - before
        stage = "publication"
        if reconstruction_provenance() != dependencies:
            raise ArtifactError("reconstruction dependencies changed during production")
        validate_source_snapshot(producer)
        metadata = {"schema": "zipper-algebra-production-evidence-v1",
                    "solver": report["solver"], "reconstruction": {"dependencies": dependencies,
                    "options": options, "counts": reconstructed.statistics},
                    "completion": completion,
                    "verification_limits": {"max_operations": max_operations,
                                            "max_nodes": max_nodes, "max_dimension": max_words},
                    "completion_options": {"requested": verify_completion,
                                           "max_reductions": max_reductions, "max_words": max_words},
                    "timings": dict(timings),
                    "dimension": len(word_basis.basis_words)}
        before = time.monotonic()
        products = _publish(root, {"presentation.json": presentation.to_json(),
                    "rewrite-candidate.json": reconstructed.certificate.to_json(),
                    "word-basis-certificate.json": word_basis.to_json(),
                    "production.json": metadata}, producer)
        timings["publication"] = time.monotonic() - before
        report.update(status="verified", verified=True, word_basis_verified=True,
                      groebner_complete=True if completion and completion.get("basis_complete") else None,
                      dimension=len(word_basis.basis_words), **products)
    except (ExecutionError, ValueError, TypeError, ArithmeticError, RuntimeError, ImportError, OSError) as error:
        report.update(status="inconclusive", reason=str(error), stage=stage)
    report["seconds"] = time.monotonic() - started
    report["timings"] = timings
    write_json_atomic(root / "result.json", report)
    return report
