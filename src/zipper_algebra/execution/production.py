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
from ..action_production import word_basis_from_groebner
from ..cofactor.artifacts import write_json_atomic
from ..production import reconstruct_groebner
from ..provenance import capture_provenance, write_bundle_provenance
from ..verification.finite import WordBasisCertificate
from ..verification.groebner import GroebnerCertificate
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
        generic = GroebnerCertificate.from_json(loaded["groebner-certificate.json"])
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
                    max_operations=10000000, max_products=10000000,
                    max_reductions=2000000, max_words=100000,
                    max_primes=5, max_degree=24, threads=1, start_prime=0,
                    shift=None) -> dict:
    """Run f4ncgb, reconstruct a generic basis and certify finite generator actions.

    Results and numerical solver data stay in the configured ignored run folder.
    A verified result requires both exact generic completion and an independently
    checked finite word-basis certificate. Failures are recorded as inconclusive.
    Numerical sampling never establishes a generic dimension. Solver timeout is
    process-level; algebraic reconstruction limits are operation/search bounds.
    """
    if not isinstance(presentation, AlgebraPresentation) or not isinstance(config, RunConfig):
        raise TypeError("production requires an AlgebraPresentation and RunConfig")
    root = config.output_root
    prepare_output_root(root)
    started = time.monotonic()
    stage = "solver"
    report = {"schema": "zipper-algebra-production-v1",
              "presentation_hash": presentation.content_hash(), "status": "running",
              "verified": False, "dimension": None, "basis_complete": False}
    # Invalidate a previous selected result before any expensive work. A killed
    # process then leaves an unfinished attempt, never a stale verified verdict.
    write_json_atomic(root / "result.json", report)
    producer = capture_provenance()
    report["producer"] = producer
    try:
        sample = run_algebra_point(presentation, point, config, refresh=refresh, producer=producer,
                                   max_operations=max_operations, max_nodes=max_nodes)
        report["solver"] = {"key": sample.key, "store": "samples.sqlite3",
                            "generation": sample.record["generation"], "hit": sample.hit,
                            "seconds": sample.seconds, "request": sample.record["request"]}
        stage = "reconstruction"
        dependencies = reconstruction_provenance()
        options = {"max_nodes": max_nodes, "max_operations": max_operations,
                   "max_products": max_products, "max_reductions": max_reductions,
                   "max_words": max_words, "max_primes": max_primes,
                   "max_degree": max_degree, "threads": threads, "start_prime": start_prime,
                   "shift": None if shift is None else list(shift)}
        reconstructed = reconstruct_groebner(presentation, solver_certificate(sample, presentation),
                                             point=point, **options)
        stage = "actions"
        word_basis = word_basis_from_groebner(presentation, reconstructed.certificate,
                         max_nodes=max_nodes, max_operations=max_operations,
                         max_reductions=max_reductions, max_words=max_words)
        stage = "publication"
        if reconstruction_provenance() != dependencies:
            raise ArtifactError("reconstruction dependencies changed during production")
        validate_source_snapshot(producer)
        metadata = {"schema": "zipper-algebra-production-evidence-v1",
                    "solver": report["solver"], "reconstruction": {"dependencies": dependencies,
                    "options": options, "counts": reconstructed.statistics},
                    "completion": reconstructed.reduction.report,
                    "dimension": len(word_basis.basis_words)}
        products = _publish(root, {"presentation.json": presentation.to_json(),
                    "groebner-certificate.json": reconstructed.certificate.to_json(),
                    "word-basis-certificate.json": word_basis.to_json(),
                    "production.json": metadata}, producer)
        report.update(status="verified", verified=True, basis_complete=True,
                      dimension=len(word_basis.basis_words), **products)
    except (ExecutionError, ValueError, TypeError, ArithmeticError, RuntimeError, ImportError, OSError) as error:
        report.update(status="inconclusive", reason=str(error), stage=stage)
    report["seconds"] = time.monotonic() - started
    write_json_atomic(root / "result.json", report)
    return report
