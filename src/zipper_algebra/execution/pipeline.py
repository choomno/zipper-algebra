"""Produce and exactly verify a cofactor from a presentation and solver proofs."""
from __future__ import annotations

from pathlib import Path
import json
import os
import tempfile
import time

from ..cofactor.artifacts import (presentation_identity, load_artifact, save_artifact, target_identity,
                                 write_json_atomic)
from ..cofactor.cofactor_io import Cofactor, resolve_target
from ..cofactor.discovery import build_system, reconstruct_cofactor
from ..cofactor.skeleton import Anchor, Skeleton, reduce
from ..cofactor.verify import check
from ..coefficients import NumericRing
from ..presentation import Presentation
from ..matrices import ExactMatrix
from ..verification.linear import KernelCertificate
from .artifacts import ArtifactError, canonical_json, digest_bytes, file_identity, prepare_output_root
from .runner import ExecutionError, validate_source_snapshot
from .sampling import SampleSet
from .dependencies import reconstruction_provenance
from .store import read as read_solver_run
from .config import RunConfig
from ..provenance import capture_provenance, write_bundle_provenance


REDUCTION_ORDER = "length-lex-index-v1"


def _skeleton(pres: Presentation, target: str, samples: SampleSet) -> Skeleton:
    basis, proof = samples.anchor_basis(), samples.anchor_proof()
    polynomial = resolve_target(pres, target, NumericRing(*samples.anchor))
    steps, remainder = reduce(polynomial, basis)
    return Skeleton.from_reduction(pres, target, samples.anchor, steps, remainder, basis, proof)



def _solver_runs(samples: SampleSet) -> list[dict]:
    """Identify all validated runs of this discovery, including unstable samples."""
    return [{"point": list(point), "key": result.key, "store": str(result.store),
             "generation": result.record["generation"], "files": result.record["files"],
             "request": result.record["request"]}
            for point, result in samples._runs.items()]


def _solver_links(root, pres, target, solver_runs, producer):
    """Bind supplied run records to retained output-root-relative solver data.

    Reading a stored generation checks its identities and bytes; it does not
    establish basis completeness or completeness of the chosen context family.
    """
    presentation = presentation_identity(pres)
    records: list[dict] = []
    links = {"schema": "zipper-solver-links-v1", "presentation": presentation,
             "target": target_identity(pres, target), "recorded": solver_runs is not None,
             "store_path_base": "output-root", "runs": records}
    if solver_runs is None:
        return links
    seen = set()
    try:
        for run in solver_runs:
            point = run["point"]
            if (not isinstance(point, list) or len(point) != 2
                    or any(type(value) is not int for value in point)):
                raise ArtifactError("solver link point must be an integer q,t pair")
            if tuple(point) in seen:
                raise ArtifactError("solver links must contain one generation per point")
            seen.add(tuple(point))
            path = Path(run["store"])
            relative = path.resolve().relative_to(Path(root).resolve()).as_posix()
            _, saved = read_solver_run(path, run["key"], generation=run["generation"])
            request = saved["request"]
            if (canonical_json(request) != canonical_json(run["request"])
                    or canonical_json(request["presentation"]) != canonical_json(presentation)
                    or canonical_json(request["point"]) != canonical_json([*point, 0])
                    or request["specialized_domain"] != "Q"
                    or canonical_json(request["producer"]) != canonical_json(producer)
                    or canonical_json(saved["files"]) != canonical_json(run["files"])):
                raise ArtifactError("solver link metadata differ from the retained computation")
            records.append({"point": list(point), "key": run["key"], "store": relative,
                            "generation": run["generation"], "files": saved["files"]})
    except (KeyError, TypeError, ValueError, OSError) as error:
        raise ArtifactError(f"invalid solver links: {error}") from error
    records.sort(key=lambda run: run["point"])
    return links


def _generation_matches(directory, files):
    """An immutable generation contains exactly its declared file names and bytes."""
    return (directory.is_dir() and {path.name for path in directory.iterdir()} == set(files)
            and all(file_identity(directory / name) == identity for name, identity in files.items()))


def _publish_verified(root, pres, target, skeleton, discovery, order,
                      producer=None, *, solver_runs=None, reconstruction=None):
    """Publish one immutable product generation; the result record selects it later."""
    destination = root / "verified"
    producer = capture_provenance() if producer is None else producer
    solver_links = _solver_links(root, pres, target, solver_runs, producer)
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".publish-", dir=destination) as temporary:
        stage = Path(temporary) / "generation"
        stage.mkdir()
        save_artifact(stage / "skeleton.json", skeleton, pres, target, order=order, producer=producer)
        documents = {"support.json": discovery.system.to_json(),
                     "matrix.json": discovery.system.matrix.to_json(),
                     "kernel.json": discovery.kernel.to_json(),
                     "reconstruction.json": reconstruction, "solver-runs.json": solver_links}
        for name, document in documents.items():
            write_json_atomic(stage / name, document)
        save_artifact(stage / "certificate.json", discovery.cofactor, pres, target, producer=producer)
        discovery.cofactor.save(stage / "certificate.txt", provenance=False)
        # Reparse published encodings and bind them to the objects whose exact
        # identities were checked. Equality avoids repeating large AN products.
        try:
            if any(json.loads((stage / name).read_text()) != document
                   for name, document in documents.items()):
                raise ValueError("saved reconstruction metadata differ from checked objects")
            bound_skeleton = load_artifact(stage / "skeleton.json", pres, target,
                                           kind="skeleton", order=order)
            bound_certificate = load_artifact(stage / "certificate.json", pres, target, kind="cofactor")
            matrix = ExactMatrix.from_json(json.loads((stage / "matrix.json").read_text()))
            kernel = KernelCertificate.from_json(json.loads((stage / "kernel.json").read_text()))
            certificate = Cofactor.load(str(stage / "certificate.txt"), pres=pres)
            if (bound_skeleton.to_json() != skeleton.to_json()
                    or bound_certificate.header != discovery.cofactor.header
                    or bound_certificate.terms != discovery.cofactor.terms
                    or matrix.to_json() != discovery.system.matrix.to_json()
                    or kernel.to_json() != discovery.kernel.to_json()
                    or certificate.header != discovery.cofactor.header
                    or certificate.terms != discovery.cofactor.terms):
                raise ValueError("saved reconstruction objects differ from checked objects")
        except (ValueError, TypeError, KeyError) as error:
            raise ArtifactError(f"invalid serialized reconstruction: {error}") from error
        write_bundle_provenance(stage, (path.name for path in stage.iterdir()), producer=producer)
        validate_source_snapshot(producer)
        files = {path.name: file_identity(path) for path in stage.iterdir()}
        generation = digest_bytes(canonical_json(files))
        final = destination / generation
        if final.exists():
            if not _generation_matches(final, files):
                raise ArtifactError("existing verified generation has different content")
        else:
            try:
                os.replace(stage, final)
            except OSError as error:
                # Another producer may have published these same immutable
                # bytes after the existence check. Only that complete winner
                # permits us to treat a failed rename as successful reuse.
                try:
                    winner_matches = _generation_matches(final, files)
                except OSError:
                    raise error
                if not winner_matches:
                    raise
    return {"generation": generation, "files": files,
            "artifacts": {name: str(final / name) for name in files},
            "certificate": str(final / "certificate.txt")}


def run(pres: Presentation, target: str, anchor: Anchor, config: RunConfig, *,
        threads: int = 1, max_primes: int = 5, reconstruction_degree: int = 24,
        max_contexts: int = 100000, max_products: int = 10000000,
        start_prime: int = 0, shift=None, refresh: bool = False,
        verbose: bool = True) -> dict:
    """Discover an exact generic cofactor using one fresh solver proof graph.

    FiniteFlow solves the contextual word equations over Q(q,t). The reconstructed
    kernel and the resulting free-algebra identity are checked independently.
    An insufficient proof support or exhausted search establishes no nonmembership.
    """
    for name, value, minimum in (("threads", threads, 1), ("max_primes", max_primes, 1),
                                 ("reconstruction_degree", reconstruction_degree, 1),
                                 ("max_contexts", max_contexts, 1),
                                 ("max_products", max_products, 1),
                                 ("start_prime", start_prime, 0)):
        if type(value) is not int or value < minimum:
            raise ValueError(f"{name} must be an integer at least {minimum}")
    if (max_primes > 2047 or start_prime + max_primes >= 2048
            or reconstruction_degree >= 2**16 or threads >= 2**32):
        raise ValueError("reconstruction options exceed FiniteFlow native limits")
    if shift is not None:
        shift = tuple(shift)
        if len(shift) != 2 or any(type(x) is not int or not 0 <= x < 2**64 for x in shift):
            raise ValueError("shift must have two unsigned exact integer coordinates")
    if (not isinstance(anchor, tuple) or len(anchor) != 3
            or any(type(value) is not int for value in anchor) or anchor[2] != 0):
        raise ValueError("anchor must be an integer (q,t,0) tuple")
    producer = capture_provenance()
    identity, polynomial = presentation_identity(pres), target_identity(pres, target)
    root = Path(config.output_root)
    prepare_output_root(root)
    order = {"solver": config.monomial_order, "reduction": REDUCTION_ORDER}
    samples = SampleSet(pres, target, anchor, config, producer=producer)
    started = time.monotonic()
    stage = "solver"

    def finish(status, reason=None, **data):
        report = {"schema": "zipper-discovery-result-v2", "status": status,
                  "verified": status == "verified", "in_ideal": True if status == "verified" else None,
                  "basis_complete": False, "method": "proof-contexts-finiteflow",
                  "producer": producer,
                  "presentation": identity, "target": polynomial, "target_name": target,
                  "anchor": list(samples.anchor), "points": len(samples.points),
                  "runs": _solver_runs(samples),
                  "seconds": time.monotonic() - started, **data}
        if reason is not None:
            report["reason"] = reason
        write_json_atomic(root / "result.json", report)
        return report

    try:
        samples.run_anchor(refresh=refresh)
        skeleton = _skeleton(pres, target, samples)
        save_artifact(root / "candidate-skeleton.json", skeleton, pres, target, order=order,
                      producer=producer)
        if not skeleton.in_ideal:
            return finish("inconclusive", "nonzero remainder against an uncertified basis",
                          remainder_terms=skeleton.remainder_terms)
        stage = "support"
        system = build_system(pres, target, skeleton, samples.anchor_basis(),
                              samples.anchor_proof(), max_contexts=max_contexts,
                              max_products=max_products)
        stage = "reconstruction"
        dependencies = reconstruction_provenance()
        reconstruction = {"schema": "zipper-reconstruction-run-v1", "dependencies": dependencies,
                          "options": {"threads": threads, "max_primes": max_primes,
                                      "max_degree": reconstruction_degree, "start_prime": start_prime,
                                      "shift": None if shift is None else list(shift),
                                      "max_contexts": max_contexts, "max_products": max_products}}
        discovery = reconstruct_cofactor(pres, target, system, anchor=samples.anchor,
                         max_primes=max_primes, max_degree=reconstruction_degree,
                         threads=threads, max_products=max_products, start_prime=start_prime,
                         shift=shift)
        # This independent free-algebra check also protects the publication boundary.
        report = check(discovery.cofactor, pres, target, verbose=verbose)
        if not report["pass"]:
            return finish("inconclusive", "candidate failed exact verification")
        stage = "publication"
        if reconstruction_provenance() != dependencies:
            raise ArtifactError("reconstruction dependencies changed during discovery")
        validate_source_snapshot(producer)
        try:
            products = _publish_verified(root, pres, target, skeleton, discovery, order, producer,
                                         solver_runs=_solver_runs(samples),
                                         reconstruction=reconstruction)
        except OSError as error:
            raise ArtifactError(f"could not publish verified generation: {error}") from error
        return finish("verified", exact=report, tuples=len(discovery.cofactor.terms),
                      reconstruction=reconstruction, support=system.statistics, **products)
    except (ExecutionError, ArtifactError) as error:
        return finish("inconclusive", str(error), execution_status="failed", stage=stage)
    except (RuntimeError, ImportError) as error:
        return finish("inconclusive", str(error), stage=stage)
