"""Produce and exactly verify a cofactor from a presentation and solver proofs."""
from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import os
import tempfile
import time

from ..cofactor.artifacts import (presentation_identity, save_artifact, target_identity,
                                 write_json_atomic)
from ..cofactor import compose as composition
from ..cofactor.atoms import AtomSet
from ..cofactor.cofactor_io import resolve_target
from ..cofactor.lift import NeedMoreSamples, lift
from ..cofactor.skeleton import Anchor, Skeleton, reduce
from ..cofactor.verify import check
from ..coefficients import NumericRing
from ..presentation import Presentation
from .artifacts import ArtifactError, canonical_json, digest_bytes, file_identity, prepare_output_root
from .runner import ExecutionError, validate_source_snapshot
from .sampling import Point, SampleSet
from .store import read as read_solver_run
from .config import RunConfig
from ..provenance import capture_provenance, write_bundle_provenance


REDUCTION_ORDER = "length-lex-index-v1"


def _skeleton(pres: Presentation, target: str, samples: SampleSet) -> Skeleton:
    basis, proof = samples.anchor_basis(), samples.anchor_proof()
    polynomial = resolve_target(pres, target, NumericRing(*samples.anchor))
    steps, remainder = reduce(polynomial, basis)
    return Skeleton.from_reduction(pres, target, samples.anchor, steps, remainder, basis, proof)


def multiplier_samples(skeleton: Skeleton, samples: SampleSet,
                       pres: Presentation) -> dict[Point, list[Fraction]]:
    """Exact multipliers for matching contextual occurrences, including repeated uses."""
    signatures = [(tuple(step["left"]), step["basis_key"], tuple(step["right"]))
                  for step in skeleton.steps]
    result: dict[Point, list[Fraction]] = {}
    for point in samples.stable(skeleton):
        basis = samples.basis(*point)
        target = resolve_target(pres, skeleton.target, NumericRing(*point, 0))
        steps, remainder = reduce(target, basis)
        observed = [(tuple(left), basis.key[index], tuple(right))
                    for left, index, right, _ in steps]
        if not remainder.terms and observed == signatures:
            result[point] = [Fraction(coef.v) for _, _, _, coef in steps]
    return result


def _solver_runs(samples: SampleSet) -> list[dict]:
    """Identify all validated runs of this discovery, including unstable samples."""
    return [{"point": list(point), "key": result.key, "store": str(result.store),
             "generation": result.record["generation"], "files": result.record["files"],
             "request": result.record["request"]}
            for point, result in samples._runs.items()]


def _solver_links(root, pres, target, solver_runs, producer):
    """Bind supplied run records to retained output-root-relative solver data.

    Reading a stored generation checks its identities and bytes; it does not
    establish basis completeness or which samples contributed to a fitted slot.
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


def _publish_verified(root, pres, target, skeleton, fitted, multipliers, certificate, order,
                      producer=None, *, solver_runs=None):
    """Publish one immutable product generation; the result record selects it later."""
    destination = root / "verified"
    producer = capture_provenance() if producer is None else producer
    solver_links = _solver_links(root, pres, target, solver_runs, producer)
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".publish-", dir=destination) as temporary:
        stage = Path(temporary) / "generation"
        stage.mkdir()
        save_artifact(stage / "skeleton.json", skeleton, pres, target, order=order, producer=producer)
        save_artifact(stage / "lift.json", fitted, pres, target, order=order, producer=producer)
        write_json_atomic(stage / "multipliers.json", composition.multipliers_to_json(
            multipliers, skeleton, fitted.atoms, pres))
        save_artifact(stage / "certificate.json", certificate, pres, target, producer=producer)
        certificate.save(stage / "certificate.txt", provenance=False)
        write_json_atomic(stage / "solver-runs.json", solver_links)
        write_bundle_provenance(stage, (path.name for path in stage.iterdir()), producer=producer)
        validate_source_snapshot(producer)
        files = {path.name: file_identity(path) for path in stage.iterdir()}
        generation = digest_bytes(canonical_json(files))
        final = destination / generation
        if final.exists():
            if any(file_identity(final / name) != identity for name, identity in files.items()):
                raise ArtifactError("existing verified generation has different content")
        else:
            try:
                os.replace(stage, final)
            except OSError as error:
                # Another producer may have published these same immutable
                # bytes after the existence check. Only that complete winner
                # permits us to treat a failed rename as successful reuse.
                try:
                    winner_matches = final.is_dir() and all(
                        file_identity(final / name) == identity for name, identity in files.items())
                except OSError:
                    raise error
                if not winner_matches:
                    raise
    return {"generation": generation, "files": files,
            "artifacts": {name: str(final / name) for name in files},
            "certificate": str(final / "certificate.txt")}


def run(pres: Presentation, target: str, anchor: Anchor, config: RunConfig, *,
        start: int = 100, step: int = 40, maxpts: int = 260, jobs: int = 1,
        atoms: AtomSet | None = None, refresh: bool = False, verbose: bool = True) -> dict:
    """Run bounded discovery; success means an exact identity, not a complete basis.

    Sampling and proof data are produced from the presentation. Saved lifts or
    historical cofactor expressions are not inputs to this procedure.
    """
    if any(type(value) is not int or value < 1 for value in (start, step, maxpts, jobs)):
        raise ValueError("start, step, maxpts and jobs must be positive integers")
    if start > maxpts:
        raise ValueError("start must not exceed maxpts")
    if (not isinstance(anchor, tuple) or len(anchor) != 3
            or any(type(value) is not int for value in anchor) or anchor[2] != 0):
        raise ValueError("anchor must be an integer (q,t,0) tuple")
    producer = capture_provenance()
    identity, polynomial = presentation_identity(pres), target_identity(pres, target)
    root = Path(config.output_root)
    prepare_output_root(root)
    order = {"solver": config.monomial_order, "reduction": REDUCTION_ORDER}
    samples = SampleSet(pres, target, anchor, config, atoms, producer=producer)
    started = time.monotonic()

    def finish(status, reason=None, **data):
        report = {"schema": "zipper-discovery-result-v1", "status": status,
                  "verified": status == "verified", "in_ideal": True if status == "verified" else None,
                  "basis_complete": False,
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
        desired = start
        while True:
            samples.grow(max(0, desired - len(samples.points)), jobs=jobs, refresh=refresh)
            modal, off = samples.modal_shapes(skeleton)
            if off:
                candidates = samples.stable(skeleton)
                if candidates:
                    samples.reanchor(candidates[0])
                    skeleton = _skeleton(pres, target, samples)
                    save_artifact(root / "candidate-skeleton.json", skeleton, pres, target, order=order,
                                  producer=producer)
                    if not skeleton.in_ideal:
                        return finish("inconclusive", "nonzero remainder at replacement anchor",
                                      remainder_terms=skeleton.remainder_terms)
                    # Rebuilding the reduction can introduce new dependencies. Its
                    # anchor must carry the modal shape for that new closure too.
                    if samples.modal_shapes(skeleton)[1]:
                        reason = "replacement anchor has incompatible proof shapes for its new closure"
                        if desired >= maxpts:
                            return finish("inconclusive", reason)
                        desired = min(desired + step, maxpts)
                        continue
                else:
                    reason = "no sample carries all modal proof shapes"
                    if desired >= maxpts:
                        return finish("inconclusive", reason)
                    desired = min(desired + step, maxpts)
                    continue
            try:
                fitted = lift(skeleton, samples, atoms=atoms, known="lines", verbose=verbose)
                values = multiplier_samples(skeleton, samples, pres)
                multipliers = composition.lift_multipliers(skeleton, values, fitted.atoms)
                candidate = composition.compose(fitted, skeleton, pres,
                                                multipliers=multipliers, verbose=verbose)
            except NeedMoreSamples as error:
                reason = str(error)
            else:
                report = check(candidate, pres, target, verbose=verbose)
                if report["pass"]:
                    validate_source_snapshot(producer)
                    products = _publish_verified(root, pres, target, skeleton, fitted,
                                                 multipliers, candidate, order, producer,
                                                 solver_runs=_solver_runs(samples))
                    return finish("verified", exact=report, tuples=len(candidate.terms), **products)
                reason = "candidate failed exact verification"
            if desired >= maxpts:
                return finish("inconclusive", reason)
            desired = min(desired + step, maxpts)
    except (ExecutionError, ArtifactError) as error:
        return finish("inconclusive", str(error), execution_status="failed")
