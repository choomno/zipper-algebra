"""Produce and exactly verify a cofactor from a presentation and solver proofs."""
from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from pathlib import Path
import os
import tempfile
import time

from ..cofactor.artifacts import (presentation_identity, save_artifact, target_identity,
                                 write_json_atomic)
from ..cofactor import compose as composition
from ..cofactor.cofactor_io import resolve_target
from ..cofactor.lift import NeedMoreSamples, lift
from ..cofactor.skeleton import Skeleton, reduce
from ..cofactor.verify import check
from ..coefficients import NumericRing
from .artifacts import ArtifactError, canonical_json, digest_bytes, file_identity
from .runner import ExecutionError
from .sampling import SampleSet


REDUCTION_ORDER = "length-lex-index-v1"


def _skeleton(pres, target, samples):
    basis, proof = samples.anchor_basis(), samples.anchor_proof()
    polynomial = resolve_target(pres, target, NumericRing(*samples.anchor))
    steps, remainder = reduce(polynomial, basis)
    return Skeleton.from_reduction(pres, target, samples.anchor, steps, remainder, basis, proof)


def multiplier_samples(skeleton, samples, pres):
    """Exact multipliers for matching contextual occurrences, including repeated uses."""
    signatures = [(tuple(step["left"]), step["basis_key"], tuple(step["right"]))
                  for step in skeleton.steps]
    result = {}
    for point in samples.stable(skeleton):
        basis = samples.basis(*point)
        target = resolve_target(pres, skeleton.target, NumericRing(*point, 0))
        steps, remainder = reduce(target, basis)
        observed = [(tuple(left), basis.key[index], tuple(right))
                    for left, index, right, _ in steps]
        if not remainder.terms and observed == signatures:
            result[point] = [Fraction(coef.v) for _, _, _, coef in steps]
    return result


def _publish_verified(root, pres, target, skeleton, fitted, multipliers, certificate, order):
    """Publish one immutable product generation; the result record selects it later."""
    destination = root / "verified"
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".publish-", dir=destination) as temporary:
        stage = Path(temporary) / "generation"
        stage.mkdir()
        save_artifact(stage / "skeleton.json", skeleton, pres, target, order=order)
        save_artifact(stage / "lift.json", fitted, pres, target, order=order)
        write_json_atomic(stage / "multipliers.json", composition.multipliers_to_json(
            multipliers, skeleton, fitted.atoms, pres))
        save_artifact(stage / "certificate.json", certificate, pres, target)
        certificate.save(stage / "certificate.txt")
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


def run(pres, target, anchor, config, *, start=100, step=40, maxpts=260,
        jobs=1, atoms=None, refresh=False, verbose=True):
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
    identity, polynomial = presentation_identity(pres), target_identity(pres, target)
    root = Path(config.output_root)
    root.mkdir(parents=True, exist_ok=True)
    order = {"solver": config.monomial_order, "reduction": REDUCTION_ORDER}
    samples = SampleSet(pres, target, anchor, replace(config, output_root=root / "samples"), atoms)
    started = time.monotonic()

    def finish(status, reason=None, **data):
        report = {"schema": "zipper-discovery-result-v1", "status": status,
                  "verified": status == "verified", "in_ideal": True if status == "verified" else None,
                  "basis_complete": False,
                  "presentation": identity, "target": polynomial, "target_name": target,
                  "anchor": list(samples.anchor), "points": len(samples.points),
                  "runs": [{"point": list(point), "key": result.key,
                            "record": str(result.directory.parent / "record.json"),
                            "generation": result.record["generation"],
                            "directory": str(result.directory), "files": result.record["files"],
                            "request": result.record["request"]}
                           for point, result in samples._runs.items()],
                  "seconds": time.monotonic() - started, **data}
        if reason is not None:
            report["reason"] = reason
        write_json_atomic(root / "result.json", report)
        return report

    try:
        samples.run_anchor(refresh=refresh)
        skeleton = _skeleton(pres, target, samples)
        save_artifact(root / "candidate-skeleton.json", skeleton, pres, target, order=order)
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
                    save_artifact(root / "candidate-skeleton.json", skeleton, pres, target, order=order)
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
                    products = _publish_verified(root, pres, target, skeleton, fitted,
                                                 multipliers, candidate, order)
                    return finish("verified", exact=report, tuples=len(candidate.terms), **products)
                reason = "candidate failed exact verification"
            if desired >= maxpts:
                return finish("inconclusive", reason)
            desired = min(desired + step, maxpts)
    except ExecutionError as error:
        return finish("inconclusive", str(error), execution_status="failed")
