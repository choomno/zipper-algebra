"""Run f4ncgb at an exact rational parameter point and validate its proof data."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile
import time

from .artifacts import (ArtifactError, canonical_json, digest_bytes, file_identity, load_record,
                        publish, validate_proofs)
from .config import RunConfig
from ..cofactor.artifacts import presentation_identity, target_identity
from ..cofactor.cofactor_io import resolve_target, target_name
from ..coefficients import NumericRing


class ExecutionError(RuntimeError):
    """The external solver did not produce usable proof output."""


@dataclass(frozen=True)
class SolverRun:
    key: str
    directory: Path
    record: dict
    hit: bool
    seconds: float

    @property
    def basis_path(self) -> Path:
        return self.directory / "basis.txt"

    @property
    def proof_path(self) -> Path:
        return self.directory / "proof.txt"


def run_point(pres, target: str, point: tuple[int, int, int], config: RunConfig,
              *, refresh: bool = False) -> SolverRun:
    """Emit the presentation, execute the solver and publish validated proof data.

    ``target`` identifies the intended computation but is never appended to the
    defining relations. Only characteristic-zero points are supported here.
    """
    if (not isinstance(point, tuple) or len(point) != 3
            or any(type(value) is not int for value in point)):
        raise ValueError("point must be an integer (q, t, characteristic) tuple")
    if point[2] != 0:
        raise ValueError("the reconstruction runner supports characteristic zero only")
    ring = NumericRing(*point)
    if not isinstance(target, str):
        raise ValueError("the runner requires a named target")
    resolve_target(pres, target, ring)
    names, inputs, _, _ = pres.build(ring)
    text = ",".join(names) + "\n0\n" + ",\n".join(
        polynomial.to_msolve_string(names) for polynomial in inputs) + "\n"
    input_bytes = text.encode()
    request = {
        "schema": "zipper-solver-request-v1",
        "presentation": presentation_identity(pres),
        "point": list(point), "specialized_domain": "Q", "target": target_name(target),
        "target_polynomial": target_identity(pres, target),
        "monomial_order": config.monomial_order, "input_sha256": digest_bytes(input_bytes),
        "solver": {"binary": str(config.binary), "identity": file_identity(config.binary),
                   "revision": config.solver_revision,
                   "prefix": [{"path": str(path), **file_identity(path)}
                              for path in config.command_prefix],
                   "flags": list(config.arguments()), "timeout_seconds": config.timeout_seconds},
    }
    key = digest_bytes(canonical_json(request))
    if not refresh:
        loaded = load_record(config.output_root, key, request, inputs, names)
        if loaded is not None:
            return SolverRun(key, *loaded, True, 0)
        if config.cache_root is not None:
            loaded = load_record(config.cache_root, key, request, inputs, names)
            if loaded is not None:
                directory, record = publish(config.output_root, key, loaded[0], request,
                                            loaded[1]["basis_entries"])
                return SolverRun(key, directory, record, True, 0)
    config.output_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".solver-", dir=config.output_root) as temporary:
        stage = Path(temporary)
        input_path = stage / "input.ms"
        input_path.write_bytes(input_bytes)
        command = [str(config.binary), *(str(path) for path in config.command_prefix),
                   "-i", str(input_path), "-o", str(stage / "basis.txt"),
                   "-p", str(stage / "proof.txt"), *config.arguments()]
        started = time.monotonic()
        try:
            with (stage / "stdout.txt").open("wb") as stdout, (stage / "stderr.txt").open("wb") as stderr:
                result = subprocess.run(command, stdout=stdout, stderr=stderr,
                                        timeout=config.timeout_seconds, cwd=stage)
        except subprocess.TimeoutExpired as error:
            raise ExecutionError(f"solver exceeded {config.timeout_seconds} seconds") from error
        except OSError as error:
            raise ExecutionError(f"solver could not start: {error}") from error
        seconds = time.monotonic() - started
        if result.returncode != 0:
            with (stage / "stderr.txt").open("rb") as stderr:
                stderr.seek(max(0, stderr.seek(0, 2) - 800))
                tail = stderr.read().decode(errors="replace")
            raise ExecutionError(f"solver exited with code {result.returncode}: "
                                 f"{tail}")
        if input_path.read_bytes() != input_bytes:
            raise ArtifactError("solver changed the requested presentation input")
        count = validate_proofs(stage, inputs, names)
        if config.cache_root is not None:
            publish(config.cache_root, key, stage, request, count)
        directory, record = publish(config.output_root, key, stage, request, count)
        return SolverRun(key, directory, record, False, seconds)
