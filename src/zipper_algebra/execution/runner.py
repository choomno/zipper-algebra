"""Run f4ncgb at an exact rational parameter point and validate its proof data."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile
import time

from .artifacts import (ARTIFACT_NAMES, ArtifactError, canonical_json, digest_bytes, file_identity,
                        prepare_output_root, validate_proofs, validate_source_snapshot)
from .config import RunConfig
from .store import load, publish
from ..cofactor.artifacts import presentation_identity, target_identity
from ..cofactor.cofactor_io import resolve_target, target_name
from ..coefficients import NumericRing
from ..provenance import capture_provenance, validate_provenance


class ExecutionError(RuntimeError):
    """The external solver did not produce usable proof output."""


@dataclass(frozen=True)
class SolverRun:
    key: str
    store: Path
    files: dict[str, bytes]
    record: dict
    hit: bool
    seconds: float

def _validate_cached(files, record, inputs, names):
    count = validate_proofs(files, inputs, names)
    if count != record["basis_entries"]:
        raise ArtifactError("saved basis entry count differs from validated proof data")


def _run_file_solver(config, input_bytes):
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
        files = {name: (stage / name).read_bytes() for name in ARTIFACT_NAMES
                 if (stage / name).is_file()}
        return files, seconds


def run_point(pres, target: str, point: tuple[int, int, int], config: RunConfig,
              *, refresh: bool = False, producer=None) -> SolverRun:
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
    producer = capture_provenance() if producer is None else validate_provenance(producer)
    validate_source_snapshot(producer)
    request = {
        "schema": "zipper-solver-request-v2",
        "presentation": presentation_identity(pres),
        "point": list(point), "specialized_domain": "Q", "target": target_name(target),
        "target_polynomial": target_identity(pres, target),
        "monomial_order": config.monomial_order, "input_sha256": digest_bytes(input_bytes),
        "producer": producer,
        "solver": {"binary": str(config.binary), "identity": file_identity(config.binary),
                   "revision": config.solver_revision, "io_mode": config.solver_io,
                   "prefix": [{"path": str(path), **file_identity(path)}
                              for path in config.command_prefix],
                   "flags": list(config.arguments()), "timeout_seconds": config.timeout_seconds},
    }
    key = digest_bytes(canonical_json(request))
    prepare_output_root(config.output_root)
    output_store = config.output_root / "samples.sqlite3"
    if config.cache_root is not None:
        prepare_output_root(config.cache_root)
    if not refresh:
        loaded = load(output_store, key, request)
        if loaded is not None:
            files, record = loaded
            _validate_cached(files, record, inputs, names)
            validate_source_snapshot(producer)
            return SolverRun(key, output_store, files, record, True, 0)
        if config.cache_root is not None:
            loaded = load(config.cache_root / "samples.sqlite3", key, request)
            if loaded is not None:
                files, record = loaded
                _validate_cached(files, record, inputs, names)
                validate_source_snapshot(producer)
                record = publish(output_store, key, files, request, record["basis_entries"])
                return SolverRun(key, output_store, files, record, True, 0)
    if config.solver_io == "sqlite":
        from .native_io import NativeExecutionError, execute
        try:
            files, seconds = execute(config, input_bytes, output_store)
        except NativeExecutionError as error:
            raise ExecutionError(str(error)) from error
    else:
        files, seconds = _run_file_solver(config, input_bytes)
    if files.get("input.ms") != input_bytes:
        raise ArtifactError("solver changed the requested presentation input")
    count = validate_proofs(files, inputs, names)
    validate_source_snapshot(producer)
    if config.cache_root is not None:
        publish(config.cache_root / "samples.sqlite3", key, files, request, count)
    record = publish(output_store, key, files, request, count)
    return SolverRun(key, output_store, files, record, False, seconds)
