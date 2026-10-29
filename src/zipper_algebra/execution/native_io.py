"""Use f4ncgb's optional SQLite transport without intermediate solver files.

Transport rows are pending computation data. Only the runner's independent proof
checks allow their bytes to become a selected, validated sample generation.
"""
from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
import subprocess
import time
import uuid

from .artifacts import ArtifactError
from .store import transaction


class NativeExecutionError(RuntimeError):
    """The native database route did not finish its solver process successfully."""


def _discard(database, job_id):
    with transaction(database, write=True) as connection:
        connection.execute("DELETE FROM f4ncgb_io WHERE job_id=?", (job_id,))


@contextmanager
def _job(database, input_bytes):
    job_id = uuid.uuid4().hex
    with transaction(database, write=True) as connection:
        connection.execute("INSERT INTO f4ncgb_io VALUES (?, ?, NULL, NULL, 'ready')",
                           (job_id, input_bytes))
    try:
        yield job_id
    except BaseException:
        # Preserve the computation failure if storage is unavailable during
        # cleanup. Unique attempt ids prevent an old process selecting another
        # attempt's output, even when abrupt termination leaves a pending row.
        try:
            _discard(database, job_id)
        except (ArtifactError, OSError):
            pass
        raise
    else:
        _discard(database, job_id)


def execute(config, input_bytes, database):
    """Return raw solver artifacts and elapsed time from a unique database job."""
    if not isinstance(input_bytes, bytes):
        raise ArtifactError("native solver input must be exact bytes")
    database = Path(database)
    with _job(database, input_bytes) as job_id:
        command = [str(config.binary), *(str(path) for path in config.command_prefix),
                   "--sqlite", str(database), "--sqlite-job", job_id, *config.arguments()]
        started = time.monotonic()
        try:
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    timeout=config.timeout_seconds, cwd=config.output_root)
        except subprocess.TimeoutExpired as error:
            raise NativeExecutionError(f"solver exceeded {config.timeout_seconds} seconds") from error
        except OSError as error:
            raise NativeExecutionError(f"solver could not start: {error}") from error
        seconds = time.monotonic() - started
        if result.returncode != 0:
            tail = result.stderr[-800:].decode(errors="replace")
            raise NativeExecutionError(f"solver exited with code {result.returncode}: {tail}")
        with transaction(database) as connection:
            row = None if connection is None else connection.execute(
                "SELECT input,basis,proof,state FROM f4ncgb_io WHERE job_id=?", (job_id,)).fetchone()
        if row is None:
            raise ArtifactError("native solver job disappeared before completion")
        original, basis, proof, state = row
        if not isinstance(original, bytes) or original != input_bytes:
            raise ArtifactError("native solver changed the requested presentation input")
        if state != "complete":
            raise ArtifactError("native solver did not publish a complete result")
        if not isinstance(basis, bytes) or not isinstance(proof, bytes):
            raise ArtifactError("native solver basis and proof must be exact bytes")
        files = {"input.ms": input_bytes, "basis.txt": basis, "proof.txt": proof,
                 "stdout.txt": result.stdout, "stderr.txt": result.stderr}
    return files, seconds
