"""Single-thread C++ execution of the GPU prototype's exact packed programs.

Admission, packing, collection, reconstruction and verification are shared with
the OpenCL prototype. Only instruction execution changes, allowing a native CPU
comparison without claiming a full native scheduler implementation.
"""

from __future__ import annotations

import ctypes
import fcntl
import hashlib
import json
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
from time import perf_counter

from zipper_algebra.zipping_algorithm.gpu_rewrite import ExpansionResult, OpenClRewriteExpander
from zipper_algebra.zipping_algorithm.interpreter import ReplayError


SOURCE = Path(__file__).with_name("native") / "zipper_rewrite_cpu.cpp"
FLAGS = ("-O3", "-std=c++17", "-shared", "-fPIC", "-Wall", "-Wextra")
ABI = 0x5A495050455201


class NativeBuildError(RuntimeError):
    """Native compiler, library, or provenance contract failed."""


def _manifest(path: Path) -> dict:
    try:
        record = json.loads(path.with_suffix(".json").read_text())
        if (record["source_sha256"] != hashlib.sha256(SOURCE.read_bytes()).hexdigest()
                or record["binary_sha256"] != hashlib.sha256(path.read_bytes()).hexdigest()
                or record["abi"] != ABI):
            raise NativeBuildError("Native library source/binary/ABI provenance differs")
        return record
    except (OSError, ValueError, KeyError) as error:
        raise NativeBuildError("Native library requires a valid matching provenance manifest") from error


def build_library(*, compiler: str = "g++", build_directory: Path | None = None,
                  flags: tuple[str, ...] = FLAGS):
    """Build once into a source/compiler/flags-bound project-local cache."""
    executable = shutil.which(compiler)
    if executable is None:
        raise NativeBuildError("C++ compiler is unavailable: " + compiler)
    try:
        version = subprocess.run([executable, "--version"], check=True,
                                 capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise NativeBuildError("Cannot identify C++ compiler") from error
    source_sha = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    identity = {"source_sha256": source_sha, "compiler": executable,
                "compiler_version": version, "flags": list(flags), "abi": ABI}
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    base = (Path(build_directory) if build_directory is not None else
            Path(__file__).resolve().parents[2] / ".local-tools" / "native-rewrite")
    directory = base / key
    directory.mkdir(parents=True, exist_ok=True)
    library = directory / "libzipper_rewrite_cpu.so"
    with (directory / "build.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if library.exists():
            return library, _manifest(library), False
        started = perf_counter()
        with tempfile.TemporaryDirectory(prefix="compile-", dir=directory) as temporary:
            candidate = Path(temporary) / library.name
            try:
                subprocess.run([executable, *flags, str(SOURCE), "-o", str(candidate)],
                               check=True, capture_output=True, text=True)
            except (OSError, subprocess.CalledProcessError) as error:
                detail = getattr(error, "stderr", "") or str(error)
                raise NativeBuildError("Native C++ compilation failed: " + detail) from error
            record = dict(identity,
                          binary_sha256=hashlib.sha256(candidate.read_bytes()).hexdigest(),
                          build_seconds=perf_counter() - started)
            manifest = Path(temporary) / "manifest.json"
            manifest.write_text(json.dumps(record, indent=2) + "\n")
            candidate.replace(library)
            manifest.replace(library.with_suffix(".json"))
    return library, record, True


class NativeRewriteExpander(OpenClRewriteExpander):
    """Execute the retained packed programs on one C++ CPU thread."""

    def __init__(self, library: Path | None = None, *, compiler: str = "g++",
                 build_directory: Path | None = None, flags: tuple[str, ...] = FLAGS):
        start = perf_counter()
        if sys.byteorder != "little":
            raise NativeBuildError("Packed native prototype requires a little-endian host")
        if library is None:
            path, record, compiled = build_library(compiler=compiler,
                                                  build_directory=build_directory, flags=flags)
        else:
            path = Path(library).resolve()
            record, compiled = _manifest(path), False
        try:
            self.library = ctypes.CDLL(str(path))
            self.library.zipper_rewrite_abi.argtypes = []
            self.library.zipper_rewrite_abi.restype = ctypes.c_uint64
            if self.library.zipper_rewrite_abi() != ABI:
                raise NativeBuildError("Native packed rewrite ABI differs")
            void = ctypes.c_void_p
            integer = ctypes.c_int64
            self.library.zipper_evaluate_registers.argtypes = [void] * 6 + [integer]
            self.library.zipper_expand_words.argtypes = [void] * 10 + [integer]
            self.library.zipper_propagate_coefficients.argtypes = [void] * 9 + [integer, integer]
            for name in ("zipper_evaluate_registers", "zipper_expand_words",
                         "zipper_propagate_coefficients"):
                getattr(self.library, name).restype = None
        except (OSError, AttributeError) as error:
            raise NativeBuildError("Cannot load native packed rewrite library") from error
        self.metadata = {
            "backend": "cpp", "name": "single-thread C++ packed rewrite interpreter",
            "threads": 1, "library": str(path), **record,
            "compiled_during_initialization": compiled,
            "wrapper_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "arithmetic": "exact signed multi-limb integers and checked i64 coordinates",
            "successor_word_generation": "native CPU symbolic program execution",
            "host_work": "shared Python guard/LHS admission, packing, duplicate collection and controls",
        }
        self._templates = []
        self._template_keys = {}
        self._template_arrays = {}
        self._template_buffers = {}
        self._program_term_indices = []
        self._term_factor_indices = []
        self._templates_dirty = False
        self.initialization_seconds = perf_counter() - start

    def expand(self, cases, *, verify=True):
        started = perf_counter()
        result = super().expand(cases, verify=verify)
        metadata = dict(result.metadata)
        metadata["native_output_bytes"] = metadata["d2h_bytes"]
        metadata["native_template_buffer_bytes"] = sum(map(len, self._template_arrays.values()))
        metadata["template_upload_happened"] = False
        metadata["template_upload_bytes"] = 0
        metadata["h2d_bytes"] = 0
        metadata["d2h_bytes"] = 0
        metadata["h2d_byte_measure"] = "native CPU execution; no host-to-device copies"
        timings = dict(result.timings)
        timings["total"] = perf_counter() - started
        return ExpansionResult(result.output, timings, metadata)

    def _execute(self, packed):
        start = perf_counter()
        # CPython bytes remain alive through these pointer wrappers. All inputs
        # use the same immutable little-endian buffers as the GPU path.
        buffers = {name: ctypes.c_char_p(payload or bytes(8))
                   for name, payload in packed.arrays.items()}
        templates = {name: ctypes.c_char_p(payload or bytes(8))
                     for name, payload in self._template_arrays.items()}
        sizes = {
            "word_output": packed.word_tokens * 6 * 8,
            "coefficient_output": packed.coefficient_records * packed.limbs * 8,
            "q_output": packed.coefficient_records * 8,
            "status": len(packed.admitted) * 4,
        }
        output = {name: (ctypes.c_uint64 * max(1, (size + 7) // 8))()
                  for name, size in sizes.items()}
        values = (ctypes.c_int64 * max(1, packed.expression_values))()
        self._templates_dirty = False
        allocated = perf_counter()
        self.library.zipper_evaluate_registers(
            templates["programs"], templates["expressions"], buffers["cases"],
            buffers["parameters"], values, output["status"], len(packed.admitted))
        evaluated = perf_counter()
        self.library.zipper_expand_words(
            buffers["cases"], values, templates["terms"], templates["primitives"],
            buffers["words"], buffers["refs"], buffers["named_words"], buffers["word_jobs"],
            output["word_output"], output["status"], len(packed.word_jobs))
        expanded = perf_counter()
        self.library.zipper_propagate_coefficients(
            buffers["cases"], values, templates["factors"], buffers["old_q"],
            buffers["coefficients"], buffers["coefficient_jobs"], output["q_output"],
            output["coefficient_output"], output["status"], packed.limbs,
            len(packed.coefficient_jobs))
        propagated = perf_counter()
        raw = {name: ctypes.string_at(ctypes.addressof(output[name]), size)
               for name, size in sizes.items()}
        statuses = struct.unpack("<" + "I" * len(packed.admitted), raw["status"])
        if any(statuses):
            raise ReplayError("Native rewrite checked arithmetic/allocation failure: " + repr(statuses))
        timings = {"native_buffer_preparation": allocated - start,
                   "evaluate": evaluated - allocated,
                   "word_expansion": expanded - evaluated,
                   "coefficient_propagation": propagated - expanded,
                   "native_output_materialization": perf_counter() - propagated}
        return raw, timings, 0
