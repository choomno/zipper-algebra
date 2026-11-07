"""Separate C++/GMP replay verifier; no producer implementation is executed.

The executable reconstructs the entire calculation and compares its certificate.
This module supplies authenticated data and process transport only. The producer's
mathematical routines, state representation, caches and witnesses are not oracles.
"""
from __future__ import annotations

from contextlib import suppress
from copy import deepcopy
import fcntl
import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from threading import Lock
from time import perf_counter
from typing import Any, TextIO

ABI = "zipper-independent-verifier-v1"
FLAGS = ("-O3", "-std=c++17", "-pthread", "-Wall", "-Wextra")
LIBRARIES = ("gmpxx", "gmp", "crypto")
_REFERENCE_NAMES = ("algorithm.py", "catalog.py", "upward.py", "framed.py", "frames.py",
                    "turnaround.py", "transport.py", "interpreter.py", "extraction.py",
                    "semilinear.py", "collected.py")
_BUILD_FIELDS = ("abi", "binary_sha256", "compiler", "compiler_version", "flags",
                 "wrapper_source_sha256", "link_libraries")
_PACKAGE = Path(__file__).resolve().parent
_CHECKOUT_SOURCES = _PACKAGE.parents[2] / "native" / "zipper_verifier"
SOURCE_DIRECTORY = (_PACKAGE / "_verifier_sources" if (_PACKAGE / "_verifier_sources").is_dir()
                    else _CHECKOUT_SOURCES)
_PRODUCER_SOURCES = (_PACKAGE / "_native_sources" if (_PACKAGE / "_native_sources").is_dir()
                     else _PACKAGE.parents[2] / "native" / "zipper_engine")


class VerifierError(ValueError):
    """Invalid certificate, failed independent replay, or process contract."""


class VerifierBuildError(VerifierError):
    """The verifier cannot be built or its source/executable binding differs."""


def _hash(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _same(a: Any, b: Any) -> bool:
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(_same(value, b[key]) for key, value in a.items())
    if isinstance(a, list):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return bool(a == b)


def _sources(directory: Path) -> dict[str, str]:
    paths = sorted(p for p in directory.rglob("*") if p.is_file() and p.suffix in (".cpp", ".hpp", ".h"))
    if not paths or not (directory / "main.cpp").is_file():
        raise VerifierBuildError("Native source directory is incomplete: " + str(directory))
    if any(not path.resolve().is_relative_to(directory.resolve()) for path in paths):
        raise VerifierBuildError("Native source path leaves its authenticated directory")
    return {path.relative_to(directory).as_posix(): _hash(path) for path in paths}


def source_identity() -> dict[str, str]:
    return _sources(SOURCE_DIRECTORY)


def reference_identity() -> dict[str, str]:
    """Reference hashes are provenance, never supplied transition computations."""
    return {name: _hash(_PACKAGE / name) for name in _REFERENCE_NAMES}


def _key(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text())
        if not isinstance(value, dict):
            raise ValueError("object required")
        return value
    except (OSError, ValueError, TypeError) as error:
        raise VerifierBuildError("Invalid build manifest: " + str(path)) from error


def _compiler_identity(compiler: str, flags: tuple[str, ...]) -> dict[str, Any]:
    executable = shutil.which(compiler)
    if executable is None:
        raise VerifierBuildError("C++ compiler is unavailable: " + compiler)
    executable = str(Path(executable).resolve())
    try:
        version = subprocess.run([executable, "--version"], capture_output=True, text=True, check=True).stdout.strip()
        dependencies = {}
        for name in LIBRARIES:
            path = Path(subprocess.run([executable, "-print-file-name=lib" + name + ".so"],
                                      capture_output=True, text=True, check=True).stdout.strip()).resolve()
            if not path.is_file():
                raise VerifierBuildError("Compiler cannot locate library: " + name)
            dependencies[name] = {"path": str(path), "sha256": _hash(path)}
    except (OSError, subprocess.CalledProcessError) as error:
        raise VerifierBuildError("Cannot identify verifier compiler/dependencies") from error
    return {"abi": ABI, "compiler": executable, "compiler_version": version,
            "flags": list(flags), "link_libraries": list(LIBRARIES), "libraries": dependencies}


def _authenticated_binary(path: Path, expected: dict[str, Any] | None = None) -> dict[str, Any]:
    record = _load(path.with_suffix(".json"))
    if (record.get("abi") != ABI or record.get("source_sha256") != source_identity()
            or record.get("binary_sha256") != _hash(path)
            or expected is not None and any(not _same(record.get(key), value) for key, value in expected.items())):
        raise VerifierBuildError("Verifier executable/source/build identity differs")
    if not isinstance(record.get("libraries"), dict) or set(record["libraries"]) != set(LIBRARIES):
        raise VerifierBuildError("Verifier library identity is missing")
    for value in record["libraries"].values():
        if not isinstance(value, dict) or set(value) != {"path", "sha256"} or _hash(Path(value["path"])) != value["sha256"]:
            raise VerifierBuildError("Verifier linked library identity differs")
    return record


def build_verifier(*, compiler: str = "g++", build_directory: Path | None = None,
                   flags: tuple[str, ...] = FLAGS) -> tuple[Path, dict[str, Any], bool]:
    identity = dict(_compiler_identity(compiler, flags), source_sha256=source_identity())
    if build_directory is not None:
        base = Path(build_directory)
    elif SOURCE_DIRECTORY == _CHECKOUT_SOURCES:
        base = SOURCE_DIRECTORY.parents[2] / ".local-tools" / "native-verifier"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache"))) / "zipper-algebra" / "native-verifier"
    directory = base / _key(identity)
    directory.mkdir(parents=True, exist_ok=True)
    binary = directory / "zipper-verifier"
    with (directory / "build.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if binary.exists() and not binary.with_suffix(".json").exists():
            binary.unlink()
        if binary.exists():
            return binary, _authenticated_binary(binary, identity), False
        started = perf_counter()
        objects = []
        compiled = 0
        for source in sorted(SOURCE_DIRECTORY.glob("*.cpp")):
            inputs = {name: digest for name, digest in identity["source_sha256"].items()
                      if name == source.name or Path(name).suffix in (".hpp", ".h")}
            object_identity = {key: identity[key] for key in ("compiler", "compiler_version", "flags", "libraries")}
            object_identity.update(source=source.name, inputs=inputs)
            object_directory = base / "objects" / _key(object_identity)
            object_directory.mkdir(parents=True, exist_ok=True)
            target = object_directory / (source.stem + ".o")
            with (object_directory / "build.lock").open("a") as object_lock:
                fcntl.flock(object_lock, fcntl.LOCK_EX)
                if target.exists() and not target.with_suffix(".json").exists():
                    target.unlink()
                if target.exists():
                    recorded = _load(target.with_suffix(".json"))
                    if not _same(recorded.get("identity"), object_identity) or recorded.get("sha256") != _hash(target):
                        raise VerifierBuildError("Verifier object cache identity differs")
                else:
                    with tempfile.TemporaryDirectory(dir=object_directory, prefix="build-") as temp:
                        candidate = Path(temp) / target.name
                        _compile([identity["compiler"], *flags, "-c", str(source), "-o", str(candidate)])
                        if any(_hash(SOURCE_DIRECTORY / name) != digest for name, digest in inputs.items()):
                            raise VerifierBuildError("Verifier sources changed during compilation")
                        object_record = {"identity": object_identity, "sha256": _hash(candidate)}
                        manifest = Path(temp) / "manifest.json"
                        manifest.write_text(json.dumps(object_record, indent=2) + "\n")
                        manifest.replace(target.with_suffix(".json"))
                        candidate.replace(target)
                        compiled += 1
                objects.append(target)
        with tempfile.TemporaryDirectory(dir=directory, prefix="link-") as temp:
            candidate = Path(temp) / binary.name
            _compile([identity["compiler"], *flags, *map(str, objects),
                      *("-l" + name for name in LIBRARIES), "-o", str(candidate)])
            if identity["source_sha256"] != source_identity():
                raise VerifierBuildError("Verifier source changed during linking")
            record = dict(identity, binary_sha256=_hash(candidate), build_seconds=perf_counter() - started,
                          compiled_translation_units=compiled, reused_translation_units=len(objects) - compiled)
            manifest = Path(temp) / "manifest.json"
            manifest.write_text(json.dumps(record, indent=2) + "\n")
            manifest.replace(binary.with_suffix(".json"))
            candidate.replace(binary)
    return binary, record, True


def _compile(command: list[str]) -> None:
    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise VerifierBuildError("Verifier compilation failed: " + (getattr(error, "stderr", "") or str(error))) from error


def _decimal(value: int) -> str:
    # Base 10**9 chunks avoid Python's process-global integer-string limit.
    negative = value < 0
    value = abs(value)
    chunks = []
    while value >= 10**9:
        value, digit = divmod(value, 10**9)
        chunks.append(f"{digit:09d}")
    return ("-" if negative else "") + str(value) + "".join(reversed(chunks))


def _validate_candidate(value: Any) -> None:
    """Preserve the native certificate's JSON types before transport encoding."""
    if isinstance(value, dict):
        if any(type(key) is not str for key in value):
            raise VerifierError("Certificate object keys must be strings")
        for key, item in value.items():
            if (key == "coefficient" and isinstance(item, list)
                    and all(isinstance(pair, list) and len(pair) == 2 for pair in item)
                    and any(type(part) is not int for pair in item for part in pair)):
                raise VerifierError("Certificate coefficient pairs require integer scalars and exponents")
            _validate_candidate(item)
    elif isinstance(value, list):
        for item in value:
            _validate_candidate(item)
    elif value is not None and type(value) not in (str, bool, int, float):
        raise VerifierError("Certificate values must retain JSON types")


def _wire(value: Any) -> Any:
    if isinstance(value, (tuple, list)):
        return [_wire(x) for x in value]
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key == "coefficient" and isinstance(item, (list, tuple)):
                result[key] = [[q, _decimal(c)] if isinstance(pair, (list, tuple)) and len(pair) == 2
                               and type(q := pair[0]) is int and type(c := pair[1]) is int else _wire(pair)
                               for pair in item]
            else:
                result[key] = _wire(item)
        return result
    return value


def _producer_identity(metadata: dict[str, Any], catalog_sha: str, references: dict[str, str]) -> dict[str, Any]:
    if not isinstance(metadata, dict) or type(metadata.get("binary")) is not str:
        raise VerifierError("An authenticated producer build metadata record is required")
    binary = Path(metadata["binary"])
    manifest = _load(binary.with_suffix(".json"))
    expected_sources = _sources(_PRODUCER_SOURCES)
    if (metadata.get("source_sha256") != expected_sources or manifest.get("source_sha256") != expected_sources
            or metadata.get("binary_sha256") != _hash(binary)
            or metadata.get("wrapper_source_sha256") != _hash(_PACKAGE / "native.py")
            or any(not _same(metadata.get(key), manifest.get(key)) for key in _BUILD_FIELDS if key != "wrapper_source_sha256")):
        raise VerifierError("Producer build/source provenance differs")
    return {"schema": "zipper-native-run-v1", "scope": "literal_word_native_collected_spanning_procedure",
            "source_sha256": expected_sources, "python_reference_source_sha256": references,
            "catalog_sha256": catalog_sha, "coefficient_ring": "Z[q,q^-1]",
            "native_build": {key: deepcopy(metadata[key]) for key in _BUILD_FIELDS}}


def _allocator_identity(library: Path) -> dict[str, str]:
    try:
        path = library.resolve(strict=True)
        if any(character.isspace() or character == ":" for character in str(path)):
            raise VerifierError("Allocator library path contains LD_PRELOAD separators")
        return {"path": str(path), "sha256": _hash(path)}
    except OSError as error:
        raise VerifierError("Allocator library is unavailable") from error


def _allocator_mapped(pid: int, identity: dict[str, str]) -> bool:
    try:
        stat = Path(identity["path"]).stat()
        mappings = Path(f"/proc/{pid}/maps").read_text().splitlines()
        for line in mappings:
            fields = line.split(maxsplit=5)
            if len(fields) < 6:
                continue
            # Overlay mounts can report a different st_dev from /proc/maps.
            # Require both the resolved requested pathname and the file inode.
            if int(fields[4]) == stat.st_ino and fields[5] == identity["path"]:
                return True
        return False
    except (OSError, ValueError) as error:
        raise VerifierError("Cannot authenticate loaded allocator mapping") from error


class CppVerifier:
    """Persistent independent C++/GMP replay process with its own bounded cache."""

    def __init__(self, *, producer_metadata: dict[str, Any], cache_bytes: int = 128 * 1024**2,
                 cache_entries: int = 262144, binary: str | Path | None = None,
                 compiler: str = "g++", build_directory: Path | None = None,
                 flags: tuple[str, ...] = FLAGS, allocator_library: Path | None = None) -> None:
        if any(type(value) is not int or not 0 <= value <= 2**63 - 1 for value in (cache_bytes, cache_entries)):
            raise VerifierError("Verifier cache budgets must be nonnegative signed 64-bit integers")
        started = perf_counter()
        self._lock = Lock()
        self._process: subprocess.Popen[str] | None = None
        self._stderr: TextIO | None = None
        self._allocator = None if allocator_library is None else _allocator_identity(allocator_library)
        self._reference = reference_identity()
        self._wrapper_sha = _hash(Path(__file__))
        payload = files(__package__).joinpath("data/full-catalog.json").read_bytes()
        self._catalog_sha = hashlib.sha256(payload).hexdigest()
        self._identity = _producer_identity(producer_metadata, self._catalog_sha, self._reference)
        if binary is None:
            path, record, compiled = build_verifier(compiler=compiler, build_directory=build_directory, flags=flags)
        else:
            path = Path(binary).resolve()
            record, compiled = _authenticated_binary(path), False
        self.metadata = dict(record, binary=str(path), wrapper_source_sha256=self._wrapper_sha,
                             compiled_during_initialization=compiled, reference_source_sha256=self._reference,
                             allocator_library=deepcopy(self._allocator))
        self.last_call_seconds = 0.0
        self.last_result: dict[str, Any] = {}
        child_environment = None
        if self._allocator is not None:
            child_environment = dict(os.environ, LD_PRELOAD=self._allocator["path"])
        try:
            self._stderr = tempfile.TemporaryFile(mode="w+t")
            self._process = subprocess.Popen([str(path)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                             stderr=self._stderr, encoding="utf-8", text=True, bufsize=1,
                                             env=child_environment)
            ready = self._request({"op": "initialize", "abi": ABI, "catalog": json.loads(payload),
                                   "catalog_sha256": self._catalog_sha, "catalog_text": payload.decode("utf-8"),
                                   "reference_sources": self._reference,
                                   "native_identity": self._identity,
                                   "cache_entries": cache_entries, "cache_bytes": cache_bytes})
            if not _same(ready, {"abi": ABI, "ready": True}):
                raise VerifierError("Verifier initialization response differs")
            if self._allocator is not None:
                if (_allocator_identity(Path(self._allocator["path"])) != self._allocator
                        or not _allocator_mapped(self._process.pid, self._allocator)):
                    raise VerifierError("Requested allocator library is not authenticated in verifier mappings")
                self.metadata["allocator_library"]["mapped_verified"] = True
        except BaseException:
            self.close()
            raise
        self.initialization_seconds = perf_counter() - started

    @property
    def pid(self) -> int | None:
        return self._process.pid if self._process is not None else None

    def _request(self, value: dict[str, Any]) -> dict[str, Any]:
        process = self._process
        if process is None:
            raise VerifierError("Verifier session is closed")
        try:
            encoded = json.dumps(value, separators=(",", ":"), allow_nan=False) + "\n"
        except (ValueError, TypeError) as error:
            raise VerifierError("Verifier request cannot be encoded as JSON") from error
        try:
            if process.stdin is None or process.stdout is None:
                raise VerifierError("Verifier pipes are unavailable")
            process.stdin.write(encoded)
            process.stdin.flush()
            line = process.stdout.readline()
            if not line:
                raise VerifierError("Verifier exited without a response")
            result = json.loads(line)
            if not isinstance(result, dict):
                raise VerifierError("Verifier protocol requires an object")
            return result
        except BaseException:
            self._close_unlocked(abort=True)
            raise

    def _assert_source(self) -> None:
        if self._allocator is not None:
            try:
                if _allocator_identity(Path(self._allocator["path"])) != self._allocator:
                    raise VerifierError("Allocator library changed")
            except VerifierError:
                self._close_unlocked(abort=True)
                raise
        if (source_identity() != self.metadata["source_sha256"] or reference_identity() != self._reference
                or _hash(Path(__file__)) != self._wrapper_sha
                or _sources(_PRODUCER_SOURCES) != self._identity["source_sha256"]
                or _hash(_PACKAGE / "native.py") != self._identity["native_build"]["wrapper_source_sha256"]
                or hashlib.sha256(files(__package__).joinpath("data/full-catalog.json").read_bytes()).hexdigest() != self._catalog_sha):
            raise VerifierError("Verifier, reference, producer, or catalog sources changed")

    def check(self, certificate: dict[str, Any]) -> bool:
        """Replay and compare in C++; no producer/Python mathematical fallback."""
        if not isinstance(certificate, dict):
            raise VerifierError("Certificate must be an object")
        _validate_candidate(certificate)
        with self._lock:
            self._assert_source()
            started = perf_counter()
            result = self._request({"op": "verify", "certificate": _wire(certificate)})
            self.last_call_seconds = perf_counter() - started
            self.last_result = result
            if result.get("ok") is not True:
                raise VerifierError(str(result.get("error", "Independent C++ replay rejected the certificate")))
            return True

    def diagnostic(self, operation: str, **parameters: Any) -> dict[str, Any]:
        """Execute independent replay/transition diagnostics for differential tests."""
        if operation not in ("replay", "transition"):
            raise VerifierError("Unknown verifier diagnostic")
        with self._lock:
            self._assert_source()
            result = self._request(dict(op=operation, **_wire(parameters)))
            if result.get("ok") is False:
                raise VerifierError(str(result.get("error", "Verifier diagnostic failed")))
            return result

    def cache_info(self) -> dict[str, Any]:
        with self._lock:
            return self._request({"op": "cache_info"})

    def clear_cache(self) -> None:
        with self._lock:
            if not _same(self._request({"op": "clear_cache"}), {"cleared": True}):
                raise VerifierError("Verifier cache reset response differs")

    def _close_unlocked(self, *, abort: bool) -> None:
        process = self._process
        stderr = None
        try:
            self._process = None
            stderr = getattr(self, "_stderr", None)
            self._stderr = None
            if process is not None:
                if abort:
                    with suppress(OSError):
                        process.kill()
                if process.stdin is not None:
                    with suppress(OSError, ValueError):
                        process.stdin.close()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=5)
        except BaseException:
            if process is not None:
                with suppress(BaseException):
                    process.kill()
                with suppress(BaseException):
                    process.wait(timeout=5)
            raise
        finally:
            # Retain the original pipe if an interrupt arrived while reading
            # its attribute, before the local handle could be assigned.
            saved_stderr = self.__dict__.pop("_stderr", None)
            if stderr is None:
                stderr = saved_stderr
            if process is not None:
                for pipe in (process.stdin, process.stdout):
                    if pipe is not None:
                        with suppress(OSError, ValueError):
                            pipe.close()
            if stderr is not None:
                with suppress(OSError, ValueError):
                    stderr.close()

    def close(self) -> None:
        with self._lock:
            self._close_unlocked(abort=False)

    def __enter__(self) -> CppVerifier:
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def __del__(self) -> None:
        if getattr(self, "_process", None) is not None:
            self.close()
