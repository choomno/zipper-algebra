"""Persistent C++ zipper state engine with an independent Python replay check.

Only complete requests and final certificates cross the process boundary. Packed
words, arbitrary-precision Laurent coefficients, transitions, collection and
scheduling remain in C++ for the entire calculation. Native builds are explicit,
content addressed and local; importing this module does not invoke a compiler.
"""
from __future__ import annotations

from contextlib import suppress
from copy import deepcopy
from dataclasses import asdict
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
from typing import Any, Iterable, TextIO, TYPE_CHECKING

from . import collected
from .catalog import Catalog
from .interpreter import Limits

if TYPE_CHECKING:
    from .independent_replay import IndependentVerifier

FLAGS = ("-O3", "-std=c++17", "-pthread", "-Wall", "-Wextra")
ABI = "zipper-native-engine-v1"
_CHECKOUT_SOURCES = Path(__file__).resolve().parents[3] / "native" / "zipper_engine"
_PACKAGED_SOURCES = Path(__file__).with_name("_native_sources")
SOURCE_DIRECTORY = (_PACKAGED_SOURCES if _PACKAGED_SOURCES.is_dir() else
                    _CHECKOUT_SOURCES)


class NativeError(collected.CollectedError):
    """The native engine, protocol, or certificate contract failed."""


class NativeBuildError(NativeError):
    """Native compiler or build provenance is unavailable or inconsistent."""


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_identity(source_directory: Path = SOURCE_DIRECTORY) -> dict[str, str]:
    """Hash all C++ sources and headers, including the pinned JSON dependency."""
    paths = sorted(p for p in source_directory.rglob("*")
                   if p.is_file() and p.suffix in (".cpp", ".hpp", ".h"))
    if not paths or not (source_directory / "main.cpp").is_file():
        raise NativeBuildError("The native engine source checkout is unavailable")
    return {str(p.relative_to(source_directory)): _digest(p) for p in paths}


def _manifest(binary: Path, identity: dict[str, Any] | None = None) -> dict[str, Any]:
    try:
        record = json.loads(binary.with_suffix(".json").read_text())
        if (not isinstance(record, dict)
                or any(type(record.get(key)) is not str for key in
                       ("abi", "binary_sha256", "compiler", "compiler_version"))
                or not isinstance(record.get("source_sha256"), dict)
                or any(type(key) is not str or type(value) is not str
                       for key, value in record["source_sha256"].items())
                or any(not isinstance(record.get(key), list)
                       or any(type(value) is not str for value in record[key])
                       for key in ("flags", "link_libraries"))):
            raise NativeBuildError("Native executable build metadata is malformed")
        if (record["abi"] != ABI or record["source_sha256"] != source_identity()
                or record["binary_sha256"] != _digest(binary)
                or identity is not None and any(record[k] != v for k, v in identity.items())):
            raise NativeBuildError("Native executable source/binary/build provenance differs")
        return record
    except NativeBuildError:
        raise
    except (OSError, ValueError, KeyError, TypeError) as error:
        raise NativeBuildError("Native executable requires a matching provenance manifest") from error


def _compile_object(source: Path, *, base: Path, identity: dict[str, Any],
                    executable: str, flags: tuple[str, ...]) -> tuple[Path, bool]:
    """Retain one object per translation unit and complete local header identity."""
    relative = str(source.relative_to(SOURCE_DIRECTORY))
    inputs = {key: value for key, value in identity["source_sha256"].items()
              if key == relative or Path(key).suffix in (".hpp", ".h")}
    object_identity = {key: identity[key] for key in ("compiler", "compiler_version", "flags")}
    object_identity.update(source=relative, inputs=inputs)
    key = hashlib.sha256(json.dumps(object_identity, sort_keys=True).encode()).hexdigest()
    directory = base / "objects" / key
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / (source.stem + ".o")
    with (directory / "build.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        # An older interrupted publication may leave an object without its
        # manifest. It is derived cache data: discard it and compile afresh.
        if target.exists() and not target.with_suffix(".json").exists():
            target.unlink()
        if target.exists():
            try:
                record = json.loads(target.with_suffix(".json").read_text())
                if (record["identity"] != object_identity or
                        record["object_sha256"] != _digest(target)):
                    raise NativeBuildError("Native object cache provenance differs")
            except (OSError, KeyError, TypeError, json.JSONDecodeError) as error:
                raise NativeBuildError("Native object cache manifest is invalid") from error
            return target, False
        with tempfile.TemporaryDirectory(prefix="compile-", dir=directory) as temporary:
            candidate = Path(temporary) / target.name
            try:
                subprocess.run([executable, *flags, "-I", str(SOURCE_DIRECTORY / "vendor"),
                                "-c", str(source), "-o", str(candidate)],
                               check=True, capture_output=True, text=True)
            except (OSError, subprocess.CalledProcessError) as error:
                detail = getattr(error, "stderr", "") or str(error)
                raise NativeBuildError("C++ translation-unit compilation failed: " + detail) from error
            if any(_digest(SOURCE_DIRECTORY / key) != value for key, value in inputs.items()):
                raise NativeBuildError("Native source or headers changed during object compilation")
            record = {"identity": object_identity, "object_sha256": _digest(candidate)}
            manifest = Path(temporary) / "manifest.json"
            manifest.write_text(json.dumps(record, indent=2) + "\n")
            # Publish the manifest first under the lock. An interruption then
            # leaves only metadata, never an executable unbound to its inputs.
            manifest.replace(target.with_suffix(".json"))
            candidate.replace(target)
    return target, True


def build_engine(*, compiler: str = "g++", build_directory: Path | None = None,
                 flags: tuple[str, ...] = FLAGS) -> tuple[Path, dict[str, Any], bool]:
    """Compile once per source/compiler/flags identity under a process lock."""
    executable = shutil.which(compiler)
    if executable is None:
        raise NativeBuildError("C++ compiler is unavailable: " + compiler)
    try:
        version = subprocess.run([executable, "--version"], check=True,
                                 capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise NativeBuildError("Cannot identify the C++ compiler") from error
    identity = {"abi": ABI, "source_sha256": source_identity(),
                "compiler": str(Path(executable).resolve()), "compiler_version": version,
                "flags": list(flags), "link_libraries": ["crypto"]}
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    if build_directory is not None:
        base = Path(build_directory)
    elif SOURCE_DIRECTORY == _CHECKOUT_SOURCES:
        base = SOURCE_DIRECTORY.parents[2] / ".local-tools" / "native-engine"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", str(Path.home() / ".cache"))) / "zipper-algebra" / "native-engine"
    directory = base / key
    directory.mkdir(parents=True, exist_ok=True)
    binary = directory / "zipper-engine"
    with (directory / "build.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if binary.exists() and not binary.with_suffix(".json").exists():
            binary.unlink()
        if binary.exists():
            return binary, _manifest(binary, identity), False
        started = perf_counter()
        objects = [_compile_object(source, base=base, identity=identity,
                                   executable=executable, flags=flags)
                   for source in sorted(SOURCE_DIRECTORY.glob("*.cpp"))]
        with tempfile.TemporaryDirectory(prefix="compile-", dir=directory) as temporary:
            candidate = Path(temporary) / binary.name
            try:
                subprocess.run([executable, *flags, *[str(path) for path, _ in objects],
                                "-lcrypto", "-o", str(candidate)], check=True, capture_output=True, text=True)
            except (OSError, subprocess.CalledProcessError) as error:
                detail = getattr(error, "stderr", "") or str(error)
                raise NativeBuildError("C++ engine compilation failed: " + detail) from error
            # A concurrent source edit must never acquire a stale build identity.
            if identity["source_sha256"] != source_identity():
                raise NativeBuildError("Native sources changed during compilation")
            record = dict(identity, binary_sha256=_digest(candidate),
                          build_seconds=perf_counter() - started,
                          compiled_translation_units=sum(compiled for _, compiled in objects),
                          reused_translation_units=sum(not compiled for _, compiled in objects))
            manifest = Path(temporary) / "manifest.json"
            manifest.write_text(json.dumps(record, indent=2) + "\n")
            manifest.replace(binary.with_suffix(".json"))
            candidate.replace(binary)
    return binary, record, True


def _decimal_integer(value: int | str) -> int:
    """Read arbitrarily wide decimal coefficients without changing Python globals."""
    if type(value) is int:
        return value
    if not isinstance(value, str) or not value or value == "-":
        raise NativeError("Malformed native integer coefficient")
    negative = value.startswith("-")
    digits = value[1:] if negative else value
    if not digits.isascii() or not digits.isdigit():
        raise NativeError("Malformed native integer coefficient")
    result = 0
    for start in range(0, len(digits), 9):
        chunk = digits[start:start + 9]
        result = result * 10**len(chunk) + int(chunk)
    return -result if negative else result


def _integer_decimal(value: int) -> str:
    """Write an unbounded integer without Python's decimal conversion limit."""
    if type(value) is not int:
        raise NativeError("Integer coefficient required")
    if value == 0:
        return "0"
    negative = value < 0
    value = abs(value)
    chunks = []
    while value:
        value, tail = divmod(value, 10**9)
        chunks.append(tail)
    return ("-" if negative else "") + str(chunks[-1]) + "".join(
        f"{chunk:09d}" for chunk in reversed(chunks[:-1]))


def _coefficient_pairs(value: Any) -> bool:
    return isinstance(value, (list, tuple)) and all(
        isinstance(pair, (list, tuple)) and len(pair) == 2 and type(pair[0]) is int
        for pair in value)


def _encode_coefficients(value: Any) -> Any:
    if isinstance(value, (list, tuple)):
        return [_encode_coefficients(x) for x in value]
    if isinstance(value, dict):
        return {key: ([[q, _integer_decimal(c) if type(c) is int else c] for q, c in item]
                      if key == "coefficient" and _coefficient_pairs(item)
                      else _encode_coefficients(item))
                for key, item in value.items()}
    return value


def _decode_coefficients(value: Any) -> Any:
    if isinstance(value, list):
        return [_decode_coefficients(x) for x in value]
    if isinstance(value, dict):
        return {key: ([[int(q), _decimal_integer(c)] for q, c in item] if key == "coefficient" and _coefficient_pairs(item)
                      else _decode_coefficients(item)) for key, item in value.items()}
    return value


def _json_arrays(value: Any) -> Any:
    """Normalize Python witness tuples as JSON arrays without encoding integers."""
    if isinstance(value, (list, tuple)):
        return [_json_arrays(item) for item in value]
    if isinstance(value, dict):
        return {key: _json_arrays(item) for key, item in value.items()}
    return value


def _strict_equal(expected: Any, actual: Any) -> bool:
    """Compare JSON values without bool/int/float coercion or decimal conversion."""
    if type(expected) is not type(actual):
        return False
    if isinstance(expected, dict):
        return (all(type(key) is str for key in expected)
                and all(type(key) is str for key in actual)
                and expected.keys() == actual.keys()
                and all(_strict_equal(item, actual[key]) for key, item in expected.items()))
    if isinstance(expected, list):
        return len(expected) == len(actual) and all(
            _strict_equal(left, right) for left, right in zip(expected, actual))
    return bool(expected == actual)


def reference_certificate(certificate: dict[str, Any]) -> dict[str, Any]:
    """Convert native provenance only; retain every mathematical/replay field."""
    result = deepcopy(certificate)
    result["schema"] = "zipper-collected-run-v1"
    result["scope"] = "literal_word_collected_spanning_procedure"
    result["source_sha256"] = result.pop("python_reference_source_sha256")
    result.pop("native_build", None)
    return result


class NativeSession:
    """One retained C++ process; each reduction has an isolated packed frontier."""

    def __init__(self, *, compiler: str = "g++", build_directory: Path | None = None,
                 flags: tuple[str, ...] = FLAGS, binary: Path | None = None,
                 cache_entries: int = 65536, cache_bytes: int = 256 * 1024**2) -> None:
        started = perf_counter()
        if any(type(x) is not int or x < 0 or x > 2**63 - 1
               for x in (cache_entries, cache_bytes)):
            raise NativeError("Cache budgets must be nonnegative signed 64-bit integers")
        self._lock = Lock()
        self._process: subprocess.Popen[str] | None = None
        self._stderr: TextIO | None = None
        catalog = Catalog.bundled()
        self._catalog_sha256 = catalog.sha256
        self._reference_source = collected._source()
        if binary is None:
            path, record, compiled = build_engine(compiler=compiler,
                                                 build_directory=build_directory, flags=flags)
        else:
            path = Path(binary).resolve()
            record, compiled = _manifest(path), False
        self.metadata = dict(record, binary=str(path), backend="cpp",
                             compiled_during_initialization=compiled,
                             arithmetic="arbitrary-precision integers over Z[q,q^-1]",
                             wrapper_source_sha256=_digest(Path(__file__)))
        self.last_call_seconds = 0.0
        self.last_runtime_info: dict[str, Any] = {}
        try:
            self._stderr = tempfile.TemporaryFile(mode="w+t")
            self._process = subprocess.Popen([str(path)], stdin=subprocess.PIPE,
                                             stdout=subprocess.PIPE, stderr=self._stderr,
                                             text=True, encoding="utf-8", bufsize=1)
            payload = files(__package__).joinpath("data/full-catalog.json").read_bytes()
            if hashlib.sha256(payload).hexdigest() != self._catalog_sha256:
                raise NativeError("Catalog changed during native initialization")
            data = json.loads(payload)
            ready = self._request({"op": "initialize", "abi": ABI, "catalog": data,
                                   "cache_entries": cache_entries, "cache_bytes": cache_bytes,
                                   "reference_sources": self._reference_source,
                                   "catalog_sha256": self._catalog_sha256})
            if ready != {"abi": ABI, "ready": True}:
                raise NativeError("Native engine initialization response differs")
        except BaseException:
            self.close()
            raise
        self.initialization_seconds = perf_counter() - started

    @property
    def pid(self) -> int | None:
        return None if self._process is None else self._process.pid

    def _request(self, value: dict[str, Any]) -> dict[str, Any]:
        process = self._process
        if process is None or process.stdin is None or process.stdout is None:
            raise NativeError("Native session is closed")
        try:
            process.stdin.write(json.dumps({key: _encode_coefficients(item) if key in ("work", "branch") else item
                                            for key, item in value.items()},
                                           separators=(",", ":"), allow_nan=False) + "\n")
            process.stdin.flush()
            line = process.stdout.readline()
            if not line:
                if self._stderr is not None:
                    self._stderr.flush()
                    self._stderr.seek(0)
                    detail = self._stderr.read()[-8000:]
                else:
                    detail = ""
                raise NativeError("Native engine terminated without a response: " + detail)
            result = json.loads(line)
            if not isinstance(result, dict):
                raise NativeError("Native engine returned a non-object response")
            if "protocol_error" in result:
                raise NativeError(str(result["protocol_error"]))
            return result
        except BaseException as error:
            # The request may already be in flight, even if its write/flush was
            # interrupted. Never let its eventual response answer another call.
            # The caller holds _lock, so cleanup must not acquire it again.
            with suppress(BaseException):
                self._close_unlocked(abort=True)
            if isinstance(error, (OSError, json.JSONDecodeError)):
                raise NativeError("Native engine protocol failed; session closed") from error
            raise

    def reduce(self, word: Iterable[int], *, n: int, limits: Limits | None = None,
               trace: bool = False, batch_size: int = 64, workers: int = 1) -> dict[str, Any]:
        if (type(n) is not int or n < 1 or type(trace) is not bool
                or type(batch_size) is not int or batch_size < 1
                or type(workers) is not int or workers != 1):
            raise NativeError("Positive rank and batch size, boolean trace and workers=1 required")
        if not isinstance(word, (list, tuple)) or any(
                type(x) is not int or not 1 <= abs(x) < n for x in word):
            raise NativeError("A literal signed-generator list or tuple in the stated rank is required")
        bounds = Limits() if limits is None else limits
        if not isinstance(bounds, Limits):
            raise NativeError("limits must be a Limits instance")
        if max(n, batch_size, *asdict(bounds).values()) > 2**63 - 1:
            raise NativeError("Rank, batch size and limits must fit signed 64-bit coordinates")
        started = perf_counter()
        with self._lock:
            body = self._request({"op": "reduce", "n": n, "word": list(word),
                                  "limits": asdict(bounds), "trace": trace,
                                  "batch_size": batch_size})
        self.last_runtime_info = body.pop("native_runtime", {})
        result = _decode_coefficients(body)
        result.update(n=n, input_word=list(word), coefficient_ring="Z[q,q^-1]",
                      trace_enabled=trace, limits=asdict(bounds),
                      schema="zipper-native-run-v1",
                      scope="literal_word_native_collected_spanning_procedure",
                      source_sha256=deepcopy(self.metadata["source_sha256"]),
                      python_reference_source_sha256=deepcopy(self._reference_source),
                      catalog_sha256=self._catalog_sha256,
                      native_build={k: deepcopy(self.metadata[k]) for k in
                                    ("abi", "binary_sha256", "compiler", "compiler_version", "flags",
                                     "wrapper_source_sha256", "link_libraries")})
        self.last_call_seconds = perf_counter() - started
        return result

    def _validate_certificate_identity(self, certificate: dict[str, Any]) -> None:
        """Authenticate native producer and independent reference identities."""
        if (not isinstance(certificate, dict)
                or certificate.get("schema") != "zipper-native-run-v1"
                or certificate.get("scope") != "literal_word_native_collected_spanning_procedure"
                or certificate.get("source_sha256") != self.metadata["source_sha256"]
                or certificate.get("python_reference_source_sha256") != collected._source()
                or certificate.get("catalog_sha256") != self._catalog_sha256
                or certificate.get("native_build") != {k: self.metadata[k] for k in
                    ("abi", "binary_sha256", "compiler", "compiler_version", "flags",
                     "wrapper_source_sha256", "link_libraries")}):
            raise NativeError("Native certificate source or catalog identity differs")
    def check(self, certificate: dict[str, Any], *,
              verifier: IndependentVerifier | None = None, cold: bool = False) -> bool:
        """Replay independently in Python; retain only Python-computed transitions.

        ``cold=True`` disables both Python cache layers for an audit. An explicit
        verifier allows a runner to budget and record verifier memory separately.
        """
        from .independent_replay import IndependentVerifier

        if type(cold) is not bool or (cold and verifier is not None):
            raise NativeError("cold must be boolean and cannot accompany a verifier")
        if not cold:
            if verifier is None:
                verifier = getattr(self, "_verifier", None)
                if verifier is None:
                    verifier = self._verifier = IndependentVerifier()
            return verifier.check_native(self, certificate)
        self._validate_certificate_identity(certificate)
        checker = collected.CollectedSession()
        checker._catalog.configure_rule_cache(max_bytes=0, max_entries=0)
        try:
            expected = checker.reduce(
                certificate["input_word"], n=certificate["n"],
                limits=Limits(**certificate["limits"]), trace=certificate["trace_enabled"],
                batch_size=certificate["execution"]["batch_size"],
            )
            actual = reference_certificate(certificate)
        except (KeyError, TypeError, ValueError) as error:
            raise NativeError("Malformed native certificate") from error
        if not _strict_equal(_json_arrays(expected), actual):
            raise NativeError("Native certificate differs from independent Python replay")
        return True

    def cache_info(self) -> dict[str, Any]:
        with self._lock:
            return self._request({"op": "cache_info"})

    def clear_cache(self) -> None:
        with self._lock:
            self._request({"op": "clear_cache"})

    def _close_unlocked(self, *, abort: bool) -> None:
        """Reap this session's child even if an interrupt arrives during detachment."""
        process = self._process
        stderr = None
        try:
            # Capture and detach under the same cleanup guard. An interrupt
            # while reading stderr must not strand an already detached child.
            stderr, self._stderr = self._stderr, None
            self._process = None
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
                    with suppress(OSError):
                        process.kill()
                    process.wait(timeout=5)
        except BaseException:
            self._process = None
            if process is not None:
                with suppress(BaseException):
                    process.kill()
                with suppress(BaseException):
                    process.wait(timeout=5)
            # If the first stderr read was interrupted, retain and close its
            # handle on this recovery path as well.
            if stderr is None:
                with suppress(BaseException):
                    stderr, self._stderr = self._stderr, None
            raise
        finally:
            if process is not None:
                for stream in (process.stdin, process.stdout):
                    if stream is not None:
                        with suppress(OSError, ValueError):
                            stream.close()
            if stderr is not None:
                with suppress(OSError, ValueError):
                    stderr.close()

    def close(self) -> None:
        with self._lock:
            self._close_unlocked(abort=False)

    def __enter__(self) -> NativeSession:
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()

    def __del__(self) -> None:
        if getattr(self, "_process", None) is not None:
            self.close()


def reduce_word(word: Iterable[int], *, n: int, limits: Limits | None = None,
                trace: bool = False, batch_size: int = 64) -> dict[str, Any]:
    with NativeSession() as session:
        return session.reduce(word, n=n, limits=limits, trace=trace, batch_size=batch_size)


def resume(certificate: dict[str, Any], *, limits: Limits | None = None) -> dict[str, Any]:
    """Independently check the saved prefix, then replay with enlarged bounds."""
    with NativeSession() as session:
        session.check(certificate)
        old = Limits(**certificate["limits"])
        following = old if limits is None else limits
        if any(getattr(following, key) < value for key, value in asdict(old).items()):
            raise NativeError("Continuation limits must not decrease")
        return session.reduce(certificate["input_word"], n=certificate["n"], limits=following,
                              trace=certificate["trace_enabled"],
                              batch_size=certificate["execution"]["batch_size"])
