"""Transport type preservation, build provenance, and interruption cleanup."""
from __future__ import annotations

from io import StringIO
from pathlib import Path
import subprocess
from threading import Lock
from typing import Any

import pytest

from zipper_algebra.zipping_algorithm import native_verifier as verifier


def test_decimal_without_changing_interpreter_limits() -> None:
    import sys
    limit = sys.get_int_max_str_digits()
    value = 10**10000 + 123
    assert verifier._decimal(value) == "1" + "0" * 9997 + "123"
    assert verifier._decimal(-value) == "-1" + "0" * 9997 + "123"
    assert verifier._decimal(0) == "0"
    assert sys.get_int_max_str_digits() == limit


@pytest.mark.parametrize("scalar", ["1", True, 1.0, None])
def test_coefficient_type_cannot_be_erased(scalar: Any) -> None:
    with pytest.raises(verifier.VerifierError, match="integer"):
        verifier._validate_candidate({"nested": [{"coefficient": [[0, scalar]]}]})


def test_wire_preserves_exact_coefficient() -> None:
    value = {"coefficient": [[-2, 10**5000]], "word": [1, -1]}
    verifier._validate_candidate(value)
    assert verifier._wire(value) == {"coefficient": [[-2, "1" + "0" * 5000]], "word": [1, -1]}
    assert type(value["coefficient"][0][1]) is int
    for malformed in ({"x": (1, 2)}, {1: "x"}):
        with pytest.raises(verifier.VerifierError):
            verifier._validate_candidate(malformed)
    assert not verifier._same({"a": [True]}, {"a": [1]})
    assert not verifier._same({"a": [1.0]}, {"a": [1]})


@pytest.fixture
def fake_build(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, list[list[str]]]:
    source = tmp_path / "sources"
    source.mkdir()
    for name in ("main.cpp", "other.cpp", "core.hpp"):
        (source / name).write_text(name)
    libraries = {}
    for name in verifier.LIBRARIES:
        path = tmp_path / (name + ".so")
        path.write_bytes(name.encode())
        libraries[name] = {"path": str(path), "sha256": verifier._hash(path)}
    identity = {"abi": verifier.ABI, "compiler": "fake", "compiler_version": "fake-v1",
                "flags": list(verifier.FLAGS), "link_libraries": list(verifier.LIBRARIES),
                "libraries": libraries}
    monkeypatch.setattr(verifier, "SOURCE_DIRECTORY", source)
    monkeypatch.setattr(verifier, "_compiler_identity", lambda *_: identity)
    calls: list[list[str]] = []

    def compile_(command: list[str]) -> None:
        calls.append(command)
        Path(command[-1]).write_bytes(repr(command).encode())

    monkeypatch.setattr(verifier, "_compile", compile_)
    return source, calls


def test_build_cache_tracks_headers_and_individual_translation_units(
    fake_build: tuple[Path, list[list[str]]], tmp_path: Path,
) -> None:
    source, calls = fake_build
    directory = tmp_path / "build"
    binary, metadata, compiled = verifier.build_verifier(build_directory=directory)
    assert compiled and metadata["compiled_translation_units"] == 2 and len(calls) == 3
    assert verifier.build_verifier(build_directory=directory)[2] is False
    assert len(calls) == 3
    (source / "other.cpp").write_text("new TU")
    _, metadata, _ = verifier.build_verifier(build_directory=directory)
    assert metadata["compiled_translation_units"] == 1 and metadata["reused_translation_units"] == 1
    (source / "core.hpp").write_text("new header")
    _, metadata, _ = verifier.build_verifier(build_directory=directory)
    assert metadata["compiled_translation_units"] == 2
    binary.write_bytes(b"corruption")
    # Restoring old source selects the old content-addressed executable.
    (source / "other.cpp").write_text("other.cpp")
    (source / "core.hpp").write_text("core.hpp")
    with pytest.raises(verifier.VerifierBuildError, match="identity differs"):
        verifier.build_verifier(build_directory=directory)


def test_build_refuses_changed_source_during_compilation(
    fake_build: tuple[Path, list[list[str]]], tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source, _ = fake_build
    original = verifier._compile

    def changing(command: list[str]) -> None:
        original(command)
        (source / "core.hpp").write_text("modified during build")

    monkeypatch.setattr(verifier, "_compile", changing)
    with pytest.raises(verifier.VerifierBuildError, match="changed during compilation"):
        verifier.build_verifier(build_directory=tmp_path / "build")
    assert not list((tmp_path / "build").rglob("*.o"))


class FakeProcess:
    def __init__(self, response: str) -> None:
        self.stdin = StringIO()
        self.stdout = StringIO(response)
        self.killed = False
        self.waited = False

    def kill(self) -> None:
        self.killed = True

    def wait(self, **_: Any) -> int:
        self.waited = True
        return 0


def session_with_process(process: Any) -> verifier.CppVerifier:
    session = object.__new__(verifier.CppVerifier)
    session._lock = Lock()
    session._allocator = None
    session._process = process
    session._stderr = StringIO()
    return session


@pytest.mark.parametrize("response", ["", "[]\n", "malformed\n"])
def test_failed_protocol_closes_and_reaps(response: str) -> None:
    process = FakeProcess(response)
    session = session_with_process(process)
    stderr = session._stderr
    with pytest.raises(ValueError):
        session._request({"op": "verify"})
    assert process.killed and process.waited
    assert process.stdin.closed and process.stdout.closed
    assert stderr is not None and stderr.closed
    with pytest.raises(verifier.VerifierError, match="closed"):
        session._request({"op": "verify"})
    session.close()


def test_certificate_rejection_does_not_poison_session() -> None:
    process = FakeProcess('{"ok":false,"error":"bad certificate"}\n{"ok":true}\n')
    session = session_with_process(process)
    try:
        assert session._request({})["ok"] is False
        assert session._request({})["ok"] is True
        assert not process.killed
    finally:
        session.close()
    assert process.waited


def test_interrupted_read_reaps_actual_subprocess() -> None:
    process = subprocess.Popen(["sleep", "60"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    session = session_with_process(process)

    class InterruptedReader:
        closed = False

        def readline(self) -> str:
            raise KeyboardInterrupt

        def close(self) -> None:
            self.closed = True

    original_stdout = process.stdout
    reader = InterruptedReader()
    process.stdout = reader  # type: ignore[assignment]
    try:
        with pytest.raises(KeyboardInterrupt):
            session._request({})
        assert process.poll() is not None and reader.closed
        assert process.stdin is not None and process.stdin.closed
    finally:
        if original_stdout is not None:
            original_stdout.close()
        if process.poll() is None:
            process.kill()
            process.wait()


def test_cleanup_interrupt_before_stderr_assignment() -> None:
    class InterruptedSession(verifier.CppVerifier):
        def __getattribute__(self, name: str) -> Any:
            if name == "_stderr":
                raise KeyboardInterrupt
            return super().__getattribute__(name)

    process = FakeProcess("")
    session = object.__new__(InterruptedSession)
    session._lock = Lock()
    session._process = process  # type: ignore[assignment]
    stderr = StringIO()
    session._stderr = stderr
    with pytest.raises(KeyboardInterrupt):
        session.close()
    assert process.killed and process.waited and stderr.closed
    assert process.stdin.closed and process.stdout.closed


def test_cleanup_timeout_forces_reaping() -> None:
    class TimeoutProcess(FakeProcess):
        def wait(self, **_: Any) -> int:
            if not self.killed:
                raise subprocess.TimeoutExpired("verifier", 5)
            return super().wait()

    process = TimeoutProcess("")
    session = session_with_process(process)
    session.close()
    assert process.killed and process.waited


@pytest.mark.parametrize("budget", [-1, True, 1.5, 2**63])
def test_constructor_rejects_invalid_cache_budgets(budget: Any) -> None:
    with pytest.raises(verifier.VerifierError, match="cache budgets"):
        verifier.CppVerifier(producer_metadata={}, cache_bytes=budget)


def test_public_operations_use_independent_protocol(monkeypatch: pytest.MonkeyPatch) -> None:
    process = FakeProcess('{"ok":true}\n{"states":2}\n{"cleared":true}\n'
                          '{"payload":1}\n{"ok":false,"error":"diagnostic rejected"}\n'
                          '{"ok":false,"error":"certificate rejected"}\n{"cleared":false}\n')
    session = session_with_process(process)
    monkeypatch.setattr(session, "_assert_source", lambda: None)
    try:
        assert session.check({"coefficient": [[0, 1]]})
        assert session.cache_info() == {"states": 2}
        session.clear_cache()
        assert session.diagnostic("replay", n=2) == {"payload": 1}
        with pytest.raises(verifier.VerifierError, match="diagnostic rejected"):
            session.diagnostic("transition")
        with pytest.raises(verifier.VerifierError, match="certificate rejected"):
            session.check({})
        with pytest.raises(verifier.VerifierError, match="reset"):
            session.clear_cache()
        with pytest.raises(verifier.VerifierError, match="Unknown"):
            session.diagnostic("producer")
        with pytest.raises(verifier.VerifierError, match="object"):
            session.check([])  # type: ignore[arg-type]
        assert '"coefficient":[[0,"1"]]' in process.stdin.getvalue()
    finally:
        session.close()
    assert session.pid is None


def test_missing_pipes_abort_process() -> None:
    process = FakeProcess("")
    process.stdin.close()
    process.stdin = None  # type: ignore[assignment]
    session = session_with_process(process)
    with pytest.raises(verifier.VerifierError, match="pipes"):
        session._request({})
    assert process.killed and process.waited and process.stdout.closed


@pytest.mark.parametrize("ready", ['{"ready":true,"abi":"zipper-independent-verifier-v1"}\n', '{}\n'])
def test_initialization_protocol_and_cleanup(
    ready: str, monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    process = FakeProcess(ready)
    process.pid = 123  # type: ignore[attr-defined]
    monkeypatch.setattr(verifier, "reference_identity", lambda: {"reference.py": "hash"})
    monkeypatch.setattr(verifier, "_producer_identity", lambda *_: {"schema": "trusted"})
    monkeypatch.setattr(verifier, "build_verifier", lambda **_: (tmp_path / "binary", {}, False))
    monkeypatch.setattr(verifier.subprocess, "Popen", lambda *_, **__: process)
    if ready == '{}\n':
        with pytest.raises(verifier.VerifierError, match="initialization"):
            verifier.CppVerifier(producer_metadata={})
        assert process.waited and process.stdin.closed and process.stdout.closed
    else:
        with verifier.CppVerifier(producer_metadata={}) as session:
            assert session.pid == 123
            assert session.initialization_seconds >= 0
            request = process.stdin.getvalue()
            assert '"catalog_text"' in request and '"native_identity":{"schema":"trusted"}' in request
        assert process.waited


def test_source_change_guard_fails_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    session = session_with_process(FakeProcess(""))
    session.metadata = {"source_sha256": {"main.cpp": "original"}}
    monkeypatch.setattr(verifier, "source_identity", lambda: {"main.cpp": "changed"})
    try:
        with pytest.raises(verifier.VerifierError, match="sources changed"):
            session.check({})
        assert session._process is not None
        assert session._process.stdin is not None
        assert session._process.stdin.getvalue() == ""  # type: ignore[attr-defined]
    finally:
        session.close()


def test_explicit_authenticated_binary_and_destructor(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    process = FakeProcess('{"ready":true,"abi":"zipper-independent-verifier-v1"}\n')
    monkeypatch.setattr(verifier, "reference_identity", lambda: {})
    monkeypatch.setattr(verifier, "_producer_identity", lambda *_: {})
    binary = tmp_path / "verified"
    authenticated = []

    def authenticate(path: Path) -> dict[str, Any]:
        authenticated.append(path)
        return {}

    monkeypatch.setattr(verifier, "_authenticated_binary", authenticate)
    monkeypatch.setattr(verifier.subprocess, "Popen", lambda *_, **__: process)
    session = verifier.CppVerifier(producer_metadata={}, binary=binary)
    assert authenticated == [binary]
    session.__del__()
    assert process.waited and process.stdin.closed
    session.__del__()


def test_interrupted_cleanup_without_process() -> None:
    class InterruptedSession(verifier.CppVerifier):
        def __getattribute__(self, name: str) -> Any:
            if name == "_stderr":
                raise KeyboardInterrupt
            return super().__getattribute__(name)

    session = object.__new__(InterruptedSession)
    session._lock = Lock()
    session._process = None
    session._stderr = StringIO()
    stderr = session.__dict__["_stderr"]
    with pytest.raises(KeyboardInterrupt):
        session.close()
    assert stderr.closed


def test_actual_native_certificate_independent_smoke() -> None:
    from zipper_algebra.zipping_algorithm.native import NativeSession
    with NativeSession() as producer:
        with verifier.CppVerifier(producer_metadata=producer.metadata) as checker:
            for trace in (False, True):
                certificate = producer.reduce([1, 2, 1, -2], n=3, trace=trace)
                assert checker.check(certificate)
            checker.clear_cache()
            assert isinstance(checker.cache_info(), dict)


def test_unserializable_request_preserves_transport() -> None:
    process = FakeProcess('{"ok":true}\n')
    session = session_with_process(process)
    try:
        with pytest.raises(verifier.VerifierError, match="encoded"):
            session._request({"value": object()})
        assert not process.killed and not process.waited
        assert process.stdin.getvalue() == ""
        assert session._request({}) == {"ok": True}
    finally:
        session.close()


def test_malformed_coordinate_preserves_healthy_verifier() -> None:
    from copy import deepcopy
    import sys

    from zipper_algebra.zipping_algorithm.native import NativeSession

    with NativeSession() as producer:
        certificate = producer.reduce([1, 2, 1, -2], n=3)
        with verifier.CppVerifier(producer_metadata=producer.metadata) as checker:
            pid = checker.pid
            bad_values = [float("nan"), float("inf"), float("-inf")]
            if sys.get_int_max_str_digits():
                bad_values.append(10 ** (sys.get_int_max_str_digits() + 1))
            for coordinate in bad_values:
                malformed = deepcopy(certificate)
                malformed["n"] = coordinate
                with pytest.raises(verifier.VerifierError, match="encoded"):
                    checker.check(malformed)
                assert checker.pid == pid
                assert checker.check(certificate)


def test_allocator_identity_errors(tmp_path: Path) -> None:
    with pytest.raises(verifier.VerifierError, match="unavailable"):
        verifier._allocator_identity(tmp_path / "missing")
    invalid = tmp_path / "with space.so"
    invalid.write_bytes(b"library")
    with pytest.raises(verifier.VerifierError, match="separators"):
        verifier._allocator_identity(invalid)


def test_allocator_mapping_authentication(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import os
    library = tmp_path / "allocator.so"
    library.write_bytes(b"library")
    identity = verifier._allocator_identity(library)
    stat = library.stat()
    device = f"{os.major(stat.st_dev):02x}:{os.minor(stat.st_dev):02x}"
    monkeypatch.setattr(Path, "read_text", lambda *_: f"short\n0-1 r-xp 0 00:00 0\n0-1 r-xp 0 {device} {stat.st_ino} {library}\n")
    assert verifier._allocator_mapped(123, identity)
    monkeypatch.setattr(Path, "read_text", lambda *_: "0-1 r-xp 0 00:00 0\n")
    assert not verifier._allocator_mapped(123, identity)
    monkeypatch.setattr(Path, "read_text", lambda *_: "0-1 r-xp 0 invalid invalid /library\n")
    with pytest.raises(verifier.VerifierError, match="mapping"):
        verifier._allocator_mapped(123, identity)
    library.unlink()
    with pytest.raises(verifier.VerifierError, match="mapping"):
        verifier._allocator_mapped(123, identity)


@pytest.mark.parametrize("mode", ["unmapped", "changed", "valid", "default"])
def test_allocator_child_environment_and_mapping_failure(
    mode: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    import os
    library = tmp_path / "allocator.so"
    library.write_bytes(b"library")
    process = FakeProcess('{"ready":true,"abi":"zipper-independent-verifier-v1"}\n')
    process.pid = 123  # type: ignore[attr-defined]
    monkeypatch.setenv("LD_PRELOAD", "original-parent-value")
    monkeypatch.setattr(verifier, "reference_identity", lambda: {})
    monkeypatch.setattr(verifier, "_producer_identity", lambda *_: {})
    monkeypatch.setattr(verifier, "build_verifier", lambda **_: (tmp_path / "binary", {}, False))
    monkeypatch.setattr(verifier, "_allocator_mapped", lambda *_: mode != "unmapped")
    environments = []

    def popen(*_: Any, **kwargs: Any) -> FakeProcess:
        environments.append(kwargs["env"])
        if mode == "changed":
            library.write_bytes(b"changed after launch")
        return process

    monkeypatch.setattr(verifier.subprocess, "Popen", popen)
    if mode in ("unmapped", "changed"):
        with pytest.raises(verifier.VerifierError, match="not authenticated"):
            verifier.CppVerifier(producer_metadata={}, allocator_library=library)
        assert process.waited and process.stdin.closed and process.stdout.closed
    else:
        with verifier.CppVerifier(producer_metadata={}, allocator_library=None if mode == "default" else library) as checker:
            if mode == "default":
                assert checker.metadata["allocator_library"] is None
                assert environments == [None]
            else:
                assert checker.metadata["allocator_library"] == dict(verifier._allocator_identity(library), mapped_verified=True)
                assert environments[0]["LD_PRELOAD"] == str(library)
    assert os.environ["LD_PRELOAD"] == "original-parent-value"


@pytest.mark.parametrize("remove", [False, True])
def test_allocator_guard_closes_on_drift(tmp_path: Path, remove: bool) -> None:
    library = tmp_path / "allocator.so"
    library.write_bytes(b"library")
    process = FakeProcess("")
    session = session_with_process(process)
    session._allocator = verifier._allocator_identity(library)
    if remove:
        library.unlink()
    else:
        library.write_bytes(b"modified")
    with pytest.raises(verifier.VerifierError, match="changed|unavailable"):
        session._assert_source()
    assert process.killed and process.waited


def test_actual_tcmalloc_verifier() -> None:
    import os
    from zipper_algebra.zipping_algorithm.native import NativeSession
    library = Path("/lib64/libtcmalloc_minimal.so.4")
    if not library.exists():
        pytest.skip("tcmalloc is not installed")
    environment = os.environ.copy()
    with NativeSession() as producer:
        certificate = producer.reduce([1, 2, 1, -2], n=3, trace=True)
        with verifier.CppVerifier(producer_metadata=producer.metadata, allocator_library=library) as checker:
            assert checker.check(certificate)
            assert checker.metadata["allocator_library"] == dict(verifier._allocator_identity(library), mapped_verified=True)
    assert os.environ == environment
