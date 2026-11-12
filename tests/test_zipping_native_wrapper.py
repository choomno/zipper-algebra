"""Boundary and build-contract tests independent of native numerical tests."""
from copy import deepcopy
import subprocess

import pytest

from zipper_algebra.zipping_algorithm import native


def test_decimal_arbitrary_width():
    digits = "9" * 10000
    assert native._decimal_integer(digits) == 10**10000 - 1
    assert native._decimal_integer("-" + digits) == 1 - 10**10000
    assert native._decimal_integer(17) == 17
    for invalid in ("", "-", "+1", " 2", "1.0", "１２", True, None):
        with pytest.raises(native.NativeError, match="Malformed"):
            native._decimal_integer(invalid)


def test_decode_coefficients_preserves_other_strings():
    value = {"coefficient": [[0, "999999999999999999999999999999"]],
             "child": [{"coefficient": [[-4, "-3"]], "name": "123"}]}
    result = native._decode_coefficients(value)
    assert result["coefficient"] == [[0, 999999999999999999999999999999]]
    assert result["child"][0] == {"coefficient": [[-4, -3]], "name": "123"}
    assert value["coefficient"][0][1] == "999999999999999999999999999999"


def test_reference_projection_preserves_calculation():
    original = {"schema": "zipper-native-run-v1", "scope": "native",
                "source_sha256": {"a.cpp": "abc"},
                "python_reference_source_sha256": {"a.py": "def"},
                "native_build": {"abi": native.ABI},
                "output": [{"word": [1, -1], "coefficient": [[0, 2**300]]}],
                "execution": {"batch_size": 7}, "transitions": [{"evidence": 3}]}
    saved = deepcopy(original)
    projected = native.reference_certificate(original)
    assert projected["schema"] == "zipper-collected-run-v1"
    assert projected["source_sha256"] == {"a.py": "def"}
    assert projected["output"] == original["output"]
    assert projected["transitions"] == original["transitions"]
    assert projected["execution"] == original["execution"]
    assert original == saved


def test_source_identity_requires_entrypoint(tmp_path):
    with pytest.raises(native.NativeBuildError, match="source checkout"):
        native.source_identity(tmp_path)
    (tmp_path / "main.cpp").write_text("int main() {}")
    assert list(native.source_identity(tmp_path)) == ["main.cpp"]


def test_build_missing_compiler(monkeypatch):
    monkeypatch.setattr(native.shutil, "which", lambda _: None)
    with pytest.raises(native.NativeBuildError, match="compiler is unavailable"):
        native.build_engine()


def test_build_compiler_identity_failure(monkeypatch):
    monkeypatch.setattr(native.shutil, "which", lambda _: "/bin/g++")
    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(1, args[0])
    monkeypatch.setattr(native.subprocess, "run", fail)
    with pytest.raises(native.NativeBuildError, match="identify"):
        native.build_engine()


def test_manifest_rejects_missing_file(tmp_path):
    with pytest.raises(native.NativeBuildError, match="matching provenance"):
        native._manifest(tmp_path / "missing")


def test_closed_protocol():
    session = object.__new__(native.NativeSession)
    session._process = None
    with pytest.raises(native.NativeError, match="closed"):
        session._request({})


def test_manifest_rejects_tampered_binary(tmp_path, monkeypatch):
    import json
    executable = tmp_path / "engine"
    executable.write_bytes(b"binary")
    monkeypatch.setattr(native, "source_identity", lambda: {"x": "y"})
    record = {"abi": native.ABI, "source_sha256": {"x": "y"},
              "binary_sha256": native._digest(executable), "flags": [],
              "compiler": "g++", "compiler_version": "test", "link_libraries": ["crypto"]}
    executable.with_suffix(".json").write_text(json.dumps(record))
    assert native._manifest(executable) == record
    with pytest.raises(native.NativeBuildError, match="provenance differs"):
        native._manifest(executable, {"flags": ["different"]})
    executable.write_bytes(b"tampered")
    with pytest.raises(native.NativeBuildError, match="provenance differs"):
        native._manifest(executable)


def test_reduce_input_validation_without_process():
    session = object.__new__(native.NativeSession)
    session._process = None
    for word, options in (([True], {"n": 2}), ([2], {"n": 2}),
                          ([], {"n": True}), ([], {"n": 2, "workers": 2}),
                          ([], {"n": 2, "trace": 1}), ([], {"n": 2, "batch_size": 0}),
                          (iter([]), {"n": 2}), ([], {"n": 2, "limits": {}}),
                          ([], {"n": 2**63})):
        with pytest.raises(native.NativeError):
            session.reduce(word, **options)


def test_decimal_round_trip_and_catalog_ast():
    for number in (0, 1, -1, 2**32770 + 123, -(2**32770 + 123)):
        assert native._decimal_integer(native._integer_decimal(number)) == number
    with pytest.raises(native.NativeError, match="Integer coefficient"):
        native._integer_decimal(True)
    scalar_ast = {"coefficient": {"integer": -1, "q_power": "i"}}
    assert native._encode_coefficients(scalar_ast) == scalar_ast
    assert native._decode_coefficients(scalar_ast) == scalar_ast


def test_protocol_encodes_only_work_coefficients():
    import io
    import json
    from types import SimpleNamespace
    session = object.__new__(native.NativeSession)
    session._process = SimpleNamespace(stdin=io.StringIO(), stdout=io.StringIO('{"ok":true}\n'))
    session._stderr = None
    scalar = 2**32770
    request = {"op": "transition", "work": {"coefficient": [[0, scalar]]},
               "catalog": {"coefficient": [[0, 1]]}}
    assert session._request(request) == {"ok": True}
    serialized = json.loads(session._process.stdin.getvalue())
    assert native._decimal_integer(serialized["work"]["coefficient"][0][1]) == scalar
    assert serialized["catalog"]["coefficient"] == [[0, 1]]
    session._process = None


@pytest.mark.parametrize("line", ['[]\n', '{"protocol_error":"bad input"}\n', 'not JSON\n'])
def test_protocol_rejects_bad_responses(line):
    import io
    from types import SimpleNamespace
    session = object.__new__(native.NativeSession)
    session._process = SimpleNamespace(stdin=io.StringIO(), stdout=io.StringIO(line))
    session._stderr = None
    try:
        with pytest.raises(native.NativeError):
            session._request({})
    finally:
        session._process = None


def test_object_cache_reuses_and_invalidates_one_translation_unit(tmp_path, monkeypatch):
    from pathlib import Path
    source_root = tmp_path / "source"
    source_root.mkdir()
    first = source_root / "one.cpp"
    second = source_root / "two.cpp"
    header = source_root / "engine.hpp"
    first.write_text("int one() {return 1;}")
    second.write_text("int two() {return 2;}")
    header.write_text("// declarations")
    monkeypatch.setattr(native, "SOURCE_DIRECTORY", source_root)
    calls = []
    def compile_mock(command, **kwargs):
        calls.append(command)
        Path(command[-1]).write_bytes(Path(command[command.index("-c") + 1]).read_bytes())
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(native.subprocess, "run", compile_mock)
    def identity():
        return {"source_sha256": {p.name: native._digest(p) for p in source_root.iterdir()},
                "compiler": "/fake/g++", "compiler_version": "test", "flags": ["-O3"]}
    def build(source):
        return native._compile_object(source, base=tmp_path / "cache", identity=identity(),
                                      executable="g++", flags=("-O3",))
    one, compiled = build(first)
    assert compiled
    two, compiled = build(second)
    assert compiled
    assert build(first) == (one, False)
    assert len(calls) == 2
    second.write_text("int two() {return 22;}")
    assert build(first) == (one, False)
    changed, compiled = build(second)
    assert compiled and changed != two and len(calls) == 3
    header.write_text("// changed declarations")
    changed, compiled = build(first)
    assert compiled and changed != one
    changed.write_bytes(b"tampered")
    with pytest.raises(native.NativeBuildError, match="provenance differs"):
        build(first)


def test_strict_certificate_comparison_rejects_scalar_type_tampering():
    expected = {"output": [{"slider": ["slider", 2, 1, False, False, False],
                             "coefficient": [[0, 2**32770]]}],
                "transitions": [{"state": {"world_phi": True, "depth": 2}}]}
    assert native._strict_equal(expected, deepcopy(expected))
    for value in (0, 0.0):
        changed = deepcopy(expected)
        changed["output"][0]["slider"][3] = value
        assert not native._strict_equal(expected, changed)
    changed = deepcopy(expected)
    changed["transitions"][0]["state"]["depth"] = 2.0
    assert not native._strict_equal(expected, changed)
    changed = deepcopy(expected)
    changed["transitions"][0]["state"]["world_phi"] = 1
    assert not native._strict_equal(expected, changed)
    changed = deepcopy(expected)
    changed["output"] = tuple(changed["output"])
    assert not native._strict_equal(expected, changed)
    assert not native._strict_equal({"x": 1}, {"x": 1, "y": 2})
    assert not native._strict_equal([1], [1, 2])
    assert not native._strict_equal({1: 1}, {True: 1})


def test_expected_python_witness_tuples_become_json_arrays():
    expected = {"witness": {"word": (1, -1), "rows": [(2**32770, True)]}}
    actual = {"witness": {"word": [1, -1], "rows": [[2**32770, True]]}}
    assert native._strict_equal(native._json_arrays(expected), actual)
    assert isinstance(expected["witness"]["word"], tuple)


def test_check_rejects_nested_boolean_and_float_tampering():
    expected = native.collected.CollectedSession().reduce([1], n=2, trace=True)
    session = object.__new__(native.NativeSession)
    session._process = None
    session.metadata = {"source_sha256": {"engine.cpp": "source"}, "abi": native.ABI,
                        "binary_sha256": "binary", "compiler": "g++",
                        "compiler_version": "test", "flags": [],
                        "wrapper_source_sha256": "wrapper", "link_libraries": ["crypto"]}
    session._catalog_sha256 = expected["catalog_sha256"]
    certificate = native._json_arrays(expected)
    certificate.update(schema="zipper-native-run-v1",
                       scope="literal_word_native_collected_spanning_procedure",
                       source_sha256=session.metadata["source_sha256"],
                       python_reference_source_sha256=expected["source_sha256"],
                       native_build={key: value for key, value in session.metadata.items()
                                     if key != "source_sha256"})
    assert session.check(certificate)
    changed = deepcopy(certificate)
    changed["metrics"]["steps"] = float(changed["metrics"]["steps"])
    with pytest.raises(native.NativeError, match="differs from independent"):
        session.check(changed)
    changed = deepcopy(certificate)
    slider = next(token for item in changed["output"] for token in item["word"]
                  if isinstance(token, list))
    assert slider[3] is False
    slider[3] = 0
    with pytest.raises(native.NativeError, match="differs from independent"):
        session.check(changed)


@pytest.mark.parametrize("stage,error", [("read", KeyboardInterrupt), ("read", OSError),
                                         ("write", KeyboardInterrupt), ("write", BrokenPipeError),
                                         ("flush", KeyboardInterrupt), ("flush", OSError)])
def test_interrupted_protocol_permanently_closes_session(stage, error):
    import io
    from threading import Lock
    class Input(io.StringIO):
        def write(self, value):
            if stage == "write":
                raise error("interrupted write")
            return super().write(value)
        def flush(self):
            if stage == "flush":
                raise error("interrupted flush")
            return super().flush()
    class Output(io.StringIO):
        def readline(self, *args):
            if stage == "read":
                raise error("interrupted read")
            return super().readline(*args)
    class Process:
        def __init__(self):
            self.stdin = Input()
            self.stdout = Output('{"stale":true}\n')
            self.killed = self.waited = False
        def kill(self):
            self.killed = True
        def wait(self, timeout):
            self.waited = True
            return -9
    session = object.__new__(native.NativeSession)
    session._lock = Lock()
    process = session._process = Process()
    stderr = session._stderr = io.StringIO()
    expected_error = KeyboardInterrupt if error is KeyboardInterrupt else native.NativeError
    # This is the same lock context as reduce; cleanup must not reacquire it.
    with session._lock:
        with pytest.raises(expected_error):
            session._request({"op": "reduce", "n": 3})
    assert session.pid is None and process.killed and process.waited
    assert process.stdin.closed and process.stdout.closed and stderr.closed
    with pytest.raises(native.NativeError, match="closed"):
        session._request({"op": "reduce", "n": 1})
    session.close()


def test_close_reaps_child_after_broken_pipe():
    import io
    from threading import Lock
    class Input(io.StringIO):
        def close(self):
            super().close()
            raise BrokenPipeError("child exited")
    class Process:
        def __init__(self):
            self.stdin = Input()
            self.stdout = io.StringIO()
            self.waited = False
        def wait(self, timeout):
            self.waited = True
            return 0
    session = object.__new__(native.NativeSession)
    session._lock = Lock()
    process = session._process = Process()
    stderr = session._stderr = io.StringIO()
    session.close()
    assert session.pid is None and process.waited
    assert process.stdin.closed and process.stdout.closed and stderr.closed
    session.close()


def test_object_cache_recovers_interrupted_publication(tmp_path, monkeypatch):
    from pathlib import Path
    source_root = tmp_path / "source"
    source_root.mkdir()
    source = source_root / "one.cpp"
    source.write_text("int one() {return 1;}")
    monkeypatch.setattr(native, "SOURCE_DIRECTORY", source_root)
    identity = {"source_sha256": {source.name: native._digest(source)},
                "compiler": "g++", "compiler_version": "test", "flags": []}
    def compile_mock(command, **kwargs):
        Path(command[-1]).write_bytes(b"object")
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(native.subprocess, "run", compile_mock)
    def build():
        return native._compile_object(source, base=tmp_path / "cache", identity=identity,
                                      executable="g++", flags=())
    replace = Path.replace
    def interrupt_publish(path, target):
        if path.suffix == ".o":
            raise KeyboardInterrupt("interrupted object publication")
        return replace(path, target)
    monkeypatch.setattr(Path, "replace", interrupt_publish)
    with pytest.raises(KeyboardInterrupt):
        build()
    assert list((tmp_path / "cache").rglob("*.json"))
    assert not list((tmp_path / "cache").rglob("*.o"))
    monkeypatch.setattr(Path, "replace", replace)
    target, compiled = build()
    assert compiled and target.read_bytes() == b"object"
    # Also recover entries written by the former binary-before-manifest order.
    target.with_suffix(".json").unlink()
    assert build() == (target, True)


@pytest.fixture
def fake_compiler(tmp_path, monkeypatch):
    from pathlib import Path
    source_root = tmp_path / "checkout" / "native" / "zipper_engine"
    source_root.mkdir(parents=True)
    (source_root / "main.cpp").write_text("int main() {return 0;}")
    (source_root / "engine.hpp").write_text("// declarations")
    monkeypatch.setattr(native, "SOURCE_DIRECTORY", source_root)
    monkeypatch.setattr(native, "_CHECKOUT_SOURCES", source_root)
    monkeypatch.setattr(native, "source_identity", lambda: {
        path.name: native._digest(path) for path in source_root.iterdir()})
    monkeypatch.setattr(native.shutil, "which", lambda _: "/fake/g++")
    state = {"compile_fail": False, "link_fail": False, "mutate_object": False,
             "mutate_link": False, "calls": [], "source": source_root}
    def run(command, **kwargs):
        state["calls"].append(command)
        if "--version" in command:
            return subprocess.CompletedProcess(command, 0, stdout="fake C++ compiler")
        compiling = "-c" in command
        if state["compile_fail" if compiling else "link_fail"]:
            raise subprocess.CalledProcessError(1, command, stderr="injected compiler failure")
        Path(command[-1]).write_bytes(b"object" if compiling else b"executable")
        if state["mutate_object" if compiling else "mutate_link"]:
            (source_root / "engine.hpp").write_text("// changed during compilation")
        return subprocess.CompletedProcess(command, 0)
    monkeypatch.setattr(native.subprocess, "run", run)
    return state


def test_build_cache_recovery_and_installed_location(fake_compiler, tmp_path, monkeypatch):
    from pathlib import Path
    path, manifest, built = native.build_engine(build_directory=tmp_path / "explicit")
    assert built and manifest["compiled_translation_units"] == 1
    assert native.build_engine(build_directory=tmp_path / "explicit") == (path, manifest, False)
    path.with_suffix(".json").unlink()
    same, rebuilt, built = native.build_engine(build_directory=tmp_path / "explicit")
    assert same == path and built and rebuilt["reused_translation_units"] == 1
    checkout, _, _ = native.build_engine()
    assert ".local-tools" in str(checkout)
    monkeypatch.setattr(native, "_CHECKOUT_SOURCES", Path("/different/checkout"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    installed, _, _ = native.build_engine()
    assert installed.is_relative_to(tmp_path / "xdg" / "zipper-algebra" / "native-engine")


@pytest.mark.parametrize("failure,message", [
    ("compile_fail", "translation-unit compilation failed"),
    ("link_fail", "engine compilation failed"),
    ("mutate_object", "changed during object compilation"),
    ("mutate_link", "sources changed during compilation"),
])
def test_build_rejects_failed_or_moving_sources(fake_compiler, tmp_path, failure, message):
    fake_compiler[failure] = True
    with pytest.raises(native.NativeBuildError, match=message):
        native.build_engine(build_directory=tmp_path / "cache")
    assert not list((tmp_path / "cache").rglob("zipper-engine"))


def test_corrupted_object_manifest_is_rejected(fake_compiler, tmp_path):
    native.build_engine(build_directory=tmp_path / "cache")
    object_manifest = next((tmp_path / "cache" / "objects").rglob("*.json"))
    object_manifest.write_text("invalid JSON")
    binary = next((tmp_path / "cache").glob("*/zipper-engine"))
    binary.unlink()
    with pytest.raises(native.NativeBuildError, match="object cache manifest is invalid"):
        native.build_engine(build_directory=tmp_path / "cache")


@pytest.mark.parametrize("bad", [[], {}, {"abi": 1}])
def test_manifest_rejects_incomplete_build_metadata(tmp_path, bad):
    import json
    path = tmp_path / "engine"
    path.with_suffix(".json").write_text(json.dumps(bad))
    with pytest.raises(native.NativeBuildError, match="metadata is malformed"):
        native._manifest(path)


@pytest.fixture
def fake_native_process(monkeypatch, tmp_path):
    import io
    import json
    from pathlib import Path
    expected = native._json_arrays(native.collected.CollectedSession().reduce([], n=1))
    metadata = {"abi": native.ABI, "source_sha256": {"engine.cpp": "native-source"},
                "binary_sha256": "binary", "compiler": "g++", "compiler_version": "test",
                "flags": [], "link_libraries": ["crypto"]}
    state = {"replies": [{"abi": native.ABI, "ready": True}, expected,
                         {"hits": 1}, {"hits": 0}], "children": []}
    class Process:
        def __init__(self):
            self.pid = 1234
            self.stdin = io.StringIO()
            self.stdout = io.StringIO("".join(json.dumps(reply) + "\n" for reply in state["replies"]))
            self.killed = self.waited = False
        def kill(self):
            self.killed = True
        def wait(self, timeout):
            self.waited = True
            return 0
    def popen(*args, **kwargs):
        process = Process()
        state["children"].append(process)
        return process
    monkeypatch.setattr(native.subprocess, "Popen", popen)
    monkeypatch.setattr(native, "build_engine", lambda **kwargs: (tmp_path / "engine", metadata, False))
    monkeypatch.setattr(native, "_manifest", lambda _: metadata)
    state["binary"] = Path(tmp_path / "engine")
    return state


@pytest.mark.parametrize("explicit_binary", [False, True])
def test_session_boundary_lifecycle_and_check(fake_native_process, explicit_binary):
    options = {"binary": fake_native_process["binary"]} if explicit_binary else {}
    with native.NativeSession(**options) as session:
        assert session.pid == 1234
        result = session.reduce([], n=1)
        assert result["input_word"] == [] and session.check(result)
        assert session.cache_info() == {"hits": 1}
        session.clear_cache()
        malformed = deepcopy(result)
        del malformed["input_word"]
        with pytest.raises(native.NativeError, match="Malformed native certificate"):
            session.check(malformed)
        altered = deepcopy(result)
        altered["native_build"]["flags"] = ["changed"]
        with pytest.raises(native.NativeError, match="identity differs"):
            session.check(altered)
    assert session.pid is None
    assert fake_native_process["children"][0].waited
    assert session.initialization_seconds >= 0 and session.last_call_seconds >= 0


@pytest.mark.parametrize("budget", [-1, True, 2**63])
def test_session_rejects_invalid_cache_budgets(budget):
    with pytest.raises(native.NativeError, match="Cache budgets"):
        native.NativeSession(cache_bytes=budget)


def test_initialization_rejects_changed_catalog(fake_native_process, monkeypatch):
    class File:
        def joinpath(self, *args):
            return self
        def read_bytes(self):
            return b"{}"
    monkeypatch.setattr(native, "files", lambda _: File())
    with pytest.raises(native.NativeError, match="Catalog changed"):
        native.NativeSession()
    assert fake_native_process["children"][0].waited


def test_initialization_rejects_wrong_handshake(fake_native_process):
    fake_native_process["replies"] = [{"ready": True, "abi": "wrong"}]
    with pytest.raises(native.NativeError, match="initialization response differs"):
        native.NativeSession()
    assert fake_native_process["children"][0].waited


def test_reduce_word_resume_and_destructor(fake_native_process):
    from dataclasses import replace
    from zipper_algebra.zipping_algorithm.interpreter import Limits
    result = native.reduce_word([], n=1)
    assert result["status"] == "complete"
    assert native.resume(result)["input_word"] == []
    assert native.resume(result, limits=replace(Limits(), max_steps=10001))["limits"]["max_steps"] == 10001
    with pytest.raises(native.NativeError, match="must not decrease"):
        native.resume(result, limits=replace(Limits(), max_steps=1))
    session = native.NativeSession()
    process = session._process
    session.__del__()
    assert session.pid is None and process.waited


@pytest.mark.parametrize("stderr_present", [False, True])
def test_protocol_eof_closes_session(stderr_present):
    import io
    from threading import Lock
    class Process:
        stdin = io.StringIO()
        stdout = io.StringIO()
        def kill(self):
            pass
        def wait(self, timeout):
            return 0
    session = object.__new__(native.NativeSession)
    session._lock = Lock()
    session._process = Process()
    session._stderr = io.StringIO("failed native call") if stderr_present else None
    with pytest.raises(native.NativeError, match="terminated without a response"):
        session._request({})
    assert session.pid is None


def test_close_kills_child_that_does_not_exit():
    from threading import Lock
    class Process:
        stdin = stdout = None
        def __init__(self):
            self.waits = 0
            self.killed = False
        def wait(self, timeout):
            self.waits += 1
            if self.waits == 1:
                raise subprocess.TimeoutExpired("engine", timeout)
            return -9
        def kill(self):
            self.killed = True
    session = object.__new__(native.NativeSession)
    session._lock = Lock()
    process = session._process = Process()
    session._stderr = None
    session.close()
    assert process.killed and process.waits == 2 and session.pid is None


def test_encode_nested_diagnostic_state_arrays():
    state = {"word": (1, -1, ["slider", 2, 1, False, False, False]),
             "contexts": [{"coefficient": [[0, 2**32770]]}]}
    encoded = native._encode_coefficients(state)
    assert encoded["word"] == [1, -1, ["slider", 2, 1, False, False, False]]
    assert native._decimal_integer(encoded["contexts"][0]["coefficient"][0][1]) == 2**32770
    assert isinstance(state["word"], tuple)
