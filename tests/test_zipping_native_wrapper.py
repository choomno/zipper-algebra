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
              "binary_sha256": native._digest(executable), "flags": []}
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
