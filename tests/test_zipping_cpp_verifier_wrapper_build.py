"""Fail-closed build cache and producer provenance tests for the verifier."""
from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
from typing import Any

import pytest

from zipper_algebra.zipping_algorithm import native_verifier as v


@pytest.fixture
def build_environment(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Path]:
    source = tmp_path / "source"
    source.mkdir()
    (source / "main.cpp").write_text("main")
    (source / "core.h").write_text("header")
    libraries = {}
    for name in v.LIBRARIES:
        path = tmp_path / (name + ".so")
        path.write_text(name)
        libraries[name] = {"path": str(path), "sha256": v._hash(path)}
    identity = {"abi": v.ABI, "compiler": "fake", "compiler_version": "v1",
                "flags": list(v.FLAGS), "link_libraries": list(v.LIBRARIES), "libraries": libraries}
    monkeypatch.setattr(v, "SOURCE_DIRECTORY", source)
    monkeypatch.setattr(v, "_compiler_identity", lambda *_: deepcopy(identity))
    monkeypatch.setattr(v, "_compile", lambda command: Path(command[-1]).write_bytes(b"compiled"))
    return source, tmp_path / "cache"


def test_source_inventory_rejects_incomplete_and_escaped_files(tmp_path: Path) -> None:
    with pytest.raises(v.VerifierBuildError, match="incomplete"):
        v._sources(tmp_path)
    (tmp_path / "main.cpp").write_text("main")
    outside = tmp_path.parent / (tmp_path.name + "-outside.h")
    outside.write_text("external")
    (tmp_path / "alias.h").symlink_to(outside)
    with pytest.raises(v.VerifierBuildError, match="leaves"):
        v._sources(tmp_path)


@pytest.mark.parametrize("content", [None, "[]", "not json"])
def test_invalid_manifest(tmp_path: Path, content: str | None) -> None:
    path = tmp_path / "manifest.json"
    if content is not None:
        path.write_text(content)
    with pytest.raises(v.VerifierBuildError, match="manifest"):
        v._load(path)


def test_partial_publication_is_rebuilt(build_environment: tuple[Path, Path]) -> None:
    _, cache = build_environment
    binary, _, _ = v.build_verifier(build_directory=cache)
    binary.with_suffix(".json").unlink()
    obj = next(cache.rglob("*.o"))
    obj.with_suffix(".json").unlink()
    binary, metadata, rebuilt = v.build_verifier(build_directory=cache)
    assert rebuilt and metadata["compiled_translation_units"] == 1
    assert v._authenticated_binary(binary) == metadata


@pytest.mark.parametrize("tamper_manifest", [False, True])
def test_object_cache_tamper_rejected(build_environment: tuple[Path, Path], tamper_manifest: bool) -> None:
    _, cache = build_environment
    binary, _, _ = v.build_verifier(build_directory=cache)
    binary.unlink()
    obj = next(cache.rglob("*.o"))
    if tamper_manifest:
        record = json.loads(obj.with_suffix(".json").read_text())
        record["identity"]["flags"] = []
        obj.with_suffix(".json").write_text(json.dumps(record))
    else:
        obj.write_bytes(b"tampered")
    with pytest.raises(v.VerifierBuildError, match="object cache"):
        v.build_verifier(build_directory=cache)


def test_link_time_source_drift_is_rejected(
    build_environment: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch,
) -> None:
    source, cache = build_environment
    original = v._compile

    def compile_(command: list[str]) -> None:
        original(command)
        if "-c" not in command:
            (source / "main.cpp").write_text("changed")

    monkeypatch.setattr(v, "_compile", compile_)
    with pytest.raises(v.VerifierBuildError, match="linking"):
        v.build_verifier(build_directory=cache)
    assert not list(cache.rglob("zipper-verifier"))


def test_installed_cache_location(
    build_environment: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch,
) -> None:
    _, cache = build_environment
    monkeypatch.setenv("XDG_CACHE_HOME", str(cache))
    binary, _, _ = v.build_verifier()
    assert binary.is_relative_to(cache / "zipper-algebra" / "native-verifier")


@pytest.mark.parametrize("mutation", ["abi", "source", "binary", "expected", "missing_libraries", "wrong_libraries", "library_record", "library_hash"])
def test_binary_authentication_rejects_mismatch(build_environment: tuple[Path, Path], mutation: str) -> None:
    _, cache = build_environment
    binary, metadata, _ = v.build_verifier(build_directory=cache)
    expected = None
    if mutation == "abi":
        metadata["abi"] = "other"
    elif mutation == "source":
        metadata["source_sha256"] = {}
    elif mutation == "binary":
        binary.write_bytes(b"tampered")
    elif mutation == "expected":
        expected = {"flags": ["unexpected"]}
    elif mutation == "missing_libraries":
        metadata.pop("libraries")
    elif mutation == "wrong_libraries":
        metadata["libraries"] = {}
    elif mutation == "library_record":
        metadata["libraries"]["gmp"] = []
    else:
        metadata["libraries"]["gmp"]["sha256"] = "bad"
    binary.with_suffix(".json").write_text(json.dumps(metadata))
    with pytest.raises(v.VerifierBuildError):
        v._authenticated_binary(binary, expected)


@pytest.mark.parametrize("error", [OSError("not executable"), subprocess.CalledProcessError(1, ["cc"], stderr="compile error")])
def test_compiler_failure_is_reported(monkeypatch: pytest.MonkeyPatch, error: Exception) -> None:
    def fail(*args: Any, **kwargs: Any) -> None:
        raise error
    monkeypatch.setattr(v.subprocess, "run", fail)
    with pytest.raises(v.VerifierBuildError, match="compilation failed"):
        v._compile(["cc"])


def test_compiler_identity_missing_compiler(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(v.shutil, "which", lambda _: None)
    with pytest.raises(v.VerifierBuildError, match="unavailable"):
        v._compiler_identity("missing", ())


@pytest.mark.parametrize("mode", ["ok", "missing_library", "failure"])
def test_compiler_identity_libraries(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str) -> None:
    compiler = tmp_path / "compiler"
    monkeypatch.setattr(v.shutil, "which", lambda _: str(compiler))
    library = tmp_path / "lib.so"
    library.write_bytes(b"library")

    def run(command: list[str], **_: Any) -> subprocess.CompletedProcess[str]:
        if mode == "failure":
            raise OSError("cannot run")
        answer = "version" if command[-1] == "--version" else str(library if mode == "ok" else tmp_path / "absent")
        return subprocess.CompletedProcess(command, 0, stdout=answer)

    monkeypatch.setattr(v.subprocess, "run", run)
    if mode == "ok":
        record = v._compiler_identity("cc", ("-O3",))
        assert record["compiler_version"] == "version"
        assert record["libraries"]["gmp"]["sha256"] == v._hash(library)
    else:
        with pytest.raises(v.VerifierBuildError):
            v._compiler_identity("cc", ())


@pytest.fixture
def producer_metadata(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    source = tmp_path / "producer"
    source.mkdir()
    (source / "main.cpp").write_text("main")
    package = tmp_path / "package"
    package.mkdir()
    (package / "native.py").write_text("wrapper")
    binary = tmp_path / "producer-bin"
    binary.write_bytes(b"producer")
    monkeypatch.setattr(v, "_PRODUCER_SOURCES", source)
    monkeypatch.setattr(v, "_PACKAGE", package)
    metadata = {key: "test" for key in v._BUILD_FIELDS}
    metadata.update(binary=str(binary), source_sha256=v._sources(source), binary_sha256=v._hash(binary),
                    wrapper_source_sha256=v._hash(package / "native.py"))
    binary.with_suffix(".json").write_text(json.dumps(metadata))
    return metadata


def test_producer_identity_is_bound_and_copied(producer_metadata: dict[str, Any]) -> None:
    identity = v._producer_identity(producer_metadata, "catalog", {"reference": "hash"})
    assert identity["catalog_sha256"] == "catalog"
    assert identity["native_build"]["compiler"] == "test"
    producer_metadata["compiler"] = "modified"
    assert identity["native_build"]["compiler"] == "test"


@pytest.mark.parametrize("field", ["source_sha256", "binary_sha256", "wrapper_source_sha256", "compiler"])
def test_producer_identity_rejects_changed_fields(producer_metadata: dict[str, Any], field: str) -> None:
    producer_metadata[field] = "tampered"
    with pytest.raises(v.VerifierError, match="provenance differs"):
        v._producer_identity(producer_metadata, "catalog", {})


@pytest.mark.parametrize("value", [None, {}, {"binary": 3}])
def test_producer_identity_requires_metadata(value: Any) -> None:
    with pytest.raises(v.VerifierError, match="metadata"):
        v._producer_identity(value, "catalog", {})


def test_reference_identity_hashes_only_declared_files(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(v, "_PACKAGE", tmp_path)
    for name in v._REFERENCE_NAMES:
        (tmp_path / name).write_text(name)
    (tmp_path / "unrelated.py").write_text("not reference provenance")
    assert v.reference_identity() == {name: v._hash(tmp_path / name) for name in v._REFERENCE_NAMES}


def test_checkout_default_cache_location(
    build_environment: tuple[Path, Path], monkeypatch: pytest.MonkeyPatch,
) -> None:
    source, _ = build_environment
    # Three nested parents reproduce the repository/native/source layout and
    # ensure the test never writes outside its temporary directory.
    nested = source / "repo" / "native" / "verifier"
    nested.mkdir(parents=True)
    (nested / "main.cpp").write_text("main")
    monkeypatch.setattr(v, "SOURCE_DIRECTORY", nested)
    monkeypatch.setattr(v, "_CHECKOUT_SOURCES", nested)
    binary, _, _ = v.build_verifier()
    assert binary.is_relative_to(source / ".local-tools" / "native-verifier")
