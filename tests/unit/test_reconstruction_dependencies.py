"""The provenance record identifies actual dependency bytes and rejects races."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from zipper_algebra.execution import dependencies


@pytest.fixture
def installed(tmp_path, monkeypatch):
    paths = {"binding": tmp_path / "fflow.py", "extension": tmp_path / "_cffi_fflow.so",
             "library": tmp_path / "libfflow.so"}
    for name, path in paths.items():
        path.write_bytes((name + " exact bytes\n").encode())
    binding = SimpleNamespace(__file__=str(paths["binding"]), __version__="2.0",
                              __fork__="example/finiteflow", __fork_api_version__=1)
    extension = SimpleNamespace(__file__=str(paths["extension"]))
    modules = {"fflow": binding, "_cffi_fflow": extension}
    monkeypatch.setattr(dependencies.importlib, "import_module", modules.__getitem__)
    maps = {str(paths[name]): {(paths[name].stat().st_dev, paths[name].stat().st_ino)}
            for name in ("extension", "library")}
    monkeypatch.setattr(dependencies, "_mapped_files", lambda: maps)
    return paths, binding, extension, maps


def test_record_hashes_actual_three_files_without_installation_paths(installed):
    paths, _, _, _ = installed
    record = dependencies.reconstruction_provenance()
    assert record["schema"] == "zipper-reconstruction-dependencies-v1"
    assert record["implementation"] == "FiniteFlow"
    assert record["version"] == "2.0"
    assert record["fork"] == "example/finiteflow"
    assert record["fork_api_version"] == 1
    for name, path in paths.items():
        assert record["files"][name] == {"name": path.name, "bytes": path.stat().st_size,
                                         "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    assert str(paths["binding"].parent) not in json.dumps(record)
    assert "revision" not in record


def test_recapture_detects_python_binding_bytes_changed_between_stages(installed):
    paths, _, _, _ = installed
    before = dependencies.reconstruction_provenance()
    paths["binding"].write_text("# a changed implementation\n")
    after = dependencies.reconstruction_provenance()
    assert before != after
    assert before["files"]["library"] == after["files"]["library"]


def test_atomic_replacement_of_loaded_library_is_rejected(installed):
    paths, _, _, _ = installed
    replacement = paths["library"].with_name("replacement.so")
    replacement.write_bytes(paths["library"].read_bytes())
    replacement.replace(paths["library"])
    with pytest.raises(RuntimeError, match="differs from the loaded"):
        dependencies.reconstruction_provenance()


@pytest.mark.parametrize("mode", ["missing", "ambiguous", "deleted", "old-inode"])
def test_unavailable_or_ambiguous_loaded_library_fails_closed(installed, mode):
    paths, _, _, maps = installed
    original = str(paths["library"])
    if mode == "missing":
        del maps[original]
    elif mode == "ambiguous":
        second = paths["library"].with_name("libfflow.so.2")
        second.write_text("other version")
        maps[str(second)] = {(second.stat().st_dev, second.stat().st_ino)}
    elif mode == "deleted":
        maps[original + " (deleted)"] = maps.pop(original)
    else:
        maps[original].add((paths["library"].stat().st_dev, 0))
    with pytest.raises(RuntimeError, match="loaded"):
        dependencies.reconstruction_provenance()


def test_versioned_native_library_filename_is_supported(installed):
    paths, _, _, maps = installed
    versioned = paths["library"].with_name("libfflow.so.2")
    paths["library"].rename(versioned)
    maps[str(versioned)] = maps.pop(str(paths["library"]))
    assert dependencies.reconstruction_provenance()["files"]["library"]["name"] == "libfflow.so.2"


def test_extension_symlink_names_resolve_to_the_loaded_file(installed):
    paths, _, extension, _ = installed
    alias = paths["extension"].with_name("alias.so")
    alias.symlink_to(paths["extension"])
    extension.__file__ = str(alias)
    assert dependencies.reconstruction_provenance()["files"]["extension"]["name"] == "_cffi_fflow.so"


def test_overlay_device_alias_preserves_matching_mapped_inode(installed):
    paths, _, _, maps = installed
    original = str(paths["library"])
    maps[original] = {(os.makedev(0, 999), paths["library"].stat().st_ino)}
    assert dependencies.reconstruction_provenance()["files"]["library"]["name"] == "libfflow.so"


def test_unmapped_binding_extension_is_rejected(installed):
    paths, _, _, maps = installed
    del maps[str(paths["extension"])]
    with pytest.raises(RuntimeError, match="CFFI extension"):
        dependencies.reconstruction_provenance()


@pytest.mark.parametrize("attribute,value", [("__version__", None), ("__version__", ""),
                                             ("__fork__", None), ("__fork__", ""),
                                             ("__fork_api_version__", None),
                                             ("__fork_api_version__", True),
                                             ("__fork_api_version__", 0)])
def test_missing_or_malformed_fork_metadata_is_not_inferred_from_required_pin(installed, attribute, value):
    _, binding, _, _ = installed
    setattr(binding, attribute, value)
    with pytest.raises(RuntimeError, match="maintained fork API"):
        dependencies.reconstruction_provenance()


@pytest.mark.parametrize("which,value", [("binding", None), ("binding", ""),
                                         ("extension", None), ("extension", "")])
def test_anonymous_binding_or_extension_cannot_receive_a_byte_identity(installed, which, value):
    _, binding, extension, _ = installed
    setattr(binding if which == "binding" else extension, "__file__", value)
    with pytest.raises(RuntimeError, match="installed files"):
        dependencies.reconstruction_provenance()


def test_absent_native_installation_has_actionable_error(monkeypatch):
    def absent(_):
        raise ImportError("not installed")
    monkeypatch.setattr(dependencies.importlib, "import_module", absent)
    with pytest.raises(ImportError, match="docs/RECONSTRUCTION.md"):
        dependencies.reconstruction_provenance()


def test_unreadable_dependency_is_reported_as_runtime_failure(installed):
    paths, _, _, _ = installed
    paths["binding"].unlink()
    with pytest.raises(RuntimeError, match="could not read") as error:
        dependencies.reconstruction_provenance()
    assert isinstance(error.value.__cause__, OSError)


def test_nonregular_dependency_is_rejected_without_reading_a_device():
    with pytest.raises(RuntimeError, match="not a regular file"):
        dependencies._file_record(Path("/dev/null"))


@pytest.mark.parametrize("mutation", ["edit", "replace"])
def test_dependency_mutation_during_hashing_is_rejected(tmp_path, monkeypatch, mutation):
    path = tmp_path / "fflow.py"
    path.write_text("original\n")
    digest = hashlib.file_digest
    def mutate(stream, algorithm):
        result = digest(stream, algorithm)
        if mutation == "edit":
            path.write_text("changed\n")
        else:
            other = path.with_name("replacement.py")
            other.write_bytes(path.read_bytes())
            other.replace(path)
        return result
    monkeypatch.setattr(dependencies.hashlib, "file_digest", mutate)
    with pytest.raises(RuntimeError, match="changed while recording"):
        dependencies._file_record(path)


def test_linux_map_reader_collects_segments_and_distinguishes_file_generations(monkeypatch):
    text = ("123-456 r-xp 0 00:24 71 /tmp/a directory/libfflow.so\n"
            "456-789 r--p 10 00:24 71 /tmp/a directory/libfflow.so\n"
            "789-abc r-xp 0 00:24 72 /tmp/a directory/libfflow.so\n"
            "abc-def rw-p 0 00:00 0 [heap]\n"
            "def-fff rw-p 0 00:00 0\n")
    monkeypatch.setattr(dependencies.Path, "read_text", lambda _: text)
    assert dependencies._mapped_files() == {"/tmp/a directory/libfflow.so":
                                            {(os.makedev(0, 0x24), 71),
                                             (os.makedev(0, 0x24), 72)}}


def test_linux_maps_read_error_is_an_actionable_failure(installed, monkeypatch):
    def missing():
        raise OSError("maps unavailable")
    monkeypatch.setattr(dependencies, "_mapped_files", missing)
    with pytest.raises(RuntimeError, match="could not read"):
        dependencies.reconstruction_provenance()
