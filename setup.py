"""Preserve source provenance in built distributions without editing the checkout."""
from __future__ import annotations

import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import tempfile

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.sdist import sdist


ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT / "src" / "zipper_algebra"
METADATA_NAME = "_build_provenance.json"
GENERATED_METADATA = {"PKG-INFO", "SOURCES.txt", "dependency_links.txt", "entry_points.txt",
                      "requires.txt", "top_level.txt", "namespace_packages.txt"}


def source_metadata():
    """Retain the original build origin when rebuilding an extracted sdist."""
    namespace = runpy.run_path(str(PACKAGE / "provenance.py"))
    recorded = PACKAGE / METADATA_NAME
    if recorded.is_file() and not (ROOT / ".git").exists():
        metadata = json.loads(recorded.read_text(encoding="utf-8"))
        if (not isinstance(metadata, dict)
                or set(metadata) != {"schema", "origin"}
                or metadata["schema"] != "zipper-build-provenance-v1"
                or not isinstance(metadata["origin"], dict)
                or metadata["origin"].get("schema") != "zipper-source-provenance-v1"):
            raise ValueError("invalid embedded build provenance")
        namespace["validate_provenance"](metadata["origin"])
        prefix = "src/zipper_algebra/"
        build_sources = namespace["BUILD_SOURCES"]
        expected = {name: identity for name, identity in metadata["origin"]["source_files"].items()
                    if name.startswith((prefix, "experiments/")) or name in build_sources}
        current = namespace["current_source_files"]()
        for name in build_sources:
            path = ROOT / name
            if path.is_file():
                if not path.resolve().is_relative_to(ROOT):
                    raise ValueError("build source path leaves its source tree")
                current[name] = namespace["file_identity"](path)
        for path in (ROOT / "experiments").rglob("*.py"):
            if not path.resolve().is_relative_to(ROOT):
                raise ValueError("research source path leaves its source tree")
            current[path.relative_to(ROOT).as_posix()] = namespace["file_identity"](path)
        if current != expected:
            raise ValueError("source archive package files differ from embedded build provenance; "
                             "use a Git checkout to record a new source origin")
        return metadata
    recipes = tuple(sorted((ROOT / "experiments").rglob("*.py")))
    origin = namespace["capture_provenance"](extra_sources=recipes)
    return {"schema": "zipper-build-provenance-v1", "origin": origin}


def write_metadata(directory, metadata):
    path = Path(directory) / "zipper_algebra" / METADATA_NAME
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def project_relative(path):
    """Validate a lexical path in the source tree without following payload links."""
    try:
        relative = Path(path).absolute().relative_to(ROOT)
    except ValueError as error:
        raise ValueError("release payload path leaves the project") from error
    if ".." in relative.parts:
        raise ValueError("release payload path leaves the project")
    return relative.as_posix()


def ignored_release_payloads(names):
    """Classify release paths without letting one symlink disable Git filtering."""
    boundaries = {}
    for name in names:
        path = Path(name)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError("release payload path leaves the project")
        partial = Path()
        for part in path.parts:
            partial /= part
            if (ROOT / partial).is_symlink():
                boundaries[name] = partial.as_posix()
                break
    candidates = sorted({boundaries.get(name, name) for name in names})
    ignored = set()
    if (ROOT / ".git").exists() and candidates:
        namespace = runpy.run_path(str(PACKAGE / "provenance.py"))
        git = namespace["_git"]
        try:
            top = Path(os.fsdecode(git(ROOT, "rev-parse", "--show-toplevel")).strip()).resolve()
            if top != ROOT:
                raise ValueError("release Git checkout does not match its source directory")
            try:
                content = git(ROOT, "check-ignore", "-z", "--stdin",
                              input=b"\0".join(os.fsencode(name) for name in candidates) + b"\0")
            except subprocess.CalledProcessError as error:
                if error.returncode != 1:
                    raise
                content = b""
        except (OSError, subprocess.CalledProcessError) as error:
            raise ValueError("Git ignore classification failed; release build stopped") from error
        ignored = {os.fsdecode(name) for name in content.split(b"\0") if name}
    if any(boundary not in ignored for boundary in boundaries.values()):
        raise ValueError("release payload crosses a symbolic link without an ignore rule")
    return {name for name in names if boundaries.get(name, name) in ignored}


class ProvenanceBuildPy(build_py):
    def source_metadata(self):
        if not hasattr(self, "_zipper_metadata"):
            self._zipper_metadata = source_metadata()
        return self._zipper_metadata

    def find_package_modules(self, package, package_dir):
        allowed = self.source_metadata()["origin"]["source_files"]
        modules = super().find_package_modules(package, package_dir)
        return [entry for entry in modules
                if project_relative(entry[2]) in allowed]

    def run(self):
        metadata = self._zipper_metadata = source_metadata()
        destination = Path(self.build_lib) / "zipper_algebra"
        if destination.resolve().is_relative_to(PACKAGE):
            raise ValueError("distribution output must not replace the source package")
        super().run()
        allowed = metadata["origin"]["source_files"]
        for path in destination.rglob("*.py"):
            name = "src/zipper_algebra/" + path.relative_to(destination).as_posix()
            if name not in allowed:
                path.unlink()
        write_metadata(self.build_lib, metadata)

    def get_outputs(self, include_bytecode=1):
        outputs = super().get_outputs(include_bytecode)
        metadata = str(Path(self.build_lib) / "zipper_algebra" / METADATA_NAME)
        return outputs if metadata in outputs else [*outputs, metadata]


class ProvenanceSdist(sdist):
    def make_release_tree(self, base_dir, files):
        metadata = source_metadata()
        allowed = metadata["origin"]["source_files"]
        files = [project_relative(name) for name in files]
        ignored = ignored_release_payloads(files)
        generated_root = project_relative(self.get_finalized_command("egg_info").egg_info)
        if (ROOT / generated_root).is_symlink():
            raise ValueError("release metadata directory must not be a symbolic link")
        selected = []
        for name in files:
            relative = name
            generated = (Path(relative).parent.as_posix() == generated_root
                         and Path(relative).name in GENERATED_METADATA
                         and not (ROOT / relative).is_symlink())
            if relative in ignored and not generated:
                continue
            computation = (relative.startswith(("src/zipper_algebra/", "experiments/"))
                           and relative.endswith(".py"))
            if not computation or relative in allowed:
                selected.append(name)
        super().make_release_tree(base_dir, selected)
        for name in selected:
            if name.endswith(".egg-info/SOURCES.txt"):
                manifest = Path(base_dir) / name
                manifest.unlink()
                manifest.write_text("\n".join(selected) + "\n", encoding="utf-8")
        write_metadata(Path(base_dir) / "src", metadata)


class ProvenanceWheel(bdist_wheel):
    def egg2dist(self, egginfo_path, distinfo_path):
        """Remove ignored local extras from the temporary wheel metadata copy."""
        source = Path(self.get_finalized_command("egg_info").egg_info).absolute()
        if source.is_symlink():
            raise ValueError("release metadata directory must not be a symbolic link")
        internal = source.resolve().is_relative_to(ROOT)
        payloads = {(project_relative(path) if internal else path.relative_to(source).as_posix()): path
                    for path in source.rglob("*") if path.is_file() or path.is_symlink()}
        # PEP 660 creates fresh metadata outside the source tree. Only the known
        # generated files are eligible there; local extras have no Git policy.
        ignored = ignored_release_payloads(list(payloads)) if internal else set(payloads)
        with tempfile.TemporaryDirectory(prefix="zipper-wheel-metadata-") as directory:
            filtered = Path(directory) / "egg-info"
            filtered.mkdir()
            for name, path in payloads.items():
                generated = (path.parent == source and path.name in GENERATED_METADATA
                             and not path.is_symlink())
                if name in ignored and not generated:
                    continue
                destination = filtered / path.relative_to(source)
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, destination)
            super().egg2dist(str(filtered), distinfo_path)
        stage = Path(egginfo_path)
        if stage.resolve() != source.resolve():
            shutil.rmtree(stage)


setup(cmdclass={"build_py": ProvenanceBuildPy, "sdist": ProvenanceSdist,
               "bdist_wheel": ProvenanceWheel})
