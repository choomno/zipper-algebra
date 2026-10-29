"""Distribution provenance must describe the implementation bytes being shipped."""
from copy import deepcopy
import hashlib
import importlib.util
import importlib.metadata
import json
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tarfile
import pytest
from unittest.mock import patch
import zipfile

from ci.check_distribution import build_provenance
from zipper_algebra import provenance


class TestDistributionProvenance:
    @pytest.fixture(autouse=True)
    def distribution(self):
        self.files = {"zipper_algebra/__init__.py": b'"""Example package."""\n',
                      "zipper_algebra/provenance.py": b"SCHEMA = 'example'\n"}
        sources = {"src/" + name: {"bytes": len(content),
                                    "sha256": hashlib.sha256(content).hexdigest()}
                   for name, content in self.files.items()}
        self.metadata = {"schema": "zipper-build-provenance-v1",
                         "origin": {"schema": "zipper-source-provenance-v1",
                                    "omitted_source_count": 0,
                                    "source_files": sources,
                                    "source_sha256": hashlib.sha256(json.dumps(
                                        sources, sort_keys=True, separators=(",", ":")
                                    ).encode()).hexdigest()}}

    def check(self, metadata=None, files=None):
        files = self.files if files is None else files
        metadata = self.metadata if metadata is None else metadata
        return build_provenance(json.dumps(metadata), list(files), files.__getitem__,
                                "zipper_algebra/")

    def test_recorded_source_matches_every_packaged_python_file(self):
        assert (self.check()) == (self.metadata)

    def test_changed_bytes_and_incomplete_source_inventories_are_rejected(self):
        changed = dict(self.files, **{"zipper_algebra/provenance.py": b"modified\n"})
        with pytest.raises(AssertionError):
            self.check(files=changed)
        extra = dict(self.files, **{"zipper_algebra/unrecorded.py": b"new\n"})
        with pytest.raises(AssertionError):
            self.check(files=extra)
        missing = dict(self.files)
        del missing["zipper_algebra/provenance.py"]
        with pytest.raises(AssertionError):
            self.check(files=missing)
        forged = deepcopy(self.metadata)
        forged["origin"]["source_files"]["src/zipper_algebra/provenance.py"]["sha256"] = "0" * 64
        with pytest.raises(AssertionError):
            self.check(forged)


def current_setuptools_available():
    return (importlib.util.find_spec("setuptools") is not None
            and int(importlib.metadata.version("setuptools").split(".", 1)[0]) >= 77)


@pytest.fixture(scope="module")
def source_build_hooks():
    root = Path(__file__).resolve().parents[2]
    with patch("setuptools.setup"):
        return runpy.run_path(str(root / "setup.py"))


@pytest.mark.skipif(not current_setuptools_available(), reason="build hook checks need setuptools>=77")
class TestSourceArchiveProvenance:
    @pytest.fixture(autouse=True)
    def source_tree(self, tmp_path, source_build_hooks):
        self.metadata_hook = source_build_hooks["source_metadata"]
        self.ignore_hook = source_build_hooks["ignored_release_payloads"]
        self.build_class = source_build_hooks["ProvenanceBuildPy"]
        self.root = tmp_path
        self.package = self.root / "src/zipper_algebra"
        self.package.mkdir(parents=True)
        (self.package / "__init__.py").write_text('"""Source archive fixture."""\n')
        (self.package / "provenance.py").write_bytes(Path(provenance.__file__).read_bytes())
        (self.root / "setup.py").write_text('# Original build hook.\n')
        (self.root / "pyproject.toml").write_text('[project]\nname = "zipper-algebra"\n')
        (self.root / "experiments").mkdir()
        (self.root / "experiments/example.py").write_text('VALUE = 1\n')
        sources = {"src/zipper_algebra/" + p.name: provenance.file_identity(p)
                   for p in self.package.glob("*.py")}
        sources.update({name: provenance.file_identity(self.root / name)
                        for name in ("setup.py", "pyproject.toml")})
        sources["experiments/example.py"] = provenance.file_identity(
            self.root / "experiments/example.py")
        origin = {"schema": "zipper-source-provenance-v1", "origin": "checkout",
                  "git": {"commit": "a" * 40, "dirty": False, "changes": []},
                  "source_files": sources,
                  "source_sha256": hashlib.sha256(json.dumps(
                      sources, sort_keys=True, separators=(",", ":")
                  ).encode()).hexdigest(),
                  "changes_since_build": [], "extra_source_git": [], "omitted_source_count": 0}
        self.metadata = {"schema": "zipper-build-provenance-v1", "origin": origin}
        (self.package / "_build_provenance.json").write_text(json.dumps(self.metadata))

    def read_metadata(self):
        with patch.dict(self.metadata_hook.__globals__, ROOT=self.root, PACKAGE=self.package):
            return self.metadata_hook()

    def test_unchanged_source_archive_keeps_its_original_origin(self):
        assert (self.read_metadata()) == (self.metadata)

    def test_invalid_checkout_does_not_disable_release_filtering(self):
        (self.root / ".git").mkdir()
        with patch.dict(self.ignore_hook.__globals__, ROOT=self.root, PACKAGE=self.package):
            with pytest.raises(ValueError, match='classification failed'):
                self.ignore_hook(["tests/private_scratch.py"])

    def test_package_module_discovery_omits_unrecorded_payload_symlinks(self):
        from setuptools import Distribution
        (self.package / "private_notes.py").symlink_to(self.package / "__init__.py")
        distribution = Distribution({"packages": ["zipper_algebra"]})
        distribution.script_name = str(self.root / "setup.py")
        command = self.build_class(distribution)
        command.ensure_finalized()
        command._zipper_metadata = self.metadata
        with patch.dict(self.metadata_hook.__globals__, ROOT=self.root, PACKAGE=self.package):
            modules = command.find_package_modules("zipper_algebra", str(self.package))
        assert ({name for _, name, _ in modules}) == ({'__init__', 'provenance'})

    @pytest.mark.parametrize("escape", ["external", "parent-traversal"])
    def test_release_paths_keep_payload_names_and_reject_source_tree_escape(self, escape):
        relative = self.metadata_hook.__globals__["project_relative"]
        link = self.package / "private_notes.py"
        link.symlink_to(self.package / "__init__.py")
        with patch.dict(relative.__globals__, ROOT=self.root):
            assert relative(link) == "src/zipper_algebra/private_notes.py"
            assert relative(self.package / "__init__.py") == "src/zipper_algebra/__init__.py"
            path = (self.root.parent / "external.py" if escape == "external"
                    else self.root / ".." / "external.py")
            with pytest.raises(ValueError, match="leaves the project"):
                relative(path)

    def test_unversioned_source_export_reports_unknown_revision(self):
        (self.package / "_build_provenance.json").unlink()
        metadata = self.read_metadata()
        origin = metadata["origin"]
        assert (origin['origin']) == ('unversioned')
        assert (origin['git']) == ({'commit': None, 'dirty': None, 'changes': []})
        assert ('experiments/example.py') in (origin['source_files'])
        assert ('setup.py') in (origin['source_files'])

    def test_changed_source_cannot_inherit_the_old_clean_origin(self):
        (self.package / "__init__.py").write_text('"""Changed source."""\n')
        with pytest.raises(ValueError, match='package files differ'):
            self.read_metadata()

    def test_deleted_source_cannot_inherit_the_old_clean_origin(self):
        (self.package / "__init__.py").unlink()
        with pytest.raises(ValueError, match='package files differ'):
            self.read_metadata()

    def test_added_source_cannot_inherit_the_old_clean_origin(self):
        (self.package / "new.py").write_text('VALUE = 1\n')
        with pytest.raises(ValueError, match='package files differ'):
            self.read_metadata()

    @pytest.mark.parametrize('name', ['setup.py', 'pyproject.toml'])
    def test_changed_build_input_cannot_inherit_the_old_clean_origin(self, name):
        path = self.root / name
        original = path.read_bytes()
        try:
            path.write_bytes(original + b"\n# Changed build input\n")
            with pytest.raises(ValueError, match='package files differ'):
                self.read_metadata()
        finally:
            path.write_bytes(original)

    @pytest.mark.parametrize('name', ['setup.py', 'pyproject.toml'])
    def test_deleted_build_input_cannot_inherit_the_old_clean_origin(self, name):
        path = self.root / name
        original = path.read_bytes()
        try:
            path.unlink()
            with pytest.raises(ValueError, match='package files differ'):
                self.read_metadata()
        finally:
            path.write_bytes(original)

    def test_added_build_input_cannot_inherit_the_old_clean_origin(self):
        (self.root / "MANIFEST.in").write_text('include README.md\n')
        with pytest.raises(ValueError, match='package files differ'):
            self.read_metadata()

    def test_changed_recipe_cannot_inherit_the_old_clean_origin(self):
        (self.root / "experiments/example.py").write_text('VALUE = 2\n')
        with pytest.raises(ValueError, match='package files differ'):
            self.read_metadata()

    def test_deleted_recipe_cannot_inherit_the_old_clean_origin(self):
        (self.root / "experiments/example.py").unlink()
        with pytest.raises(ValueError, match='package files differ'):
            self.read_metadata()

    def test_added_recipe_cannot_inherit_the_old_clean_origin(self):
        (self.root / "experiments/added.py").write_text('VALUE = 3\n')
        with pytest.raises(ValueError, match='package files differ'):
            self.read_metadata()

    def test_legacy_fixture_without_build_inputs_remains_explicit(self):
        sources = self.metadata["origin"]["source_files"]
        for name in ("setup.py", "pyproject.toml"):
            (self.root / name).unlink()
            del sources[name]
        self.metadata["origin"]["source_sha256"] = hashlib.sha256(json.dumps(
            sources, sort_keys=True, separators=(",", ":")
        ).encode()).hexdigest()
        (self.package / "_build_provenance.json").write_text(json.dumps(self.metadata))
        assert (self.read_metadata()) == (self.metadata)
        (self.root / "setup.py").write_text('# Newly added build input.\n')
        with pytest.raises(ValueError, match='package files differ'):
            self.read_metadata()


def release_build_available():
    return (importlib.util.find_spec("build") is not None
            and current_setuptools_available()
            and shutil.which("git") is not None)


@pytest.mark.skipif(not (current_setuptools_available()), reason="editable metadata checks need setuptools>=77")
class TestEditableMetadata:
    def test_pep660_external_metadata_omits_local_extras(self, tmp_path):
        project = Path(__file__).resolve().parents[2]
        script = '''
from pathlib import Path
import json, sys
from setuptools.command.egg_info import egg_info
original = egg_info.run
def run(self):
    original(self)
    folder = Path(self.egg_info)
    assert not folder.resolve().is_relative_to(Path.cwd().resolve())
    (folder / 'private_scratch.py').write_text('UNPUBLISHED_EXTERNAL_METADATA')
egg_info.run = run
from setuptools.build_meta import prepare_metadata_for_build_editable
name = prepare_metadata_for_build_editable(sys.argv[1])
print('ZIPPER_METADATA_RESULT=' + json.dumps(name))
'''
        directory = tmp_path
        output = Path(directory)
        result = subprocess.run([sys.executable, "-I", "-c", script, str(output)],
                                cwd=project, text=True, capture_output=True, timeout=60)
        assert (result.returncode) == (0), result.stdout + result.stderr
        line = next(line for line in result.stdout.splitlines()
                    if line.startswith("ZIPPER_METADATA_RESULT="))
        name = json.loads(line.split("=", 1)[1])
        metadata = output / name
        assert (metadata / 'METADATA').is_file()
        assert not ((metadata / 'private_scratch.py').exists())
        for path in metadata.rglob("*"):
            if path.is_file():
                assert (b'UNPUBLISHED_EXTERNAL_METADATA') not in (path.read_bytes())
        original = next(output.glob("*.egg-info")) / "private_scratch.py"
        assert (original.read_text()) == ('UNPUBLISHED_EXTERNAL_METADATA')


@pytest.mark.skipif(not (release_build_available()), reason="release fixture needs current build tools and Git")
class TestIgnoredPythonRelease:
    def test_ignored_local_python_files_are_excluded_from_wheel_and_source_archive(self, tmp_path):
        self.check_release(tmp_path)

    def test_ignored_symlink_does_not_disable_release_filtering(self, tmp_path):
        self.check_release(tmp_path, symlinks=True)

    def check_release(self, directory, *, symlinks=False):
        project = Path(__file__).resolve().parents[2]
        root = Path(directory)
        package = root / "src/zipper_algebra"
        package.mkdir(parents=True)
        shutil.copyfile(project / "setup.py", root / "setup.py")
        shutil.copyfile(Path(provenance.__file__), package / "provenance.py")
        (package / "__init__.py").write_text('"""Release fixture."""\n')
        (root / "pyproject.toml").write_text('''[build-system]
requires = ["setuptools>=77"]
build-backend = "setuptools.build_meta"
[project]
name = "zipper-algebra"
version = "0.0.0"
[tool.setuptools]
package-dir = {"" = "src"}
include-package-data = false
[tool.setuptools.packages.find]
where = ["src"]
[tool.setuptools.package-data]
zipper_algebra = ["_build_provenance.json"]
''')
        (root / "MANIFEST.in").write_text(
            "include setup.py pyproject.toml .gitignore\n"
            "recursive-include src *.py\nrecursive-include experiments *.py\n"
            "recursive-include tests *.py\nrecursive-include ci *.py\n"
            "recursive-include examples *.py\nrecursive-include zipping_algorithm *.py\n")
        (root / "README.md").write_text("Release build fixture.\n")
        (root / "experiments").mkdir()
        (root / "experiments/example.py").write_text("VALUE = 1\n")
        other_trees = ("tests", "ci", "examples", "zipping_algorithm")
        for name in other_trees:
            (root / name).mkdir()
            (root / name / "tracked_reference.py").write_text("VALUE = 1\n")
        (root / ".gitignore").write_text(
            "private_notes.py\nbuild/\ndist/\n*.egg-info/\n__pycache__/\n")
        subprocess.run(["git", "init", "-q", str(root)], check=True)
        subprocess.run(["git", "-C", str(root), "add", "."], check=True)
        subprocess.run(["git", "-C", str(root), "-c", "user.name=Release fixture",
                        "-c", "user.email=release-fixture@example.invalid", "commit", "-q",
                        "-m", "Public release fixture"], check=True)
        commit = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"],
                                         text=True).strip()
        sentinel = b"UNPUBLISHED_LOCAL_PYTHON_CONTENT"
        with (root / ".git/info/exclude").open("a") as ignore:
            ignore.write("\nprivate_scratch.py\ntracked_reference.py\nlocal\n")
        links = []
        if symlinks:
            target = root / "private-link-target"
            target.mkdir()
            (target / "linked_private.py").write_bytes(sentinel)
            for parent in (*other_trees, "src/zipper_algebra.egg-info"):
                link = root / parent / "local"
                link.parent.mkdir(parents=True, exist_ok=True)
                try:
                    link.symlink_to(target, target_is_directory=True)
                except OSError as error:
                    pytest.skip(f'directory symlinks unavailable: {error}')
                links.append(link)
        local_private = (package / "private_notes.py", root / "experiments/private_notes.py",
                         root / "src/zipper_algebra.egg-info/private_scratch.py",
                         *(root / name / "private_scratch.py" for name in other_trees))
        for path in (*local_private, root / "build/lib/zipper_algebra/private_notes.py"):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(sentinel)
        result = subprocess.run([sys.executable, "-I", "-m", "build", "--no-isolation",
                                 "--wheel", "--sdist"], cwd=root, text=True,
                                capture_output=True, timeout=60)
        assert (result.returncode) == (0), result.stdout + result.stderr
        assert not ((package / '_build_provenance.json').exists())
        for path in local_private:
            assert (path.read_bytes()) == (sentinel)
        for link in links:
            assert link.is_symlink()
            assert ((link / 'linked_private.py').read_bytes()) == (sentinel)
        with zipfile.ZipFile(next((root / "dist").glob("*.whl"))) as wheel:
            assert not (any((name.endswith('private_notes.py') for name in wheel.namelist())))
            assert not (any((name.endswith('private_scratch.py') for name in wheel.namelist())))
            assert not (any((name.endswith('linked_private.py') for name in wheel.namelist())))
            assert not (any((sentinel in wheel.read(name) for name in wheel.namelist())))
            origin = json.loads(wheel.read("zipper_algebra/_build_provenance.json"))["origin"]
            assert (origin['git']['commit']) == (commit)
            assert not (origin['git']['dirty'])
            assert (origin['omitted_source_count']) == (2)
            assert ('experiments/example.py') in (origin['source_files'])
            assert not (any((name.endswith('private_notes.py') for name in origin['source_files'])))
        with tarfile.open(next((root / "dist").glob("*.tar.gz"))) as source:
            assert not (any((name.endswith('private_notes.py') for name in source.getnames())))
            assert not (any((name.endswith('private_scratch.py') for name in source.getnames())))
            assert not (any((name.endswith('linked_private.py') for name in source.getnames())))
            for name in other_trees:
                assert any((path.endswith(f'/{name}/tracked_reference.py') for path in source.getnames()))
            for member in source:
                if member.isfile():
                    content = source.extractfile(member).read()
                    assert (sentinel) not in (content)
                    if member.name.endswith("SOURCES.txt"):
                        assert (b'private_notes.py') not in (content)
                        assert (b'private_scratch.py') not in (content)
                        assert (b'linked_private.py') not in (content)
