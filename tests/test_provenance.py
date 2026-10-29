"""Revision and byte provenance through real repositories and exported sources."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import pytest
from unittest.mock import patch

from zipper_algebra import provenance as source


class ProvenanceFixture:
    @pytest.fixture(autouse=True)
    def resources(self, request):
        self.temporary = tempfile.TemporaryDirectory(prefix='zipper-provenance-')
        request.addfinalizer(self.temporary.cleanup)
        self.base = Path(self.temporary.name).resolve()
        self.root, self.package = self.project(self.base / 'project')

    @staticmethod
    def project(root):
        package = root / 'src/zipper_algebra'
        package.mkdir(parents=True)
        (root / 'pyproject.toml').write_text('[project]\nname="zipper-algebra"\n', newline="\n")
        (root / '.gitignore').write_text('CONTRIBUTING.md\n.claude/\nruns/\n', newline="\n")
        (root / 'README.md').write_text('Public documentation\n', newline="\n")
        (package / '__init__.py').write_text('VERSION = "test"\n', newline="\n")
        (package / 'alpha.py').write_text('VALUE = 7\n', newline="\n")
        (package / 'provenance.py').write_text('# provenance module fixture\n', newline="\n")
        return root, package

    @staticmethod
    def git(root, *arguments):
        return subprocess.check_output(['git', '-C', str(root), *arguments], stderr=subprocess.DEVNULL).decode().strip()

    def initialize(self, root=None, commit=True):
        root = root or self.root
        self.git(root, 'init', '-q')
        self.git(root, 'config', 'user.name', 'Test author')
        self.git(root, 'config', 'user.email', 'test@example.invalid')
        self.git(root, 'config', 'core.autocrlf', 'false')
        if commit:
            self.git(root, 'add', '.')
            self.git(root, 'commit', '-q', '-m', 'Initial public sources')
            return self.git(root, 'rev-parse', 'HEAD')

    def remove_git(self):
        shutil.rmtree(self.root / '.git')

    def capture(self, package=None, **kwargs):
        with patch.object(source, '__file__', str((package or self.package) / 'provenance.py')):
            return source.capture_provenance(**kwargs)

    def installed(self, baseline):
        package = self.base / 'environment/site-packages/zipper_algebra'
        package.parent.mkdir(parents=True)
        shutil.copytree(self.package, package)
        (package / '_build_provenance.json').write_text(json.dumps(
            {'schema': source.BUILD_SCHEMA, 'origin': baseline}))
        return package


class TestSourceSnapshots(ProvenanceFixture):
    def test_clean_checkout_records_exact_revision_and_sources(self):
        revision = self.initialize()
        result = self.capture()
        assert result['origin'] == 'checkout'
        assert result['git'] == {'commit': revision, 'dirty': False, 'changes': []}
        assert set(result['source_files']) == {'src/zipper_algebra/' + name for name in ('alpha.py', '__init__.py', 'provenance.py')} | {'pyproject.toml'}
        assert result['source_files']['src/zipper_algebra/alpha.py'] == {'bytes': 10, 'sha256': hashlib.sha256(b'VALUE = 7\n').hexdigest()}
        assert str(self.root) not in json.dumps(result)
        assert result == self.capture()

    def test_source_only_rehash_detects_additions_modifications_and_deletions_without_git(self):
        self.initialize()
        with patch.object(source, '__file__', str(self.package / 'provenance.py')), \
             patch.object(source.subprocess, 'check_output', side_effect=AssertionError('Git must not run')):
            original = source.current_source_files()
            (self.package / 'alpha.py').unlink()
            (self.package / '__init__.py').write_text('VERSION = "updated"\n')
            (self.package / 'new.py').write_text('VALUE = 12\n')
            current = source.current_source_files()
        assert 'src/zipper_algebra/alpha.py' not in current
        assert 'src/zipper_algebra/new.py' not in original
        assert 'src/zipper_algebra/new.py' in current
        assert original['src/zipper_algebra/__init__.py'] != current['src/zipper_algebra/__init__.py']
        assert current['src/zipper_algebra/new.py'] == source.file_identity(self.package / 'new.py')

    def test_git_queries_disable_optional_index_writes(self):
        self.initialize()
        index = self.root / '.git/index'
        original = index.read_bytes()
        (self.package / 'alpha.py').write_text('VALUE = 123\n')
        with patch.object(source.subprocess, 'check_output',
                          wraps=subprocess.check_output) as execute:
            assert self.capture()['git']['dirty']
        assert execute.call_args_list
        for call in execute.call_args_list:
            assert call.kwargs['env']['GIT_OPTIONAL_LOCKS'] == '0'
        assert index.read_bytes() == original

    def test_build_sources_include_untracked_hooks_and_exact_configuration_hashes(self):
        self.initialize()
        setup = self.root / 'setup.py'
        setup.write_text('# local build hook\n')
        manifest = self.root / 'MANIFEST.in'
        manifest.write_text('include README.md\n')
        result = self.capture()
        assert result['git']['dirty']
        for path in (setup, manifest, self.root / 'pyproject.toml'):
            assert result['source_files'][path.name] == source.file_identity(path)
        changes = {entry['path']: entry for entry in result['git']['changes']}
        assert changes['setup.py']['index'] == '?'
        assert changes['MANIFEST.in']['index'] == '?'
        with patch.object(source, '__file__', str(self.package / 'provenance.py')):
            package_only = source.current_source_files()
        assert not set(source.BUILD_SOURCES) & package_only.keys()

    def test_build_source_symlink_cannot_escape_the_declared_project(self):
        external = self.base / 'outside_setup.py'
        external.write_text('# not project source\n')
        try:
            (self.root / 'setup.py').symlink_to(external)
        except OSError:
            pytest.skip('symlinks are unavailable')
        with pytest.raises(ValueError, match='build source path leaves'):
            self.capture()

    def test_git_change_names_preserve_legal_colons_and_newlines(self, subtests):
        self.initialize()
        names = ('calculation:notes.txt', 'calculation\nnotes.txt')
        for name in names:
            try:
                (self.root / name).write_text('original\n')
            except OSError:
                pytest.skip('filesystem does not support legal POSIX names')
        self.git(self.root, 'add', '.')
        self.git(self.root, 'commit', '-q', '-m', 'Track legal names')
        for name in names:
            (self.root / name).write_text('changed\n')
        result = self.capture()
        assert {entry['path'] for entry in result['git']['changes']} == set(names)
        encoded = json.dumps(result)
        assert 'calculation\\nnotes.txt' in encoded
        assert 'changed' not in encoded
        assert json.loads(encoded) == result
        for invalid in ('/absolute', '../outside', 'null\0name', None):
            with subtests.test(invalid=invalid), pytest.raises(ValueError):
                source._git_relative(invalid)

    def test_ignored_untracked_python_is_excluded_but_tracked_matches_remain(self):
        recipe = self.root / 'experiments/public.py'
        recipe.parent.mkdir()
        recipe.write_text('PUBLIC = 1\n')
        ignored = self.package / 'secret_notes.py'
        ignored_recipe = recipe.with_name('secret_notes.py')
        tracked = self.package / 'tracked_allowed.py'
        tracked.write_text('TRACKED = 1\n')
        with (self.root / '.gitignore').open('a') as stream:
            stream.write('**/secret_notes.py\n**/tracked_allowed.py\n')
        self.initialize(commit=False)
        self.git(self.root, 'add', '.')
        self.git(self.root, 'add', '-f', 'src/zipper_algebra/tracked_allowed.py')
        self.git(self.root, 'commit', '-q', '-m', 'Public source and ignore rules')
        before = self.capture(extra_sources=[recipe])
        ignored.write_text('PRIVATE_NOTE = 1\n')
        ignored_recipe.write_text('PRIVATE_NOTE = 2\n')
        after = self.capture(extra_sources=[recipe, ignored_recipe, ignored_recipe])
        assert before['omitted_source_count'] == 0
        assert after['omitted_source_count'] == 2
        assert {key: value for key, value in before.items() if key != 'omitted_source_count'} == {key: value for key, value in after.items() if key != 'omitted_source_count'}
        assert before['source_sha256'] == after['source_sha256']
        assert not after['git']['dirty']
        assert 'secret_notes' not in json.dumps(after)
        assert 'src/zipper_algebra/tracked_allowed.py' in after['source_files']
        public_new = self.package / 'public_new.py'
        public_new.write_text('PUBLIC = 2\n')
        current = self.capture(extra_sources=[recipe])
        assert current['git']['dirty']
        assert current['omitted_source_count'] == 1
        assert 'src/zipper_algebra/public_new.py' in current['source_files']
        assert 'secret_notes' not in json.dumps(current)

    def test_declared_source_rehash_ignores_new_unused_files_and_detects_known_drift(self, subtests):
        self.initialize()
        producer = self.capture()
        expected = {name: identity for name, identity in producer['source_files'].items()
                    if name.startswith(source.PREFIX)}
        with patch.object(source, '__file__', str(self.package / 'provenance.py')), \
             patch.object(source.subprocess, 'check_output', side_effect=AssertionError('Git must not run')):
            (self.package / 'unused_new.py').write_text('UNUSED = 1\n')
            assert source.current_source_files(expected) == expected
            (self.package / 'alpha.py').write_text('VALUE = 22\n')
            assert source.current_source_files(expected) != expected
            (self.package / 'alpha.py').unlink()
            assert 'src/zipper_algebra/alpha.py' not in source.current_source_files(expected)
            for invalid in ('private.py', 'src/zipper_algebra/file.txt', '../outside.py'):
                with subtests.test(invalid=invalid), pytest.raises(ValueError):
                    source.current_source_files([invalid])

    def test_private_ignored_notes_and_outputs_do_not_change_snapshot(self):
        self.initialize()
        before = self.capture()
        (self.root / 'CONTRIBUTING.md').write_text('Private instruction text')
        (self.root / '.claude').mkdir()
        (self.root / '.claude/memory.md').write_text('Hidden session text')
        (self.root / 'runs').mkdir()
        (self.root / 'runs/output.txt').write_text('generated')
        (self.root / 'untracked-note.md').write_text('not a declared source')
        after = self.capture()
        assert before == after
        for hidden in ('Private instruction', 'Hidden session', 'CONTRIBUTING', 'memory.md', 'untracked-note'):
            assert hidden not in json.dumps(after)

    def test_staged_unstaged_and_declared_untracked_sources(self):
        self.initialize()
        alpha = self.package / 'alpha.py'
        alpha.write_text('VALUE = 8\n')
        self.git(self.root, 'add', 'src/zipper_algebra/alpha.py')
        alpha.write_text('VALUE = 9\n')
        new = self.package / 'new.py'
        new.write_text('VALUE = 2\n')
        (self.root / 'README.md').write_text('A local documentation change\n')
        result = self.capture()
        changes = {entry['path']: entry for entry in result['git']['changes']}
        assert result['git']['dirty']
        assert (changes['src/zipper_algebra/alpha.py']['index'], changes['src/zipper_algebra/alpha.py']['worktree']) == ('M', 'M')
        assert (changes['src/zipper_algebra/new.py']['index'], changes['src/zipper_algebra/new.py']['worktree']) == ('?', '?')
        assert changes['README.md']['sha256'] == source.file_identity(self.root / 'README.md')['sha256']
        assert result['source_files']['src/zipper_algebra/alpha.py'] == source.file_identity(alpha)
        assert 'A local documentation change' not in json.dumps(result)

    def test_deleted_and_renamed_files_keep_status_and_current_hash(self):
        self.initialize()
        self.git(self.root, 'mv', 'src/zipper_algebra/alpha.py', 'src/zipper_algebra/renamed.py')
        (self.root / 'README.md').unlink()
        changes = {entry['path']: entry for entry in self.capture()['git']['changes']}
        renamed = changes['src/zipper_algebra/renamed.py']
        assert renamed['index'] == 'R'
        assert renamed['previous_path'] == 'src/zipper_algebra/alpha.py'
        assert renamed['sha256'] == source.file_identity(self.package / 'renamed.py')['sha256']
        assert (changes['README.md']['bytes'], changes['README.md']['sha256']) == (None, None)

    def test_no_git_and_unborn_git_are_unknown(self):
        result = self.capture()
        assert result['origin'] == 'unversioned'
        assert result['git'] == {'commit': None, 'dirty': None, 'changes': []}
        self.initialize(commit=False)
        assert self.capture()['git'] == result['git']

    def test_missing_git_and_mismatched_toplevel_are_unknown(self):
        self.initialize()
        with patch.object(source.subprocess, 'check_output', side_effect=FileNotFoundError('git')):
            assert self.capture()['git']['commit'] is None
        with patch.object(source, '_git', return_value=str(self.base).encode()):
            assert self.capture()['git']['commit'] is None

    def test_source_export_does_not_inherit_an_enclosing_repository(self):
        self.initialize(self.base)
        assert self.capture()['git']['commit'] is None

    def test_package_source_symlink_cannot_leave_source_tree(self):
        external = self.base / 'outside.py'
        external.write_text('PRIVATE = "not source"\n')
        link = self.package / 'escape.py'
        try:
            link.symlink_to(external)
        except OSError:
            pytest.skip('symlinks are unavailable')
        with pytest.raises(ValueError, match='leaves its source tree'):
            self.capture()

    def test_missing_or_other_project_marker_is_not_a_checkout(self):
        self.initialize()
        marker = self.root / 'pyproject.toml'
        marker.unlink()
        assert self.capture()['git']['commit'] is None
        marker.write_text('[project]\nname="another-project"\n')
        assert self.capture()['git']['commit'] is None

    def test_empty_package_has_no_attributable_source(self):
        for path in self.package.glob('*.py'):
            path.unlink()
        with pytest.raises(ValueError, match='needs source files'):
            self.capture()


class TestBuildAndRecipeOrigins(ProvenanceFixture):
    def test_installed_build_keeps_original_commit_and_rehashes_current_files(self):
        revision = self.initialize()
        baseline = self.capture()
        package = self.installed(baseline)
        unchanged = self.capture(package)
        assert unchanged['git']['commit'] == revision
        assert unchanged['origin'] == 'build'
        assert unchanged['changes_since_build'] == []
        (package / 'alpha.py').unlink()
        (package / '__init__.py').write_text('VERSION = "modified"\n')
        (package / 'new.py').write_text('VALUE = 3\n')
        result = self.capture(package)
        assert result['git'] == baseline['git']
        assert result['changes_since_build'] == [{'path': 'src/zipper_algebra/__init__.py', 'status': 'modified'}, {'path': 'src/zipper_algebra/alpha.py', 'status': 'deleted'}, {'path': 'src/zipper_algebra/new.py', 'status': 'added'}]
        assert result['source_sha256'] != baseline['source_sha256']

    def test_own_unborn_checkout_never_reuses_embedded_old_commit(self):
        self.initialize()
        baseline = self.capture()
        (self.package / '_build_provenance.json').write_text(json.dumps(
            {'schema': source.BUILD_SCHEMA, 'origin': baseline}))
        self.remove_git()
        self.initialize(commit=False)
        assert self.capture()['origin'] == 'unversioned'
        assert self.capture()['git']['commit'] is None

    def test_export_recipe_and_build_changes_are_compared_with_their_baseline(self):
        recipes = self.root / 'experiments'
        recipes.mkdir()
        changed, deleted, unused = [recipes / (name + '.py') for name in ('changed', 'deleted', 'unused')]
        for path in (changed, deleted, unused):
            path.write_text('VALUE = 1\n')
        setup = self.root / 'setup.py'
        setup.write_text('# build setup\n')
        manifest = self.root / 'MANIFEST.in'
        manifest.write_text('include README.md\n')
        self.initialize()
        baseline = self.capture(extra_sources=[changed, deleted, unused])
        (self.package / '_build_provenance.json').write_text(json.dumps(
            {'schema': source.BUILD_SCHEMA, 'origin': baseline}))
        self.remove_git()
        changed.write_text('VALUE = 2\n')
        deleted.unlink()
        setup.unlink()
        manifest.write_text('include LICENSE\n')
        added = recipes / 'added.py'
        added.write_text('VALUE = 3\n')
        result = self.capture(extra_sources=[changed, deleted, added])
        assert result['git'] == baseline['git']
        assert result['changes_since_build'] == [{'path': 'MANIFEST.in', 'status': 'modified'}, {'path': 'experiments/added.py', 'status': 'added'}, {'path': 'experiments/changed.py', 'status': 'modified'}, {'path': 'experiments/deleted.py', 'status': 'deleted'}, {'path': 'setup.py', 'status': 'deleted'}]
        assert 'experiments/unused.py' not in result['source_files']
        assert result['extra_source_git'] == [{'paths': ['experiments/added.py', 'experiments/changed.py', 'experiments/deleted.py'], 'git': {'commit': None, 'dirty': None, 'changes': []}}]

    def test_recipe_without_build_baseline_has_explicit_unknown_current_origin(self):
        self.initialize()
        baseline = self.capture()
        (self.package / '_build_provenance.json').write_text(json.dumps(
            {'schema': source.BUILD_SCHEMA, 'origin': baseline}))
        self.remove_git()
        recipe = self.root / 'experiments/late.py'
        recipe.parent.mkdir()
        recipe.write_text('VALUE = 5\n')
        result = self.capture(extra_sources=[recipe])
        assert result['git'] == baseline['git']
        assert result['changes_since_build'] == [{'path': 'experiments/late.py', 'status': 'added'}]
        assert result['extra_source_git'] == [{'paths': ['experiments/late.py'], 'git': {'commit': None, 'dirty': None, 'changes': []}}]
        (self.package / '_build_provenance.json').unlink()
        unversioned = self.capture(extra_sources=[recipe])
        assert unversioned['origin'] == 'unversioned'
        assert unversioned['extra_source_git'][0]['git']['dirty'] == None

    def test_invalid_embedded_build_metadata_is_rejected(self, subtests):
        self.initialize()
        package = self.installed(self.capture())
        path = package / '_build_provenance.json'
        for bad in ([], {}, {'schema': 'other', 'origin': {}},
                    {'schema': source.BUILD_SCHEMA, 'origin': {}}):
            with subtests.test(bad=bad):
                path.write_text(json.dumps(bad))
                with pytest.raises(ValueError):
                    self.capture(package)
        path.write_text('not JSON')
        with pytest.raises(ValueError):
            self.capture(package)

    def test_declared_recipe_sources_are_hashed_and_track_untracked_files(self):
        self.initialize()
        recipe = self.root / 'experiments/recipe.py'
        recipe.parent.mkdir()
        recipe.write_text('RESULT = 1\n')
        result = self.capture(extra_sources=[recipe, recipe])
        assert result['source_files']['experiments/recipe.py'] == source.file_identity(recipe)
        assert result['extra_source_git'] == []
        assert 'experiments/recipe.py' in [item['path'] for item in result['git']['changes']]

    def test_wheel_and_recipe_checkout_have_separate_origins(self):
        self.initialize()
        baseline = self.capture()
        package = self.installed(baseline)
        other, _ = self.project(self.base / 'other')
        recipe = other / 'experiments/recipe.py'
        recipe.parent.mkdir()
        recipe.write_text('RESULT = 2\n')
        recipe_revision = self.initialize(other)
        result = self.capture(package, extra_sources=[recipe])
        assert result['git'] == baseline['git']
        assert result['extra_source_git'] == [{'paths': ['experiments/recipe.py'], 'git': {'commit': recipe_revision, 'dirty': False, 'changes': []}}]

    def test_additional_sources_are_confined_to_declared_public_python_files(self, subtests):
        private = self.root / 'private.py'
        private.write_text('PRIVATE = 1\n')
        outside = self.base / 'outside.py'
        outside.write_text('VALUE = 1\n')
        wrong_type = self.root / 'experiments/data.txt'
        wrong_type.parent.mkdir()
        wrong_type.write_text('not Python')
        for path in (private, outside, wrong_type, self.root / 'missing.py', wrong_type.parent / 'missing.py'):
            with subtests.test(path=path):
                with pytest.raises(ValueError):
                    self.capture(extra_sources=[path])
        baseline = self.capture()
        package = self.installed(baseline)
        with pytest.raises(ValueError, match='outside the declared project'):
            self.capture(package, extra_sources=[outside])
        a, _ = self.project(self.base / 'a')
        b, _ = self.project(self.base / 'b')
        paths = []
        for root in (a, b):
            (root / 'experiments').mkdir()
            path = root / 'experiments/recipe.py'
            path.write_text('VALUE = 1\n')
            paths.append(path)
        with pytest.raises(ValueError, match='one project'):
            self.capture(package, extra_sources=paths)


class TestMetadataValidation(ProvenanceFixture):
    def test_invalid_schema_git_source_and_build_change_metadata(self, subtests):
        value = self.capture()
        invalid = [None, {}, {**value, 'schema': 'other'}, {**value, 'origin': 'cwd'},
                   {**value, 'git': []}, {**value, 'git': {'commit': 'wrong', 'dirty': False, 'changes': []}},
                   {**value, 'git': {'commit': None, 'dirty': 1, 'changes': []}},
                   {**value, 'git': {'commit': None, 'dirty': None, 'changes': None}},
                   {**value, 'source_files': []}, {**value, 'source_files': {}},
                   {**value, 'source_sha256': '0' * 64}, {**value, 'changes_since_build': {}},
                   {**value, 'changes_since_build': [{}]}, {**value, 'extra_source_git': {}},
                   {**value, 'extra_source_git': [{}]},
                   {**value, 'omitted_source_count': -1}, {**value, 'omitted_source_count': True},
                   {**value, 'omitted_source_count': '0'}]
        for bad in invalid:
            with subtests.test(bad=bad):
                with pytest.raises(ValueError):
                    source.validate_provenance(bad)

    def test_invalid_source_paths_and_hash_shapes(self, subtests):
        value = self.capture()
        for name in ('/absolute.py', '../parent.py', 'x\\y.py', 'drive:x.py', 'bad\n.py', '', 'data.txt', 'private.py', 'private.toml'):
            bad = deepcopy(value)
            bad['source_files'] = {name: {'bytes': 0, 'sha256': '0' * 64}}
            with subtests.test(name=name), pytest.raises(ValueError):
                source.validate_provenance(bad)
        for identity in ({}, {'bytes': True, 'sha256': '0' * 64},
                         {'bytes': -1, 'sha256': '0' * 64}, {'bytes': 0, 'sha256': 'g' * 64},
                         {'bytes': 0, 'sha256': None}):
            bad = deepcopy(value)
            bad['source_files'] = {'src/zipper_algebra/a.py': identity}
            with subtests.test(identity=identity), pytest.raises(ValueError):
                source.validate_provenance(bad)
        with pytest.raises(ValueError):
            source._relative(None)

    def test_invalid_change_metadata_and_extra_origin_references(self, subtests):
        value = self.capture()
        valid = {'path': 'README.md', 'index': 'M', 'worktree': ' ', 'bytes': 0, 'sha256': '0' * 64}
        for change in ({}, {**valid, 'extra': 1}, {**valid, 'index': 'MM'},
                       {**valid, 'sha256': 'bad'}, {**valid, 'previous_path': '../old'}):
            bad = deepcopy(value)
            bad['git']['changes'] = [change]
            with subtests.test(change=change), pytest.raises(ValueError):
                source.validate_provenance(bad)
        for entry in ({'paths': ['missing.py'], 'git': value['git']},
                      {'paths': None, 'git': value['git']}):
            bad = deepcopy(value)
            bad['extra_source_git'] = [entry]
            with pytest.raises(ValueError):
                source.validate_provenance(bad)


class TestOutputTags(ProvenanceFixture):
    def test_sidecar_binds_output_and_prehashed_input_without_paths_or_contents(self):
        self.initialize()
        output = self.base / 'result.txt'
        output.write_bytes(b'certificate bytes\r\n')
        input_path = self.base / 'input.txt'
        input_path.write_bytes(b'original input')
        input_tag = {'name': input_path.name, **source.file_identity(input_path)}
        producer = self.capture()
        input_path.write_bytes(b'changed later')
        record = source.write_output_provenance(output, producer=producer, inputs=[input_tag])
        assert record['producer'] == producer
        assert record['inputs'] == [input_tag]
        assert record['outputs'] == {'result.txt': source.file_identity(output)}
        assert json.loads(output.with_name('result.txt.provenance.json').read_text()) == record
        assert output.read_bytes() == b'certificate bytes\r\n'
        assert str(self.base) not in json.dumps(record)
        assert 'original input' not in json.dumps(record)
        second = source.write_output_provenance(output, producer=producer, inputs=[input_path])
        assert second['inputs'][0]['sha256'] == source.file_identity(input_path)['sha256']

    def test_bundle_and_implicit_capture(self):
        directory = self.base / 'outputs'
        directory.mkdir()
        (directory / 'a.txt').write_text('a')
        (directory / 'b.json').write_text('{}')
        with patch.object(source, '__file__', str(self.package / 'provenance.py')):
            record = source.write_bundle_provenance(directory, ['a.txt', 'b.json'])
            single = source.write_output_provenance(directory / 'a.txt')
        assert record['schema'] == source.OUTPUT_SCHEMA
        assert set(record['outputs']) == {'a.txt', 'b.json'}
        assert record['producer'] == single['producer']
        assert json.loads((directory / 'provenance.json').read_text()) == record

    def test_unsafe_duplicate_missing_and_empty_outputs_do_not_publish(self, subtests):
        (self.base / 'ok.txt').write_text('ok')
        producer = self.capture()
        for names in ([], ['../outside'], ['/absolute'], [''], ['.'], ['..'],
                      ['provenance.json'], ['ok.txt', 'ok.txt'], ['missing.txt'], [None], ['windows\\file.txt'], ['drive:x']):
            with subtests.test(names=names), pytest.raises((ValueError, FileNotFoundError)):
                source.write_bundle_provenance(self.base, names, producer=producer)
        assert not (self.base / 'provenance.json').exists()

    def test_bad_prehashed_inputs_do_not_replace_existing_metadata(self, subtests):
        output = self.base / 'result.txt'
        output.write_text('result')
        producer = self.capture()
        original = source.write_output_provenance(output, producer=producer)
        for tag in ({}, {'name': '../x', 'bytes': 0, 'sha256': '0' * 64},
                    {'name': 'x', 'bytes': True, 'sha256': '0' * 64},
                    {'name': '..', 'bytes': 0, 'sha256': '0' * 64},
                    {'name': 'C:\\private', 'bytes': 0, 'sha256': '0' * 64},
                    {'name': None, 'bytes': 0, 'sha256': '0' * 64}):
            with subtests.test(tag=tag), pytest.raises(ValueError):
                source.write_output_provenance(output, producer=producer, inputs=[tag])
        assert json.loads(output.with_name('result.txt.provenance.json').read_text()) == original

    def test_failed_atomic_publication_preserves_previous_sidecar(self):
        output = self.base / 'result.txt'
        output.write_text('result')
        producer = self.capture()
        source.write_output_provenance(output, producer=producer)
        sidecar = output.with_name('result.txt.provenance.json')
        previous = sidecar.read_bytes()
        with patch.object(source.os, 'replace', side_effect=OSError('interrupted')):
            with pytest.raises(OSError):
                source.write_output_provenance(output, producer=producer)
        assert sidecar.read_bytes() == previous
        assert sorted((p.name for p in self.base.iterdir())) == ['project', 'result.txt', 'result.txt.provenance.json']
