"""Revision and byte provenance through real repositories and exported sources."""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from zipper_algebra import provenance as source


class ProvenanceFixture(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='zipper-provenance-')
        self.addCleanup(self.temporary.cleanup)
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
        def retry_readonly(function, path, error_info):
            error = error_info[1]
            if not isinstance(error, PermissionError):
                raise error
            os.chmod(path, os.stat(path).st_mode | stat.S_IWRITE)
            function(path)
        shutil.rmtree(self.root / '.git', onerror=retry_readonly)

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
        self.assertEqual(result['origin'], 'checkout')
        self.assertEqual(result['git'], {'commit': revision, 'dirty': False, 'changes': []})
        self.assertEqual(set(result['source_files']),
                         {'src/zipper_algebra/' + name for name in ('alpha.py', '__init__.py', 'provenance.py')}
                         | {'pyproject.toml'})
        self.assertEqual(result['source_files']['src/zipper_algebra/alpha.py'],
                         {'bytes': 10, 'sha256': hashlib.sha256(b'VALUE = 7\n').hexdigest()})
        self.assertNotIn(str(self.root), json.dumps(result))
        self.assertEqual(result, self.capture())

    def test_source_only_rehash_detects_additions_modifications_and_deletions_without_git(self):
        self.initialize()
        with patch.object(source, '__file__', str(self.package / 'provenance.py')), \
             patch.object(source.subprocess, 'check_output', side_effect=AssertionError('Git must not run')):
            original = source.current_source_files()
            (self.package / 'alpha.py').unlink()
            (self.package / '__init__.py').write_text('VERSION = "updated"\n')
            (self.package / 'new.py').write_text('VALUE = 12\n')
            current = source.current_source_files()
        self.assertNotIn('src/zipper_algebra/alpha.py', current)
        self.assertNotIn('src/zipper_algebra/new.py', original)
        self.assertIn('src/zipper_algebra/new.py', current)
        self.assertNotEqual(original['src/zipper_algebra/__init__.py'],
                            current['src/zipper_algebra/__init__.py'])
        self.assertEqual(current['src/zipper_algebra/new.py'],
                         source.file_identity(self.package / 'new.py'))

    def test_git_queries_disable_optional_index_writes(self):
        self.initialize()
        index = self.root / '.git/index'
        original = index.read_bytes()
        (self.package / 'alpha.py').write_text('VALUE = 123\n')
        with patch.object(source.subprocess, 'check_output',
                          wraps=subprocess.check_output) as execute:
            self.assertTrue(self.capture()['git']['dirty'])
        self.assertTrue(execute.call_args_list)
        for call in execute.call_args_list:
            self.assertEqual(call.kwargs['env']['GIT_OPTIONAL_LOCKS'], '0')
        self.assertEqual(index.read_bytes(), original)

    def test_build_sources_include_untracked_hooks_and_exact_configuration_hashes(self):
        self.initialize()
        setup = self.root / 'setup.py'
        setup.write_text('# local build hook\n')
        manifest = self.root / 'MANIFEST.in'
        manifest.write_text('include README.md\n')
        result = self.capture()
        self.assertTrue(result['git']['dirty'])
        for path in (setup, manifest, self.root / 'pyproject.toml'):
            self.assertEqual(result['source_files'][path.name], source.file_identity(path))
        changes = {entry['path']: entry for entry in result['git']['changes']}
        self.assertEqual(changes['setup.py']['index'], '?')
        self.assertEqual(changes['MANIFEST.in']['index'], '?')
        with patch.object(source, '__file__', str(self.package / 'provenance.py')):
            package_only = source.current_source_files()
        self.assertFalse(set(source.BUILD_SOURCES) & package_only.keys())

    def test_build_source_symlink_cannot_escape_the_declared_project(self):
        external = self.base / 'outside_setup.py'
        external.write_text('# not project source\n')
        try:
            (self.root / 'setup.py').symlink_to(external)
        except OSError:
            self.skipTest('symlinks are unavailable')
        with self.assertRaisesRegex(ValueError, 'build source path leaves'):
            self.capture()

    def test_git_change_names_preserve_legal_colons_and_newlines(self):
        self.initialize()
        names = ('calculation:notes.txt', 'calculation\nnotes.txt')
        for name in names:
            try:
                (self.root / name).write_text('original\n')
            except OSError:
                self.skipTest('filesystem does not support legal POSIX names')
        self.git(self.root, 'add', '.')
        self.git(self.root, 'commit', '-q', '-m', 'Track legal names')
        for name in names:
            (self.root / name).write_text('changed\n')
        result = self.capture()
        self.assertEqual({entry['path'] for entry in result['git']['changes']}, set(names))
        encoded = json.dumps(result)
        self.assertIn('calculation\\nnotes.txt', encoded)
        self.assertNotIn('changed', encoded)
        self.assertEqual(json.loads(encoded), result)
        for invalid in ('/absolute', '../outside', 'null\0name', None):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
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
        self.assertEqual(before['omitted_source_count'], 0)
        self.assertEqual(after['omitted_source_count'], 2)
        self.assertEqual({key: value for key, value in before.items() if key != 'omitted_source_count'},
                         {key: value for key, value in after.items() if key != 'omitted_source_count'})
        self.assertEqual(before['source_sha256'], after['source_sha256'])
        self.assertFalse(after['git']['dirty'])
        self.assertNotIn('secret_notes', json.dumps(after))
        self.assertIn('src/zipper_algebra/tracked_allowed.py', after['source_files'])
        public_new = self.package / 'public_new.py'
        public_new.write_text('PUBLIC = 2\n')
        current = self.capture(extra_sources=[recipe])
        self.assertTrue(current['git']['dirty'])
        self.assertEqual(current['omitted_source_count'], 1)
        self.assertIn('src/zipper_algebra/public_new.py', current['source_files'])
        self.assertNotIn('secret_notes', json.dumps(current))

    def test_declared_source_rehash_ignores_new_unused_files_and_detects_known_drift(self):
        self.initialize()
        producer = self.capture()
        expected = {name: identity for name, identity in producer['source_files'].items()
                    if name.startswith(source.PREFIX)}
        with patch.object(source, '__file__', str(self.package / 'provenance.py')), \
             patch.object(source.subprocess, 'check_output', side_effect=AssertionError('Git must not run')):
            (self.package / 'unused_new.py').write_text('UNUSED = 1\n')
            self.assertEqual(source.current_source_files(expected), expected)
            (self.package / 'alpha.py').write_text('VALUE = 22\n')
            self.assertNotEqual(source.current_source_files(expected), expected)
            (self.package / 'alpha.py').unlink()
            self.assertNotIn('src/zipper_algebra/alpha.py', source.current_source_files(expected))
            for invalid in ('private.py', 'src/zipper_algebra/file.txt', '../outside.py'):
                with self.subTest(invalid=invalid), self.assertRaises(ValueError):
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
        self.assertEqual(before, after)
        for hidden in ('Private instruction', 'Hidden session', 'CONTRIBUTING', 'memory.md', 'untracked-note'):
            self.assertNotIn(hidden, json.dumps(after))

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
        self.assertTrue(result['git']['dirty'])
        self.assertEqual((changes['src/zipper_algebra/alpha.py']['index'],
                          changes['src/zipper_algebra/alpha.py']['worktree']), ('M', 'M'))
        self.assertEqual((changes['src/zipper_algebra/new.py']['index'],
                          changes['src/zipper_algebra/new.py']['worktree']), ('?', '?'))
        self.assertEqual(changes['README.md']['sha256'], source.file_identity(self.root / 'README.md')['sha256'])
        self.assertEqual(result['source_files']['src/zipper_algebra/alpha.py'], source.file_identity(alpha))
        self.assertNotIn('A local documentation change', json.dumps(result))

    def test_deleted_and_renamed_files_keep_status_and_current_hash(self):
        self.initialize()
        self.git(self.root, 'mv', 'src/zipper_algebra/alpha.py', 'src/zipper_algebra/renamed.py')
        (self.root / 'README.md').unlink()
        changes = {entry['path']: entry for entry in self.capture()['git']['changes']}
        renamed = changes['src/zipper_algebra/renamed.py']
        self.assertEqual(renamed['index'], 'R')
        self.assertEqual(renamed['previous_path'], 'src/zipper_algebra/alpha.py')
        self.assertEqual(renamed['sha256'], source.file_identity(self.package / 'renamed.py')['sha256'])
        self.assertEqual((changes['README.md']['bytes'], changes['README.md']['sha256']), (None, None))

    def test_no_git_and_unborn_git_are_unknown(self):
        result = self.capture()
        self.assertEqual(result['origin'], 'unversioned')
        self.assertEqual(result['git'], {'commit': None, 'dirty': None, 'changes': []})
        self.initialize(commit=False)
        self.assertEqual(self.capture()['git'], result['git'])

    def test_missing_git_and_mismatched_toplevel_are_unknown(self):
        self.initialize()
        with patch.object(source.subprocess, 'check_output', side_effect=FileNotFoundError('git')):
            self.assertIsNone(self.capture()['git']['commit'])
        with patch.object(source, '_git', return_value=str(self.base).encode()):
            self.assertIsNone(self.capture()['git']['commit'])

    def test_source_export_does_not_inherit_an_enclosing_repository(self):
        self.initialize(self.base)
        self.assertIsNone(self.capture()['git']['commit'])

    def test_package_source_symlink_cannot_leave_source_tree(self):
        external = self.base / 'outside.py'
        external.write_text('PRIVATE = "not source"\n')
        link = self.package / 'escape.py'
        try:
            link.symlink_to(external)
        except OSError:
            self.skipTest('symlinks are unavailable')
        with self.assertRaisesRegex(ValueError, 'leaves its source tree'):
            self.capture()

    def test_missing_or_other_project_marker_is_not_a_checkout(self):
        self.initialize()
        marker = self.root / 'pyproject.toml'
        marker.unlink()
        self.assertIsNone(self.capture()['git']['commit'])
        marker.write_text('[project]\nname="another-project"\n')
        self.assertIsNone(self.capture()['git']['commit'])

    def test_empty_package_has_no_attributable_source(self):
        for path in self.package.glob('*.py'):
            path.unlink()
        with self.assertRaisesRegex(ValueError, 'needs source files'):
            self.capture()


class TestBuildAndRecipeOrigins(ProvenanceFixture):
    def test_installed_build_keeps_original_commit_and_rehashes_current_files(self):
        revision = self.initialize()
        baseline = self.capture()
        package = self.installed(baseline)
        unchanged = self.capture(package)
        self.assertEqual(unchanged['git']['commit'], revision)
        self.assertEqual(unchanged['origin'], 'build')
        self.assertEqual(unchanged['changes_since_build'], [])
        (package / 'alpha.py').unlink()
        (package / '__init__.py').write_text('VERSION = "modified"\n')
        (package / 'new.py').write_text('VALUE = 3\n')
        result = self.capture(package)
        self.assertEqual(result['git'], baseline['git'])
        self.assertEqual(result['changes_since_build'], [
            {'path': 'src/zipper_algebra/__init__.py', 'status': 'modified'},
            {'path': 'src/zipper_algebra/alpha.py', 'status': 'deleted'},
            {'path': 'src/zipper_algebra/new.py', 'status': 'added'},
        ])
        self.assertNotEqual(result['source_sha256'], baseline['source_sha256'])

    def test_own_unborn_checkout_never_reuses_embedded_old_commit(self):
        self.initialize()
        baseline = self.capture()
        (self.package / '_build_provenance.json').write_text(json.dumps(
            {'schema': source.BUILD_SCHEMA, 'origin': baseline}))
        self.remove_git()
        self.initialize(commit=False)
        self.assertEqual(self.capture()['origin'], 'unversioned')
        self.assertIsNone(self.capture()['git']['commit'])

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
        self.assertEqual(result['git'], baseline['git'])
        self.assertEqual(result['changes_since_build'], [
            {'path': 'MANIFEST.in', 'status': 'modified'},
            {'path': 'experiments/added.py', 'status': 'added'},
            {'path': 'experiments/changed.py', 'status': 'modified'},
            {'path': 'experiments/deleted.py', 'status': 'deleted'},
            {'path': 'setup.py', 'status': 'deleted'},
        ])
        self.assertNotIn('experiments/unused.py', result['source_files'])
        self.assertEqual(result['extra_source_git'], [{
            'paths': ['experiments/added.py', 'experiments/changed.py', 'experiments/deleted.py'],
            'git': {'commit': None, 'dirty': None, 'changes': []},
        }])

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
        self.assertEqual(result['git'], baseline['git'])
        self.assertEqual(result['changes_since_build'], [{'path': 'experiments/late.py', 'status': 'added'}])
        self.assertEqual(result['extra_source_git'], [{
            'paths': ['experiments/late.py'], 'git': {'commit': None, 'dirty': None, 'changes': []},
        }])
        (self.package / '_build_provenance.json').unlink()
        unversioned = self.capture(extra_sources=[recipe])
        self.assertEqual(unversioned['origin'], 'unversioned')
        self.assertEqual(unversioned['extra_source_git'][0]['git']['dirty'], None)

    def test_invalid_embedded_build_metadata_is_rejected(self):
        self.initialize()
        package = self.installed(self.capture())
        path = package / '_build_provenance.json'
        for bad in ([], {}, {'schema': 'other', 'origin': {}},
                    {'schema': source.BUILD_SCHEMA, 'origin': {}}):
            with self.subTest(bad=bad):
                path.write_text(json.dumps(bad))
                with self.assertRaises(ValueError):
                    self.capture(package)
        path.write_text('not JSON')
        with self.assertRaises(ValueError):
            self.capture(package)

    def test_declared_recipe_sources_are_hashed_and_track_untracked_files(self):
        self.initialize()
        recipe = self.root / 'experiments/recipe.py'
        recipe.parent.mkdir()
        recipe.write_text('RESULT = 1\n')
        result = self.capture(extra_sources=[recipe, recipe])
        self.assertEqual(result['source_files']['experiments/recipe.py'], source.file_identity(recipe))
        self.assertEqual(result['extra_source_git'], [])
        self.assertIn('experiments/recipe.py', [item['path'] for item in result['git']['changes']])

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
        self.assertEqual(result['git'], baseline['git'])
        self.assertEqual(result['extra_source_git'], [
            {'paths': ['experiments/recipe.py'],
             'git': {'commit': recipe_revision, 'dirty': False, 'changes': []}}])

    def test_additional_sources_are_confined_to_declared_public_python_files(self):
        private = self.root / 'private.py'
        private.write_text('PRIVATE = 1\n')
        outside = self.base / 'outside.py'
        outside.write_text('VALUE = 1\n')
        wrong_type = self.root / 'experiments/data.txt'
        wrong_type.parent.mkdir()
        wrong_type.write_text('not Python')
        for path in (private, outside, wrong_type, self.root / 'missing.py', wrong_type.parent / 'missing.py'):
            with self.subTest(path=path):
                with self.assertRaises(ValueError):
                    self.capture(extra_sources=[path])
        baseline = self.capture()
        package = self.installed(baseline)
        with self.assertRaisesRegex(ValueError, 'outside the declared project'):
            self.capture(package, extra_sources=[outside])
        a, _ = self.project(self.base / 'a')
        b, _ = self.project(self.base / 'b')
        paths = []
        for root in (a, b):
            (root / 'experiments').mkdir()
            path = root / 'experiments/recipe.py'
            path.write_text('VALUE = 1\n')
            paths.append(path)
        with self.assertRaisesRegex(ValueError, 'one project'):
            self.capture(package, extra_sources=paths)


class TestMetadataValidation(ProvenanceFixture):
    def test_invalid_schema_git_source_and_build_change_metadata(self):
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
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    source.validate_provenance(bad)

    def test_invalid_source_paths_and_hash_shapes(self):
        value = self.capture()
        for name in ('/absolute.py', '../parent.py', 'x\\y.py', 'drive:x.py', 'bad\n.py', '', 'data.txt', 'private.py', 'private.toml'):
            bad = deepcopy(value)
            bad['source_files'] = {name: {'bytes': 0, 'sha256': '0' * 64}}
            with self.subTest(name=name), self.assertRaises(ValueError):
                source.validate_provenance(bad)
        for identity in ({}, {'bytes': True, 'sha256': '0' * 64},
                         {'bytes': -1, 'sha256': '0' * 64}, {'bytes': 0, 'sha256': 'g' * 64},
                         {'bytes': 0, 'sha256': None}):
            bad = deepcopy(value)
            bad['source_files'] = {'src/zipper_algebra/a.py': identity}
            with self.subTest(identity=identity), self.assertRaises(ValueError):
                source.validate_provenance(bad)
        with self.assertRaises(ValueError):
            source._relative(None)

    def test_invalid_change_metadata_and_extra_origin_references(self):
        value = self.capture()
        valid = {'path': 'README.md', 'index': 'M', 'worktree': ' ', 'bytes': 0, 'sha256': '0' * 64}
        for change in ({}, {**valid, 'extra': 1}, {**valid, 'index': 'MM'},
                       {**valid, 'sha256': 'bad'}, {**valid, 'previous_path': '../old'}):
            bad = deepcopy(value)
            bad['git']['changes'] = [change]
            with self.subTest(change=change), self.assertRaises(ValueError):
                source.validate_provenance(bad)
        for entry in ({'paths': ['missing.py'], 'git': value['git']},
                      {'paths': None, 'git': value['git']}):
            bad = deepcopy(value)
            bad['extra_source_git'] = [entry]
            with self.assertRaises(ValueError):
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
        self.assertEqual(record['producer'], producer)
        self.assertEqual(record['inputs'], [input_tag])
        self.assertEqual(record['outputs'], {'result.txt': source.file_identity(output)})
        self.assertEqual(json.loads(output.with_name('result.txt.provenance.json').read_text()), record)
        self.assertEqual(output.read_bytes(), b'certificate bytes\r\n')
        self.assertNotIn(str(self.base), json.dumps(record))
        self.assertNotIn('original input', json.dumps(record))
        second = source.write_output_provenance(output, producer=producer, inputs=[input_path])
        self.assertEqual(second['inputs'][0]['sha256'], source.file_identity(input_path)['sha256'])

    def test_bundle_and_implicit_capture(self):
        directory = self.base / 'outputs'
        directory.mkdir()
        (directory / 'a.txt').write_text('a')
        (directory / 'b.json').write_text('{}')
        with patch.object(source, '__file__', str(self.package / 'provenance.py')):
            record = source.write_bundle_provenance(directory, ['a.txt', 'b.json'])
            single = source.write_output_provenance(directory / 'a.txt')
        self.assertEqual(record['schema'], source.OUTPUT_SCHEMA)
        self.assertEqual(set(record['outputs']), {'a.txt', 'b.json'})
        self.assertEqual(record['producer'], single['producer'])
        self.assertEqual(json.loads((directory / 'provenance.json').read_text()), record)

    def test_unsafe_duplicate_missing_and_empty_outputs_do_not_publish(self):
        (self.base / 'ok.txt').write_text('ok')
        producer = self.capture()
        for names in ([], ['../outside'], ['/absolute'], [''], ['.'], ['..'],
                      ['provenance.json'], ['ok.txt', 'ok.txt'], ['missing.txt'], [None], ['windows\\file.txt'], ['drive:x']):
            with self.subTest(names=names), self.assertRaises((ValueError, FileNotFoundError)):
                source.write_bundle_provenance(self.base, names, producer=producer)
        self.assertFalse((self.base / 'provenance.json').exists())

    def test_bad_prehashed_inputs_do_not_replace_existing_metadata(self):
        output = self.base / 'result.txt'
        output.write_text('result')
        producer = self.capture()
        original = source.write_output_provenance(output, producer=producer)
        for tag in ({}, {'name': '../x', 'bytes': 0, 'sha256': '0' * 64},
                    {'name': 'x', 'bytes': True, 'sha256': '0' * 64},
                    {'name': '..', 'bytes': 0, 'sha256': '0' * 64},
                    {'name': 'C:\\private', 'bytes': 0, 'sha256': '0' * 64},
                    {'name': None, 'bytes': 0, 'sha256': '0' * 64}):
            with self.subTest(tag=tag), self.assertRaises(ValueError):
                source.write_output_provenance(output, producer=producer, inputs=[tag])
        self.assertEqual(json.loads(output.with_name('result.txt.provenance.json').read_text()), original)

    def test_failed_atomic_publication_preserves_previous_sidecar(self):
        output = self.base / 'result.txt'
        output.write_text('result')
        producer = self.capture()
        source.write_output_provenance(output, producer=producer)
        sidecar = output.with_name('result.txt.provenance.json')
        previous = sidecar.read_bytes()
        with patch.object(source.os, 'replace', side_effect=OSError('interrupted')):
            with self.assertRaises(OSError):
                source.write_output_provenance(output, producer=producer)
        self.assertEqual(sidecar.read_bytes(), previous)
        self.assertEqual(sorted(p.name for p in self.base.iterdir()),
                         ['project', 'result.txt', 'result.txt.provenance.json'])


if __name__ == '__main__':
    unittest.main()
