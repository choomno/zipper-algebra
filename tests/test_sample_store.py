"""Byte identity, transactions, corruption and export of consolidated samples."""
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from copy import deepcopy
import json
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from threading import Event
import unittest
from unittest.mock import patch

from zipper_algebra.execution import store
from zipper_algebra.execution.artifacts import ArtifactError, canonical_json, digest_bytes


class TestSampleStore(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.database = self.root / 'samples.sqlite3'
        sources = {'src/zipper_algebra/example.py': {'bytes': 3, 'sha256': digest_bytes(b'abc')}}
        self.producer = {'schema': 'zipper-source-provenance-v1', 'origin': 'unversioned',
                         'git': {'commit': None, 'dirty': None, 'changes': []},
                         'source_files': sources, 'source_sha256': digest_bytes(canonical_json(sources)),
                         'changes_since_build': [], 'extra_source_git': [], 'omitted_source_count': 0}
        self.files = {'input.ms': b'x\n0\nx\n', 'basis.txt': b'x\n', 'proof.txt': b'[i0]\n',
                      'stdout.txt': b'completed\r\n', 'stderr.txt': b''}
        self.request = {'schema': 'zipper-solver-request-v2', 'producer': self.producer,
                        'input_sha256': digest_bytes(self.files['input.ms']), 'point': [5, 7, 0]}
        self.key = digest_bytes(canonical_json(self.request))
        # Source acquisition/guard behaviour is tested separately; this unit
        # exercises storage of a worker's already captured source snapshot.
        self.guard = self.enterContext(patch.object(store, 'validate_source_snapshot'))

    def publish(self, files=None, **kwargs):
        return store.publish(kwargs.pop('path', self.database), kwargs.pop('key', self.key),
                             self.files if files is None else files,
                             kwargs.pop('request', self.request), kwargs.pop('basis_entries', 1), **kwargs)

    def sql(self, statement, arguments=()):
        with closing(sqlite3.connect(self.database)) as connection:
            with connection:
                return connection.execute(statement, arguments).fetchall()

    def record_bytes(self, content):
        self.sql('UPDATE generations SET record=?', (content,))

    def test_missing_and_initializing_databases_are_cache_misses_without_side_effects(self):
        self.assertIsNone(store.load(self.database, self.key, self.request))
        self.assertFalse(self.database.exists())
        with self.assertRaisesRegex(ArtifactError, 'not present'):
            store.read(self.database, self.key)
        self.database.touch()
        self.assertIsNone(store.load(self.database, self.key, self.request))
        self.publish()
        absent = 'a' * 64
        request = dict(self.request, point=[7, 11, 0])
        key = digest_bytes(canonical_json(request))
        self.assertIsNone(store.load(self.database, key, request))
        with self.assertRaisesRegex(ArtifactError, 'not present'):
            store.read(self.database, absent)

    def test_roundtrip_is_byte_exact_and_one_database_remains_after_close(self):
        record = self.publish()
        self.guard.assert_called_once_with(self.producer)
        files, loaded = store.load(self.database, self.key, self.request)
        self.assertEqual(files, self.files)
        self.assertEqual(loaded, record)
        self.assertEqual(store.read(self.database, self.key), (files, record))
        self.assertEqual(record['request']['producer'], self.producer)
        self.assertFalse(record['basis_complete'])
        self.assertEqual(set(p.name for p in self.root.iterdir()), {'.gitignore', 'samples.sqlite3'})
        before = self.database.read_bytes()
        store.read(self.database, self.key)
        self.assertEqual(self.database.read_bytes(), before)

    def test_read_and_export_work_on_a_read_only_copied_database(self):
        record = self.publish()
        copy = self.root / 'copied.sqlite3'
        shutil.copyfile(self.database, copy)
        copy.chmod(0o444)
        self.addCleanup(copy.chmod, 0o600)
        self.assertEqual(store.read(copy, self.key), (self.files, record))
        self.assertEqual(store.export(copy, self.key, self.root / 'exported'), record)
        self.assertEqual((self.root / 'exported' / 'input.ms').read_bytes(), self.files['input.ms'])

    def test_schema_and_unrelated_database_rejection_preserve_existing_content(self):
        self.database.write_bytes(b'not a SQLite database')
        with self.assertRaisesRegex(ArtifactError, 'database'):
            store.load(self.database, self.key, self.request)
        self.assertEqual(self.database.read_bytes(), b'not a SQLite database')
        self.database.unlink()
        self.sql('CREATE TABLE unrelated (payload TEXT)')
        with self.assertRaisesRegex(ArtifactError, 'supported'):
            self.publish()
        self.assertEqual(self.sql("SELECT name FROM sqlite_master WHERE type='table'"), [('unrelated',)])

    def test_empty_database_with_foreign_marker_is_not_initialized(self):
        for marker in ('application_id', 'user_version'):
            with self.subTest(marker=marker):
                self.database.unlink(missing_ok=True)
                self.sql(f'PRAGMA {marker}=99')
                with self.assertRaisesRegex(ArtifactError, 'unrelated schema marker'):
                    self.publish()

    def test_directory_instead_of_database_is_rejected_without_changes(self):
        self.database.mkdir()
        with self.assertRaisesRegex(ArtifactError, 'database'):
            self.publish()
        self.assertTrue(self.database.is_dir())

    def test_database_markers_objects_and_table_definitions_are_checked(self):
        corruptions = [
            'PRAGMA application_id=0', 'PRAGMA user_version=99',
            'CREATE TABLE extra (value TEXT)',
            'CREATE TABLE sqliteXextra (value TEXT)',
            'CREATE VIEW sqliteXview AS SELECT * FROM metadata',
            'CREATE TRIGGER sqliteXtrigger BEFORE INSERT ON selections BEGIN SELECT 1; END',
            'ALTER TABLE metadata RENAME COLUMN schema TO other',
            "UPDATE metadata SET schema='other'",
            'INSERT INTO metadata SELECT * FROM metadata',
        ]
        for command in corruptions:
            with self.subTest(command=command):
                self.database.unlink(missing_ok=True)
                self.publish()
                self.sql(command)
                with self.assertRaises(ArtifactError):
                    store.read(self.database, self.key)

    def test_corrupt_artifact_content_types_names_and_rows_are_rejected(self):
        corruptions = [
            ("UPDATE artifacts SET content=? WHERE name='basis.txt'", (b'changed',)),
            ("UPDATE artifacts SET content=? WHERE name='stderr.txt'", ('text not a blob',)),
            ("DELETE FROM artifacts WHERE name='proof.txt'", ()),
            ("UPDATE artifacts SET name='unexpected.txt' WHERE name='stderr.txt'", ()),
        ]
        for command, arguments in corruptions:
            with self.subTest(command=command):
                self.database.unlink(missing_ok=True)
                self.publish()
                self.sql(command, arguments)
                with self.assertRaisesRegex(ArtifactError, 'invalid stored sample'):
                    store.load(self.database, self.key, self.request)

    def test_corrupt_record_fields_and_non_json_record_are_rejected(self):
        record = self.publish()
        mutations = [lambda data: data.update(extra=True), lambda data: data.update(schema='other'),
                     lambda data: data.update(status='complete'), lambda data: data.update(basis_complete=True),
                     lambda data: data.update(basis_entries=True), lambda data: data.update(basis_entries=0),
                     lambda data: data.update(generation='0' * 64), lambda data: data.update(request={}),
                     lambda data: data['request']['producer'].update(source_sha256='0' * 64),
                     lambda data: data.update(files={})]
        for mutate in mutations:
            data = deepcopy(record)
            mutate(data)
            self.record_bytes(canonical_json(data))
            with self.subTest(mutate=mutate), self.assertRaisesRegex(ArtifactError, 'invalid stored sample'):
                store.read(self.database, self.key)
        for content in (b'not JSON', b'[]', 'text instead of blob'):
            self.record_bytes(content)
            with self.subTest(content=content), self.assertRaisesRegex(ArtifactError, 'invalid stored sample'):
                store.read(self.database, self.key)

    def test_selected_generation_and_request_key_are_checked(self):
        self.publish()
        self.sql("UPDATE selections SET generation='broken'")
        with self.assertRaisesRegex(ArtifactError, 'generation'):
            store.read(self.database, self.key)
        self.sql('UPDATE selections SET generation=?', ('0' * 64,))
        with self.assertRaisesRegex(ArtifactError, 'generation is missing'):
            store.read(self.database, self.key)
        self.publish()
        data = store.read(self.database, self.key)[1]
        data['request']['point'] = [7, 13, 0]
        self.record_bytes(canonical_json(data))
        with self.assertRaisesRegex(ArtifactError, 'identity'):
            store.read(self.database, self.key)

    def test_request_keys_input_bytes_counts_and_metadata_are_validated_before_publication(self):
        for key in ('../unsafe', '0' * 64, True):
            with self.subTest(key=key), self.assertRaises(ArtifactError):
                self.publish(key=key)
        for count in (0, True, 1.5):
            with self.subTest(count=count), self.assertRaises(ArtifactError):
                self.publish(basis_entries=count)
        for request in ({}, dict(self.request, producer={}), dict(self.request, point=float('nan'))):
            with self.subTest(request=request), self.assertRaises(ArtifactError):
                self.publish(request=request)
        for files in ([], {}, dict(self.files, **{'basis.txt': 'not bytes'}),
                      dict(self.files, **{'input.ms': b'another presentation'})):
            with self.subTest(files=files), self.assertRaises(ArtifactError):
                self.publish(files=files)
        with self.assertRaises(ArtifactError):
            store.load(self.database, 'a' * 64, self.request)
        with self.assertRaises(ArtifactError):
            store.read(self.database, self.key, 'bad')

    def test_first_failed_publication_leaves_a_resumable_empty_store(self):
        self.guard.side_effect = ArtifactError('source changed during computation')
        with self.assertRaisesRegex(ArtifactError, 'source changed'):
            self.publish()
        self.guard.side_effect = None
        self.assertIsNone(store.load(self.database, self.key, self.request))
        self.assertEqual(self.sql('SELECT count(*) FROM generations'), [(0,)])
        self.publish()
        self.assertIsNotNone(store.load(self.database, self.key, self.request))

    def test_interrupted_refresh_rolls_back_bytes_history_and_selection(self):
        first = self.publish()
        changed = dict(self.files, **{'stdout.txt': b'another valid run'})
        self.guard.side_effect = ArtifactError('source changed during computation')
        with self.assertRaises(ArtifactError):
            self.publish(changed)
        self.guard.side_effect = None
        self.assertEqual(store.read(self.database, self.key), (self.files, first))
        self.assertEqual(self.sql('SELECT count(*) FROM generations'), [(1,)])
        self.assertEqual(self.sql('SELECT count(*) FROM artifacts'), [(5,)])

    def test_refresh_retains_history_and_repairs_corrupt_same_generation(self):
        first = self.publish()
        changed = dict(self.files, **{'stdout.txt': b'another valid run'})
        second = self.publish(changed)
        self.assertNotEqual(first['generation'], second['generation'])
        self.assertEqual(store.read(self.database, self.key), (changed, second))
        self.assertEqual(store.read(self.database, self.key, first['generation']), (self.files, first))
        self.sql("UPDATE artifacts SET content=? WHERE generation=? AND name='basis.txt'",
                 (b'corrupt', second['generation']))
        self.assertEqual(self.publish(changed), second)
        self.assertEqual(store.read(self.database, self.key), (changed, second))

    def test_concurrent_identical_and_different_generations_publish_without_partial_records(self):
        with ThreadPoolExecutor(max_workers=8) as executor:
            records = list(executor.map(lambda _: self.publish(), range(8)))
        self.assertTrue(all(record == records[0] for record in records))
        self.assertEqual(self.sql('SELECT count(*) FROM generations'), [(1,)])
        variants = [dict(self.files, **{'stdout.txt': str(i).encode()}) for i in range(6)]
        with ThreadPoolExecutor(max_workers=6) as executor:
            records = list(executor.map(self.publish, variants))
        for files, record in zip(variants, records):
            self.assertEqual(store.read(self.database, self.key, record['generation']), (files, record))
        current, record = store.read(self.database, self.key)
        self.assertIn(current, variants)
        self.assertIn(record, records)
        self.assertFalse(self.database.with_name('samples.sqlite3-journal').exists())

    def test_reader_observes_committed_selection_while_refresh_is_uncommitted(self):
        first = self.publish()
        changed = dict(self.files, **{'stdout.txt': b'new'})
        written, release = Event(), Event()
        def pause(_producer):
            written.set()
            self.assertTrue(release.wait(10))
        self.guard.side_effect = pause
        with ThreadPoolExecutor(max_workers=1) as executor:
            pending = executor.submit(self.publish, changed)
            self.assertTrue(written.wait(10))
            try:
                self.assertEqual(store.read(self.database, self.key), (self.files, first))
            finally:
                release.set()
            second = pending.result(timeout=10)
        self.assertEqual(store.read(self.database, self.key), (changed, second))

    def test_export_preserves_exact_original_provenance_and_reuses_identical_files(self):
        first = self.publish()
        output = self.root / 'export'
        self.assertEqual(store.export(self.database, self.key, output), first)
        self.assertEqual(store.export(self.database, self.key, output, first['generation']), first)
        for name, content in self.files.items():
            self.assertEqual((output / name).read_bytes(), content)
        self.assertEqual(json.loads((output / 'record.json').read_text()), first)
        self.assertEqual(set(path.name for path in output.iterdir()),
                         set(self.files) | {'record.json', '.gitignore'})

    def test_export_refuses_conflicting_files_directories_and_git_roots(self):
        self.publish()
        output = self.root / 'export'
        output.mkdir()
        blocked = output / 'proof.txt'
        blocked.write_bytes(b'precious')
        with self.assertRaisesRegex(ArtifactError, 'overwrite'):
            store.export(self.database, self.key, output)
        self.assertEqual(blocked.read_bytes(), b'precious')
        self.assertFalse((output / 'input.ms').exists())
        blocked.unlink()
        blocked.mkdir()
        with self.assertRaisesRegex(ArtifactError, 'overwrite'):
            store.export(self.database, self.key, output)
        checkout = self.root / 'checkout'
        (checkout / '.git').mkdir(parents=True)
        with self.assertRaisesRegex(ArtifactError, 'Git checkout'):
            store.export(self.database, self.key, checkout)
        self.assertFalse((checkout / '.gitignore').exists())

    def test_export_link_races_preserve_a_differing_winner_and_accept_identical_bytes(self):
        self.publish()
        link = store.os.link
        for equal in (True, False):
            output = self.root / str(equal)
            def raced(source, destination):
                Path(destination).write_bytes(Path(source).read_bytes() if equal else b'competing bytes')
                return link(source, destination)
            with patch.object(store.os, 'link', side_effect=raced):
                if equal:
                    store.export(self.database, self.key, output)
                else:
                    with self.assertRaisesRegex(ArtifactError, 'differing file'):
                        store.export(self.database, self.key, output)
            self.assertFalse(list(output.glob('.export-*')))

    def test_export_write_failure_cleans_temporary_files(self):
        self.publish()
        output = self.root / 'export'
        with patch.object(store.os, 'fsync', side_effect=OSError('disk failure')):
            with self.assertRaisesRegex(OSError, 'disk failure'):
                store.export(self.database, self.key, output)
        self.assertFalse(list(output.glob('.export-*')))
        self.assertFalse((output / 'input.ms').exists())

    def test_database_symlinks_are_rejected_without_following_them(self):
        self.publish()
        alias = self.root / 'alias.sqlite3'
        try:
            alias.symlink_to(self.database)
        except OSError:
            self.skipTest('creating symlinks is unavailable')
        with self.assertRaisesRegex(ArtifactError, 'symbolic link'):
            store.read(alias, self.key)
        with self.assertRaisesRegex(ArtifactError, 'symbolic link'):
            self.publish(path=alias)

    def test_native_transport_transaction_preserves_raw_blobs_and_read_only_snapshots(self):
        with store.transaction(self.database) as connection:
            self.assertIsNone(connection)
        self.assertFalse(self.database.exists())
        with store.transaction(self.database, write=True) as connection:
            connection.execute('INSERT INTO f4ncgb_io (job_id,input,state) VALUES (?,?,?)',
                               ('native-job', b'x\r\n\x00', 'ready'))
        with store.transaction(self.database) as connection:
            self.assertEqual(connection.execute('SELECT * FROM f4ncgb_io').fetchall(),
                             [('native-job', b'x\r\n\x00', None, None, 'ready')])
            with self.assertRaises(sqlite3.OperationalError):
                connection.execute('DELETE FROM f4ncgb_io')
        with store.transaction(self.database, write=True) as connection:
            connection.execute('UPDATE f4ncgb_io SET basis=?,proof=?,state=? WHERE job_id=?',
                               (b'x\n', b'[i0]\n', 'complete', 'native-job'))
        with store.transaction(self.database) as connection:
            self.assertEqual(connection.execute('SELECT basis,proof,state FROM f4ncgb_io').fetchall(),
                             [(b'x\n', b'[i0]\n', 'complete')])
        with store.transaction(self.database, write=True) as connection:
            connection.execute('DELETE FROM f4ncgb_io WHERE job_id=?', ('native-job',))
        self.assertFalse(self.database.with_name('samples.sqlite3-journal').exists())

    def test_native_transport_write_failure_rolls_back_its_job(self):
        with self.assertRaisesRegex(RuntimeError, 'interrupted'):
            with store.transaction(self.database, write=True) as connection:
                connection.execute('INSERT INTO f4ncgb_io (job_id,input,state) VALUES (?,?,?)',
                                   ('native-job', b'x\n', 'ready'))
                raise RuntimeError('interrupted native job enqueue')
        with store.transaction(self.database) as connection:
            self.assertEqual(connection.execute('SELECT * FROM f4ncgb_io').fetchall(), [])
        with self.assertRaises(ArtifactError):
            with store.transaction(self.database, write=True) as connection:
                connection.execute('INSERT INTO f4ncgb_io (job_id,input,state) VALUES (?,?,?)',
                                   ('native-job', b'x\n', 'unrecognized'))

    def crash_spilled_writer(self):
        self.publish()
        source = '''import sqlite3,sys,time
connection=sqlite3.connect(sys.argv[1],isolation_level=None)
connection.execute('PRAGMA cache_size=1')
connection.execute('PRAGMA cache_spill=ON')
connection.execute('BEGIN IMMEDIATE')
connection.execute("UPDATE artifacts SET content=? WHERE name='stdout.txt'",(b'uncommitted'*262144,))
print('spilled',flush=True)
time.sleep(60)
'''
        process = subprocess.Popen([sys.executable, '-u', '-c', source, str(self.database)],
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            self.assertEqual(process.stdout.readline().strip(), 'spilled')
            journal = self.database.with_name(self.database.name + '-journal')
            with journal.open('rb') as stream:
                self.assertEqual(stream.read(8), bytes.fromhex('d9d505f920a163d7'))
            self.assertGreater(self.database.stat().st_size, 1024 * 1024)
            process.kill()
            process.wait(timeout=10)
        finally:
            if process.poll() is None:
                process.kill()
                process.wait(timeout=10)
            process.stdout.close()
            process.stderr.close()
        return journal

    def test_killed_spilled_writer_recovers_original_bytes_and_can_resume(self):
        journal = self.crash_spilled_writer()
        recover = store._recover_hot_journal
        with patch.object(store, '_recover_hot_journal', wraps=recover) as recovery:
            files, record = store.load(self.database, self.key, self.request)
        recovery.assert_called_once_with(self.database)
        self.assertEqual(files, self.files)
        self.assertFalse(journal.exists())
        self.assertEqual(store.read(self.database, self.key)[1], record)
        changed = dict(self.files, **{'stdout.txt': b'valid resumed output'})
        self.publish(changed)
        self.assertEqual(store.read(self.database, self.key)[0], changed)

    def test_hot_journal_that_cannot_be_recovered_remains_an_explicit_failure(self):
        journal = self.crash_spilled_writer()
        connect = sqlite3.connect
        def forbid_recovery(database, *args, **kwargs):
            # Real SQLite read-only mode reproduces a genuine inability to
            # restore journaled pages, independently of the test user's UID.
            return connect(str(database).replace('?mode=rw', '?mode=ro'), *args, **kwargs)
        with patch.object(store.sqlite3, 'connect', side_effect=forbid_recovery):
            with self.assertRaisesRegex(ArtifactError, 'needs writable recovery'):
                store.load(self.database, self.key, self.request)
        self.assertTrue(journal.exists())
        self.assertEqual(store.load(self.database, self.key, self.request)[0], self.files)
