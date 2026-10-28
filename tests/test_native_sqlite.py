"""Actual process-boundary tests of the optional native database route."""
from dataclasses import replace
import json
from pathlib import Path
import sqlite3
import sys
import tempfile
import unittest
from unittest.mock import patch

from zipper_algebra import Presentation
from zipper_algebra.execution import ArtifactError, ExecutionError, RunConfig, run_point
from zipper_algebra.execution import native_io, pipeline, store
from zipper_algebra.execution.config import DEFAULT_FLAGS
from zipper_algebra.cofactor.atoms import AtomSet


PRODUCER = r'''
import argparse
import json
from pathlib import Path
import re
import sqlite3
import sys
import time

p = argparse.ArgumentParser()
p.add_argument('--sqlite')
p.add_argument('--sqlite-job')
p.add_argument('-i')
p.add_argument('-o')
p.add_argument('-p')
p.add_argument('-m', type=int)
p.add_argument('-d', type=int)
p.add_argument('--mode', default='success')
p.add_argument('--counter')
a, _ = p.parse_known_args()
if a.sqlite:
    assert a.sqlite_job and not any((a.i, a.o, a.p))
    db = sqlite3.connect(a.sqlite)
    original = db.execute('SELECT input FROM f4ncgb_io WHERE job_id=?',
                          (a.sqlite_job,)).fetchone()[0]
    db.close()
else:
    assert all((a.i, a.o, a.p))
    original = Path(a.i).read_bytes()
if a.counter:
    counter = Path(a.counter)
    counter.write_text(str(int(counter.read_text()) + 1) if counter.exists() else '1')
if a.mode == 'timeout':
    time.sleep(10)
if a.mode == 'exit':
    print('synthetic solver failure', file=sys.stderr)
    sys.exit(7)
relations = [line.rstrip(',') for line in original.decode().splitlines()[2:]]
spaced = [re.sub(r'(?<!^)([+-])', r' \1 ', line) for line in relations]
basis = (spaced[0] + (' + 1' if a.mode == 'false-proof' else '') + '\n').encode()
proof = b'[i0]\n'
if a.mode == 'all':
    basis = ('\n'.join(spaced) + '\n').encode()
    proof = ('\n'.join('[i%d]' % i for i in range(len(spaced))) + '\n').encode()
if a.sqlite:
    db = sqlite3.connect(a.sqlite)
    if a.mode == 'delete-row':
        db.execute('DELETE FROM f4ncgb_io WHERE job_id=?', (a.sqlite_job,))
    elif a.mode not in ('missing', 'delete-db'):
        if a.mode == 'text-basis': basis = basis.decode()
        if a.mode == 'text-proof': proof = proof.decode()
        if a.mode == 'partial': proof = None
        if a.mode == 'bad-bytes': basis = b'\xff'
        incoming = (b'changed input' if a.mode == 'mutate-input' else
                    original.decode() if a.mode == 'text-input' else original)
        db.execute("UPDATE f4ncgb_io SET input=?,basis=?,proof=?,state='complete' WHERE job_id=?",
                   (incoming, basis, proof, a.sqlite_job))
        if a.mode == 'transaction-timeout':
            Path(str(a.counter) + '.transaction').write_text('uncommitted result')
            time.sleep(10)
    db.commit()
    db.close()
    if a.mode == 'delete-db': Path(a.sqlite).unlink()
else:
    Path(a.o).write_bytes(basis)
    Path(a.p).write_bytes(proof)
print(json.dumps({'iterations': a.m, 'degree': a.d}))
'''


class TestNativeSQLite(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix='zipper-native-')
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.script = self.root / 'producer.py'
        self.script.write_text(PRODUCER)
        self.counter = self.root / 'counter'
        self.pres = Presentation(2)
        self.point = (5, 7, 0)

    def configuration(self, mode='success', **kwargs):
        values = dict(binary=Path(sys.executable), output_root=self.root / 'output',
                      command_prefix=(self.script,), solver_io='sqlite', timeout_seconds=3,
                      max_iterations=13, max_degree=7,
                      flags=(*DEFAULT_FLAGS, '--mode', mode, '--counter', str(self.counter)),
                      solver_revision='synthetic SQLite process producer')
        values.update(kwargs)
        return RunConfig(**values)

    def execute(self, mode='success', **kwargs):
        return run_point(self.pres, 'X_2', self.point, self.configuration(mode, **kwargs))

    def assert_no_pending_jobs(self, path=None):
        path = path or self.root / 'output/samples.sqlite3'
        with sqlite3.connect(path) as connection:
            self.assertEqual(connection.execute('SELECT count(*) FROM f4ncgb_io').fetchone()[0], 0)

    def test_native_and_file_routes_retain_identical_bytes_with_distinct_request_modes(self):
        native = self.execute()
        file = run_point(self.pres, 'X_2', self.point, self.configuration(solver_io='files'))
        self.assertEqual(native.files, file.files)
        self.assertNotEqual(native.key, file.key)
        self.assertEqual(native.record['request']['solver']['io_mode'], 'sqlite')
        self.assertEqual(file.record['request']['solver']['io_mode'], 'files')
        self.assert_no_pending_jobs()

    def test_native_route_creates_no_solver_files_and_reuses_validated_cache(self):
        config = self.configuration(cache_root=self.root / 'cache')
        first = run_point(self.pres, 'X_2', self.point, config)
        self.assertFalse(first.hit)
        self.assertEqual(self.counter.read_text(), '1')
        self.assertEqual({path.name for path in config.output_root.iterdir()},
                         {'.gitignore', 'samples.sqlite3'})
        again = run_point(self.pres, 'X_2', self.point, config)
        self.assertTrue(again.hit)
        self.assertEqual(self.counter.read_text(), '1')
        shared = run_point(self.pres, 'X_2', self.point,
                           replace(config, output_root=self.root / 'another'))
        self.assertTrue(shared.hit)
        self.assertEqual(shared.files, first.files)
        self.assertEqual(self.counter.read_text(), '1')
        refreshed = run_point(self.pres, 'X_2', self.point, config, refresh=True)
        self.assertFalse(refreshed.hit)
        self.assertEqual(self.counter.read_text(), '2')
        self.assert_no_pending_jobs()

    def test_nonzero_exit_and_timeout_preserve_process_errors_and_remove_pending_jobs(self):
        for mode, kwargs, message in [('exit', {}, 'code 7'),
                                      ('timeout', {'timeout_seconds':0.1}, 'exceeded')]:
            with self.subTest(mode=mode), self.assertRaisesRegex(ExecutionError, message):
                self.execute(mode, **kwargs)
            self.assert_no_pending_jobs()
        self.assert_no_pending_jobs()

    def test_invalid_native_database_outputs_cannot_become_validated_samples(self):
        for mode, message in [('missing', 'complete result'), ('delete-row', 'disappeared'),
                              ('delete-db', 'disappeared'), ('mutate-input', 'changed'),
                              ('text-input', 'changed'), ('text-basis', 'exact bytes'),
                              ('text-proof', 'exact bytes'), ('partial', 'exact bytes'),
                              ('bad-bytes', 'invalid solver proof'),
                              ('false-proof', 'input-ideal identity')]:
            with self.subTest(mode=mode), self.assertRaisesRegex(ArtifactError, message):
                self.execute(mode)
            self.assert_no_pending_jobs()
            with sqlite3.connect(self.root / 'output/samples.sqlite3') as connection:
                self.assertEqual(connection.execute('SELECT count(*) FROM selections').fetchone()[0], 0)

    def test_timeout_rolls_back_a_native_transaction_and_the_store_remains_resumable(self):
        with self.assertRaisesRegex(ExecutionError, 'exceeded'):
            self.execute('transaction-timeout', timeout_seconds=2)
        self.assertTrue(self.counter.with_name(self.counter.name + '.transaction').is_file())
        self.assert_no_pending_jobs()
        resumed = self.execute()
        self.assertFalse(resumed.hit)
        self.assert_no_pending_jobs()

    def test_configuration_rejects_unknown_routes_and_native_option_hijacking(self):
        for mode in ('unknown', None, True, []):
            with self.subTest(mode=mode), self.assertRaisesRegex(ValueError, 'solver_io'):
                self.configuration(solver_io=mode)
        for option in ('--sqlite', '--sqlite=elsewhere', '--sqlite-job=other', '--features'):
            with self.subTest(option=option), self.assertRaisesRegex(ValueError, 'owned by RunConfig'):
                self.configuration(flags=(option,))

    def test_native_input_rejects_nonbyte_values_before_creating_any_database(self):
        config = self.configuration()
        path = self.root / 'other.sqlite3'
        with self.assertRaisesRegex(ArtifactError, 'exact bytes'):
            native_io.execute(config, 'not bytes', path)
        self.assertFalse(path.exists())

    def test_start_failure_removes_its_pending_row(self):
        binary = self.root / 'not-an-executable'
        binary.write_bytes(b'not an executable')
        config = self.configuration(binary=binary, command_prefix=())
        with self.assertRaisesRegex(native_io.NativeExecutionError, 'could not start'):
            native_io.execute(config, b'x\n0\nx\n', config.output_root / 'samples.sqlite3')
        self.assert_no_pending_jobs()

    def test_cleanup_failure_does_not_hide_the_original_process_error(self):
        with patch.object(native_io, '_discard', side_effect=ArtifactError('cleanup unavailable')):
            with self.assertRaisesRegex(ExecutionError, 'code 7'):
                self.execute('exit')

    def test_cleanup_failure_after_success_is_not_reported_as_a_validated_sample(self):
        with patch.object(native_io, '_discard', side_effect=ArtifactError('cleanup unavailable')):
            with self.assertRaisesRegex(ArtifactError, 'cleanup unavailable'):
                self.execute()
        with sqlite3.connect(self.root / 'output/samples.sqlite3') as connection:
            self.assertEqual(connection.execute('SELECT count(*) FROM selections').fetchone()[0], 0)

    def test_full_discovery_and_export_use_the_same_native_byte_validation(self):
        pres = Presentation(2, zero_from=2)
        report = pipeline.run(pres, 'X_2', self.point, self.configuration('all'),
                              start=7, step=2, maxpts=9, atoms=AtomSet([]), verbose=False)
        self.assertEqual(report['status'], 'verified')
        self.assert_no_pending_jobs()
        self.assertEqual(len(report['runs']), 7)
        sample = report['runs'][0]
        files, record = store.read(sample['store'], sample['key'], sample['generation'])
        output = self.root / 'export'
        self.assertEqual(store.export(sample['store'], sample['key'], output), record)
        for name, content in files.items():
            self.assertEqual((output / name).read_bytes(), content)
        self.assertEqual(json.loads((output / 'record.json').read_text()), record)


if __name__ == '__main__':
    unittest.main()
