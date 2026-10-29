"""Producer identities survive discovery, cache reuse and artifact migration."""
from copy import deepcopy
import json
from pathlib import Path
from unittest.mock import patch

from tests.test_discovery_pipeline import DiscoveryFixture
from zipper_algebra.cofactor import artifacts as mathematical
from zipper_algebra.execution import artifacts, pipeline, runner
from zipper_algebra.execution import store
from zipper_algebra.provenance import capture_provenance, file_identity, validate_provenance


class TestProducedProvenance(DiscoveryFixture):
    def test_discovery_binds_one_source_snapshot_to_every_product_and_sample(self):
        result = self.discover()
        self.assertEqual(result['status'], 'verified')
        producer = validate_provenance(result['producer'])
        self.assertIn('src/zipper_algebra/provenance.py', producer['source_files'])
        self.assertEqual(len(result['files']), 6)
        directory = Path(result['certificate']).parent
        self.assertEqual(set(path.name for path in directory.iterdir()), set(result['files']))
        self.assertFalse(list(directory.glob('*.provenance.json')))
        bundle = json.loads((directory / 'provenance.json').read_text())
        self.assertEqual(bundle['producer'], producer)
        self.assertEqual(set(bundle['outputs']), set(result['files']) - {'provenance.json'})
        for name, identity in bundle['outputs'].items():
            self.assertEqual(identity, file_identity(directory / name))
        for name in ('lift.json', 'skeleton.json', 'certificate.json'):
            document = json.loads((directory / name).read_text())
            self.assertEqual(document['schema'], mathematical.SCHEMA)
            self.assertEqual(document['provenance'], producer)
        for run in result['runs']:
            _, record = store.read(run['store'], run['key'], generation=run['generation'])
            self.assertEqual(record['request']['schema'], 'zipper-solver-request-v2')
            self.assertEqual(record['request']['producer'], producer)
            self.assertEqual(run['request']['producer'], producer)

    def test_cached_sample_retains_its_producer_and_new_source_state_has_a_new_key(self):
        first = runner.run_point(self.pres, self.target, (5, 7, 0), self.config)
        again = runner.run_point(self.pres, self.target, (5, 7, 0), self.config)
        self.assertTrue(again.hit)
        self.assertEqual(first.record, again.record)
        other = deepcopy(first.record['request']['producer'])
        # A different commit with byte-identical sources is a distinct declared
        # revision. It must not silently reuse a differently attributed record.
        other['git']['commit'] = 'f' * 40
        changed = runner.run_point(self.pres, self.target, (5, 7, 0), self.config,
                                   producer=other)
        self.assertNotEqual(first.key, changed.key)
        self.assertFalse(changed.hit)
        self.assertEqual(changed.record['request']['producer'], other)

    def test_changed_sources_cannot_publish_or_leave_a_stale_verified_selection(self):
        previous = self.discover()
        before = {name: Path(path).read_bytes() for name, path in previous['artifacts'].items()}
        with patch.object(artifacts, 'current_source_files', return_value={}):
            failed = self.discover()
        self.assertEqual(failed['status'], 'inconclusive')
        self.assertEqual(failed['execution_status'], 'failed')
        self.assertIn('source files changed', failed['reason'])
        self.assertNotIn('certificate', failed)
        for name, path in previous['artifacts'].items():
            self.assertEqual(Path(path).read_bytes(), before[name])

    def test_unreadable_recorded_sources_are_operational_publication_failures(self):
        producer = capture_provenance()
        for error in (PermissionError('unreadable source'), ValueError('source symlink changed')):
            with patch.object(artifacts, 'current_source_files', side_effect=error):
                with self.assertRaisesRegex(artifacts.ArtifactError, 'source files changed'):
                    artifacts.validate_source_snapshot(producer)

    def test_publication_rechecks_sources_after_solver_and_before_selecting_records(self):
        first = runner.run_point(self.pres, self.target, (5, 7, 0), self.config)
        before = store.read(first.store, first.key)
        with patch.object(artifacts, 'current_source_files', return_value={}):
            with self.assertRaisesRegex(artifacts.ArtifactError, 'source files changed'):
                store.publish(first.store, first.key, first.files,
                              first.record['request'], first.record['basis_entries'])
        self.assertEqual(store.read(first.store, first.key), before)

    def test_proof_validation_normalizes_missing_and_nonbyte_store_artifacts(self):
        from zipper_algebra.coefficients import NumericRing
        sample = runner.run_point(self.pres, self.target, (5, 7, 0), self.config)
        names, relations, _, _ = self.pres.build(NumericRing(5, 7, 0))
        for malformed in ({'basis.txt':sample.files['basis.txt']},
                          {**sample.files, 'proof.txt':'not raw bytes'}):
            with self.subTest(files=malformed), self.assertRaises(artifacts.ArtifactError):
                artifacts.validate_proofs(malformed, relations, names)

    def test_export_command_restores_exact_sample_bytes_with_their_original_producer(self):
        import io
        from contextlib import redirect_stdout
        from zipper_algebra.cli import main
        sample = runner.run_point(self.pres, self.target, (5, 7, 0), self.config)
        output = self.root / 'exported'
        for machine_readable in (False, True):
            stream = io.StringIO()
            arguments = ['export-sample', '--store', str(sample.store), '--key', sample.key,
                         '--generation', sample.record['generation'], '--output', str(output)]
            with redirect_stdout(stream):
                self.assertEqual(main(arguments + (['--json'] if machine_readable else [])), 0)
            if machine_readable:
                report = json.loads(stream.getvalue())
                self.assertEqual(report['original_producer'], sample.record['request']['producer'])
                validate_provenance(report['producer'])
            else:
                self.assertIn('Exported sample', stream.getvalue())
        for name, content in sample.files.items():
            self.assertEqual((output / name).read_bytes(), content)
        self.assertEqual(json.loads((output / 'record.json').read_text()), sample.record)

    def test_verified_generation_rechecks_sources_before_publication(self):
        result = self.discover()
        directory = Path(result['certificate']).parent
        order = {'solver': 'degree-lex', 'reduction': pipeline.REDUCTION_ORDER}
        skeleton = mathematical.load_artifact(directory / 'skeleton.json', self.pres, self.target,
                                              kind='skeleton', order=order)
        fitted = mathematical.load_artifact(directory / 'lift.json', self.pres, self.target,
                                            kind='lift', order=order)
        certificate = mathematical.load_artifact(directory / 'certificate.json', self.pres, self.target,
                                                 kind='cofactor')
        from zipper_algebra.cofactor.compose import load_multipliers
        multipliers = load_multipliers(directory / 'multipliers.json', skeleton, fitted.atoms, self.pres)
        with patch.object(artifacts, 'current_source_files', return_value={}):
            with self.assertRaisesRegex(artifacts.ArtifactError, 'source files changed'):
                pipeline._publish_verified(self.root / 'other', self.pres, self.target, skeleton,
                                           fitted, multipliers, certificate, order)
        self.assertFalse(any(path.name == 'provenance.json'
                             for path in (self.root / 'other').rglob('*')))

    def test_retired_schema_and_invalid_current_source_metadata_are_rejected(self):
        result = self.discover()
        document = json.loads(Path(result['artifacts']['certificate.json']).read_text())
        old = self.root / 'old.json'
        legacy = deepcopy(document)
        legacy['schema'] = 'zipper-mathematical-artifact-v1'
        del legacy['provenance']
        old.write_text(json.dumps(legacy))
        with self.assertRaisesRegex(ValueError, 'schema'):
            mathematical.load_artifact(old, self.pres, self.target)
        document['provenance']['source_sha256'] = '0' * 64
        old.write_text(json.dumps(document))
        with self.assertRaisesRegex(ValueError, 'digest'):
            mathematical.load_artifact(old, self.pres, self.target)
