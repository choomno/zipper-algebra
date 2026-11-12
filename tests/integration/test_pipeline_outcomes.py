"""Bounded discovery outcomes preserve exact acceptance and operational failures."""
from __future__ import annotations

from contextlib import redirect_stdout
from dataclasses import replace
import io
import sys
from unittest.mock import patch

import pytest

from tests.test_discovery_pipeline import DiscoveryFixture
from zipper_algebra import cli
from zipper_algebra.coefficients import NumericRing
from zipper_algebra.execution import pipeline


class TestPipelineReconstructionOutcome(DiscoveryFixture):
    @pytest.mark.parametrize('error',[RuntimeError('degree budget exhausted'),
                                      ImportError('FiniteFlow bindings are not installed')])
    def test_failed_reconstruction_is_inconclusive_and_never_selects_old_certificate(self,error):
        prior=self.discover()
        certificate=prior['certificate']
        from pathlib import Path
        before=Path(certificate).read_bytes()
        with patch.object(pipeline,'reconstruct_cofactor',side_effect=error):
            result=self.discover()
        assert result['status']=='inconclusive'
        assert result['in_ideal'] is None and not result['verified']
        assert result['stage']=='reconstruction'
        assert result['reason']==str(error)
        assert 'certificate' not in result
        assert Path(certificate).read_bytes()==before

    def test_publication_independently_rejects_a_false_reconstructed_cofactor(self):
        from zipper_algebra.cofactor.cofactor_io import Cofactor
        _,arguments=self.material()
        discovery=arguments[3]
        false=Cofactor.from_presentation(self.pres,self.target,self.atoms,[],
                                         'false candidate',anchor=(5,7,0))
        with patch.object(pipeline,'reconstruct_cofactor',return_value=replace(discovery,cofactor=false)):
            result=self.discover()
        assert result['status']=='inconclusive'
        assert result['reason']=='candidate failed exact verification'
        assert 'certificate' not in result

    def test_native_dependency_changes_prevent_publication(self):
        with patch.object(pipeline,'reconstruction_provenance',side_effect=[{'version':1},{'version':2}]):
            result=self.discover()
        assert result['status']=='inconclusive'
        assert result['stage']=='publication'
        assert result['execution_status']=='failed'
        assert result['reason']=='reconstruction dependencies changed during discovery'
        assert 'certificate' not in result

    @pytest.mark.parametrize('limits',[{'max_contexts':1},{'max_products':1}])
    def test_support_resource_exhaustion_is_inconclusive(self,limits):
        result=self.discover(**limits)
        assert result['status']=='inconclusive'
        assert result['stage']=='support'
        assert 'budget exhausted' in result['reason']
        assert result['in_ideal'] is None

    def test_programmer_errors_are_not_reported_as_computational_inconclusiveness(self):
        with patch.object(pipeline,'reconstruct_cofactor',side_effect=TypeError('wrong interface')):
            with pytest.raises(TypeError,match='wrong interface'):
                self.discover()

    def test_cli_returns_nonzero_for_inconclusive_reconstruction(self):
        with patch.object(pipeline,'run',return_value={'status':'inconclusive','verified':False,
                                                     'reason':'degree budget'}),redirect_stdout(io.StringIO()):
            assert cli.main(['discover','--n','2','--target','X_2','--solver',sys.executable,
                             '--output',str(self.root)])==3

    def test_defining_relation_nonzero_against_truncated_basis_is_inconclusive(self):
        relation=self.pres.target(2,1,NumericRing(5,7,0))
        result=self.discover(config=self.configuration('omit'))
        assert result['status']=='inconclusive'
        assert result['in_ideal'] is None
        assert result['remainder_terms']==len(relation.terms)
        assert not result['basis_complete']

    @pytest.mark.parametrize('mode',['schema','matrix','kernel','certificate','representative',
                                      'support','reconstruction','solver-links'])
    def test_serialization_corruption_cannot_publish_checked_metadata_for_different_objects(self,mode):
        from pathlib import Path
        from zipper_algebra.cofactor.cofactor_io import Cofactor
        prior=self.discover()
        previous=Path(prior['certificate']).read_bytes()
        writer=pipeline.write_json_atomic
        saver=Cofactor.save
        def altered(path,data):
            name=Path(path).name
            if name=={'support':'support.json','reconstruction':'reconstruction.json',
                       'solver-links':'solver-runs.json'}.get(mode):
                data={}
            if name=='matrix.json' and mode in ('schema','matrix'):
                data={**data,'schema':'invalid'} if mode=='schema' else {**data,'entries':[]}
            if name=='matrix.json' and mode=='representative':
                from zipper_algebra.fields import FieldElement
                from zipper_algebra.matrices import ExactMatrix
                matrix=ExactMatrix.from_json(data)
                first=matrix.entries[0][2]
                factor=first.field.parameter('q')-11
                top=first.field.polynomial(dict(first.numerator))*factor
                bottom=first.field.polynomial(dict(first.denominator))*factor
                changed=FieldElement(first.field,dict(top.numerator),dict(bottom.numerator))
                assert changed==first  # Same rational function, different locus.
                encoded=matrix.to_json()
                encoded['entries'][0][2]=changed.to_json()
                assert ExactMatrix.from_json(encoded).content_hash()!=matrix.content_hash()
                data=encoded
            if name=='kernel.json' and mode=='kernel':
                data={**data,'matrix_hash':'0'*64}
            return writer(path,data)
        def alter_cofactor(cofactor,path,**options):
            saver(cofactor,path,**options)
            if mode=='certificate':
                Path(path).write_text(Path(path).read_text().replace('1*q^0*t^0','2*q^0*t^0'))
        with patch.object(pipeline,'write_json_atomic',side_effect=altered), \
             patch.object(Cofactor,'save',alter_cofactor):
            result=self.discover()
        assert result['status']=='inconclusive', result
        assert result['execution_status']=='failed'
        assert result['stage']=='publication'
        assert 'invalid serialized reconstruction' in result['reason']
        assert 'certificate' not in result
        assert Path(prior['certificate']).read_bytes()==previous

    def test_saved_reconstruction_options_match_explicit_controls(self):
        import json
        from pathlib import Path
        report=self.discover(threads=2,max_primes=3,reconstruction_degree=19,
                             start_prime=1,shift=[3,11],max_contexts=123,max_products=777)
        assert report['verified']
        assert report['reconstruction']['options']=={
            'threads':2,'max_primes':3,'max_degree':19,'start_prime':1,
            'shift':[3,11],'max_contexts':123,'max_products':777}
        assert json.loads(Path(report['artifacts']['reconstruction.json']).read_text())==report['reconstruction']

    @pytest.mark.parametrize('kind',['skeleton','cofactor-header','cofactor-terms'])
    def test_bound_mathematical_artifact_must_match_checked_discovery(self,kind):
        from pathlib import Path
        from copy import deepcopy
        prior=self.discover()
        before=Path(prior['certificate']).read_bytes()
        save=pipeline.save_artifact
        def changed(path,value,*arguments,**options):
            name=Path(path).name
            if name=='skeleton.json' and kind=='skeleton':
                value=deepcopy(value)
                value.anchor=(11,13,0)
            if name=='certificate.json' and kind.startswith('cofactor'):
                value=deepcopy(value)
                if kind=='cofactor-header':
                    value.producer='other technical producer'
                else:
                    value.terms=[]
            return save(path,value,*arguments,**options)
        with patch.object(pipeline,'save_artifact',side_effect=changed):
            report=self.discover()
        assert report['status']=='inconclusive'
        assert 'invalid serialized reconstruction' in report['reason']
        assert Path(prior['certificate']).read_bytes()==before

    @pytest.mark.parametrize('kind',['extra','ordinary-file'])
    def test_preexisting_generation_must_have_exact_file_set(self,kind):
        from pathlib import Path
        import shutil
        prior=self.discover()
        directory=Path(prior['certificate']).parent
        if kind=='extra':
            (directory/'unexpected.txt').write_text('not part of this generation')
        else:
            shutil.rmtree(directory)
            directory.write_text('a directory was required')
        report=self.discover()
        assert report['status']=='inconclusive'
        assert report['execution_status']=='failed'
        assert 'existing verified generation has different content' in report['reason']

    def test_existing_generation_entry_changed_into_directory_is_failed_publication(self):
        from pathlib import Path
        prior=self.discover()
        certificate=Path(prior['certificate'])
        certificate.unlink()
        certificate.mkdir()
        result=self.discover()
        assert result['status']=='inconclusive'
        assert result['execution_status']=='failed'
        assert result['stage']=='publication'
        assert 'could not publish verified generation' in result['reason']
        assert 'certificate' not in result
        assert certificate.is_dir()

    def test_denied_generation_rename_is_a_failed_publication(self):
        from pathlib import Path
        replace=pipeline.os.replace
        def denied(source,destination):
            if Path(source).name=='generation' and Path(destination).parent.name=='verified':
                raise PermissionError('publication rename denied')
            return replace(source,destination)
        with patch.object(pipeline.os,'replace',side_effect=denied):
            result=self.discover()
        assert result['status']=='inconclusive'
        assert result['execution_status']=='failed'
        assert result['stage']=='publication'
        assert 'publication rename denied' in result['reason']
        assert 'certificate' not in result
