"""Portable discovery using exact, explicitly synthetic process outputs."""
from __future__ import annotations

from contextlib import closing, redirect_stdout, redirect_stderr
from dataclasses import replace
from fractions import Fraction
from itertools import islice
import io
import json
from pathlib import Path
import sys
import tempfile
import pytest
from unittest.mock import patch

from zipper_algebra import Presentation, Atom, AtomSet, Cofactor, cli
from zipper_algebra.coefficients import NumericRing
from zipper_algebra.cofactor.artifacts import load_artifact
from zipper_algebra.cofactor.compose import main as compose_command
from zipper_algebra.execution.config import RunConfig
from zipper_algebra.execution import pipeline
from zipper_algebra.execution.sampling import SampleSet, admissible, spiral


# This is a process-control fixture, not a Gröbner-basis implementation. It
# returns one explicitly imposed input relation, optionally multiplied by q.
PRODUCER = r'''
import argparse
from fractions import Fraction
from pathlib import Path
import re
import sys
import time
p = argparse.ArgumentParser()
p.add_argument('-i');p.add_argument('-o');p.add_argument('-p')
p.add_argument('--mode',default='plain')
a,_=p.parse_known_args()
if a.mode=='exit':sys.exit(9)
if a.mode=='logged':print(time.time_ns())
relations=Path(a.i).read_text().splitlines()[2:]
if a.mode=='reanchor-remainder':
    # Fixed three-strand fixture; the parent checks this input-label position.
    t2=relations[11]
    q=int(re.search(r'([+-]?\d+)\*x1\*y1',t2).group(1))
    terms={}
    for factor,source in ((1,relations[-1].rstrip(',')),(q-5,relations[0].rstrip(','))):
        for piece in re.findall(r'[+-]?[^+-]+',source):
            sign,coefficient,word=re.match(r'([+-]?)(?:(\d+(?:/\d+)?)\*?)?(.*)',piece).groups()
            coefficient=Fraction(coefficient or 1)*(-1 if sign=='-' else 1)*factor
            terms[word]=terms.get(word,0)+coefficient
    pieces=[]
    for word,coefficient in terms.items():
        if coefficient:
            text=str(abs(coefficient))+('*'+word if word else '')
            pieces.append(('-' if coefficient<0 else '+' if pieces else '')+text)
    Path(a.o).write_text(re.sub(r'(?<!^)([+-])',r' \1 ',''.join(pieces))+'\n', newline='\n')
    proof='[i%d]'%(len(relations)-1)
    if q!=5:proof+=' '+('+' if q>5 else '-')+str(abs(q-5))+'*[i0]'
    Path(a.p).write_text(proof+'\n', newline='\n')
    sys.exit(0)
q=int(re.search(r'([+-]?\d+)\*y1',relations[-1]).group(1))
index=0 if a.mode=='omit' or (a.mode=='vary' and q!=5) else len(relations)-1
source=relations[index].rstrip(',')
scale=q if a.mode=='scaled' else Fraction(1,q**5) if a.mode=='pole' else 1
parts=[]
for piece in re.findall(r'[+-]?[^+-]+',source):
    match=re.match(r'([+-]?)(?:(\d+(?:/\d+)?)\*?)?(.*)',piece)
    sign,coefficient,word=match.groups()
    coefficient=Fraction(coefficient or 1)*(-1 if sign=='-' else 1)*scale
    text=(str(abs(coefficient))+('*'+word if word else ''))
    parts.append(('-' if coefficient<0 else '+' if parts else '')+text)
basis=re.sub(r'(?<!^)([+-])',r' \1 ',''.join(parts))
if a.mode=='closure':
    h=re.sub(r'(?<!^)([+-])',r' \1 ',relations[0].rstrip(','))
    hproof='[i0] +[i1] -[i1]' if q==3 else '[i0]'
    gproof='[i3]' if q==5 else '[i3] +[0] -[i0]'
    Path(a.o).write_text(h+'\n'+basis+'\n', newline='\n')
    Path(a.p).write_text(hproof+'\n'+gproof+'\n', newline='\n')
    sys.exit(0)
Path(a.o).write_text(basis+'\n', newline='\n')
proof=(str(scale)+'*' if scale!=1 else '')+'[i%d]'%index
if a.mode=='bad-proof':proof='[i0]'
if a.mode=='shape' and q!=5:proof+=' +[i0] -[i0]'
Path(a.p).write_text(proof+'\n', newline='\n')
'''


class DiscoveryFixture:
    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        self.root = tmp_path
        script = self.root / 'producer.py'
        script.write_text(PRODUCER)
        self.pres = Presentation(2, zero_from=2)
        self.target = 'X_2'
        self.config = RunConfig(Path(sys.executable), self.root / 'result',
                                self.root / 'cache', command_prefix=(script,))
        self.atoms = AtomSet([])

    def configuration(self, mode):
        return replace(self.config, flags=('--mode', mode))

    def discover(self, **kwargs):
        return pipeline.run(self.pres, self.target, (5, 7, 0),
                            kwargs.pop('config', self.config), start=kwargs.pop('start', 7),
                            step=kwargs.pop('step', 2), maxpts=kwargs.pop('maxpts', 9),
                            atoms=kwargs.pop('atoms', self.atoms), verbose=False, **kwargs)


class TestPortableSamples(DiscoveryFixture):
    def test_spiral_preserves_the_historical_prime_pair_order(self):
        assert (list(islice(spiral(), 10))) == ([(3, 2), (2, 3), (5, 2), (5, 3), (3, 5),
                          (2, 5), (7, 2), (7, 3), (7, 5), (5, 7)])
        atoms = AtomSet([Atom('t-q', {(0, 1): 1, (1, 0): -1})])
        assert not (admissible(2, 3, atoms))
        assert not (admissible(5, 5, atoms))
        assert (admissible(5, 7, atoms))

    def test_selection_is_pure_and_parallel_growth_is_deterministically_ordered(self):
        samples = SampleSet(self.pres, self.target, (5, 7, 0), self.config, self.atoms)
        assert (samples.next_points(0)) == ([])
        assert (samples.grow(0)) == ([])
        assert not (self.config.output_root.exists())
        with pytest.raises(ValueError, match='nonnegative'):
            samples.next_points(-1)
        with pytest.raises(ValueError, match='positive'):
            samples.grow(0, jobs=0)
        selected = samples.next_points(3)
        assert (samples.grow(3, jobs=2)) == (selected)
        assert (samples.points) == (selected)
        samples.run_anchor()
        assert (samples.run_anchor().hit)
        assert (samples.char) == (0)
        assert (samples.anchor_basis()) == (samples.basis(5, 7))
        assert (samples.anchor_proof()) is (samples.proof(5, 7))
        samples.reanchor(selected[0])
        assert (samples.anchor) == ((*selected[0], 0))
        with pytest.raises(KeyError):
            samples.reanchor((101, 103))

    def test_census_names_actual_proof_shapes(self):
        samples = SampleSet(self.pres, self.target, (5, 7, 0), self.config, self.atoms)
        samples.run_anchor()
        samples.grow(2)
        skeleton = pipeline._skeleton(self.pres, self.target, samples)
        modal, off = samples.modal_shapes(skeleton)
        assert (off) == ([])
        assert (set(modal)) == (set(skeleton.basis_keys))
        assert (samples.stable(skeleton)) == (samples.points)
        varied = SampleSet(self.pres, self.target, (5, 7, 0), self.configuration('vary'), self.atoms)
        varied.run_anchor()
        varied.grow(7)
        skeleton = pipeline._skeleton(self.pres, self.target, varied)
        census = varied.shape_census(skeleton)
        assert (census['per_point'][(3,2)][skeleton.basis_keys[0]]) is None
        # Most truncated outputs omit the required key. None cannot become a
        # proof shape; the anchor's actual shape remains the comparison reference.
        assert (varied.modal_shapes(skeleton)[1]) == ([])
        assert (varied.stable(skeleton)) == ([(5,7),(5,2),(5,3)])



class TestFreshDiscovery(DiscoveryFixture):
    def test_process_to_bound_artifacts_and_exact_recomposition(self):
        report = self.discover()
        assert (report['status']) == ('verified')
        assert (report['verified'])
        assert (report['in_ideal'])
        assert not (report['basis_complete'])
        output = Path(report['certificate']).parent
        assert (json.loads((self.config.output_root/'result.json').read_text())['status']) == ('verified')
        assert (len(report['runs'])) == (7)
        for run in report['runs']:
            assert (Path(run['store']).is_file())
        assert (len({run['store'] for run in report['runs']})) == (1)
        order = {'solver': 'degree-lex', 'reduction': pipeline.REDUCTION_ORDER}
        skeleton = load_artifact(output/'skeleton.json', self.pres, self.target,
                                 kind='skeleton', order=order)
        fitted = load_artifact(output/'lift.json', self.pres, self.target, kind='lift', order=order)
        assert (len(skeleton.steps)) == (1)
        assert (fitted.samples_used) == (7)
        certificate = Cofactor.load(output/'certificate.txt', pres=self.pres)
        assert (certificate.terms[0][1]) == ('ZERO_2_1')
        recomposed = self.root/'recomposed.txt'
        with redirect_stdout(io.StringIO()):
            code = compose_command(['--n','2','--zero-from','2','--target','X_2',
                                    '--lift',str(output/'lift.json'), '--skeleton',str(output/'skeleton.json'),
                                    '--multipliers',str(output/'multipliers.json'),'--out',str(recomposed)])
        assert (code) == (0)
        assert (Cofactor.load(recomposed, pres=self.pres).terms) == (certificate.terms)
        assert (self.discover()['status']) == ('verified')  # Validated cache reuse.

    def test_symbolic_multiplier_is_reconstructed_instead_of_using_the_anchor(self):
        atoms = AtomSet([Atom('q', {(1,0): 1})], [1])
        report = self.discover(config=self.configuration('scaled'), atoms=atoms,
                               start=9, step=10, maxpts=49)
        assert (report['verified']), (report)
        multipliers = json.loads(Path(report['artifacts']['multipliers.json']).read_text())
        # X2 = (1/q)*(q X2). Its coefficient is not the anchor constant 1/5.
        assert (multipliers['steps'][0]['num']) == ([[[0,0], '1']])
        assert (multipliers['steps'][0]['denom']) == ([1])

    def test_high_order_coefficient_pole_uses_lines_and_composes_exactly(self):
        atoms = AtomSet([Atom('q', {(1,0): 1})], [6])
        report = self.discover(config=self.configuration('pole'), atoms=atoms,
                               start=100, maxpts=100)
        assert (report['verified']), (report)
        order = {'solver': 'degree-lex', 'reduction': pipeline.REDUCTION_ORDER}
        fitted = load_artifact(report['artifacts']['lift.json'], self.pres, self.target,
                               kind='lift', order=order)
        coefficient, = next(iter(fitted.entries.values()))
        assert (coefficient[0]) == ((5,))
        assert (coefficient[1]) == ({(0,0): Fraction(1)})
        assert (fitted.timing['routes']) == ({'known': 1, 'enum': 0})
        # The multiplier q^5 cancels the reconstructed basis coefficient 1/q^5.
        multipliers = json.loads(Path(report['artifacts']['multipliers.json']).read_text())
        numerator = {tuple(m): Fraction(c) for m, c in multipliers['steps'][0]['num']
                     if Fraction(c)}
        assert (numerator) == ({(5,0): Fraction(1)})
        assert (multipliers['steps'][0]['denom']) == ([0])

    def test_truncated_output_cannot_establish_nonmembership(self):
        report = self.discover(config=self.configuration('omit'))
        assert (report['status']) == ('inconclusive')
        assert (report['in_ideal']) is None
        assert (report['remainder_terms']) > (0)
        assert not ((self.config.output_root/'certificate.txt').exists())

    def test_failed_process_is_an_inconclusive_computation(self):
        report = self.discover(config=self.configuration('exit'))
        assert (report['execution_status']) == ('failed')
        assert not (report['verified'])
        assert ('code 9') in (report['reason'])

    def test_sampling_caps_grow_and_remain_inconclusive_without_enough_data(self):
        report = self.discover(start=3, step=2, maxpts=5)
        assert (report['points']) == (5)
        assert (report['status']) == ('inconclusive')
        assert ('held-out') in (report['reason'])
        assert not ((self.config.output_root/'certificate.json').exists())

    def test_degenerate_anchor_is_replaced_using_valid_equivalent_proof_shapes(self):
        report = self.discover(config=self.configuration('shape'), start=7, maxpts=11)
        assert (report['verified']), (report)
        assert (report['anchor']) == ([3,2,0])
        assert (report['points']) > (7)

    def test_false_symbolic_fit_at_the_cap_remains_inconclusive(self):
        atoms = AtomSet([Atom('q', {(1,0):1})],[1])
        report = self.discover(config=self.configuration('scaled'), atoms=atoms,
                               start=9, maxpts=11)
        assert not (report['verified'])
        assert ('failed exact verification') in (report['reason'])
        assert not ((self.config.output_root/'certificate.txt').exists())

    def test_replacement_anchor_is_checked_against_its_new_proof_closure(self):
        result = self.discover(config=self.configuration('closure'), start=11, maxpts=11)
        assert (result['status']) == ('inconclusive')
        assert ('new closure') in (result['reason'])
        assert (result['anchor']) == ([3,2,0])
        assert not ((self.config.output_root/'certificate.txt').exists())

    def test_changed_closure_can_grow_and_reanchor_without_dropping_slots(self):
        result = self.discover(config=self.configuration('closure'), start=11, step=5, maxpts=26)
        assert (result['verified']), (result)
        assert (result['anchor'][:2]) != ([3,2])
        assert (result['points']) > (11)

    def test_interrupted_final_record_preserves_prior_verified_generation(self):
        from zipper_algebra.execution.artifacts import file_identity
        prior = self.discover()
        marker = self.config.output_root/'result.json'
        prior_record = marker.read_bytes()
        prior_files = {name: Path(path).read_bytes() for name,path in prior['artifacts'].items()}
        original_write = pipeline.write_json_atomic
        def interrupted(path, data):
            if Path(path).name == 'result.json':
                raise OSError('interrupted final record promotion')
            return original_write(path,data)
        other = Presentation(2, 'square-zero', zero_from=2)
        with patch.object(pipeline,'write_json_atomic',side_effect=interrupted), \
             pytest.raises(OSError, match='interrupted'):
            pipeline.run(other,self.target,(5,7,0),self.config,start=7,step=2,maxpts=9,
                         atoms=self.atoms,verbose=False)
        assert (marker.read_bytes()) == (prior_record)
        for name,path in prior['artifacts'].items():
            assert (Path(path).read_bytes()) == (prior_files[name])
            assert (file_identity(Path(path))) == (prior['files'][name])

    def test_new_context_inconclusive_report_does_not_select_older_certificate(self):
        prior = self.discover()
        old_bytes = Path(prior['certificate']).read_bytes()
        other = Presentation(2,'square-zero',zero_from=2)
        result = pipeline.run(other,self.target,(5,7,0),self.configuration('omit'),
                              start=7,maxpts=9,atoms=self.atoms,verbose=False)
        assert (result['status']) == ('inconclusive')
        assert ('certificate') not in (result)
        assert (result['presentation']['sha256']) != (prior['presentation']['sha256'])
        assert (Path(prior['certificate']).read_bytes()) == (old_bytes)

    def test_old_verified_result_binds_immutable_solver_generations_after_refresh(self):
        from zipper_algebra.execution.store import read
        from zipper_algebra.execution.artifacts import digest_bytes
        configuration = self.configuration('logged')
        prior = self.discover(config=configuration)
        refreshed = self.discover(config=configuration, refresh=True)
        assert (prior['runs'][0]['key']) == (refreshed['runs'][0]['key'])
        assert (prior['runs'][0]['generation']) != (refreshed['runs'][0]['generation'])
        for run in prior['runs']:
            _, selected = read(run['store'], run['key'])
            assert (selected['generation']) != (run['generation'])
            files, _ = read(run['store'], run['key'], generation=run['generation'])
            for name, identity in run['files'].items():
                assert ({'bytes':len(files[name]), 'sha256':digest_bytes(files[name])}) == (identity)
            assert (run['request']) == (selected['request'])

    def test_existing_verified_generation_is_checked_before_reuse(self):
        from zipper_algebra.execution.artifacts import ArtifactError
        from zipper_algebra.cofactor.compose import load_multipliers
        report = self.discover()
        paths = {name:Path(path) for name,path in report['artifacts'].items()}
        order={'solver':'degree-lex','reduction':pipeline.REDUCTION_ORDER}
        skeleton=load_artifact(paths['skeleton.json'],self.pres,self.target,kind='skeleton',order=order)
        fitted=load_artifact(paths['lift.json'],self.pres,self.target,kind='lift',order=order)
        certificate=load_artifact(paths['certificate.json'],self.pres,self.target,kind='cofactor')
        multipliers=load_multipliers(paths['multipliers.json'],skeleton,fitted.atoms,self.pres)
        repeated=pipeline._publish_verified(self.config.output_root,self.pres,self.target,skeleton,
                                           fitted,multipliers,certificate,order)
        assert (repeated['generation']) == (report['generation'])
        paths['certificate.txt'].write_text('corrupt')
        with pytest.raises(ArtifactError, match='different content'):
            pipeline._publish_verified(self.config.output_root,self.pres,self.target,skeleton,
                                       fitted,multipliers,certificate,order)

    def test_nonzero_remainder_at_new_anchor_is_inconclusive(self):
        # G(q)=X3+(q-5)*INV_1_p belongs to the defining ideal and equals X3
        # at the original anchor. At q=3, reducing X3 by this truncated output
        # leaves 2*INV_1_p; its absence from the supplied basis proves nothing.
        pres=Presentation(3,zero_from=3)
        assert (pres.labels()[11]) == ('T_2')
        report=pipeline.run(pres,'X_3',(5,7,0),self.configuration('reanchor-remainder'),
                            start=7,maxpts=9,atoms=self.atoms,verbose=False)
        assert (report['status']) == ('inconclusive')
        assert ('replacement anchor') in (report['reason'])
        assert (report['remainder_terms']) > (0)
        assert (report['in_ideal']) is None

    def test_invalid_limits_and_anchor_fail_before_creating_outputs(self, subtests):
        for arguments in ({'start':0}, {'step':0}, {'jobs':0}, {'maxpts':0}, {'start':10,'maxpts':9}):
            with subtests.test(arguments=arguments), pytest.raises(ValueError):
                self.discover(**arguments)
        for anchor in ((1,2,5), (True,2,0), (1,2), [1,2,0]):
            with subtests.test(anchor=anchor), pytest.raises(ValueError):
                pipeline.run(self.pres,self.target,anchor,self.config)
        assert not (self.config.output_root.exists())


class TestInstalledDiscoveryCommands(DiscoveryFixture):
    def invoke(self, arguments):
        output, errors = io.StringIO(), io.StringIO()
        with redirect_stdout(output), redirect_stderr(errors):
            code = cli.main(arguments)
        return code, output.getvalue(), errors.getvalue()

    def invoke_with_producer(self, arguments):
        # Supply only the interpreter prefix missing from the installed CLI;
        # discovery, cache/proof validation and subprocess execution all run.
        def configuration(*args, **kwargs):
            return RunConfig(*args, command_prefix=(self.root/'producer.py',), **kwargs)
        with patch('zipper_algebra.execution.config.RunConfig', side_effect=configuration):
            return self.invoke(arguments)

    def test_emit_json_and_text_preserve_bare_presentation(self):
        path = self.root/'inputs'/'bare.ms'
        code, output, errors = self.invoke(['emit','--n','2','--out',str(path),'--json'])
        assert (code) == (0), (errors)
        assert (json.loads(output)['relations']) == (3)
        assert (len(path.read_text().splitlines())-2) == (3)
        assert (self.invoke(['emit','--n','2','--out',str(path)])[0]) == (0)
        assert ('Wrote') in (self.invoke(['emit','--n','2','--out',str(path)])[1])

    def test_failed_rerun_replaces_verified_selection_and_preserves_immutable_products(self, subtests):
        script = self.root/'producer.py'
        atoms_path = self.root/'atoms.json'
        self.atoms.save(atoms_path)
        for failure in ('output-record', 'cache-record', 'bad-proof', 'exit'):
            with subtests.test(failure=failure):
                script.write_text(PRODUCER)
                output, cache = self.root/failure/'output', self.root/failure/'cache'
                for directory in (output, cache):
                    directory.mkdir(parents=True)
                    (directory/'.gitignore').write_text('existing output policy\n')
                arguments = ['discover', '--n', '2', '--zero-from', '2', '--target', 'X_2',
                             '--solver', sys.executable, '--output', str(output), '--cache', str(cache),
                             '--atoms', str(atoms_path), '--start', '7', '--step', '2',
                             '--max', '9', '--json']
                code, text, errors = self.invoke_with_producer(arguments)
                assert (code) == (0), (errors or text)
                prior = json.loads(text)
                assert (prior['verified'])
                prior_files = {Path(path): Path(path).read_bytes()
                               for path in prior['artifacts'].values()}
                import sqlite3
                databases = (output/'samples.sqlite3', cache/'samples.sqlite3')
                def raw_samples(database):
                    with closing(sqlite3.connect(database)) as connection, connection:
                        return connection.execute(
                            'SELECT key,generation,name,content FROM artifacts ORDER BY key,generation,name').fetchall()
                before_samples = {database:raw_samples(database) for database in databases}
                first = prior['runs'][0]
                if failure == 'output-record':
                    with closing(sqlite3.connect(databases[0])) as connection, connection:
                        connection.execute('UPDATE generations SET record=? WHERE key=? AND generation=?',
                                           (b'{}', first['key'], first['generation']))
                elif failure == 'cache-record':
                    with closing(sqlite3.connect(databases[0])) as connection, connection:
                        connection.execute('DELETE FROM selections WHERE key=?', (first['key'],))
                    with closing(sqlite3.connect(databases[1])) as connection, connection:
                        connection.execute('UPDATE generations SET record=? WHERE key=? AND generation=?',
                                           (b'{}', first['key'], first['generation']))
                else:
                    script.write_text(PRODUCER.replace(
                        "default='plain'", "default='"+failure+"'"))
                    arguments.append('--refresh')
                code, text, errors = self.invoke_with_producer(arguments)
                assert (code) == (3), (errors or text)
                result = json.loads(text)
                assert (result) == (json.loads((output/'result.json').read_text()))
                assert (result['status']) == ('inconclusive')
                assert (result['execution_status']) == ('failed')
                assert not (result['verified'])
                assert (result['in_ideal']) is None
                assert not (result['basis_complete'])
                for field in ('certificate', 'artifacts', 'generation', 'files', 'exact'):
                    assert (field) not in (result)
                expected = ('invalid stored sample' if failure.endswith('record') else
                            'input-ideal identity' if failure == 'bad-proof' else 'code 9')
                assert (expected) in (result['reason'])
                for path, contents in prior_files.items():
                    assert (path.read_bytes()) == (contents)
                for database, contents in before_samples.items():
                    assert (raw_samples(database)) == (contents)
                for directory in (output, cache):
                    assert ((directory/'.gitignore').read_text()) == ('existing output policy\n')

    def test_invalid_discovery_configuration_preserves_prior_selection(self):
        prior = self.discover()
        marker = self.config.output_root/'result.json'
        contents = marker.read_bytes()
        atoms_path = self.root/'atoms.json'
        self.atoms.save(atoms_path)
        arguments = ['discover', '--n', '2', '--zero-from', '2', '--target', 'X_2',
                     '--solver', sys.executable,
                     '--output', str(self.config.output_root), '--atoms', str(atoms_path),
                     '--start', '10', '--max', '9', '--json']
        code, text, errors = self.invoke_with_producer(arguments)
        assert (code) == (2), (errors)
        assert (json.loads(text)['status']) == ('invalid-input')
        assert (marker.read_bytes()) == (contents)
        assert (Path(prior['certificate']).is_file())

    def test_discover_dispatches_explicit_configuration_and_exit_statuses(self, subtests):
        args=['discover','--n','2','--zero-from','2','--target','X_2',
              '--solver',sys.executable,'--output',str(self.root/'discovery'),
              '--cache',str(self.root/'shared-cache'),'--anchor','3,7',
              '--start','7','--step','2','--max','9','--jobs','2','--refresh','--json']
        atoms_path=self.root/'atoms.json'
        self.atoms.save(atoms_path)
        for verified in (True,False):
            with subtests.test(verified=verified), patch.object(pipeline,'run',return_value={
                    'verified':verified,'status':'verified' if verified else 'inconclusive'}) as producer:
                code, output, errors=self.invoke(args+['--atoms',str(atoms_path)])
                assert (code) == (0 if verified else 3), (errors)
                assert (json.loads(output)['verified']) == (verified)
                assert (producer.call_args.args[2]) == ((3,7,0))
                assert (producer.call_args.kwargs['jobs']) == (2)
        with patch.object(pipeline,'run',return_value={'verified':True,'status':'verified'}):
            assert (self.invoke(args[:-1])[1].strip()) == ('verified')
        for value in ('bad','1,2,3','x,2'):
            with subtests.test(value=value),redirect_stderr(io.StringIO()),pytest.raises(SystemExit):
                cli.main(args+['--anchor',value])

    def test_installed_compose_has_the_shared_exact_verification_entrypoint(self):
        report=self.discover()
        output=Path(report['certificate']).parent
        code, text, errors=self.invoke(['compose','--n','2','--zero-from','2',
             '--target','X_2','--skeleton',str(output/'skeleton.json'),
             '--lift',str(output/'lift.json'),'--multipliers',str(output/'multipliers.json'),
             '--out',str(self.root/'installed.txt'),'--json'])
        assert (code) == (0), (errors)
        assert (json.loads(text)['pass'])
        assert (report['verified'])

class TestCompositionArtifactBoundaries(DiscoveryFixture):
    def test_default_sampling_factors_are_local_to_each_run(self):
        left = SampleSet(self.pres,self.target,(5,7,0),self.config)
        right = SampleSet(self.pres,self.target,(5,7,0),self.config)
        assert (left.atoms) is not (right.atoms)
        old = right.atoms.to_json()
        left.atoms.max_exp = tuple(0 for _ in left.atoms.max_exp)
        assert (right.atoms.to_json()) == (old)

    def test_bound_composition_requires_explicit_target_and_rejects_bad_documents(self):
        report=self.discover()
        output=Path(report['certificate']).parent
        destination=self.root/'protected.txt'
        destination.write_text('prior valid result\n')
        arguments=['--n','2','--zero-from','2','--lift',str(output/'lift.json'),
                   '--skeleton',str(output/'skeleton.json'),'--multipliers',str(output/'multipliers.json'),
                   '--out',str(destination),'--json']
        with redirect_stdout(io.StringIO()):
            assert (compose_command(arguments)) == (2)
        assert (destination.read_text()) == ('prior valid result\n')
        (output/'skeleton.json').write_text('[]')
        with redirect_stdout(io.StringIO()):
            assert (compose_command(arguments+['--target','X_2'])) == (2)
        assert (destination.read_text()) == ('prior valid result\n')

    def test_legacy_inputs_still_reject_an_explicit_different_target(self):
        from zipper_algebra.cofactor.skeleton import Skeleton
        from zipper_algebra.cofactor.lift import Lift
        pres=Presentation(3)
        skeleton=Skeleton(pres.slug(),pres.content_hash(),'X_2',(5,7,0),remainder_terms=3)
        fitted=Lift(self.atoms,{},0,0,(5,7,0),{})
        skeleton.save(self.root/'skeleton.json')
        fitted.save(self.root/'lift.json')
        with redirect_stdout(io.StringIO()):
            code=compose_command(['--n','3','--target','X_3','--skeleton',str(self.root/'skeleton.json'),
                 '--lift',str(self.root/'lift.json'),'--multipliers',str(self.root/'unused.json'),
                 '--out',str(self.root/'bad.txt'),'--json'])
        assert (code) == (2)
        assert not ((self.root/'bad.txt').exists())


class TestConcurrentVerifiedPublication(DiscoveryFixture):
    def material(self):
        from zipper_algebra.cofactor.compose import load_multipliers
        prior = self.discover()
        paths = {name:Path(path) for name,path in prior['artifacts'].items()}
        order = {'solver':'degree-lex','reduction':pipeline.REDUCTION_ORDER}
        skeleton = load_artifact(paths['skeleton.json'],self.pres,self.target,kind='skeleton',order=order)
        fitted = load_artifact(paths['lift.json'],self.pres,self.target,kind='lift',order=order)
        certificate = load_artifact(paths['certificate.json'],self.pres,self.target,kind='cofactor')
        multipliers = load_multipliers(paths['multipliers.json'],skeleton,fitted.atoms,self.pres)
        return prior, (self.pres,self.target,skeleton,fitted,multipliers,certificate,order)

    def test_identical_generation_publishers_reuse_the_complete_winner(self):
        from concurrent.futures import ThreadPoolExecutor
        import os
        import threading
        from zipper_algebra.execution.artifacts import file_identity
        prior, arguments = self.material()
        prior_record = (self.config.output_root/'result.json').read_bytes()
        barrier = threading.Barrier(2)
        replace = os.replace
        rename_failures = []
        def concurrent_replace(source,destination):
            if Path(source).name == 'generation' and Path(destination).parent.name == 'verified':
                barrier.wait(timeout=10)
                try:
                    return replace(source,destination)
                except OSError as error:
                    rename_failures.append(error)
                    raise
            return replace(source,destination)
        with patch.object(pipeline.os,'replace',side_effect=concurrent_replace):
            with ThreadPoolExecutor(max_workers=2) as executor:
                first = executor.submit(pipeline._publish_verified,self.root/'concurrent',*arguments)
                second = executor.submit(pipeline._publish_verified,self.root/'concurrent',*arguments)
                results = [first.result(timeout=15),second.result(timeout=15)]
        assert (len(rename_failures)) == (1)
        assert (results[0]['generation']) == (results[1]['generation'])
        for name,path in results[0]['artifacts'].items():
            assert (file_identity(Path(path))) == (results[0]['files'][name])
        assert ((self.config.output_root/'result.json').read_bytes()) == (prior_record)
        assert (Path(prior['certificate']).is_file())

    def test_unrelated_rename_failure_is_propagated_without_a_winner(self):
        import os
        prior, arguments = self.material()
        marker = self.config.output_root/'result.json'
        prior_record = marker.read_bytes()
        replace = os.replace
        error = PermissionError('generation rename denied')
        def denied(source,destination):
            if Path(source).name == 'generation' and Path(destination).parent.name == 'verified':
                raise error
            return replace(source,destination)
        with patch.object(pipeline.os,'replace',side_effect=denied), pytest.raises(PermissionError) as raised:
            pipeline._publish_verified(self.root/'denied',*arguments)
        assert (raised.value) is (error)
        assert (marker.read_bytes()) == (prior_record)
        assert (Path(prior['certificate']).is_file())

    def test_corrupt_or_incomplete_concurrent_winner_cannot_hide_a_rename_error(self, subtests):
        import os
        prior, arguments = self.material()
        marker = self.config.output_root/'result.json'
        prior_record = marker.read_bytes()
        replace = os.replace
        for mode in ('corrupt','incomplete'):
            error = OSError('raced generation rename')
            def raced(source,destination):
                if Path(source).name == 'generation' and Path(destination).parent.name == 'verified':
                    replace(source,destination)
                    product = Path(destination)/'certificate.txt'
                    if mode == 'corrupt':
                        product.write_text('corrupted concurrent product')
                    else:
                        product.unlink()
                    raise error
                return replace(source,destination)
            with subtests.test(mode=mode),patch.object(pipeline.os,'replace',side_effect=raced), \
                 pytest.raises(OSError) as raised:
                pipeline._publish_verified(self.root/mode,*arguments)
            assert (raised.value) is (error)
            assert (marker.read_bytes()) == (prior_record)
            assert (Path(prior['certificate']).is_file())
