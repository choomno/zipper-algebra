"""Entry points distinguish augmented presentations and verified symbolic identities."""
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import pytest
from unittest.mock import patch

from zipper_algebra import Presentation
from zipper_algebra import msolve as emitter
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor import compose as composition
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.cofactor.skeleton import Skeleton


class TestExplicitAugmentedPresentation:
    @staticmethod
    def invoke(arguments):
        output, errors = io.StringIO(), io.StringIO()
        with patch.object(sys, "argv", ["msolve", *arguments]):
            with redirect_stdout(output), redirect_stderr(errors):
                code = emitter.main()
        return code, output.getvalue(), errors.getvalue()

    def test_default_emission_keeps_only_defining_relations(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "bare.ms"
            assert self.invoke(['--n', '2', '--out', str(output)])[0] == 0
            lines = output.read_text().splitlines()
            assert lines[:2] == ['x1,y1', '0']
            assert len(lines) - 2 == 3
            assert lines[-1] == '-x1*x1+5*x1*y1+3*x1-35*y1+28'

    def test_append_target_is_an_explicit_extra_defining_relation(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "augmented.ms"
            code, text, errors = self.invoke(
                ["--n", "2", "--append-target", "2", "--out", str(output)])
            assert code == 0
            assert errors == ''
            lines = output.read_text().splitlines()
            assert len(lines) - 2 == 4
            assert lines[-1] == '-x1+5*y1-4'
            assert 'Augmented presentation' in text
            assert 'look at line' not in text

    def test_removed_target_flag_rejects_before_overwriting_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "existing.ms"
            output.write_text("keep the existing input\n")
            errors = io.StringIO()
            with patch.object(sys, "argv", ["msolve", "--n", "3", "--target", "3",
                                            "--out", str(output)]):
                with redirect_stderr(errors), pytest.raises(SystemExit) as raised:
                    emitter.main()
            assert raised.value.code == 2
            assert '--append-target' in errors.getvalue()
            assert 'does not prove membership' in errors.getvalue()
            assert output.read_text() == 'keep the existing input\n'

    def test_module_help_describes_actual_relations_and_no_structural_claim(self):
        process = subprocess.run([sys.executable, "-m", "zipper_algebra.msolve", "--help"],
                                 capture_output=True, text=True, timeout=20)
        assert process.returncode == 0, process.stderr
        assert '--append-target' in process.stdout
        assert 'Higher levels retain' in process.stdout
        assert 'INDUCES NON-SEMISIMPLICITY' not in process.stdout
        assert 'strictly larger' not in process.stdout


class TestSymbolicCompositionEntryPoint:
    @pytest.fixture(autouse=True)
    def resources(self, request):
        self.directory = tempfile.TemporaryDirectory()
        request.addfinalizer(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.pres = Presentation(2, zero_from=2)
        self.atoms = AtomSet([Atom("q", {(1, 0): 1})], [1])
        # g = X2/q and X2 = q*g. The anchor multiplier is 5, not the function q.
        self.lift = Lift(self.atoms,
                         {"g": [((1,), {(0, 0): F(1)}, (), "input",
                                  self.pres.labels().index("ZERO_2_1"), ())]},
                         3, 1, (5, 7, 0), {"g": 0})
        self.skeleton = Skeleton(self.pres.slug(), self.pres.content_hash(), "X_2", (5, 7, 0),
                                 [{"basis_key": "g", "left": (), "right": (), "mult": F(5)}])
        self.correct = [({(1, 0): F(1)}, (0,))]
        self.lift_path = self.root / "lift.json"
        self.skeleton_path = self.root / "skeleton.json"
        self.multiplier_path = self.root / "multipliers.json"
        self.output = self.root / "certificate.txt"
        self.lift.save(self.lift_path)
        self.skeleton.save(self.skeleton_path)
        self.write_multipliers(self.correct)

    def write_multipliers(self, multipliers):
        data = composition.multipliers_to_json(multipliers, self.skeleton, self.atoms, self.pres)
        self.multiplier_path.write_text(json.dumps(data))
        return data

    def arguments(self, *, machine=True):
        return ["--lift", str(self.lift_path), "--skeleton", str(self.skeleton_path),
                "--multipliers", str(self.multiplier_path), "--out", str(self.output),
                "--n", "2", "--zero-from", "2", *(["--json"] if machine else [])]

    def invoke(self, *, machine=True):
        output, errors = io.StringIO(), io.StringIO()
        with redirect_stdout(output), redirect_stderr(errors):
            code = composition.main(self.arguments(machine=machine))
        return code, output.getvalue(), errors.getvalue()

    def test_versioned_multiplier_data_binds_full_presentation_and_occurrences(self):
        data = json.loads(self.multiplier_path.read_text())
        assert data['schema'] == 'zipper-cofactor-multipliers-v1'
        assert data['presentation_sha256'] == hashlib.sha256(self.pres.canonical_text().encode()).hexdigest()
        assert len(data['presentation_sha256']) == 64
        assert composition.load_multipliers(self.multiplier_path, self.skeleton, self.atoms, self.pres) == self.correct

    def test_nonconstant_symbolic_multiplier_is_exactly_verified_before_output(self):
        code, output, errors = self.invoke()
        assert code == 0, (output, errors)
        assert errors == ''
        assert json.loads(output)['status'] == 'verified'
        assert json.loads(output)['exact']['pass']
        assert 'ZERO_2_1' in self.output.read_text()
        assert not list(self.root.glob('*.tmp'))

    def test_valid_text_mode_reports_exact_success(self):
        code, output, errors = self.invoke(machine=False)
        assert code == 0, errors
        assert 'Verified exact identity' in output
        assert self.output.exists()

    def test_anchor_constant_is_rejected_and_preserves_existing_output(self):
        self.write_multipliers([({(0, 0): F(5)}, (0,))])
        self.output.write_text("previous verified certificate\n")
        code, output, errors = self.invoke()
        assert code == 1, errors
        assert json.loads(output)['status'] == 'false-identity'
        assert self.output.read_text() == 'previous verified certificate\n'
        code, output, errors = self.invoke(machine=False)
        assert code == 1, errors
        assert 'Exact identity verification failed' in output
        assert self.output.read_text() == 'previous verified certificate\n'

    def test_matching_all_numeric_check_points_still_does_not_certify_identity(self):
        # P vanishes at every q used by numeric_check; q+P is nevertheless not q.
        coefficients = [F(1)]
        for root in (5, 3, 7, 11, 10):
            nxt = [F(0)] * (len(coefficients) + 1)
            for exponent, coefficient in enumerate(coefficients):
                nxt[exponent] -= root * coefficient
                nxt[exponent+1] += coefficient
            coefficients = nxt
        numerator = {(i, 0): value for i, value in enumerate(coefficients)}
        numerator[(1, 0)] += 1
        multiplier = [(numerator, (0,))]
        candidate = composition.compose(self.lift, self.skeleton, self.pres,
                                          multipliers=multiplier, verbose=False)
        assert composition.numeric_check(candidate, self.pres, verbose=False)
        self.write_multipliers(multiplier)
        code, output, errors = self.invoke()
        assert code == 1, errors
        assert not json.loads(output)['exact']['pass']
        assert not self.output.exists()

    def test_schema_presentation_target_atoms_and_context_mismatches_are_rejected(self, subtests):
        good = json.loads(self.multiplier_path.read_text())
        cases = [[], {**good, "schema": "future-schema"},
                 {**good, "presentation_sha256": "0" * 64},
                 {**good, "target": "X_3"}, {**good, "atoms": {}},
                 {**good, "steps": None}, {**good, "steps": []}]
        for change in (None, {**good["steps"][0], "basis_key": "other"},
                       {**good["steps"][0], "left": [1]},
                       {**good["steps"][0], "right": [2]}):
            cases.append({**good, "steps": [change]})
        for data in cases:
            with subtests.test(data=data):
                self.multiplier_path.write_text(json.dumps(data))
                code, output, errors = self.invoke()
                assert code == 2, (output, errors)
                assert json.loads(output)['status'] == 'invalid-input'
                assert not self.output.exists()

    def test_numeric_atom_metadata_cannot_be_replaced_by_equal_boolean_or_float(self, subtests):
        good = json.loads(self.multiplier_path.read_text())
        # Every original numeric entry is an integer, including the atom's powers,
        # coefficient and search bound. Python equality alone conflates these types.
        for replacement in (True, 1.0):
            for location in ("power", "coefficient", "bound"):
                data = deepcopy(good)
                if location == "bound":
                    data["atoms"]["max_exp"][0] = replacement
                else:
                    index = 0 if location == "power" else 2
                    data["atoms"]["atoms"][0]["coefs"][0][index] = replacement
                self.multiplier_path.write_text(json.dumps(data))
                with subtests.test(replacement=replacement, location=location):
                    with pytest.raises(ValueError, match='denominator atoms'):
                        composition.load_multipliers(self.multiplier_path, self.skeleton,
                                                     self.atoms, self.pres)

    def test_context_metadata_cannot_use_boolean_or_float_generator_indices(self, subtests):
        for side in ("left", "right"):
            self.skeleton.steps[0][side] = (1,)
            good = composition.multipliers_to_json(self.correct, self.skeleton,
                                                   self.atoms, self.pres)
            for replacement in (True, 1.0):
                data = deepcopy(good)
                data["steps"][0][side] = [replacement]
                self.multiplier_path.write_text(json.dumps(data))
                with subtests.test(side=side, replacement=replacement):
                    with pytest.raises(ValueError, match='occurrence'):
                        composition.load_multipliers(self.multiplier_path, self.skeleton,
                                                     self.atoms, self.pres)
            self.skeleton.steps[0][side] = ()

    def test_export_rejects_noninteger_context_indices(self, subtests):
        for side in ("left", "right"):
            for replacement in (True, 1.0):
                self.skeleton.steps[0][side] = (replacement,)
                with subtests.test(side=side, replacement=replacement):
                    with pytest.raises(ValueError, match='word'):
                        composition.multipliers_to_json(self.correct, self.skeleton,
                                                         self.atoms, self.pres)
            self.skeleton.steps[0][side] = ()

    def test_malformed_exact_coefficients_and_denominators_are_rejected(self, subtests):
        good = json.loads(self.multiplier_path.read_text())
        cases = [[[[1, 0], 1]], [[[1, 0], "1"], [[1, 0], "2"]],
                 [[[-1, 0], "1"]], [[[1, 0], "1/0"]], [[[1, 0], "broken"]]]
        for numerator in cases:
            data = deepcopy(good)
            data["steps"][0]["num"] = numerator
            self.multiplier_path.write_text(json.dumps(data))
            with subtests.test(numerator=numerator):
                assert self.invoke()[0] == 2
        for den in ([], [-1], [True]):
            data = deepcopy(good)
            data["steps"][0]["denom"] = den
            self.multiplier_path.write_text(json.dumps(data))
            with subtests.test(den=den):
                assert self.invoke()[0] == 2
        self.multiplier_path.write_text('{"schema":')
        code, output, errors = self.invoke(machine=False)
        assert code == 2
        assert output == ''
        assert 'compose:' in errors
        assert not self.output.exists()

    def test_export_rejects_invalid_or_missing_symbolic_multiplier_values(self, subtests):
        for values in ([], self.correct * 2,
                       [({(1, 0): 1.5}, (0,))], [({(1, 0): F(1)}, (-1,))]):
            with subtests.test(values=values), pytest.raises((TypeError, ValueError)):
                composition.multipliers_to_json(values, self.skeleton, self.atoms, self.pres)

    def test_explicit_empty_multiplier_is_zero_not_a_missing_function(self):
        data = self.write_multipliers([({}, (0,))])
        assert data['steps'][0]['num'] == []
        assert composition.load_multipliers(self.multiplier_path, self.skeleton, self.atoms, self.pres) == [({}, (0,))]
        assert self.invoke()[0] == 1
        assert not self.output.exists()

    def test_atomic_output_failure_cleans_temporary_files_and_keeps_directory(self):
        self.output.mkdir()
        marker = self.output / "keep.txt"
        marker.write_text("keep")
        code, output, errors = self.invoke()
        assert code == 2, errors
        assert marker.read_text() == 'keep'
        assert not list(self.root.glob('*.tmp'))
        self.output = self.root / "missing-parent" / "certificate.txt"
        assert self.invoke()[0] == 2
        assert not self.output.exists()

    def test_cli_requires_explicit_multiplier_file_before_loading_other_inputs(self):
        arguments = self.arguments()
        at = arguments.index("--multipliers")
        del arguments[at:at+2]
        errors = io.StringIO()
        with redirect_stderr(errors), pytest.raises(SystemExit) as raised:
            composition.main(arguments)
        assert raised.value.code == 2
        assert '--multipliers' in errors.getvalue()
        assert not self.output.exists()

    def test_module_entrypoint_uses_the_same_exact_check(self):
        result = subprocess.run([sys.executable, "-m", "zipper_algebra.cofactor.compose",
                                 *self.arguments()], capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stderr
        assert json.loads(result.stdout)['pass']
        assert self.output.exists()
