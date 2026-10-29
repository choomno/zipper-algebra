"""Denominator data and heuristic inference have independent import boundaries."""
import json
import subprocess
import sys


class TestReconstructionImportBoundaries:
    def run_python(self, source):
        result = subprocess.run([sys.executable, "-I", "-c", source],
                                capture_output=True, text=True, timeout=20)
        assert result.returncode == 0, result.stderr
        return json.loads(result.stdout)

    def test_denominator_data_does_not_import_inference_or_sympy(self):
        report = self.run_python('''
import json, sys
from zipper_algebra import Atom, AtomSet
family = AtomSet([Atom("q+t", {(1, 0): 1, (0, 1): 1})], [2])
print(json.dumps({"value": family.eval_denom((2,), 2, 3),
                  "inference": "zipper_algebra.cofactor.inference" in sys.modules,
                  "sympy": "sympy" in sys.modules}))
''')
        assert report == {'value': 25, 'inference': False, 'sympy': False}

    def test_inference_can_be_imported_before_the_factor_data(self):
        report = self.run_python('''
import json, sys
from zipper_algebra.cofactor.inference import infer
from zipper_algebra.cofactor.atoms import AtomSet
family, report = infer({}, verbose=False)
print(json.dumps({"family": family.to_json(),
                  "is_atom_set": isinstance(family, AtomSet),
                  "sympy": "sympy" in sys.modules}))
''')
        assert report['family']['atoms'] == []
        assert report['is_atom_set']
        assert not report['sympy']
