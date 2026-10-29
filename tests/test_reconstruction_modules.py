"""Denominator data and heuristic inference have independent import boundaries."""
import json
import subprocess
import sys
import unittest


class TestReconstructionImportBoundaries(unittest.TestCase):
    def run_python(self, source):
        result = subprocess.run([sys.executable, "-I", "-c", source],
                                capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stderr)
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
        self.assertEqual(report, {"value": 25, "inference": False, "sympy": False})

    def test_inference_can_be_imported_first_and_compatibility_entry_agrees(self):
        report = self.run_python('''
import json, sys
from zipper_algebra.cofactor.inference import infer
from zipper_algebra.cofactor.atoms import AtomSet, infer as compatibility_infer
direct, report = infer({}, verbose=False)
compatibility, other_report = compatibility_infer({}, verbose=False)
print(json.dumps({"direct": direct.to_json(), "compatibility": compatibility.to_json(),
                  "same_report": report == other_report,
                  "sympy": "sympy" in sys.modules}))
''')
        self.assertEqual(report["direct"], report["compatibility"])
        self.assertEqual(report["direct"]["atoms"], [])
        self.assertTrue(report["same_report"])
        self.assertFalse(report["sympy"])
