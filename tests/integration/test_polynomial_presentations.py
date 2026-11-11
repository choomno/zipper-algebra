"""Explicit polynomial presentations bind and verify exact certificates."""
import json
import subprocess
import sys

import pytest

from zipper_algebra import AtomSet, Cofactor, Presentation, UntwistingFamily, UntwistingPolynomial
from zipper_algebra.cofactor.verify import check, pit_check


def test_constant_untwisting_certificate_round_trip(tmp_path):
    presentation = Presentation.from_untwistings(2, {2: UntwistingPolynomial([1])})
    certificate = Cofactor.from_presentation(
        presentation, 'X_2', AtomSet([]),
        [((), 'T_2', (), {(0, 0): 1}, ())], 'polynomial-input-test')
    path = tmp_path / 'certificate.txt'
    certificate.save(path)
    loaded = Cofactor.load(path, pres=presentation)
    assert check(loaded, presentation, 'X_2', verbose=False)['pass']
    assert pit_check(loaded, presentation, 'X_2', verbose=False)['pass']
    assert loaded.presentation_hash == presentation.content_hash()
    # The same abstract quotient with a different literal relation is not the
    # same certificate input: T_2 is now 2*X_2.
    different = Presentation.from_untwistings(2, {2: UntwistingPolynomial([2])})
    with pytest.raises(ValueError, match='hash'):
        check(Cofactor.load(path), different, 'X_2', verbose=False)


def test_public_polynomial_imports_and_immutable_family():
    polynomial = UntwistingPolynomial([2, -3, 1])
    assert UntwistingFamily({2: polynomial}).to_json() == {'2': [2, -3, 1]}


def test_installed_polynomial_interface_outside_checkout(tmp_path):
    result = subprocess.run([sys.executable, '-I', '-c', '''
import json, sys
from zipper_algebra import Presentation, UntwistingPolynomial, UntwistingFamily
family = UntwistingFamily.from_json({'2': [2, -3, 1]})
p = Presentation.from_untwistings(2, dict(family.items))
assert 'sympy' not in sys.modules
print(json.dumps({'labels': p.labels(), 'coefficients': family.to_json()}))
'''], cwd=tmp_path, text=True, capture_output=True, check=True, timeout=15)
    assert json.loads(result.stdout) == {
        'labels': ['INV_1_p', 'INV_1_m', 'T_2'],
        'coefficients': {'2': [2, -3, 1]},
    }
