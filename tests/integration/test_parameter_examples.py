"""Known-answer examples for generic fields and exceptional specializations."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

from zipper_algebra.algebra import AlgebraElement, AlgebraPresentation, CofactorIdentity
from zipper_algebra.fields import FieldElement, RationalFunctionField


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'examples' / 'polynomial_quotients.py'


def load_example():
    spec = importlib.util.spec_from_file_location('polynomial_quotient_examples', SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_installed_parameter_api_outside_checkout_without_sympy(tmp_path):
    code = '''
import json, sys
from zipper_algebra import (RationalFunctionField, FieldElement, FreeAlgebra,
                            AlgebraElement, AlgebraPresentation, CofactorIdentity)
field = RationalFunctionField(('q', 'a', 'b'))
q, a, b = (field.parameter(name) for name in field.parameters)
presentation = AlgebraPresentation.zipper(2, field=field, q=q, untwistings={2: [b, a, 1]})
algebra = presentation.algebra
x = algebra.generator('x1')
p = x*x+a*x+b
identity = CofactorIdentity(presentation.content_hash(), p*(x-1)*(x+q),
                           ((q*p, 0, algebra.one), (-x, 2, algebra.one)))
assert isinstance(q, FieldElement) and isinstance(x, AlgebraElement)
assert isinstance(algebra, FreeAlgebra) and 'sympy' not in sys.modules
print(json.dumps({'verified': identity.verify(presentation), 'labels': presentation.labels}))
'''
    result = subprocess.run([sys.executable, '-I', '-c', code], cwd=tmp_path,
                            capture_output=True, text=True, check=True, timeout=15)
    assert json.loads(result.stdout) == {
        'verified': True, 'labels': ['INV_1_p', 'INV_1_m', 'T_2'],
    }


def test_exact_generic_and_special_fibre_dimensions():
    report = load_example().run_examples()
    assert {name: result['dimension'] for name, result in report.items()} == {
        'three_parameters': 4,
        'prime_field_specialization': 4,
        'repeated_root': 4,
        'leading_coefficient_generic': 4,
        'leading_coefficient_zero': 3,
        'constant_coefficient_generic': 4,
        'constant_coefficient_zero': 2,
    }
    assert report['repeated_root']['nonzero_square_zero_element_verified']
    assert report['constant_coefficient_zero']['removed_invertible_x_factors'] == 2
    assert report['leading_coefficient_zero']['removed_invertible_x_factors'] == 0
    assert all(result['literal_relations_verified'] == 3 for result in report.values())
    assert report['three_parameters']['field']['parameters'] == ['q', 'a', 'b']
    assert report['prime_field_specialization']['field']['characteristic'] == 5
    assert report['prime_field_specialization']['field']['parameters'] == []
    module = load_example()
    for result in report.values():
        certificate = CofactorIdentity.from_json(result['cofactor'])
        assert certificate.presentation_hash == result['presentation_hash']
        presentation = AlgebraPresentation.from_json(result['presentation'])
        assert certificate.verify(presentation)
        witness = result['matrix_witness']
        size = result['dimension']
        assert witness['dimension'] == size
        assert witness['basis_words'] == [[1]*power for power in range(size)]
        matrices = []
        for name in presentation.algebra.generators:
            matrix = module._zero(presentation.algebra.field, size)
            for row, col, value in witness['generator_actions'][name]:
                matrix[row][col] = FieldElement.from_json(value)
            matrices.append(matrix)
        assert all(module._evaluate(relation, matrices) == module._zero(presentation.algebra.field, size)
                   for relation in presentation.relations)
        field = presentation.algebra.field
        assert witness['cyclic_vector'] == [[0, field.one.to_json()]]
        # Replay both inverse equations and the lower-bound witness from saved data.
        # Word matrices act on columns, so the rightmost letter acts first.
        assert module._multiply(field, matrices[0], matrices[1]) == module._identity(field, size)
        assert module._multiply(field, matrices[1], matrices[0]) == module._identity(field, size)
        vector = [field.zero for _ in range(size)]
        for row, entry in witness['cyclic_vector']:
            vector[row] = FieldElement.from_json(entry)
        for basis_index, word in enumerate(witness['basis_words']):
            image = list(vector)
            for letter in reversed(word):
                matrix = matrices[letter-1]
                image = [sum((matrix[row][col]*image[col] for col in range(size)), field.zero)
                         for row in range(size)]
            assert image == [field.scalar(int(row == basis_index)) for row in range(size)]
        if 'square_zero_cofactor' in result:
            nilpotent = AlgebraElement.from_json(result['square_zero_element'])
            square_zero = CofactorIdentity.from_json(result['square_zero_cofactor'])
            assert square_zero.target == nilpotent*nilpotent
            assert square_zero.verify(presentation)
            action = module._evaluate(nilpotent, matrices)
            assert action != module._zero(field, size)
            assert module._multiply(field, action, action) == module._zero(field, size)


def test_example_cli_saves_replayable_identities_with_byte_bound_provenance(tmp_path):
    subprocess.run(['git', 'init', '-q', str(tmp_path)], check=True)
    path = tmp_path/'examples'/'quotients.json'
    result = subprocess.run([sys.executable, str(SCRIPT), '--output', str(path)],
                            cwd=tmp_path, capture_output=True, text=True, check=True, timeout=30)
    assert 'three_parameters: dimension 4' in result.stdout
    assert 'constant_coefficient_zero: dimension 2' in result.stdout
    report = json.loads(path.read_text())
    assert report['leading_coefficient_zero']['dimension'] == 3
    provenance = json.loads(path.with_name(path.name+'.provenance.json').read_text())
    assert set(provenance['outputs']) == {'quotients.json'}
    assert [item['name'] for item in provenance['inputs']] == ['polynomial_quotients.py']
    assert any(name.endswith('/fields.py') for name in provenance['producer']['source_files'])
    assert (path.parent/'.gitignore').read_text() == '*\n!.gitignore\n'
    ignored = subprocess.run(['git', '-C', str(tmp_path), 'check-ignore', str(path)],
                             capture_output=True, text=True, check=True)
    assert ignored.stdout.strip() == str(path)


@pytest.mark.parametrize('expanded', [[], [0], [0, 0, 0], [1]])
def test_zero_and_constant_expanded_products_are_rejected(expanded):
    module = load_example()
    field = RationalFunctionField()
    presentation = AlgebraPresentation.zipper(2, field=field, q=2, untwistings={2: [1, 1]})
    with pytest.raises(ValueError):
        module.certify_two_strand(presentation, 2, [1, 1], expanded)


def test_wrong_product_and_wrong_polynomial_cannot_claim_a_dimension():
    module = load_example()
    field = RationalFunctionField()
    presentation = AlgebraPresentation.zipper(2, field=field, q=2, untwistings={2: [1, 1]})
    # Correct F is x^3+2*x^2-x-2; changing its constant makes the relation fail.
    with pytest.raises(ArithmeticError, match='defining relation'):
        module.certify_two_strand(presentation, 2, [1, 1], [-3, -1, 2, 1])
    # The companion representation remains valid but the proposed spanning identity fails.
    with pytest.raises(ArithmeticError, match='cofactor'):
        module.certify_two_strand(presentation, 2, [2, 1], [-2, -1, 2, 1])
    with pytest.raises(ArithmeticError, match='cofactor'):
        module.certify_two_strand(presentation, 3, [1, 1], [-2, -1, 2, 1])
    # Multiplying the supplied product by x leaves the companion representation
    # unchanged after Laurent cancellation, but no longer matches p*(x-1)*(x+q).
    with pytest.raises(ArithmeticError, match='cofactor'):
        module.certify_two_strand(presentation, 2, [1, 1], [0, -2, -1, 2, 1])
    with pytest.raises(TypeError, match='field'):
        module.certify_two_strand(presentation, RationalFunctionField(('q',)).one,
                                 [1, 1], [-2, -1, 2, 1])


def test_missing_defining_relations_and_wrong_alphabets_are_rejected():
    module = load_example()
    field = RationalFunctionField()
    presentation = AlgebraPresentation.zipper(2, field=field, q=2, untwistings={2: [1, 1]})
    missing = AlgebraPresentation(presentation.algebra, presentation.relations[1:], presentation.labels[1:])
    with pytest.raises(ValueError, match='inverse relations'):
        module.certify_two_strand(missing, 2, [1, 1], [-2, -1, 2, 1])
    three_strands = AlgebraPresentation.zipper(3, field=field, q=2, untwistings={2: [1, 1]})
    with pytest.raises(ValueError, match='alphabet'):
        module.certify_two_strand(three_strands, 2, [1, 1], [-2, -1, 2, 1])


def test_checks_survive_optimized_python(tmp_path):
    code = f"""
import importlib.util
from zipper_algebra.algebra import AlgebraPresentation
from zipper_algebra.fields import RationalFunctionField
spec = importlib.util.spec_from_file_location('example', {str(SCRIPT)!r})
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
presentation = AlgebraPresentation.zipper(2, field=RationalFunctionField(), q=2,
                                           untwistings={{2: [1, 1]}})
try:
    module.certify_two_strand(presentation, 2, [2, 1], [-2, -1, 2, 1])
except ArithmeticError:
    print('false cofactor rejected')
else:
    raise SystemExit('false cofactor accepted')
"""
    result = subprocess.run([sys.executable, '-O', '-c', code], cwd=tmp_path,
                            capture_output=True, text=True, check=True, timeout=15)
    assert result.stdout.strip() == 'false cofactor rejected'


@pytest.mark.parametrize('inverse_label', ['INV_1_p', 'INV_1_m'])
def test_forged_inverse_labels_are_rejected_before_the_laurent_argument(inverse_label):
    module = load_example()
    field = RationalFunctionField()
    genuine = AlgebraPresentation.zipper(2, field=field, q=2, untwistings={2: [1]})
    algebra = genuine.algebra
    x = algebra.generator('x1')
    quartic_factor = (x-algebra.one)*(x+algebra.scalar(2))
    if inverse_label == 'INV_1_p':
        # The companion matrices satisfy these relations and F=2*R_0, but
        # y is free, so the actual algebra is infinite-dimensional.
        relations = (algebra.scalar(field.scalar(1)/2)*quartic_factor,
                     algebra.zero, algebra.zero)
    else:
        # The helper requires both explicit inverse polynomials rather than
        # deriving a missing inverse relation from other defining relations.
        relations = (genuine.relations[0], algebra.zero, genuine.relations[2])
    forged = AlgebraPresentation(algebra, relations, genuine.labels)
    with pytest.raises(ValueError, match=inverse_label+' must be the literal inverse relation'):
        module.certify_two_strand(forged, 2, [1], [-2, 1, 1])



def test_quadratic_relations_with_inverse_labels_do_not_prove_finite_dimension():
    module = load_example()
    field = RationalFunctionField()
    genuine = AlgebraPresentation.zipper(2, field=field, q=1, untwistings={2: [1]})
    algebra = genuine.algebra
    x = algebra.generator('x1')
    relation = x*x-algebra.one
    # The relations constrain only x; y remains free. A two-dimensional
    # representation and the true identity x²-1=R_0 cannot prove spanning.
    forged = AlgebraPresentation(algebra, (relation, relation, algebra.zero), genuine.labels)
    with pytest.raises(ValueError, match='INV_1_p must be the literal inverse relation'):
        module.certify_two_strand(forged, 1, [1], [-1, 0, 1])
