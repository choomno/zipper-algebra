"""Independent free-word identities and adversarial generic certificate inputs."""
from copy import deepcopy
from dataclasses import FrozenInstanceError
from fractions import Fraction
import hashlib
import json
import operator

import pytest

from zipper_algebra.algebra import AlgebraElement, AlgebraPresentation, CofactorIdentity, FreeAlgebra
from zipper_algebra.coefficients import NumericRing
from zipper_algebra.fields import RationalFunctionField
from zipper_algebra.presentation import Presentation
from zipper_algebra.untwisting import UntwistingPolynomial


@pytest.fixture
def algebra():
    return FreeAlgebra(RationalFunctionField(('q', 'a', 'b')), ('x', 'y'))


def coefficients(element):
    return dict(element.terms)


def quartic_witness():
    field = RationalFunctionField(('q', 'a', 'b'))
    q, a, b = (field.parameter(name) for name in field.parameters)
    presentation = AlgebraPresentation.zipper(2, field=field, q=q, untwistings={2: [b, a, 1]})
    algebra = presentation.algebra
    x, y = (algebra.generator(name) for name in algebra.generators)
    polynomial = x*x + a*x + b
    # An independently expanded ordinary polynomial, rather than evaluation
    # of the same product used on the right of the cofactor identity.
    target = (x**4 + (a+q-1)*x**3 + (b+a*(q-1)-q)*x**2
              + (b*(q-1)-a*q)*x - b*q)
    witness = CofactorIdentity(presentation.content_hash(), target,
                               ((q*polynomial, 0, algebra.one), (-x, 2, algebra.one)))
    return presentation, witness, x, y


def test_free_algebra_domain_and_immutable_snapshots(algebra):
    source = {(1,): 3, (): Fraction(1, 2), (2,): 0}
    polynomial = algebra.element(source)
    source[(1,)] = 99
    assert coefficients(polynomial) == {(): algebra.field.scalar(Fraction(1, 2)), (1,): algebra.field.scalar(3)}
    assert not polynomial.is_zero()
    assert algebra.zero.is_zero()
    assert algebra.zero.algebra == algebra
    assert algebra.scalar(0) == algebra.zero
    assert algebra.one == algebra.scalar(1)
    with pytest.raises(FrozenInstanceError):
        polynomial.terms = ()
    with pytest.raises(TypeError):
        hash(polynomial)
    assert FreeAlgebra(RationalFunctionField(), ()).one.terms[0][0] == ()


@pytest.mark.parametrize('generators', [['x'], ('x', 'x'), ('',), ('1x',), ('λ',), (1,)])
def test_invalid_alphabets_are_rejected(generators):
    with pytest.raises(ValueError):
        FreeAlgebra(RationalFunctionField(), generators)


def test_invalid_parents_and_term_containers(algebra):
    with pytest.raises(TypeError):
        FreeAlgebra('Q', ('x',))
    with pytest.raises(TypeError):
        AlgebraElement('Q', {})
    with pytest.raises(TypeError):
        algebra.element([])
    with pytest.raises(ValueError):
        algebra.generator('z')
    for word in ['x', (0,), (3,), (True,), (1.0,), (-1,)]:
        with pytest.raises(ValueError):
            algebra.element({word: 1})
    with pytest.raises(TypeError):
        algebra.scalar(0.5)
    with pytest.raises(TypeError):
        algebra.scalar(True)


@pytest.mark.parametrize('foreign', [
    FreeAlgebra(RationalFunctionField(('a', 'q', 'b')), ('x', 'y')),
    FreeAlgebra(RationalFunctionField(('q', 'a', 'b'), 3), ('x', 'y')),
    FreeAlgebra(RationalFunctionField(('q', 'a', 'b')), ('y', 'x')),
])
@pytest.mark.parametrize('operation', [operator.add, operator.sub, operator.mul])
def test_full_parent_is_checked_even_for_zero_operands(algebra, foreign, operation):
    assert algebra.zero != foreign.zero
    with pytest.raises(TypeError, match='fields or alphabets'):
        operation(algebra.zero, foreign.zero)
    with pytest.raises(TypeError, match='fields or alphabets'):
        operation(foreign.zero, algebra.one)


def test_word_concatenation_and_independently_expanded_powers(algebra):
    x, y = algebra.generator('x'), algebra.generator('y')
    assert x*y != y*x
    assert coefficients((x+y)**3) == {
        (1, 1, 1): 1, (1, 1, 2): 1, (1, 2, 1): 1, (1, 2, 2): 1,
        (2, 1, 1): 1, (2, 1, 2): 1, (2, 2, 1): 1, (2, 2, 2): 1,
    }
    assert coefficients((x+2*y)*(x-y)) == {(1, 1): 1, (1, 2): -1, (2, 1): 2, (2, 2): -2}
    assert (x+y)**0 == algebra.one
    assert x**4 == algebra.element({(1, 1, 1, 1): 1})
    assert x-x == algebra.zero
    assert coefficients(x+1) == {(): 1, (1,): 1}
    assert 1+x == x+1
    assert coefficients(1-x) == {(): 1, (1,): -1}
    assert Fraction(1, 2)*x == x*Fraction(1, 2)
    q = algebra.field.parameter('q')
    assert x*q == q*x
    assert (-x) == -1*x
    assert x != 7
    with pytest.raises(TypeError):
        x+True
    with pytest.raises(TypeError):
        x*RationalFunctionField(('q',)).zero
    for exponent in [-1, True, Fraction(1, 2)]:
        with pytest.raises(ValueError):
            x**exponent


@pytest.mark.parametrize('operation', [
    lambda x: x+object(), lambda x: x-object(), lambda x: object()-x,
    lambda x: x*object(),
])
def test_unsupported_operands_return_python_type_errors(algebra, operation):
    with pytest.raises(TypeError):
        operation(algebra.generator('x'))


def test_free_algebra_and_element_json_round_trips(algebra):
    assert FreeAlgebra.from_json(json.loads(json.dumps(algebra.to_json()))) == algebra
    q = algebra.field.parameter('q')
    polynomial = algebra.element({(1, 2): q/(q-1), (): Fraction(1, 3)})
    assert AlgebraElement.from_json(json.loads(json.dumps(polynomial.to_json()))) == polynomial
    assert AlgebraElement.from_json(algebra.zero.to_json()) == algebra.zero


@pytest.mark.parametrize('value', [None, [], {},
    {'schema': 'other', 'field': {}, 'generators': []},
    {'schema': 'zipper-free-algebra-v1', 'field': {}, 'generators': ()},
])
def test_malformed_free_algebra_json(value):
    with pytest.raises(ValueError):
        FreeAlgebra.from_json(value)


@pytest.mark.parametrize('mutation', [
    lambda value: value.update(schema='other'),
    lambda value: value.update(extra='ignored?'),
    lambda value: value.pop('algebra'),
    lambda value: value.update(terms={}),
    lambda value: value.update(terms=[None]),
    lambda value: value.update(terms=[[]]),
    lambda value: value['terms'][0].__setitem__(0, (1,)),
    lambda value: value['terms'][0].__setitem__(0, [True]),
    lambda value: value['terms'][0].__setitem__(0, [0]),
    lambda value: value['terms'].append(deepcopy(value['terms'][0])),
    lambda value: value['terms'][0].__setitem__(1, RationalFunctionField(('s',)).one.to_json()),
])
def test_corrupt_free_polynomial_data_are_rejected(algebra, mutation):
    data = algebra.generator('x').to_json()
    mutation(data)
    with pytest.raises((TypeError, ValueError)):
        AlgebraElement.from_json(data)


@pytest.mark.parametrize('value', [None, []])
def test_nonobject_polynomial_json(value):
    with pytest.raises(ValueError):
        AlgebraElement.from_json(value)


def test_n2_parameter_relation_has_independent_word_coefficients():
    field = RationalFunctionField(('q', 'a', 'b'))
    q, a, b = (field.parameter(name) for name in field.parameters)
    presentation = AlgebraPresentation.zipper(2, field=field, q=q, untwistings={2: [b, a, 1]})
    assert presentation.algebra.generators == ('x1', 'y1')
    assert presentation.labels == ('INV_1_p', 'INV_1_m', 'T_2')
    assert coefficients(presentation.relations[-1]) == {
        (): b*(1-q), (2,): b*q, (1,): a*(1-q)-b,
        (1, 2): a*q, (1, 1): 1-q-a, (1, 1, 2): q, (1, 1, 1): -1,
    }


def test_shifted_lasso_and_missing_levels_are_literal():
    field = RationalFunctionField()
    presentation = AlgebraPresentation.zipper(4, field=field, q=2,
        untwistings={3: [2, 1, 0]}, max_zipper_level=2)
    assert [label for label in presentation.labels if label.startswith('T_')] == ['T_3', 'T_3_j2']
    assert 'Z_3' not in presentation.labels
    relation = presentation.relations[presentation.labels.index('T_3_j2')]
    theta = (2, 3, 3, 2)
    slider = {(5, 6, 5): 8, (5, 6): -4, (5, 6, 2): -4,
              (2, 3, 5): -2, (2, 3): 1, (2, 3, 2): 1}
    expected = {theta+word: coefficient for word, coefficient in slider.items()}
    expected.update({word: 2*coefficient for word, coefficient in slider.items()})
    assert coefficients(relation) == expected
    assert AlgebraPresentation.zipper(2, field=field, q=2, untwistings={}).labels == ('INV_1_p', 'INV_1_m')


@pytest.mark.parametrize('characteristic', [0, 7])
def test_specialized_generic_relations_match_exact_numeric_constructor(characteristic):
    field = RationalFunctionField(('q', 'a', 'b'))
    q, a, b = (field.parameter(name) for name in field.parameters)
    generic = AlgebraPresentation.zipper(3, field=field, q=q,
        untwistings={2: [b, a, 1], 3: [Fraction(1, 2), -3, 2]})
    specialized = generic.specialize({'q': 5, 'a': -3, 'b': 2}, characteristic=characteristic)
    legacy = Presentation.from_untwistings(3, {
        2: UntwistingPolynomial([2, -3, 1]), 3: UntwistingPolynomial([Fraction(1, 2), -3, 2])})
    assert specialized.labels == tuple(legacy.labels())
    for new, old in zip(specialized.relations, legacy.relations(NumericRing(5, 7, characteristic)), strict=True):
        assert coefficients(new) == {word: coefficient.v for word, coefficient in old.terms.items()}
    assert Presentation(4).content_hash() == '091bbb92271ab870'


def test_raw_nonmonic_and_zero_constant_polynomials_are_not_normalized():
    field = RationalFunctionField()
    first = AlgebraPresentation.zipper(2, field=field, q=2, untwistings={2: [0, 1]})
    second = AlgebraPresentation.zipper(2, field=field, q=2, untwistings={2: [0, 2]})
    assert second.relations[-1] == 2*first.relations[-1]
    assert second.content_hash() != first.content_hash()
    assert coefficients(first.relations[-1]) == {(1, 2): 2, (1,): -1, (1, 1): -1}


@pytest.mark.parametrize('untwistings', [[], {1: [1]}, {3: [1]}, {True: [1]},
    {'2': [1]}, {2: []}, {2: [0, 0]}, {2: [0.5, 1]}, {2: [True, 1]}])
def test_bad_polynomial_family_inputs_are_rejected(untwistings):
    with pytest.raises((TypeError, ValueError)):
        AlgebraPresentation.zipper(2, field=RationalFunctionField(), q=2, untwistings=untwistings)


def test_presentation_validates_literal_data_and_hashes(algebra):
    x, y = algebra.generator('x'), algebra.generator('y')
    original = AlgebraPresentation(algebra, (x, y), ('r', 's'))
    restored = AlgebraPresentation.from_json(json.loads(json.dumps(original.to_json())))
    assert restored == original
    assert restored.content_hash() == original.content_hash()
    assert original.content_hash() != AlgebraPresentation(algebra, (y, x), ('r', 's')).content_hash()
    assert original.content_hash() != AlgebraPresentation(algebra, (x, y), ('s', 'r')).content_hash()
    foreign = FreeAlgebra(RationalFunctionField(('a', 'q', 'b')), ('x', 'y'))
    assert original.content_hash() != AlgebraPresentation(foreign,
        (foreign.generator('x'), foreign.generator('y')), ('r', 's')).content_hash()
    q = algebra.field.parameter('q')
    factored = algebra.scalar((q*q-1)/(q-1))
    reduced = algebra.scalar(q+1)
    assert factored == reduced
    assert AlgebraPresentation(algebra, (factored,), ('r',)).content_hash() != \
        AlgebraPresentation(algebra, (reduced,), ('r',)).content_hash()
    with pytest.raises(TypeError):
        AlgebraPresentation(None, (), ())
    for relations in [[x], (7,), (foreign.one,)]:
        with pytest.raises(TypeError):
            AlgebraPresentation(algebra, relations, ('r',))
    for labels in [['r', 's'], ('r',), ('', 's'), ('r', 2), ('r', 'r')]:
        with pytest.raises(ValueError):
            AlgebraPresentation(algebra, (x, y), labels)


@pytest.mark.parametrize('value', [None, [], {},
    {'schema': 'other', 'algebra': {}, 'relations': [], 'labels': []},
    {'schema': 'zipper-algebra-presentation-v1', 'algebra': {}, 'relations': (), 'labels': []},
    {'schema': 'zipper-algebra-presentation-v1', 'algebra': {}, 'relations': [], 'labels': ()},
])
def test_malformed_presentation_json(value):
    with pytest.raises(ValueError):
        AlgebraPresentation.from_json(value)


def test_specialization_keeps_a_zero_relation_and_rejects_stored_poles():
    field = RationalFunctionField(('s',))
    s = field.parameter('s')
    algebra = FreeAlgebra(field, ('x',))
    presentation = AlgebraPresentation(algebra, (s*algebra.generator('x'),), ('s*x',))
    fiber = presentation.specialize({'s': 0})
    assert fiber.relations == (fiber.algebra.zero,)
    assert fiber.labels == ('s*x',)
    pole = AlgebraPresentation(algebra, (algebra.scalar((s*s-1)/(s-1)),), ('r',))
    with pytest.raises(ZeroDivisionError):
        pole.specialize({'s': 1})
    with pytest.raises(ValueError):
        AlgebraPresentation(algebra, (), ()).specialize({})


def test_cofactor_is_an_exact_three_parameter_identity_and_round_trips(tmp_path):
    presentation, witness, _, _ = quartic_witness()
    assert witness.verify(presentation)
    assert CofactorIdentity.from_json(json.loads(json.dumps(witness.to_json()))).verify(presentation)
    output = tmp_path/'cofactor.json'
    witness.save(output)
    assert CofactorIdentity.load(output).verify(presentation)
    record = json.loads(output.with_name(output.name+'.provenance.json').read_text())
    assert record['outputs'][output.name]['sha256'] == hashlib.sha256(output.read_bytes()).hexdigest()
    assert 'src/zipper_algebra/algebra.py' in record['producer']['source_files']
    assert 'src/zipper_algebra/fields.py' in record['producer']['source_files']


def test_rational_function_cofactors_preserve_literal_data_and_detect_false_fractions():
    field = RationalFunctionField(('q', 'a', 'b'))
    q, a, b = (field.parameter(name) for name in field.parameters)
    algebra = FreeAlgebra(field, ('x',))
    x = algebra.generator('x')
    presentation = AlgebraPresentation(algebra, (a/(q-1)*x,), ('r',))
    left = algebra.scalar((q*q-1)/(a*b))
    right = algebra.scalar(b/(q+1))
    identity = CofactorIdentity(presentation.content_hash(), x, ((left, 0, right),))
    assert identity.verify(presentation)
    encoded = json.loads(json.dumps(identity.to_json()))
    restored = CofactorIdentity.from_json(encoded)
    assert restored.to_json() == encoded
    assert restored.presentation_hash == presentation.content_hash()
    assert restored.verify(presentation)
    # The coefficient's numerator and denominator both matter. Neither
    # corrupted expression can be rescued by ideal membership alone: this
    # checker establishes the supplied literal cofactor equality.
    wrong_denominator = algebra.scalar(b/(q+2))
    wrong_numerator = algebra.scalar((b+1)/(q+1))
    for false_right in (wrong_denominator, wrong_numerator):
        assert not CofactorIdentity(presentation.content_hash(), x,
                                    ((left, 0, false_right),)).verify(presentation)
    changed = AlgebraPresentation(algebra, (a/(q+1)*x,), ('r',))
    with pytest.raises(ValueError, match='hash'):
        restored.verify(changed)


def test_symbolic_prime_field_untwisting_has_an_exact_quartic_cofactor():
    field = RationalFunctionField(('q', 'a', 'b'), characteristic=7)
    q, a, b = (field.parameter(name) for name in field.parameters)
    presentation = AlgebraPresentation.zipper(2, field=field, q=q, untwistings={2: [b, a, 1]})
    algebra = presentation.algebra
    x = algebra.generator('x1')
    polynomial = x*x+a*x+b
    # Coefficients of (x²+a*x+b)(x-1)(x+q), expanded in GF(7)(q,a,b).
    quartic = algebra.element({
        (): -b*q,
        (1,): b*(q-1)-a*q,
        (1, 1): b+a*(q-1)-q,
        (1, 1, 1): a+q-1,
        (1, 1, 1, 1): 1,
    })
    identity = CofactorIdentity(presentation.content_hash(), quartic,
                               ((q*polynomial, 0, algebra.one), (-x, 2, algebra.one)))
    assert identity.verify(presentation)
    assert CofactorIdentity.from_json(identity.to_json()).verify(presentation)
    assert all(coefficient.field == field for _, coefficient in quartic.terms)


def test_three_strand_zipper_has_independently_expanded_symbolic_word_coefficients():
    field = RationalFunctionField(('q',))
    q = field.parameter('q')
    presentation = AlgebraPresentation.zipper(3, field=field, q=q, untwistings={})
    relation = presentation.relations[presentation.labels.index('Z_2')]
    # Alphabet x1,x2,y1,y2. Expand
    # (q*y2+1-q-x2-q*y1*y2+x1*x2)(q*y1+1-q-x1)
    # as fifteen distinct free words, without using production multiplication.
    assert coefficients(relation) == {
        (4, 3): q*q, (4,): q*(1-q), (4, 1): -q,
        (3,): q*(1-q), (): (1-q)*(1-q), (1,): q-1,
        (2, 3): -q, (2,): q-1, (2, 1): 1,
        (3, 4, 3): -q*q, (3, 4): -q*(1-q), (3, 4, 1): q,
        (1, 2, 3): q, (1, 2): 1-q, (1, 2, 1): -1,
    }
    assert not any(label.startswith('T_') for label in presentation.labels)


def test_wrong_cofactor_sign_scalar_order_and_presentation_are_rejected():
    presentation, witness, x, y = quartic_witness()
    algebra = presentation.algebra
    first, second = witness.terms
    assert not CofactorIdentity(witness.presentation_hash, witness.target,
                               (first, (x, second[1], algebra.one))).verify(presentation)
    assert not CofactorIdentity(witness.presentation_hash, witness.target,
                               ((2*first[0], first[1], first[2]), second)).verify(presentation)
    assert not CofactorIdentity(witness.presentation_hash, witness.target+x*y-y*x,
                               witness.terms).verify(presentation)
    with pytest.raises(ValueError, match='hash'):
        witness.verify(AlgebraPresentation(algebra, tuple(reversed(presentation.relations)),
                                          tuple(reversed(presentation.labels))))
    with pytest.raises(TypeError, match='AlgebraPresentation'):
        witness.verify(Presentation(2))
    with pytest.raises(ValueError, match='out of range'):
        CofactorIdentity(witness.presentation_hash, algebra.zero,
                         ((algebra.zero, 3, algebra.zero),)).verify(presentation)
    other_algebra = FreeAlgebra(RationalFunctionField(), ('z',))
    with pytest.raises(TypeError, match='different free algebras'):
        CofactorIdentity(witness.presentation_hash, other_algebra.zero, ()).verify(presentation)
    assert CofactorIdentity(witness.presentation_hash, algebra.zero, ()).verify(presentation)


@pytest.mark.parametrize('presentation_hash', [None, 'abcd', 'F'*64, 'g'*64])
def test_invalid_cofactor_hash(presentation_hash, algebra):
    with pytest.raises(ValueError):
        CofactorIdentity(presentation_hash, algebra.zero, ())


def test_cofactor_constructor_rejects_malformed_or_foreign_factors(algebra):
    digest = 'a'*64
    with pytest.raises(TypeError):
        CofactorIdentity(digest, 0, ())
    with pytest.raises(TypeError):
        CofactorIdentity(digest, algebra.zero, [])
    for term in [None, (), (algebra.one, -1, algebra.one),
                 (algebra.one, True, algebra.one), [algebra.one, 0, algebra.one]]:
        with pytest.raises(ValueError):
            CofactorIdentity(digest, algebra.zero, (term,))
    for factor in [7, FreeAlgebra(RationalFunctionField(), ('x',)).one]:
        with pytest.raises(TypeError):
            CofactorIdentity(digest, algebra.zero, ((factor, 0, algebra.one),))


@pytest.mark.parametrize('value', [None, [], {},
    {'schema': 'other', 'presentation_hash': 'a'*64, 'target': {}, 'terms': []},
    {'schema': 'zipper-cofactor-identity-v1', 'presentation_hash': 'a'*64, 'target': {}, 'terms': ()},
])
def test_malformed_cofactor_json(value):
    with pytest.raises(ValueError):
        CofactorIdentity.from_json(value)


@pytest.mark.parametrize('term', [None, [], [1, 2]])
def test_malformed_cofactor_terms_json(algebra, term):
    data = CofactorIdentity('a'*64, algebra.zero, ()).to_json()
    data['terms'] = [term]
    with pytest.raises(ValueError):
        CofactorIdentity.from_json(data)
