"""Exact independent oracles for polynomial untwisting input and relations."""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F

import pytest

from zipper_algebra.coefficients import NumericRing, SymbolicRing
from zipper_algebra.linear import matrix_rank
from zipper_algebra.polynomial import Poly
from zipper_algebra.presentation import Presentation
from zipper_algebra.untwisting import UntwistingFamily, UntwistingPolynomial


def numeric_terms(polynomial):
    return {word: coefficient.v for word, coefficient in polynomial.terms.items()}


def convolution(left, right):
    """Ordinary coefficient convolution, independently of free-polynomial arithmetic."""
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def identity(size):
    return [[F(i == j) for j in range(size)] for i in range(size)]


def matrix_product(left, right):
    size = len(left)
    return [[sum((left[i][k] * right[k][j] for k in range(size)), F(0))
             for j in range(size)] for i in range(size)]


def matrix_polynomial(coefficients, matrix):
    """Power sum, independently of the production polynomial's Horner evaluation."""
    size = len(matrix)
    result = [[F(0)] * size for _ in range(size)]
    power = identity(size)
    for coefficient in coefficients:
        for i in range(size):
            for j in range(size):
                result[i][j] += coefficient * power[i][j]
        power = matrix_product(power, matrix)
    return result


def evaluate_free_polynomial(polynomial, generators):
    size = len(generators[0])
    result = [[F(0)] * size for _ in range(size)]
    for word, coefficient in polynomial.terms.items():
        term = identity(size)
        for generator in word:
            term = matrix_product(term, generators[generator - 1])
        for i in range(size):
            for j in range(size):
                result[i][j] += coefficient.v * term[i][j]
    return result


class TestUntwistingPolynomial:
    def test_exact_normalization_preserves_nonmonic_and_zero_constant_inputs(self):
        source = [0, F(3, 2), 2, 0, 0]
        polynomial = UntwistingPolynomial(source)
        source[1] = 99
        assert polynomial.coefficients == (F(0), F(3, 2), F(2))
        assert polynomial.degree == 2
        assert polynomial.constant_term == 0
        assert not polynomial.is_monic
        assert UntwistingPolynomial([1, 1]).is_monic
        assert UntwistingPolynomial([7]).degree == 0
        assert hash(polynomial) == hash(UntwistingPolynomial([0, F(3, 2), 2]))
        assert polynomial == UntwistingPolynomial([0, F(3, 2), 2])
        with pytest.raises(FrozenInstanceError):
            polynomial.coefficients = (F(1),)

    @pytest.mark.parametrize('coefficients', [[], [0], [0, 0, 0]])
    def test_zero_polynomial_is_not_an_untwisting(self, coefficients):
        with pytest.raises(ValueError):
            UntwistingPolynomial(coefficients)

    @pytest.mark.parametrize('coefficients', [
        [True, 1], [1.0, 1], ['1/2', 1], [None, 1], [complex(1), 1],
    ])
    def test_constructor_does_not_guess_exact_coefficients(self, coefficients):
        with pytest.raises((TypeError, ValueError)):
            UntwistingPolynomial(coefficients)

    def test_generator_input_is_consumed_once(self):
        polynomial = UntwistingPolynomial(value for value in [F(1, 3), -2, 1, 0])
        assert polynomial.coefficients == (F(1, 3), F(-2), F(1))
        assert polynomial.to_json() == ['1/3', -2, 1]

    def test_json_round_trip_and_integer_strings(self):
        polynomial = UntwistingPolynomial.from_json(['-3/2', '0', 1, 0])
        assert polynomial == UntwistingPolynomial([F(-3, 2), 0, 1])
        assert polynomial.to_json() == ['-3/2', 0, 1]
        assert UntwistingPolynomial.from_json(polynomial.to_json()) == polynomial

    @pytest.mark.parametrize('value', [
        None, {}, (), '1,1', [], [0], [True, 1], [0.5, 1], ['0.5', 1],
        ['1e2', 1], ['1/0', 1], ['1/2/3', 1], ['', 1], [None, 1],
    ])
    def test_malformed_json_does_not_become_another_polynomial(self, value):
        with pytest.raises((TypeError, ValueError, ZeroDivisionError)):
            UntwistingPolynomial.from_json(value)

    def test_horner_evaluation_has_literal_exact_coefficients(self):
        ring = NumericRing(5, 7, 0)
        operator = Poly.generator_over(1, ring.one)
        result = UntwistingPolynomial([F(1, 3), -2, 1]).evaluate(operator, ring)
        assert numeric_terms(result) == {(): F(1, 3), (1,): F(-2), (1, 1): F(1)}
        constant = UntwistingPolynomial([F(2, 3)]).evaluate(operator, ring)
        assert numeric_terms(constant) == {(): F(2, 3)}

    def test_symbolic_evaluation_keeps_the_supplied_operator(self):
        ring = SymbolicRing()
        operator = Poly.generator_over(1, ring.one).scalar_mul(ring.q)
        result = UntwistingPolynomial([-2, 1]).evaluate(operator, ring)
        assert {word: coefficient.v for word, coefficient in result.terms.items()} == {
            (): {(0, 0): F(-2)}, (1,): {(1, 0): F(1)},
        }

    def test_evaluation_rejects_wrong_argument_and_characteristic(self):
        polynomial = UntwistingPolynomial([1, 1])
        with pytest.raises(TypeError, match='argument'):
            polynomial.evaluate(7, NumericRing(5, 7, 0))
        with pytest.raises(TypeError, match='characteristic'):
            polynomial.evaluate(Poly.generator(1, 0), NumericRing(5, 7, 3))

    def test_rational_coefficients_require_a_valid_residue_specialization(self):
        polynomial = UntwistingPolynomial([F(1, 2), 1])
        ring = NumericRing(5, 7, 3)
        operator = Poly.generator_over(1, ring.one)
        assert numeric_terms(polynomial.evaluate(operator, ring)) == {(): 2, (1,): 1}
        bad_ring = NumericRing(5, 7, 2)
        with pytest.raises(ZeroDivisionError):
            polynomial.evaluate(Poly.generator_over(1, bad_ring.one), bad_ring)


class TestUntwistingFamily:
    def test_family_is_an_immutable_sorted_snapshot(self):
        quadratic = UntwistingPolynomial([2, -3, 1])
        linear = UntwistingPolynomial([-7, 1])
        source = {4: linear, 2: quadratic}
        family = UntwistingFamily(source)
        source[3] = linear
        assert family.items == ((2, quadratic), (4, linear))
        assert family == UntwistingFamily({2: quadratic, 4: linear})
        assert hash(family) == hash(UntwistingFamily({4: linear, 2: quadratic}))
        with pytest.raises(FrozenInstanceError):
            family.items = ()
        assert UntwistingFamily({}).items == ()

    def test_family_json_round_trip_preserves_levels_and_exact_coefficients(self):
        family = UntwistingFamily.from_json({'4': [-1, 1], '2': ['1/3', -2, 1]})
        assert family.to_json() == {'2': ['1/3', -2, 1], '4': [-1, 1]}
        assert UntwistingFamily.from_json(family.to_json()) == family

    @pytest.mark.parametrize('level', [True, 1, 0, -2, 2.0, '2'])
    def test_family_levels_are_exact_supported_integers(self, level):
        with pytest.raises((TypeError, ValueError)):
            UntwistingFamily({level: UntwistingPolynomial([1])})

    @pytest.mark.parametrize('value', [None, [], '2', {2: [1, 1]}])
    def test_family_requires_a_mapping_of_polynomial_objects(self, value):
        with pytest.raises((TypeError, ValueError)):
            UntwistingFamily(value)

    @pytest.mark.parametrize('value', [None, [], {'02': [1]}, {'-2': [1]},
                                             {'2.0': [1]}, {2: [1]}, {'1': [1]}])
    def test_family_json_rejects_ambiguous_or_invalid_levels(self, value):
        with pytest.raises((TypeError, ValueError)):
            UntwistingFamily.from_json(value)


class TestPolynomialPresentation:
    def test_two_strand_free_relator_has_independently_expanded_coefficients(self):
        presentation = Presentation.from_untwistings(
            2, {2: UntwistingPolynomial([2, -3, 1])})
        names, relations, _, tags = presentation.build(NumericRing(5, 7, 0))
        assert names == ['x1', 'y1']
        assert [tag.label(2) for tag in tags] == ['INV_1_p', 'INV_1_m', 'T_2']
        # (x^2-3x+2)(5y-4-x), with no cancellation modulo inverse relations.
        assert numeric_terms(relations[-1]) == {
            (1, 1, 2): F(5), (1, 1, 1): F(-1), (1, 1): F(-1),
            (1, 2): F(-15), (1,): F(10), (2,): F(10), (): F(-8),
        }

    def test_shifted_three_strand_relator_uses_the_correct_lasso_word(self):
        presentation = Presentation.from_untwistings(4, {3: UntwistingPolynomial([2, 1])})
        _, relations, _, tags = presentation.build(NumericRing(2, 7, 0))
        relator = next(rel for rel, tag in zip(relations, tags) if tag.label(4) == 'T_3_j2')
        theta = (2, 3, 3, 2)
        # X_{3,2}=(4*y2*y3-x2*x3)(2*y2-1-x2).
        slider = {(5, 6, 5): 8, (5, 6): -4, (5, 6, 2): -4,
                  (2, 3, 5): -2, (2, 3): 1, (2, 3, 2): 1}
        expected = {theta + word: F(coefficient) for word, coefficient in slider.items()}
        expected.update({word: F(2 * coefficient) for word, coefficient in slider.items()})
        assert numeric_terms(relator) == expected
        assert all(not label.startswith('T_2') for label in presentation.labels())
        assert 'T_3' in presentation.labels() and 'T_3_j2' in presentation.labels()
        assert 'T_4' not in presentation.labels()

    def test_empty_family_and_missing_levels_are_not_filled_in(self):
        bare = Presentation(4, untwist='none')
        empty = Presentation.from_untwistings(4, {})
        assert empty.canonical_text() == bare.canonical_text()
        assert empty.content_hash() == bare.content_hash()
        partial = Presentation.from_untwistings(5, {4: UntwistingPolynomial([1, 1])})
        assert [label for label in partial.labels() if label.startswith('T_')] == ['T_4', 'T_4_j2']

    def test_nonmonic_input_is_not_implicitly_normalized(self):
        one = Presentation.from_untwistings(2, {2: UntwistingPolynomial([1, 1])})
        two = Presentation.from_untwistings(2, {2: UntwistingPolynomial([2, 2])})
        ring = NumericRing(5, 7, 0)
        first = numeric_terms(one.relations(ring)[-1])
        assert numeric_terms(two.relations(ring)[-1]) == {word: 2 * c for word, c in first.items()}
        assert one.content_hash() != two.content_hash()
        assert one.slug() != two.slug()

    def test_raw_zero_constant_input_is_retained(self):
        presentation = Presentation.from_untwistings(2, {2: UntwistingPolynomial([0, 1])})
        relation = numeric_terms(presentation.relations(NumericRing(5, 7, 0))[-1])
        assert relation == {(1, 2): F(5), (1,): F(-4), (1, 1): F(-1)}

    def test_levels_above_the_strand_count_are_rejected(self):
        with pytest.raises(ValueError):
            Presentation.from_untwistings(3, {4: UntwistingPolynomial([1, 1])})

    @pytest.mark.parametrize('options', [
        {'untwist': 'bigelow'}, {'framing_power': 2}, {'theta_power': 2},
        {'max_untwist_level': 2},
    ])
    def test_custom_family_cannot_mix_with_preset_controls(self, options):
        with pytest.raises((TypeError, ValueError)):
            Presentation(3, untwistings=UntwistingFamily({2: UntwistingPolynomial([1, 1])}),
                         **({'untwist': 'none'} | options))

    def test_constructor_rejects_an_unvalidated_family_object(self):
        with pytest.raises(TypeError, match='UntwistingFamily'):
            Presentation(3, untwist='none', untwistings={2: UntwistingPolynomial([1, 1])})

    def test_custom_relations_precede_explicit_zero_and_extra_relations(self):
        presentation = Presentation.from_untwistings(
            4, {3: UntwistingPolynomial([1, 1]), 2: UntwistingPolynomial([-3, 1])},
            max_zipper_level=2, zero_from=4, extra=('X_4',))
        bare_count = Presentation(4, untwist='none', max_zipper_level=2).n_relations()
        assert presentation.labels()[bare_count:] == [
            'T_2', 'T_2_j2', 'T_2_j3', 'T_3', 'T_3_j2', 'ZERO_4_1', 'EXTRA[X_4]',
        ]
        assert 'Z_3' not in presentation.labels()

    def test_custom_presentation_emits_the_literal_solver_input(self, tmp_path):
        presentation = Presentation.from_untwistings(
            2, {2: UntwistingPolynomial([2, -3, 1])})
        output = tmp_path / 'custom.ms'
        assert presentation.emit_ms(str(output), 5, 7, 0, provenance=False) == 3
        assert output.read_bytes() == (
            b'x1,y1\n0\nx1*y1-1,\ny1*x1-1,\n'
            b'-x1*x1*x1+5*x1*x1*y1-x1*x1-15*x1*y1+10*x1+10*y1-8\n')

    def test_existing_preset_identity_is_unchanged(self):
        assert Presentation(4).content_hash() == '091bbb92271ab870'
        assert Presentation(4).n_relations() == 31

    @pytest.mark.parametrize('coefficients', [
        [1], [-3, 1], [4, -4, 1], [F(3, 2), F(-7, 3), 1], [-5, 2, 0, 1],
    ])
    def test_two_strand_companion_representation_and_laurent_upper_bound(self, coefficients):
        q = F(5)
        polynomial = UntwistingPolynomial(coefficients)
        presentation = Presentation.from_untwistings(2, {2: polynomial})
        _, relations, _, _ = presentation.build(NumericRing(int(q), 7, 0))
        # The quotient polynomial is f=p(x)(x-1)(x+q), monic with unit constant.
        f = convolution(coefficients, [-q, q - 1, F(1)])
        size = polynomial.degree + 2
        assert len(f) == size + 1 and f[-1] == 1 and f[0] != 0
        x = [[F(0)] * size for _ in range(size)]
        for column in range(size - 1):
            x[column + 1][column] = F(1)
        for row in range(size):
            x[row][-1] = -f[row]
        # f(x)=0 implies x^{-1}=-(f_1+...+f_N*x^{N-1})/f_0.
        y = matrix_polynomial([-c / f[0] for c in f[1:]], x)
        assert matrix_product(x, y) == identity(size)
        assert matrix_product(y, x) == identity(size)
        for relation in relations:
            assert evaluate_free_polynomial(relation, [x, y]) == [[F(0)] * size for _ in range(size)]
        # Powers 1,x,...,x^{N-1} applied to e_0 are the coordinate basis: lower bound N.
        unit = [F(1)] + [F(0)] * (size - 1)
        orbit = []
        for _ in range(size):
            orbit.append(unit)
            unit = [sum((x[i][j] * unit[j] for j in range(size)), F(0)) for i in range(size)]
        assert orbit == identity(size)
        assert matrix_rank(orbit) == size
        # In the Laurent quotient by xy=yx=1, x*T_2=-f(x): the monic upper bound N.
        cleared = {}
        for word, coefficient in relations[-1].terms.items():
            degree = 1 + word.count(1) - word.count(2)
            cleared[degree] = cleared.get(degree, F(0)) + coefficient.v
        assert {degree: c for degree, c in cleared.items() if c} == {
            degree: -c for degree, c in enumerate(f) if c}
