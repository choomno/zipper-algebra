"""Exact denominator arithmetic and small functional reconstruction examples."""
from fractions import Fraction as F
import importlib.util
import random
import runpy
import sys
from pathlib import Path

import pytest

from zipper_algebra.cofactor import atoms as a
from zipper_algebra.cofactor.inference import infer
from zipper_algebra.cofactor import inference as i
from zipper_algebra.cofactor import univariate as u


HAS_SYMPY = importlib.util.find_spec('sympy') is not None


def grid_values(denominator, qs=range(2, 10), ts=range(11, 19)):
    return [((q, t), F(1, denominator(q, t))) for q in qs for t in ts]


class TestAtomArithmetic:

    def test_integer_polynomial_validation(self):
        for name in (None, '', ' ', '1'):
            with pytest.raises(ValueError):
                a.Atom(name, {(0, 0): 1})
        for coefs in ({(1,): 1}, {(1, -1): 1}, {(True, 0): 1},
                      {('q', 0): 1}, {(1, 0): True}, {(1, 0): 1.0},
                      {(1, 0): F(1, 2)}, {}):
            with pytest.raises(ValueError):
                a.Atom('bad', coefs)
        atom = a.Atom('q', {(1, 0): F(1), (0, 0): 0})
        assert atom.coefs == {(1, 0): 1}
        atom.coefs = []
        with pytest.raises(ValueError):
            atom.validate()

    def test_atom_labels_reject_unbalanced_parentheses(self):
        for name in ('q)', 'q(', ')('):
            with pytest.raises(ValueError):
                a.Atom(name, {(1, 0): 1})

    def test_balanced_atom_labels_round_trip_as_denominator_factors(self):
        for name in ('g(t)', '(q*t+1)', 'q^2', 'q*t'):
            family = a.AtomSet([a.Atom(name, {(1, 0): 1, (0, 1): 1})])
            for exponent in (1, 2):
                written = family.denom_str((exponent,))
                assert family.parse_denom_str(written) == (exponent,)
                assert family.eval_denom((exponent,), 2, 3) == 5 ** exponent

    def test_polynomial_values_persistence_and_sign_units(self):
        atom = a.Atom('difference', {(1, 0): 1, (0, 1): -1})
        negative = a.Atom('opposite', {(1, 0): -1, (0, 1): 1})
        assert atom == negative
        assert hash(atom) == hash(negative)
        assert atom != 'q-t'
        assert atom.eval(F(3, 2), F(1, 2)) == 1
        assert atom.bidegree() == (1, 1)
        assert a.Atom.from_json(atom.to_json()) == atom
        assert repr(atom) == 'Atom(difference)'
        assert a.poly_str({(2, 3): 2, (1, 1): -1, (0, 1): 3, (0, 0): -4}) == '2*q^2*t^3-q*t+3*t-4'
        assert a.poly_str({}) == '0'
        assert a.poly_str({(1, 0): 1, (0, 2): -1}) == 'q-t^2'

    def test_atom_sets_preserve_denominator_positions(self, tmp_path):
        q = a.Atom('q', {(1, 0): 1})
        t = a.Atom('t', {(0, 1): 1})
        compound = a.Atom('q*t+1', {(1, 1): 1, (0, 0): 1})
        family = a.AtomSet([q, t, compound], [2, 1, 2], 'initial')
        assert len(family) == 3
        assert family.zero() == (0, 0, 0)
        assert family.eval_denom((2, 0, 1), 2, 3) == 28
        assert family.denom_str((2, 0, 1)) == 'q^2*(q*t+1)'
        assert family.denom_str((0, 1, 2)) == 't*(q*t+1)^2'
        assert family.denom_str((0, 0, 0)) == '1'
        assert family.parse_denom_str('q*q*(q*t+1)') == (2, 0, 1)
        assert family.parse_denom_str('q^2*(q*t+1)^2') == (2, 0, 2)
        assert family.parse_denom_str('q^0') == (0, 0, 0)
        assert family.parse_denom_str('') == (0, 0, 0)
        assert family.parse_denom_str('1') == (0, 0, 0)
        assert a.AtomSet.lcm((2, 0, 1), (1, 3, 0)) == (2, 3, 1)
        assert a.AtomSet.mul((2, 0, 1), (1, 3, 0)) == (3, 3, 1)
        assert a.AtomSet([q], [2]).candidates() == [(0,), (1,), (2,)]
        extended = family.extended(
            [a.Atom('minus_q', {(1, 0): -1}),
             a.Atom('q+1', {(1, 0): 1, (0, 0): 1})], 4, 'extended')
        assert extended.names() == ['q', 't', 'q*t+1', 'q+1']
        assert extended.max_exp == [2, 1, 2, 4]
        assert extended.provenance == 'extended'
        assert family.extended([]).provenance == 'initial'
        assert 'box=[2, 1, 2]' in repr(family)
        path = tmp_path / 'atoms.json'
        family.save(path)
        restored = a.AtomSet.load(path)
        assert restored.to_json() == family.to_json()
        assert a.AtomSet.from_json({'atoms': [q.to_json()]}).max_exp == [3]

    def test_polynomial_display_uses_coefficients_and_preserves_factor_powers(self):
        declared = a.AtomSet([
            a.Atom('q', {(1, 0): 1}),
            a.Atom('t-q', {(0, 1): 1, (1, 0): -1}),
            a.Atom('t+q^2', {(0, 1): 1, (2, 0): 1}),
        ])
        exponents = (2, 1, 2)
        expected = 'q^2*(-q+t)*(q^2+t)^2'
        assert declared.denom_polynomial_str(exponents) == expected
        assert declared.denom_str(exponents) == 'q^2*(t-q)*(t+q^2)^2'
        assert declared.denom_polynomial_str(declared.zero()) == '1'
        assert a.AtomSet([]).denom_polynomial_str(()) == '1'
        family = a.AtomSet([
            a.Atom('t+5', {(1, 0): 1, (0, 0): -1}),
            a.Atom('q', {(1, 1): -2, (0, 0): 3}),
            a.Atom('constant', {(0, 0): -2}),
        ])
        assert family.denom_polynomial_str((1, 2, 1)) == '(q-1)*(-2*q*t+3)^2*(-2)'
        assert family.denom_str((1, 2, 1)) == '(t+5)*q^2*constant'
        assert family.eval_denom((1, 2, 1), 2, 3) == -162
        with pytest.raises(ValueError):
            family.denom_polynomial_str((1,))

    def test_malformed_families_and_denominators_are_rejected(self):
        q = a.Atom('q', {(1, 0): 1})
        for values in (['q'], [q, q]):
            with pytest.raises(ValueError):
                a.AtomSet(values)
        family = a.AtomSet([q])
        for exponents in ((-1,), (True,), (1.0,), (), '1'):
            with pytest.raises(ValueError):
                family.validate_exponents(exponents)
        for pair in [((1,), (1, 2)), ((-1,), (0,)), ('1', (1,))]:
            with pytest.raises(ValueError):
                a.AtomSet.lcm(*pair)
        for text in ('q)', '(q', '(q)^-1', 'q^-1', 'unknown', 'q**q'):
            with pytest.raises(ValueError):
                family.parse_denom_str(text)


class TestUnivariateDenominatorArithmetic:

    def test_division_exact_and_with_remainder(self):
        quotient, remainder = u._poly_divmod([-1, 0, 1, 0], [-1, 1, 0])
        assert quotient == [1, 1, 0]
        assert remainder == []
        assert u._poly_divmod([1], [1, 1]) == ([], [1])
        assert u._poly_divmod([1, 0, 1], [1, 1]) == ([-1, 1], [2])
        with pytest.raises(ZeroDivisionError):
            u._poly_divmod([1], [0, 0])
        assert u._poly_eval([F(1, 2), 2, 3], 2) == F(33, 2)

    def test_linear_systems_detect_missing_information(self):
        assert u._solve([[F(0), F(1), F(3)], [F(1), F(0), F(2)]], 2) == [2, 3]
        assert u._solve([[F(0), F(1)]], 1) is None
        assert u._solve([[F(0), F(0)]], 1) is None
        assert u._solve([[F(1), F(2)], [F(2), F(4)]], 1) == [2]

    def test_samples_are_exact_and_distinct(self):
        for value in (True, 0.5, '1'):
            with pytest.raises(ValueError):
                u.exact_fraction(value)
        assert u.normalise_points([(2, 4), (1, 1), (2, 4)], 1) == [(F(1), F(1)), (F(2), F(4))]
        assert u.normalise_points([((1, 2), 3)], 2) == [((F(1), F(2)), F(3))]
        for points, dimension in [([1], 1), ([(1, 2, 3)], 1),
                                  ([(1, 2)], 2), ([((1,), 2)], 2),
                                  ([(1, 2), (1, 3)], 1)]:
            with pytest.raises(ValueError):
                u.normalise_points(points, dimension)

    def test_fraction_reduction_and_rational_reconstruction(self):
        assert u._reduce_univariate_fraction([0, 0], [2, 2]) == ([0], [1])
        assert u._reduce_univariate_fraction([-2, 0, 2], [-2, 2]) == ([1, 1], [1])
        with pytest.raises(ValueError):
            u._reduce_univariate_fraction([1], [0])
        for holdout in (0, -1, True, 1.5):
            with pytest.raises(ValueError):
                u.reconstruct_univariate([], holdout)
        assert u.reconstruct_univariate([(1, 1)]) is None
        points = [(x, F(x + 2, x - 1)) for x in range(2, 11)]
        numerator, denominator = u.reconstruct_univariate(points)
        assert numerator == [2, 1]
        assert denominator == [-1, 1]
        assert u.reconstruct_univariate([(x, 0) for x in range(7)]) == ([0], [1])
        assert u.reconstruct_univariate([(x, int(x == 6)) for x in range(7)]) is None

    def test_specialization_and_candidate_order(self):
        polynomial = {(2, 1): 1, (1, 0): 2, (0, 2): -1}
        assert u._specialise(polynomial, 't', 3) == [-9, 2, 3]
        assert u._specialise(polynomial, 'q', 2) == [4, 4, -1]
        assert u._specialise({(1, 1): 1, (1, 0): -1}, 't', 1) == [0]
        assert u._specialise({}, 'q', 2) == [0]
        candidates = i._candidate_polys((1, 0), 2)
        expected = {'q', 'q-1', 'q+1', 'q-2', 'q+2', '2*q-1', '2*q+1'}
        assert {a.poly_str(p) for p in candidates} == expected
        sizes = [sum(abs(c) for c in p.values()) for p in candidates]
        assert sizes == sorted(sizes)
        assert i._candidate_polys((0, 0), 2) == []

    def test_primitive_and_readable_univariate_polynomials(self):
        assert i._primitive([]) == ()
        assert i._primitive([0, 0]) == (0, 0)
        assert i._primitive([F(1, 2), -1]) == (-1, 2)
        assert i._fmt_univariate([0, -1, 2, -3], 'q') == '-q+2*q^2-3*q^3'
        assert i._fmt_univariate([0], 'q') == '0'
        assert i._fmt_univariate([2, 1], 't') == '2+t'

    def test_polynomial_line_fits_have_independent_holdouts(self):
        assert i._fit_univariate_exact([(x, x * x - 1) for x in range(7)], 3) == {0: F(-1), 2: F(1)}
        assert i._fit_univariate_exact([(0, 1), (1, 2)], 3) is None
        assert i._fit_univariate_exact([(1, 0), (1, 0), (1, 0), (2, 1), (3, 2)], 4) is None
        assert i._fit_univariate_exact([(x, x * x) for x in range(5)], 1) is None


class TestAtomInference:

    def test_empty_and_unreconstructible_lines(self, capsys):
        family, report = infer({}, max_bideg=(1, 0), coef_bound=1)
        assert family.names() == []
        assert report['lines_tried'] == 0
        assert report['residual_factors'] == []
        assert 'every reconstructed denominator' in capsys.readouterr().out
        points = [((x, 1), int(x == 6)) for x in range(1, 7)]
        family, report = infer({'bad': points}, max_bideg=(1, 0), coef_bound=1, verbose=False)
        assert report['lines_tried'] == 1
        assert report['lines_reconstructed'] == 0

    def test_polynomial_coefficients_need_no_denominator_atoms(self):
        points = [((q, t), q + t) for q in range(2, 9) for t in range(11, 18)]
        family, report = infer({0: points}, max_bideg=(1, 0), coef_bound=1, verbose=False)
        assert family.names() == []
        assert report['lines_reconstructed'] == 14
        assert report['unexplained'] == 0

    def test_quiet_search_has_the_same_exact_factor_result(self):
        family, report = infer({0: grid_values(lambda q, t: q + t)},
                               max_bideg=(1, 1), coef_bound=1, verbose=False)
        assert family.atoms == [a.Atom('q+t', {(1, 0): 1, (0, 1): 1})]
        assert report['unexplained'] == 0

    def test_seed_divisibility_and_multiplicity(self):
        seed = a.AtomSet([a.Atom('q+t', {(1, 0): 1, (0, 1): 1})], [0])
        family, report = infer({0: grid_values(lambda q, t: (q + t) ** 2)},
                               seed, max_bideg=(1, 1), coef_bound=1, verbose=False)
        assert family.names() == ['q+t']
        assert report['multiplicity'] == {'q+t': 2}
        assert report['residual_factors'] == []
        assert family.max_exp == [3]
        assert report['unexplained'] == 0

    def test_greedy_search_finds_independent_factors(self, capsys):
        family, report = infer({0: grid_values(lambda q, t: (q + t) * (q - t))},
                               max_bideg=(1, 1), coef_bound=1)
        assert set(family.atoms) == {
            a.Atom('q+t', {(1, 0): 1, (0, 1): 1}),
            a.Atom('q-t', {(1, 0): 1, (0, 1): -1}),
        }
        assert report['unexplained'] == 0
        assert 'found' in capsys.readouterr().out

    def test_one_line_cannot_support_a_bivariate_factor(self, monkeypatch):
        monkeypatch.setitem(sys.modules, 'sympy', None)
        family, report = infer({0: grid_values(lambda q, t: q, range(1, 8), [0])},
                               max_bideg=(1, 1), coef_bound=1, min_lines=1, verbose=False)
        assert family.names() == []
        assert report['unexplained'] == 1

    def test_unavailable_factorization_remains_inconclusive(self, capsys, monkeypatch):
        points = grid_values(lambda q, t: q + t + 5, range(1, 8), range(1, 4))
        monkeypatch.setitem(sys.modules, 'sympy', None)
        family, report = infer({0: points}, max_bideg=(0, 0), coef_bound=1)
        assert family.names() == []
        assert report['unexplained'] == 3
        assert report['residual_factors'] == []
        output = capsys.readouterr().out
        assert 'factorisation unavailable' in output
        assert 'UNRESOLVED' in output

    def test_bounded_search_remains_available_without_optional_factorization(self, monkeypatch):
        monkeypatch.setitem(sys.modules, 'sympy', None)
        family, report = infer({0: grid_values(lambda q, t: q + t)},
                               max_bideg=(1, 1), coef_bound=1, verbose=False)
        assert family.atoms == [a.Atom('q+t', {(1, 0): 1, (0, 1): 1})]
        assert report['unexplained'] == 0

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_one_parallel_line_does_not_establish_a_repeated_factor(self):
        family, report = infer({0: grid_values(lambda q, t: q, range(1, 8), [0])},
                               max_bideg=(0, 0), coef_bound=1, min_lines=1, verbose=False)
        assert family.names() == []
        assert report['unexplained'] == 1

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_factor_first_separates_repeated_factors_before_bivariate_lifting(self):
        points = grid_values(
            lambda q, t: (q + 5) ** 2 * (q - 3) * (t * t + 7) * (q * q + t + 11),
            range(4, 17), range(17, 30))
        family, report = infer({0: points}, max_bideg=(0, 0), coef_bound=1, verbose=False)
        assert set(family.atoms) == {
            a.Atom('q+5', {(1, 0): 1, (0, 0): 5}),
            a.Atom('q-3', {(1, 0): 1, (0, 0): -3}),
            a.Atom('t^2+7', {(0, 2): 1, (0, 0): 7}),
            a.Atom('q^2+t+11', {(2, 0): 1, (0, 1): 1, (0, 0): 11}),
        }
        assert report['multiplicity']['q+5'] == 2
        assert report['unexplained'] == 0
        assert report['seed'] == []

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_bounded_fallback_keeps_factors_found_before_factorization_failure(self, monkeypatch):
        factor = i.residual_factors
        calls = 0

        def fail_second_pass(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise ArithmeticError('factorization backend failed')
            return factor(*args, **kwargs)

        monkeypatch.setattr(i, 'residual_factors', fail_second_pass)
        family, report = infer({0: grid_values(lambda q, t: q * (q + t))},
                               max_bideg=(1, 1), coef_bound=1, verbose=False)
        assert set(family.atoms) == {
            a.Atom('q', {(1, 0): 1}),
            a.Atom('q+t', {(1, 0): 1, (0, 1): 1}),
        }
        assert report['unexplained'] == 0

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_factor_lifting_repeats_after_resolving_ambiguous_degree_matches(self, monkeypatch):
        points = grid_values(
            lambda q, t: (q + t + 3) * (q + t * t + 5) * (q * q + t * t + 7),
            range(2, 15), range(17, 30))
        monkeypatch.setattr(i, '_candidate_polys',
                            lambda *_: pytest.fail('resolved factors should avoid enumeration'))
        family, report = infer({0: points}, max_bideg=(0, 0), coef_bound=1, verbose=False)
        assert set(family.atoms) == {
            a.Atom('q+t+3', {(1, 0): 1, (0, 1): 1, (0, 0): 3}),
            a.Atom('q+t^2+5', {(1, 0): 1, (0, 2): 1, (0, 0): 5}),
            a.Atom('q^2+t^2+7', {(2, 0): 1, (0, 2): 1, (0, 0): 7}),
        }
        assert report['unexplained'] == 0

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_residual_lifting_respects_requested_minimum_lines(self):
        points = grid_values(lambda q, t: q + t + 5, range(2, 8), range(11, 17))
        family, report = infer({0: points}, max_bideg=(0, 0), coef_bound=1,
                               min_lines=5, verbose=False)
        assert family.atoms == []
        assert report['unexplained'] == 12

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_three_specializations_do_not_establish_a_constant_factor(self):
        points = grid_values(lambda q, t: q + (t - 1) * (t - 2) * (t - 3),
                             range(2, 15), range(1, 10))
        family, report = infer({0: points}, max_bideg=(0, 0), coef_bound=1, verbose=False)
        assert a.Atom('q', {(1, 0): 1}) not in family.atoms
        assert a.Atom('q+t^3-6*t^2+11*t-6', {
            (1, 0): 1, (0, 3): 1, (0, 2): -6, (0, 1): 11, (0, 0): -6}) in family.atoms
        assert report['unexplained'] == 0

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_stripping_arithmetic_failure_propagates(self, monkeypatch):
        def fail_division(*_args):
            raise ArithmeticError('division failed')

        monkeypatch.setattr(i, '_poly_divmod', fail_division)
        with pytest.raises(ArithmeticError, match='division failed'):
            infer({0: grid_values(lambda q, t: q)}, max_bideg=(0, 0), verbose=False)

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_bounded_factor_can_disambiguate_a_larger_residual_factor(self):
        points = grid_values(lambda q, t: (q + t) * (q + t + 5))
        family, report = infer({0: points}, max_bideg=(1, 1), coef_bound=1, verbose=False)
        assert set(family.atoms) == {
            a.Atom('q+t', {(1, 0): 1, (0, 1): 1}),
            a.Atom('q+t+5', {(1, 0): 1, (0, 1): 1, (0, 0): 5}),
        }
        assert report['unexplained'] == 0

    def test_optional_factorization_fallback_handles_different_slot_denominators(self, monkeypatch):
        monkeypatch.setitem(sys.modules, 'sympy', None)
        family, report = infer({0: grid_values(lambda q, t: q),
                                1: grid_values(lambda q, t: q + t)},
                               max_bideg=(1, 1), coef_bound=1, verbose=False)
        assert set(family.atoms) == {
            a.Atom('q', {(1, 0): 1}),
            a.Atom('q+t', {(1, 0): 1, (0, 1): 1}),
        }
        assert report['unexplained'] == 0

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    @pytest.mark.parametrize('minimum,coordinates', [(1, 3), (2, 4)])
    def test_small_minimum_does_not_bypass_degree_bound_in_varying_lift(
            self, minimum, coordinates):
        def denominator(q, t):
            coincident = (t - 1) * (t - 2) * (t - 3)
            return q + coincident if coordinates == 3 else q + t + coincident * (t - 4)

        points = grid_values(denominator, range(2, 15), range(1, coordinates + 1))
        family, report = infer({0: points}, max_bideg=(0, 0), coef_bound=1,
                               min_lines=minimum, verbose=False)
        assert family.atoms == []
        assert report['unexplained'] == coordinates

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_residual_lifting_repeats_after_bounded_factors_remove_ambiguity(self):
        points = grid_values(
            lambda q, t: ((q + t + 1) * (q + t - 1) * (q + t * t + 5)
                          * (q * q + t * t + 6) * (q * q + t + 7)),
            range(2, 15), range(17, 30))
        family, report = infer({0: points}, max_bideg=(1, 1), coef_bound=1, verbose=False)
        assert set(family.atoms) == {
            a.Atom('q+t+1', {(1, 0): 1, (0, 1): 1, (0, 0): 1}),
            a.Atom('q+t-1', {(1, 0): 1, (0, 1): 1, (0, 0): -1}),
            a.Atom('q+t^2+5', {(1, 0): 1, (0, 2): 1, (0, 0): 5}),
            a.Atom('q^2+t^2+6', {(2, 0): 1, (0, 2): 1, (0, 0): 6}),
            a.Atom('q^2+t+7', {(2, 0): 1, (0, 1): 1, (0, 0): 7}),
        }
        assert report['unexplained'] == 0

    def test_line_factor_lifting_and_ambiguity(self, capsys):
        factors = {('t', (x, 1)): [x] for x in range(1, 8)}
        factors.update({('q', (x, 1)): [x] for x in range(1, 8)})
        result = i.atoms_from_residual_factors(factors)
        assert result == [a.Atom('q+t', {(1, 0): 1, (0, 1): 1})]
        assert 'lifted to a bivariate atom' in capsys.readouterr().out
        odd = range(1, 12, 2)
        rational = {('t', ((x + 1) // 2, 1)): [x] for x in odd}
        assert i.atoms_from_residual_factors(rational, verbose=False) == [
            a.Atom('2*q+t+1', {(1, 0): 2, (0, 1): 1, (0, 0): 1}),
        ]
        assert i.atoms_from_residual_factors({('t', (1, 1)): [1]}, verbose=False) == []
        conflicting = {('t', (1, 1)): list(range(7)), ('t', (2, 1)): list(range(7))}
        assert i.atoms_from_residual_factors(conflicting, verbose=False) == []
        nonpolynomial = {('t', (2 ** x, 1)): [x] for x in range(7)}
        assert i.atoms_from_residual_factors(nonpolynomial, max_other_deg=1, verbose=False) == []
        constant = {('t', (1,)): list(range(7))}
        assert i.atoms_from_residual_factors(constant, verbose=False) == []

    def test_standalone_example_explains_missing_optional_factorization(self, monkeypatch):
        namespace = runpy.run_path(
            Path(__file__).resolve().parents[2] / 'examples/infer_denominators.py')
        monkeypatch.setattr(importlib.util, 'find_spec', lambda _name: None)
        with pytest.raises(SystemExit, match='Install the inference extra'):
            namespace['main']()

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_standalone_example_discovers_factors_outside_its_initial_search_box(
            self, capsys, request):
        namespace = runpy.run_path(
            Path(__file__).resolve().parents[2] / 'examples/infer_denominators.py')
        state = random.getstate()
        request.addfinalizer(lambda: random.setstate(state))
        family, report = namespace['main']()
        assert random.getstate() == state
        assert set(family.atoms) == {
            a.Atom('q+t+5', {(1, 0): 1, (0, 1): 1, (0, 0): 5}),
            a.Atom('q^2+t+7', {(2, 0): 1, (0, 1): 1, (0, 0): 7}),
        }
        assert report['seed'] == []
        assert report['search'] == {'max_bideg': [1, 1], 'coef_bound': 1}
        assert report['unexplained'] == 0
        assert report['residual_factors']
        assert 'No denominator family or seed was supplied.' in capsys.readouterr().out

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_actual_residual_factorization(self, capsys):
        # (x-1)^2*(x^2+1) has these two irreducible rational factors.
        lines = [('t', 2, [1, -2, 2, -2, 1]), ('q', 3, [-1, 1]), ('t', 4, [1])]
        factors = i.residual_factors(lines)
        assert factors == {('t', (-1, 1)): [2], ('t', (1, 0, 1)): [2], ('q', (-1, 1)): [3]}
        assert 'residual irreducible factors' in capsys.readouterr().out
        assert i.residual_factors([], verbose=False) == {}

    @pytest.mark.skipif(not HAS_SYMPY, reason='optional Sympy inference dependency')
    def test_missing_large_constant_is_identified_by_residual_lifting(self):
        points = grid_values(lambda q, t: q + t + 5)
        family, report = infer({0: points}, max_bideg=(1, 1), coef_bound=1)
        assert family.atoms == [a.Atom('q+t+5', {(1, 0): 1, (0, 1): 1, (0, 0): 5})]
        assert report['unexplained'] == 0
        assert report['residual_factors']
