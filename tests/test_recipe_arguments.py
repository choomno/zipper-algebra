"""Explicit experiment inputs and the actual linear Markov relaxation."""
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction as F
import io
import importlib.util
import pytest
from unittest.mock import patch

from experiments.variants import markov_search as markov
from zipper_algebra import core


class TestMarkovRecipeArguments:
    def test_default_selection_requires_all_four_basis_inputs(self):
        with redirect_stderr(io.StringIO()) as output, pytest.raises(SystemExit) as error:
            markov.main([])
        assert error.value.code == 2
        for flag in ('--framing-small', '--framing-big', '--square-zero-small', '--square-zero-big'):
            assert flag in output.getvalue()

    def test_one_selected_family_does_not_require_other_family_inputs(self):
        with patch.object(markov, 'analyse_variant', return_value={'observed': True}) as analyse, \
                redirect_stdout(io.StringIO()):
            assert markov.main([
                '--only', 'framing', '--framing-small', 'three.basis',
                '--framing-big', 'four.basis',
            ]) == 0
        analyse.assert_called_once_with('FRAMING a=2', 'three.basis', 'four.basis')

    def test_partial_family_is_rejected_before_analysis(self):
        with patch.object(markov, 'analyse_variant') as analyse, \
                redirect_stderr(io.StringIO()) as output, pytest.raises(SystemExit):
            markov.main(['--only', 'square-zero', '--square-zero-small', 'three.basis'])
        analyse.assert_not_called()
        message = output.getvalue().rsplit('error:', 1)[1]
        assert '--square-zero-big' in message
        assert '--framing-small' not in message

    def test_linear_relaxation_rows_have_independent_auxiliary_blocks(self):
        # The larger model is Q<x,y>/(x-2,y-3), with one-dimensional basis (1).
        # Its trace restricts to 1, 2 and 3 on 1, x and y respectively.
        delta, match, patterns, _ = core.build_aho_corasick([('x',), ('y',)], ['x', 'y'])
        larger = {'N': 1, 'mon_idx': {(): 0},
                  'lm_to_tail': {('x',): [(F(2), ())], ('y',): [(F(3), ())]},
                  'ac_d': delta, 'ac_m': match, 'ac_p': patterns}
        rows, smaller_dim, larger_dim = markov.build_markov_system(
            {'N': 1, 'standard': [()]}, larger, [[F(1)]], [[F(1)]], 'x', 'y')
        assert (smaller_dim, larger_dim) == (1, 1)
        assert rows == [[1, -1, 0, 0], [2, 0, -1, 0], [3, 0, 0, -1]]
        rank, kernel = markov.gaussian_rank_nullity(rows, 4)
        assert rank == 3
        assert len(kernel) == 1
        assert all(sum(a*b for a, b in zip(row, kernel[0])) == 0 for row in rows)


class TestIdealClosure:
    @pytest.mark.parametrize("generator, dimension", (({}, 0), ({0: F(1)}, 2), ({1: F(1)}, 1)))
    def test_two_sided_closure_in_the_dual_number_algebra(self, generator, dimension):
        from experiments.variants.square_zero_ideal_nilpotence import two_sided_ideal_basis
        # In Q[x]/(x^2), <0>=0, <1>=A and <x>=Q*x. No cached operators are used.
        action = {'x': {0: {1: F(1)}}}
        ideal = two_sided_ideal_basis(generator, ['x'], action, action)
        assert ideal.dim() == dimension


@pytest.mark.skipif(not importlib.util.find_spec('sympy'), reason='optional symbolic representation recipe')
class TestRepresentationSampleValidation:
    def test_symbolic_and_constant_entries_are_rechecked_against_exact_samples(self):
        import sympy as sp
        from experiments.variants import square_zero_regular_rep_qt as representation
        matrix = sp.Matrix([[representation.q + 2*representation.t]])
        representation._validate_matrix('x', matrix, [(3, 5, {'x': [[F(13)]]})])
        representation._validate_matrix('x', sp.Matrix([[sp.Rational(2, 3)]]),
                                        [(3, 5, {'x': [[F(2, 3)]]})])
        with pytest.raises(RuntimeError, match='VALIDATION FAIL'):
            representation._validate_matrix('x', matrix, [(3, 5, {'x': [[F(12)]]})])
        with pytest.raises(RuntimeError, match='denominator vanishes'):
            representation._validate_matrix(
                'x', sp.Matrix([[1/(representation.q-representation.t)]]),
                [(3, 3, {'x': [[F(0)]]})])


class TestMinimalPolynomialRecipe:
    def test_specialised_basis_determines_report_without_filename_guesses(self, tmp_path, capsys):
        from experiments.variants import sigma_min_poly as sigma
        # F_5<x,y>/(x^2-1,y-x) has basis (1,x), so multiplication by x
        # has minimal polynomial x^2-1 and two distinct eigenvalues.
        basis = tmp_path / 'framing-n2-k17.basis'
        basis.write_text('x1*x1-1\ny1-x1\n')
        basis.with_suffix('.ms').write_text('x1,y1\n5\n')
        result = sigma.process(str(basis), n=2)
        assert result == {
            'path': str(basis), 'n': 2, 'dim': 2,
            'min_poly_deg': 2, 'n_distinct_eigenvalues': 2,
            'factors': [('-1 + x^2', 1)],
        }
        report = capsys.readouterr().out
        assert 'prime p = 5' in report
        assert 'min poly        = -1 + x^2' in report
        assert 'predicted' not in report
        assert 'conjecture' not in report
        assert 'annihilator' not in report
