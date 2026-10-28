"""Small regression checks for exploratory recipe frontends and their domains."""
from fractions import Fraction
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from experiments.research import theta3_min_poly
from experiments.variants.compute_dim_aho import hilbert_series
from experiments.variants.trace_space import rank_Q
from experiments.variants.tower_check import _gauss_rank_Q
from zipper_algebra import Presentation, core


class TestExperimentalRecipeBoundaries(unittest.TestCase):
    def test_unit_ideal_has_no_standard_word_including_the_unit(self):
        delta, matches, _, count = core.build_aho_corasick([()], ['x1', 'y1'])
        self.assertEqual(hilbert_series(delta, matches, count, ['x1', 'y1'], 4), [0] * 5)

    def test_avoidance_counts_keep_generators_absent_from_the_leading_words(self):
        delta, matches, _, count = core.build_aho_corasick([('x1', 'x1')], ['x1', 'y1'])
        self.assertEqual(hilbert_series(delta, matches, count, ['x1', 'y1'], 3), [1, 2, 3, 5])
        delta, matches, _, count = core.build_aho_corasick([], ['x1', 'y1'])
        self.assertEqual(hilbert_series(delta, matches, count, ['x1', 'y1'], 3), [1, 2, 4, 8])

    def test_count_cli_requires_a_complete_declared_alphabet(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'basis.txt'
            source.write_text('x1*x1\n')
            command = [sys.executable, '-m', 'experiments.variants.compute_dim_aho', str(source)]
            result = subprocess.run(command + ['--n', '2', '--max-deg', '3'],
                                    text=True, capture_output=True, timeout=20)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('[1, 2, 3, 5]', result.stdout)
            for arguments in ([], ['--vars', 'x1,x1'], ['--vars', ''], ['--n', '1'],
                              ['--vars', 'y1']):
                with self.subTest(arguments=arguments):
                    result = subprocess.run(command + arguments, text=True, capture_output=True,
                                            timeout=20)
                    self.assertEqual(result.returncode, 2, result.stderr)

    def test_rank_recipes_share_exact_rational_arithmetic(self):
        cases = [([], 3, 0), ([[Fraction(2147483647 * 2147483629)]], 1, 1),
                 ([[Fraction(1, 2), 1], [1, 2]], 2, 1), ([[0, 1], [1, 0]], 2, 2)]
        for rows, columns, expected in cases:
            with self.subTest(rows=rows):
                self.assertEqual(rank_Q(rows, columns), expected)
                self.assertEqual(_gauss_rank_Q(rows, columns), expected)

    def test_theta_recipe_emits_its_declared_presentation_before_the_solver_boundary(self):
        # A controlled unit-ideal basis exercises emission/reduction, not a real solver run.
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(theta3_min_poly.subprocess, 'run') as execute, \
                 patch.object(theta3_min_poly.core, 'parse_basis_Fp',
                              return_value=[((), [(1, ())])]):
                self.assertEqual(theta3_min_poly.compute('square-zero', 5, 7, 11,
                                                         'test-producer', directory), (0, None))
            execute.assert_called_once()
            emitted = Path(directory) / 'theta3-square-zero-5-7.ms'
            expected = Path(directory) / 'expected.ms'
            Presentation(3, untwist='square-zero', max_untwist_level=2).emit_ms(expected, 5, 7, 11)
            self.assertEqual(emitted.read_bytes(), expected.read_bytes())
            self.assertEqual(execute.call_args.args[0][0], 'test-producer')

    @unittest.skipUnless(importlib.util.find_spec('sympy'), 'symbolic recipe needs the inference extra')
    def test_symbolic_representation_recipe_uses_the_maintained_presentation(self):
        from experiments.variants.square_zero_regular_rep_qt import write_square_zero_ms
        with tempfile.TemporaryDirectory() as directory:
            actual, expected = Path(directory) / 'actual.ms', Path(directory) / 'expected.ms'
            write_square_zero_ms(str(actual), -2, 3)
            Presentation(3, untwist='square-zero').emit_ms(expected, -2, 3, 0)
            self.assertEqual(actual.read_bytes(), expected.read_bytes())


@unittest.skipUnless(importlib.util.find_spec('sympy'), 'symbolic recipe needs the inference extra')
class TestAllWordCyclicBasis(unittest.TestCase):
    @staticmethod
    def example():
        import sympy as sp
        # a*b acts differently from b*a. Both send the unit to their named
        # single-generator basis vectors, while only a*b reaches e_(a*b).
        a, b = sp.zeros(4), sp.zeros(4)
        a[1, 0] = 1
        a[3, 2] = 1
        b[2, 0] = 1
        return {'a': a, 'b': b}, [(), ('a',), ('b',), ('a', 'b')]

    def test_all_words_with_noncommuting_actions_are_checked_in_order(self):
        from experiments.variants.square_zero_regular_rep_qt import check_cyclic_basis
        matrices, words = self.example()
        self.assertTrue(check_cyclic_basis(matrices, words))

    def test_single_generator_checks_cannot_conceal_a_bad_longer_word(self):
        import sympy as sp
        from experiments.variants.square_zero_regular_rep_qt import check_cyclic_basis
        matrices, words = self.example()
        unit = sp.eye(4)[:, 0]
        for index, generator in enumerate(('a', 'b'), start=1):
            self.assertEqual(matrices[generator] * unit, sp.eye(4)[:, index])
        matrices['a'][3, 2] = 0
        self.assertFalse(check_cyclic_basis(matrices, words))

    def test_unit_word_may_have_a_nonzero_position(self):
        import sympy as sp
        from experiments.variants.square_zero_regular_rep_qt import check_cyclic_basis
        a = sp.zeros(3)
        a[0, 1] = 1
        a[2, 0] = 1
        self.assertTrue(check_cyclic_basis({'a': a}, [('a',), (), ('a', 'a')]))

    def test_missing_unit_and_incompatible_sizes_are_invalid(self):
        import sympy as sp
        from experiments.variants.square_zero_regular_rep_qt import check_cyclic_basis
        with self.assertRaisesRegex(ValueError, 'empty word'):
            check_cyclic_basis({'a': sp.eye(1)}, [('a',)])
        with self.assertRaisesRegex(ValueError, 'matrix size'):
            check_cyclic_basis({'a': sp.eye(2)}, [()])


if __name__ == '__main__':
    unittest.main()
