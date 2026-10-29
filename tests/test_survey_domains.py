"""Small exact field/rank regressions and qualified experimental observations."""
from __future__ import annotations

import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from zipper_algebra.linear import matrix_rank, validate_field_characteristic
from experiments.research import theta_min_poly as theta
from experiments.survey import pointalg, semisimple, structural, survey
from experiments.survey.run import run_spec, SURVEY_SCHEMA
from experiments.survey.spec import Spec


class TestExactMatrixRank(unittest.TestCase):
    def test_field_characteristics_are_exactly_checked(self):
        for characteristic in (0, 2, 3, 5, 31, 2147483647, 2147483629):
            with self.subTest(characteristic=characteristic):
                self.assertEqual(validate_field_characteristic(characteristic), characteristic)
        for characteristic in (-1, 1, 4, 9, 25, 341, True, False, 5.0, "5", None):
            with self.subTest(characteristic=characteristic):
                with self.assertRaises(ValueError):
                    validate_field_characteristic(characteristic)

    def test_empty_zero_and_rectangular_matrices(self):
        cases = (([], 0), ([[], []], 0), ([[0, 0], [0, 0]], 0),
                 ([[0, 1, 2], [0, 2, 4]], 1), ([[1, 2, 3]], 1),
                 ([[0, 1], [1, 0], [0, 0], [1, 1]], 2))
        for characteristic in (0, 2, 5):
            for rows, expected in cases:
                with self.subTest(rows=rows, characteristic=characteristic):
                    before = [row[:] for row in rows]
                    self.assertEqual(matrix_rank(rows, characteristic=characteristic), expected)
                    self.assertEqual(rows, before)

    def test_rational_rank_does_not_equal_two_agreeing_modular_ranks(self):
        product = 2147483647 * 2147483629
        self.assertEqual(matrix_rank([[F(product)]]), 1)
        self.assertEqual(semisimple.rank_of([[F(product)]], 0), 1)
        for prime in (2147483647, 2147483629):
            self.assertEqual(matrix_rank([[product]], characteristic=prime), 0)
        self.assertEqual(matrix_rank([[F(1, 2147483647)]]), 1)

    def test_rational_coefficients_and_modular_specialization(self):
        rows = [[F(1, 2), 1], [3, 1]]
        self.assertEqual(matrix_rank(rows), 2)
        self.assertEqual(matrix_rank(rows, characteristic=5), 1)
        self.assertEqual(matrix_rank([[F(1, 2), F(1, 3)], [F(3, 2), 1]]), 1)
        self.assertEqual(semisimple._rank_Fp([[5, 1], [1, 0]], 5), 2)

    def test_invalid_entries_shape_and_field_are_rejected(self):
        for value in (True, 0.5, complex(1), "1", None):
            with self.subTest(value=value):
                with self.assertRaises(TypeError):
                    matrix_rank([[value]])
        with self.assertRaises(ValueError):
            matrix_rank([[1], []])
        with self.assertRaises(ZeroDivisionError):
            matrix_rank([[F(1, 5)]], characteristic=5)
        with self.assertRaises(ValueError):
            matrix_rank([], characteristic=6)
        with self.assertRaises(ValueError):
            semisimple._rank_Fp([[1]], 0)


class TestSurveyCoefficientDomains(unittest.TestCase):
    @staticmethod
    def basis(char):
        scalar = int if char else F
        return [(('x1', 'x1'), [(scalar(1), ('x1', 'x1')), (scalar(4), ())]),
                (('y1',), [(scalar(1), ('y1',)), (scalar(-1), ())])]

    def test_point_reducer_uses_declared_prime_and_full_alphabet(self):
        for char, expected in ((5, 1), (7, 3), (0, F(-4))):
            with self.subTest(char=char):
                with patch.object(pointalg, 'point_gb', return_value=self.basis(char)):
                    model = pointalg.PointAlgebra(Spec('linear', 2), 5, 7, '/unused', char)
                self.assertEqual(model.nf(('x1', 'x1')), {0: expected})
                self.assertEqual(model.standard, [(), ('x1',)])
                self.assertEqual(model.nf(('y1', 'x1')), {1: 1})
                self.assertFalse(model.basis_complete)
                self.assertEqual(model.trace_form(), [[2 % char if char else F(2), 0],
                                                      [0, (2 * expected) % char if char else 2 * expected]])
        with patch.object(pointalg, 'point_gb') as solver:
            with self.assertRaises(ValueError):
                pointalg.PointAlgebra(Spec('linear', 2), 5, 7, '/unused', 6)
            solver.assert_not_called()

    def test_reducer_checks_field_alphabet_and_enumeration_budget(self):
        basis = self.basis(5)
        for ring, char in (('wrong', None), ('Fp', 0), ('Fp', 9), ('Q', 5)):
            with self.subTest(ring=ring, char=char):
                with self.assertRaises(ValueError):
                    theta.make_reducer(basis, ring, char=char)
        for variables in (['x1'], ['x1', 'x1', 'y1']):
            with self.assertRaises(ValueError):
                theta.make_reducer(basis, 'Fp', char=5, variables=variables)
        with self.assertRaisesRegex(RuntimeError, 'budget'):
            theta.make_reducer(basis[:1], 'Fp', char=5,
                               variables=['x1', 'y1'], max_total=5)
        # The default historical modulus remains explicit compatibility behaviour.
        _, nf, _, _, _ = theta.make_reducer(basis, 'Fp')
        self.assertEqual(nf(('x1', 'x1')), {0: theta.P31 - 4})

    def test_krylov_inverse_uses_same_prime_as_reduction(self):
        for ring, char, expected in (('Fp', 5, [4, 0, 1]), ('Q', 0, [F(4), F(0), F(1)])):
            _, nf, standard, zero, norm = theta.make_reducer(self.basis(char), ring, char=char)
            result = theta.min_poly(nf, standard, [(1, ())], ('x1',), ring, zero, norm,
                                    char=char)
            self.assertEqual(result, expected)

    def test_survey_presets_delegate_to_presentation(self):
        for name, untwist, power in (('linear', 'bigelow', 1),
                                     ('square-zero', 'square-zero', 1),
                                     ('framing', 'bigelow', 2)):
            spec = Spec(name, 3)
            self.assertEqual(spec.presentation().untwist, untwist)
            self.assertEqual(spec.presentation().framing_power, power)
        with self.assertRaises(ValueError):
            Spec('unknown', 3)
        with self.assertRaises(ValueError):
            Spec('linear', 1)

    def test_solver_boundary_emits_recorded_characteristic(self):
        real_run, solver_calls = pointalg.subprocess.run, []
        with tempfile.TemporaryDirectory() as directory:
            def produce(args, **kwargs):
                if args[0] not in {pointalg.F4, theta.F4}:
                    return real_run(args, **kwargs)
                solver_calls.append(args)
                output = args[args.index('-o') + 1]
                Path(output).write_text('x1*x1+4\ny1-1\n')
            with patch.object(pointalg.subprocess, 'run', side_effect=produce):
                basis = pointalg.point_gb(Spec('linear', 2), 5, 7, directory, 5)
            self.assertEqual(basis, self.basis(5))
            self.assertIn('5\n', next(Path(directory).glob('*.ms')).read_text())
            self.assertEqual(len(solver_calls), 1)
            with patch.object(theta.subprocess, 'run', side_effect=produce):
                path = theta.build_basis('linear', 2, 5, 7, 5, directory)
                self.assertTrue(Path(path).is_file())


class TestQualifiedSurveyObservations(unittest.TestCase):
    def test_point_records_do_not_merge_different_dimensions_or_claim_structure(self):
        models = [SimpleNamespace(dim=2, trace_form=lambda: [[1, 0], [0, 0]]),
                  SimpleNamespace(dim=3, trace_form=lambda: [[1, 0, 0], [0, 1, 0], [0, 0, 1]])]
        with patch.object(semisimple, 'PointAlgebra', side_effect=models):
            result = semisimple.semisimplicity(Spec('linear', 3), '/unused', char=0)
        self.assertFalse(result['standard_word_counts_agree'])
        self.assertEqual([(p['standard_word_count'], p['trace_rank']) for p in result['points']],
                         [(2, 1), (3, 3)])
        self.assertEqual(result['coefficient_domain'], 'Q')
        self.assertFalse(result['basis_complete'])
        for key in ('dim', 'rank', 'generic_rank', 'radical_dim', 'semisimple'):
            self.assertNotIn(key, result)
        with self.assertRaises(ValueError):
            semisimple.semisimplicity(Spec('linear', 3), '/unused', points=[])

    def test_positive_characteristic_degeneracy_remains_an_observation(self):
        # M_2(F_2) is semisimple but its regular trace form is zero.
        units = [(i, j) for i in range(2) for j in range(2)]
        actions = []
        for i, j in units:
            action = [[0] * 4 for _ in range(4)]
            for column, (k, ell) in enumerate(units):
                if j == k:
                    action[units.index((i, ell))][column] = 1
            actions.append(action)
        gram = [[sum(left[i][j] * right[j][i]
                     for i in range(4) for j in range(4)) % 2
                 for right in actions] for left in actions]
        self.assertEqual(gram, [[0] * 4 for _ in range(4)])
        model = SimpleNamespace(dim=4, trace_form=lambda: gram)
        with patch.object(semisimple, 'PointAlgebra', return_value=model):
            result = semisimple.semisimplicity(Spec('linear', 3), '/unused', points=[(5, 7)], char=2)
        self.assertEqual(result['points'][0]['trace_rank'], 0)
        self.assertEqual(result['coefficient_domain'], 'GF(2)')
        self.assertNotIn('semisimple', result)
        with patch.object(semisimple, 'PointAlgebra',
                          return_value=SimpleNamespace(dim=2, trace_form=lambda: [[0, 1], [0, 0]])):
            result = semisimple.semisimplicity(Spec('linear', 3), '/unused', points=[(5, 7)])
        self.assertFalse(result['points'][0]['trace_form_symmetric'])

    def test_structural_nullspace_respects_prime_and_returns_actual_kernel(self):
        matrix = {(0, 0): 1, (0, 1): 2}
        basis = structural._nullspace_cols_Fp(matrix, 3, 5)
        self.assertEqual(len(basis), 2)
        for vector in basis:
            self.assertEqual((vector.get(0, 0) + 2 * vector.get(1, 0)) % 5, 0)
        self.assertEqual(structural._nullspace_cols_Fp({}, 2, 5), [{0: 1}, {1: 1}])
        for char in (0, 4):
            with self.assertRaises(ValueError):
                structural._nullspace_cols_Fp(matrix, 3, char)
            with self.assertRaises(ValueError):
                structural.pi_kernel_nilpotency(2, 2, 3, '/unused', char)

    def test_structural_reductions_do_not_claim_radical_identification(self):
        def normal_form(word):
            value = 1
            for letter in word:
                value = value * (3 if letter == 'x1' else 2) % 5
            return {0: value} if value else {}
        # Evaluation x=3,y=2 in F5 is an actual one-dimensional algebra model.
        model = SimpleNamespace(dim=1, standard=[()], nf=normal_form)
        messages = []
        with patch.object(structural, 'PointAlgebra', return_value=model):
            result = structural.pi_kernel_nilpotency(2, 2, 3, '/unused', 5, messages.append)
        self.assertTrue(result['observed_square_zero_products'])
        self.assertFalse(result['basis_complete'])
        self.assertEqual(result['reduction_map_nullity'], 0)
        self.assertNotIn('radical_dim', result)
        self.assertNotIn('certifies_rad_eq_ker', result)
        self.assertIn('observed reductions all vanish', messages[-1])

    def test_summary_recomputes_old_unqualified_records_and_aggregation_labels_them(self):
        with tempfile.TemporaryDirectory() as directory:
            spec = Spec('linear', 3)
            path = Path(directory) / f'{spec.name}.json'
            legacy = {'spec': spec.name, 'semisimple': True, 'rank': 123}
            path.write_text(json.dumps(legacy))
            trace = {'points': [{'q': 5, 't': 7, 'standard_word_count': 2, 'trace_rank': 1}],
                     'coefficient_domain': 'Q'}
            with patch('experiments.survey.run.ss.semisimplicity', return_value=trace) as compute:
                result = run_spec(spec, directory, directory, do_basis=False, char=0, log=lambda _: None)
            self.assertEqual(compute.call_count, 1)
            self.assertEqual(result['schema'], SURVEY_SCHEMA)
            self.assertEqual(json.loads(path.read_text()), result)
            display = survey._fmt([legacy, result])
            self.assertIn('legacy-unqualified', display)
            self.assertIn('observed', display)
            self.assertNotIn('123', display)
            self.assertNotIn('semisimple', display)

    @unittest.skipUnless(importlib.util.find_spec('sympy'), 'optional candidate reconstruction uses Sympy')
    def test_candidate_sample_agreement_is_not_basis_certification(self):
        import sympy as sp
        from experiments.survey import reconstruct
        candidate = {'LMS': [('x1',)], 'atoms': [], 'D': sp.Integer(1),
                     'GEN': {('x1',): {(): sp.Integer(1)}}}
        trace = {'points': [], 'coefficient_domain': 'Q'}
        with tempfile.TemporaryDirectory() as directory:
            with patch('experiments.survey.run.ss.semisimplicity', return_value=trace), \
                 patch.object(reconstruct, 'reconstruct_gb', return_value=candidate), \
                 patch.object(reconstruct, 'validate', return_value=(True, 'samples agree')):
                result = run_spec(Spec('linear', 3), directory, directory, log=lambda _: None)
            self.assertTrue(result['basis_candidate']['sample_check_passed'])
            self.assertFalse(result['basis_candidate']['basis_complete'])
            self.assertEqual(result['basis_candidate']['status'], 'candidate')
            self.assertTrue(Path(result['basis_candidate']['artifact']).is_file())


if __name__ == '__main__':
    unittest.main()
