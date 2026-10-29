"""Exact proof graphs with forward references and contextual polynomial reconstruction."""
from fractions import Fraction
import unittest

from zipper_algebra.coefficients import Coef
from zipper_algebra.polynomial import Poly
from zipper_algebra.cofactor.proof import Basis, Proof


def polynomial(terms):
    result = Poly(0)
    for word, coefficient in terms:
        result.add_term(word, Coef(coefficient, 0))
    return result


class TestForwardProofReferences(unittest.TestCase):
    def make_proof(self):
        # g0 = g2; g1 = 2*a*r0; g2 = 3*g1*b + r1.
        # The graph visits indices 0, 2, 1 and expands to 6*a*r0*b + r1.
        return Proof([
            [(Fraction(1), (), "basis", 2, ())],
            [(Fraction(2), (1,), "input", 0, ())],
            [(Fraction(3), (), "basis", 1, (2,)),
             (Fraction(1), (), "input", 1, ())],
        ])

    def test_forward_references_preserve_coefficients_and_ordered_contexts(self):
        proof = self.make_proof()
        expected = {((1,), 0, (2,)): Fraction(6), ((), 1, ()): Fraction(1)}
        self.assertEqual(proof.expand_dict(0), expected)
        self.assertEqual(proof.expand_dict(2), expected)
        self.assertEqual(proof.deps(0), {1, 2})
        self.assertEqual(proof.deps(1), set())
        self.assertEqual(proof.deps(2), {1})

    def test_forward_graph_reconstructs_independently_specified_polynomials(self):
        # r0 = a-1, r1 = b-1, so 6*a*r0*b+r1 = 6*a*a*b-6*a*b+b-1.
        relations = [polynomial([((1,), 1), ((), -1)]),
                     polynomial([((2,), 1), ((), -1)])]
        target = polynomial([((1, 1, 2), 6), ((1, 2), -6), ((2,), 1), ((), -1)])
        middle = polynomial([((1, 1), 2), ((1,), -2)])
        basis = Basis([target, middle, target], 0, ["a", "b"])
        proof = self.make_proof()
        for index in range(3):
            with self.subTest(index=index):
                self.assertTrue(proof.reconstructs(index, basis, relations))
        altered = Basis([target + polynomial([((), 1)]), middle, target], 0, ["a", "b"])
        self.assertFalse(proof.reconstructs(0, altered, relations))


if __name__ == "__main__":
    unittest.main()
