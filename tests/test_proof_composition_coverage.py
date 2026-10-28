"""Small exact examples for proof graphs, reduction and contextual substitution."""
from __future__ import annotations

from contextlib import redirect_stdout
from fractions import Fraction as F
import io
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest

from zipper_algebra import Presentation
from zipper_algebra.zipper_in_ms import Coef, NumericRing, Poly
from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor import compose as composition
from zipper_algebra.cofactor.lift import Lift, NeedMoreSamples
from zipper_algebra.cofactor.proof import (
    Basis, Proof, _parse_coef_str, _parse_noex_term, _split_terms,
    leading_term, parse_basis, parse_noex_proof, word_text,
)
from zipper_algebra.cofactor.skeleton import (
    Skeleton, compose_at, evaluate, find_substring, reduce,
)


def poly(terms, characteristic=0):
    value = Poly(characteristic)
    for word, coefficient in terms:
        value.add_term(word, Coef(coefficient, characteristic))
    return value


def values(p):
    return {word: coefficient.v for word, coefficient in p.terms.items()}


class TestProofParsingExamples(unittest.TestCase):
    def test_degree_lex_order_and_zero_basis_rejection(self):
        p = poly([((2,), 9), ((1, 1), -3), ((), 7)])
        self.assertEqual(leading_term(p)[0], (1, 1))
        self.assertEqual(leading_term(p)[1].v, -3)
        self.assertEqual(word_text((), ["a", "b"]), "1")
        self.assertEqual(word_text((2, 1), ["a", "b"]), "b*a")
        for operation in (lambda: leading_term(Poly(0)),
                          lambda: Basis([Poly(0)], 0, ["a", "b"])):
            with self.assertRaises(ValueError):
                operation()

    def test_signed_tokens_preserve_exact_rational_coefficients(self):
        self.assertEqual(_parse_coef_str("-5/7"), F(-5, 7))
        self.assertEqual(_parse_coef_str("8"), F(8))
        self.assertEqual(_split_terms(" , "), [])
        self.assertEqual(_split_terms("a + -2*b - -3*a + 4*b,"),
                         ["+a", "-2*b", "+3*a", "+4*b"])
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/"basis.txt"
            p.write_text("\na - 2*b + -1/2*a,\n1,\n")
            basis = Basis.load(str(p), 0, ["a", "b"])
            self.assertEqual(len(basis), 2)
            self.assertEqual(values(basis.polys[0]), {(1,): F(1, 2), (2,): F(-2)})
            self.assertEqual(values(basis.polys[1]), {(): F(1)})
            p.write_text("a + c\n")
            with self.assertRaisesRegex(ValueError, "unknown variable"):
                parse_basis(str(p), 0, ["a", "b"])

    def test_leading_word_collisions_and_identical_repeats_have_distinct_keys(self):
        p = poly([((2,), 1), ((1,), -1)])
        q = poly([((2,), 1), ((), -1)])
        basis = Basis([p, q, p], 0, ["a", "b"])
        self.assertEqual(len(set(basis.key)), 3)
        self.assertTrue(basis.key[2].endswith("#2"))
        self.assertEqual(basis.index(basis.key[1]), 1)
        with self.assertRaisesRegex(KeyError, "absent"):
            basis.index("absent")

    def test_proof_terms_have_one_reference_and_known_context_letters(self):
        self.assertEqual(_parse_noex_term("-3/2*a*[i2]*b", {"a": 1, "b": 2}),
                         (F(-3, 2), (1,), "input", 2, (2,)))
        for text in ("+a*[0]*[i0]", "+a*[i0]*[0]", "+a", "+c*[i0]"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                _parse_noex_term(text, {"a": 1})
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/"proof.txt"
            path.write_text("\n+2*a*[i0]*b\n")
            proof = Proof.load(str(path), ["a", "b"])
            self.assertEqual(len(proof), 2)
            self.assertEqual(proof.expand(0), [])
            self.assertEqual(parse_noex_proof(str(path), ["a", "b"])[1],
                             [(F(2), (1,), "input", 0, (2,))])


    def test_negative_references_are_invalid_in_the_saved_proof_format(self):
        for reference in ("[i-1]", "[-1]"):
            with self.subTest(reference=reference):
                with self.assertRaises(ValueError):
                    _parse_noex_term("+a*"+reference+"*b", {"a": 1, "b": 2})
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory)/"invalid-proof.txt"
                    path.write_text("+"+reference+"\n")
                    with self.assertRaises(ValueError):
                        Proof.load(str(path), ["a", "b"])


class TestExactProofGraphs(unittest.TestCase):
    def test_diamond_dependencies_cached_expansion_and_term_cancellation(self):
        # g0 = a*r*b; g1 = 2*g0; g2 = g0 + g1; g3 = g2 - 3*g0 = 0.
        proof = Proof([
            [(F(1), (1,), "input", 0, (2,))],
            [(F(2), (), "basis", 0, ())],
            [(F(1), (), "basis", 0, ()), (F(1), (), "basis", 1, ())],
            [(F(1), (), "basis", 2, ()), (F(-3), (), "basis", 0, ())],
        ])
        self.assertEqual(proof.deps(1), {0})
        self.assertEqual(proof.deps(3), {0, 1, 2})
        self.assertEqual(proof.deps(3), {0, 1, 2})
        self.assertEqual(proof.expand(0), [(F(1), (1,), 0, (2,))])
        self.assertEqual(proof.expand(3), [])
        self.assertEqual(proof.expand_all(), [
            [(F(1), (1,), 0, (2,))], [(F(2), (1,), 0, (2,))],
            [(F(3), (1,), 0, (2,))], [],
        ])
        basis = Basis([poly([((1, 2), 1)]), poly([((2, 1), 1)]),
                       poly([((1, 1), 1)]), poly([((2, 2), 1)])], 0, ["a", "b"])
        self.assertEqual(proof.shape(2), (((), "basis", 0, ()), ((), "basis", 1, ())))
        self.assertEqual(proof.shape(2, basis),
                         (((), "basis", basis.key[0], ()), ((), "basis", basis.key[1], ())))

    def test_out_of_range_references_and_cycles_are_rejected(self):
        for index in (-1, 1):
            with self.subTest(index=index), self.assertRaisesRegex(ValueError, "references"):
                Proof([[(F(1), (), "basis", index, ())]]).expand(0)
        cyclic = Proof([[(F(1), (), "basis", 1, ())], [(F(1), (), "basis", 0, ())]])
        self.assertEqual(cyclic.deps(0), {0, 1})
        with self.assertRaisesRegex(ValueError, "cyclic proof"):
            cyclic.expand(0)

    def test_reconstruction_checks_the_actual_polynomial_not_only_the_graph(self):
        relation = poly([((1,), 1), ((), -1)])
        g = poly([((2, 1, 2), 2), ((2, 2), -2)])
        proof = Proof([[(F(2), (2,), "input", 0, (2,))]])
        self.assertTrue(proof.reconstructs(0, Basis([g], 0, ["a", "b"]), [relation]))
        altered = g + poly([((), 1)])
        self.assertFalse(proof.reconstructs(0, Basis([altered], 0, ["a", "b"]), [relation]))


    def test_invalid_input_indices_cannot_alias_the_last_relation(self):
        relation = poly([((1,), 1), ((), -1)])
        basis = Basis([relation], 0, ["a", "b"])
        # In particular -1 would select the final input through Python indexing
        # and incorrectly report reconstruction of this very basis element.
        for index in (-1, 1, True, F(1, 2)):
            proof = Proof([[(F(1), (), "input", index, ())]])
            with self.subTest(index=index), self.assertRaises(ValueError):
                proof.reconstructs(0, basis, [relation])


class TestSmallReductionExamples(unittest.TestCase):
    def test_leftmost_subword_and_characteristic_validation(self):
        self.assertEqual(find_substring((1, 2, 1, 2), (1, 2)), ((), (1, 2)))
        self.assertEqual(find_substring((1, 2), ()), ((), (1, 2)))
        self.assertIsNone(find_substring((1,), (2, 2)))
        self.assertIsNone(find_substring((1, 1), (2,)))
        with self.assertRaisesRegex(ValueError, "char"):
            reduce(poly([((1,), 1)], 5), Basis([poly([((1,), 1)])], 0, ["a"]))

    def test_finite_field_reduction_and_an_irreducible_remainder(self):
        basis = Basis([poly([((1,), 2), ((), -2)], 5)], 5, ["a", "b"])
        steps, remainder = reduce(poly([((2, 1, 2), 3), ((2,), 4)], 5), basis)
        self.assertEqual(len(steps), 1)
        self.assertEqual(steps[0][:3], ((2,), 0, (2,)))
        self.assertEqual(steps[0][3].v, 4)  # 3/2 = 4 in F_5.
        self.assertEqual(values(remainder), {(2, 2): 3, (2,): 4})

    def test_unit_ideal_reduces_every_word(self):
        basis = Basis([poly([((), 2)])], 0, ["a", "b"])
        target = poly([((1, 2), 3), ((2,), -4), ((), 5)])
        steps, remainder = reduce(target, basis)
        self.assertFalse(remainder.terms)
        self.assertEqual(evaluate([(s[3].v, s[0], s[1], s[2]) for s in steps],
                                  basis.polys, 0).terms, target.terms)

    def test_skeleton_closure_reidentification_and_cancelling_occurrences(self):
        p = Presentation(2)
        relation = poly([((1,), 1), ((), -1)])
        g = poly([((1, 1), 1), ((1,), -1)])
        basis = Basis([relation, g], 0, p.var_names())
        proof = Proof([[(F(1), (), "input", 0, ())], [(F(1), (1,), "basis", 0, ())]])
        steps = [((), 1, (2,), Coef(3, 0)), ((), 1, (2,), Coef(-3, 0))]
        skel = Skeleton.from_reduction(p, "X_2", (5, 7, 0), steps, Poly(0), basis, proof)
        self.assertEqual(skel.basis_keys, basis.key)
        self.assertTrue(skel.in_ideal)
        self.assertEqual(skel.step_keys(), [basis.key[1], basis.key[1]])
        self.assertEqual(compose_at(skel, basis, proof), [])
        reordered = Basis([g, relation], 0, p.var_names())
        self.assertEqual(skel.step_indices(reordered), [0, 0])
        self.assertEqual(skel.closure_indices(reordered), [1, 0])
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/"skeleton.json"
            skel.save(str(path))
            self.assertEqual(Skeleton.load(str(path)), skel)
        self.assertNotEqual(skel, object())
        skel.remainder_terms = 1
        self.assertFalse(skel.in_ideal)
        with self.assertRaisesRegex(ValueError, "proof has"):
            Skeleton.from_reduction(p, "X_2", (5, 7, 0), [], Poly(0), basis, Proof([]))


class TestRationalContextComposition(unittest.TestCase):
    def setUp(self):
        self.atoms = AtomSet([Atom("q", {(1, 0): 1})], max_exp=[1])
        self.pres = Presentation(2)

    def lift(self, entries):
        return Lift(self.atoms, entries, 0, 0, (5, 7, 0), {})

    def term(self, coefficient=1, kind="input", ref=0, left=(), right=()):
        return ((0,), {(0, 0): F(coefficient)}, left, kind, ref, right)

    def skel(self, steps=None):
        return Skeleton(self.pres.slug(), self.pres.content_hash(), "X_2", (5, 7, 0),
                        steps or [{"basis_key": "g", "left": (), "right": (), "mult": F(1)}])

    def test_rational_arithmetic_obeys_addition_and_multiplication(self):
        arithmetic = composition.Rat(self.atoms)
        a = ({(0, 0): F(1)}, (1,))
        b = arithmetic.const(F(2, 3))
        for q in (2, 5):
            self.assertEqual(arithmetic.eval(arithmetic.add(a, b), q, 7), F(1, q)+F(2, 3))
            self.assertEqual(arithmetic.eval(arithmetic.mul(a, b), q, 7), F(2, 3*q))
        self.assertEqual(arithmetic.const(0), arithmetic.zero())
        self.assertEqual(arithmetic.eval(arithmetic.zero(), 0, 7), 0)
        self.assertEqual(composition.poly_add({(0, 0): F(1)}, {(0, 0): F(-1)}), {})
        self.assertEqual(composition.poly_add({}, {(0, 0): F(0)}), {})
        self.assertEqual(composition.poly_mul({(1, 0): F(2)}, {(0, 1): F(3)}), {(1, 1): F(6)})

    def test_filtered_expansion_keeps_dependencies_and_drops_unrelated_entries(self):
        lift = self.lift({"outer": [self.term(2, "basis", "inner", (1,), ())],
                          "unused": [self.term(4)], "inner": [self.term(3)]})
        expanded = composition.expand_over_inputs(lift, keys=["outer"])
        self.assertEqual(set(expanded), {"outer", "inner"})
        self.assertEqual(expanded["outer"], {((1,), 0, ()): ({(0, 0): F(6)}, (0,))})

    def test_validation_precedes_zero_coefficient_cancellation(self):
        invalid = [
            (self.lift([]), ValueError, "mapping"),
            (self.lift({"g": None}), ValueError, "list"),
            (self.lift({"g": [()]}), ValueError, "malformed"),
            (self.lift({"g": [self.term(kind="basis", ref=[])]}), ValueError, "basis key"),
            (self.lift({"g": [self.term(left=(0,))]}), ValueError, "positive integer"),
        ]
        for lift, error, message in invalid:
            with self.subTest(message=message), self.assertRaisesRegex(error, message):
                composition.expand_over_inputs(lift)

    def test_skeleton_and_multiplier_boundaries_and_verbose_results(self):
        lift = self.lift({"g": [self.term()]})
        for steps, error in (([{}], ValueError),
                             ([{"basis_key": [], "left": (), "right": ()}], ValueError),
                             ([{"basis_key": "absent", "left": (), "right": ()}], KeyError)):
            with self.subTest(steps=steps), self.assertRaises(error):
                composition.compose(lift, self.skel(steps), self.pres,
                                    multipliers=[({(0, 0): F(1)}, (0,))], verbose=False)
        with self.assertRaisesRegex(ValueError, "multipliers="):
            composition.compose(lift, self.skel(), self.pres, verbose=False)
        with self.assertRaisesRegex(ValueError, "numerator"):
            composition.compose(lift, self.skel(), self.pres, multipliers=[None], verbose=False)
        stream = io.StringIO()
        with redirect_stdout(stream):
            certificate = composition.compose(lift, self.skel(), self.pres,
                       samples={(2, 3): [F(2)], (3, 5): [F(2)], (5, 7): [F(2)]})
        self.assertEqual(certificate.terms[0][3], {(0, 0): F(2)})
        self.assertIn("1 R-input terms", stream.getvalue())
        ring_target = self.pres.relations(NumericRing(2, 3, 0))[0]
        # The certificate is twice its input relation, so exactly that direct target matches.
        with redirect_stdout(stream):
            self.assertTrue(composition.numeric_check(certificate, self.pres, points=[(2, 3)],
                                                      target=ring_target.scalar_mul(Coef(2, 0))))
            self.assertFalse(composition.numeric_check(certificate, self.pres, points=[(2, 3)]))
        self.assertIn("MISMATCH", stream.getvalue())

    def test_insufficient_and_unmodellable_multiplier_samples(self):
        skel = self.skel()
        with self.assertRaises(NeedMoreSamples):
            composition.lift_multipliers(skel, {(2, 3): [F(1)]}, self.atoms)
        with self.assertRaises(NeedMoreSamples):
            composition.lift_multipliers(skel, {(2, 3): [F(1)], (3, 5): [F(2)], (5, 7): [F(7)]},
                                         AtomSet([], []), max_pq=0, max_pt=0)

    def test_zero_and_nonzero_anchor_multipliers(self):
        skel = self.skel([{"basis_key": "g", "left": (), "right": (), "mult": F(0)},
                          {"basis_key": "g", "left": (), "right": (), "mult": F(-2, 3)}])
        self.assertEqual(composition.multipliers_at_anchor(skel, self.atoms),
                         [({}, (0,)), ({(0, 0): F(-2, 3)}, (0,))])


class TestCompositionCommandLine(unittest.TestCase):
    def run_command(self, target, entries):
        pres = Presentation(3)
        atoms = AtomSet([], [])
        lift = Lift(atoms, entries, 0, 0, (5, 7, 0), {"g": 0})
        skel = Skeleton(pres.slug(), pres.content_hash(), target, (5, 7, 0),
                        [{"basis_key": "g", "left": (), "right": (), "mult": F(1)}])
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)
            lift.save(str(path/"lift.json")); skel.save(str(path/"skeleton.json"))
            args = [sys.executable, "-m", "zipper_algebra.cofactor.compose",
                    "--lift", str(path/"lift.json"), "--skeleton", str(path/"skeleton.json"),
                    "--out", str(path/"certificate.txt"), "--n", "3"]
            result = subprocess.run(args, cwd=Path(__file__).resolve().parents[1],
                                    text=True, capture_output=True, timeout=30)
            self.assertTrue((path/"certificate.txt").exists(), result.stderr)
            return result.returncode, result.stdout

    def test_cli_records_and_rejects_a_false_empty_expression(self):
        code, output = self.run_command("X_2", {"g": []})
        self.assertEqual(code, 1)
        self.assertIn("MISMATCH", output)

    def test_cli_verifies_a_contextual_untwisting_identity(self):
        # U*T3 + q²*(x1*(x2*y2-1)*y1 + (x1*y1-1))*S3
        # is (Theta3-q²)*S3 as a literal free-polynomial identity.
        pres = Presentation(3)
        terms = [((), {(0, 0): F(1)}, (1, 2), "input", pres.labels().index("T_3"), ())]
        slider = [((3, 4, 3), {(3, 0): F(1)}),
                  ((3, 4), {(2, 0): F(1), (3, 0): F(-1)}),
                  ((3, 4, 1), {(2, 0): F(-1)}),
                  ((1, 2, 3), {(1, 0): F(-1)}),
                  ((1, 2), {(0, 0): F(-1), (1, 0): F(1)}),
                  ((1, 2, 1), {(0, 0): F(1)})]
        for word, numerator in slider:
            shifted = {(q+2, t): c for (q, t), c in numerator.items()}
            terms.append(((), shifted, (1,), "input", pres.labels().index("INV_2_p"), (3,)+word))
            terms.append(((), shifted, (), "input", pres.labels().index("INV_1_p"), word))
        code, output = self.run_command("(Theta_3-q^2)X_3", {"g": terms})
        self.assertEqual(code, 0, output)
        self.assertNotIn("MISMATCH", output)


if __name__ == "__main__":
    unittest.main()
