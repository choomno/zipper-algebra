"""Exact basis parsing, language recognition and presentation-view regressions."""
from fractions import Fraction as F
from itertools import product
import contextlib
import io
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch
import warnings

from zipper_algebra import core
from zipper_algebra.presentation import Presentation
from zipper_algebra.coefficients import NumericRing, SymbolicRing
from zipper_algebra.polynomial import Poly


class TestBasisReadersAndReduction(unittest.TestCase):
    def parse(self, text, reader=core.parse_basis_Q):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "basis.txt"
            path.write_text(text)
            return reader(str(path))

    def test_ordered_terms_and_zero_polynomials(self):
        basis = self.parse("\n0,\n-0,\na-a,\n2*a + -3*b + 1,\n-a*b + 2/3*c - 5\n")
        self.assertEqual(basis, [
            (("a",), [(F(2), ("a",)), (F(-3), ("b",)), (F(1), ())]),
            (("a", "b"), [(F(-1), ("a", "b")), (F(2, 3), ("c",)), (F(-5), ())]),
        ])
        self.assertEqual(core.collect_variables(basis), ["a", "b", "c"])
        self.assertEqual(core.collect_variables_from_lms([("b", "a"), (), ("a", "c")]),
                         ["a", "b", "c"])
        self.assertEqual(core.collect_variables([]), [])
        self.assertEqual(core.collect_variables_from_lms([]), [])
        residues = self.parse("-x*x+4*x-1,\n3*y\n", core.parse_basis_Fp)
        self.assertEqual(residues[0], (("x", "x"), [(-1, ("x", "x")), (4, ("x",)), (-1, ())]))
        self.assertEqual(self.parse("x*x-1\ny\n", core.parse_leading_monomials),
                         [("x", "x"), ("y",)])
        self.assertEqual(self.parse("0\nx\n-3\ny\n", core.parse_leading_monomials), [()])

    def test_parser_rejects_missing_operators_names_and_coefficients(self):
        for text in (",\n", "a,,\n", "a,b\n", "a b\n", "a++b\n", "+\n", "a*\n", "1/0*a\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.parse(text)
        with self.assertRaisesRegex(ValueError, "integer residue"):
            self.parse("1/2*x\n", core.parse_basis_Fp)
        with self.assertRaisesRegex(ValueError, "invalid word"):
            core._parse_mon("x**y")
        with self.assertRaisesRegex(ValueError, "invalid polynomial term"):
            core._parse_term("@")
        self.assertEqual(core._parse_term("+3/2"), (F(3, 2), ()))
        self.assertEqual(core._parse_term("-x"), (F(-1), ("x",)))
        self.assertEqual(list(core._line_terms("")), [])

    def test_automaton_matches_contiguous_words_by_direct_search(self):
        alphabet = ["a", "b", "c"]
        patterns = [("a", "b", "a", "c"), ("b", "a", "b"), ("b", "a"),
                    ("c",), ("a", "b", "c"), ("a", "a", "b")]
        delta, matched, pattern, states = core.build_aho_corasick(patterns, alphabet)
        self.assertEqual(states, len(matched))
        for length in range(6):
            for word in product(alphabet, repeat=length):
                occurrences = [(start, start + len(lm), lm)
                               for start in range(len(word)) for lm in patterns
                               if word[start:start + len(lm)] == lm]
                actual = core.find_leading_match(word, delta, matched, pattern)
                if not occurrences:
                    self.assertIsNone(actual)
                else:
                    self.assertIn(actual, occurrences)
                    self.assertEqual(actual[1], min(end for _, end, _ in occurrences))
        self.assertIsNone(core.find_leading_match(("outside",), delta, matched, pattern))
        # A suffix that is itself forbidden must be recognized inside a longer trie prefix.
        d, m, p, _ = core.build_aho_corasick([("a", "b", "c"), ("b",)], alphabet)
        self.assertEqual(core.find_leading_match(("a", "b"), d, m, p), (1, 2, ("b",)))

    def test_standard_language_budget_and_unit_ideal(self):
        alphabet = ["a", "b"]
        all_pairs = list(product(alphabet, repeat=2))
        delta, matched, _, _ = core.build_aho_corasick(all_pairs, alphabet)
        self.assertEqual(core.enumerate_standard_monomials(delta, matched, alphabet),
                         [(), ("a",), ("b",)])
        self.assertIsNone(core.enumerate_standard_monomials(delta, matched, alphabet, max_total=2))
        d, m, p, _ = core.build_aho_corasick([()], alphabet)
        self.assertEqual(core.enumerate_standard_monomials(d, m, alphabet), [])
        self.assertEqual(core.find_leading_match(("a",), d, m, p), (0, 0, ()))
        d, m, _, _ = core.build_aho_corasick([], [])
        self.assertEqual(core.enumerate_standard_monomials(d, m, []), [()])

    def test_zero_budget_counts_the_unit_word_and_respects_the_unit_ideal(self):
        delta, matched, _, _ = core.build_aho_corasick([("x",)], ["x"])
        # The quotient has just the unit word: a zero-element budget is exceeded.
        self.assertEqual(core.enumerate_standard_monomials(delta, matched, ["x"]), [()])
        self.assertIsNone(core.enumerate_standard_monomials(delta, matched, ["x"], max_total=0))
        delta, matched, _, _ = core.build_aho_corasick([()], ["x"])
        self.assertEqual(core.enumerate_standard_monomials(delta, matched, ["x"], max_total=0), [])

    def test_exact_and_modular_reductions_preserve_polynomial_identities(self):
        # In Q<x>/(x^2-2x-1), x^3=5x+2 and x^4=12x+5.
        x = ("x",)
        polynomial = [(x + x, [(F(1), x + x), (F(-2), x), (F(-1), ())])]
        delta, matched, pattern, _ = core.build_aho_corasick([x + x], ["x"])
        basis = {(): 0, x: 1}
        tails = core.lm_to_tail_Q(polynomial)
        self.assertEqual(tails, {x + x: [(F(2), x), (F(1), ())]})
        self.assertEqual(core.reduce_to_basis_Q([(F(1), x * 4)], basis, tails,
                                                delta, matched, pattern), {1: F(12), 0: F(5)})
        self.assertEqual(core.reduce_to_basis_Q([(F(1, 2), x * 2), (F(-1), x), (F(-1, 2), ()),
                                                (F(0), x * 9)], basis, tails,
                                                delta, matched, pattern), {})
        modular_polynomial = [(lm, [(int(c), mon) for c, mon in terms]) for lm, terms in polynomial]
        modular_tails = core.lm_to_tail_Fp(modular_polynomial, 5)
        self.assertEqual(core.reduce_to_basis_Fp([(1, x * 3), (-2, ()), (0, x * 7)],
                                                 basis, modular_tails, delta, matched, pattern, 5), {})
        self.assertEqual(core.reduce_to_basis_Fp([(1, x * 4)], basis, modular_tails,
                                                 delta, matched, pattern, 5), {1: 2})
        with self.assertRaisesRegex(ValueError, "not invertible"):
            core.lm_to_tail_Fp([(x, [(2, x), (1, ())])], 4)
        # A nonmonic relation must divide by its leading coefficient.
        self.assertEqual(core.lm_to_tail_Q([(x, [(F(2), x), (F(-3), ())])]), {x: [(F(3, 2), ())]})
        self.assertEqual(core.lm_to_tail_Fp([(x, [(2, x), (-3, ())])], 5), {x: [(4, ())]})

    def test_noncommuting_contexts_and_incomplete_basis_diagnostics(self):
        # y(xy)=yx+2yy for the rule xy -> x+2y.
        xy = ("x", "y")
        delta, matched, pattern, _ = core.build_aho_corasick([xy], ["x", "y"])
        basis = {("y", "x"): 0, ("y", "y"): 1}
        tails = {xy: [(F(1), ("x",)), (F(2), ("y",))]}
        self.assertEqual(core.reduce_to_basis_Q([(F(3), ("y", "x", "y"))], basis,
                                                tails, delta, matched, pattern), {0: F(3), 1: F(6)})
        self.assertEqual(core.reduce_to_basis_Fp([(3, ("y", "x", "y"))], basis,
                                                 {xy: [(1, ("x",)), (2, ("y",))]}, delta, matched, pattern, 5), {0: 3, 1: 1})
        for reducer, args in ((core.reduce_to_basis_Q, ()), (core.reduce_to_basis_Fp, (5,))):
            with self.subTest(reducer=reducer.__name__), self.assertRaisesRegex(RuntimeError, "incomplete"):
                reducer([(1, ("z",))], {}, tails, delta, matched, pattern, *args)


class TestPresentationViews(unittest.TestCase):
    def test_two_strand_symbolic_relations_have_the_declared_parameters(self):
        pres = Presentation(2)
        self.assertEqual(pres.var_names(), ["x1", "y1"])
        self.assertEqual(pres.slug(), "n2-bigelow")
        expected_t2 = {
            (1, 2): {(1, 0): F(1)},
            (1,): {(0, 0): F(1), (1, 0): F(-1), (0, 1): F(1)},
            (1, 1): {(0, 0): F(-1)},
            (2,): {(1, 1): F(-1)},
            (): {(1, 1): F(1), (0, 1): F(-1)},
        }
        fitted = pres.symbolic_fitted()
        self.assertEqual(fitted[2], expected_t2)
        relations, sliders = pres.symbolic_ncpolys()
        self.assertEqual(relations[2], expected_t2)
        self.assertEqual(sliders[(2, 1)], {(2,): {(1, 0): F(1)},
                                           (): {(0, 0): F(1), (1, 0): F(-1)},
                                           (1,): {(0, 0): F(-1)}})
        self.assertEqual(pres.labels(), ["INV_1_p", "INV_1_m", "T_2"])
        self.assertEqual(pres.n_relations(), 3)
        ring = NumericRing(2, 3, 0)
        self.assertEqual(len(pres.relations(ring)), 3)
        self.assertEqual({word: coef.v for word, coef in pres.target(2, 1, ring).terms.items()},
                         {(2,): F(2), (): F(-1), (1,): F(-1)})

    def test_zero_sliders_and_extra_relations_retain_order_and_labels(self):
        extras = ("X_4", "X_5", "(Theta_3-q^2)X_3", "(Theta_4+t)X_4", "(Theta_5+t)X_5")
        pres = Presentation(5, untwist="none", zero_from=4, extra=extras)
        ring = NumericRing(2, 3, 0)
        names, relations, sliders, tags = pres.build(ring)
        self.assertEqual(names, ["x1", "x2", "x3", "x4", "y1", "y2", "y3", "y4"])
        self.assertEqual([tag.kind for tag in tags[-8:]], ["ZERO"] * 3 + ["EXTRA"] * 5)
        self.assertEqual([rel.terms for rel in relations[-8:-5]],
                         [sliders[(4, 1)].terms, sliders[(4, 2)].terms, sliders[(5, 1)].terms])
        self.assertEqual(relations[-5].terms, sliders[(4, 1)].terms)
        self.assertEqual(relations[-4].terms, sliders[(5, 1)].terms)
        for offset, level, constant in ((-3, 3, -4), (-2, 4, 3), (-1, 5, 3)):
            twist = Poly(0)
            word = tuple(range(1, level)) + tuple(range(level - 1, 0, -1))
            twist.add_term(word, ring.one)
            twist.add_term((), ring.scalar(constant))
            self.assertEqual(relations[offset].terms, (twist * sliders[(level, 1)]).terms)
        self.assertEqual(pres.labels()[-5:], [f"EXTRA[{name}]" for name in extras])
        self.assertIn("ZERO_4_1", pres.labels())
        text = pres.canonical_text()
        self.assertIn("(Theta_5+t)X_5 = ", text)
        self.assertEqual(text, pres.canonical_text())
        self.assertNotEqual(pres.content_hash(), Presentation(5, untwist="none").content_hash())

    def test_slug_records_every_presentation_choice(self):
        pres = Presentation(5, "none", framing_power=2, theta_power=3,
                            max_zipper_level=2, max_untwist_level=3, zero_from=4,
                            extra=("(Theta_5+t)X_5",))
        self.assertEqual(pres.slug(), "n5-none-a2-k3-Z2-L3-z4-x_Theta_5_t_X_5")

    def test_emitted_input_encodes_exact_two_strand_relations(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.ms"
            self.assertEqual(Presentation(2).emit_ms(str(path), 2, 3, 0), 3)
            lines = path.read_text().splitlines()
        self.assertEqual(lines[:2], ["x1,y1", "0"])
        self.assertTrue(all(line.endswith(",") for line in lines[2:-1]))
        terms = [core._parse_term(term) for term in core._line_terms(lines[-1])]
        self.assertEqual({word: coef for coef, word in terms},
                         {("x1", "x1"): F(-1), ("x1", "y1"): F(2),
                          ("x1",): F(2), ("y1",): F(-6), (): F(3)})

    def test_module_cli_reports_requested_presentation_and_labels(self):
        output = io.StringIO()
        with patch("sys.argv", ["presentation", "--n", "3", "--untwist", "none",
                                 "--max-untwist-level", "2", "--max-zipper-level", "2",
                                 "--extra", "(Theta_3-q^2)X_3"]), contextlib.redirect_stdout(output), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            runpy.run_module("zipper_algebra.presentation", run_name="__main__")
        text = output.getvalue()
        self.assertIn("n3-none-Z2-L2-x_Theta_3_q_2_X_3", text)
        self.assertIn("EXTRA[(Theta_3-q^2)X_3]", text)
        self.assertIn("relations", text)
