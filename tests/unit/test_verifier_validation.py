"""Adversarial certificate checks independent of solver binaries and cached samples."""
from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.cofactor_io import Cofactor, symbolic_target
from zipper_algebra.presentation import Presentation
from zipper_algebra.cofactor.verify import check, pit_check
from zipper_algebra.coefficients import Coef, SymbolicRing
from zipper_algebra.polynomial import Poly


class TestVerifierValidation(unittest.TestCase):
    def setUp(self):
        self.pres = Presentation(4)

    def certificate(self, terms=(), atoms=None):
        return Cofactor.from_presentation(self.pres, "zero", atoms or AtomSet([]), terms,
                                          "tests_verifier_validation")

    def reject_both(self, cof, target=None, pres=None):
        for verifier in (check, pit_check):
            with self.subTest(verifier=verifier.__name__):
                with self.assertRaises((ValueError, TypeError)):
                    verifier(cof, pres or self.pres, {} if target is None else target,
                             verbose=False)

    def test_zero_denominator_is_rejected_at_construction(self):
        for coeffs in ({}, {(0, 0): 0}, {(0, 0): Fraction(1, 2)}):
            with self.subTest(coeffs=coeffs), self.assertRaises(ValueError):
                Atom("bad", coeffs)

    def test_mutated_zero_atom_cannot_clear_a_false_identity(self):
        atoms = AtomSet([Atom("a", {(1, 0): 1})])
        cof = self.certificate([((), self.pres.labels()[0], (), {}, (1,))], atoms)
        atoms.atoms[0].coefs.clear()
        self.reject_both(cof, {(): {(0, 0): 1}})

    def test_atom_and_box_schema(self):
        for coeffs in ({(-1, 0): 1}, {(0.0, 0): 1}, {(0, 0): 1.2}):
            with self.subTest(coeffs=coeffs), self.assertRaises(ValueError):
                Atom("a", coeffs)
        with self.assertRaises(ValueError):
            AtomSet([Atom("q", {(1, 0): 1})], [])
        with self.assertRaises(ValueError):
            AtomSet([Atom("q", {(1, 0): 1}), Atom("q", {(0, 1): 1})])

    def test_invalid_terms_are_rejected_even_with_zero_numerator(self):
        base = ((), self.pres.labels()[0], (), {}, ())
        bad_terms = [((0,), *base[1:]), ((7,), *base[1:]), ((True,), *base[1:]),
                     ((), "unknown", (), {}, ()), ((), base[1], (), {}, (1,)),
                     ((), base[1], (), {(-1, 0): 1}, ()),
                     ((), base[1], (), {(0, 0): 0.5}, ())]
        for term in bad_terms:
            with self.subTest(term=term):
                cof = self.certificate()
                cof.terms.append(term)
                self.reject_both(cof)

    def test_denominator_exponents_have_exact_arity_and_domain(self):
        atoms = AtomSet([Atom("q", {(1, 0): 1})])
        for den in ((), (0, 0), (-1,), (0.5,), (True,)):
            with self.subTest(den=den):
                cof = self.certificate(atoms=atoms)
                cof.terms.append(((), self.pres.labels()[0], (), {}, den))
                self.reject_both(cof)

    def test_hash_and_strand_count_checked_by_both_verifiers(self):
        cof = self.certificate()
        self.reject_both(cof, pres=Presentation(4, "square-zero"))
        cof = self.certificate()
        cof.n = 3
        self.reject_both(cof)

    def test_untagged_certificate_requires_explicit_presentation(self):
        cof = self.certificate()
        cof.presentation_hash = ""
        for verifier in (check, pit_check):
            with self.subTest(verifier=verifier.__name__):
                with self.assertRaisesRegex(ValueError, "explicit pres"):
                    verifier(cof, target={}, verbose=False)
                self.assertTrue(verifier(cof, self.pres, target={}, verbose=False)["pass"])

    def test_valid_zero_and_nonzero_identities(self):
        for verifier in (check, pit_check):
            self.assertTrue(verifier(self.certificate(), self.pres, {}, verbose=False)["pass"])
        atoms = AtomSet([Atom("q", {(1, 0): 1})], max_exp=[0])
        cof = self.certificate([((), self.pres.labels()[0], (), {(2, 0): 1}, (2,))], atoms)
        target = self.pres.symbolic_ncpolys()[0][0]
        for verifier in (check, pit_check):
            # Search exponent boxes are not certificate validity restrictions.
            self.assertTrue(verifier(cof, self.pres, target, verbose=False)["pass"])

    def test_pit_evaluates_cleared_identity_at_poles(self):
        atoms = AtomSet([Atom("t-q", {(0, 1): 1, (1, 0): -1})])
        rel = self.pres.labels()[0]
        cof = self.certificate([((), rel, (), {(0, 0): 1}, (1,)),
                                ((), rel, (), {(0, 0): -1}, (1,))], atoms)
        self.assertTrue(check(cof, self.pres, {}, verbose=False)["pass"])
        self.assertTrue(pit_check(cof, self.pres, {}, verbose=False)["pass"])

    def test_verification_displays_actual_polynomial_for_misleading_atom_name(self):
        polynomial = {(1, 0): 1, (0, 0): -1}
        atoms = AtomSet([Atom("t+5", polynomial)])
        # The two coefficients 1/(q-1) and (q-2)/(q-1) each have a pole at q=1;
        # their sum is 1, so the nonzero relation is an exact target.
        rel = self.pres.labels()[0]
        cof = self.certificate([((), rel, (), {(0, 0): 1}, (1,)),
                                ((), rel, (), {(1, 0): 1, (0, 0): -2}, (1,))], atoms)
        target = self.pres.symbolic_ncpolys()[0][0]
        for verifier in (check, pit_check):
            with self.subTest(verifier=verifier.__name__), redirect_stdout(StringIO()) as output:
                report = verifier(cof, self.pres, target)
                self.assertTrue(report["pass"])
                self.assertEqual(report["lcm"], "(q-1)")
                self.assertEqual(report["lcm_exps"], [1])
                self.assertIn("M = (q-1)", output.getvalue())
                self.assertNotIn("t+5", output.getvalue())
        with self.assertRaises(ZeroDivisionError):
            cof.evaluate(1, 2, self.pres)
        self.assertEqual(atoms.denom_str((1,)), "(t+5)")
        self.assertEqual(cof.header["atoms"]["atoms"][0]["name"], "t+5")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "misleading-name.txt"
            cof.save(path)
            restored = Cofactor.load(path, pres=self.pres)
            self.assertEqual(restored.header, cof.header)
            self.assertEqual(restored.terms, cof.terms)
            self.assertIn("(t+5)", Cofactor.strip_header(path.read_text()))
            described = subprocess.run([sys.executable, "-m", "zipper_algebra.cofactor.cofactor_io",
                                        str(path)], check=True, capture_output=True, text=True)
            self.assertIn("LCM denominator: (q-1)", described.stdout)

    def test_laurent_grid_forgery_is_rejected(self):
        # The old PIT bound ignored negative powers. This nonzero Laurent polynomial
        # vanishes at every q in the old n=4 grid, giving a false positive.
        num = {(0, 0): Fraction(1)}
        for k in range(2, 12):
            out = {}
            for (i, j), c in num.items():
                out[(i - 1, j)] = out.get((i - 1, j), 0) + c
                out[(i, j)] = out.get((i, j), 0) - c / k
            num = out
        cof = self.certificate()
        cof.terms.append(((), self.pres.labels()[0], (), num, ()))
        self.reject_both(cof)

    def test_target_degree_is_included_in_pit_bound(self):
        num = {(0, 0): Fraction(1)}
        for k in range(2, 12):
            out = {}
            for (i, j), c in num.items():
                out[(i + 1, j)] = out.get((i + 1, j), 0) + c
                out[(i, j)] = out.get((i, j), 0) - k * c
            num = out
        target = {(): num}
        self.assertFalse(check(self.certificate(), self.pres, target, verbose=False)["pass"])
        report = pit_check(self.certificate(), self.pres, target, verbose=False)
        self.assertFalse(report["pass"])
        self.assertGreaterEqual(report["bidegree"][0], 10)

    def test_extra_and_poly_targets_use_the_same_resolution(self):
        pres = Presentation(4, extra=("X_4",))
        cof = Cofactor.from_presentation(pres, "EXTRA[X_4]", AtomSet([]),
                 [((), "EXTRA[X_4]", (), {(0, 0): 1}, ())], "tests_verifier_validation")
        poly_target = pres.target(4, 1, SymbolicRing())
        self.assertEqual(symbolic_target(pres, "EXTRA[X_4]"), symbolic_target(pres, poly_target))
        for verifier in (check, pit_check):
            for target in ("EXTRA[X_4]", poly_target):
                self.assertTrue(verifier(cof, pres, target, verbose=False)["pass"])

    def test_malformed_targets_are_rejected(self):
        for target in ({(7,): {(0, 0): 1}}, {(): {(-1, 0): 1}},
                       {(): {(0, 0): 0.1}}, (2.0, 1), "S[3][1]trailing"):
            with self.subTest(target=target):
                self.reject_both(self.certificate(), target)
        finite = Poly.scalar(Coef.one(5))
        self.reject_both(self.certificate(), finite)

    def test_atom_names_containing_products_round_trip(self):
        atoms = AtomSet([Atom("2*t+1", {(0, 1): 2, (0, 0): 1})])
        self.assertEqual(atoms.parse_denom_str(atoms.denom_str((2,))), (2,))
        monomial = AtomSet([Atom("q*t", {(1, 1): 1})])
        self.assertEqual(monomial.parse_denom_str(monomial.denom_str((2,))), (2,))
        with self.assertRaises(ValueError):
            atoms.parse_denom_str("(2*t+1)^-1")

    def test_legacy_extraction_rejects_invalid_relation_indices(self):
        for labels in (None, self.pres.labels()):
            with self.subTest(labels=labels), self.assertRaisesRegex(ValueError, "relation index"):
                Cofactor._parse_extract_line("+1 * () * [-1] * ()", {}, AtomSet([]), labels)
        with self.assertRaisesRegex(ValueError, "relation index"):
            Cofactor._parse_extract_line("+1 * () * [999] * ()", {}, AtomSet([]),
                                         self.pres.labels())

    def test_current_example_certificate_body_round_trips(self):
        source = Path(__file__).resolve().parents[2] / "examples/x3/certificate.txt"
        cof = Cofactor.load(str(source), pres=self.pres)
        self.assertTrue(check(cof, self.pres, "X_3", verbose=False)["pass"])
        with tempfile.TemporaryDirectory() as directory:
            dest = Path(directory) / "cofactor.txt"
            cof.save(str(dest))
            self.assertEqual(Cofactor.strip_header(source.read_text()),
                             Cofactor.strip_header(dest.read_text()))


if __name__ == "__main__":
    unittest.main(verbosity=2)
