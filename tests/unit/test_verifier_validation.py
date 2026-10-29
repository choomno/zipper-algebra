"""Adversarial certificate checks independent of solver binaries and cached samples."""
from __future__ import annotations

from contextlib import redirect_stdout
from io import StringIO
import subprocess
import sys
import tempfile
import pytest
from fractions import Fraction
from pathlib import Path

from zipper_algebra.cofactor.atoms import Atom, AtomSet
from zipper_algebra.cofactor.cofactor_io import Cofactor, symbolic_target
from zipper_algebra.presentation import Presentation
from zipper_algebra.cofactor.verify import check, pit_check
from zipper_algebra.coefficients import Coef, SymbolicRing
from zipper_algebra.polynomial import Poly


class TestVerifierValidation:
    @pytest.fixture(autouse=True)
    def case_groups(self, subtests):
        self.subtests = subtests

    @pytest.fixture(autouse=True)
    def setup_algebra(self):
        self.pres = Presentation(4)

    def certificate(self, terms=(), atoms=None):
        return Cofactor.from_presentation(self.pres, "zero", atoms or AtomSet([]), terms,
                                          "tests_verifier_validation")

    def reject_both(self, cof, target=None, pres=None):
        for verifier in (check, pit_check):
            with self.subtests.test(verifier=verifier.__name__):
                with pytest.raises((ValueError, TypeError)):
                    verifier(cof, pres or self.pres, {} if target is None else target,
                             verbose=False)

    def test_zero_denominator_is_rejected_at_construction(self):
        for coeffs in ({}, {(0, 0): 0}, {(0, 0): Fraction(1, 2)}):
            with self.subtests.test(coeffs=coeffs):
                with pytest.raises(ValueError):
                    Atom("bad", coeffs)

    def test_mutated_zero_atom_cannot_clear_a_false_identity(self):
        atoms = AtomSet([Atom("a", {(1, 0): 1})])
        cof = self.certificate([((), self.pres.labels()[0], (), {}, (1,))], atoms)
        atoms.atoms[0].coefs.clear()
        self.reject_both(cof, {(): {(0, 0): 1}})

    def test_atom_and_box_schema(self):
        for coeffs in ({(-1, 0): 1}, {(0.0, 0): 1}, {(0, 0): 1.2}):
            with self.subtests.test(coeffs=coeffs):
                with pytest.raises(ValueError):
                    Atom("a", coeffs)
        with pytest.raises(ValueError):
            AtomSet([Atom("q", {(1, 0): 1})], [])
        with pytest.raises(ValueError):
            AtomSet([Atom("q", {(1, 0): 1}), Atom("q", {(0, 1): 1})])

    def test_invalid_terms_are_rejected_even_with_zero_numerator(self):
        base = ((), self.pres.labels()[0], (), {}, ())
        bad_terms = [((0,), *base[1:]), ((7,), *base[1:]), ((True,), *base[1:]),
                     ((), "unknown", (), {}, ()), ((), base[1], (), {}, (1,)),
                     ((), base[1], (), {(-1, 0): 1}, ()),
                     ((), base[1], (), {(0, 0): 0.5}, ())]
        for term in bad_terms:
            with self.subtests.test(term=term):
                cof = self.certificate()
                cof.terms.append(term)
                self.reject_both(cof)

    def test_denominator_exponents_have_exact_arity_and_domain(self):
        atoms = AtomSet([Atom("q", {(1, 0): 1})])
        for den in ((), (0, 0), (-1,), (0.5,), (True,)):
            with self.subtests.test(den=den):
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
            with self.subtests.test(verifier=verifier.__name__):
                with pytest.raises(ValueError, match='explicit pres'):
                    verifier(cof, target={}, verbose=False)
                assert verifier(cof, self.pres, target={}, verbose=False)['pass']

    def test_valid_zero_and_nonzero_identities(self):
        for verifier in (check, pit_check):
            assert verifier(self.certificate(), self.pres, {}, verbose=False)['pass']
        atoms = AtomSet([Atom("q", {(1, 0): 1})], max_exp=[0])
        cof = self.certificate([((), self.pres.labels()[0], (), {(2, 0): 1}, (2,))], atoms)
        target = self.pres.symbolic_ncpolys()[0][0]
        for verifier in (check, pit_check):
            # Search exponent boxes are not certificate validity restrictions.
            assert verifier(cof, self.pres, target, verbose=False)['pass']

    def test_pit_evaluates_cleared_identity_at_poles(self):
        atoms = AtomSet([Atom("t-q", {(0, 1): 1, (1, 0): -1})])
        rel = self.pres.labels()[0]
        cof = self.certificate([((), rel, (), {(0, 0): 1}, (1,)),
                                ((), rel, (), {(0, 0): -1}, (1,))], atoms)
        assert check(cof, self.pres, {}, verbose=False)['pass']
        assert pit_check(cof, self.pres, {}, verbose=False)['pass']

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
            with self.subtests.test(verifier=verifier.__name__):
                with redirect_stdout(StringIO()) as output:
                    report = verifier(cof, self.pres, target)
                    assert report['pass']
                    assert report['lcm'] == '(q-1)'
                    assert report['lcm_exps'] == [1]
                    assert 'M = (q-1)' in output.getvalue()
                    assert 't+5' not in output.getvalue()
        with pytest.raises(ZeroDivisionError):
            cof.evaluate(1, 2, self.pres)
        assert atoms.denom_str((1,)) == '(t+5)'
        assert cof.header['atoms']['atoms'][0]['name'] == 't+5'
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "misleading-name.txt"
            cof.save(path)
            restored = Cofactor.load(path, pres=self.pres)
            assert restored.header == cof.header
            assert restored.terms == cof.terms
            assert '(t+5)' in Cofactor.strip_header(path.read_text())
            described = subprocess.run([sys.executable, "-m", "zipper_algebra.cofactor.cofactor_io",
                                        str(path)], check=True, capture_output=True, text=True)
            assert 'LCM denominator: (q-1)' in described.stdout

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
        assert not check(self.certificate(), self.pres, target, verbose=False)['pass']
        report = pit_check(self.certificate(), self.pres, target, verbose=False)
        assert not report['pass']
        assert report['bidegree'][0] >= 10

    @pytest.mark.parametrize('degree_source', ['target', 'numerator', 'denominator'])
    def test_t_degree_bound_cannot_accept_a_false_identity(self, degree_source):
        # f(t)=(t-2)(t-3)(t-4)(t-5).  An undersized grid containing only t=2
        # would mistake this nonzero polynomial for zero.  These coefficients
        # are stated independently of the polynomial multiplication routines.
        polynomial = {(0, 4): 1, (0, 3): -14, (0, 2): 71,
                      (0, 1): -154, (0, 0): 120}
        pres = Presentation(2, 'none')
        label = pres.labels()[0]  # x1*y1-1, with constant coefficients.
        atoms, terms, target = AtomSet([]), [], {}
        if degree_source == 'target':
            target = {(): polynomial}
        elif degree_source == 'numerator':
            terms = [((), label, (), polynomial, ())]
        else:
            atoms = AtomSet([Atom('f', polynomial)])
            # The zero term still declares a denominator. Clearing it makes
            # the first term f(t)*(x1*y1-1), whose t-degree must be counted.
            terms = [((), label, (), {(0, 0): 1}, (0,)),
                     ((), label, (), {}, (1,))]
        cof = Cofactor.from_presentation(pres, 'zero', atoms, terms,
                                         'tests_verifier_validation')
        assert not check(cof, pres, target, verbose=False)['pass']
        report = pit_check(cof, pres, target, verbose=False)
        assert not report['pass']
        assert report['bidegree'][1] >= 4
        assert report['grid'][1] >= 5

    @pytest.mark.parametrize('perturbation', [0, 1])
    def test_mixed_denominators_require_their_clearing_complements(self, perturbation):
        pres = Presentation(2, 'none')
        label = pres.labels()[0]  # R = x1*y1-1.
        atoms = AtomSet([Atom('A', {(1, 0): 1, (0, 0): -2}),
                         Atom('B', {(0, 1): 1, (0, 0): -3})])
        # R/A + R/B - (A+B)*R/(A*B) = 0.  The literal last numerator
        # is -q-t+5; adding one gives the nonzero residual R/(A*B).
        terms = [((), label, (), {(0, 0): 1}, (1, 0)),
                 ((), label, (), {(0, 0): 1}, (0, 1)),
                 ((), label, (), {(1, 0): -1, (0, 1): -1,
                                   (0, 0): 5 + perturbation}, (1, 1))]
        cof = Cofactor.from_presentation(pres, 'zero', atoms, terms,
                                         'tests_verifier_validation')
        for verifier in (check, pit_check):
            report = verifier(cof, pres, {}, verbose=False)
            assert report['pass'] is (perturbation == 0)
            assert report['lcm_exps'] == [1, 1]

    def test_extra_and_poly_targets_use_the_same_resolution(self):
        pres = Presentation(4, extra=("X_4",))
        cof = Cofactor.from_presentation(pres, "EXTRA[X_4]", AtomSet([]),
                 [((), "EXTRA[X_4]", (), {(0, 0): 1}, ())], "tests_verifier_validation")
        poly_target = pres.target(4, 1, SymbolicRing())
        assert symbolic_target(pres, 'EXTRA[X_4]') == symbolic_target(pres, poly_target)
        for verifier in (check, pit_check):
            for target in ("EXTRA[X_4]", poly_target):
                assert verifier(cof, pres, target, verbose=False)['pass']

    def test_malformed_targets_are_rejected(self):
        for target in ({(7,): {(0, 0): 1}}, {(): {(-1, 0): 1}},
                       {(): {(0, 0): 0.1}}, (2.0, 1), "S[3][1]trailing"):
            with self.subtests.test(target=target):
                self.reject_both(self.certificate(), target)
        finite = Poly.scalar(Coef.one(5))
        self.reject_both(self.certificate(), finite)

    def test_atom_names_containing_products_round_trip(self):
        atoms = AtomSet([Atom("2*t+1", {(0, 1): 2, (0, 0): 1})])
        assert atoms.parse_denom_str(atoms.denom_str((2,))) == (2,)
        monomial = AtomSet([Atom("q*t", {(1, 1): 1})])
        assert monomial.parse_denom_str(monomial.denom_str((2,))) == (2,)
        with pytest.raises(ValueError):
            atoms.parse_denom_str("(2*t+1)^-1")

    def test_legacy_extraction_rejects_invalid_relation_indices(self):
        for labels in (None, self.pres.labels()):
            with self.subtests.test(labels=labels):
                with pytest.raises(ValueError, match='relation index'):
                    Cofactor._parse_extract_line("+1 * () * [-1] * ()", {}, AtomSet([]), labels)
        with pytest.raises(ValueError, match='relation index'):
            Cofactor._parse_extract_line("+1 * () * [999] * ()", {}, AtomSet([]),
                                         self.pres.labels())

    def test_current_example_certificate_body_round_trips(self):
        source = Path(__file__).resolve().parents[2] / "examples/x3/certificate.txt"
        cof = Cofactor.load(str(source), pres=self.pres)
        assert check(cof, self.pres, 'X_3', verbose=False)['pass']
        with tempfile.TemporaryDirectory() as directory:
            dest = Path(directory) / "cofactor.txt"
            cof.save(str(dest))
            assert Cofactor.strip_header(source.read_text()) == Cofactor.strip_header(dest.read_text())
