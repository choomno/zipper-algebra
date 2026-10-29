#!/usr/bin/env python3
"""Equivalence tests for cofactor_io / compose / verify against the shipped X_3 artifacts.

The claim under test is not "the new modules run" but "they produce the same bytes and the same
verdict as the scripts they replace".  Everything here reads `tests/fixtures/legacy_x3/`; nothing
writes there.
"""
from __future__ import annotations

import json
from copy import deepcopy
import os
from pathlib import Path
import tempfile
import unittest
from fractions import Fraction

from tests.support import LEGACY_X3

from zipper_algebra.cofactor.atoms import AtomSet
from zipper_algebra.cofactor.cofactor_io import Cofactor
from zipper_algebra.cofactor.compose import compose, expand_over_inputs, numeric_check
from zipper_algebra.cofactor.lift import Lift
from zipper_algebra.presentation import Presentation
from zipper_algebra.cofactor import verify


# Immutable pre-consolidation inputs are carried by this repository.
FROZEN = str(LEGACY_X3)


def frozen(name: str) -> str:
    return os.path.join(FROZEN, name)


SYMBOLIC = os.path.join(FROZEN, "cofactor-X3-n4-Qqt-symbolic.txt")
LIFT_JSON = frozen("basis-symbolic-Qqt.json")
SKEL_JSON = frozen("skeleton-n4-Q.json")
_OUTPUT = None
OUT = None


def setUpModule():
    global _OUTPUT, OUT
    _OUTPUT = tempfile.TemporaryDirectory(prefix="zipper-legacy-tests-")
    OUT = _OUTPUT.name


def tearDownModule():
    if _OUTPUT is not None:
        _OUTPUT.cleanup()


PRES = Presentation(4)
TARGET = "S[3][1]"


# --- adapters from the legacy on-disk formats -------------------------------

def legacy_lift() -> Lift:
    """basis-symbolic-Qqt.json (integer-string keys, 'idx', words as variable names) -> Lift."""
    d = json.loads(Path(LIFT_JSON).read_text())
    # the frozen snapshot predates the atom feature and carries no "atoms" key; the six
    # historical atoms are exactly the class the legacy script enumerated.
    atoms = AtomSet.from_json(d["atoms"]) if "atoms" in d else AtomSet.DEFAULT
    names = PRES.var_names()
    index = {v: i + 1 for i, v in enumerate(names)}
    entries, anchor_index = {}, {}
    for k, terms in d["lifted"].items():
        out = []
        for t in terms:
            if t is None:
                continue
            out.append((tuple(t["denom"]),
                        {tuple(m): Fraction(c) for m, c in t["num"]},
                        tuple(index[w] for w in t["left"]),
                        t["kind"],
                        str(t["idx"]) if t["kind"] == "basis" else t["idx"],
                        tuple(index[w] for w in t["right"])))
        entries[k] = out
        anchor_index[k] = int(k)
    return Lift(atoms, entries, samples_used=d.get("samples_used", 0), holdout=0,
                anchor=(5, 7, 0), anchor_index=anchor_index)


class _Skel:
    """skeleton-n4-Q.json as compose expects it (steps carry a basis_key)."""

    def __init__(self, path):
        d = json.loads(Path(path).read_text())
        self.target = d["target"]
        a = d["anchor"]
        self.anchor = (a["q"], a["t"], a["char"])
        self.steps = [{"basis_key": str(s["basis_index"]), "left": tuple(s["left"]),
                       "right": tuple(s["right"]), "mult": Fraction(s["mult"])}
                      for s in d["steps"]]
        self.remainder_terms = d["remainder_terms"]


def legacy_multipliers(atoms: AtomSet):
    """q^3, -q^2, q(q-1) -- the multipliers the legacy compose typed in, for the byte test."""
    z = atoms.zero()
    return [({(3, 0): Fraction(1)}, z),
            ({(2, 0): Fraction(-1)}, z),
            ({(2, 0): Fraction(1), (1, 0): Fraction(-1)}, z)]


# --- cofactor_io ------------------------------------------------------------

class TestCofactorIO(unittest.TestCase):
    def test_load_save_round_trip_is_byte_identical_in_the_body(self):
        cof = Cofactor.load(SYMBOLIC, pres=PRES)
        self.assertEqual(len(cof.terms), 695)
        path = os.path.join(OUT, "roundtrip.txt")
        cof.save(path)
        original = Cofactor.strip_header(Path(SYMBOLIC).read_text())
        written = Cofactor.strip_header(Path(path).read_text())
        self.assertEqual(written, original)

    def test_header_round_trips_through_a_reload(self):
        cof = Cofactor.load(SYMBOLIC, pres=PRES)
        cof.target, cof.presentation = TARGET, PRES.slug()
        cof.presentation_hash = PRES.content_hash()
        path = os.path.join(OUT, "withheader.txt")
        cof.save(path)
        back = Cofactor.load(path)
        self.assertEqual(back.target, TARGET)
        self.assertEqual(back.presentation_hash, PRES.content_hash())
        self.assertEqual(back.atoms.names(), cof.atoms.names())
        self.assertEqual(Cofactor.strip_header(Path(path).read_text()),
                         Cofactor.strip_header(Path(SYMBOLIC).read_text()))

    def test_denominators_parse_through_the_atom_set(self):
        cof = Cofactor.load(SYMBOLIC, pres=PRES)
        for _l, _lab, _r, _num, den in cof.terms:
            self.assertEqual(len(den), len(cof.atoms))
            self.assertEqual(cof.atoms.parse_denom_str(cof.atoms.denom_str(den)), tuple(den))


# --- compose ----------------------------------------------------------------

class TestCompose(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lift = legacy_lift()
        cls.skel = _Skel(SKEL_JSON)
        cls.expanded = expand_over_inputs(cls.lift)
        cls.certificate = compose(cls.lift, cls.skel, PRES,
                                  multipliers=legacy_multipliers(cls.lift.atoms), verbose=False)

    def test_expansion_sizes_match_the_legacy_transcript(self):
        exp = self.expanded
        # the driver's phase 4 log records basis[68] -> 378 and basis[102] -> 695
        self.assertEqual(len(exp["102"]), 695)
        self.assertEqual(len(exp["68"]), 378)
        self.assertEqual(set(exp), set(self.lift.entries))

    def test_composed_body_is_byte_identical_to_the_shipped_certificate(self):
        cof = deepcopy(self.certificate)
        self.assertEqual(len(cof.terms), 695)
        path = os.path.join(OUT, "composed.txt")
        cof.save(path)
        self.assertEqual(Cofactor.strip_header(Path(path).read_text()),
                         Cofactor.strip_header(Path(SYMBOLIC).read_text()))

    def test_composed_cofactor_checks_out_numerically(self):
        cof = deepcopy(self.certificate)
        self.assertTrue(numeric_check(cof, PRES, target=TARGET, verbose=False))


# --- verify -----------------------------------------------------------------

class TestVerify(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cof = Cofactor.load(SYMBOLIC, pres=PRES)

    def test_check_passes_with_the_transcript_numbers(self):
        rep = verify.check(self.cof, PRES, TARGET, verbose=False)
        self.assertTrue(rep["pass"])
        self.assertEqual(rep["tuples"], 695)
        self.assertEqual(rep["relations_used"],
                         [0, 1, 2, 3, 6, 7, 14, 15, 16, 17, 22, 25, 26, 28])
        self.assertEqual(rep["lcm"], "q^5*t^4*(q-1)*(t+q^2)")
        with open(os.path.join(OUT, "verify-report.json"), "w") as f:
            json.dump(rep, f, indent=1)

    def test_a_tampered_coefficient_fails(self):
        bad = Cofactor.load(SYMBOLIC, pres=PRES)
        left, lab, right, num, den = bad.terms[0]
        tampered = dict(num)
        k = next(iter(tampered))
        tampered[k] = tampered[k] + 1
        bad.terms[0] = (left, lab, right, tampered, den)
        self.assertFalse(verify.check(bad, PRES, TARGET, verbose=False)["pass"])

    def test_a_cofactor_for_another_algebra_is_refused(self):
        """Labels coincide between bigelow and square-zero, so only the hash catches this."""
        tagged = Cofactor.load(SYMBOLIC, pres=PRES)
        tagged.presentation, tagged.presentation_hash = PRES.slug(), PRES.content_hash()
        other = Presentation(4, "square-zero")
        self.assertEqual(PRES.labels(), other.labels())          # the trap being guarded against
        with self.assertRaises(ValueError):
            verify.check(tagged, other, TARGET, verbose=False)
        self.assertTrue(verify.check(tagged, PRES, TARGET, verbose=False)["pass"])

    def test_pit_check_agrees_independently(self):
        rep = verify.pit_check(self.cof, PRES, TARGET, verbose=False)
        self.assertTrue(rep["pass"])
        self.assertEqual(rep["failures"], [])
        with open(os.path.join(OUT, "pit-report.json"), "w") as f:
            json.dump(rep, f, indent=1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
