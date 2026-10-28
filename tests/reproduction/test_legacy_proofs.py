#!/usr/bin/env python3
"""Equivalence tests for proof.py and skeleton.py against extract_cofactor.py on the X3 artifacts.

Every assertion is an EQUALITY between the legacy script's output and the new modules' on
logs_and_traces/f4ncgb/basis-n4-Q.txt and proof-n4-Q-noex.txt (the (5,7) anchor over Q), plus
synthetic cases for what the legacy could not do (forward references; leading-word collisions).
Inputs are immutable; outputs use a temporary directory.

    cd zipper-algebra && python3 -m unittest tests.reproduction.test_legacy_proofs -v
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys
import tempfile
import unittest
from fractions import Fraction

from tests.support import ROOT, LOGS as LOG_DIR, LEGACY_X3
from tests.reference import legacy_x3 as E
from zipper_algebra.cofactor.proof import Basis, Proof, mon_sort_key, word_text
from zipper_algebra.cofactor.skeleton import Skeleton, reduce, compose_at, evaluate, closure
from zipper_algebra.presentation import Presentation
from zipper_algebra.zipper_in_ms import NumericRing, Poly, Coef

LOGS = str(LOG_DIR)
FROZEN = str(LEGACY_X3)


def frozen(name: str) -> str:
    p = os.path.join(FROZEN, name)
    return p if os.path.exists(p) else os.path.join(LOGS, name)
_OUTPUT = tempfile.TemporaryDirectory(prefix="zipper-proof-tests-")
OUT = _OUTPUT.name


def tearDownModule():
    _OUTPUT.cleanup()
BASIS_PATH = os.path.join(LOGS, "basis-n4-Q.txt")
PROOF_PATH = os.path.join(FROZEN, "proof-n4-Q-noex.txt")
ANCHOR = (5, 7, 0)

def _terms_of(p: Poly) -> dict:
    return {m: c.v for m, c in p.terms.items()}


def _json(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


class _Shared:
    """Parsed once for the whole run (the legacy expansion is the slow part)."""
    pres = basis = proof = legacy_basis = legacy_noex = legacy_expanded = None
    target = inputs = None

    @classmethod
    def setup(cls):
        if cls.basis is not None:
            return
        os.makedirs(OUT, exist_ok=True)
        cls.pres = Presentation(4)
        names = cls.pres.var_names()
        E.configure(4)
        cls.basis = Basis.load(BASIS_PATH, 0, names)
        cls.proof = Proof.load(PROOF_PATH, names)
        cls.legacy_basis = E.parse_basis(BASIS_PATH, 0)
        cls.legacy_noex = E.parse_noex_proof(PROOF_PATH)
        cls.legacy_expanded = E.expand_proofs(cls.legacy_noex)
        ring = NumericRing(*ANCHOR)
        _, cls.inputs, S, _ = cls.pres.build(ring)
        cls.target = S[(3, 1)]


class TestBasisAndProofParsers(unittest.TestCase):
    def setUp(self):
        _Shared.setup()

    def test_1_basis_polys_equal_legacy_term_for_term(self):
        b, lb = _Shared.basis, _Shared.legacy_basis
        self.assertEqual(len(b), len(lb))
        self.assertEqual(len(b), 137)
        for k, (p, lp) in enumerate(zip(b.polys, lb)):
            self.assertEqual(_terms_of(p), _terms_of(lp), f"basis[{k}]")
            self.assertEqual(b.lm[k], E.leading_term(lp)[0], f"lm of basis[{k}]")

    def test_1b_keys_are_leading_words_and_distinct(self):
        b = _Shared.basis
        self.assertEqual(len(set(b.key)), len(b))
        for k in range(len(b)):
            self.assertEqual(b.key[k], word_text(b.lm[k], b.var_names))
            self.assertEqual(b.index_by_key[b.key[k]], k)
        self.assertEqual([b.key[k] for k in (16, 17, 102)],
                         ["y1*y2*x1", "y1*y2*y1", "x1*x2*y1"])

    def test_2_proof_terms_equal_legacy(self):
        pr, lp = _Shared.proof, _Shared.legacy_noex
        self.assertEqual(len(pr), len(lp))
        for k in range(len(pr)):
            want = [(c, tuple(E.VAR_INDEX[n] for n in l), kind, idx,
                     tuple(E.VAR_INDEX[n] for n in r)) for (c, l, kind, idx, r) in lp[k]]
            self.assertEqual(pr.terms[k], want, f"proof line {k}")


class TestExpansion(unittest.TestCase):
    def setUp(self):
        _Shared.setup()

    def test_3_expand_equals_legacy_for_every_k(self):
        pr, le = _Shared.proof, _Shared.legacy_expanded
        pr.expand(len(pr) - 1)                        # topological path over the whole graph
        for k in range(len(pr)):
            mine = pr.expand(k)
            self.assertEqual({(l, i, r): c for c, l, i, r in mine},
                             {(l, i, r): c for c, l, i, r in le[k]}, f"expand({k})")
            self.assertEqual(mine, le[k], f"expand({k}) term order")

    def test_3b_every_expansion_reconstructs_its_basis_element(self):
        b, pr = _Shared.basis, _Shared.proof
        for k in range(len(b)):
            self.assertTrue(pr.reconstructs(k, b, _Shared.inputs), f"basis[{k}]")

    def test_3c_deps_match_legacy_transitive_closure(self):
        pr = _Shared.proof
        needed = _json(frozen("basis-symbolic-Qqt.json"))["needed"]
        self.assertEqual(sorted(closure([16, 17, 102], pr)), needed)
        for k in range(len(pr)):
            refs = {idx for _c, _l, kind, idx, _r in pr.terms[k] if kind == "basis"}
            self.assertTrue(refs <= pr.deps(k))
            for j in refs:
                self.assertTrue(pr.deps(j) <= pr.deps(k))

    def test_forward_reference_is_expanded_topologically(self):
        # basis[0] -> basis[2] -> basis[1]: the legacy raised on any forward reference
        terms = [[(Fraction(1), (), "basis", 2, ())],
                 [(Fraction(2), (1,), "input", 0, ())],
                 [(Fraction(3), (), "basis", 1, (2,)), (Fraction(1), (), "input", 1, ())]]
        pr = Proof(terms)
        self.assertEqual(pr.expand_dict(0), {((1,), 0, (2,)): Fraction(6), ((), 1, ()): Fraction(1)})
        self.assertEqual(pr.deps(0), {1, 2})
        self.assertEqual(pr.deps(1), set())
        names = ["x1", "x2"]
        legacy = [[(c, [names[i - 1] for i in l], k, i_, [names[i - 1] for i in r])
                   for c, l, k, i_, r in line] for line in terms]
        with self.assertRaises(NotImplementedError):
            E.expand_proofs(legacy)

    def test_cycle_is_an_error(self):
        pr = Proof([[(Fraction(1), (), "basis", 1, ())], [(Fraction(1), (), "basis", 0, ())]])
        with self.assertRaises(ValueError):
            pr.expand(0)
        self.assertEqual(pr.deps(0), {0, 1})


class TestSkeleton(unittest.TestCase):
    def setUp(self):
        _Shared.setup()
        self.steps, self.rem = reduce(_Shared.target, _Shared.basis)
        self.skel = Skeleton.from_reduction(_Shared.pres, "S[3][1]", ANCHOR, self.steps,
                                            self.rem, _Shared.basis, _Shared.proof)

    def test_4_reduce_equals_legacy(self):
        rem_l, steps_l = E.reduce_against_basis(_Shared.target, _Shared.legacy_basis, 0)
        self.assertEqual([(u, k, v, m.v) for u, k, v, m in self.steps],
                         [(u, k, v, m.v) for u, k, v, m in steps_l])
        self.assertEqual(_terms_of(self.rem), _terms_of(rem_l))
        self.assertEqual(len(self.steps), 3)
        self.assertEqual(self.rem.terms, {})

    def test_4b_x3_skeleton_is_the_historical_one(self):
        b = _Shared.basis
        s = self.skel
        self.assertEqual(s.step_keys(), [b.key[17], b.key[16], b.key[102]])
        self.assertEqual({st["basis_key"]: st["mult"] for st in s.steps},
                         {b.key[16]: Fraction(-25), b.key[17]: Fraction(125), b.key[102]: Fraction(20)})
        self.assertTrue(all(st["left"] == () and st["right"] == () for st in s.steps))
        self.assertEqual(s.remainder_terms, 0)
        self.assertTrue(s.in_ideal)
        self.assertEqual(s.step_indices(b), [17, 16, 102])
        legacy = _json(frozen("skeleton-n4-Q.json"))
        self.assertEqual([b.key[st["basis_index"]] for st in legacy["steps"]], s.step_keys())
        self.assertEqual([str(st["mult"]) for st in s.steps], [st["mult"] for st in legacy["steps"]])

    def test_4c_closure_is_the_25_needed_elements(self):
        b = _Shared.basis
        needed = _json(frozen("basis-symbolic-Qqt.json"))["needed"]
        self.assertEqual(len(needed), 25)
        self.assertEqual(self.skel.basis_keys, [b.key[k] for k in needed])
        self.assertEqual(self.skel.closure_indices(b), needed)

    def test_4d_json_round_trip_carries_the_required_fields(self):
        path = os.path.join(OUT, "skeleton-S3_1-n4-bigelow.json")
        self.skel.save(path)
        d = _json(path)
        self.assertEqual(d["presentation"], "n4-bigelow")
        self.assertEqual(d["presentation_hash"], _Shared.pres.content_hash())
        self.assertEqual(d["target"], "S[3][1]")
        self.assertEqual(d["anchor"], {"q": 5, "t": 7, "char": 0})
        self.assertEqual(d["remainder_terms"], 0)
        self.assertEqual(len(d["steps"]), 3)
        self.assertEqual(len(d["basis_keys"]), 25)
        self.assertNotIn("basis_index", json.dumps(d))          # indices are never stored
        back = Skeleton.load(path)
        self.assertEqual(back, self.skel)
        self.assertEqual(back.steps[0]["mult"], Fraction(125))
        self.assertEqual(back.anchor, ANCHOR)

    def test_5_composition_reconstructs_x3_with_695_tuples(self):
        terms = compose_at(self.skel, _Shared.basis, _Shared.proof)
        self.assertEqual(len(terms), 695)
        self.assertTrue(E.verify_cofactor(terms, _Shared.inputs, _Shared.target, 0))
        got = evaluate(terms, _Shared.inputs, 0)
        self.assertEqual(_terms_of(got), _terms_of(_Shared.target))
        _rem, steps_l = E.reduce_against_basis(_Shared.target, _Shared.legacy_basis, 0)
        legacy = E.compose_cofactor(steps_l, _Shared.legacy_expanded, 0)
        self.assertEqual({(l, i, r): c for c, l, i, r in terms},
                         {(l, i, r): c for c, l, i, r in legacy})

    def test_reduce_identity_for_an_arbitrary_target(self):
        # Target is ANY Poly: here a word times X_3 times a word plus a slider. The reduction
        # identity target = sum mult * left * basis[k] * right + remainder must hold exactly,
        # and the composed cofactor must rebuild target - remainder from the inputs.
        b, pr = _Shared.basis, _Shared.proof
        _, _, S, _ = _Shared.pres.build(NumericRing(*ANCHOR))
        w1 = Poly(0); w1.add_term((1, 2, 1, 4), Coef.one(0))
        w2 = Poly(0); w2.add_term((5, 3), Coef.one(0))
        target = w1 * S[(3, 1)] * w2 + S[(2, 2)].scalar_mul(Coef(7, 0))
        steps, rem = reduce(target, b)
        self.assertGreater(len(steps), 3)
        acc = Poly(0)
        for u, k, v, mult in steps:
            lp = Poly(0); lp.add_term(u, Coef.one(0))
            rp = Poly(0); rp.add_term(v, Coef.one(0))
            acc = acc + (lp * b.polys[k] * rp).scalar_mul(mult)
        self.assertEqual(_terms_of(acc + rem), _terms_of(target))
        skel = Skeleton.from_reduction(_Shared.pres, "x1*x2*x1*y1 * S[3][1] * y2*x3 + 7 S[2][2]",
                                       ANCHOR, steps, rem, b, pr)
        got = evaluate(compose_at(skel, b, pr), _Shared.inputs, 0)
        self.assertEqual(_terms_of(got + rem), _terms_of(target))
        self.assertTrue(set(skel.step_keys()) <= set(skel.basis_keys))

    def test_keys_reidentify_at_other_points(self):
        # Reordering the actual anchor basis changes numeric references but not
        # leading-word identities. Per-point data below supply a separate census.
        names = _Shared.pres.var_names()
        reordered = Basis(list(reversed(_Shared.basis.polys)), 0, names)
        last = len(_Shared.basis) - 1
        self.assertEqual(self.skel.step_indices(reordered),
                         [last - index for index in (17, 16, 102)])
        self.assertEqual(self.skel.closure_indices(reordered),
                         [last - index for index in self.skel.closure_indices(_Shared.basis)])
        ok, missing = 0, []
        for path in sorted(glob.glob(os.path.join(LOGS, "basis-n4-Q-q*-t*.txt"))):
            q, t = map(int, re.search(r"q(\d+)-t(\d+)", path).groups())
            bp = Basis.load(path, 0, names)
            try:
                self.skel.step_indices(bp); self.skel.closure_indices(bp)
                ok += 1
            except KeyError as e:
                missing.append(((q, t), str(e)))
        with open(os.path.join(OUT, "reidentification.txt"), "w") as f:
            f.write(f"points whose basis carries every skeleton key: {ok}/{ok + len(missing)}\n")
            for pt, why in missing:
                f.write(f"  {pt}: {why}\n")
        self.assertGreater(ok, 0)


class TestKeyCollisions(unittest.TestCase):
    def test_colliding_leading_words_get_fingerprints(self):
        names = ["x1", "x2"]
        def poly(terms):
            p = Poly(0)
            for m, c in terms:
                p.add_term(m, Coef(c, 0))
            return p
        a = poly([((1, 2), 1), ((), -1)])          # x1*x2 - 1
        b = poly([((1, 2), 1), ((2,), -1)])        # x1*x2 - x2
        c = poly([((2, 2), 1)])                    # x2*x2
        basis = Basis([a, b, c, a], 0, names)
        self.assertEqual(basis.key, ["x1*x2 | x1*x2-1", "x1*x2 | x1*x2-x2", "x2*x2",
                                     "x1*x2 | x1*x2-1 #2"])
        self.assertEqual(len(set(basis.key)), 4)
        with self.assertRaises(KeyError):
            basis.index("x1*x1")


if __name__ == "__main__":
    unittest.main(verbosity=2)
