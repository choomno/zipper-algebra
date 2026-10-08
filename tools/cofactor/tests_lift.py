#!/usr/bin/env python3
"""Equivalence tests for lift.py against symbolic_lift_x3.py on the X_3 artifacts on disk.

Inputs (read only): logs_and_traces/f4ncgb/basis-n4-Q.txt, proof-n4-Q-noex.txt, the 100
per-point basis-n4-Q-q*-t*.txt / proof-n4-Q-q*-t*-noex.txt, skeleton-n4-Q.json, atoms-n4-Q.json,
and basis-symbolic-Qqt.json (the legacy script's output).  Outputs (lift JSONs, timing) go to
notes/refactor/out/lift/.

Run from zipper-algebra/::

    python3 tools/cofactor/tests_lift.py
    python3 -m tools.cofactor.tests_lift

The points are loaded through a small adapter over the legacy parsers (extract_cofactor
.parse_basis / .parse_noex_proof) shaped like SPEC.md's Basis / Proof / SampleSet, to be
replaced by proof.py / samples.py when they land.  Keys here are the leading word as
'x1*y2*y1'; lift.py never interprets the key text.
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys
import time
import unittest
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py and core.py live
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# The scripts these tests compare against were archived on 2026-09-17 once these modules replaced
# them (notes/REMOVED.md). They stay importable from the archive precisely so the comparison keeps
# running: an equivalence test that can no longer reach the thing it claims equivalence with is
# not a test.
_LEGACY = os.path.join(os.path.dirname(os.path.dirname(HERE)), "archive", "cofactor-pre-consolidation")
if _LEGACY not in sys.path:
    sys.path.insert(0, _LEGACY)

import extract_cofactor as EC            # noqa: E402  legacy parsers + leading_term
import symbolic_lift_x3 as LEGACY        # noqa: E402  the script being replaced
import lift as LIFT                      # noqa: E402
from atoms import AtomSet                # noqa: E402

ZA = os.path.dirname(PARENT)                              # zipper-algebra/
REPO = os.path.dirname(ZA)
LOGS = os.path.join(ZA, "logs_and_traces", "f4ncgb")     # read only

# Legacy-format fixtures come from the frozen snapshot; the pipeline rewrites the live ones.
FROZEN = os.path.join(str(REPO), "notes", "pipeline_test", "pre-atoms-20260917")


def frozen(name: str) -> str:
    p = os.path.join(FROZEN, name)
    return p if os.path.exists(p) else os.path.join(str(LOGS), name)
OUT = os.path.join(REPO, "notes", "refactor", "out", "lift")
VAR_NAME = {v: k for k, v in EC.VAR_INDEX.items()}


# ----------------------------------------------------------------------------
# adapter: legacy parsers -> SPEC-shaped Basis / Proof / SampleSet (swap for proof.py/samples.py)
# ----------------------------------------------------------------------------

class _Basis:
    def __init__(self, lms):
        self.lm = list(lms)
        self.key = ["*".join(VAR_NAME[v] for v in lm) if lm else "1" for lm in self.lm]
        assert len(set(self.key)) == len(self.key), "LM collision: keys need a fingerprint"
        self.index_by_key = {k: i for i, k in enumerate(self.key)}


class _Proof:
    def __init__(self, legacy_terms):
        # legacy words are variable-name lists; SPEC words are 1-based index tuples
        self.terms = [[(c, tuple(EC.VAR_INDEX[v] for v in left), kind, idx,
                        tuple(EC.VAR_INDEX[v] for v in right))
                       for c, left, kind, idx, right in terms] for terms in legacy_terms]


class LegacyX3Samples:
    """Every per-point proof on disk, in the legacy glob (string) order; the anchor (5,7)
    reads the anchor files, which are byte-identical to its per-point files."""
    QT_RE = re.compile(r"proof-n4-Q-q(\d+)-t(\d+)-noex\.txt$")

    def __init__(self, anchor=(5, 7)):
        self.anchor = anchor
        self.points = []
        for pf in sorted(glob.glob(os.path.join(LOGS, "proof-n4-Q-q*-noex.txt"))):
            m = self.QT_RE.search(pf)
            if m:
                self.points.append((int(m.group(1)), int(m.group(2))))
        self._proof, self._basis = {}, {}

    def _paths(self, q, t):
        if (q, t) == self.anchor:
            return (os.path.join(LOGS, "basis-n4-Q.txt"), os.path.join(LOGS, "proof-n4-Q-noex.txt"))
        return (os.path.join(LOGS, f"basis-n4-Q-q{q}-t{t}.txt"),
                os.path.join(LOGS, f"proof-n4-Q-q{q}-t{t}-noex.txt"))

    def proof(self, q, t):
        if (q, t) not in self._proof:
            self._proof[(q, t)] = _Proof(EC.parse_noex_proof(self._paths(q, t)[1]))
        return self._proof[(q, t)]

    def basis(self, q, t):
        if (q, t) not in self._basis:
            polys = EC.parse_basis(self._paths(q, t)[0], 0)
            self._basis[(q, t)] = _Basis(EC.leading_term(p)[0] for p in polys)
        return self._basis[(q, t)]


class _Skeleton:
    def __init__(self, anchor, steps, basis_keys=None):
        self.anchor = anchor
        self.steps = steps
        self.basis_keys = basis_keys


def skeleton_from_disk(samples):
    with open(frozen("skeleton-n4-Q.json")) as f:
        s = json.load(f)
    a = (s["anchor"]["q"], s["anchor"]["t"], s["anchor"]["char"])
    basis = samples.basis(a[0], a[1])
    steps = [{"left": tuple(st["left"]), "basis_key": basis.key[st["basis_index"]],
              "right": tuple(st["right"]), "mult": Fraction(st["mult"])} for st in s["steps"]]
    return _Skeleton(a, steps), s


# ----------------------------------------------------------------------------
# shared fixtures (parsed once: ~4 s for the 100 points)
# ----------------------------------------------------------------------------

class _F:
    samples = None
    skeleton = None
    skel_json = None
    atoms = None
    anchor_legacy = None
    legacy_json = None


def setUpModule():
    os.makedirs(OUT, exist_ok=True)
    _F.samples = LegacyX3Samples()
    _F.skeleton, _F.skel_json = skeleton_from_disk(_F.samples)
    _F.atoms = AtomSet.load(os.path.join(LOGS, "atoms-n4-Q.json"))
    _F.anchor_legacy = LEGACY.parse_noex_proof(os.path.join(LOGS, "proof-n4-Q-noex.txt"))
    with open(frozen("basis-symbolic-Qqt.json")) as f:
        _F.legacy_json = json.load(f)


def _legacy_needed():
    return LEGACY.transitive_deps(_F.anchor_legacy, _F.skel_json["basis_indices"])


def _legacy_stable():
    anchor_skel = [LEGACY.skel(t) for t in _F.anchor_legacy]
    return LEGACY.load_stable_samples(_F.anchor_legacy, anchor_skel, _legacy_needed())


def _num_dict(items):
    return {tuple(m): Fraction(c) for m, c in items}


# ----------------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------------

class TestClosureAndStability(unittest.TestCase):
    def test_closure_equals_legacy_transitive_deps(self):
        s = _F.samples
        a = s.basis(5, 7)
        keys = LIFT.closure([st["basis_key"] for st in _F.skeleton.steps], s.proof(5, 7), a)
        self.assertEqual([a.index_by_key[k] for k in keys], sorted(_legacy_needed()))
        self.assertEqual(len(keys), 25)

    def test_stable_points_equal_legacy_filter(self):
        s = _F.samples
        a = s.basis(5, 7)
        keys = LIFT.closure([st["basis_key"] for st in _F.skeleton.steps], s.proof(5, 7), a)
        mine = LIFT.stable_points(keys, s.proof(5, 7), a, s)
        legacy = [(q, t) for q, t, _ in _legacy_stable()]
        self.assertEqual(mine, legacy)
        self.assertEqual(len(mine), 99)
        self.assertNotIn((8, 7), mine)   # loses a term of anchor[31], [32], [68] at that point

    def test_stability_survives_renumbering(self):
        """A basis listed in another order re-identifies by leading word; the legacy filter,
        comparing indices, would reject it."""
        s = _F.samples
        a, p = s.basis(5, 7), s.proof(5, 7)
        n = len(a.key)
        perm = [(i * 7 + 3) % n for i in range(n)]          # a permutation of 0..n-1
        self.assertEqual(sorted(perm), list(range(n)))
        inv = {perm[i]: i for i in range(n)}
        b2 = _Basis([a.lm[inv[j]] for j in range(n)])         # element j of b2 is old inv[j]
        t2 = [None] * n
        for k in range(n):
            t2[perm[k]] = [(c, l, kind, perm[idx] if kind == "basis" else idx, r)
                           for c, l, kind, idx, r in p.terms[k]]
        p2 = _Proof([])
        p2.terms = t2

        class _S:
            points = [(5, 7), (101, 103)]

            def proof(self, q, t):
                return p if (q, t) == (5, 7) else p2

            def basis(self, q, t):
                return a if (q, t) == (5, 7) else b2

        keys = LIFT.closure([st["basis_key"] for st in _F.skeleton.steps], p, a)
        self.assertEqual(LIFT.stable_points(keys, p, a, _S()), [(5, 7), (101, 103)])
        anchor_skel = [LEGACY.skel(t) for t in _F.anchor_legacy]
        legacy_t2 = [[(c, [VAR_NAME[v] for v in l], kind, idx, [VAR_NAME[v] for v in r])
                      for c, l, kind, idx, r in terms] for terms in t2]
        self.assertTrue(any(LEGACY.skel(legacy_t2[k]) != anchor_skel[k] for k in _legacy_needed()))


class TestFitter(unittest.TestCase):
    def test_fit_rational_equals_legacy_on_every_slot(self):
        """Same first-found denominator and the same numerator dict (zeros included) as
        symbolic_lift_x3.fit_rational, slot by slot, on the legacy fit set."""
        samples = _legacy_stable()
        h = max(4, len(samples) // 10)
        fit = samples[:-h]
        needed = sorted(_legacy_needed())
        cands = _F.atoms.candidates()
        self.assertEqual(cands, LEGACY.DENOMS)
        n_slots = 0
        for k in needed:
            for i in range(len(_F.anchor_legacy[k])):
                pts = [((q, t), pr[k][i][0]) for q, t, pr in fit]
                self.assertEqual(LIFT.fit_rational(pts, _F.atoms, candidates=cands),
                                 LEGACY.fit_rational(pts), (k, i))
                n_slots += 1
        self.assertEqual(n_slots, 231)

    def test_divisibility(self):
        q, t = {(1, 0): 1}, {(0, 1): 1}
        qm1, tq2 = {(1, 0): 1, (0, 0): -1}, {(0, 1): 1, (2, 0): 1}
        prod = {(1, 0): -1, (0, 1): 1, (3, 0): 1, (2, 0): -1, (1, 1): 1, (0, 0): 0}  # (q-1)(t+q^2) - t + t
        prod = {(3, 0): 1, (2, 0): -1, (1, 1): 1, (0, 1): -1}                         # (q-1)(t+q^2)
        self.assertTrue(LIFT._divisible(prod, qm1))
        self.assertTrue(LIFT._divisible(prod, tq2))
        self.assertFalse(LIFT._divisible(prod, q))
        self.assertFalse(LIFT._divisible(prod, t))
        self.assertTrue(LIFT._divisible({(2, 1): 3}, q))
        self.assertFalse(LIFT._divisible({(2, 1): 3, (0, 0): 1}, q))
        self.assertTrue(LIFT._divisible({}, q))
        six = AtomSet.DEFAULT
        self.assertTrue(LIFT.is_reduced(six, (1, 0, 1, 0, 0, 0), {(0, 1): 1, (0, 0): 2}))
        self.assertFalse(LIFT.is_reduced(six, (1, 0, 1, 0, 0, 0), {(1, 1): 1, (1, 0): 2}))

    def test_known_denominator_route_agrees_with_enumeration_per_slot(self):
        samples = _legacy_stable()
        h = max(4, len(samples) // 10)
        fit = samples[:-h]
        cands = _F.atoms.candidates()
        n_known, n_slots, below, other = 0, 0, [], []
        for k in sorted(_legacy_needed()):
            for i in range(len(_F.anchor_legacy[k])):
                n_slots += 1
                pts = [((q, t), pr[k][i][0]) for q, t, pr in fit]
                enum = LIFT.fit_rational(pts, _F.atoms, candidates=cands)
                exps = LIFT.line_denominator(pts, _F.atoms)
                self.assertIsNotNone(exps, (k, i))
                known = LIFT.fit_rational_known(pts, _F.atoms, exps)
                if known is not None:
                    n_known += 1
                    self.assertEqual(known, enum, (k, i))
                    self.assertEqual(exps, enum[0], (k, i))   # the line denominator IS the reduced one
                elif all(e <= f for e, f in zip(exps, enum[0])):
                    below.append((k, i, exps, enum[0]))        # a proper lower bound: a line lost a factor
                else:
                    other.append((k, i, exps, enum[0]))
        self.assertEqual(n_slots, 231)
        self.assertGreater(n_known, 0)
        self.assertEqual(other, [])       # never above the reduced denominator (reducedness check aside)
        with open(os.path.join(OUT, "line-denominator-diagnostics.txt"), "w") as f:
            f.write(f"line_denominator reconstructed the reduced denominator for {n_known}/231 slots;\n"
                    f"{len(below)} were proper lower bounds (fell back to the enumeration):\n")
            for k, i, exps, true in below:
                f.write(f"  basis[{k}] term {i}: lines gave {exps}, enumeration {true}\n")


class TestLift(unittest.TestCase):
    """The whole lift versus the legacy script's output on disk, entry for entry."""

    def _compare_with_legacy_json(self, L):
        s = _F.samples
        a = s.basis(5, 7)
        legacy = _F.legacy_json["lifted"]
        self.assertEqual(sorted(a.index_by_key[k] for k in L.entries), _F.legacy_json["needed"])
        want = ([x["name"] for x in _F.legacy_json["atoms"]["atoms"]]
                if "atoms" in _F.legacy_json else AtomSet.DEFAULT.names())
        self.assertEqual(L.atoms.names(), want)
        n = 0
        for k_str, terms in legacy.items():
            k = int(k_str)
            key = a.key[k]
            self.assertEqual(L.anchor_index[key], k)
            mine = L.entries[key]
            self.assertEqual(len(mine), len(terms), key)
            for i, (it, m) in enumerate(zip(terms, mine)):
                den, num, left, kind, ref, right = m
                self.assertEqual(den, tuple(it["denom"]), (k, i))
                self.assertEqual(num, _num_dict(it["num"]), (k, i))
                self.assertEqual(left, tuple(EC.VAR_INDEX[v] for v in it["left"]), (k, i))
                self.assertEqual(right, tuple(EC.VAR_INDEX[v] for v in it["right"]), (k, i))
                self.assertEqual(kind, it["kind"], (k, i))
                self.assertEqual(ref, a.key[it["idx"]] if kind == "basis" else it["idx"], (k, i))
                if kind == "basis":
                    self.assertEqual(L.anchor_index[ref], it["idx"])
                n += 1
        self.assertEqual(n, 231)

    def test_enumeration_route_equals_legacy_json(self):
        t0 = time.time()
        L = LIFT.lift(_F.skeleton, _F.samples, _F.atoms, verbose=False)
        dt = time.time() - t0
        self.assertEqual(L.samples_used, 99)
        self.assertEqual(L.holdout, 9)
        self.assertEqual(L.timing["routes"], {"known": 0, "enum": 231})
        self._compare_with_legacy_json(L)
        L.save(os.path.join(OUT, "lift-X3-enumeration.json"))
        with open(os.path.join(OUT, "timing.txt"), "a") as f:
            f.write(f"enumeration route: {dt:.2f}s for 231 slots, fit loop {L.timing['seconds']}s\n")

    def test_known_denominator_route_equals_enumeration(self):
        t0 = time.time()
        L2 = LIFT.lift(_F.skeleton, _F.samples, _F.atoms, known="lines", verbose=False)
        dt = time.time() - t0
        routes = L2.timing["routes"]
        self.assertEqual(routes["known"] + routes["enum"], 231)
        self.assertGreater(routes["known"], 0)
        self._compare_with_legacy_json(L2)
        L1 = LIFT.lift(_F.skeleton, _F.samples, _F.atoms, verbose=False)
        self.assertEqual(L1, L2)
        L2.save(os.path.join(OUT, "lift-X3-known.json"))
        with open(os.path.join(OUT, "timing.txt"), "a") as f:
            f.write(f"known-denominator route: {dt:.2f}s for 231 slots, fit loop "
                    f"{L2.timing['seconds']}s of which line reconstruction + one-candidate fit "
                    f"{L2.timing['known_route_seconds']}s; routes {routes}\n")

    def test_json_round_trip(self):
        L = LIFT.lift(_F.skeleton, _F.samples, _F.atoms, verbose=False)
        path = os.path.join(OUT, "lift-X3-roundtrip.json")
        L.save(path)
        L2 = LIFT.Lift.load(path)
        self.assertEqual(L, L2)
        self.assertEqual(L.anchor_index, L2.anchor_index)
        self.assertEqual(L.anchor, L2.anchor)
        with open(path) as f:
            d = json.load(f)
        self.assertEqual(d["entries"]["x1*y2*y1"]["anchor_index"], 68)

    def test_skeleton_closure_mismatch_is_an_error(self):
        sk = _Skeleton(_F.skeleton.anchor, _F.skeleton.steps, basis_keys=["x1*y1"])
        with self.assertRaises(ValueError):
            LIFT.lift(sk, _F.samples, _F.atoms, verbose=False)


class _Synthetic:
    """One basis element, one input term, coefficient f(q,t): samples shaped like the SPEC."""

    def __init__(self, points, f):
        self.points = list(points)
        self.f = f

    def proof(self, q, t):
        p = _Proof([])
        p.terms = [[(self.f(q, t), (), "input", 0, ())]]
        return p

    def basis(self, q, t):
        return _Basis([(1,)])


def _synthetic_skeleton():
    return _Skeleton((5, 7, 0), [{"left": (), "basis_key": "x1", "right": (), "mult": Fraction(1)}])


class TestResultTypes(unittest.TestCase):
    GRID = [(q, t) for q in (5, 3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41) for t in (7, 3, 5, 11, 13)]

    def test_too_few_points_for_holdout(self):
        s = _Synthetic(self.GRID[:6], lambda q, t: Fraction(q * t + 1))
        with self.assertRaises(LIFT.NeedMoreSamples):
            LIFT.lift(_synthetic_skeleton(), s, AtomSet.DEFAULT, verbose=False)

    def test_no_model_in_class(self):
        s = _Synthetic(self.GRID, lambda q, t: Fraction(1, q * q + 1))    # q^2+1 is not an atom
        with self.assertRaises(LIFT.NeedMoreSamples) as cm:
            LIFT.lift(_synthetic_skeleton(), s, AtomSet.DEFAULT, verbose=False)
        self.assertEqual(cm.exception.slots, [("x1", 0)])

    def test_holdout_failure_is_wrong_class(self):
        # 14 points, HOLDOUT 4: the fit points lie on t in {3, 5}, the held-out ones on t = 7,
        # and (t-3)(t-5) q^9 + q is the polynomial q on every fit point and not beyond them.
        pts = [(q, t) for q in (5, 3, 7, 11, 13) for t in (3, 5)] + [(q, 7) for q in (5, 3, 7, 11)]
        s = _Synthetic(pts, lambda q, t: Fraction((t - 3) * (t - 5) * q ** 9 + q))
        with self.assertRaises(LIFT.WrongClass) as cm:
            LIFT.lift(_synthetic_skeleton(), s, AtomSet.DEFAULT, verbose=False)
        self.assertEqual(cm.exception.slots, [("x1", 0, "holdout")])

    def test_enough_points_lift_exactly(self):
        s = _Synthetic(self.GRID, lambda q, t: Fraction(q ** 3 - t, q * (q - 1) ** 2))
        L = LIFT.lift(_synthetic_skeleton(), s, AtomSet.DEFAULT, verbose=False)
        den, num, left, kind, ref, right = L.entries["x1"][0]
        self.assertEqual(den, (1, 0, 2, 0, 0, 0))
        self.assertEqual(num, {(3, 0): Fraction(1), (0, 1): Fraction(-1)})
        self.assertEqual((left, kind, ref, right), ((), "input", 0, ()))
        L2 = LIFT.lift(_synthetic_skeleton(), s, AtomSet.DEFAULT, known="lines", verbose=False)
        self.assertEqual(L, L2)
        self.assertEqual(L2.timing["routes"], {"known": 1, "enum": 0})

    def test_atoms_inferred_when_not_supplied(self):
        s = _Synthetic(self.GRID, lambda q, t: Fraction(q ** 3 - t, q * (q - 1) ** 2))
        L = LIFT.lift(_synthetic_skeleton(), s, None, verbose=False)
        self.assertEqual(L.atoms.names()[:6], AtomSet.DEFAULT.names())
        self.assertIsNotNone(L.atom_report)
        self.assertEqual(L.entries["x1"][0][0][:3], (1, 0, 2))


if __name__ == "__main__":
    unittest.main()
