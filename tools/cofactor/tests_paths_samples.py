#!/usr/bin/env python3
"""Equivalence tests for paths.py and samples.py against the legacy X_3 scripts.

Every assertion compares the new modules with what the legacy scripts compute (their own
constants and functions, imported, not re-typed) or with the artifacts on disk under
logs_and_traces/f4ncgb/.  Nothing is written there; outputs go to notes/refactor/out/paths_samples/.
f4ncgb is never run: the anchor test requires a cache hit and fails before running otherwise.

    cd zipper-algebra && python3 tools/cofactor/tests_paths_samples.py
"""
from __future__ import annotations

import glob
import inspect
import json
import os
import re
import shutil
import subprocess
import sys
import unittest
from itertools import islice
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent            # tools/cofactor
TOOLS = HERE.parent
ZA = TOOLS.parent
REPO = ZA.parent
for p in (str(HERE), str(TOOLS), str(REPO), str(REPO / "papers" / "generic")):
    if p not in sys.path:
        sys.path.insert(0, p)

# The scripts these tests compare against were archived on 2026-09-17 once this module replaced
# them (notes/REMOVED.md). They stay importable from the archive precisely so the comparison
# keeps running: an equivalence test that can no longer reach the thing it claims equivalence
# with is not a test.
_LEGACY = os.path.join(ZA, "archive", "cofactor-pre-consolidation")
if _LEGACY not in sys.path:
    sys.path.insert(0, _LEGACY)

OUT = REPO / "notes" / "refactor" / "out" / "paths_samples"
LOGS = ZA / "logs_and_traces" / "f4ncgb"

# Legacy-format fixtures come from the frozen snapshot; the pipeline rewrites the live ones.
FROZEN = os.path.join(str(REPO), "notes", "pipeline_test", "pre-atoms-20260917")


def frozen(name: str) -> str:
    p = os.path.join(FROZEN, name)
    return p if os.path.exists(p) else os.path.join(str(LOGS), name)

from presentation import Presentation  # noqa: E402
from cache import Cache  # noqa: E402
import paths as P_  # noqa: E402
import samples as S_  # noqa: E402
from paths import Paths  # noqa: E402
from samples import SampleSet, run_f4ncgb, f4ncgb_key, DEFAULT_FLAGS  # noqa: E402

# legacy scripts, imported for their constants and functions
import lift_until_verified as luv  # noqa: E402
import symbolic_lift_x3 as slx  # noqa: E402
import compose_x3_cofactor_qt as cxq  # noqa: E402
import extract_cofactor_qt as ecq  # noqa: E402
import extract_cofactor as ec  # noqa: E402
import fit_symbolic as fs  # noqa: E402   papers/generic


def _mtimes(d: Path) -> dict:
    return {f.name: f.stat().st_mtime_ns for f in d.iterdir()}


class _NoWritesToLogs(unittest.TestCase):
    """Every test class snapshots logs_and_traces/f4ncgb and checks nothing changed."""

    @classmethod
    def setUpClass(cls):
        cls._before = _mtimes(LOGS)
        OUT.mkdir(parents=True, exist_ok=True)

    @classmethod
    def tearDownClass(cls):
        assert _mtimes(LOGS) == cls._before, "a test wrote into logs_and_traces/f4ncgb/"


class TestLegacyPaths(_NoWritesToLogs):

    def setUp(self):
        self.P = Paths.legacy_x3()

    def test_directories_match_legacy_constants(self):
        P = self.P
        self.assertEqual(str(P.dir), str(luv.LOGS))
        self.assertEqual(str(P.dir), slx.PROOF_DIR)
        self.assertEqual(str(P.dir), ecq.LOGS)
        self.assertEqual(str(P.dir), os.path.join(ecq.REPO, "logs_and_traces", "f4ncgb"))
        self.assertEqual(str(P.ms_dir), os.path.join(ecq.REPO, "ms_inputs", "qt_sweep"))

    def test_anchor_names_match_driver_and_extract_cofactor(self):
        P = self.P
        # reproduce_X3_zero.sh phase 1: INPUT_MS, BASIS_OUT, PROOF_OUT; phase 2: SKEL
        self.assertEqual(P.anchor_ms, LOGS / "zipper-n4-Q-q5-t7.ms")
        self.assertEqual(P.anchor_basis, LOGS / "basis-n4-Q.txt")
        self.assertEqual(P.anchor_proof, LOGS / "proof-n4-Q-noex.txt")
        self.assertEqual(P.skeleton, LOGS / "skeleton-n4-Q.json")
        # extract_cofactor.py main(): out_path for (3,1) and for another slider
        self.assertEqual(P.cofactor_anchor, LOGS / "cofactor-X3-n4-Q.txt")
        P41 = Paths.legacy_x3(target="4,1")
        self.assertEqual(P41.target, "S[4][1]")
        self.assertEqual(P41.skeleton, LOGS / "skeleton-n4-S4_1-Q.json")
        self.assertEqual(P41.cofactor_anchor, LOGS / "cofactor-S4_1-n4-Q.txt")
        self.assertTrue(P41.skeleton.is_file() and P41.cofactor_anchor.is_file())

    def test_lift_and_compose_names_match_legacy_constants(self):
        P = self.P
        self.assertEqual(str(P.atoms), slx.ATOMS_PATH)
        self.assertEqual(str(P.atoms_report), slx.ATOMS_PATH[:-5] + "-report.json")
        self.assertEqual(str(P.lift), cxq.JSON_PATH)
        self.assertEqual(str(P.lift), os.path.join(slx.PROOF_DIR, "basis-symbolic-Qqt.json"))
        self.assertEqual(str(P.lift_txt), os.path.join(slx.PROOF_DIR, "basis-symbolic-Qqt.txt"))
        # compose_x3_cofactor_qt.py out_path; reproduce_X3_zero.sh SYM_OUT/COMP_OUT/SIMP_OUT/THREE_OUT/ZT_OUT
        self.assertEqual(P.cofactor(), LOGS / "cofactor-X3-n4-Qqt-symbolic.txt")
        for ph in ("compressed", "simplified", "three-rels", "ZT"):
            self.assertEqual(P.cofactor(ph), LOGS / f"cofactor-X3-n4-Qqt-{ph}.txt")
        self.assertEqual(P.transcript, ZA / "logs_and_traces" / "reproduce-X3-zero-transcript.txt")

    def test_per_point_names_match_extract_cofactor_qt(self):
        P = self.P
        for (q, t) in ((5, 7), (11, 13), (31, 37)):
            # extract_cofactor_qt.py main(): ms_path, basis_path, proof_path, cofactor_path
            self.assertEqual(str(P.ms(q, t)), os.path.join(ecq.REPO, "ms_inputs", "qt_sweep",
                                                           f"zipper-n4-Q-q{q}-t{t}.ms"))
            self.assertEqual(str(P.basis(q, t)), os.path.join(ecq.LOGS, f"basis-n4-Q-q{q}-t{t}.txt"))
            self.assertEqual(str(P.proof(q, t)), os.path.join(ecq.LOGS, f"proof-n4-Q-q{q}-t{t}-noex.txt"))
            self.assertEqual(str(P.cofactor_at(q, t)),
                             os.path.join(ecq.LOGS, f"cofactor-X3-n4-Q-q{q}-t{t}.txt"))
        self.assertEqual(P.proof_re.pattern, luv.PROOF_RE.pattern)
        self.assertEqual(P.proof_glob, "proof-n4-Q-q*-t*-noex.txt")     # lift_until_verified.have_points

    def test_every_legacy_artifact_exists_on_disk(self):
        P = self.P
        for name in ("anchor_ms", "anchor_basis", "anchor_proof", "skeleton", "atoms",
                     "atoms_report", "lift", "lift_txt", "cofactor_anchor", "transcript"):
            self.assertTrue(getattr(P, name).is_file(), name)
        for ph in ("symbolic", "compressed", "simplified", "three-rels", "ZT"):
            self.assertTrue(P.cofactor(ph).is_file(), ph)
        pts = P.points_on_disk()
        self.assertEqual(len(pts), 100)
        for (q, t) in pts:
            self.assertTrue(P.ms(q, t).is_file(), (q, t))
            self.assertTrue(P.basis(q, t).is_file(), (q, t))
            self.assertTrue(P.proof(q, t).is_file(), (q, t))

    def test_legacy_is_char0_only(self):
        with self.assertRaises(ValueError):
            Paths.legacy_x3(anchor=(5, 7, 65521))


class TestNewLayout(_NoWritesToLogs):

    def test_names_and_no_side_effects(self):
        root = OUT / "root"
        shutil.rmtree(root, ignore_errors=True)
        P = Paths(Presentation.DEFAULT, "S[3][1]", (5, 7, 0), root=root)
        self.assertFalse(root.exists())                     # nothing created by construction
        d = root / "n4-bigelow"
        self.assertEqual(P.dir, d)
        self.assertEqual(P.anchor_ms, d / "input-5-7-c0.ms")
        self.assertEqual(P.ms(11, 13), d / "input-11-13-c0.ms")
        self.assertEqual(P.anchor_basis, P.basis(5, 7))
        self.assertEqual(P.basis(11, 13), d / "basis-11-13.txt")
        self.assertEqual(P.proof(11, 13), d / "proof-11-13-noex.txt")
        self.assertEqual(P.skeleton, d / "skeleton-S3_1.json")
        self.assertEqual(P.atoms, d / "atoms.json")
        self.assertEqual(P.lift, d / "lift-S3_1.json")
        self.assertEqual(P.cofactor(), d / "cofactor-S3_1-symbolic.txt")
        self.assertEqual(P.points_on_disk(), [])
        P.ensure_dirs()
        (d / "proof-11-13-noex.txt").write_text("")
        (d / "proof-3-5-noex.txt").write_text("")
        self.assertEqual(P.points_on_disk(), [(11, 13), (3, 5)])   # file-name order, as legacy
        # characteristic p lives in the directory name; the .ms carries it too
        Pp = Paths(Presentation.DEFAULT, "S[3][1]", (5, 7, 65521), root=root)
        self.assertEqual(Pp.dir, root / "n4-bigelow-c65521")
        self.assertEqual(Pp.anchor_ms.name, "input-5-7-c65521.ms")
        # other targets and presentations
        Px = Paths(Presentation(4, extra=("(Theta_3-q^2)X_3",)), "EXTRA[(Theta_3-q^2)X_3]",
                   (5, 7, 0), root=root)
        self.assertEqual(Px.skeleton.name, "skeleton-EXTRA_Theta_3_q_2_X_3.json")
        self.assertEqual(Px.dir.name, Presentation(4, extra=("(Theta_3-q^2)X_3",)).slug())
        self.assertEqual(P_.normalise_target((4, 1)), "S[4][1]")
        self.assertEqual(P_.normalise_target("4,1"), "S[4][1]")

    def test_default_root_from_settings(self):
        from repo_settings import settings
        self.assertEqual(P_.default_root(), settings().path("output") / "cofactor")


class TestSampleSet(_NoWritesToLogs):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.S = SampleSet(Paths.legacy_x3())

    def test_points_equal_legacy_glob(self):
        S = self.S
        self.assertEqual(set(S.points), luv.have_points())
        self.assertEqual(len(S.points), 100)
        # order: symbolic_lift_x3.load_stable_samples iterates sorted(glob(...))
        qt_re = re.compile(r"proof-n4-Q-q(\d+)-t(\d+)-noex.txt")
        legacy_order = []
        for pf in sorted(glob.glob(os.path.join(slx.PROOF_DIR, "proof-n4-Q-q*-noex.txt"))):
            m = qt_re.search(pf)
            if m:
                legacy_order.append((int(m.group(1)), int(m.group(2))))
        self.assertEqual(S.points, legacy_order)

    def test_anchor_run_hits_cache_with_identical_bytes(self):
        S = self.S
        P = S.paths
        # phase 1a into OUT: same basename as the driver's INPUT_MS, so the key is the driver's
        ms = OUT / P.anchor_ms.name
        Presentation.DEFAULT.emit_ms(str(ms), 5, 7, 0)
        self.assertEqual(ms.read_bytes(), P.anchor_ms.read_bytes())
        key = f4ncgb_key(ms)
        self.assertEqual(key, f4ncgb_key(P.anchor_ms))
        # the driver's own key command
        from repo_settings import settings
        sh = subprocess.run([sys.executable, str(TOOLS / "cache.py"), "key", "--file", str(ms),
                             "--tool", settings().binary("f4ncgb"), "--tag", "f4ncgb", "--",
                             *DEFAULT_FLAGS], capture_output=True, text=True, check=True)
        self.assertEqual(sh.stdout.strip(), key)
        slot = Cache().root / key / "manifest.json"
        self.assertTrue(slot.is_file(),
                        f"no cache slot {key}: refusing to run f4ncgb inside a test "
                        f"(run reproduce_X3_zero.sh phase 1b once)")
        basis, proof = OUT / P.anchor_basis.name, OUT / P.anchor_proof.name
        for f in (basis, proof):
            if f.exists():
                f.unlink()
        r = run_f4ncgb(ms, basis, proof)
        self.assertTrue(r.hit, r)
        self.assertEqual(r.key, key)
        self.assertEqual(basis.read_bytes(), P.anchor_basis.read_bytes())
        self.assertEqual(proof.read_bytes(), P.anchor_proof.read_bytes())

    def test_check_ms_accepts_emitter_output_and_rejects_malformed(self):
        S_.check_ms(self.S.paths.anchor_ms)
        bad = OUT / "malformed.ms"
        bad.write_text("x1,y1\n0\nx1*y1 + -1,\nx1,\n")     # sign run and trailing comma
        with self.assertRaises(ValueError):
            S_.check_ms(bad)

    def test_spiral_is_verbatim(self):
        self.assertEqual(list(islice(S_.spiral(), 30)), list(islice(fs.spiral(), 30)))
        for name in ("primes", "spiral", "_spiral_indices"):
            self.assertEqual(inspect.getsource(getattr(S_, name)),
                             inspect.getsource(getattr(fs, name)), name)
        self.assertEqual(S_.PRIMES, fs.PRIMES)

    def test_grow_selection(self):
        S = self.S
        pts = S.next_points(12)
        self.assertEqual(len(pts), 12)
        have = luv.have_points()
        expected = []
        for (q, t) in fs.spiral():
            if (q, t) in have or q in (2, 4) or any(a.eval(q, t) == 0 for a in S.atoms.atoms):
                continue
            expected.append((q, t))
            if len(expected) == 12:
                break
        self.assertEqual(pts, expected)
        # exclusion rules, stated independently
        self.assertFalse(S_.admissible(2, 3, S.atoms))          # q in {2,4}
        self.assertFalse(S_.admissible(7, 7, S.atoms))          # t - q
        self.assertFalse(S_.admissible(3, -9, S.atoms))         # t + q^2
        self.assertTrue(S_.admissible(3, 5, S.atoms))
        # grow of nothing runs nothing
        self.assertEqual(S.grow(0), [])

    def test_stable_equals_legacy_load_stable_samples(self):
        S = self.S
        P = S.paths
        ec.configure(4)
        anchor_main = ec.parse_noex_proof(str(P.anchor_proof))
        anchor_skel = [slx.skel(t) for t in anchor_main]
        # the frozen legacy-format skeleton: the live one is rewritten by the pipeline
        seeds = json.loads(open(frozen("skeleton-n4-Q.json")).read())["basis_indices"]
        self.assertEqual(seeds, [16, 17, 102])
        needed = slx.transitive_deps(anchor_main, seeds)
        self.assertEqual(len(needed), 25)                       # transcript: "lifted entries for 25"
        legacy = [(q, t) for q, t, _ in slx.load_stable_samples(anchor_main, anchor_skel, needed)]
        ab = S.anchor_basis()
        self.assertEqual([ab.key[k] for k in seeds], ["y1*y2*x1", "y1*y2*y1", "x1*x2*y1"])
        keys = [ab.key[k] for k in sorted(needed)]
        mine = S.stable(SimpleNamespace(basis_keys=keys))
        self.assertEqual(mine, legacy)
        self.assertEqual(S.stable({"basis_keys": keys}), legacy)
        # transitive closure through the Proof type agrees with symbolic_lift_x3.transitive_deps
        ap = S.anchor_proof()
        closure = set(seeds).union(*(ap.deps(k) for k in seeds))
        self.assertEqual(closure, needed)
        # a key with a different structure at every point is stable nowhere but the anchor copy
        with self.assertRaises(KeyError):
            S.stable(SimpleNamespace(basis_keys=["no*such*word"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
