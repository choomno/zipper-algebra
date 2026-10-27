#!/usr/bin/env python3
"""The cofactor pipeline: one driver, any presentation, any target.

    python3 -m tools.cofactor.pipeline --n 4 --target 3,1
    python3 -m tools.cofactor.pipeline --n 4 --untwist square-zero --max-untwist-level 3 --target 3,1
    python3 -m tools.cofactor.pipeline --legacy-x3          # the X_3 run, phases 1-10, in place

Phases 1-5 produce a certificate: a finite list of tuples (L, relation, R, c(q,t)) with

    sum_j  c_j(q,t) * L_j * rel_{i_j} * R_j  =  target      in Q(q,t)<x_1..y_{n-1}>,

verified with q,t transcendental.  Phases 6-9 reduce a certificate that already exists; they are
still the X_3-at-n=4 scripts and run only under --legacy-x3.

What replaced what.  This driver supersedes `reproduce_X3_zero.sh`, which was a bash script whose
ten phases called nine scripts that shared state through fixed file names in one flat directory,
hard-coded `n=4`, the anchor `(5,7)`, the basis indices `16/17/102`, the multipliers `q^3, -q^2,
q(q-1)`, and a six-atom denominator class.  All of those are now arguments or computed data:
`presentation.Presentation`, `paths.Paths`, `samples.SampleSet`, `skeleton.Skeleton`,
`atoms.AtomSet`, `lift.Lift`, `cofactor_io.Cofactor`.

The pass criterion is phase 5 and nothing else.  A lift that fits every sampled point can still be
wrong -- that is how a cofactor fitted on 37 points came to disagree with the exact answer in 695
of 695 slots -- so the loop grows the grid and re-lifts until the symbolic identity holds, and
reports failed reconstruction as inconclusive. A failed candidate can require more independent
data, a revised hypothesis class, or a correction to the producer.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from fractions import Fraction

_HERE = os.path.dirname(os.path.abspath(__file__))
_TOOLS = os.path.dirname(_HERE)
_ZA = os.path.dirname(_TOOLS)
_REPO = os.path.dirname(_ZA)
if __package__:
    from . import compose as compose_mod, verify as verify_mod
    from .atoms import AtomSet
    from .lift import NeedMoreSamples, WrongClass, lift as do_lift
    from .paths import Paths
    from ..presentation import Presentation
    from .skeleton import Skeleton, reduce as reduce_target
    from ..zipper_in_ms import NumericRing
else:
    for _p in (_HERE, _TOOLS, _REPO):
        if _p not in sys.path:
            sys.path.insert(0, _p)
    import compose as compose_mod
    import verify as verify_mod
    from atoms import AtomSet
    from lift import NeedMoreSamples, WrongClass, lift as do_lift
    from paths import Paths
    from presentation import Presentation
    from skeleton import Skeleton, reduce as reduce_target
    from zipper_in_ms import NumericRing


BAR = "=" * 64


def phase(n, title):
    print(f"\n{BAR}\nPHASE {n}: {title}\n{BAR}", flush=True)


class NotInIdeal(Exception):
    """Compatibility exception for an unresolved nonzero reduction remainder.

    A nonzero remainder proves non-membership only against a certified complete basis.
    The experimental pipeline does not supply such a certificate, so this outcome is
    inconclusive, including at the specialized anchor.
    """

    def __init__(self, remainder, skeleton):
        super().__init__(f"uncertified target reduction: remainder has "
                         f"{len(remainder.terms)} terms")
        self.remainder = remainder
        self.skeleton = skeleton


def _inconclusive(pres, paths, anchor, maxpts, reason, cofactor=None):
    report = {"status": "inconclusive", "verified": False, "in_ideal": None,
              "presentation": pres.slug(), "presentation_hash": pres.content_hash(),
              "target": paths.target, "anchor": list(anchor), "max_points": maxpts,
              "reason": reason}
    if cofactor is not None:
        report["candidate_cofactor"] = str(cofactor)
    return report


def _unresolved_remainder(pres, paths, anchor, maxpts, skeleton):
    report = _inconclusive(pres, paths, anchor, maxpts,
                           "nonzero remainder against an uncertified basis")
    report.update({"remainder_terms": skeleton.remainder_terms,
                   "basis_completeness": "uncertified", "skeleton": str(paths.skeleton)})
    return report


def run(pres: Presentation, target, anchor, *, start=100, step=40, maxpts=260,
        legacy=False, refresh=False, jobs=None, atoms_path=None, verbose=True) -> dict:
    if any(type(value) is not int or value < 1 for value in (start, step, maxpts)):
        raise ValueError("start, step and maxpts must be positive integers")
    if start > maxpts:
        raise ValueError("start must not exceed maxpts")
    if len(anchor) != 3 or type(anchor[2]) is not int or anchor[2] != 0:
        raise ValueError("rational reconstruction requires a characteristic-zero anchor")
    if __package__:
        from .samples import SampleSet
    else:
        from samples import SampleSet
    t0 = time.time()
    paths = Paths.legacy_x3(pres, target, anchor) if legacy else Paths(pres, target, anchor)
    paths.ensure_dirs()
    print(f"presentation : {pres.slug()}  hash {pres.content_hash()}  "
          f"{pres.n_relations()} relations")
    print(f"target       : {paths.target}")
    print(f"anchor       : (q,t,char) = {anchor}")
    print(f"artifacts    : {paths.root}")

    # --- 1. emit ------------------------------------------------------------
    phase(1, "Emit the presentation at the anchor")
    ms = paths.anchor_ms
    n_rel = pres.emit_ms(str(ms), anchor[0], anchor[1], anchor[2])
    print(f"{n_rel} relations -> {ms.name}")

    # --- 1b. Groebner basis + proof ----------------------------------------
    phase("1b", "Groebner basis with proof logging (cached on the input bytes)")
    samples = SampleSet(paths, atoms=AtomSet.load(atoms_path) if atoms_path else None)
    res = samples.run_anchor(refresh=refresh)
    print(f"{'[cache] HIT ' if res.hit else '[cache] MISS'}  {res.key}"
          f"{'' if res.hit else f'  ({res.seconds:.1f}s)'}")
    basis = samples.anchor_basis()
    proof = samples.anchor_proof()
    print(f"[size] basis: {len(basis)} elements; proof: {len(proof)} lines", flush=True)

    # --- 2. skeleton --------------------------------------------------------
    phase(2, "Skeleton: reduce the target against the anchor basis")
    ring = NumericRing(anchor[0], anchor[1], anchor[2])
    tgt = compose_mod.resolve_target(pres, paths.target, ring)
    steps, remainder = reduce_target(tgt, basis)
    skel = Skeleton.from_reduction(pres, paths.target, anchor, steps, remainder, basis, proof)
    skel.save(str(paths.skeleton))
    print(f"[skeleton] {paths.target} = {len(skel.steps)} step(s) on {skel.step_keys()}, "
          f"remainder {skel.remainder_terms} terms")
    print(f"[skeleton] closure: {len(skel.basis_keys)} basis elements -> {paths.skeleton.name}", flush=True)
    if not skel.in_ideal:
        return _unresolved_remainder(pres, paths, anchor, maxpts, skel)

    # --- 3-5. lift, compose, verify; grow the grid until phase 5 passes -----
    target_pts = min(maxpts, max(start, len(paths.points_on_disk())))
    attempt = 0
    while True:
        attempt += 1
        phase(3, f"Multi-point lift (attempt {attempt}, target {target_pts} points)")
        have = paths.points_on_disk()
        if len(have) < target_pts:
            print(f"[sample] have {len(have)}, generating {target_pts - len(have)} more ...",
                  flush=True)
            samples.grow(target_pts - len(have), jobs=jobs)
        print(f"[sample] {len(paths.points_on_disk())} points on disk", flush=True)

        # The anchor must be GENERIC. A coefficient that is generically non-zero can vanish at
        # the anchor, collapsing a proof term; then the anchor's term structure is the exception,
        # every other point looks "unstable", and -- worse -- the lift would size each slot list
        # by the anchor's term count and silently drop the collapsed one. Measured 2026-09-17 at
        # n=3 square-zero: (5,7) is degenerate for `x2*x1*y2*x1` (109 proof terms against 110
        # everywhere else) and 2 of 160 points survived. So: census the shapes, and if the anchor
        # is not modal, move to a point that is and redo phase 2 there.
        modal, off_modal = samples.modal_shapes(skel)
        if off_modal:
            candidates = samples.stable(skel)
            if not candidates:
                reason = ("no sample carries the modal proof shape for every skeleton "
                          "element; no stable reconstruction is available")
                print(f"[anchor] INCONCLUSIVE: {reason}", flush=True)
                if target_pts >= maxpts:
                    return _inconclusive(pres, paths, anchor, maxpts, reason)
                target_pts = min(target_pts + step, maxpts)
                continue
            new_anchor = (candidates[0][0], candidates[0][1], anchor[2])
            print(f"[anchor] ({anchor[0]},{anchor[1]}) is DEGENERATE for {off_modal} -- a "
                  f"generically non-zero coefficient vanishes there.")
            print(f"[anchor] re-anchoring on {new_anchor[:2]}, which carries the modal shape, "
                  f"and recomputing the skeleton.")
            anchor = new_anchor
            paths = paths.rebase(anchor) if hasattr(paths, "rebase") else type(paths)(
                pres, target, anchor)
            paths.ensure_dirs()
            samples = SampleSet(paths, atoms=samples._atoms)
            samples.run_anchor()
            basis, proof = samples.anchor_basis(), samples.anchor_proof()
            ring = NumericRing(anchor[0], anchor[1], anchor[2])
            tgt = compose_mod.resolve_target(pres, paths.target, ring)
            steps, remainder = reduce_target(tgt, basis)
            skel = Skeleton.from_reduction(pres, paths.target, anchor, steps, remainder,
                                           basis, proof)
            skel.save(str(paths.skeleton))
            print(f"[anchor] new skeleton: {len(skel.steps)} step(s) on {skel.step_keys()}, "
                  f"remainder {skel.remainder_terms} terms")
            if not skel.in_ideal:
                return _unresolved_remainder(pres, paths, anchor, maxpts, skel)
        try:
            lf = do_lift(skel, samples, atoms=samples.atoms if atoms_path else None)
        except NeedMoreSamples as e:
            print(f"[lift] INCONCLUSIVE: {e}", flush=True)
            if target_pts >= maxpts:
                return _inconclusive(pres, paths, anchor, maxpts, str(e))
            target_pts = min(target_pts + step, maxpts)
            continue
        lf.save(str(paths.lift))
        print(f"[lift] {len(lf.entries)} basis elements lifted from {lf.samples_used} "
              f"stable points -> {paths.lift.name}")
        if lf.atom_report:
            print(f"[atoms] {lf.atoms.names()}  box {lf.atoms.max_exp}", flush=True)

        phase(4, "Symbolic composition over the input relations")
        try:
            mult_samples = _multiplier_samples(skel, samples, pres)
            cof = compose_mod.compose(lf, skel, pres, samples=mult_samples)
        except NeedMoreSamples as e:
            print(f"[compose] INCONCLUSIVE: {e}", flush=True)
            if target_pts >= maxpts:
                return _inconclusive(pres, paths, anchor, maxpts, str(e))
            target_pts = min(target_pts + step, maxpts)
            continue
        cof.save(str(paths.cofactor("symbolic")))
        print(f"[size] {len(cof.terms)} tuples -> {paths.cofactor('symbolic').name}", flush=True)

        phase(5, "Verification: the free-algebra identity over Q(q,t)")
        rep = verify_mod.check(cof, pres, paths.target)
        if rep["pass"]:
            break
        print("\n[verify] Candidate identity failed exact verification; reconstruction is "
              "inconclusive.")
        if target_pts >= maxpts:
            return _inconclusive(pres, paths, anchor, maxpts,
                                 "candidate failed exact verification",
                                 paths.cofactor("symbolic"))
        target_pts = min(target_pts + step, maxpts)

    out = {"status": "verified", "verified": True,
           "presentation": pres.slug(), "presentation_hash": pres.content_hash(),
           "target": paths.target, "anchor": list(anchor), "in_ideal": True,
           "tuples": len(cof.terms), "verify": rep, "points": lf.samples_used,
           "atoms": lf.atoms.names(), "seconds": time.time() - t0,
           "cofactor": str(paths.cofactor("symbolic"))}
    if legacy:
        _legacy_tail(paths)
    print(f"\n{BAR}\nCOMPLETE in {out['seconds']:.1f}s -- {paths.cofactor('symbolic')}\n{BAR}")
    return out


def _multiplier_samples(skel, samples, pres) -> dict:
    """{(q,t): [multiplier, ...]} in skeleton-step order, at compatible stable points.

    The multipliers are not constants: each is a ratio of leading coefficients, so it is a
    rational function to be lifted like any coefficient (WORKFLOW.md notes that their being
    monomials in q alone is a property of this algebra and this target). They are recovered by
    re-reducing the target against each point's own basis.  A sample is usable only when its
    ordered (left word, basis key, right word) signatures match the anchor's.  Basis keys
    identify the elements across changed basis numbering; step positions distinguish repeated
    uses, including repeated uses in the same context.  A different reduction schedule is
    skipped rather than assigning its multipliers to unrelated anchor steps.
    """
    out = {}
    signatures = [(tuple(s["left"]), s["basis_key"], tuple(s["right"]))
                  for s in skel.steps]
    for (q, t) in samples.stable(skel):
        try:
            b = samples.basis(q, t)
            ring = NumericRing(q, t, samples.char)
            tgt = compose_mod.resolve_target(pres, skel.target, ring)
            steps, rem = reduce_target(tgt, b)
        except Exception:
            continue
        if rem.terms:
            continue
        point_signatures = [(tuple(u), b.key[k], tuple(v)) for u, k, v, _mult in steps]
        if point_signatures == signatures:
            out[(q, t)] = [Fraction(mult.v) for _u, _k, _v, mult in steps]
    return out


def _legacy_tail(paths) -> None:
    """Phases 6-10 of the X_3 run: still the original scripts, still B_3-only."""
    scripts = [(6, "Compression: substitute derivable relations", "compress_cofactor.py"),
               (7, "Simplification: Q(q,t)[B_3] via Garside normalization", "simplify_cofactor.py"),
               (8, "Conjugacy reduction to three relations", "reduce_to_three_rels.py"),
               (9, "Relabeling: Z_2, T_2, T_3 convention", "rename_to_ZT.py")]
    for n, title, script in scripts:
        phase(n, title)
        r = subprocess.run([sys.executable, os.path.join(_HERE, script)],
                           capture_output=True, text=True)
        print("\n".join(r.stdout.strip().splitlines()[-10:]))
        if r.returncode != 0:
            print(r.stderr[-600:])
            raise SystemExit(f"phase {n} failed ({script})")
    phase(10, "Final verification: ZT file integrity + ideal-equality logic")
    r = subprocess.run([sys.executable, os.path.join(_HERE, "verify_ZT_and_ideal_equality.py")],
                       capture_output=True, text=True)
    print("\n".join(r.stdout.strip().splitlines()[-25:]))
    if r.returncode != 0:
        raise SystemExit("phase 10 failed")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--untwist", default="bigelow",
                    choices=["bigelow", "square-zero", "silly_question", "none"])
    ap.add_argument("--framing-power", type=int, default=1)
    ap.add_argument("--theta-power", type=int, default=1)
    ap.add_argument("--max-zipper-level", type=int, default=None)
    ap.add_argument("--max-untwist-level", type=int, default=None)
    ap.add_argument("--extra", action="append", default=[])
    ap.add_argument("--target", default="3,1",
                    help="'i,j' for the slider S[i][j], or a name like 'X_3'")
    ap.add_argument("--anchor", default="5,7", help="q,t of the anchor point")
    ap.add_argument("--char", type=int, default=0)
    ap.add_argument("--start", type=int, default=100)
    ap.add_argument("--step", type=int, default=40)
    ap.add_argument("--max", type=int, default=260)
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--atoms", default=None, help="an AtomSet JSON; default is to infer")
    ap.add_argument("--refresh", action="store_true", help="recompute the anchor basis")
    ap.add_argument("--legacy-x3", action="store_true",
                    help="write into the historical X_3 artifact names and run phases 6-10")
    a = ap.parse_args()

    pres = Presentation(a.n, a.untwist, a.framing_power, a.theta_power,
                        a.max_zipper_level, a.max_untwist_level, extra=tuple(a.extra))
    q, t = (int(v) for v in a.anchor.split(","))
    target = tuple(int(v) for v in a.target.split(",")) if "," in a.target else a.target
    try:
        result = run(pres, target, (q, t, a.char), start=a.start, step=a.step, maxpts=a.max,
            legacy=a.legacy_x3, refresh=a.refresh, jobs=a.jobs, atoms_path=a.atoms)
    except NotInIdeal as e:
        print(f"INCONCLUSIVE: {e}; the basis has no completeness certificate")
        return 3
    if result.get("status") == "inconclusive":
        print(f"INCONCLUSIVE: {result['reason']}")
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
