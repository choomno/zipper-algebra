#!/usr/bin/env python3
"""Sample points of a presentation: the anchor, the grid on disk, and how the grid grows.

One f4ncgb run per point (q, t) over the anchor's characteristic gives a basis and a
non-expanded proof.  The lift (lift.py) reads those proofs at every point whose basis
re-identifies the skeleton's elements, and fits each coefficient as a rational function in
Q(q,t).  This module owns the points:

    run_f4ncgb(ms, basis_out, proof_out, flags)   f4ncgb through cache.py, memory-capped
    spiral()                                      (q,t) in PRIMES x PRIMES, spiralling out,
                                                  verbatim from papers/generic/fit_symbolic.py
    admissible(q, t, atoms)                       no atom vanishes, q not in {2, 4}
    SampleSet(paths)                              points on disk, proof/basis loaders,
                                                  stable(skeleton), grow(n_more)

Where the pieces come from.  The runner is reproduce_X3_zero.sh phase 1b (same command line,
same cache key: .ms bytes + flags + the binary's content).  The point loaders are proof.py's
`Basis.load` / `Proof.load`; the stability test is symbolic_lift_x3.load_stable_samples with
basis INDICES replaced by leading-word keys (`Proof.shape`), so that a renumbered basis at some
point is recognised rather than discarded.  grow() is
lift_until_verified.generate (a thread pool of per-point runs) with the points taken from the
spiral instead of that script's fixed candidate list.

Grid exclusions.  A point is unusable if any denominator atom vanishes there: the pipeline's
coefficients are N(q,t)/prod a_i^{e_i}, and a fit through such a point is meaningless.  The
atoms come from the AtomSet (atoms.py), not from a constant.  On top of that the historical
rule q not in {2, 4} stays as a default: run_lift_grid.py recorded a coefficient collapse at
those q (679 / 674 tuples instead of 695) without an atom to explain it, and the 226-point grid
the shipped cofactor was verified on was drawn under it.  Only q = 2 can occur on the prime
spiral; the rule is kept in the form it was documented.

Coefficients are Q(q,t) throughout; a point is an interior node of the lift, never a result.
"""
from __future__ import annotations

import concurrent.futures as cf
import os
import re
import resource
from collections import Counter
import subprocess
import sys
import time
from fractions import Fraction
from itertools import count
from pathlib import Path
from typing import NamedTuple

_HERE = Path(__file__).resolve().parent          # tools/cofactor
_TOOLS = _HERE.parent                            # tools
_REPO = _TOOLS.parent.parent
for _p in (str(_HERE), str(_TOOLS), str(_REPO)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from atoms import AtomSet  # noqa: E402
from cache import Cache  # noqa: E402
from paths import Anchor, Paths  # noqa: E402
from proof import Basis, Proof  # noqa: E402
from repo_settings import settings  # noqa: E402

__all__ = ["Anchor", "DEFAULT_FLAGS", "EXCLUDED_Q", "RunResult", "run_f4ncgb", "f4ncgb_key",
           "spiral", "PRIMES", "admissible", "SampleSet"]

# reproduce_X3_zero.sh phase 1b: F4_ARGS="-e 0 -v 2 -r 0" -- non-expanded proofs, which the
# pipeline re-expands itself (extract_cofactor.py: f4ncgb 1.0.1's expanded proofs agree with
# the basis only modulo the ideal).
DEFAULT_FLAGS = ("-e", "0", "-v", "2", "-r", "0")

# run_lift_grid.py / lift_until_verified.candidate_points: "q in {2,4}: observed coefficient
# collapse" -- see the module docstring.
EXCLUDED_Q = frozenset({2, 4})


# ---------------------------------------------------------------------------
# f4ncgb through the cache
# ---------------------------------------------------------------------------

class RunResult(NamedTuple):
    key: str            # cache key (16 hex chars), as the driver prints it
    hit: bool           # True: artifacts came from the cache, f4ncgb did not run
    seconds: float      # wall time of the f4ncgb run (0.0 on a hit)
    returncode: int     # f4ncgb's exit status (0 on a hit)


def _f4ncgb_binary() -> str:
    return settings().binary("f4ncgb")


def f4ncgb_key(ms_path, flags=DEFAULT_FLAGS, binary: str | None = None,
               cache: Cache | None = None) -> str:
    """The driver's key: `cache.py key --file <ms> --tool <f4ncgb> --tag f4ncgb -- <flags>`.

    NB cache.key hashes the .ms basename together with its bytes, so the same relations under a
    different file name are a different key."""
    return (cache or Cache()).key(files=[str(ms_path)], tag="f4ncgb", args=list(flags),
                                  tools=[binary or _f4ncgb_binary()])


_IDENT_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def check_ms(ms_path) -> None:
    """Refuse input the release f4ncgb would crash on.  A malformed .ms segfaults the release
    build and the WSL kernel writes a 40-60 GB core dump (settings.toml).  The emitter's output
    is grammar-valid; this catches a hand-edited or truncated file before it reaches the binary.
    Checks: a variable line, an integer characteristic, every relation line but the last ending
    in a comma and the last not, no sign runs ('+ -', '- -'), no empty relation."""
    lines = Path(ms_path).read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{ms_path}: fewer than 3 lines; not an msolve input")
    if not all(_IDENT_RE.match(v) for v in lines[0].split(",")):
        raise ValueError(f"{ms_path}: line 1 is not a variable list: {lines[0]!r}")
    if not re.match(r"^\d+$", lines[1].strip()):
        raise ValueError(f"{ms_path}: line 2 is not a characteristic: {lines[1]!r}")
    rels = [l.rstrip() for l in lines[2:] if l.strip()]
    if not rels:
        raise ValueError(f"{ms_path}: no relations")
    for i, l in enumerate(rels):
        last = i == len(rels) - 1
        if last and l.endswith(","):
            raise ValueError(f"{ms_path}: trailing comma on the last relation")
        if not last and not l.endswith(","):
            raise ValueError(f"{ms_path}: relation {i} does not end in a comma")
        body = l[:-1] if l.endswith(",") else l
        if not body.strip():
            raise ValueError(f"{ms_path}: relation {i} is empty")
        if re.search(r"[+-]\s*[+-]", body):
            raise ValueError(f"{ms_path}: sign run in relation {i}: {body[:60]!r}")


def run_f4ncgb(ms_path, basis_out, proof_out, flags=DEFAULT_FLAGS, *, binary: str | None = None,
               cache: Cache | None = None, refresh: bool = False,
               memcap_gb: float | None = None, timeout: float | None = None) -> RunResult:
    """f4ncgb -i <ms> -o <basis_out> -p <proof_out> <flags>, through the cache.

    Exactly reproduce_X3_zero.sh phase 1b: the key is the .ms content + the flags + the binary;
    on a hit the cached basis and proof are copied to the two output paths and f4ncgb does not
    run; on a miss it runs under a memory cap and a timeout ([limits] in settings.toml, or the
    arguments), and the outputs are stored under the key only if it exits 0.  `refresh=True` is
    the driver's REFRESH=1: recompute and overwrite the slot."""
    ms_path, basis_out, proof_out = Path(ms_path), Path(basis_out), Path(proof_out)
    if not ms_path.is_file():
        raise FileNotFoundError(ms_path)
    binary = binary or _f4ncgb_binary()
    cache = cache or Cache()
    key = f4ncgb_key(ms_path, flags, binary, cache)
    if not refresh and cache.get(key, [str(basis_out), str(proof_out)]):
        return RunResult(key, True, 0.0, 0)

    # The safety rule that actually applies to f4ncgb is INPUT VALIDATION, not a memory cap:
    # malformed input makes a release build segfault and the WSL kernel then writes a 40-60 GB
    # core dump to the Windows disk. check_ms is that guard.
    check_ms(ms_path)
    s = settings()
    timeout = float(timeout if timeout is not None else s.limit("timeout_seconds"))

    # No address-space rlimit by default. MEASURED 2026-09-17 on this machine: f4ncgb reserves a
    # fixed arena up front, independent of problem size -- a 14-relation n=3 input aborts with
    # std::bad_alloc/SIGABRT under `ulimit -v 12G` AND under `ulimit -d 12G`, and runs fine at
    # 24 GB or uncapped. So any cap below the reservation kills every run, and any cap above it
    # exceeds this machine's 23 GB of RAM and caps nothing. (It is not mimalloc: ldd shows no
    # alternative allocator.) `ulimit -v` stays correct for Singular, which is the solver that
    # took down WSL uncapped. Pass memcap_gb explicitly to override, knowing the above.
    preexec = None
    if memcap_gb is not None:
        cap = int(float(memcap_gb) * 1024 ** 3)

        def preexec():                                      # noqa: F811
            resource.setrlimit(resource.RLIMIT_AS, (cap, cap))

    basis_out.parent.mkdir(parents=True, exist_ok=True)
    proof_out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "-i", str(ms_path), "-o", str(basis_out), "-p", str(proof_out), *flags]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                           preexec_fn=preexec)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"f4ncgb exceeded {timeout:.0f} s on {ms_path.name}") from None
    secs = time.time() - t0
    if r.returncode != 0 or not (basis_out.is_file() and proof_out.is_file()):
        tail = (r.stdout[-800:] + r.stderr[-800:]).strip()
        raise RuntimeError(f"f4ncgb rc={r.returncode} on {ms_path.name} after {secs:.1f} s "
                           f"(timeout {timeout:.0f} s)\n{tail}")
    cache.put(key, [str(basis_out), str(proof_out)])
    return RunResult(key, False, secs, r.returncode)


# ---------------------------------------------------------------------------
# sampling: a spiral through Z x Z  (papers/generic/fit_symbolic.py, verbatim)
# ---------------------------------------------------------------------------

def primes(limit=100000):
    """Primes, ascending."""
    sieve = bytearray([1]) * limit
    sieve[0] = sieve[1] = 0
    for i in range(2, int(limit ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i, b in enumerate(sieve) if b]


PRIMES = primes()


def spiral():
    """(q,t) drawn from PRIMES x PRIMES, spiralling out from the small primes,
    with the diagonal removed.

    Primes rather than all of Z: the theory's degeneracy loci are small algebraic
    conditions on (q,t) -- the denominator atoms of the X_3 cofactor are
    {q, t, q-1, q+1, t-q, t+q^2} -- and prime coordinates avoid the accidental
    coincidences that small composites fall into. q = 0, 1 are excluded outright
    (q and 1-q vanish there).

    The diagonal q = t is dropped: t - q is itself one of those atoms, so the
    diagonal is a degeneracy locus of the pipeline.
    """
    for (a, b) in _spiral_indices():
        if a == b:                      # diagonal q = t: t - q vanishes
            continue
        yield (PRIMES[a], PRIMES[b])


def _spiral_indices():
    """Non-negative lattice points, spiralling out from (0,0)."""
    x = y = 0
    dx, dy = 1, 0
    for leg in count(1):
        for _ in range(2):
            for _ in range(leg):
                x, y = x + dx, y + dy
                if x >= 0 and y >= 0:
                    yield (x, y)
            dx, dy = -dy, dx
        if leg > 200:
            raise RuntimeError("spiral exhausted without a fit")


# ---------------------------------------------------------------------------
# grid exclusions
# ---------------------------------------------------------------------------

def admissible(q: int, t: int, atoms: AtomSet, excluded_q=EXCLUDED_Q) -> bool:
    """False where a denominator atom vanishes, or q is in the documented exclusion set."""
    if q in excluded_q:
        return False
    return all(a.eval(q, t) != 0 for a in atoms.atoms)


# ---------------------------------------------------------------------------
# proof / basis loaders (proof.py's types; the variable names come from the presentation)
# ---------------------------------------------------------------------------

def load_proof(path, pres) -> Proof:
    return Proof.load(str(path), pres.var_names())


def load_basis(path, pres, char: int) -> Basis:
    return Basis.load(str(path), char, pres.var_names())


# ---------------------------------------------------------------------------
# SampleSet
# ---------------------------------------------------------------------------

class SampleSet:
    """The points of one presentation at the anchor's characteristic, over a Paths.

        S = SampleSet(Paths.legacy_x3())
        S.points                      -> the 100 (q,t) with a proof on disk
        S.run_anchor()                -> phase 1a + 1b (emit .ms, f4ncgb through the cache)
        S.proof(q, t) / S.basis(q, t) -> Proof / Basis at a point (proof.py's types)
        S.stable(skeleton)            -> points whose proofs re-identify every skeleton key
        S.next_points(n)              -> the next n admissible spiral points not on disk
        S.grow(n)                     -> run f4ncgb (cached) at those, return them
    """

    def __init__(self, paths: Paths, atoms: AtomSet | None = None, flags=DEFAULT_FLAGS):
        self.paths = paths
        self.pres = paths.pres
        self.anchor: Anchor = paths.anchor
        self.flags = tuple(flags)
        self._atoms = atoms
        self._census_memo = (None, None)
        self._proofs: dict = {}
        self._bases: dict = {}

    # --- atoms -------------------------------------------------------------
    @property
    def atoms(self) -> AtomSet:
        """The atom file for this presentation if written, else the historical six as a prior
        (the same rule as symbolic_lift_x3.py)."""
        if self._atoms is None:
            p = self.paths.atoms
            self._atoms = AtomSet.load(str(p)) if p.is_file() else AtomSet.DEFAULT
        return self._atoms

    # --- points on disk ----------------------------------------------------------
    @property
    def points(self) -> list[tuple[int, int]]:
        return self.paths.points_on_disk()

    @property
    def char(self) -> int:
        return self.anchor[2]

    # --- loaders ---------------------------------------------------------------------
    def proof(self, q: int, t: int):
        if (q, t) not in self._proofs:
            self._proofs[(q, t)] = load_proof(self.paths.proof(q, t), self.pres)
        return self._proofs[(q, t)]

    def basis(self, q: int, t: int):
        if (q, t) not in self._bases:
            self._bases[(q, t)] = load_basis(self.paths.basis(q, t), self.pres, self.char)
        return self._bases[(q, t)]

    def anchor_proof(self):
        if "anchor" not in self._proofs:
            self._proofs["anchor"] = load_proof(self.paths.anchor_proof, self.pres)
        return self._proofs["anchor"]

    def anchor_basis(self):
        if "anchor" not in self._bases:
            self._bases["anchor"] = load_basis(self.paths.anchor_basis, self.pres, self.char)
        return self._bases["anchor"]

    # --- running f4ncgb --------------------------------------------------------------
    def emit(self, q: int, t: int) -> Path:
        """The presentation's msolve input at (q, t), byte-identical to zipper_in_ms.py's CLI."""
        path = self.paths.anchor_ms if (q, t) == self.anchor[:2] else self.paths.ms(q, t)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.pres.emit_ms(str(path), q, t, self.char)
        return path

    def run_anchor(self, refresh: bool = False, **kw) -> RunResult:
        """reproduce_X3_zero.sh phases 1a + 1b."""
        self.paths.ensure_dirs()
        ms = self.emit(*self.anchor[:2])
        return run_f4ncgb(ms, self.paths.anchor_basis, self.paths.anchor_proof, self.flags,
                          refresh=refresh, **kw)

    def run_point(self, q: int, t: int, refresh: bool = False, **kw) -> RunResult:
        self.paths.ensure_dirs()
        ms = self.emit(q, t)
        return run_f4ncgb(ms, self.paths.basis(q, t), self.paths.proof(q, t), self.flags,
                          refresh=refresh, **kw)

    def next_points(self, n_more: int) -> list[tuple[int, int]]:
        """The next n_more spiral points that are admissible and not yet on disk, in spiral
        order.  Pure selection: nothing is run or written."""
        have = set(self.points)
        atoms = self.atoms
        out: list[tuple[int, int]] = []
        for (q, t) in spiral():
            if len(out) >= n_more:
                break
            if (q, t) in have or (q, t) in out or not admissible(q, t, atoms):
                continue
            out.append((q, t))
        return out

    def grow(self, n_more: int, jobs: int | None = None, refresh: bool = False,
             **kw) -> list[tuple[int, int]]:
        """Run f4ncgb (through the cache) at next_points(n_more), `jobs` at a time.

        Threads, as in lift_until_verified.generate: every task is a subprocess.  The memory
        cap is per f4ncgb process, so jobs x memcap_gb is the worst case; `jobs` defaults to
        [limits].jobs.  Every run must succeed; the first failure raises with f4ncgb's tail."""
        pts = self.next_points(n_more)
        if not pts:
            return []
        self.paths.ensure_dirs()
        jobs = int(jobs or settings().limit("jobs"))

        def one(qt):
            return qt, self.run_point(*qt, refresh=refresh, **kw)

        with cf.ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
            results = list(ex.map(one, pts))
        return [qt for qt, _ in results]

    # --- stability ---------------------------------------------------------------------
    def shape_census(self, skeleton) -> dict:
        """{key: Counter(shape -> points)} over every point on disk, plus the anchor's shape.

        Why this exists: a proof's term structure is what must agree across points for their
        coefficients to be specializations of the same rational functions, and the obvious
        reference is the anchor's.  But the anchor can be DEGENERATE -- a coefficient that is
        generically non-zero can vanish there, collapsing a term, and then the anchor's shape is
        the exception rather than the rule and almost every point is rejected.  Measured
        2026-09-17 on n=3 square-zero: element `x2*x1*y2*x1` has 109 proof terms at (5,7) and 110
        at every other point, so 2 of 160 points were "stable".  The reference must therefore be
        the MODAL shape, and a non-modal anchor is a fact to report, not to absorb.
        """
        keys = getattr(skeleton, "basis_keys", None)
        if keys is None:
            keys = skeleton["basis_keys"]
        keys = list(keys)
        # Memoised on (keys, points): the census is asked for three times in one run -- once by
        # the driver to test the anchor, twice more inside stable() via modal_shapes -- and each
        # pass is O(points x keys) shape constructions over every parsed proof.
        memo_key = (tuple(keys), tuple(self.points), self.paths.anchor)
        if getattr(self, "_census_memo", (None,))[0] == memo_key:
            return self._census_memo[1]
        ab, ap = self.anchor_basis(), self.anchor_proof()
        missing = [k for k in keys if k not in ab.index_by_key]
        if missing:
            raise KeyError(f"skeleton keys not in the anchor basis: {missing[:5]}")
        anchor_shape = {k: ap.shape(ab.index_by_key[k], ab) for k in keys}
        census = {k: Counter() for k in keys}
        per_point = {}
        for (q, t) in self.points:
            b, p = self.basis(q, t), self.proof(q, t)
            shapes = {}
            for k in keys:
                i = b.index_by_key.get(k)
                shapes[k] = p.shape(i, b) if (i is not None and i < len(p.terms)) else None
                census[k][shapes[k]] += 1
            per_point[(q, t)] = shapes
        out = {"keys": keys, "census": census, "anchor_shape": anchor_shape,
               "per_point": per_point}
        self._census_memo = (memo_key, out)
        return out

    def modal_shapes(self, skeleton) -> tuple[dict, list]:
        """(modal shape per key, keys where the anchor is not modal)."""
        c = self.shape_census(skeleton)
        modal, off = {}, []
        for k in c["keys"]:
            counts = c["census"][k]
            best = max(counts.items(), key=lambda kv: (kv[1], kv[0] is not None))[0] \
                if counts else c["anchor_shape"][k]
            modal[k] = best if best is not None else c["anchor_shape"][k]
            if modal[k] != c["anchor_shape"][k]:
                off.append(k)
        return modal, off

    def stable(self, skeleton, reference="modal") -> list[tuple[int, int]]:
        """Points whose proofs re-identify every key of `skeleton.basis_keys` with the reference
        term structure (`Proof.shape`: the (left, kind, key-or-input-idx, right) sequence).

        reference="modal" (default) compares against the most common shape per key, so a
        degenerate anchor costs one point rather than all of them; reference="anchor" is the
        historical criterion (symbolic_lift_x3.skel with basis indices replaced by keys).  The
        two agree whenever the anchor is generic, which it is for n=4 bigelow at (5,7): 99 of
        100 points either way.  Order is the order of `self.points`.
        """
        c = self.shape_census(skeleton)
        keys, per_point = c["keys"], c["per_point"]
        want = c["anchor_shape"] if reference == "anchor" else self.modal_shapes(skeleton)[0]
        return [pt for pt in self.points
                if all(per_point[pt][k] == want[k] for k in keys)]

    def __repr__(self) -> str:
        return f"SampleSet({self.paths!r}, {len(self.points)} points on disk)"


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="inspect the sample grid of the legacy X_3 layout")
    ap.add_argument("--next", type=int, default=10, help="show the next N admissible spiral points")
    a = ap.parse_args()
    S = SampleSet(Paths.legacy_x3())
    print(S)
    print("atoms:", S.atoms)
    print(f"next {a.next} points not on disk:", S.next_points(a.next))
