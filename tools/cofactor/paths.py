#!/usr/bin/env python3
"""Every artifact path of the cofactor pipeline, from (presentation, target, anchor).

Before this module each phase spelled its own file names: the driver, `extract_cofactor.py`,
`extract_cofactor_qt.py`, `symbolic_lift_x3.py`, `compose_x3_cofactor_qt.py` and
`lift_until_verified.py` each carried `logs_and_traces/f4ncgb/` and `n4`/`Q`/`X3` baked into
strings, and agreed only because nobody changed one of them.  Here the names are derived once
from the three things they depend on:

    pres    Presentation      -> the directory (pres.slug())
    target  "S[3][1]" etc.    -> skeleton / lift / cofactor names
    anchor  (q, t, char)      -> the anchor's .ms / basis / proof

New layout (default root = settings().path("output") / "cofactor"):

    <root>/<pres.slug()>/                         char 0; "<slug>-c<char>/" for char p
        input-<q>-<t>-c<char>.ms                  ONE name per point, the anchor included (*)
        basis-<q>-<t>.txt                         f4ncgb basis at the point
        proof-<q>-<t>-noex.txt                    f4ncgb non-expanded proof at the point
        skeleton-<target>.json                    phase 2
        atoms.json, atoms-report.json             AtomSet inferred for this presentation
        lift-<target>.json, lift-<target>.txt     phase 3
        cofactor-<target>-<phase>.txt             phase 4 ("symbolic") and later phases
        cofactor-<target>-<q>-<t>.txt             numeric cofactor at one point (a by-product)

(*) The SPEC sketched `anchor-<q>-<t>-c<char>.ms` for the anchor.  `cache.Cache.key` hashes the
.ms BASENAME as well as its bytes, so if the anchor and the grid spelled the same point
differently the anchor's f4ncgb run would be repeated as a grid run and the two copies of the
basis would live under different names, as they do today (basis-n4-Q.txt == basis-n4-Q-q5-t7.txt).
One name per point removes that; `anchor_ms` is simply `ms(*anchor)`.

Compatibility: `Paths.legacy_x3()` returns the names the X_3 pipeline uses today under
`zipper-algebra/logs_and_traces/f4ncgb/` (and `ms_inputs/qt_sweep/` for the per-point inputs),
so the byte-identity test keeps its artifact locations until the switch-over.  Every legacy
spelling lives in `_LegacyX3Paths`, read off the scripts named above; `tests_paths_samples.py`
checks them against those scripts' own constants and against the files on disk.

Nothing here creates a directory unless `ensure_dirs()` is called.
"""
from __future__ import annotations

import glob
import re
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent          # tools/cofactor
_TOOLS = _HERE.parent                            # tools
_ZA = _TOOLS.parent                              # zipper-algebra
_REPO = _ZA.parent
if __package__:
    from ..presentation import Presentation
else:
    for _p in (str(_HERE), str(_TOOLS), str(_REPO)):
        if _p not in sys.path:
            sys.path.insert(0, _p)
    from presentation import Presentation


Anchor = tuple[int, int, int]            # (q, t, char); char 0 = Q

_SLIDER_RE = re.compile(r"^S\[(\d+)\]\[(\d+)\]$")


def normalise_target(target) -> str:
    """"S[3][1]" as written by extract_cofactor's skeleton; a pair (i, j) or the driver's
    "i,j" spelling (reproduce_X3_zero.sh: TARGET=3,1) both mean S[i][j]."""
    if isinstance(target, (tuple, list)) and len(target) == 2:
        return f"S[{int(target[0])}][{int(target[1])}]"
    m = re.match(r"^\s*(\d+)\s*,\s*(\d+)\s*$", str(target))
    if m:
        return f"S[{m.group(1)}][{m.group(2)}]"
    return str(target)


def slider_indices(target: str) -> tuple[int, int] | None:
    m = _SLIDER_RE.match(target)
    return (int(m.group(1)), int(m.group(2))) if m else None


def target_slug(target: str) -> str:
    """File-name form of a target: "S[3][1]" -> "S3_1"; any other label keeps its alphanumerics
    with runs of everything else collapsed to "_" ("(Theta_3-q^2)X_3" -> "Theta_3_q_2_X_3")."""
    ij = slider_indices(target)
    if ij:
        return f"S{ij[0]}_{ij[1]}"
    s = re.sub(r"[^A-Za-z0-9]+", "_", target).strip("_")
    return s or "target"


def default_root() -> Path:
    from repo_settings import settings
    return settings().path("output") / "cofactor"


class Paths:
    """Artifact names for one (presentation, target, anchor).  See the module docstring."""

    legacy = False

    def __init__(self, pres: Presentation, target, anchor: Anchor, root=None):
        self.pres = pres
        self.target = normalise_target(target)
        self.anchor: Anchor = tuple(int(v) for v in anchor)
        if len(self.anchor) != 3:
            raise ValueError(f"anchor must be (q, t, char), got {anchor!r}")
        self.root = Path(root) if root is not None else default_root()
        self.dir = self.root / (pres.slug() + ("" if self.char == 0 else f"-c{self.char}"))
        self.tslug = target_slug(self.target)
        # how the sample set finds its points: every proof file in dir matching the glob,
        # the (q, t) read off with the regex
        self.proof_glob = "proof-*-*-noex.txt"
        self.proof_re = re.compile(r"proof-(\d+)-(\d+)-noex\.txt$")

    # --- the anchor ---------------------------------------------------------
    @property
    def q(self) -> int:
        return self.anchor[0]

    @property
    def t(self) -> int:
        return self.anchor[1]

    @property
    def char(self) -> int:
        return self.anchor[2]

    # --- per-point artifacts ------------------------------------------------
    def ms(self, q: int, t: int) -> Path:
        return self.dir / f"input-{q}-{t}-c{self.char}.ms"

    def basis(self, q: int, t: int) -> Path:
        return self.dir / f"basis-{q}-{t}.txt"

    def proof(self, q: int, t: int) -> Path:
        return self.dir / f"proof-{q}-{t}-noex.txt"

    def cofactor_at(self, q: int, t: int) -> Path:
        """Numeric cofactor of the target at one point -- what extract_cofactor_qt.py wrote per
        grid point.  No phase reads it; kept because the lift's consistency probes did."""
        return self.dir / f"cofactor-{self.tslug}-{q}-{t}.txt"

    @property
    def anchor_ms(self) -> Path:
        return self.ms(self.q, self.t)

    @property
    def anchor_basis(self) -> Path:
        return self.basis(self.q, self.t)

    @property
    def anchor_proof(self) -> Path:
        return self.proof(self.q, self.t)

    @property
    def cofactor_anchor(self) -> Path:
        """extract_cofactor.py's numeric cofactor at the anchor (phase 2's by-product)."""
        return self.cofactor_at(self.q, self.t)

    # --- per-(presentation, target) artifacts --------------------------------
    @property
    def skeleton(self) -> Path:
        return self.dir / f"skeleton-{self.tslug}.json"

    @property
    def atoms(self) -> Path:
        return self.dir / "atoms.json"

    @property
    def atoms_report(self) -> Path:
        return self.dir / "atoms-report.json"

    @property
    def lift(self) -> Path:
        return self.dir / f"lift-{self.tslug}.json"

    @property
    def lift_txt(self) -> Path:
        return self.dir / f"lift-{self.tslug}.txt"

    def cofactor(self, phase: str = "symbolic") -> Path:
        """phase 4 writes "symbolic"; phases 6-9 "compressed", "simplified", "three-rels", "ZT"."""
        return self.dir / f"cofactor-{self.tslug}-{phase}.txt"

    # --- discovery -------------------------------------------------------------
    def points_on_disk(self) -> list[tuple[int, int]]:
        """(q, t) of every proof file in dir, in lexicographic order of the FILE NAME.

        That is the order `symbolic_lift_x3.load_stable_samples` iterated (`sorted(glob(...))`),
        and the lift's hold-out split takes the last points of the list, so the order is part of
        what the lift computes; it is kept rather than replaced by numeric order."""
        out = []
        for f in sorted(glob.glob(str(self.dir / self.proof_glob))):
            m = self.proof_re.search(f)
            if m:
                out.append((int(m.group(1)), int(m.group(2))))
        return out

    def ensure_dirs(self) -> None:
        self.dir.mkdir(parents=True, exist_ok=True)

    # --- the X_3 pipeline's names ------------------------------------------------
    @classmethod
    def legacy_x3(cls, pres: Presentation | None = None, target="S[3][1]",
                  anchor: Anchor = (5, 7, 0)) -> "Paths":
        return _LegacyX3Paths(pres or Presentation.DEFAULT, target, anchor)

    def __repr__(self) -> str:
        return (f"{type(self).__name__}({self.pres.slug()}, {self.target!r}, "
                f"anchor={self.anchor}, dir={self.dir})")


class _LegacyX3Paths(Paths):
    """Today's spellings, one method each, with the script they come from.

    Only char 0 ("Q" in every name) was ever run, and only target S[3][1] ("X3") through phases
    4-9; the S{i}_{j} rule for other sliders is extract_cofactor.py's own (phase 2 only)."""

    legacy = True

    def __init__(self, pres: Presentation, target, anchor: Anchor):
        super().__init__(pres, target, anchor, root=_ZA / "logs_and_traces")
        if self.char != 0:
            raise ValueError("the legacy X_3 layout is char 0 only ('Q' in every name)")
        self.dir = _ZA / "logs_and_traces" / "f4ncgb"            # reproduce_X3_zero.sh: LOG_DIR
        self.ms_dir = _ZA / "ms_inputs" / "qt_sweep"              # extract_cofactor_qt.py
        self.transcript = _ZA / "logs_and_traces" / "reproduce-X3-zero-transcript.txt"
        n = pres.n
        ij = slider_indices(self.target)
        # extract_cofactor.py main(): tag = f"n{n}" + ("" if (ti,tj)==(3,1) else f"-S{ti}_{tj}")
        self._n = n
        self._tag = f"n{n}" + ("" if ij == (3, 1) else f"-S{ij[0]}_{ij[1]}" if ij
                               else f"-{self.tslug}")
        self._x = "X3" if ij == (3, 1) else (f"S{ij[0]}_{ij[1]}" if ij else self.tslug)
        # lift_until_verified.py: PROOF_RE and the glob of have_points()
        self.proof_glob = f"proof-n{n}-Q-q*-t*-noex.txt"
        self.proof_re = re.compile(rf"proof-n{n}-Q-q(\d+)-t(\d+)-noex\.txt$")

    # per point: extract_cofactor_qt.py main()
    def ms(self, q, t):
        return self.ms_dir / f"zipper-n{self._n}-Q-q{q}-t{t}.ms"

    def basis(self, q, t):
        return self.dir / f"basis-n{self._n}-Q-q{q}-t{t}.txt"

    def proof(self, q, t):
        return self.dir / f"proof-n{self._n}-Q-q{q}-t{t}-noex.txt"

    def cofactor_at(self, q, t):
        return self.dir / f"cofactor-{self._x}-n{self._n}-Q-q{q}-t{t}.txt"

    # the anchor: reproduce_X3_zero.sh phase 1 (INPUT_MS, BASIS_OUT, PROOF_OUT) and
    # extract_cofactor.py main() (basis_path, noex_proof_path, out_path)
    @property
    def anchor_ms(self):
        return self.dir / f"zipper-n{self._n}-Q-q{self.q}-t{self.t}.ms"

    @property
    def anchor_basis(self):
        return self.dir / f"basis-n{self._n}-Q.txt"

    @property
    def anchor_proof(self):
        return self.dir / f"proof-n{self._n}-Q-noex.txt"

    @property
    def cofactor_anchor(self):
        return self.dir / f"cofactor-{self._x}-n{self._n}-Q.txt"

    # per (presentation, target)
    @property
    def skeleton(self):                       # extract_cofactor.py: f"skeleton-{tag}-Q.json"
        return self.dir / f"skeleton-{self._tag}-Q.json"

    @property
    def atoms(self):                          # symbolic_lift_x3.py: ATOMS_PATH
        return self.dir / f"atoms-n{self._n}-Q.json"

    @property
    def atoms_report(self):                   # symbolic_lift_x3.py: ATOMS_PATH[:-5] + '-report.json'
        return self.dir / f"atoms-n{self._n}-Q-report.json"

    @property
    def lift(self):                           # symbolic_lift_x3.py out_json; compose: JSON_PATH
        return self.dir / "basis-symbolic-Qqt.json"

    @property
    def lift_txt(self):                       # symbolic_lift_x3.py out_path
        return self.dir / "basis-symbolic-Qqt.txt"

    def cofactor(self, phase="symbolic"):     # compose_x3_cofactor_qt.py out_path; phases 6-9
        return self.dir / f"cofactor-{self._x}-n{self._n}-Qqt-{phase}.txt"

    def ensure_dirs(self):
        # The legacy directories exist and are never written by the new modules (SPEC rule).
        if not self.dir.is_dir():
            raise FileNotFoundError(f"legacy artifact directory missing: {self.dir}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="print the artifact paths for a presentation")
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--untwist", default="bigelow")
    ap.add_argument("--target", default="S[3][1]")
    ap.add_argument("--anchor", default="5,7,0", metavar="q,t,char")
    ap.add_argument("--legacy-x3", action="store_true")
    a = ap.parse_args()
    pres = Presentation(a.n, a.untwist)
    anchor = tuple(int(v) for v in a.anchor.split(","))
    P = Paths.legacy_x3(pres, a.target, anchor) if a.legacy_x3 else Paths(pres, a.target, anchor)
    print(P)
    for name in ("anchor_ms", "anchor_basis", "anchor_proof", "skeleton", "atoms", "lift",
                 "cofactor_anchor"):
        print(f"  {name:16s} {getattr(P, name)}")
    print(f"  {'cofactor':16s} {P.cofactor()}")
    print(f"  {'ms(q,t)':16s} {P.ms(11, 13)}")
    print(f"  {'points_on_disk':16s} {len(P.points_on_disk())}")
