#!/usr/bin/env python3
"""Cofactor: the membership certificate, in ONE text format with ONE reader and ONE writer.

A cofactor of a target f in the two-sided ideal of a presentation is the finite sum

    f = sum_j c_j(q,t) * L_j * R_{i_j} * R'_j

with L_j, R'_j words in the generators, R_i the input relations of the presentation, and
c_j in Q(q,t) written as num / den, where den is a product of powers of the presentation's
denominator atoms (tools/cofactor/atoms.py) and num is a polynomial in q, t.  One term per
text line, exactly the line compose_x3_cofactor_qt.py has always written:

    x1*x2 * T_2 * y1  ::  (1*q^2*t^0 + -1*q^3*t^1) / (q^2*t*(q-1))

The empty word is "()"; the numerator lists its monomials in ascending (deg_q, deg_t) with
every coefficient and both exponents written out; the denominator is AtomSet.denom_str of
the exponent tuple.  The header is '#' comment lines, one key per line:

    # presentation: n4-bigelow                 Presentation.slug()
    # presentation_hash: 0123456789abcdef      Presentation.content_hash()
    # target: S[3][1]                          what the identity proves membership of
    # atoms: {"atoms": [...], ...}             AtomSet.to_json(), on one line
    # n_terms: 695
    # producer: tools/cofactor/compose.py
    # anchor: 5,7,0                            only for a cofactor specialised at a point
    # note: ...                                free text carried over from a legacy file

so a cofactor file says which algebra it is over, what it proves, and over which atoms its
denominators are written.  Until 2026-09-17 there were three writers (extract_cofactor.py
with its "+c * L * [idx] * R" lines, compose_x3_cofactor_qt.py, and the
compress/simplify/three-rels family) and none recorded the atoms or the presentation; `load`
still reads all of them, given the presentation the caller knows they were written over.

Words are tuples of 1-based variable indices in the presentation's variable order
(x1..x_{n-1}, y1..y_{n-1}), as in zipper_in_ms.Poly.  Relation labels are the strings
Presentation.labels() derives from the construction (INV_1_p, BR_1, Z_2, T_3, ...);
`relation_indices(pres)` maps them back to positions in the relation list, and the legacy
"R_<idx>" form is accepted too.

A cofactor specialised at a point (the extract format) is stored with constant numerators
and a trivial denominator; it is an interior node of the reconstruction, never a result.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py and presentation.py live
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import json  # noqa: E402
import re  # noqa: E402
from fractions import Fraction  # noqa: E402

from atoms import AtomSet  # noqa: E402
from zipper_in_ms import Coef, Poly, NumericRing  # noqa: E402
from presentation import Presentation, EXTRA_CATALOGUE, _theta_word  # noqa: E402

HEADER_KEYS = ("presentation", "presentation_hash", "target", "atoms", "n_terms", "producer")
_HEADER_RE = re.compile(r"^#\s*([a-z_]+):\s*(.*)$")
_LEGACY_ANCHOR_RE = re.compile(r"\(q,\s*t\)\s*=\s*\((-?\d+),\s*(-?\d+)\)")
_S_TARGET_RE = re.compile(r"^S\[(\d+)\]\[(\d+)\]$")
_X_TARGET_RE = re.compile(r"^X_(\d+)$")
_R_LABEL_RE = re.compile(r"^R_(\d+)$")


# ---------------------------------------------------------------------------
# Words, numerators, denominators: text <-> data
# ---------------------------------------------------------------------------

def var_names(n: int) -> list[str]:
    """The presentation's variable order: x1..x_{n-1}, y1..y_{n-1}."""
    return [f"x{i}" for i in range(1, n)] + [f"y{i}" for i in range(1, n)]


def word_str(word: tuple, names: list[str]) -> str:
    return "*".join(names[i - 1] for i in word) if word else "()"


def parse_word(s: str, var_index: dict) -> tuple:
    s = s.strip()
    if s in ("", "()"):
        return ()
    try:
        return tuple(var_index[v.strip()] for v in s.split("*"))
    except KeyError as e:
        raise ValueError(f"unknown variable {e.args[0]!r} in word {s!r}") from None


def word_from_json(w, var_index: dict) -> tuple:
    """A word from JSON: a list of 1-based indices, or a list of variable names."""
    if not w:
        return ()
    if all(isinstance(v, int) for v in w):
        return tuple(w)
    return tuple(var_index[v] for v in w)


def num_str(num: dict) -> str:
    """{(dq, dt): Fraction} -> '1*q^2*t^0 + -1*q^3*t^1' (ascending monomials); '0' if empty."""
    return " + ".join(f"{c}*q^{i}*t^{j}" for (i, j), c in sorted(num.items())) if num else "0"


def _strip_outer(s: str) -> str:
    """Remove one pair of enclosing parentheses if they enclose the whole string."""
    s = s.strip()
    if not (s.startswith("(") and s.endswith(")")):
        return s
    depth = 0
    for k, ch in enumerate(s):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0 and k != len(s) - 1:
                return s
    return s[1:-1].strip()


def parse_num(s: str) -> dict:
    """Inverse of num_str.  Also accepts the looser forms other writers produced: bare 'q^2',
    '-t', '3', fractions '1/2*q^1*t^0', and ' - ' between terms."""
    s = _strip_outer(s)
    out: dict = {}
    for term in s.replace(" - ", " + -").split(" + "):
        term = term.strip()
        if not term:
            continue
        sign = 1
        if term[0] in "+-":
            sign = -1 if term[0] == "-" else 1
            term = term[1:].strip()
        coef, i, j = Fraction(1), 0, 0
        for piece in term.split("*"):
            piece = piece.strip()
            if not piece:
                continue
            if piece[0].isdigit():
                coef *= Fraction(piece)
            elif piece[0] in "qt":
                e = int(piece[2:]) if len(piece) > 1 else 1
                if len(piece) > 1 and piece[1] != "^":
                    raise ValueError(f"cannot parse numerator factor {piece!r} in {s!r}")
                if piece[0] == "q":
                    i += e
                else:
                    j += e
            else:
                raise ValueError(f"cannot parse numerator factor {piece!r} in {s!r}")
        c = sign * coef
        if c == 0:
            continue
        out[(i, j)] = out.get((i, j), Fraction(0)) + c
        if out[(i, j)] == 0:
            del out[(i, j)]
    return out


# ---------------------------------------------------------------------------
# Targets: the string a skeleton / cofactor header names, resolved over a ring
# ---------------------------------------------------------------------------

def resolve_target(pres: Presentation, target, ring) -> Poly:
    """The target polynomial over `ring` (NumericRing or SymbolicRing).

    `target` is a Poly (returned as is), a pair (i, j) meaning the slider S[i][j], or a
    string: "S[i][j]", "X_k" (= S[k][1] by definition), or the name of an EXTRA relation
    from presentation.EXTRA_CATALOGUE (optionally written "EXTRA[name]")."""
    if isinstance(target, Poly):
        return target
    if isinstance(target, tuple) and len(target) == 2:
        return pres.target(target[0], target[1], ring)
    if not isinstance(target, str):
        raise TypeError(f"target must be a Poly, an (i, j) pair or a string, not {type(target).__name__}")
    s = target.replace(" ", "")
    m = _S_TARGET_RE.match(s)
    if m:
        return pres.target(int(m.group(1)), int(m.group(2)), ring)
    m = _X_TARGET_RE.match(s)
    if m:
        return pres.target(int(m.group(1)), 1, ring)
    name = s[6:-1] if s.startswith("EXTRA[") and s.endswith("]") else target
    if name in EXTRA_CATALOGUE:
        _, _, S, _ = pres.build(ring)
        return EXTRA_CATALOGUE[name](S, ring, lambda i, j: _theta_word(i, j, ring, pres.n))
    raise ValueError(f"cannot resolve target {target!r}: not S[i][j], X_k or an EXTRA name "
                     f"{sorted(EXTRA_CATALOGUE)}")


def target_name(target: str) -> str:
    """Display name: S[k][1] is X_k throughout this project; anything else is itself."""
    m = _S_TARGET_RE.match(str(target).replace(" ", ""))
    if m and m.group(2) == "1":
        return f"X_{m.group(1)}"
    return str(target)


# ---------------------------------------------------------------------------
# Cofactor
# ---------------------------------------------------------------------------

class Cofactor:
    """terms: list of (left_word, relation_label, right_word, num, den_exps).

    left_word / right_word: tuples of 1-based variable indices; relation_label: str;
    num: {(deg_q, deg_t): Fraction}; den_exps: exponent tuple over `atoms`."""

    def __init__(self, terms, n: int, atoms: AtomSet, presentation: str = "",
                 presentation_hash: str = "", target: str = "", producer: str = "",
                 anchor=None, notes=()):
        self.terms = list(terms)
        self.n = int(n)
        self.atoms = atoms
        self.presentation = presentation
        self.presentation_hash = presentation_hash
        self.target = target
        self.producer = producer
        self.anchor = tuple(anchor) if anchor is not None else None
        self.notes = list(notes)

    @classmethod
    def from_presentation(cls, pres: Presentation, target: str, atoms: AtomSet, terms,
                          producer: str, anchor=None, notes=()) -> "Cofactor":
        return cls(terms, pres.n, atoms, pres.slug(), pres.content_hash(), str(target),
                   producer, anchor, notes)

    # --- views -----------------------------------------------------------------
    @property
    def header(self) -> dict:
        h = {"presentation": self.presentation, "presentation_hash": self.presentation_hash,
             "target": self.target, "atoms": self.atoms.to_json(), "n_terms": len(self.terms),
             "producer": self.producer}
        if self.anchor is not None:
            h["anchor"] = self.anchor
        if self.notes:
            h["notes"] = list(self.notes)
        return h

    @property
    def var_names(self) -> list[str]:
        return var_names(self.n)

    @property
    def n_terms(self) -> int:
        return len(self.terms)

    def labels(self) -> list[str]:
        return sorted({t[1] for t in self.terms})

    def relation_indices(self, pres: Presentation) -> list[int]:
        """The 0-based relation index of every term, in term order."""
        by_label = {lab: i for i, lab in enumerate(pres.labels())}
        out = []
        for _, lab, _, _, _ in self.terms:
            if lab in by_label:
                out.append(by_label[lab])
                continue
            m = _R_LABEL_RE.match(lab)
            if m:
                out.append(int(m.group(1)))
                continue
            raise ValueError(f"relation label {lab!r} is not one of {pres.slug()}'s "
                             f"{len(by_label)} labels")
        return out

    def lcm_denominator(self) -> tuple:
        """Componentwise max of the denominator exponents: the M that clears every term."""
        exps = self.atoms.zero()
        for _, _, _, _, den in self.terms:
            exps = AtomSet.lcm(exps, den)
        return exps

    def evaluate(self, q: int, t: int, pres: Presentation) -> Poly:
        """sum_j c_j(q,t) * L_j * R_{i_j} * R'_j over Q at the integer point (q, t).
        Raises ZeroDivisionError if a denominator atom vanishes there."""
        ring = NumericRing(q, t, 0)
        _, relations, _, _ = pres.build(ring)
        idxs = self.relation_indices(pres)
        acc = Poly(0)
        for (left, _, right, num, den), idx in zip(self.terms, idxs):
            d_val = self.atoms.eval_denom(den, q, t)
            if d_val == 0:
                raise ZeroDivisionError(f"a denominator atom vanishes at (q,t)=({q},{t})")
            n_val = sum(c * Fraction(q) ** i * Fraction(t) ** j for (i, j), c in num.items())
            if n_val == 0:
                continue
            lp = Poly(0); lp.add_term(left, Coef.one(0))
            rp = Poly(0); rp.add_term(right, Coef.one(0))
            acc = acc + (lp * relations[idx] * rp).scalar_mul(Coef(n_val / d_val, 0))
        return acc

    def numeric_check(self, q: int, t: int, pres: Presentation, target=None) -> tuple[bool, int]:
        """(identity holds at (q,t), number of residual monomials)."""
        tgt = resolve_target(pres, self.target if target is None else target, NumericRing(q, t, 0))
        diff = self.evaluate(q, t, pres) - tgt
        return (not diff.terms), len(diff.terms)

    # --- text ------------------------------------------------------------------
    def header_lines(self) -> list[str]:
        lines = [f"# presentation: {self.presentation}",
                 f"# presentation_hash: {self.presentation_hash}",
                 f"# target: {self.target}",
                 f"# atoms: {json.dumps(self.atoms.to_json())}",
                 f"# n_terms: {len(self.terms)}",
                 f"# producer: {self.producer}"]
        if self.anchor is not None:
            lines.append("# anchor: " + ",".join(str(v) for v in self.anchor))
        lines += [f"# note: {s}" for s in self.notes]
        return lines

    def body_lines(self) -> list[str]:
        names = self.var_names
        out = []
        for left, lab, right, num, den in self.terms:
            out.append(f"{word_str(left, names)} * {lab} * {word_str(right, names)}  ::  "
                       f"({num_str(num)}) / ({self.atoms.denom_str(den)})")
        return out

    def save(self, path: str) -> None:
        with open(path, "w") as f:
            for ln in self.header_lines():
                f.write(ln + "\n")
            for ln in self.body_lines():
                f.write(ln + "\n")

    @staticmethod
    def strip_header(text: str) -> str:
        """The body alone: every line that is not a '#' comment (for byte comparisons)."""
        return "".join(ln for ln in text.splitlines(keepends=True) if not ln.startswith("#"))

    @classmethod
    def load(cls, path: str, pres: Presentation | None = None, atoms: AtomSet | None = None,
             target: str | None = None) -> "Cofactor":
        """Read a cofactor file in the current format or any of the three legacy formats.

        A legacy file names neither its presentation nor its atoms: `pres` is then required
        (variable order, relation labels for the extract format) and `atoms` defaults to the
        historical six, AtomSet.DEFAULT, over which every legacy denominator was written."""
        header: dict = {}
        notes: list[str] = []
        body: list[str] = []
        with open(path) as f:
            for raw in f:
                ln = raw.rstrip("\n")
                s = ln.strip()
                if not s:
                    continue
                if s.startswith("#"):
                    m = _HEADER_RE.match(s)
                    if m and m.group(1) in HEADER_KEYS + ("anchor",):
                        header[m.group(1)] = m.group(2).strip()
                    elif m and m.group(1) == "note":
                        notes.append(m.group(2).strip())
                    else:
                        text = s.lstrip("#").strip()
                        if text:
                            notes.append(text)
                    continue
                body.append(ln)

        legacy = "presentation" not in header
        if legacy and pres is None:
            raise ValueError(f"{os.path.basename(path)} has no presentation header; pass pres=")
        if "atoms" in header:
            atoms = AtomSet.from_json(json.loads(header["atoms"]))
        elif atoms is None:
            atoms = AtomSet.DEFAULT
        if pres is not None:
            n = pres.n
        else:
            m = re.match(r"^n(\d+)", header["presentation"])
            if not m:
                raise ValueError(f"cannot read n from presentation slug {header['presentation']!r}")
            n = int(m.group(1))
        names = var_names(n)
        var_index = {v: i + 1 for i, v in enumerate(names)}

        anchor = None
        if "anchor" in header:
            anchor = tuple(int(v) for v in header["anchor"].split(","))
        elif legacy:
            for s in notes:
                m = _LEGACY_ANCHOR_RE.search(s)
                if m:
                    anchor = (int(m.group(1)), int(m.group(2)), 0)
                    break

        labels = pres.labels() if pres is not None else None
        terms = []
        for ln in body:
            if "::" in ln:
                terms.append(cls._parse_symbolic_line(ln, var_index, atoms))
            elif ln.lstrip()[0] in "+-" and "[" in ln:
                terms.append(cls._parse_extract_line(ln, var_index, atoms, labels))
            else:
                raise ValueError(f"cannot parse cofactor line: {ln!r}")

        if legacy:
            producer = f"cofactor_io.load from legacy {os.path.basename(path)}"
            return cls(terms, n, atoms, pres.slug(), pres.content_hash(),
                       target if target is not None else "", producer, anchor, notes)
        return cls(terms, n, atoms, header["presentation"], header.get("presentation_hash", ""),
                   header.get("target", "") if target is None else target,
                   header.get("producer", ""), anchor, notes)

    @staticmethod
    def _parse_symbolic_line(ln: str, var_index: dict, atoms: AtomSet):
        lhs, rhs = ln.split("::", 1)
        parts = [p.strip() for p in lhs.strip().split(" * ")]
        if len(parts) != 3:
            raise ValueError(f"expected 'left * LABEL * right' before '::' in {ln!r}")
        left, lab, right = parts
        # the denominator never contains '/', so the last one separates num from den
        num_s, den_s = rhs.rsplit("/", 1)
        return (parse_word(left, var_index), lab, parse_word(right, var_index),
                parse_num(num_s), atoms.parse_denom_str(_strip_outer(den_s)))

    @staticmethod
    def _parse_extract_line(ln: str, var_index: dict, atoms: AtomSet, labels):
        """'+15735/19208 * x1 * [0] * y2': a cofactor specialised at a point."""
        pieces = [p.strip() for p in ln.strip().split(" * ")]
        coef = Fraction(pieces[0])
        left, right, idx = [], [], None
        for p in pieces[1:]:
            if p.startswith("[") and p.endswith("]"):
                if idx is not None:
                    raise ValueError(f"two relation indices in {ln!r}")
                idx = int(p[1:-1])
            elif idx is None:
                left.append(p)
            else:
                right.append(p)
        if idx is None:
            raise ValueError(f"no [index] in {ln!r}")
        lab = labels[idx] if labels is not None else f"R_{idx}"
        num = {(0, 0): coef} if coef != 0 else {}
        return (parse_word("*".join(left), var_index), lab, parse_word("*".join(right), var_index),
                num, atoms.zero())

    def __repr__(self):
        return (f"Cofactor({self.presentation}, target={self.target!r}, {len(self.terms)} terms, "
                f"atoms={self.atoms.names()})")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="read a cofactor file (any format) and rewrite it in "
                                             "the current format, or just describe it")
    ap.add_argument("path")
    ap.add_argument("--n", type=int, default=None, help="strands (legacy files only)")
    ap.add_argument("--untwist", default="bigelow")
    ap.add_argument("--max-untwist-level", type=int, default=None)
    ap.add_argument("--target", default=None)
    ap.add_argument("--out", default=None, help="write the current format here")
    a = ap.parse_args()
    pres = Presentation(a.n, a.untwist, max_untwist_level=a.max_untwist_level) if a.n else None
    cof = Cofactor.load(a.path, pres=pres, target=a.target)
    print(cof)
    print(f"  labels: {cof.labels()}")
    print(f"  LCM denominator: {cof.atoms.denom_str(cof.lcm_denominator())}")
    if a.out:
        cof.save(a.out)
        print(f"  wrote {a.out}")
