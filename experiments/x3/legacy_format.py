"""Historical X3 text readers and adapters to maintained exact arithmetic.

The compact expressions use braid-algebra equality; their normalization is in
simplify_cofactor. The parser is retained for the recorded four-factor format.
Arithmetic and presentation construction are shared with zipper_algebra.
"""
from __future__ import annotations
import os
import re
from pathlib import Path
from fractions import Fraction
from typing import Dict, List, Tuple
from zipper_algebra import Presentation
from zipper_algebra.cofactor.ncpoly import (
    QTPoly, NCPoly, NCMon,
    qt_zero, qt_const, qt_q, qt_t, qt_monomial, qt_add, qt_sub, qt_neg,
    qt_mul, qt_pow, qt_eq, qt_is_zero, qt_str,
    nc_zero, nc_const, nc_var, nc_add, nc_sub, nc_neg, nc_mul,
    nc_scalar, nc_is_zero, nc_size,
)
ROOT = str(Path(__file__).resolve().parents[2])
COF_FILE = os.path.join(ROOT, "logs_and_traces/f4ncgb/cofactor-X3-n4-Qqt-symbolic.txt")
N_STRANDS = 4
VAR_NAMES = Presentation(N_STRANDS).var_names()
VAR_INDEX = {name: i + 1 for i, name in enumerate(VAR_NAMES)}


def nc_var_named(name):
    return nc_var(VAR_INDEX[name])


def build_inputs_and_X3(presentation=None, target=(3, 1)):
    """Use the specified presentation's exact symbolic relations and target."""
    pres = presentation or Presentation(N_STRANDS)
    relations, sliders = pres.symbolic_ncpolys()
    if target not in sliders:
        raise ValueError(f"target S{list(target)} not available at n={pres.n}")
    return relations, sliders[target]


_NUM_TERM_RE = re.compile(
    r"\s*(?P<sign>[+\-])?\s*"
    r"(?P<coef>-?\d+)\s*\*\s*q\^(?P<aq>\d+)\s*\*\s*t\^(?P<bt>\d+)\s*"
)

def parse_numerator(s: str) -> QTPoly:
    s = s.strip()
    # strip outer parens if any
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()
    out: QTPoly = {}
    # tokenize by splitting at top-level ' + ' (file has explicit signs)
    # but signs may be '+' or '-'.  We scan with the regex above sequentially.
    pos = 0
    first = True
    while pos < len(s):
        m = _NUM_TERM_RE.match(s, pos)
        if not m:
            # skip a stray whitespace or sign char if present
            if s[pos].isspace() or s[pos] == "+":
                pos += 1
                continue
            raise ValueError(
                f"Cannot parse numerator at position {pos}: "
                f"{s[pos:pos+40]!r} in full string {s!r}"
            )
        sign_g = m.group("sign")
        coef = int(m.group("coef"))
        if sign_g == "-":
            coef = -coef
        aq = int(m.group("aq"))
        bt = int(m.group("bt"))
        s_new = out.get((aq, bt), Fraction(0)) + coef
        if s_new == 0:
            out.pop((aq, bt), None)
        else:
            out[(aq, bt)] = s_new
        pos = m.end()
        first = False
    return out

def parse_denominator(s: str) -> Tuple[int, int, int, int]:
    """Return (a_q, b_t, k_qm1, k_tpq2) for D = q^a * t^b * (q-1)^k_qm1 * (t+q^2)^k_tpq2."""
    s = s.strip()
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()
    a_q = 0
    b_t = 0
    k_qm1 = 0
    k_tpq2 = 0
    parts = [p.strip() for p in s.split("*") if p.strip()]
    for p in parts:
        if p == "q":
            a_q += 1
        elif p == "t":
            b_t += 1
        elif p.startswith("q^"):
            a_q += int(p[2:])
        elif p.startswith("t^"):
            b_t += int(p[2:])
        elif p == "(q-1)":
            k_qm1 += 1
        elif p == "(t+q^2)":
            k_tpq2 += 1
        else:
            raise ValueError(f"Unknown denominator factor: {p!r} in {s!r}")
    return a_q, b_t, k_qm1, k_tpq2

def parse_monomial(s: str) -> NCMon:
    s = s.strip()
    if s == "" or s == "()":
        return ()
    pieces = [p.strip() for p in s.split("*") if p.strip()]
    out: List[int] = []
    for p in pieces:
        if p not in VAR_INDEX:
            raise ValueError(f"Unknown var in monomial: {p!r}")
        out.append(VAR_INDEX[p])
    return tuple(out)

_LINE_RE = re.compile(
    r"^\s*(?P<left>.+?)\s+\*\s+(?P<lbl>[A-Z][A-Za-z0-9_]*)\s+\*\s+(?P<right>.+?)\s*::\s*"
    r"(?P<num>\([^()]*\))\s*/\s*(?P<den>\(.+\))\s*$"
)

def parse_cofactor_line(ln: str):
    from .relations_naming import parse_label
    m = _LINE_RE.match(ln)
    if not m:
        raise ValueError(f"Cannot parse cofactor line: {ln!r}")
    left = parse_monomial(m.group("left"))
    idx = parse_label(m.group("lbl"))
    right = parse_monomial(m.group("right"))
    num = parse_numerator(m.group("num"))
    den = parse_denominator(m.group("den"))
    return left, idx, right, num, den
