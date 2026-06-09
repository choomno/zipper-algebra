"""Shared utilities for f4ncgb post-processing.

The pipeline takes a noncommutative Gröbner basis (a ``.basis`` file from
f4ncgb) and uses an Aho--Corasick automaton on its leading monomials to
enumerate standard monomials, reduce arbitrary noncommutative polynomials
modulo the ideal, and answer downstream questions about the quotient algebra
(Hilbert series, trace form, cocenter dimension, etc.).

This module collects the parsing, automaton construction, and reduction
routines used by every consumer in :mod:`f4ncgb_tools`. Two coefficient
backends are supported:

* :math:`\\mathbb{F}_p` (integer arithmetic mod a prime), with parser
  :func:`parse_basis_Fp` and reducer :func:`reduce_to_basis_Fp`.
* :math:`\\mathbb{Q}` (exact rationals via :class:`fractions.Fraction`),
  with parser :func:`parse_basis_Q` and reducer :func:`reduce_to_basis_Q`.

The Aho--Corasick construction and standard-monomial enumeration are
characteristic-independent.
"""
from __future__ import annotations

import re
from collections import deque
from fractions import Fraction
from typing import Any

Mon = tuple[str, ...]
PolyFp = list[tuple[int, Mon]]
PolyQ = list[tuple[Fraction, Mon]]


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

_FP_TERM = re.compile(r"^(-?\d+)(?:\*(.+))?$")
_Q_TERM = re.compile(r"^(\d+)(?:/(\d+))?(?:\*(.+))?$")


def _parse_mon(s: str) -> Mon:
    return tuple(v.strip() for v in s.split("*"))


def parse_basis_Fp(path: str) -> list[tuple[Mon, PolyFp]]:
    """Parse an ``f4ncgb`` ``.basis`` file with integer (mod-p) coefficients.

    Each line is one polynomial in the basis, terms listed in monomial-order
    decreasing (so the first term carries the leading monomial). Returns a
    list of ``(leading_monomial, terms)`` pairs.
    """
    polys: list[tuple[Mon, PolyFp]] = []
    with open(path) as f:
        for raw in f:
            line = raw.strip().rstrip(",")
            if not line:
                continue
            terms: PolyFp = []
            for term in line.replace(" - ", " + -").split(" + "):
                term = term.strip()
                if not term:
                    continue
                m = _FP_TERM.match(term)
                if m:
                    coef = int(m.group(1))
                    mon = _parse_mon(m.group(2)) if m.group(2) else ()
                else:
                    coef = 1
                    mon = _parse_mon(term)
                terms.append((coef, mon))
            if terms:
                polys.append((terms[0][1], terms))
    return polys


def parse_basis_Q(path: str) -> list[tuple[Mon, PolyQ]]:
    """Parse an ``f4ncgb`` ``.basis`` file with rational coefficients.

    Coefficients can be integers or fractions ``n/d``. Returns
    ``(leading_monomial, terms)`` pairs where each coefficient is a
    :class:`Fraction`.
    """
    polys: list[tuple[Mon, PolyQ]] = []
    with open(path) as f:
        for raw in f:
            line = raw.strip().rstrip(",")
            if not line:
                continue
            terms: PolyQ = []
            for term in line.replace(" - ", " + -").split(" + "):
                term = term.strip()
                if not term:
                    continue
                sign = 1
                rest = term
                if rest.startswith("-"):
                    sign = -1
                    rest = rest[1:]
                m = _Q_TERM.match(rest)
                if m:
                    num = int(m.group(1))
                    den = int(m.group(2)) if m.group(2) else 1
                    coef = Fraction(sign * num, den)
                    mon = _parse_mon(m.group(3)) if m.group(3) else ()
                else:
                    coef = Fraction(sign)
                    mon = _parse_mon(rest)
                terms.append((coef, mon))
            if terms:
                polys.append((terms[0][1], terms))
    return polys


def parse_leading_monomials(path: str) -> list[Mon]:
    """Parse only the leading monomials from a ``.basis`` file (Hilbert-series
    computations need nothing else)."""
    lms: list[Mon] = []
    with open(path) as f:
        for raw in f:
            line = raw.strip().rstrip(",")
            if not line:
                continue
            # Leading term ends at the first " + " or " - " separator after
            # column 0; everything before the first " + "/" - " is the leading
            # term (coefficient optional).
            best = -1
            for sep in (" + ", " - "):
                idx = line.find(sep, 1)
                if idx > 0 and (best < 0 or idx < best):
                    best = idx
            head = line[:best] if best > 0 else line
            m = re.match(r"^([0-9]+/[0-9]+|[0-9]+)\*?", head)
            if m and (head[m.end():].lstrip("*") or "").strip() != "":
                mon_str = head[m.end():].lstrip("*").strip()
            else:
                mon_str = head.strip()
            if mon_str == "" or mon_str.isdigit():
                # Pure-scalar basis element ⇒ ideal is the whole algebra
                return [()]
            lms.append(_parse_mon(mon_str))
    return lms


def collect_variables(polys: list[tuple[Mon, Any]]) -> list[str]:
    """Sorted alphabet of variable names appearing in the basis."""
    s: set[str] = set()
    for _, terms in polys:
        for _, mon in terms:
            s.update(mon)
    return sorted(s)


def collect_variables_from_lms(lms: list[Mon]) -> list[str]:
    """Sorted alphabet, given just the leading-monomial list."""
    s: set[str] = set()
    for lm in lms:
        s.update(lm)
    return sorted(s)


# ---------------------------------------------------------------------------
# Aho--Corasick automaton on leading monomials
# ---------------------------------------------------------------------------

# An automaton is a quadruple (delta, is_match, match_pat, nstates) where:
#   delta[(state, letter)] = next state
#   is_match[state] = True iff arriving at state means a leading monomial
#       has just been matched as a contiguous substring
#   match_pat[state] = the leading-monomial tuple matched at state, or None


def build_aho_corasick(
    lms: list[Mon], variables: list[str]
) -> tuple[dict[tuple[int, str], int], list[bool], list[Mon | None], int]:
    """Construct the Aho--Corasick automaton recognising the language of words
    that contain some element of ``lms`` as a contiguous substring.

    Returns ``(delta, is_match, match_pat, nstates)``; see module docstring.
    """
    children: list[dict[str, int]] = [dict()]
    is_match: list[bool] = [False]
    match_pat: list[Mon | None] = [None]

    for lm in lms:
        s = 0
        for v in lm:
            if v not in children[s]:
                children.append(dict())
                is_match.append(False)
                match_pat.append(None)
                children[s][v] = len(children) - 1
            s = children[s][v]
        is_match[s] = True
        match_pat[s] = lm

    nstates = len(children)
    fail = [0] * nstates
    delta: dict[tuple[int, str], int] = {}
    queue: deque[int] = deque()
    for v in variables:
        if v in children[0]:
            c = children[0][v]
            fail[c] = 0
            delta[(0, v)] = c
            queue.append(c)
        else:
            delta[(0, v)] = 0
    while queue:
        s = queue.popleft()
        for v in variables:
            if v in children[s]:
                c = children[s][v]
                fs = fail[s]
                while fs != 0 and v not in children[fs]:
                    fs = fail[fs]
                if v in children[fs] and children[fs][v] != c:
                    fail[c] = children[fs][v]
                else:
                    fail[c] = 0
                if is_match[fail[c]] and not is_match[c]:
                    is_match[c] = True
                    match_pat[c] = match_pat[fail[c]]
                delta[(s, v)] = c
                queue.append(c)
            else:
                delta[(s, v)] = delta.get((fail[s], v), 0)
    return delta, is_match, match_pat, nstates


def enumerate_standard_monomials(
    delta: dict[tuple[int, str], int],
    is_match: list[bool],
    variables: list[str],
    max_total: int | None = None,
) -> list[Mon] | None:
    """All standard monomials (words containing no leading monomial as a
    substring), ordered by (length, lex). The empty monomial is included.
    Returns ``None`` if the count exceeds ``max_total`` (used as a budget
    to abort on possibly-infinite quotients).
    """
    result: list[tuple[int, Mon]] = [(0, ())]
    frontier: deque[tuple[int, Mon]] = deque([(0, ())])
    while frontier:
        state, mon = frontier.popleft()
        for v in variables:
            t = delta[(state, v)]
            if is_match[t]:
                continue
            new_mon = mon + (v,)
            result.append((t, new_mon))
            frontier.append((t, new_mon))
            if max_total is not None and len(result) > max_total:
                return None
    result.sort(key=lambda sm: (len(sm[1]), sm[1]))
    return [mon for _, mon in result]


def find_leading_match(
    word: Mon,
    delta: dict[tuple[int, str], int],
    is_match: list[bool],
    match_pat: list[Mon | None],
) -> tuple[int, int, Mon] | None:
    """The leftmost substring of ``word`` that equals some leading monomial,
    returned as ``(start, end, lm)``. ``None`` if ``word`` is standard.
    """
    state = 0
    for i, v in enumerate(word):
        if (state, v) not in delta:
            return None
        state = delta[(state, v)]
        if is_match[state]:
            lm = match_pat[state]
            assert lm is not None
            end = i + 1
            return (end - len(lm), end, lm)
    return None


# ---------------------------------------------------------------------------
# Reduction of arbitrary polynomials modulo the Gröbner basis
# ---------------------------------------------------------------------------

# A reduction sends a polynomial p ∈ k⟨x⟩ to its unique standard normal form:
# the term-by-term remainder after rewriting every occurrence of every leading
# monomial via its corresponding tail.
#
# We represent the polynomial during reduction as a worklist of (coef, word)
# pairs (the same word may appear multiple times with different coefficients;
# duplicates accumulate into the result sparse vector at the end).
#
# The caller supplies ``lm_to_tail``: for each leading monomial L, the list
# of tail terms ``[(c_1, m_1), ..., (c_k, m_k)]`` such that the substitution
# ``L ↦ c_1 m_1 + ... + c_k m_k`` reduces the basis polynomial whose leading
# monomial is L. The coefficients in the tail are already normalised so that
# the leading coefficient is 1 (i.e. the tail equals -(non-leading part) /
# leading_coef).


def lm_to_tail_Fp(
    polys: list[tuple[Mon, PolyFp]], p: int
) -> dict[Mon, PolyFp]:
    """For each polynomial in the basis, compute its reduction-substitution
    tail in :math:`\\mathbb{F}_p`.
    """
    out: dict[Mon, PolyFp] = {}
    for lm, terms in polys:
        lead = terms[0][0]
        inv = pow(lead, p - 2, p)
        out[lm] = [((-c * inv) % p, mon) for c, mon in terms[1:]]
    return out


def lm_to_tail_Q(
    polys: list[tuple[Mon, PolyQ]]
) -> dict[Mon, PolyQ]:
    """For each polynomial in the basis, compute its reduction-substitution
    tail in :math:`\\mathbb{Q}`.
    """
    out: dict[Mon, PolyQ] = {}
    for lm, terms in polys:
        inv = Fraction(1) / terms[0][0]
        out[lm] = [(-c * inv, mon) for c, mon in terms[1:]]
    return out


def reduce_to_basis_Fp(
    word_terms: PolyFp,
    mon_idx: dict[Mon, int],
    lm_to_tail: dict[Mon, PolyFp],
    delta: dict[tuple[int, str], int],
    is_match: list[bool],
    match_pat: list[Mon | None],
    p: int,
) -> dict[int, int]:
    """Reduce a noncommutative polynomial (a list of ``(coef, word)`` terms)
    to its standard normal form over :math:`\\mathbb{F}_p`, returned as a
    sparse vector ``{basis_idx: coef}``.

    Each step finds the leftmost leading-monomial substring (via the AC
    automaton) and rewrites it via its precomputed tail. The reduction
    terminates because the term order is a well-ordering on monomials.
    """
    out: dict[int, int] = {}
    worklist = list(word_terms)
    while worklist:
        coef, word = worklist.pop()
        if coef == 0:
            continue
        idx = mon_idx.get(word)
        if idx is not None:
            new = (out.get(idx, 0) + coef) % p
            if new == 0 and idx in out:
                del out[idx]
            else:
                out[idx] = new
            continue
        match = find_leading_match(word, delta, is_match, match_pat)
        if match is None:
            raise RuntimeError(
                f"Standard word {word!r} not in the enumerated basis -- "
                "the Gröbner basis may be incomplete."
            )
        start, end, lm = match
        u, v = word[:start], word[end:]
        for c, tail_mon in lm_to_tail[lm]:
            worklist.append(((coef * c) % p, u + tail_mon + v))
    return out


def reduce_to_basis_Q(
    word_terms: PolyQ,
    mon_idx: dict[Mon, int],
    lm_to_tail: dict[Mon, PolyQ],
    delta: dict[tuple[int, str], int],
    is_match: list[bool],
    match_pat: list[Mon | None],
) -> dict[int, Fraction]:
    """Reduction over :math:`\\mathbb{Q}` using :class:`Fraction` arithmetic.
    See :func:`reduce_to_basis_Fp` for the algorithm; the only difference
    is the coefficient ring.
    """
    out: dict[int, Fraction] = {}
    worklist = list(word_terms)
    while worklist:
        coef, word = worklist.pop()
        if coef == 0:
            continue
        idx = mon_idx.get(word)
        if idx is not None:
            new = out.get(idx, Fraction(0)) + coef
            if new == 0 and idx in out:
                del out[idx]
            else:
                out[idx] = new
            continue
        match = find_leading_match(word, delta, is_match, match_pat)
        if match is None:
            raise RuntimeError(
                f"Standard word {word!r} not in the enumerated basis -- "
                "the Gröbner basis may be incomplete."
            )
        start, end, lm = match
        u, v = word[:start], word[end:]
        for c, tail_mon in lm_to_tail[lm]:
            worklist.append((coef * c, u + tail_mon + v))
    return out
