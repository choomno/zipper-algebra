"""Shared utilities for f4ncgb post-processing.

The pipeline takes a noncommutative Gröbner basis (a ``.basis`` file from
f4ncgb) and uses an Aho--Corasick automaton on its leading monomials to
enumerate standard monomials, reduce arbitrary noncommutative polynomials
modulo the ideal, and answer downstream questions about the quotient algebra
(Hilbert series, trace form, cocenter dimension, etc.).

This module collects the parsing, automaton construction, and reduction
routines used by every consumer in :mod:`tools`. Two coefficient
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

_TERM = re.compile(r"^(?P<sign>[+-]?)(?:(?P<coef>\d+(?:/\d+)?)(?:\*(?P<word>.+))?|(?P<bare>[A-Za-z_].*))$")
_VARIABLE = re.compile(r"^[A-Za-z_][A-Za-z_0-9]*$")


def _parse_mon(s: str) -> Mon:
    words = tuple(v.strip() for v in s.split("*"))
    if any(not _VARIABLE.fullmatch(v) for v in words):
        raise ValueError(f"invalid word {s!r}")
    return words


def _line_terms(line: str):
    # Both pretty printed bases and the emitter's compact +/- syntax are accepted.
    if re.search(r"[A-Za-z_0-9]\s+[A-Za-z_0-9]", line):
        raise ValueError(f"missing operator in polynomial {line!r}")
    line = re.sub(r"\+\s*-", "-", line)
    compact = re.sub(r"\s+", "", line)
    position = 0
    for match in re.finditer(r"[+-]?[^+-]+", compact):
        if match.start() != position:
            raise ValueError(f"invalid polynomial {line!r}")
        position = match.end()
        yield match.group()
    if position != len(compact):
        raise ValueError(f"invalid polynomial {line!r}")


def _parse_term(term: str, *, integers=False):
    match = _TERM.fullmatch(term)
    if match is None:
        raise ValueError(f"invalid polynomial term {term!r}")
    sign = -1 if match.group("sign") == "-" else 1
    coefficient = sign * Fraction(match.group("coef") or 1)
    if integers and coefficient.denominator != 1:
        raise ValueError("a finite-field basis must use integer residue coefficients")
    word = match.group("word") or match.group("bare")
    monomial = _parse_mon(word) if word else ()
    return (int(coefficient) if integers else coefficient), monomial


def _basis_rows(path: str, *, integers=False):
    with open(path) as stream:
        for line_number, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line:
                continue
            if line.endswith(","):
                line = line[:-1].strip()
            if not line or "," in line:
                raise ValueError(f"{path}:{line_number}: invalid polynomial separator")
            coefficients = {}
            try:
                for term in _line_terms(line):
                    coefficient, word = _parse_term(term, integers=integers)
                    coefficients[word] = coefficients.get(word, 0) + coefficient
            except (ValueError, ZeroDivisionError) as error:
                raise ValueError(f"{path}:{line_number}: {error}") from error
            terms = [(coefficient, word) for word, coefficient in coefficients.items()
                     if coefficient]
            if terms:
                yield terms[0][1], terms


def parse_basis_Fp(path: str) -> list[tuple[Mon, PolyFp]]:
    """Read integer residue coefficients and ordered words from an f4ncgb basis.

    The modulus is supplied later to the reducer. Fractions are rejected here,
    since this reader cannot convert them without knowing that modulus.
    """
    return list(_basis_rows(path, integers=True))


def parse_basis_Q(path: str) -> list[tuple[Mon, PolyQ]]:
    """Read exact rational coefficients and ordered words from an f4ncgb basis.

    Signed bare words, rational leading coefficients and scalar relations are
    supported. Zero polynomials are omitted, and repeated words are collected.
    The input order must be the engine's decreasing monomial order.
    """
    return list(_basis_rows(path))


def parse_leading_monomials(path: str) -> list[Mon]:
    """Read leading words, using the same signed-term grammar as both full readers.

    A nonzero scalar relation gives the unit ideal, represented by ``[()]``.
    An explicit zero polynomial contributes no leading word.
    """
    lms = []
    for word, _terms in _basis_rows(path):
        if not word:
            return [()]
        lms.append(word)
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
    if is_match[0]:
        return []  # The empty leading word divides every word, including the unit.
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
    if is_match[0]:
        return (0, 0, ())
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
        try:
            inv = pow(lead, -1, p)
        except ValueError:
            raise ValueError(f"leading coefficient {lead} is not invertible modulo {p}") from None
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
        coef %= p
        if coef == 0:
            continue
        idx = mon_idx.get(word)
        if idx is not None:
            new = (out.get(idx, 0) + coef) % p
            if new == 0:
                out.pop(idx, None)
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
