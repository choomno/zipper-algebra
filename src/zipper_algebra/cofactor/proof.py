#!/usr/bin/env python3
"""f4ncgb output at one point: the basis and its non-expanded proof.

`notes/refactor/SPEC.md`, data types `Basis` and `Proof`.  The parsers are the ones of
`extract_cofactor.py` (`parse_basis`, `parse_noex_proof`), returning objects instead of bare
lists, with the variable names passed in rather than read from module globals.

    var_names = Presentation(4).var_names()                       # ['x1', 'x2', 'x3', 'y1', 'y2', 'y3']
    basis = Basis.load("basis-n4-Q.txt", char=0, var_names=var_names)
    proof = Proof.load("proof-n4-Q-noex.txt", var_names=var_names)
    proof.expand(k)   # basis[k] = sum coef * left * input[idx] * right, exact
    proof.deps(k)     # basis indices basis[k]'s proof transitively references

Monomial order.  f4ncgb orders words degree first, then lexicographically in the variable order
of the .ms header (x1 < x2 < ... < y_{n-1}); `mon_sort_key` is that order on 1-based index tuples
and the leading word is the maximum.  Nothing here assumes n = 4.

Keys.  A basis element is identified across points by its leading word as text ("y1*y2*x1"),
never by its index: indices are an accident of the run, leading words are a property of the
ideal and the order.  When two elements share a leading word (possible with -r 0, never in a
reduced basis) every element with that word is keyed by the word plus its full text instead.

Expansion.  The non-expanded proof writes basis[k] as a combination of inputs and of OTHER
basis elements; `expand` substitutes until only inputs remain.  It walks the reference graph in
topological order, so a proof line may reference a later basis element; only a genuine cycle
(basis[k] needed to express basis[k]) is an error.  `extract_cofactor.expand_proofs` raised
NotImplementedError on any forward reference; on the X3 artifacts none occur, and the two agree
term for term there (tests_proof_skeleton.py).
"""
from __future__ import annotations

import re
from fractions import Fraction

from ..zipper_in_ms import Coef, Poly


Word = tuple[int, ...]


# ---------------------------------------------------------------------------
# Monomial order
# ---------------------------------------------------------------------------

def mon_sort_key(mon: Word) -> tuple:
    """Degree-lex key in the variable order of the .ms header (max = leading)."""
    return (len(mon),) + tuple(mon)


def leading_term(p: Poly, order=mon_sort_key) -> tuple[Word, Coef]:
    if not p.terms:
        raise ValueError("leading_term of the zero polynomial")
    lm = max(p.terms.keys(), key=order)
    return lm, p.terms[lm]


def word_text(mon: Word, var_names: list[str]) -> str:
    """'x1*y2' for (1, 5); '1' for the empty word."""
    return "*".join(var_names[i - 1] for i in mon) if mon else "1"


# ---------------------------------------------------------------------------
# Tokenizers (as in extract_cofactor.py)
# ---------------------------------------------------------------------------

_FRAC_RE = re.compile(r"^(-?\d+(?:/\d+)?)$")


def _parse_coef_str(s: str) -> Fraction:
    if "/" in s:
        num, den = s.split("/")
        return Fraction(int(num), int(den))
    return Fraction(int(s))


def _split_terms(line: str) -> list[str]:
    """Split a polynomial line into signed term tokens.

    Tokens are sometimes joined to a leading sign ('-4/35*[0]') and sometimes separated by a
    bare ' + ' or ' - '; whitespace-split and re-attach pending signs."""
    line = line.strip().rstrip(",").strip()
    if not line:
        return []
    out: list[str] = []
    pending: str | None = None
    for tok in line.split():
        if tok == "+":
            pending = "+"
            continue
        if tok == "-":
            pending = "-"
            continue
        if pending is None:
            out.append(tok if tok[0] in "+-" else "+" + tok)
        else:
            if tok[0] in "+-":
                combined = pending + tok[0]
                sign = "+" if combined.count("-") % 2 == 0 else "-"
                out.append(sign + tok[1:])
            else:
                out.append(pending + tok)
            pending = None
    return out


def _split_coef(tok: str) -> tuple[Fraction, list[str]]:
    """'-4/35*x1*[0]*y1' -> (-4/35, ['x1', '[0]', 'y1']); a missing coefficient is 1."""
    sign_is_neg = (tok[0] == "-")
    pieces = tok[1:].split("*")
    coef = Fraction(1)
    if _FRAC_RE.match(pieces[0]):
        coef = _parse_coef_str(pieces[0])
        rest = pieces[1:]
    else:
        rest = pieces
    return (-coef if sign_is_neg else coef), rest


# ---------------------------------------------------------------------------
# Basis
# ---------------------------------------------------------------------------

class Basis:
    """Parsed f4ncgb basis at a point.

    polys        : list[Poly] over Q (char 0) or GF(char)
    lm           : leading word of each element under `order`
    key          : stable text key per element (leading word; fingerprint on collision)
    index_by_key : key -> index in this file
    """

    def __init__(self, polys: list[Poly], char: int, var_names: list[str], order=mon_sort_key):
        self.polys = polys
        self.char = char
        self.var_names = list(var_names)
        self.order = order
        self.lm: list[Word] = []
        for k, p in enumerate(polys):
            if not p.terms:
                raise ValueError(f"basis[{k}] is the zero polynomial")
            self.lm.append(leading_term(p, order)[0])
        self.key = self._keys()
        self.index_by_key = {key: k for k, key in enumerate(self.key)}

    def _keys(self) -> list[str]:
        lm_text = [word_text(m, self.var_names) for m in self.lm]
        count: dict[str, int] = {}
        for s in lm_text:
            count[s] = count.get(s, 0) + 1
        keys = []
        for s, p in zip(lm_text, self.polys):
            keys.append(s if count[s] == 1 else f"{s} | {p.to_msolve_string(self.var_names)}")
        # identical polynomials listed twice: number the repeats so keys stay unique
        seen: dict[str, int] = {}
        for i, s in enumerate(keys):
            seen[s] = seen.get(s, 0) + 1
            if seen[s] > 1:
                keys[i] = f"{s} #{seen[s]}"
        return keys

    def __len__(self) -> int:
        return len(self.polys)

    def index(self, key: str) -> int:
        """Index of the element with this key here; KeyError names the key if absent."""
        if key not in self.index_by_key:
            raise KeyError(f"no basis element with leading word {key!r} at this point")
        return self.index_by_key[key]

    @classmethod
    def load(cls, path: str, char: int, var_names: list[str], order=mon_sort_key) -> "Basis":
        return cls(parse_basis(path, char, var_names), char, var_names, order)


def parse_basis(path: str, char: int, var_names: list[str]) -> list[Poly]:
    """One Poly per non-empty line of an f4ncgb basis file."""
    var_index = {name: i + 1 for i, name in enumerate(var_names)}
    polys: list[Poly] = []
    with open(path) as f:
        for raw in f:
            line = raw.strip().rstrip(",")
            if not line:
                continue
            p = Poly(char)
            for tok in _split_terms(line):
                coef, names = _split_coef(tok)
                for v in names:
                    if v not in var_index:
                        raise ValueError(f"unknown variable {v!r} in basis term {tok!r} of {path}")
                p.add_term(tuple(var_index[v] for v in names), Coef(coef, char))
            polys.append(p)
    return polys


# ---------------------------------------------------------------------------
# Proof
# ---------------------------------------------------------------------------

Term = tuple[Fraction, Word, str, int, Word]        # (coef, left, kind, idx, right)
Expanded = tuple[Fraction, Word, int, Word]         # (coef, left, input_idx, right)


def _parse_noex_term(tok: str, var_index: dict[str, int]) -> Term:
    """One term of the non-expanded proof: coef * left * ([ik] | [k]) * right."""
    coef, pieces = _split_coef(tok)
    left: list[int] = []
    right: list[int] = []
    idx: int | None = None
    kind: str | None = None
    for p in pieces:
        if p.startswith("[i") and p.endswith("]"):
            if idx is not None:
                raise ValueError(f"two indices in {tok!r}")
            idx, kind = int(p[2:-1]), "input"
        elif p.startswith("[") and p.endswith("]"):
            if idx is not None:
                raise ValueError(f"two indices in {tok!r}")
            idx, kind = int(p[1:-1]), "basis"
        else:
            if p not in var_index:
                raise ValueError(f"unknown variable {p!r} in {tok!r}")
            (left if idx is None else right).append(var_index[p])
    if idx is None or kind is None:
        raise ValueError(f"no [..] index in {tok!r}")
    return coef, tuple(left), kind, idx, tuple(right)


def parse_noex_proof(path: str, var_names: list[str]) -> list[list[Term]]:
    """One term list per line; an empty line is an element with no proof terms."""
    var_index = {name: i + 1 for i, name in enumerate(var_names)}
    out: list[list[Term]] = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            out.append([_parse_noex_term(tok, var_index) for tok in _split_terms(line)]
                       if line else [])
    return out


def _accumulate(bucket: dict, key, value: Fraction) -> None:
    """bucket[key] += value, deleting the key when the sum cancels (legacy behaviour, so that
    the term order of the result matches extract_cofactor.expand_proofs exactly)."""
    nv = bucket.get(key, Fraction(0)) + value
    if nv == 0:
        bucket.pop(key, None)
    else:
        bucket[key] = nv


class Proof:
    """Non-expanded proof: terms[k] writes basis[k] over inputs and other basis elements."""

    def __init__(self, terms: list[list[Term]]):
        self.terms = terms
        self._expanded: dict[int, dict[tuple[Word, int, Word], Fraction]] = {}
        self._deps: dict[int, set[int]] = {}

    def __len__(self) -> int:
        return len(self.terms)

    @classmethod
    def load(cls, path: str, var_names: list[str]) -> "Proof":
        return cls(parse_noex_proof(path, var_names))

    # --- the reference graph -------------------------------------------------
    def refs(self, k: int) -> list[int]:
        """Basis indices basis[k]'s own line references, in line order, without repeats."""
        seen: list[int] = []
        for _c, _l, kind, idx, _r in self.terms[k]:
            if kind == "basis" and idx not in seen:
                seen.append(idx)
        return seen

    def _topological(self, k: int) -> list[int]:
        """Indices reachable from k, each after everything it references (iterative DFS).
        Raises ValueError on a cycle, naming it."""
        order: list[int] = []
        state: dict[int, int] = {}          # 1 = on the stack, 2 = finished
        stack: list[tuple[int, int]] = [(k, 0)]
        state[k] = 1
        path: list[int] = [k]
        while stack:
            node, i = stack[-1]
            if node in self._expanded:
                state[node] = 2
                stack.pop(); path.pop()
                continue
            refs = self.refs(node)
            if i < len(refs):
                stack[-1] = (node, i + 1)
                j = refs[i]
                if j >= len(self.terms) or j < 0:
                    raise ValueError(f"basis[{node}] references basis[{j}]; the proof has "
                                     f"{len(self.terms)} lines")
                st = state.get(j, 0)
                if st == 1:
                    cyc = path[path.index(j):] + [j]
                    raise ValueError("cyclic proof: " + " -> ".join(f"basis[{x}]" for x in cyc))
                if st == 0:
                    state[j] = 1
                    stack.append((j, 0)); path.append(j)
            else:
                state[node] = 2
                order.append(node)
                stack.pop(); path.pop()
        return order

    # --- expansion -------------------------------------------------------------
    def _expand_one(self, k: int) -> None:
        """Expand basis[k] assuming everything it references is already expanded."""
        bucket: dict[tuple[Word, int, Word], Fraction] = {}
        for coef, left, kind, idx, right in self.terms[k]:
            if kind == "input":
                _accumulate(bucket, (left, idx, right), coef)
            else:
                for (sub_l, sub_i, sub_r), sub_coef in self._expanded[idx].items():
                    _accumulate(bucket, (left + sub_l, sub_i, sub_r + right), coef * sub_coef)
        self._expanded[k] = bucket

    def expand_dict(self, k: int) -> dict[tuple[Word, int, Word], Fraction]:
        """basis[k] over the inputs as {(left, input_idx, right): coef}; memoised."""
        if k not in self._expanded:
            for j in self._topological(k):
                if j not in self._expanded:
                    self._expand_one(j)
        return self._expanded[k]

    def expand(self, k: int) -> list[Expanded]:
        """basis[k] = sum coef * left * input[idx] * right, as (coef, left, idx, right) tuples."""
        return [(c, l, i, r) for (l, i, r), c in self.expand_dict(k).items()]

    def expand_all(self) -> list[list[Expanded]]:
        return [self.expand(k) for k in range(len(self.terms))]

    def deps(self, k: int) -> set[int]:
        """Basis indices basis[k] transitively references (k itself only through a cycle)."""
        if k in self._deps:
            return self._deps[k]
        out: set[int] = set()
        work = list(self.refs(k))
        while work:
            j = work.pop()
            if j in out:
                continue
            out.add(j)
            if j in self._deps:
                out |= self._deps[j]
            else:
                work.extend(self.refs(j))
        self._deps[k] = out
        return out

    def shape(self, k: int, basis: "Basis | None" = None) -> tuple:
        """Term structure of basis[k]'s proof without coefficients: (left, kind, ref, right)
        per term, where ref is the input index or, given a basis, the referenced element's
        key.  Two points whose proofs have the same shape for every element in a skeleton's
        closure carry the same lift slots."""
        out = []
        for _c, left, kind, idx, right in self.terms[k]:
            ref = basis.key[idx] if (basis is not None and kind == "basis") else idx
            out.append((left, kind, ref, right))
        return tuple(out)

    # --- checks ------------------------------------------------------------------
    def reconstructs(self, k: int, basis: Basis, inputs: list[Poly]) -> bool:
        """Does expand(k) rebuild basis[k] exactly from the input polynomials?"""
        char = basis.char
        acc = Poly(char)
        for coef, u, idx, v in self.expand(k):
            lp = Poly(char); lp.add_term(u, Coef.one(char))
            rp = Poly(char); rp.add_term(v, Coef.one(char))
            acc = acc + (lp * inputs[idx] * rp).scalar_mul(Coef(coef, char))
        return not (acc - basis.polys[k]).terms
