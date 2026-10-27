#!/usr/bin/env python3
"""Reduce a target against a basis at a point and record the steps as a Skeleton.

`notes/refactor/SPEC.md`, data type `Skeleton`.  The reduction is
`extract_cofactor.reduce_against_basis` with the basis passed as a `Basis`; the target is ANY
`Poly` at the anchor (a slider from `Presentation.target(i, j, ring)`, an EXTRA, a product of
words) -- this module does not know what it is proving.

    steps, remainder = reduce(target, basis)
    skel = Skeleton.from_reduction(pres, "S[3][1]", (5, 7, 0), steps, remainder, basis, proof)
    skel.save(path); Skeleton.load(path)
    skel.step_indices(other_basis)        # the same steps at another point, or KeyError

A step records one subtraction  p -> p - mult * left * basis[k] * right  and names basis[k] by
its KEY (leading word), never by k: the index changes from run to run, the leading word does
not.  `basis_keys` is the transitive closure of the steps' keys under the proof's references,
i.e. every basis element whose lift the target's cofactor needs (25 for X3 at n=4).  The
multipliers are stored at the anchor; lifting them to Q(q,t) is compose.py's job.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from fractions import Fraction

from ..zipper_in_ms import Coef, Poly
from .proof import Basis, Proof, Word, Expanded, leading_term, mon_sort_key


Anchor = tuple[int, int, int]                       # (q, t, char); char 0 = Q
Step = tuple[Word, int, Word, Coef]                 # (left, basis_idx, right, mult)


# ---------------------------------------------------------------------------
# Leading-word reduction
# ---------------------------------------------------------------------------

def find_substring(mu: Word, lam: Word) -> tuple[Word, Word] | None:
    """(u, v) with mu = u + lam + v at the leftmost occurrence, or None."""
    n, m = len(mu), len(lam)
    if m == 0 or m > n:
        return None
    for i in range(n - m + 1):
        if mu[i:i + m] == lam:
            return mu[:i], mu[i + m:]
    return None


def reduce(target: Poly, basis: Basis, order=mon_sort_key) -> tuple[list[Step], Poly]:
    """Leading-word-driven non-commutative reduction of `target` by `basis`.

    Returns (steps, remainder).  Each step is (left, k, right, mult) meaning the running
    polynomial lost mult * left * basis[k] * right.  Divisor choice: among the basis elements
    whose leading word occurs in the current leading word, the one with the shortest leading
    word, then lex-smallest, then lowest index -- the same canonical choice as the legacy
    reducer, so the X3 skeleton is the historical one (basis[17], basis[16], basis[102]).
    """
    if target.char != basis.char:
        raise ValueError(f"target over char {target.char}, basis over char {basis.char}")
    char = basis.char
    lt_basis = sorted(((basis.lm[k], basis.polys[k].terms[basis.lm[k]], k)
                       for k in range(len(basis))),
                      key=lambda t: (len(t[0]), t[0], t[2]))

    steps: list[Step] = []
    p = Poly(char)
    for mon, coef in target.terms.items():
        p.add_term(mon, coef)
    remainder = Poly(char)

    while p.terms:
        lm_p, lc_p = leading_term(p, order)
        for lm_g, lc_g, k in lt_basis:
            split = find_substring(lm_p, lm_g)
            if split is None:
                continue
            u, v = split
            if char == 0:
                mult = Coef(Fraction(lc_p.v) / Fraction(lc_g.v), char)
            else:
                mult = Coef(int(lc_p.v) * pow(int(lc_g.v), -1, char), char)
            for mg, cg in basis.polys[k].terms.items():
                p.add_term(u + mg + v, -(mult * cg))
            steps.append((u, k, v, mult))
            break
        else:
            remainder.add_term(lm_p, lc_p)
            del p.terms[lm_p]
    return steps, remainder


# ---------------------------------------------------------------------------
# Skeleton
# ---------------------------------------------------------------------------

def closure(indices, proof: Proof) -> set[int]:
    """The given basis indices and everything their proofs transitively reference."""
    out: set[int] = set()
    for k in indices:
        out.add(k)
        out |= proof.deps(k)
    return out


@dataclass
class Skeleton:
    presentation: str                       # Presentation.slug()
    presentation_hash: str                  # Presentation.content_hash()
    target: str                             # "S[3][1]", an EXTRA label, a word expression
    anchor: Anchor
    steps: list[dict] = field(default_factory=list)   # {left, basis_key, right, mult: Fraction}
    basis_keys: list[str] = field(default_factory=list)
    remainder_terms: int = 0

    # --- construction --------------------------------------------------------
    @classmethod
    def from_reduction(cls, pres, target: str, anchor: Anchor, steps: list[Step],
                       remainder: Poly, basis: Basis, proof: Proof) -> "Skeleton":
        """Name each step's basis element by key and close the step set under the proof's
        references.  `pres` is a Presentation (slug and content_hash are taken from it)."""
        if len(proof) != len(basis):
            raise ValueError(f"proof has {len(proof)} lines, basis {len(basis)} elements")
        step_dicts = [{"left": u, "basis_key": basis.key[k], "right": v,
                       "mult": Fraction(mult.v)} for (u, k, v, mult) in steps]
        needed = closure((k for _u, k, _v, _m in steps), proof)
        keys = [basis.key[k] for k in sorted(needed)]     # anchor index order
        return cls(pres.slug(), pres.content_hash(), target, tuple(anchor), step_dicts,
                   keys, len(remainder.terms))

    @property
    def in_ideal(self) -> bool:
        return self.remainder_terms == 0

    def step_keys(self) -> list[str]:
        return [s["basis_key"] for s in self.steps]

    # --- re-identification at a point --------------------------------------
    def step_indices(self, basis: Basis) -> list[int]:
        """The steps' basis elements at this basis, by key.  KeyError if one is absent
        (the point does not carry this skeleton)."""
        return [basis.index(s["basis_key"]) for s in self.steps]

    def closure_indices(self, basis: Basis) -> list[int]:
        return [basis.index(key) for key in self.basis_keys]

    # --- persistence ---------------------------------------------------------
    def to_json(self) -> dict:
        q, t, char = self.anchor
        return {
            "presentation": self.presentation,
            "presentation_hash": self.presentation_hash,
            "target": self.target,
            "anchor": {"q": q, "t": t, "char": char},
            "n_steps": len(self.steps),
            "remainder_terms": self.remainder_terms,
            "steps": [{"left": list(s["left"]), "basis_key": s["basis_key"],
                       "right": list(s["right"]), "mult": str(s["mult"])} for s in self.steps],
            "basis_keys": list(self.basis_keys),
        }

    @classmethod
    def from_json(cls, d: dict) -> "Skeleton":
        a = d["anchor"]
        steps = [{"left": tuple(s["left"]), "basis_key": s["basis_key"],
                  "right": tuple(s["right"]), "mult": Fraction(s["mult"])} for s in d["steps"]]
        return cls(d["presentation"], d["presentation_hash"], d["target"],
                   (int(a["q"]), int(a["t"]), int(a["char"])), steps,
                   list(d["basis_keys"]), int(d["remainder_terms"]))

    def save(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_json(), f, indent=1)

    @classmethod
    def load(cls, path: str) -> "Skeleton":
        with open(path) as f:
            return cls.from_json(json.load(f))

    def __eq__(self, other) -> bool:
        return isinstance(other, Skeleton) and self.to_json() == other.to_json()


# ---------------------------------------------------------------------------
# Numeric composition at a point: the check that steps + expansion give the target
# ---------------------------------------------------------------------------

def compose_at(skel: Skeleton, basis: Basis, proof: Proof) -> list[Expanded]:
    """target = sum coef * left * input[idx] * right at the point `basis`/`proof` live at,
    obtained by pushing each step's multiplier through that element's expanded proof.  The
    steps are re-identified by key, so this works at the anchor and at any other point
    carrying the skeleton.  The multipliers are the anchor's; at another point they are only
    right if the lifted multipliers specialise to them (compose.py's concern, not this one's)."""
    bucket: dict[tuple[Word, int, Word], Fraction] = {}
    for s, k in zip(skel.steps, skel.step_indices(basis)):
        u, v, mult = s["left"], s["right"], s["mult"]
        for (left, idx, right), coef in proof.expand_dict(k).items():
            key = (u + left, idx, right + v)
            nv = bucket.get(key, Fraction(0)) + mult * coef
            if nv == 0:
                bucket.pop(key, None)
            else:
                bucket[key] = nv
    return [(c, l, i, r) for (l, i, r), c in bucket.items()]


def evaluate(terms: list[Expanded], inputs: list[Poly], char: int) -> Poly:
    """sum coef * left * inputs[idx] * right as a Poly -- compare with the target."""
    acc = Poly(char)
    for coef, u, idx, v in terms:
        lp = Poly(char); lp.add_term(u, Coef.one(char))
        rp = Poly(char); rp.add_term(v, Coef.one(char))
        acc = acc + (lp * inputs[idx] * rp).scalar_mul(Coef(coef, char))
    return acc
