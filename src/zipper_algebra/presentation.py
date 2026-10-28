#!/usr/bin/env python3
"""Ordered defining data for the implemented zipper-algebra quotients.

``Presentation`` selects a strand count, untwisting preset and additional relations.
Its numeric and symbolic views use one relation builder. The content hash identifies
that ordered symbolic presentation, independently of parameter specialization.
Computational order, solver configuration and completion evidence belong to run records.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

from .coefficients import NumericRing, SymbolicRing
from .polynomial import Poly
from .relations import RelTag, build_relations_ring


# ---------------------------------------------------------------------------
# Extra relations that can be adjoined to a presentation, by name.
# Moved from papers/generic/generic_dim.py (where they were numeric-only).
# Each takes (S, ring, theta) and returns a Poly over the ring.
# ---------------------------------------------------------------------------

def _theta_word(i: int, j: int, ring, n: int) -> Poly:
    """Full twist on an i-strand block: thesis Theta_{i-1,j}, represented by the positive word."""
    def x(k):
        return Poly.generator_over(k, ring.one)
    asc = x(j)
    for kk in range(j + 1, j + i - 1):
        asc = asc * x(kk)
    desc = x(j + i - 2)
    for kk in range(j + i - 3, j - 1, -1):
        desc = desc * x(kk)
    return asc * desc


EXTRA_CATALOGUE = {
    "X_4":              lambda S, ring, theta: S[(4, 1)],
    "X_5":              lambda S, ring, theta: S[(5, 1)],
    "(Theta_3-q^2)X_3": lambda S, ring, theta: (theta(3, 1) - Poly.scalar(ring.q * ring.q)) * S[(3, 1)],
    "(Theta_4+t)X_4":   lambda S, ring, theta: (theta(4, 1) + Poly.scalar(ring.t)) * S[(4, 1)],
    "(Theta_5+t)X_5":   lambda S, ring, theta: (theta(5, 1) + Poly.scalar(ring.t)) * S[(5, 1)],
}


@dataclass(frozen=True)
class Presentation:
    n: int
    untwist: str = "bigelow"
    framing_power: int = 1
    theta_power: int = 1
    max_zipper_level: int | None = None
    max_untwist_level: int | None = None
    zero_from: int | None = None
    extra: tuple = field(default_factory=tuple)    # names from EXTRA_CATALOGUE, in order

    # --- construction ------------------------------------------------------
    def __post_init__(self):
        if type(self.n) is not int or self.n < 2:
            raise ValueError("n must be an integer at least 2")
        if self.untwist not in ("bigelow", "square-zero", "silly_question", "none"):
            raise ValueError(f"unknown untwisting choice {self.untwist!r}")
        for name in ("framing_power", "theta_power"):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        for name in ("max_zipper_level", "max_untwist_level", "zero_from"):
            value = getattr(self, name)
            if value is not None and (type(value) is not int or not 2 <= value <= self.n):
                raise ValueError(f"{name} must be an integer from 2 to n, or None")
        if not isinstance(self.extra, tuple) or any(not isinstance(name, str) for name in self.extra):
            raise ValueError("extra must be a tuple of relation names")
        if len(set(self.extra)) != len(self.extra):
            raise ValueError("extra relation names must be distinct")
        required_strands = {
            "X_4": 4, "X_5": 5, "(Theta_3-q^2)X_3": 3,
            "(Theta_4+t)X_4": 4, "(Theta_5+t)X_5": 5,
        }
        for name in self.extra:
            if name not in EXTRA_CATALOGUE:
                raise ValueError(f"unknown extra relation {name!r}")
            if self.n < required_strands[name]:
                raise ValueError(f"extra relation {name!r} requires at least {required_strands[name]} strands")

    def var_names(self) -> list[str]:
        return [f"x{i}" for i in range(1, self.n)] + [f"y{i}" for i in range(1, self.n)]

    def build(self, ring):
        """(var_names, relations, S, tags) over the given ring. Zero-slider relations and
        extras are appended after the core presentation, each with a tag."""
        names, rels, S, tags = build_relations_ring(
            self.n, ring, self.untwist, self.framing_power, self.theta_power,
            self.max_zipper_level, self.max_untwist_level)
        if self.zero_from is not None:
            for i in range(self.zero_from, self.n + 1):
                for j in range(1, self.n - i + 2):
                    if (i, j) in S:
                        rels.append(S[(i, j)])
                        tags.append(RelTag("ZERO", i, j))
        for name in self.extra:
            if name not in EXTRA_CATALOGUE:
                raise KeyError(f"unknown extra relation {name!r}; known: {sorted(EXTRA_CATALOGUE)}")
            theta = lambda i, j: _theta_word(i, j, ring, self.n)  # noqa: E731
            rels.append(EXTRA_CATALOGUE[name](S, ring, theta))
            tags.append(RelTag("EXTRA", name, None))
        return names, rels, S, tags

    def relations(self, ring):
        return self.build(ring)[1]

    def sliders(self, ring):
        return self.build(ring)[2]

    def target(self, i: int, j: int, ring) -> Poly:
        if (type(i) is not int or type(j) is not int
                or not 2 <= i <= self.n or not 1 <= j <= self.n - i + 1):
            raise ValueError("slider indices must be integers with 2 <= i <= n and 1 <= j <= n-i+1")
        S = self.sliders(ring)
        if (i, j) not in S:
            raise ValueError(f"S[{i}][{j}] does not exist at n={self.n}; have {sorted(S)}")
        return S[(i, j)]

    # --- identity ------------------------------------------------------------
    def labels(self) -> list[str]:
        """Readable relation labels in emission order, derived from the construction."""
        _, _, _, tags = self.build(NumericRing(2, 3, 0))
        return [t.label(self.n) if t.kind != "EXTRA" else f"EXTRA[{t.level}]" for t in tags]

    def n_relations(self) -> int:
        return len(self.labels())

    def slug(self) -> str:
        s = f"n{self.n}-{self.untwist}"
        if self.framing_power != 1:
            s += f"-a{self.framing_power}"
        if self.theta_power != 1:
            s += f"-k{self.theta_power}"
        if self.max_zipper_level is not None:
            s += f"-Z{self.max_zipper_level}"
        if self.max_untwist_level is not None:
            s += f"-L{self.max_untwist_level}"
        if self.zero_from is not None:
            s += f"-z{self.zero_from}"
        for e in self.extra:
            s += "-x" + "".join(ch if ch.isalnum() else "_" for ch in e)
        return s

    def canonical_text(self) -> str:
        """The symbolic relations rendered deterministically: the presentation's identity,
        independent of any specialization."""
        names, rels, _, tags = self.build(SymbolicRing())
        lines = []
        for rel, tag in zip(rels, tags):
            terms = sorted(((m, sorted(c.v.items())) for m, c in rel.terms.items()),
                           key=lambda kv: (len(kv[0]), kv[0]))
            body = " + ".join("*".join(names[i - 1] for i in m) + ":" +
                              ",".join(f"{i},{j}:{c}" for (i, j), c in cv)
                              for m, cv in terms)
            lines.append(f"{tag.label(self.n) if tag.kind != 'EXTRA' else tag.level} = {body}")
        return "\n".join(lines) + "\n"

    def content_hash(self) -> str:
        return hashlib.sha256(self.canonical_text().encode()).hexdigest()[:16]

    # --- views for the consumers ---------------------------------------------
    def emit_ms(self, path: str, q: int, t: int, char: int, *, provenance=True) -> int:
        """Write the msolve input exactly as zipper_in_ms.py's CLI does (same bytes, so the
        driver's content-addressed cache keys are unchanged). Returns the relation count."""
        if provenance:
            from .provenance import capture_provenance, write_output_provenance
            producer = capture_provenance()
        names, rels, _, _ = self.build(NumericRing(q, t, char))
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(",".join(names) + "\n")
            f.write(str(char) + "\n")
            for i, rel in enumerate(rels):
                sep = "," if i < len(rels) - 1 else ""
                f.write(rel.to_msolve_string(names) + sep + "\n")
        if provenance:
            write_output_provenance(path, producer=producer)
        return len(rels)

    def symbolic_fitted(self) -> list[dict]:
        """[{monomial: {(deg_q, deg_t): Fraction}}] per relation -- the format
        papers/generic/generic_dim.py consumed from the sampling fitter, now exact."""
        _, rels, _, _ = self.build(SymbolicRing())
        return [{m: dict(c.v) for m, c in rel.terms.items()} for rel in rels]

    def symbolic_ncpolys(self):
        """(relations, S) as symbolic_identity_check.NCPoly dicts {mon: {(i,j): Fraction}}."""
        _, rels, S, _ = self.build(SymbolicRing())
        conv = lambda p: {m: dict(c.v) for m, c in p.terms.items()}  # noqa: E731
        return [conv(r) for r in rels], {k: conv(v) for k, v in S.items()}


# The algebra the X_3 pipeline was built on: n=4, Bigelow untwisting, all levels (31 relations).
Presentation.DEFAULT = Presentation(4)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="describe a presentation")
    ap.add_argument("--n", type=int, default=4)
    ap.add_argument("--untwist", default="bigelow")
    ap.add_argument("--max-untwist-level", type=int, default=None)
    ap.add_argument("--max-zipper-level", type=int, default=None)
    ap.add_argument("--extra", action="append", default=[])
    a = ap.parse_args()
    p = Presentation(a.n, a.untwist, max_zipper_level=a.max_zipper_level,
                     max_untwist_level=a.max_untwist_level, extra=tuple(a.extra))
    print(p.slug(), p.content_hash(), p.n_relations(), "relations")
    for i, l in enumerate(p.labels()):
        print(f"  R_{i:<3} {l}")
