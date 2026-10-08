#!/usr/bin/env python3
"""
Extract a cofactor representation of X_3 = S[3][1] over the n=4 Bigelow
zipper algebra at (q,t) = (5,7) over Q.

Reads:
  logs_and_traces/f4ncgb/basis-n4-Q.txt        (137-element GB of the 31-input ideal)
  logs_and_traces/f4ncgb/proof-n4-Q-noex.txt   (non-expanded cofactors:
                                                each B_k = sum c * left * (B_j | I_j) * right)

The expanded-proof option (-e on) of f4ncgb 1.0.1 produces cofactors that
agree with the basis only modulo the ideal in many cases; the
non-expanded variant is exact, so we re-expand the proofs ourselves.

Builds X_3 with the same Poly machinery used in zipper_in_ms.py, reduces
it modulo the Gröbner basis tracking every (left, basis_idx, right, coef)
reduction step, then composes the steps with the (self-expanded)
non-expanded proof to express X_3 as a finite Z-linear combination of
words

    sum_t  coef_t * (left_mon_t) * (input_{i_t}) * (right_mon_t).

Verifies by substituting the 31 input polynomials back in and checking the
result equals X_3 exactly.

Writes a human-readable file to
  logs_and_traces/f4ncgb/cofactor-X3-n4-Q.txt

Run from the repo root:
    python3 tools/extract_cofactor.py
"""

from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py lives
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)


import re
import time
from fractions import Fraction


# Make the sibling module importable when invoked as a script.
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PARENT_DIR = os.path.dirname(_SCRIPT_DIR)  # tools/, where zipper_in_ms.py lives
for _p in (_SCRIPT_DIR, _PARENT_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from zipper_in_ms import Coef, Poly, build_relations  # noqa: E402


# ---------------------------------------------------------------------------
# Variable / monomial order matching f4ncgb's output:
#     Monomial order:    x1 < x2 < x3 < y1 < y2 < y3
# Our 1-based encoding (x1=1, x2=2, x3=3, y1=4, y2=5, y3=6) lines up.
# ---------------------------------------------------------------------------

N = 4                       # default; configure(n) rebinds these for other n
NUM_VARS = 2 * (N - 1)
VAR_NAMES = [f"x{i}" for i in range(1, N)] + [f"y{i}" for i in range(1, N)]
VAR_INDEX = {name: idx + 1 for idx, name in enumerate(VAR_NAMES)}   # 1-based


def configure(n: int) -> None:
    """Rebind the strand-dependent globals. Call before anything else."""
    global N, NUM_VARS, VAR_NAMES, VAR_INDEX
    N = n
    NUM_VARS = 2 * (N - 1)
    VAR_NAMES = [f"x{i}" for i in range(1, N)] + [f"y{i}" for i in range(1, N)]
    VAR_INDEX = {name: idx + 1 for idx, name in enumerate(VAR_NAMES)}


def mon_sort_key(mon: tuple[int, ...]) -> tuple:
    """Degree-lex key with x1 < x2 < x3 < y1 < y2 < y3 (max = leading)."""
    return (len(mon),) + tuple(mon)


def leading_term(p: Poly) -> tuple[tuple[int, ...], Coef]:
    assert p.terms, "leading_term called on zero polynomial"
    lm = max(p.terms.keys(), key=mon_sort_key)
    return lm, p.terms[lm]


# ---------------------------------------------------------------------------
# Tokenizers / parsers
# ---------------------------------------------------------------------------

_FRAC_RE = re.compile(r"^(-?\d+(?:/\d+)?)$")


def _parse_coef_str(s: str) -> Fraction:
    if "/" in s:
        num, den = s.split("/")
        return Fraction(int(num), int(den))
    return Fraction(int(s))


def _split_terms(line: str) -> list[str]:
    """Split a polynomial line into signed term tokens.

    Tokens are sometimes joined to a leading sign (e.g. '-4/35*[0]') and
    sometimes separated by a bare ' + ' or ' - '. Strategy: whitespace-split
    and re-attach pending signs.
    """
    line = line.strip().rstrip(",").strip()
    if not line:
        return []
    raw = line.split()
    out: list[str] = []
    pending: str | None = None
    for tok in raw:
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
                inner = tok[1:]
                sign = "+" if combined.count("-") % 2 == 0 else "-"
                out.append(sign + inner)
            else:
                out.append(pending + tok)
            pending = None
    return out


def _parse_basis_term(tok: str) -> tuple[Fraction, list[str]]:
    """Parse a basis term (no [..] brackets). Returns (coef, var_names)."""
    sign_is_neg = (tok[0] == "-")
    body = tok[1:]
    pieces = body.split("*")
    coef = Fraction(1)
    if _FRAC_RE.match(pieces[0]):
        coef = _parse_coef_str(pieces[0])
        rest = pieces[1:]
    else:
        rest = pieces
    if sign_is_neg:
        coef = -coef
    for p in rest:
        if p not in VAR_INDEX:
            raise ValueError(f"unknown variable {p!r} in basis term {tok!r}")
    return coef, rest


def parse_basis(path: str, char: int) -> list[Poly]:
    polys: list[Poly] = []
    with open(path) as f:
        for raw in f:
            line = raw.strip().rstrip(",")
            if not line:
                continue
            p = Poly(char)
            for tok in _split_terms(line):
                coef, factors = _parse_basis_term(tok)
                mon = tuple(VAR_INDEX[v] for v in factors)
                p.add_term(mon, Coef(coef, char))
            polys.append(p)
    return polys


def _parse_noex_term(tok: str) -> tuple[Fraction, list[str], str, int, list[str]]:
    """Parse one term of the non-expanded proof. Returns
       (coef, left_names, kind, idx, right_names) where kind is
       'input' (the [ik] form) or 'basis' (the [k] form)."""
    sign_is_neg = (tok[0] == "-")
    body = tok[1:]
    pieces = body.split("*")
    coef = Fraction(1)
    if _FRAC_RE.match(pieces[0]):
        coef = _parse_coef_str(pieces[0])
        rest = pieces[1:]
    else:
        rest = pieces
    if sign_is_neg:
        coef = -coef
    left: list[str] = []
    right: list[str] = []
    idx: int | None = None
    kind: str | None = None
    for p in rest:
        if p.startswith("[i") and p.endswith("]"):
            if idx is not None:
                raise ValueError(f"two indices in {tok!r}")
            idx = int(p[2:-1])
            kind = "input"
        elif p.startswith("[") and p.endswith("]"):
            if idx is not None:
                raise ValueError(f"two indices in {tok!r}")
            idx = int(p[1:-1])
            kind = "basis"
        else:
            if p not in VAR_INDEX:
                raise ValueError(f"unknown variable {p!r} in {tok!r}")
            if idx is None:
                left.append(p)
            else:
                right.append(p)
    if idx is None:
        raise ValueError(f"no [..] index in {tok!r}")
    return coef, left, kind, idx, right


def parse_noex_proof(path: str
                     ) -> list[list[tuple[Fraction, list[str], str, int, list[str]]]]:
    out: list[list[tuple[Fraction, list[str], str, int, list[str]]]] = []
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line:
                out.append([])
                continue
            terms: list[tuple[Fraction, list[str], str, int, list[str]]] = []
            for tok in _split_terms(line):
                terms.append(_parse_noex_term(tok))
            out.append(terms)
    return out


# ---------------------------------------------------------------------------
# Recursively expand non-expanded proofs into pure-input cofactors.
# Each entry of the expanded result is (coef, left_mon_tuple,
# input_idx, right_mon_tuple). The expansion is done iteratively (smallest
# basis index first), so when we process basis[k] each referenced basis
# index j < k has already been fully expanded.
#
# Empirically, basis[k]'s non-expanded proof only references basis indices
# j < k (and inputs). If that ever fails we fall back to recursion with
# memoisation.
# ---------------------------------------------------------------------------

def _names_to_mon(names: list[str]) -> tuple[int, ...]:
    return tuple(VAR_INDEX[n] for n in names)


def expand_proofs(noex: list[list[tuple[Fraction, list[str], str, int, list[str]]]]
                   ) -> list[list[tuple[Fraction, tuple[int, ...], int, tuple[int, ...]]]]:
    """Convert non-expanded proofs into pure-input cofactors.

    Returns expanded[k] = list of (coef, left_mon, input_idx, right_mon)
    such that basis[k] = sum_t coef_t * left_mon_t * input_{input_idx_t}
                          * right_mon_t (exact polynomial identity).
    """
    expanded: list[list[tuple[Fraction, tuple[int, ...], int, tuple[int, ...]]]] = []
    for k in range(len(noex)):
        bucket: dict[tuple[tuple[int, ...], int, tuple[int, ...]], Fraction] = {}
        for coef, left, kind, idx, right in noex[k]:
            l_mon = _names_to_mon(left)
            r_mon = _names_to_mon(right)
            if kind == "input":
                key = (l_mon, idx, r_mon)
                cur = bucket.get(key, Fraction(0))
                nv = cur + coef
                if nv == 0:
                    if key in bucket:
                        del bucket[key]
                else:
                    bucket[key] = nv
            else:  # 'basis'
                if idx >= k:
                    raise NotImplementedError(
                        f"basis[{k}] references basis[{idx}] forward — "
                        f"add recursion with memoisation")
                # Substitute basis[idx]'s already-expanded cofactor.
                for sub_coef, sub_l, sub_i, sub_r in expanded[idx]:
                    new_left = l_mon + sub_l
                    new_right = sub_r + r_mon
                    new_coef = coef * sub_coef
                    key = (new_left, sub_i, new_right)
                    cur = bucket.get(key, Fraction(0))
                    nv = cur + new_coef
                    if nv == 0:
                        if key in bucket:
                            del bucket[key]
                    else:
                        bucket[key] = nv
        expanded.append([(c, l, i, r) for (l, i, r), c in bucket.items()])
    return expanded


# ---------------------------------------------------------------------------
# Non-commutative leading-monomial reduction
# ---------------------------------------------------------------------------

def find_substring(mu: tuple[int, ...], lam: tuple[int, ...]
                   ) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    n, m = len(mu), len(lam)
    if m == 0 or m > n:
        return None
    for i in range(n - m + 1):
        if mu[i:i + m] == lam:
            return mu[:i], mu[i + m:]
    return None


def reduce_against_basis(target: Poly, basis: list[Poly], char: int
                          ) -> tuple[Poly, list[tuple[tuple[int, ...], int,
                                                      tuple[int, ...], Coef]]]:
    """Leading-monomial-driven non-commutative reduction. Returns
    (remainder, steps). Steps record subtractions of the form
    mult * u * B_k * v from the running poly."""
    lt_basis = []
    for k, g in enumerate(basis):
        if not g.terms:
            continue
        lm, lc = leading_term(g)
        lt_basis.append((lm, lc, k))
    # Canonical tiebreak: lex on (length, lm, basis_index). For n=4 every
    # leading monomial in the 137-element basis is distinct, so the k field
    # never breaks a tie. We still spell it out so the canonical choice is
    # documented and stable under future bases where LM collisions occur.
    lt_basis.sort(key=lambda t: (len(t[0]), t[0], t[2]))

    steps: list[tuple[tuple[int, ...], int, tuple[int, ...], Coef]] = []
    p = Poly(char)
    for mon, coef in target.terms.items():
        p.add_term(mon, coef)
    remainder = Poly(char)

    while p.terms:
        lm_p, lc_p = leading_term(p)
        reduced = False
        for lm_g, lc_g, k in lt_basis:
            split = find_substring(lm_p, lm_g)
            if split is None:
                continue
            u, v = split
            if char == 0:
                mult_v = Fraction(lc_p.v) / Fraction(lc_g.v)
                mult = Coef(mult_v, char)
            else:
                inv = pow(int(lc_g.v), -1, char)
                mult = Coef((int(lc_p.v) * inv) % char, char)
            g = basis[k]
            for mg, cg in g.terms.items():
                new_mon = u + mg + v
                new_coef = mult * cg
                p.add_term(new_mon, -new_coef)
            steps.append((u, k, v, mult))
            reduced = True
            break
        if not reduced:
            remainder.add_term(lm_p, lc_p)
            del p.terms[lm_p]
    return remainder, steps


# ---------------------------------------------------------------------------
# Compose reduction steps with the input-level cofactor of each basis
# polynomial.
# ---------------------------------------------------------------------------

def compose_cofactor(steps: list[tuple[tuple[int, ...], int, tuple[int, ...], Coef]],
                      expanded_proofs: list[list[tuple[Fraction,
                                                       tuple[int, ...],
                                                       int,
                                                       tuple[int, ...]]]],
                      char: int
                      ) -> list[tuple[Fraction, tuple[int, ...], int,
                                       tuple[int, ...]]]:
    bucket: dict[tuple[tuple[int, ...], int, tuple[int, ...]], Fraction] = {}
    for u, k, v, mult in steps:
        mult_q = Fraction(mult.v) if char == 0 else Fraction(int(mult.v))
        for coef_t, left_mon, idx, right_mon in expanded_proofs[k]:
            new_left = u + left_mon
            new_right = right_mon + v
            new_coef = mult_q * coef_t
            key = (new_left, idx, new_right)
            cur = bucket.get(key, Fraction(0))
            nv = cur + new_coef
            if nv == 0:
                if key in bucket:
                    del bucket[key]
            else:
                bucket[key] = nv
    return [(c, l, i, r) for (l, i, r), c in bucket.items()]


# ---------------------------------------------------------------------------
# Verification: substitute originals back in and check we get X_3.
# ---------------------------------------------------------------------------

def build_originals_and_target(char: int = 0, q_int: int = 5, t_int: int = 7,
                               target: tuple[int, int] = (3, 1),
                               untw: str = "bigelow",
                               max_untwist_level: int | None = None
                               ) -> tuple[list[Poly], Poly]:
    """Inputs of the presentation, and the target slider S[i][j].

    The target was X_3 = S[3][1] at n=4 until 2026-09-16; it is a parameter now,
    because "which element are we proving lies in the ideal" is exactly what
    changes between jobs.
    """
    from presentation import Presentation
    from zipper_in_ms import NumericRing
    var_names, relations, S, _tags = Presentation(
        N, untw, max_untwist_level=max_untwist_level).build(NumericRing(q_int, t_int, char))
    assert var_names == VAR_NAMES, (var_names, VAR_NAMES)
    if target not in S:
        raise SystemExit(f"target S{list(target)} not available at n={N}; "
                         f"have {sorted(S)[:8]}...")
    return relations, S[target]


# kept so older callers do not break
def build_originals_and_x3(char: int = 0, q_int: int = 5, t_int: int = 7):
    return build_originals_and_target(char, q_int, t_int, (3, 1))


def verify_cofactor(cofactor_terms: list[tuple[Fraction, tuple[int, ...],
                                                int, tuple[int, ...]]],
                     originals: list[Poly], x3: Poly, char: int) -> bool:
    acc = Poly(char)
    for coef, u, idx, v in cofactor_terms:
        coef_c = Coef(coef, char)
        u_p = Poly(char)
        u_p.add_term(u, Coef.one(char))
        v_p = Poly(char)
        v_p.add_term(v, Coef.one(char))
        prod = u_p * originals[idx] * v_p
        prod = prod.scalar_mul(coef_c)
        acc = acc + prod
    diff = acc - x3
    return not diff.terms


# ---------------------------------------------------------------------------
# Pretty output
# ---------------------------------------------------------------------------

def _mon_str(mon: tuple[int, ...]) -> str:
    if not mon:
        return ""
    return "*".join(VAR_NAMES[i - 1] for i in mon)


def _frac_str(f: Fraction) -> str:
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"


def write_cofactor(path: str,
                   cofactor_terms: list[tuple[Fraction, tuple[int, ...],
                                              int, tuple[int, ...]]],
                   originals: list[Poly]) -> None:
    terms = sorted(cofactor_terms,
                   key=lambda t: (t[2], len(t[1]), t[1], len(t[3]), t[3]))
    used = sorted({t[2] for t in terms})
    with open(path, "w") as f:
        f.write("# Cofactor representation of X_3 = S[3][1] over the n=4 "
                "Bigelow zipper algebra,\n")
        f.write("# specialised at (q, t) = (5, 7) over Q.\n")
        f.write("#\n")
        f.write("# X_3 = sum_t  c_t * (left_mon_t) * R_{i_t} * (right_mon_t)\n")
        f.write("# where R_0..R_30 are the 31 original input relations\n")
        f.write("# (output of zipper_in_ms.py with --n 4 --q 5 --t 7\n")
        f.write("# --char 0 --untwist bigelow, no --target).\n")
        f.write(f"#\n# Number of terms: {len(terms)}\n")
        f.write(f"# Distinct inputs used (0-based): {used}\n")
        f.write("#\n# Each line is:  signed_coef * left_mon * [input_index] * right_mon\n")
        f.write("# An empty left/right monomial is omitted; coef +/-1 is\n")
        f.write("# printed explicitly.\n")
        f.write("#\n")
        f.write("# --- Original input relations R_0 .. R_30 ---\n")
        for i, g in enumerate(originals):
            s = g.to_msolve_string(VAR_NAMES)
            f.write(f"# [{i}] = {s}\n")
        f.write("# --- End of input relations ---\n#\n")
        for coef, u, idx, v in terms:
            cs = _frac_str(coef)
            left = _mon_str(u)
            right = _mon_str(v)
            pieces: list[str] = []
            if coef > 0:
                pieces.append("+" + cs)
            else:
                pieces.append(cs)  # already includes '-'
            if left:
                pieces.append(left)
            pieces.append(f"[{idx}]")
            if right:
                pieces.append(right)
            f.write(" * ".join(pieces) + "\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Phase 2: reduce a target against the "
                                             "anchor basis and emit the skeleton.")
    ap.add_argument("--n", type=int, default=4, help="strand count")
    ap.add_argument("--target", default="3,1", metavar="i,j",
                    help="the slider S[i][j] to prove lies in the ideal (default 3,1 = X_3)")
    ap.add_argument("--q", type=int, default=5)
    ap.add_argument("--t", type=int, default=7)
    ap.add_argument("--char", type=int, default=0)
    ap.add_argument("--untwist", default="bigelow")
    ap.add_argument("--max-untwist-level", type=int, default=None)
    a = ap.parse_args()

    configure(a.n)
    char = a.char
    q_int, t_int = a.q, a.t
    ti, tj = (int(v) for v in a.target.split(","))
    tag = f"n{a.n}" + ("" if (ti, tj) == (3, 1) else f"-S{ti}_{tj}")

    # cofactor/ -> tools/ -> zipper-algebra/ . Was ".." (one level short), a
    # leftover from when this script lived one directory higher.
    repo_root = os.path.abspath(os.path.join(_SCRIPT_DIR, "..", ".."))
    # Defaults reproduce the historical filenames exactly at n=4, target (3,1).
    _L = os.path.join(repo_root, "logs_and_traces", "f4ncgb")
    basis_path = os.path.join(_L, f"basis-n{a.n}-Q.txt")
    noex_proof_path = os.path.join(_L, f"proof-n{a.n}-Q-noex.txt")
    out_path = os.path.join(_L, (f"cofactor-X3-n{a.n}-Q.txt" if (ti, tj) == (3, 1)
                                 else f"cofactor-S{ti}_{tj}-n{a.n}-Q.txt"))

    print(f"Parsing basis: {basis_path}")
    t0 = time.time()
    basis = parse_basis(basis_path, char)
    print(f"  parsed {len(basis)} basis polynomials in {time.time()-t0:.2f}s")

    print(f"Parsing non-expanded proof: {noex_proof_path}")
    t0 = time.time()
    noex = parse_noex_proof(noex_proof_path)
    print(f"  parsed {len(noex)} proof lines in {time.time()-t0:.2f}s")
    assert len(basis) == len(noex), (len(basis), len(noex))

    print("Recursively expanding non-expanded proofs into input cofactors...")
    t0 = time.time()
    expanded = expand_proofs(noex)
    expanded_sizes = [len(p) for p in expanded]
    print(f"  expansion took {time.time()-t0:.2f}s; "
          f"max expanded length = {max(expanded_sizes)}, "
          f"mean = {sum(expanded_sizes)/len(expanded_sizes):.1f}")

    print(f"Building originals and S[{ti}][{tj}] via zipper_in_ms.py...")
    originals, x3 = build_originals_and_target(char, q_int, t_int, (ti, tj),
                                               a.untwist, a.max_untwist_level)
    print(f"  S[{ti}][{tj}] has {len(x3.terms)} terms; "
          f"{len(originals)} original inputs")

    # Sanity check: verify each expanded[k] reconstructs basis[k].
    print("Sanity-checking that expanded cofactors reconstruct each basis poly...")
    t0 = time.time()
    bad = 0
    for k, exp in enumerate(expanded):
        acc = Poly(char)
        for coef, u, idx, v in exp:
            lp = Poly(char)
            lp.add_term(u, Coef.one(char))
            rp = Poly(char)
            rp.add_term(v, Coef.one(char))
            prod = (lp * originals[idx] * rp).scalar_mul(Coef(coef, char))
            acc = acc + prod
        if (acc - basis[k]).terms:
            bad += 1
            if bad <= 3:
                print(f"  basis[{k}] mismatch")
    print(f"  sanity check took {time.time()-t0:.2f}s; bad: {bad}/{len(basis)}")
    if bad:
        print("  FATAL: expansion produced cofactors that do not match basis polys.")
        return 1

    print(f"Reducing S[{ti}][{tj}] against the Gröbner basis...")
    t0 = time.time()
    remainder, steps = reduce_against_basis(x3, basis, char)
    print(f"  reduction took {time.time()-t0:.2f}s; "
          f"{len(steps)} steps; remainder terms: {len(remainder.terms)}")
    if remainder.terms:
        print("  ERROR: nonzero remainder; X_3 not in the ideal?!")

    # Emit the SKELETON as data, so downstream phases stop hard-coding it.
    # symbolic_lift_x3.py previously carried `transitive_deps(anchor, [16,17,102])`
    # as constants; those indices are a property of THIS target against THIS
    # anchor, and change with either. Writing them out makes the dependency
    # explicit and is a prerequisite for targets other than X_3.
    skel_path = os.path.join(_L, f"skeleton-{tag}-Q.json")
    try:
        import json as _json
        _skel = {
            "target": f"S[{ti}][{tj}]",
            "n": a.n,
            "anchor": {"q": q_int, "t": t_int, "char": char},
            "basis_size": len(basis),
            "n_steps": len(steps),
            "remainder_terms": len(remainder.terms),
            # one entry per one-step reduction: X -> X - mult * left * basis[k] * right
            "steps": [{"basis_index": k, "left": list(u), "right": list(v),
                       "mult": str(mult)} for (u, k, v, mult) in steps],
            "basis_indices": sorted({k for (_u, k, _v, _m) in steps}),
        }
        with open(skel_path, "w") as _f:
            _json.dump(_skel, _f, indent=1)
        print(f"  wrote skeleton: {skel_path} "
              f"(indices {_skel['basis_indices']})")
    except Exception as _e:
        print(f"  WARNING: could not write skeleton: {_e}")
        return 1

    print(f"Composing cofactor of S[{ti}][{tj}] over the {len(originals)} original inputs...")
    t0 = time.time()
    cofactor_terms = compose_cofactor(steps, expanded, char)
    print(f"  composition took {time.time()-t0:.2f}s; "
          f"{len(cofactor_terms)} collected terms")

    print("Verifying by substituting originals back...")
    t0 = time.time()
    ok = verify_cofactor(cofactor_terms, originals, x3, char)
    print(f"  verification took {time.time()-t0:.2f}s -> "
          f"{'OK' if ok else 'FAILED'}")
    if not ok:
        return 1

    print(f"Writing cofactor file: {out_path}")
    write_cofactor(out_path, cofactor_terms, originals)
    size = os.path.getsize(out_path)
    print(f"  wrote {size} bytes; "
          f"{len(cofactor_terms)} (left, input, right, coef) tuples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
