"""Frozen historical X3 implementation used for comparison tests.

Source: archive/cofactor-pre-consolidation/extract_cofactor.py at 0db1dd9.
The original file SHA-256 is 515c5f22a2b35ddca22d6fb463560a5adffd03da9708c2c2be794d321d79a3f2.
The parser, proof expansion, reduction and composition bodies are retained;
script path setup and the obsolete command-line driver are omitted.
This implementation shares Coef/Poly arithmetic with the maintained package.
It is a historical implementation comparison, not an independent arithmetic checker.
"""
from __future__ import annotations
import re
from fractions import Fraction
from zipper_algebra.zipper_in_ms import Coef, Poly, build_relations

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

