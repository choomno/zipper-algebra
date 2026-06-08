#!/usr/bin/env python3
"""
Generate an f4ncgb input file (msolve format) for the zipper presentation
specialized at chosen integer (q, t) over either Q (characteristic 0) or GF(p).

Usage:
    python3 zipper_to_ms.py --n 3 --q 5 --t 7 --char 2147483647 \
        --untwist bigelow --out zipper-3-2147483647-5-7.ms

    # Or over Q with char=0:
    python3 zipper_to_ms.py --n 3 --q 5 --t 7 --char 0 --out zipper-3-Q-5-7.ms

The resulting .ms file feeds directly into f4ncgb:
    f4ncgb -p -o basis.txt zipper-3-...ms      # -p enables cofactor proof
"""

from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction
from typing import Iterable


# ---------------------------------------------------------------------------
# Coefficient arithmetic: either Fraction (char=0) or int mod p
# ---------------------------------------------------------------------------

class Coef:
    """Lightweight scalar wrapper supporting both Q and GF(p) backends."""
    __slots__ = ("v", "char")

    def __init__(self, v, char: int):
        self.char = char
        if char == 0:
            self.v = Fraction(v)
        else:
            self.v = int(v) % char

    @staticmethod
    def zero(char: int) -> "Coef":
        return Coef(0, char)

    @staticmethod
    def one(char: int) -> "Coef":
        return Coef(1, char)

    def __add__(self, other: "Coef") -> "Coef":
        return Coef(self.v + other.v, self.char)

    def __sub__(self, other: "Coef") -> "Coef":
        return Coef(self.v - other.v, self.char)

    def __mul__(self, other: "Coef") -> "Coef":
        return Coef(self.v * other.v, self.char)

    def __neg__(self) -> "Coef":
        return Coef(-self.v, self.char)

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Coef) and self.v == other.v and self.char == other.char

    def is_zero(self) -> bool:
        return self.v == 0

    def __str__(self) -> str:
        if self.char == 0:
            assert isinstance(self.v, Fraction)
            return (str(self.v.numerator) if self.v.denominator == 1
                    else f"{self.v.numerator}/{self.v.denominator}")
        return str(self.v)


# ---------------------------------------------------------------------------
# Non-commutative polynomial: dict mapping monomial-tuple -> Coef
# ---------------------------------------------------------------------------

class Poly:
    """A non-commutative polynomial. monomials are tuples of variable indices
    (1-based). The empty tuple () is the identity. Coefficients drop on collision."""
    __slots__ = ("terms", "char")

    def __init__(self, char: int):
        self.terms: dict[tuple[int, ...], Coef] = {}
        self.char = char

    def add_term(self, mon: tuple[int, ...], coef: Coef) -> None:
        if coef.is_zero():
            return
        existing = self.terms.get(mon)
        if existing is None:
            self.terms[mon] = coef
        else:
            new = existing + coef
            if new.is_zero():
                del self.terms[mon]
            else:
                self.terms[mon] = new

    def __add__(self, other: "Poly") -> "Poly":
        result = Poly(self.char)
        for mon, c in self.terms.items():
            result.add_term(mon, c)
        for mon, c in other.terms.items():
            result.add_term(mon, c)
        return result

    def __sub__(self, other: "Poly") -> "Poly":
        result = Poly(self.char)
        for mon, c in self.terms.items():
            result.add_term(mon, c)
        for mon, c in other.terms.items():
            result.add_term(mon, -c)
        return result

    def __neg__(self) -> "Poly":
        result = Poly(self.char)
        for mon, c in self.terms.items():
            result.add_term(mon, -c)
        return result

    def __mul__(self, other: "Poly") -> "Poly":
        result = Poly(self.char)
        for m1, c1 in self.terms.items():
            for m2, c2 in other.terms.items():
                result.add_term(m1 + m2, c1 * c2)
        return result

    def scalar_mul(self, c: Coef) -> "Poly":
        result = Poly(self.char)
        for mon, coef in self.terms.items():
            result.add_term(mon, coef * c)
        return result

    @staticmethod
    def scalar(c: Coef) -> "Poly":
        p = Poly(c.char)
        p.add_term((), c)
        return p

    @staticmethod
    def generator(idx: int, char: int) -> "Poly":
        p = Poly(char)
        p.add_term((idx,), Coef.one(char))
        return p

    def to_msolve_string(self, var_names: list[str]) -> str:
        """Render in msolve syntax. Variable names are the strings the user
        chose (e.g. 'x1', 'y2'); 'mon' tuples are 1-based indices into them."""
        if not self.terms:
            return "0"
        # Sort terms: longer monomials first, then lex.
        items = sorted(self.terms.items(), key=lambda kv: (-len(kv[0]), kv[0]))
        parts: list[str] = []
        for mon, coef in items:
            cs = str(coef)
            # mon is a tuple of 1-based var indices
            mon_str = "*".join(var_names[i - 1] for i in mon) if mon else ""
            if mon_str == "":
                term = cs
            elif cs == "1":
                term = mon_str
            elif cs == "-1":
                term = "-" + mon_str
            else:
                term = f"{cs}*{mon_str}"
            # Connect with sign
            if parts and not term.startswith("-"):
                parts.append("+" + term)
            else:
                parts.append(term)
        return "".join(parts)


# ---------------------------------------------------------------------------
# build_relations: mirrors zipper.gap's function of the same name
# ---------------------------------------------------------------------------

def build_relations(n: int, q_int: int, t_int: int, char: int, untw: str,
                    framing_power: int = 1,
                    theta_power: int = 1
                    ) -> tuple[list[str], list[Poly],
                               dict[tuple[int, int], Poly]]:
    """Generate variable names and the full list of input relations.

    Variable order: x1, x2, ..., x_{n-1}, y1, ..., y_{n-1}. 1-based indices
    in monomials map to this order.

    The ``framing_power`` parameter only affects the level-2 untwisting
    relation in ``bigelow`` mode: with framing_power = k >= 1, the level-2
    relation becomes (x_j^k - t) * S[2][j], i.e. the eigenvalue equation
    sigma_j^k * X_{2,j} = t * X_{2,j} (with t fixed, not t^k). k=1 is the
    canonical Bigelow framing. For k > 1 the slider need not be a sigma_j
    eigenvector and the resulting algebra is strictly larger than BMW.
    """
    var_names = [f"x{i}" for i in range(1, n)] + [f"y{i}" for i in range(1, n)]

    # Helpers
    one_p = Poly.scalar(Coef.one(char))

    def x(i: int) -> Poly:
        return Poly.generator(i, char)  # x1 is index 1, x2 is 2, ...

    def y(i: int) -> Poly:
        return Poly.generator((n - 1) + i, char)  # y1 is index n, y2 is n+1, ...

    q = Poly.scalar(Coef(q_int, char))
    t = Poly.scalar(Coef(t_int, char))

    relations: list[Poly] = []

    # Invertibility: x_i * y_i = 1, y_i * x_i = 1
    for i in range(1, n):
        relations.append(x(i) * y(i) - one_p)
        relations.append(y(i) * x(i) - one_p)

    # Braid + far-commutativity + mixed-sign length-3
    for i in range(1, n - 1):
        for j in range(i + 1, n):
            if j - i >= 2:
                relations.append(x(i) * x(j) - x(j) * x(i))
                relations.append(y(i) * y(j) - y(j) * y(i))
                relations.append(y(i) * x(j) - x(j) * y(i))
                relations.append(x(i) * y(j) - y(j) * x(i))
            else:  # j == i+1
                relations.append(x(i) * x(j) * x(i) - x(j) * x(i) * x(j))
                relations.append(y(i) * y(j) * y(i) - y(j) * y(i) * y(j))
    for i in range(1, n - 1):
        relations.append(x(i) * x(i + 1) * y(i) - y(i + 1) * x(i) * x(i + 1))
        relations.append(y(i) * y(i + 1) * x(i) - x(i + 1) * y(i) * y(i + 1))
        relations.append(x(i) * y(i + 1) * y(i) - y(i + 1) * y(i) * x(i + 1))
        relations.append(y(i) * x(i + 1) * x(i) - x(i + 1) * x(i) * y(i + 1))

    # Sliders S[i][j], for i in [2..n], j in [1..n-i+1]
    S: dict[tuple[int, int], Poly] = {}
    for j in range(1, n):
        S[(2, j)] = q * y(j) + (one_p - q) * one_p - x(j)
        # Note: (one_p - q) * one_p folds (1-q)*1; cleaner just to compute (1-q) as scalar
    # Simpler: redo S[2][j] with scalar (1-q)
    for j in range(1, n):
        S[(2, j)] = q * y(j) + Poly.scalar(Coef.one(char) - Coef(q_int, char)) - x(j)
    for i in range(3, n + 1):
        for j in range(1, n - i + 2):
            yProd = y(j)
            for k in range(j + 1, j + i - 1):
                yProd = yProd * y(k)
            xProd = x(j)
            for k in range(j + 1, j + i - 1):
                xProd = xProd * x(k)
            q_pow_v = Coef(1, char)
            for _ in range(i - 1):
                q_pow_v = q_pow_v * Coef(q_int, char)
            q_pow = Poly.scalar(q_pow_v)
            S[(i, j)] = (q_pow * yProd - xProd) * S[(i - 1, j)]

    # Zipper relations
    for j in range(1, n - 1):
        rel = (q * y(j + 1) + Poly.scalar(Coef.one(char) - Coef(q_int, char)) - x(j + 1)
               - q * y(j) * y(j + 1) + x(j) * x(j + 1)) * S[(2, j)]
        relations.append(rel)
    for i in range(3, n):
        for j in range(1, n - i + 1):
            yProd_after = y(j + 1)
            for k in range(j + 2, j + i):
                yProd_after = yProd_after * y(k)
            yProd_full = y(j) * yProd_after
            xProd_after = x(j + 1)
            for k in range(j + 2, j + i):
                xProd_after = xProd_after * x(k)
            xProd_full = x(j) * xProd_after
            q_pow_v = Coef(1, char)
            for _ in range(i - 1):
                q_pow_v = q_pow_v * Coef(q_int, char)
            q_pow = Poly.scalar(q_pow_v)
            rel = (q_pow * yProd_after - xProd_after
                   - q_pow * yProd_full + xProd_full) * S[(i, j)]
            relations.append(rel)

    # Untwisting relations
    if untw != "none":
        # 2-strand
        for j in range(1, n):
            if untw == "silly_question":
                # silly: (sigma_j^2 - t^2) X_{2,j} = 0 (squares-of-both-sides
                # form). Distinct from the framing_power variant below.
                rel = (x(j) * x(j) - t * t) * S[(2, j)]
            elif untw == "jordan":
                # (sigma_j - t)^2 X_{2,j} = (sigma_j^2 - 2t sigma_j + t^2) X_{2,j} = 0.
                # This makes sigma_j have eigenvalue t with multiplicity 2 on the
                # slider line; a Jordan block forms (eigenvector + generalized
                # eigenvector), forcing the algebra to be NON-semisimple.
                two_t = t + t
                rel = (x(j) * x(j) - two_t * x(j) + t * t) * S[(2, j)]
            else:
                # bigelow with framing_power k:
                #   (sigma_j^k - t) X_{2,j} = 0  i.e. sigma_j^k * X = t * X.
                # k=1 is the canonical Bigelow framing; k>=2 is a weaker
                # eigenvalue equation (k-th-power eigenvalue is t, but the
                # k-th roots of t are all candidate sigma_j eigenvalues).
                if framing_power < 1:
                    raise ValueError("framing_power must be >= 1")
                lhs = x(j)
                for _ in range(framing_power - 1):
                    lhs = lhs * x(j)
                rel = (lhs - t) * S[(2, j)]
            relations.append(rel)
        # k-strand (k >= 3)
        for i in range(3, n + 1):
            for j in range(1, n - i + 2):
                xProd_rev = x(j + i - 2)
                yProd_rev = y(j + i - 2)
                for k in range(2, i):
                    xProd_rev = xProd_rev * x(j + i - 1 - k)
                    yProd_rev = yProd_rev * y(j + i - 1 - k)
                q_pow_v = Coef(1, char)
                for _ in range(i - 1):
                    q_pow_v = q_pow_v * Coef(q_int, char)
                q_pow = Poly.scalar(q_pow_v)
                if untw == "indeterminate":
                    raise NotImplementedError("indeterminate mode needs an extra "
                                              "parameter s; not wired up here")
                # theta_power = m: replace canonical (Theta_i - q^{i-1}) X_i = 0
                # with (Theta_i^m - q^{(i-1)m}) X_i = 0.
                # Per sliders_and_full_twists.tex:
                #   - k_i = 1 at i = 3 is FORCED by rigidity of the slider module M
                #     when paired with canonical framing k_2 = 1.
                #   - For higher framing k_2 = a, the natural pairing at i = 3 is
                #     k_3 = a, giving (Theta_3^a - q^{2a}) X_3 = 0.
                #   - For i >= 4 the "right" k_i is conjectural; we keep canonical
                #     k_i = 1 (per the paper's note: "forced at i = 3 by rigidity
                #     but only conjectural for i > 3").
                # We therefore apply theta_power ONLY at i = 3; at i >= 4 the
                # relation is canonical regardless of theta_power.
                effective_power = theta_power if i == 3 else 1
                if effective_power < 1:
                    raise ValueError("theta_power must be >= 1")
                if effective_power == 1:
                    rel = (xProd_rev - q_pow * yProd_rev) * S[(i, j)]
                else:
                    # Theta_{i,j} as a positive word in x's:
                    #   Theta = (ascending Coxeter on band {j..j+i-2})
                    #         * (descending Coxeter on band {j..j+i-2})
                    #         = (x_j x_{j+1} ... x_{j+i-2}) * (x_{j+i-2} ... x_j)
                    # Verified equivalent to canonical k=1 relation:
                    #   in the algebra y_desc^{-1} = x_asc, so y_desc^{-1} * x_desc
                    #   = Theta, and (x_desc - q^{i-1} y_desc) S = 0 multiplies on
                    #   the left by y_desc^{-1} to give (Theta - q^{i-1}) S = 0.
                    # For k>=2 the (x_desc^k - q^{...} y_desc^k) form does NOT
                    # represent Theta^k = (x_asc x_desc)^k S, since x_desc and
                    # y_desc do not commute. So we build the genuine Theta word.
                    xProd_asc = x(j)
                    for kk in range(j + 1, j + i - 1):
                        xProd_asc = xProd_asc * x(kk)
                    theta_word = xProd_asc * xProd_rev
                    thetaPow = theta_word
                    q_pow_acc = q_pow
                    for _ in range(effective_power - 1):
                        thetaPow = thetaPow * theta_word
                        q_pow_acc = q_pow_acc * q_pow
                    rel = (thetaPow - q_pow_acc) * S[(i, j)]
                relations.append(rel)

    return var_names, relations, S


def append_zero_sliders(relations: list[Poly], S: dict[tuple[int, int], Poly],
                        zero_from: int, n: int) -> int:
    """For each i in [zero_from..n] and each j in [1..n-i+1], add the relation
    S[i][j] = 0 to `relations` (forcing the i-strand sliders to vanish).
    Returns the number of relations appended."""
    added = 0
    for i in range(zero_from, n + 1):
        for j in range(1, n - i + 2):
            if (i, j) in S:
                relations.append(S[(i, j)])
                added += 1
    return added


def append_target(relations: list[Poly], S: dict[tuple[int, int], Poly],
                  target_k: int, target_j: int = 1) -> int:
    """Append S[target_k][target_j] as the LAST input relation. Returns the
    0-based index of the appended target. After running f4ncgb with -e on,
    the corresponding proof line will encode the target's cofactor
    representation over the original input relations."""
    if (target_k, target_j) not in S:
        raise ValueError(f"S[{target_k}][{target_j}] not in slider dict; "
                         f"check that target_k >= 2 and "
                         f"1 <= target_j <= n - target_k + 1.")
    relations.append(S[(target_k, target_j)])
    return len(relations) - 1


# ---------------------------------------------------------------------------
# Main: emit msolve format
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, required=True, help="number of strands")
    ap.add_argument("--q", type=int, default=5, help="q specialization (integer)")
    ap.add_argument("--t", type=int, default=7, help="t specialization (integer)")
    ap.add_argument("--char", type=int, default=0,
                    help="coefficient characteristic (0 for Q, prime for GF(p))")
    ap.add_argument("--untwist", choices=["bigelow", "silly_question", "jordan", "none"],
                    default="bigelow",
                    help="'bigelow': (sigma^a - t)X_2 = 0 (the higher-framing variant). "
                         "'silly_question': (sigma^2 - t^2)X_2 = 0. "
                         "'jordan': (sigma - t)^2 X_2 = 0, i.e. sigma has eigenvalue t with multiplicity 2 "
                         "on the slider; INDUCES NON-SEMISIMPLICITY (a deliberate non-rigid framing). "
                         "'none': no level-2 framing relation at all.")
    ap.add_argument("--framing-power", type=int, default=1, metavar="K",
                    help="In 'bigelow' mode only: use the level-2 relation "
                         "(sigma_j^K - t) * S[2][j] = 0 instead of the "
                         "canonical K=1 case. K=1 is canonical; K>=2 gives "
                         "an algebra strictly larger than BMW.")
    ap.add_argument("--theta-power", type=int, default=1, metavar="M",
                    help="For i >= 3 untwisting (any 'untwist' mode): use "
                         "(Theta_i^M - q^{(i-1)M}) * S[i][j] = 0 in place of "
                         "the canonical M=1 relation. M=1 is canonical. NOTE: "
                         "k_3 = 1 is forced by rigidity per sliders_and_full_twists.tex "
                         "(rigidity of the slider module); M >= 2 at i = 3 gives a "
                         "weaker relation that is not consistent with the natural "
                         "presentation. Recommend leaving at 1.")
    ap.add_argument("--zero-from", type=int, default=None, metavar="K",
                    help="Force X_i = 0 for all i >= K by appending S[i][j] = 0 "
                         "as extra relations for every valid (i,j). Use to test "
                         "whether forcing a level-K vanishing makes the algebra "
                         "finite-dimensional (e.g., --zero-from 4 forces X_4 = 0).")
    ap.add_argument("--out", required=True, help="output .ms path")
    ap.add_argument("--target", type=int, default=None, metavar="K",
                    help="append the K-strand slider S[K][1] as the last input. "
                         "Its proof line gives the cofactor of X_K = S[K][1] in "
                         "terms of the other inputs. K must satisfy "
                         "2 <= K <= n. Common choices: K=3 (slider on first "
                         "three strands of an n-strand algebra).")
    args = ap.parse_args()

    var_names, relations, S = build_relations(args.n, args.q, args.t, args.char,
                                              args.untwist,
                                              framing_power=args.framing_power,
                                              theta_power=args.theta_power)
    if args.zero_from is not None:
        added = append_zero_sliders(relations, S, args.zero_from, args.n)
        print(f"Appended {added} zero-slider relations for i >= {args.zero_from}")
    if args.target is not None:
        target_idx = append_target(relations, S, args.target)
        print(f"Appended S[{args.target}][1] as input #{target_idx} "
              f"(0-based); look at line {target_idx + 1} of the proof.")

    with open(args.out, "w") as f:
        f.write(",".join(var_names) + "\n")
        f.write(str(args.char) + "\n")
        for i, rel in enumerate(relations):
            sep = "," if i < len(relations) - 1 else ""
            f.write(rel.to_msolve_string(var_names) + sep + "\n")

    print(f"Wrote {len(relations)} relations in {len(var_names)} variables to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
