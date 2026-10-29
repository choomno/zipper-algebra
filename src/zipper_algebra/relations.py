"""Zipper presentation relations, sliders and their ordered relation labels."""
from __future__ import annotations

from .coefficients import NumericRing
from .polynomial import Poly

__all__ = ["RelTag", "build_relations", "build_relations_ring",
           "append_zero_sliders", "append_target"]


class RelTag:
    """What a relation IS: kind in {INV, BR, FC, MIX, Z, T}, level, anchor j, variant.
    Recorded by build_relations_ring at every append, so the readable label is derived from
    the construction rather than maintained in a second place (relations_naming.py used to
    hold a hand-typed n=4 table, which silently mislabels anything after an omitted level)."""
    __slots__ = ("kind", "level", "anchor", "variant")

    def __init__(self, kind: str, level, anchor, variant: str = ""):
        self.kind, self.level, self.anchor, self.variant = kind, level, anchor, variant

    def label(self, n: int) -> str:
        k, i, j, v = self.kind, self.level, self.anchor, self.variant
        if k == "INV":
            return f"INV_{i}_{v}"
        if k == "BR":
            return f"BR_{i}" + ("_y" if v == "y" else "")
        if k == "FC":
            # n=4 has a single far-commuting pair (1,3); its historical labels carry no indices.
            return f"FC_{v}" if n <= 4 else f"FC_{i}{j}_{v}"
        if k == "MIX":
            return f"MIX_{i}_{v}"
        if k in ("Z", "T"):
            return f"{k}_{i}" + (f"_j{j}" if j != 1 else "")
        return f"{k}_{i}_{j}{('_' + v) if v else ''}"

    def __repr__(self):
        return f"RelTag({self.kind},{self.level},{self.anchor},{self.variant!r})"


def build_relations(n: int, q_int: int, t_int: int, char: int, untw: str,
                    framing_power: int = 1,
                    theta_power: int = 1,
                    max_zipper_level: int | None = None,
                    max_untwist_level: int | None = None
                    ) -> tuple[list[str], list[Poly],
                               dict[tuple[int, int], Poly]]:
    """Classic entry point: q, t specialised to integers. See build_relations_ring."""
    names, rels, S, _tags = build_relations_ring(n, NumericRing(q_int, t_int, char), untw,
                                                 framing_power, theta_power,
                                                 max_zipper_level, max_untwist_level)
    return names, rels, S


def build_relations_ring(n: int, ring, untw: str,
                         framing_power: int = 1,
                         theta_power: int = 1,
                         max_zipper_level: int | None = None,
                         max_untwist_level: int | None = None
                         ) -> tuple[list[str], list[Poly],
                                    dict[tuple[int, int], Poly], list[RelTag]]:
    """Generate variable names, the input relations, the sliders S, and a RelTag per relation.

    `ring` is NumericRing (integers q,t over Q or GF(p)) or SymbolicRing (q,t transcendental).

    Variable order: x1, x2, ..., x_{n-1}, y1, ..., y_{n-1}. 1-based indices
    in monomials map to this order.

    The ``framing_power`` parameter only affects the level-2 untwisting
    relation in ``bigelow`` mode: with framing_power = k >= 1, the level-2
    relation becomes (x_j^k - t) * S[2][j], i.e. the eigenvalue equation
    sigma_j^k * X_{2,j} = t * X_{2,j} (with t fixed, not t^k). k=1 is the
    canonical Bigelow framing. For k > 1 this specifies polynomial annihilation
    by sigma_j^k - t; strict quotient comparisons require separate evidence.
    """
    var_names = [f"x{i}" for i in range(1, n)] + [f"y{i}" for i in range(1, n)]

    # Helpers
    one_p = Poly.scalar(ring.one)

    def x(i: int) -> Poly:
        return Poly.generator_over(i, ring.one)  # x1 is index 1, x2 is 2, ...

    def y(i: int) -> Poly:
        return Poly.generator_over((n - 1) + i, ring.one)  # y1 is index n, y2 is n+1, ...

    q = Poly.scalar(ring.q)
    t = Poly.scalar(ring.t)

    relations: list[Poly] = []
    tags: list[RelTag] = []

    def add(rel: Poly, tag: RelTag) -> None:
        relations.append(rel)
        tags.append(tag)

    # Invertibility: x_i * y_i = 1, y_i * x_i = 1
    for i in range(1, n):
        add(x(i) * y(i) - one_p, RelTag("INV", i, None, "p"))
        add(y(i) * x(i) - one_p, RelTag("INV", i, None, "m"))

    # Braid + far-commutativity + mixed-sign length-3
    for i in range(1, n - 1):
        for j in range(i + 1, n):
            if j - i >= 2:
                add(x(i) * x(j) - x(j) * x(i), RelTag("FC", i, j, "xx"))
                add(y(i) * y(j) - y(j) * y(i), RelTag("FC", i, j, "yy"))
                add(y(i) * x(j) - x(j) * y(i), RelTag("FC", i, j, "yx"))
                add(x(i) * y(j) - y(j) * x(i), RelTag("FC", i, j, "xy"))
            else:  # j == i+1
                add(x(i) * x(j) * x(i) - x(j) * x(i) * x(j), RelTag("BR", i, j, "x"))
                add(y(i) * y(j) * y(i) - y(j) * y(i) * y(j), RelTag("BR", i, j, "y"))
    for i in range(1, n - 1):
        add(x(i) * x(i + 1) * y(i) - y(i + 1) * x(i) * x(i + 1), RelTag("MIX", i, i + 1, "a"))
        add(y(i) * y(i + 1) * x(i) - x(i + 1) * y(i) * y(i + 1), RelTag("MIX", i, i + 1, "b"))
        add(x(i) * y(i + 1) * y(i) - y(i + 1) * y(i) * x(i + 1), RelTag("MIX", i, i + 1, "c"))
        add(y(i) * x(i + 1) * x(i) - x(i + 1) * x(i) * y(i + 1), RelTag("MIX", i, i + 1, "d"))

    # Sliders S[i][j], for i in [2..n], j in [1..n-i+1]
    S: dict[tuple[int, int], Poly] = {}
    for j in range(1, n):
        S[(2, j)] = q * y(j) + Poly.scalar(ring.one - ring.q) - x(j)
    for i in range(3, n + 1):
        for j in range(1, n - i + 2):
            yProd = y(j)
            for k in range(j + 1, j + i - 1):
                yProd = yProd * y(k)
            xProd = x(j)
            for k in range(j + 1, j + i - 1):
                xProd = xProd * x(k)
            q_pow_v = ring.one
            for _ in range(i - 1):
                q_pow_v = q_pow_v * ring.q
            q_pow = Poly.scalar(q_pow_v)
            S[(i, j)] = (q_pow * yProd - xProd) * S[(i - 1, j)]

    # Zipper relations
    for j in range(1, n - 1):
        rel = (q * y(j + 1) + Poly.scalar(ring.one - ring.q) - x(j + 1)
               - q * y(j) * y(j + 1) + x(j) * x(j + 1)) * S[(2, j)]
        add(rel, RelTag("Z", 2, j))
    zip_high = n - 1 if max_zipper_level is None else min(n - 1, max_zipper_level)
    for i in range(3, zip_high + 1):
        for j in range(1, n - i + 1):
            yProd_after = y(j + 1)
            for k in range(j + 2, j + i):
                yProd_after = yProd_after * y(k)
            yProd_full = y(j) * yProd_after
            xProd_after = x(j + 1)
            for k in range(j + 2, j + i):
                xProd_after = xProd_after * x(k)
            xProd_full = x(j) * xProd_after
            q_pow_v = ring.one
            for _ in range(i - 1):
                q_pow_v = q_pow_v * ring.q
            q_pow = Poly.scalar(q_pow_v)
            rel = (q_pow * yProd_after - xProd_after
                   - q_pow * yProd_full + xProd_full) * S[(i, j)]
            add(rel, RelTag("Z", i, j))

    # Untwisting relations
    if untw != "none":
        # 2-strand
        for j in range(1, n):
            if untw == "silly_question":
                # silly: (sigma_j^2 - t^2) X_{2,j} = 0 (squares-of-both-sides
                # form). Distinct from the framing_power variant below.
                rel = (x(j) * x(j) - t * t) * S[(2, j)]
            elif untw == "square-zero":
                # (sigma_j - t)^2 X_{2,j} = (sigma_j^2 - 2t sigma_j + t^2) X_{2,j} = 0.
                # This is a squared annihilating relation. It does not by itself
                # establish a nonzero Jordan block or nonsemisimplicity.
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
            add(rel, RelTag("T", 2, j, untw))
        # k-strand (k >= 3)
        untw_high = n if max_untwist_level is None else min(n, max_untwist_level)
        for i in range(3, untw_high + 1):
            for j in range(1, n - i + 2):
                xProd_rev = x(j + i - 2)
                yProd_rev = y(j + i - 2)
                for k in range(2, i):
                    xProd_rev = xProd_rev * x(j + i - 1 - k)
                    yProd_rev = yProd_rev * y(j + i - 1 - k)
                q_pow_v = ring.one
                for _ in range(i - 1):
                    q_pow_v = q_pow_v * ring.q
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
                add(rel, RelTag("T", i, j))

    return var_names, relations, S, tags


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
    """Append the slider as a defining relation of an augmented presentation.

    Return its zero-based input index. Membership in the augmented ideal does
    not establish membership in the ideal generated by the original inputs.
    """
    if (target_k, target_j) not in S:
        raise ValueError(f"S[{target_k}][{target_j}] not in slider dict; "
                         f"check that target_k >= 2 and "
                         f"1 <= target_j <= n - target_k + 1.")
    relations.append(S[(target_k, target_j)])
    return len(relations) - 1
