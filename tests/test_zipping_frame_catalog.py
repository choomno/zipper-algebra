"""Every higher-rank catalog branch preserves the two marked flype frames.

The finite table is a regression for the uniform tape formulas; it does not
assert scheduler reachability or replace the frame identities' braid proof.
"""

from zipper_algebra.zipping_algorithm.catalog import (
    Catalog,
    expand_polynomial,
    condition,
)
from zipper_algebra.zipping_algorithm.frames import transport_frame, check_frame


def inverse(w):
    return tuple(-x for x in reversed(w))


def half(t, j):
    return tuple(x for length in range(t - 1, 0, -1) for x in range(j, j + length))


def friend(t, j, k):
    return tuple(x for p in range(j + k, j, -1) for x in range(p, p + t))


def data(r, j, k):
    name = r["id"].split(".I")[0]
    mir = False
    a = b = 0
    if name.startswith(("zip.far.", "zip.close.")):
        return [(j, k, mir, a, b)]
    if name.startswith("zip.flype."):
        a = -1 if ".positive" in name else 1
        return [(j, 0, True, a, -a)]
    if name == "zip.nearest.positive":
        return [(j + 1, k - 1, False, 0, 0)]
    if name == "zip.nearest.negative.i_ge_3":
        return [(j, k, False, 0, 0), (j, k, False, 0, 0), (j + 1, k - 1, False, 0, 0)]
    if name == "zip.left.negative":
        return [(j - 1, k + 1, False, 0, 0)]
    if name == "zip.left.positive.i_ge_3":
        return [(j, k, False, 0, 0), (j, k, False, 0, 0), (j - 1, k + 1, False, 0, 0)]
    if name == "zip.farthest.negative":
        return [(j, k - 1, False, 0, 0)]
    if name == "zip.outer.positive":
        return [(j, k + 1, False, 0, 0)]
    if name == "zip.farthest.positive.i_ge_3":
        return [
            (j + k, 0, True, -1, 1),
            (j + k, 0, True, -1, 1),
            (j, k - 1, False, 0, 0),
        ]
    if name == "zip.outer.negative.i_ge_3":
        return [
            (j + k + 1, 0, False, 0, 0),
            (j + k + 1, 0, False, 0, 0),
            (j, k + 1, False, 0, 0),
        ]
    return None


def test_all_higher_catalog_branches_both_orientations_preserve_frames():
    cat = Catalog.bundled()
    counts = {}
    failures = []
    for n in range(3, 8):
        for i in range(3, n + 1):
            t = i - 1
            for j in range(1, n - i + 2):
                for k in range(n - i - j + 2):
                    for l in range(1, n):
                        env = dict(n=n, i=i, j=j, k=k, l=l)
                        for r in cat.rules:
                            if (
                                not r["id"].startswith("zip.")
                                or not r["id"].endswith(".I")
                                or not condition(r["condition"], env)
                            ):
                                continue
                            desc = data(r, j, k)
                            if desc is None:
                                continue
                            for reflected in (False, True):
                                rr = cat.get(r["id"][:-1] + "M") if reflected else r

                                def reflect(w):
                                    return (
                                        tuple(
                                            (1 if x > 0 else -1) * (n - abs(x))
                                            for x in w
                                        )
                                        if reflected
                                        else w
                                    )

                                def anchor(a):
                                    return n - a - t + 2 if reflected else a

                                J = n - i - j + 2 if reflected else j
                                lhs = expand_polynomial(rr["lhs"], env)
                                (((q, word), co),) = lhs.items()
                                core = reflect(friend(t, j, k)) + (
                                    ("slider", i, J, reflected, False, False),
                                )
                                if word[-len(core) :] != core:
                                    raise Exception("badcore", rr["id"], env)
                                W = word[: -len(core)]
                                try:
                                    before = transport_frame(
                                        W,
                                        n=n,
                                        size=t,
                                        start=anchor(j + k + 1),
                                        side="left",
                                    )
                                    check_frame(before)
                                    if before.transformed != W:
                                        raise Exception("changed incoming")
                                    for term, (jj, kk, mir, aa, bb) in zip(
                                        rr["rhs"], desc, strict=True
                                    ):
                                        poly = expand_polynomial(
                                            [
                                                {
                                                    "coefficient": [
                                                        {"integer": 1, "q_power": 0}
                                                    ],
                                                    "word": term["word"],
                                                }
                                            ],
                                            env,
                                        )
                                        (((qq, ww), cc),) = poly.items()
                                        actualj = n - i - jj + 2 if reflected else jj
                                        actualmir = not mir if reflected else mir
                                        slider = (
                                            "slider",
                                            i,
                                            actualj,
                                            actualmir,
                                            False,
                                            False,
                                        )
                                        ix = ww.index(slider)
                                        prefix, suffix = ww[:ix], ww[ix + 1 :]
                                        outanchor = anchor(jj if mir else jj + kk + 1)
                                        inanchor = anchor(jj + 1 if mir else jj)
                                        h = half(t, outanchor)
                                        leftframe = (
                                            inverse(h) if aa < 0 else h if aa else ()
                                        )
                                        ff = () if mir else reflect(friend(t, jj, kk))
                                        if prefix != leftframe + ff:
                                            raise Exception(
                                                "badprefix", prefix, leftframe + ff
                                            )
                                        if outanchor != before.end:
                                            raise Exception(
                                                "output tape mismatch",
                                                outanchor,
                                                before.end,
                                            )
                                        h = half(t, inanchor)
                                        rightframe = (
                                            inverse(h) if bb < 0 else h if bb else ()
                                        )
                                        if suffix[: len(rightframe)] != rightframe:
                                            raise Exception("bad right frame")
                                        T = suffix[len(rightframe) :]
                                        after = transport_frame(
                                            T, n=n, size=t, start=anchor(j), side="left"
                                        )
                                        check_frame(after)
                                        if (
                                            after.transformed != T
                                            or after.end != inanchor
                                        ):
                                            raise Exception(
                                                "input tape mismatch",
                                                after.end,
                                                inanchor,
                                                after.transformed,
                                                T,
                                            )
                                        counts[rr["id"]] = counts.get(rr["id"], 0) + 1
                                except Exception as e:
                                    failures.append((rr["id"], env, str(e)))

    assert sum(counts.values()) == 7000
    assert failures == []
