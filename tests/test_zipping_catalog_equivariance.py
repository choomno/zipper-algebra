"""Uniform symbolic catalog identities; no assertion of raw guard-list invariance."""

from collections import Counter, defaultdict

import pytest

from zipper_algebra.zipping_algorithm.catalog import Catalog
from zipper_algebra.zipping_algorithm.semilinear import semilinear_word

V = ("n", "i", "j", "k", "l", "eps")
Z = (0,) * 6


def clean(p):
    return {k: v for k, v in p.items() if v}


def add(a, b):
    c = Counter(a)
    c.update(b)
    return clean(c)


def mul(a, b):
    c = defaultdict(int)
    for x, a1 in a.items():
        for y, b1 in b.items():
            c[tuple(i + j for i, j in zip(x, y))] += a1 * b1
    return clean(c)


def poly(x):
    if type(x) is int:
        return {Z: x} if x else {}
    if isinstance(x, str):
        if x not in V:
            raise ValueError("unknown variable")
        return {tuple(int(v == x) for v in V): 1}
    if not isinstance(x, dict) or set(x) != {"op", "args"}:
        raise ValueError("invalid expression")
    if (
        x["op"] not in ("add", "sub", "mul")
        or not isinstance(x["args"], list)
        or len(x["args"]) < 2
        or (x["op"] == "sub" and len(x["args"]) != 2)
    ):
        raise ValueError("invalid operation")
    args = [poly(a) for a in x["args"]]
    if x["op"] == "sub":
        return add(args[0], {k: -v for k, v in args[1].items()})
    r = {Z: 1} if x["op"] == "mul" else {}
    for a in args:
        r = mul(r, a) if x["op"] == "mul" else add(r, a)
    return r


def deriv(p, var):
    c = {}
    v = V.index(var)
    for m, a in p.items():
        if m[v]:
            mm = list(m)
            mm[v] -= 1
            c[tuple(mm)] = a * m[v]
    return c


def signature(x):
    p = poly(x)
    return deriv(p, "n"), add(deriv(p, "j"), deriv(p, "l"))


def invariant(x):
    return signature(x) == ({}, {})


def dependencies(x):
    return sorted({v for m in poly(x) for v, k in zip(V, m) if k})


def sub(a, b):
    return {"op": "sub", "args": [a, b]}


def sum_expr(*a):
    return {"op": "add", "args": list(a)}


def inv(w):
    return {"type": "inverse", "word": w}


def mirror(w):
    t = w["type"]
    r = dict(w)
    if t == "identity":
        return r
    if t == "sigma":
        r["index"] = sub("n", w["index"])
        return r
    if t in ("fat", "block"):
        r["start"] = (
            sum_expr("n", sub(0, w["start"]), sub(0, w["length"]), 1)
            if t == "fat"
            else sum_expr(
                "n", sub(0, w["start"]), sub(0, w["length"]), sub(0, w["count"]), 2
            )
        )
        r["sign"] = sub(0, w.get("sign", 1))
        return inv(r)
    if t in ("half_twist", "slider"):
        r["start"] = sum_expr("n", sub(0, w["start"]), sub(0, w["size"]), 2)
        if t == "slider":
            r["orientation"] = (
                "normal" if w.get("orientation", "normal") == "mirrored" else "mirrored"
            )
        return r
    if t == "concat":
        r["items"] = [mirror(v) for v in w["items"]]
        return r
    if t in ("inverse", "reverse", "bar"):
        r["word"] = mirror(w["word"])
        return r
    raise ValueError(t)


def norm(w):
    r = dict(w)
    for k in ("start", "index", "length", "size", "count", "sign"):
        if k in r:
            r[k] = sorted(poly(r[k]).items())
    if "word" in r:
        r["word"] = norm(r["word"])
    if "items" in r:
        r["items"] = [norm(v) for v in r["items"]]
    return r


@pytest.mark.parametrize(
    "expression",
    [
        True,
        False,
        1.0,
        None,
        "unknown",
        {"op": "divide", "args": [1, 2]},
        {"op": "add", "args": [1]},
        {"op": "sub", "args": [1, 2, 3]},
        {"op": "add", "args": (1, 2)},
        {"op": "add", "args": [1, 2], "extra": 3},
    ],
)
def test_polynomial_parser_rejects_invalid_syntax(expression):
    with pytest.raises(ValueError):
        poly(expression)


def nodes(word):
    assert word["type"] in {
        "identity",
        "sigma",
        "slider",
        "fat",
        "block",
        "half_twist",
        "concat",
        "inverse",
        "bar",
        "reverse",
        "reflect",
        "word_ref",
    }
    yield word
    if "word" in word:
        yield from nodes(word["word"])
    for child in word.get("items", []):
        yield from nodes(child)


def test_uniform_catalog_translation_and_coefficient_identities():
    counts = Counter()
    for rule in Catalog.bundled().rules:
        suffix = rule["id"].rsplit(".", 1)[1]
        assert suffix in ("I", "M", "phi", "Mphi")
        expected = ({Z: 1}, {Z: -1}) if suffix in ("M", "phi") else ({}, {Z: 1})
        for side in ("lhs", "rhs"):
            for term in rule[side]:
                for node in nodes(term["word"]):
                    for field in ("start", "index"):
                        if field in node:
                            assert signature(node[field]) == expected, (
                                rule["id"],
                                node,
                            )
                            counts["coordinates"] += 1
                    for field in ("length", "size", "count", "sign"):
                        if field in node:
                            assert invariant(node[field]), (rule["id"], node)
                            counts["dimensions"] += 1
                for coefficient in term["coefficient"]:
                    assert type(coefficient["integer"]) is int
                    for field in ("integer", "q_power"):
                        assert set(dependencies(coefficient[field])) <= {"i"}, (
                            rule["id"],
                            coefficient,
                        )
                        counts["coefficients"] += 1
        counts["rules"] += 1
    assert counts == {
        "rules": 128,
        "coordinates": 952,
        "dimensions": 1532,
        "coefficients": 816,
    }


def test_all_zip_mirrors_equal_symbolic_canonical_frame_transform():
    catalog = Catalog.bundled()
    rules = terms = 0
    for rule in catalog.rules:
        if not (rule["id"].startswith("zip.") and rule["id"].endswith(".I")):
            continue
        reflected = catalog.get(rule["id"][:-2] + ".M")
        assert rule["condition"] == reflected["condition"]
        for side in ("lhs", "rhs"):
            for source, target in zip(rule[side], reflected[side], strict=True):
                assert source["coefficient"] == target["coefficient"]
                assert norm(mirror(source["word"])) == norm(target["word"])
                terms += 1
        rules += 1
    assert (rules, terms) == (20, 60)


def test_bm_preserves_slider_flags_and_packet_dimensions():
    for rule in Catalog.bundled().rules:
        if not (rule["id"].startswith("zip.") and rule["id"].endswith(".M")):
            continue
        for side in ("lhs", "rhs"):
            for term in rule[side]:
                before = list(nodes(term["word"]))
                after = list(nodes(semilinear_word(term["word"])))
                assert [x for x in before if x["type"] == "slider"] == [
                    x for x in after if x["type"] == "slider"
                ]
                for kind in ("sigma", "fat", "block"):
                    left = [x for x in before if x["type"] == kind]
                    right = [x for x in after if x["type"] == kind]
                    for source, target in zip(left, right, strict=True):
                        assert poly(target.get("sign", 1)) == {
                            k: -v for k, v in poly(source.get("sign", 1)).items()
                        }
                        assert {k: v for k, v in source.items() if k != "sign"} == {
                            k: v for k, v in target.items() if k != "sign"
                        }
