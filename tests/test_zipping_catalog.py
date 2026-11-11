"""Exact syntax, orientations and refusal boundaries of the frozen catalog."""

from copy import deepcopy
import hashlib
import json
from itertools import product

import pytest

from zipper_algebra.zipping_algorithm import catalog as mod
from zipper_algebra.zipping_algorithm.catalog import (
    Catalog,
    CatalogError,
    condition,
    expand_polynomial,
    expand_word,
    integer,
    instantiate_rule,
    serialize,
    validate_catalog,
)

ENV = {"n": 6, "i": 3, "j": 1, "k": 2, "l": 1, "eps": -1}
ONE = [{"coefficient": [{"integer": 1, "q_power": 0}], "word": {"type": "identity"}}]


def test_frozen_snapshot_and_independent_copies():
    catalog = Catalog.bundled()
    assert catalog.sha256 == mod.CATALOG_SHA256
    assert catalog.source_sha256 == mod.CATALOG_SHA256
    assert catalog.version == "2026-11-05.1"
    assert catalog.author_review == "partial"
    assert len(catalog.rules) == 128
    assert len({r["base_id"] for r in catalog.rules}) == 32
    assert all(r["execution"]["algebra_verified"] for r in catalog.rules)
    records = catalog.rules
    records[0]["id"] = "changed"
    assert catalog.rules[0]["id"] != "changed"
    record = catalog.get("single.far.positive.I")
    record["lhs"].clear()
    assert catalog.get("single.far.positive.I")["lhs"]
    with pytest.raises(CatalogError, match="Unknown rule"):
        catalog.get("absent")
    manifest = json.loads(mod._data("manifest.json"))
    for name, digest in manifest["files"].items():
        assert hashlib.sha256(mod._data(name)).hexdigest() == digest
    assert b"/home/" not in mod._data("full-catalog.json")


def test_filtered_rules_match_full_catalog_filter_in_the_same_order():
    catalog = Catalog.bundled()
    records = catalog.rules
    for size, friends, index, sign in product((2, 3, 4), (0, 1, 2), range(1, 6), (-1, 1)):
        env = dict(ENV, i=size, k=friends, l=index, eps=sign)
        for suffix in (".I", ".M"):
            excluded = ("zip.flype.",) if size == 2 else ()
            expected = tuple(
                rule for rule in records
                if rule["id"].startswith("zip.") and rule["id"].endswith(suffix)
                and not rule["id"].startswith(excluded)
                and condition(rule["condition"], env)
            )
            assert catalog.select_rules(
                env, prefix="zip.", suffix=suffix, exclude_prefixes=excluded
            ) == expected


def test_filtered_rules_copy_only_matches_and_preserve_mutation_isolation(monkeypatch):
    catalog = Catalog.bundled()
    records = catalog.rules
    expected = tuple(
        rule for rule in records
        if rule["id"].startswith("zip.") and rule["id"].endswith(".I")
        and condition(rule["condition"], ENV)
    )
    copied = []
    original_copy = mod.deepcopy

    def observed_copy(value):
        copied.append(value["id"])
        return original_copy(value)

    monkeypatch.setattr(mod, "deepcopy", observed_copy)
    selected = catalog.select_rules(ENV, prefix="zip.", suffix=".I")
    assert copied == [rule["id"] for rule in expected]
    assert selected == expected and selected
    copied.clear()
    # No condition evaluation or copying is needed when the ID cannot match.
    assert catalog.select_rules({}, prefix="absent.") == ()
    assert copied == []
    with pytest.raises(CatalogError, match="Missing integer parameter"):
        catalog.select_rules({}, prefix="zip.", suffix=".I")
    assert copied == []
    selected[0]["condition"].clear()
    selected[0]["rhs"][0]["coefficient"][0]["integer"] = 999
    assert catalog.select_rules(ENV, prefix="zip.", suffix=".I") == expected
    assert catalog.get(expected[0]["id"]) == expected[0]


@pytest.mark.parametrize("name", ["full-catalog.json", "schema.json"])
def test_tampered_snapshot_rejected(monkeypatch, name):
    original = mod._data
    monkeypatch.setattr(
        mod,
        "_data",
        lambda filename: b"bad" if filename == name else original(filename),
    )
    with pytest.raises(CatalogError, match="hash mismatch"):
        Catalog.bundled()


def test_order_signs_and_half_twists():
    block = {
        "type": "block",
        "length": 2,
        "start": 2,
        "count": 2,
        "order": "descending",
    }
    assert expand_word(block, ENV) == {(0, (3, 4, 2, 3)): 1}
    block["order"] = "ascending"
    assert expand_word(block, ENV) == {(0, (2, 3, 3, 4)): 1}
    fat = {"type": "fat", "length": 3, "start": 1}
    assert expand_word({"type": "bar", "word": fat}, ENV) == {(0, (-1, -2, -3)): 1}
    assert expand_word({"type": "inverse", "word": fat}, ENV) == {(0, (-3, -2, -1)): 1}
    assert expand_word({"type": "reverse", "word": fat}, ENV) == {(0, (3, 2, 1)): 1}
    assert expand_word({"type": "reflect", "word": fat}, ENV) == {(0, (5, 4, 3)): 1}
    assert expand_word({"type": "half_twist", "size": 4, "start": 1}, ENV) == {
        (0, (1, 2, 3, 1, 2, 1)): 1
    }
    assert expand_word(
        {
            "type": "block",
            "length": 0,
            "start": 1,
            "count": 10**30,
            "order": "descending",
        },
        ENV,
    ) == {(0, ()): 1}


def test_slider_exact_definition_and_transform_flags():
    slider = {"type": "slider", "size": 2, "start": 2}
    expected = {(1, (-2,)): 1, (0, ()): 1, (1, ()): -1, (0, (2,)): -1}
    assert expand_word(slider, ENV, expand_sliders=True) == expected
    assert expand_word(
        {"type": "slider", "size": 1, "start": 1}, ENV, expand_sliders=True
    ) == {(0, ()): 1}
    x3 = {"type": "slider", "size": 3, "start": 1}
    normal = expand_word(x3, ENV, expand_sliders=True)
    reflected = {
        (q, tuple((1 if a > 0 else -1) * (3 - abs(a)) for a in w)): c
        for (q, w), c in normal.items()
    }
    assert (
        expand_word(dict(x3, orientation="mirrored"), ENV, expand_sliders=True)
        == reflected
    )
    assert expand_word({"type": "reflect", "word": x3}, ENV) == {
        (0, (("slider", 3, 4, True, False, False),)): 1
    }
    assert expand_word({"type": "reverse", "word": x3}, ENV) == {
        (0, (("slider", 3, 1, False, True, False),)): 1
    }
    assert expand_word({"type": "bar", "word": x3}, ENV) == {
        (0, (("slider", 3, 1, False, False, True),)): 1
    }
    with pytest.raises(CatalogError, match="inverse is only"):
        expand_word({"type": "inverse", "word": slider}, ENV)
    with pytest.raises(CatalogError, match="inverse is only"):
        expand_word({"type": "inverse", "word": slider}, ENV, expand_sliders=True)


def test_context_collection_and_serialization():
    node = {
        "type": "concat",
        "items": [
            {"type": "word_ref", "name": "P"},
            {"type": "slider", "size": 2, "start": 1},
        ],
    }
    assert expand_word(node, ENV, words={"P": [2, -3]}) == {
        (0, (2, -3, ("slider", 2, 1, False, False, False))): 1
    }
    # Literal inverse pairs are retained; this layer performs no braid reduction.
    terms = [
        {
            "word": {"type": "word_ref", "name": "P"},
            "coefficient": [
                {"integer": 2, "q_power": 1},
                {"integer": -1, "q_power": 1},
            ],
        }
    ]
    assert expand_polynomial(terms, ENV, words={"P": [1, -1]}) == {(1, (1, -1)): 1}
    assert serialize({(0, (("slider", 2, 1, False, False, False),)): -2}) == [
        {"integer": -2, "q_power": 0, "word": [["slider", 2, 1, False, False, False]]}
    ]
    cancelled = deepcopy(ONE)
    cancelled[0]["coefficient"].append({"integer": -1, "q_power": 0})
    assert expand_polynomial(cancelled, ENV) == {}


def test_original_symmetry_parameters_and_full_context():
    catalog = Catalog.bundled()
    env = {"n": 6, "i": 1, "j": 2, "l": 4}
    normal = catalog.instantiate("single.far.positive.I", env)
    mirror = catalog.instantiate("single.far.positive.M", env)
    assert normal["lhs"][0]["word"] == [4, ["slider", 1, 2, False, False, False]]
    assert mirror["lhs"][0]["word"] == [2, ["slider", 1, 5, True, False, False]]
    assert mirror["parameters"] == env
    assert mirror["algebra_verified"] is True
    assert (
        mirror["scope"]
        == "Formal instantiation only; control preconditions are not checked."
    )
    nearest = catalog.instantiate(
        "zip.nearest.negative.i_ge_3.I", {"n": 5, "i": 3, "j": 1, "k": 2, "l": 1}
    )
    assert len(nearest["rhs"]) >= 3
    with pytest.raises(CatalogError, match="condition is false"):
        catalog.instantiate("single.far.positive.I", dict(env, l=2))
    source = {"id": "source", "condition": True, "lhs": ONE, "rhs": ONE}
    assert instantiate_rule(source, ENV)["algebra_verified"] is False
    with pytest.raises(CatalogError, match="acceptance"):
        instantiate_rule(source, ENV, require_verified=True)


@pytest.mark.parametrize(
    "expr,expected",
    [
        (3, 3),
        ("i", 3),
        ({"op": "add", "args": [1, 2, 3]}, 6),
        ({"op": "sub", "args": ["n", "i"]}, 3),
        ({"op": "mul", "args": [2, 3, 4]}, 24),
    ],
)
def test_integer_expressions(expr, expected):
    assert integer(expr, ENV) == expected


@pytest.mark.parametrize(
    "expr",
    [
        True,
        '__import__("os")',
        "missing",
        {"op": "eval", "args": [1, 2]},
        {"op": "add", "args": [1]},
        {"op": "sub", "args": [1, 2, 3]},
        {"op": "add", "args": "12"},
    ],
)
def test_bad_integer_expressions(expr):
    with pytest.raises(CatalogError):
        integer(expr, ENV)


def test_missing_parameter():
    with pytest.raises(CatalogError, match="Missing integer"):
        integer("n", {})


@pytest.mark.parametrize(
    "op,expected",
    [
        ("eq", False),
        ("ne", True),
        ("lt", True),
        ("le", True),
        ("gt", False),
        ("ge", False),
    ],
)
def test_comparisons(op, expected):
    assert condition({"op": op, "left": 1, "right": 2}, ENV) is expected


def test_boolean_conditions():
    assert condition(True, ENV)
    assert not condition({"op": "not", "arg": True}, ENV)
    assert condition({"op": "and", "args": [True, True]}, ENV)
    assert condition({"op": "or", "args": [True, False]}, ENV)
    for expr in [None, {"op": "eval"}, {"op": "and", "args": "bad"}]:
        with pytest.raises(CatalogError):
            condition(expr, ENV)


@pytest.mark.parametrize(
    "node",
    [
        None,
        {"type": "unknown"},
        {"type": "identity", "extra": True},
        {"type": "sigma", "index": 6},
        {"type": "sigma", "index": 0},
        {"type": "sigma", "index": 1, "sign": 2},
        {"type": "fat", "length": -1, "start": 1},
        {"type": "block", "length": 1, "start": 1, "count": 1, "order": "bad"},
        {"type": "block", "length": 1, "start": 1, "count": 1},
        {"type": "half_twist", "size": 8, "start": 1},
        {"type": "slider", "size": 8, "start": 1},
        {"type": "slider", "size": 2, "start": 1, "orientation": "bad"},
        {"type": "concat", "items": "bad"},
        {"type": "word_ref", "name": "unbound"},
        {"type": "reflect", "n": 5, "word": {"type": "identity"}},
    ],
)
def test_bad_nodes(node):
    with pytest.raises(CatalogError):
        expand_word(node, ENV)


def test_bad_context_words_and_rank():
    for word in [[0], [True], [6], ["1"]]:
        with pytest.raises(CatalogError):
            expand_word({"type": "word_ref", "name": "P"}, ENV, words={"P": word})
    for env in [{}, {"n": 0}, {"n": True}]:
        with pytest.raises(CatalogError):
            expand_word({"type": "identity"}, env)


@pytest.mark.parametrize(
    "node,budget,expanded",
    [
        ({"type": "identity"}, 0, False),
        ({"type": "fat", "length": 4, "start": 1}, 3, False),
        ({"type": "half_twist", "size": 4, "start": 1}, 3, False),
        ({"type": "slider", "size": 6, "start": 1}, 5, True),
        ({"type": "slider", "size": 3, "start": 1}, 3, True),
        (
            {"type": "concat", "items": [{"type": "fat", "length": 2, "start": 1}] * 2},
            3,
            False,
        ),
    ],
)
def test_expansion_budgets(node, budget, expanded):
    with pytest.raises(CatalogError, match="budget"):
        expand_word(node, ENV, max_terms=budget, expand_sliders=expanded)
    with pytest.raises(CatalogError, match="word budget"):
        expand_word(
            {"type": "word_ref", "name": "P"}, ENV, words={"P": [1, 1]}, max_terms=1
        )


@pytest.mark.parametrize(
    "terms",
    [
        None,
        [None],
        [{"word": {"type": "identity"}, "coefficient": 2}],
        [
            {
                "word": {"type": "identity"},
                "coefficient": [{"integer": True, "q_power": 0}],
            }
        ],
    ],
)
def test_bad_polynomials(terms):
    with pytest.raises(CatalogError):
        expand_polynomial(terms, ENV)


def test_polynomial_budget():
    terms = deepcopy(ONE)
    terms[0]["coefficient"] = [
        {"integer": 1, "q_power": 1},
        {"integer": 1, "q_power": 2},
    ]
    with pytest.raises(CatalogError, match="term budget"):
        expand_polynomial(terms, ENV, max_terms=1)


def test_catalog_schema_validation():
    rule = {"id": "identity", "condition": True, "lhs": ONE, "rhs": ONE}
    data = {"schema_version": "1.0", "rules": [rule]}
    assert validate_catalog(data)
    assert Catalog(data).sha256 is None
    broken = [
        {},
        dict(data, schema_version="2.0"),
        dict(data, rules="bad"),
        dict(data, rules=[rule, rule]),
        dict(data, rules=[dict(rule, id="")]),
        dict(
            data,
            rules=[dict(rule, lhs=[{"word": {"type": "garbage"}, "coefficient": []}])],
        ),
        dict(data, rules=[dict(rule, condition={"op": "add", "args": [1, 2]})]),
        dict(
            data,
            rules=[
                dict(
                    rule,
                    lhs=[{"word": {"type": "sigma", "index": True}, "coefficient": []}],
                )
            ],
        ),
        dict(
            data,
            rules=[
                dict(
                    rule,
                    lhs=[{"word": {"type": "identity", "extra": 1}, "coefficient": []}],
                )
            ],
        ),
        dict(
            data,
            rules=[
                dict(
                    rule,
                    lhs=[
                        {
                            "word": {
                                "type": "fat",
                                "length": {"op": "sub", "args": [1, 2, 3]},
                                "start": 1,
                            },
                            "coefficient": [],
                        }
                    ],
                )
            ],
        ),
    ]
    for item in broken:
        with pytest.raises(CatalogError):
            validate_catalog(item)


def test_malformed_coefficient_entry():
    with pytest.raises(CatalogError, match="fields"):
        expand_polynomial([{"word": {"type": "identity"}, "coefficient": [None]}], ENV)
