"""Bounded derived-rule sharing retains guards, widths and independent records."""
from copy import deepcopy
from itertools import product

import pytest

from zipper_algebra.zipping_algorithm import semilinear
from zipper_algebra.zipping_algorithm.catalog import Catalog, CatalogError, condition, instantiate_rule
from zipper_algebra.zipping_algorithm.framed import rule_for_event

ENV = {"n": 3, "i": 3, "j": 1, "k": 0, "l": 1, "eps": -1}
RID = "zip.flype.negative.M"


def test_repeated_derivation_is_reused_but_returned_records_are_independent(monkeypatch):
    catalog = Catalog.bundled()
    derive = semilinear._semilinear_template
    expected = semilinear.semilinear_zip_rule(catalog.get(RID), ENV)
    calls = []

    def counted(rule):
        calls.append(rule["id"])
        return derive(rule)

    monkeypatch.setattr(semilinear, "_semilinear_template", counted)
    first = catalog.semilinear_rule(RID, ENV)
    first["rhs"].clear()
    first["semilinear"]["parameters"]["n"] = 99
    second = rule_for_event(catalog, {"rule_id": RID[:-2] + ".BM", "parameters": ENV})
    assert second == expected
    assert calls == [RID]
    second["condition"].clear()
    assert catalog.semilinear_rule(RID, ENV) == expected
    assert catalog.rule_cache_info()["hits"] == 2
    assert catalog.rule_cache_info()["entries"] == 1


def test_guards_and_integer_types_are_checked_after_cache_warmup():
    catalog = Catalog.bundled()
    catalog.semilinear_rule(RID, ENV)
    for env in (dict(ENV, i=True), dict(ENV, k=99), {"i": 3}):
        with pytest.raises(CatalogError):
            catalog.semilinear_rule(RID, env)
    with pytest.raises(CatalogError, match="Unknown rule"):
        catalog.semilinear_rule("absent", ENV)
    with pytest.raises(CatalogError, match="accepted frozen mirrored"):
        catalog.semilinear_rule("single.far.positive.I", ENV)
    assert catalog.rule_cache_info()["hits"] == 0


def test_family_reuse_keeps_fresh_typed_bindings_and_copies_inputs():
    catalog = Catalog.bundled()
    env = dict(ENV, unused=True)
    result = catalog.semilinear_rule(RID, env)
    env["unused"] = 1
    other = catalog.semilinear_rule(RID, env)
    assert type(result["semilinear"]["parameters"]["unused"]) is bool
    assert type(other["semilinear"]["parameters"]["unused"]) is int
    assert catalog.rule_cache_info()["misses"] == 1
    assert catalog.rule_cache_info()["hits"] == 1
    assert catalog.semilinear_rule(RID, dict(ENV, unused=True)) == result


def test_cache_budget_eviction_disable_oversize_and_reset():
    catalog = Catalog.bundled()
    catalog.configure_rule_cache(max_bytes=1000000, max_entries=1)
    for source in (RID, "zip.flype.positive.M", RID):
        catalog.semilinear_rule(source, ENV)
    info = catalog.rule_cache_info()
    assert info["misses"] == 3 and info["evictions"] == 2 and info["entries"] == 1
    assert info["retained_bytes"] <= info["max_bytes"]
    assert info["resident_estimate_bytes"] > info["retained_bytes"]
    catalog.configure_rule_cache(max_bytes=0, max_entries=0)
    assert catalog.rule_cache_info()["entries"] == 0
    catalog.semilinear_rule(RID, ENV)
    catalog.configure_rule_cache(max_bytes=1, max_entries=100)
    catalog.semilinear_rule(RID, ENV)
    assert catalog.rule_cache_info()["entries"] == 0
    catalog.rule_cache_clear()
    info = catalog.rule_cache_info()
    assert all(info[key] == 0 for key in ("entries", "retained_bytes", "hits", "misses", "evictions"))
    with pytest.raises(CatalogError, match="nonnegative"):
        catalog.configure_rule_cache(max_bytes=True, max_entries=1)


def test_family_template_preserves_nonstandard_extra_bindings():
    catalog = Catalog.bundled()
    for extra in ((1, 2), 1 << 2000):
        env = dict(ENV, unused=extra)
        expected = semilinear.semilinear_zip_rule(catalog.get(RID), env)
        assert catalog.semilinear_rule(RID, env) == expected
    assert catalog.rule_cache_info()["entries"] == 1
    assert catalog.rule_cache_info()["hits"] == 1


def test_every_bm_family_matches_uncached_derivation_across_recursive_widths():
    catalog = Catalog.bundled()
    seen = set()
    for rule in catalog.rules:
        if not rule["id"].startswith("zip.") or not rule["id"].endswith(".M"):
            continue
        for size in (2, 3, 4, 8):
            n = size + 3
            bindings = ({"n": n, "i": size, "j": j, "k": k, "l": l, "eps": eps}
                        for j, k, l, eps in product(range(1, n-size+2), range(n-size+1), range(1, n), (-1, 1)))
            env = next((env for env in bindings if condition(rule["condition"], env)), None)
            if env is None:
                continue
            original = deepcopy(rule)
            expected = semilinear.semilinear_zip_rule(rule, env)
            assert catalog.semilinear_rule(rule["id"], env) == expected
            assert catalog.semilinear_rule(rule["id"], env) == expected
            assert rule == original
            seen.add(rule["id"])
    assert seen == {rule["id"] for rule in catalog.rules
                    if rule["id"].startswith("zip.") and rule["id"].endswith(".M")}


def test_concurrent_same_catalog_derivation_has_one_retained_record():
    from concurrent.futures import ThreadPoolExecutor

    catalog = Catalog.bundled()
    with ThreadPoolExecutor(max_workers=8) as pool:
        values = list(pool.map(lambda _: catalog.semilinear_rule(RID, ENV), range(32)))
    assert all(value == values[0] for value in values)
    assert len({id(value) for value in values}) == 32
    assert catalog.rule_cache_info()["misses"] == 1
    assert catalog.rule_cache_info()["hits"] == 31


def test_cached_ast_does_not_cache_or_bypass_expansion_limits():
    catalog = Catalog.bundled()
    for _ in range(2):
        rule = catalog.semilinear_rule(RID, ENV)
        with pytest.raises(CatalogError, match="budget"):
            instantiate_rule(rule, ENV, max_terms=1)
        assert instantiate_rule(rule, ENV, max_terms=100000)["rhs"]
    assert catalog.rule_cache_info()["hits"] == 1


def test_programmatic_catalog_aliases_have_json_value_semantics():
    from zipper_algebra.zipping_algorithm.catalog import _data
    import json

    raw = json.loads(_data("full-catalog.json"))
    rule = next(rule for rule in raw["rules"] if rule["id"] == RID)
    shared_coefficient = {"integer": 1, "q_power": 7}
    rule["lhs"][0]["coefficient"] = [shared_coefficient]
    rule["rhs"][0]["coefficient"] = [shared_coefficient]
    # A repeated term exercises the alias in both its word and coefficient.
    rule["rhs"].append(rule["rhs"][0])
    plain = json.loads(json.dumps(raw))
    plain_rule = next(rule for rule in plain["rules"] if rule["id"] == RID)
    expected = semilinear.semilinear_zip_rule(plain_rule, ENV)
    assert expected["lhs"][0]["coefficient"][0]["q_power"] == -7
    assert semilinear.semilinear_zip_rule(rule, ENV) == expected
    for enabled in (False, True):
        catalog = Catalog(raw)
        catalog.configure_rule_cache(max_bytes=1000000 if enabled else 0, max_entries=10)
        assert catalog.semilinear_rule(RID, ENV) == expected
        assert catalog.semilinear_rule(RID, ENV) == expected
        assert catalog.get(RID)["lhs"][0]["coefficient"][0]["q_power"] == 7


def test_binding_snapshot_prevents_dict_subclass_validation_bind_disagreement():
    class AdversarialBindings(dict):
        def __getitem__(self, key):
            return 99 if key == "i" else super().__getitem__(key)

    hostile = AdversarialBindings(ENV)
    catalog = Catalog.bundled()
    expected = semilinear.semilinear_zip_rule(catalog.get(RID), ENV)
    assert semilinear.semilinear_zip_rule(catalog.get(RID), hostile) == expected
    for enabled in (False, True):
        catalog = Catalog.bundled()
        catalog.configure_rule_cache(max_bytes=1000000 if enabled else 0, max_entries=10)
        for _ in range(2):
            actual = catalog.semilinear_rule(RID, hostile)
            assert actual == expected
            assert actual["semilinear"]["parameters"]["i"] == 3
            assert actual["semilinear"]["common_slider_unit"]["q_power"] == -3
