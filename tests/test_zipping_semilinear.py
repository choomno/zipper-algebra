"""Exact generic checks of the semilinear mirrored-rule derivation."""
from collections import defaultdict
from copy import deepcopy
from itertools import product

import pytest

from zipper_algebra.zipping_algorithm.catalog import (
    Catalog, CatalogError, condition, expand_polynomial, expand_word,
)
from zipper_algebra.zipping_algorithm.semilinear import (
    PROOF_VERSION, semilinear_word, semilinear_zip_rule,
)
from zipper_algebra.zipping_algorithm.verification import (
    VerificationLimits, _Budget, _relator, _transform, _zero, certify_rule,
)


def _bindings(rule, i):
    n = i + 3
    for j, k, l in product(range(1, n - i + 2), range(n - i + 1), range(1, n)):
        env = {"n": n, "i": i, "j": j, "k": k, "l": l, "eps": 1}
        if condition(rule["condition"], env):
            return env
    return None


def _bar_poly(poly):
    return {(-q, tuple(-x for x in word)): c for (q, word), c in poly.items()}


def test_slider_and_original_relator_scale_in_literal_laurent_words():
    for i in range(2, 8):
        env = {"n": i + 2, "i": i, "j": 2, "k": 0, "l": 1, "eps": 1}
        unit_sign, unit_power = (-1) ** (i - 1), -i * (i - 1) // 2
        for orientation in ("normal", "mirrored"):
            node = {"type": "slider", "size": i, "start": 2, "orientation": orientation}
            original = expand_word(node, env, expand_sliders=True)
            assert _bar_poly(original) == {(q + unit_power, word): unit_sign * c
                                          for (q, word), c in original.items()}
        original_z = _relator(i + 2, i, 2, _Budget(VerificationLimits()))
        assert _bar_poly(original_z) == {(q + unit_power - (i - 1), word): -unit_sign * c
                                        for (q, word), c in original_z.items()}


def test_every_derived_bm_family_has_an_exact_original_relator_cofactor():
    catalog = Catalog.bundled()
    seen = set()
    for source in catalog.rules:
        if not source["id"].startswith("zip.") or not source["id"].endswith(".M"):
            continue
        for i in (2, 3, 4):
            env = _bindings(source, i)
            if env is None:
                continue
            seen.add(source["id"])
            saved_source = deepcopy(source)
            derived = semilinear_zip_rule(source, env)
            assert source == saved_source
            assert derived["id"] == source["id"][:-2] + ".BM"
            assert derived["source_record_id"] == source["id"]
            assert derived["semilinear"]["proof_version"] == PROOF_VERSION
            assert derived["semilinear"]["common_slider_unit"] == {
                "integer": (-1) ** (i - 1), "q_power": -i * (i - 1) // 2}
            assert derived["symmetry"]["invert_q"] is True
            assert derived["symmetry"]["invert_crossing_sign"] is True
            certificate = certify_rule(source["id"], env, catalog=catalog)
            residual = defaultdict(int)
            for side, sign in (("lhs", 1), ("rhs", -1)):
                for key, value in expand_polynomial(derived[side], env, expand_sliders=True).items():
                    residual[key] += sign * value
            for cofactor in certificate["cofactors"]:
                tag = cofactor["relator"]
                original = _relator(env["n"], tag["level"], tag["start"], _Budget(VerificationLimits()))
                left = tuple(-x for x in cofactor["left"])
                right = tuple(-x for x in cofactor["right"])
                q = -cofactor["q_power"] - (i - 1)
                c = -cofactor["integer"]
                for (power, word), value in original.items():
                    key = power + q, left + _transform(word, env["n"], tag["symmetry"]) + right
                    residual[key] -= c * value
            assert _zero(dict(residual), env["n"], _Budget(VerificationLimits()))
    expected = {r["id"] for r in catalog.rules if r["id"].startswith("zip.") and r["id"].endswith(".M")}
    assert seen == expected


def test_transformed_friends_are_negative_and_whole_frames_are_retained():
    catalog = Catalog.bundled()
    source = catalog.get("zip.outer.positive.M")
    env = {"n": 6, "i": 3, "j": 2, "k": 1, "l": 1}
    rule = semilinear_zip_rule(source, env)
    rhs = expand_polynomial(rule["rhs"], env)
    word = next(iter(rhs))[1]
    assert all(x < 0 for x in word if type(x) is int)
    flype = semilinear_zip_rule(catalog.get("zip.flype.positive.M"), dict(env, k=0))
    assert any(term["word"] for term in flype["rhs"])
    assert all(value["integer"] for term in flype["rhs"] for value in term["coefficient"])


def test_recursive_word_nodes_preserve_their_explicit_semantics():
    env = {"n": 5, "i": 3, "j": 1, "eps": 1}
    pure = {"type": "concat", "items": [
        {"type": "identity"}, {"type": "sigma", "index": 1, "sign": "eps"},
        {"type": "fat", "length": 2, "start": 2},
        {"type": "block", "length": 1, "start": 1, "count": 2, "order": "ascending"},
        {"type": "half_twist", "size": 3, "start": 1},
        {"type": "inverse", "word": {"type": "half_twist", "size": 2, "start": 1}},
    ]}
    for wrapper in (None, "bar", "reverse", "reflect", "inverse"):
        node = pure if wrapper is None else {"type": wrapper, "word": pure}
        original = expand_word(node, env)
        actual = expand_word(semilinear_word(node), env)
        residual = defaultdict(int)
        for (q, w), c in actual.items():
            residual[q, w] += c
        for (q, w), c in _bar_poly(original).items():
            residual[q, w] -= c
        assert _zero(dict(residual), env["n"], _Budget(VerificationLimits()))
    named = {"type": "word_ref", "name": "P"}
    assert expand_word(semilinear_word(named), env, words={"P": [1, -2]}) == {(0, (-1, 2)): 1}
    with pytest.raises(CatalogError):
        semilinear_word({"type": "inverse", "word": {"type": "slider", "size": 3, "start": 1}})
    with pytest.raises(CatalogError):
        semilinear_word({"type": "unknown"})


def test_invalid_source_guards_and_slider_degrees_are_refused():
    catalog = Catalog.bundled()
    source = catalog.get("zip.outer.positive.M")
    env = {"n": 6, "i": 3, "j": 2, "k": 1, "l": 1}
    wrong = deepcopy(source)
    wrong["execution"]["algebra_verified"] = False
    candidates = [wrong, dict(source, id="layer.positive.M"), dict(source, id="zip.outer.positive.I"),
                  dict(source, id=3), dict(source, lhs=[])]
    missing = deepcopy(source)
    missing["rhs"][0]["word"] = {"type": "identity"}
    candidates.append(missing)
    other = deepcopy(source)
    other["rhs"][0]["word"] = {"type": "slider", "size": 2, "start": 1}
    candidates.append(other)
    inverse = deepcopy(source)
    inverse["rhs"][0]["word"] = {"type": "inverse", "word": {"type": "slider", "size": 3, "start": 1}}
    candidates.append(inverse)
    for rule in candidates:
        with pytest.raises(CatalogError):
            semilinear_zip_rule(rule, env)
    for bad_env in (dict(env, i=True), dict(env, i=1), dict(env, k=-1)):
        with pytest.raises(CatalogError):
            semilinear_zip_rule(source, bad_env)
