"""Packed compiler semantics; emulation here is not GPU hardware evidence."""

from collections import defaultdict
from copy import deepcopy
from itertools import product

import pytest

from zipper_algebra.zipping_algorithm import gpu_rewrite_program as gp
from zipper_algebra.zipping_algorithm.catalog import (
    Catalog, CatalogError, condition, expand_polynomial,
)


def rule(word, coefficient=None):
    return {"rhs": [{"word": word, "coefficient": coefficient or [
        {"integer": 1, "q_power": 0}]}]}


def sigma(index=1):
    return {"type": "sigma", "index": index}


def emit(program, env, words=None):
    """Independent readable descriptor interpreter, used only by these tests."""
    values = gp.evaluate(program, env)
    result = defaultdict(int)
    for term in program.terms:
        output = []
        for kind, a, b, c, d, flags in term.primitives:
            if kind == gp.WORD_REF:
                word = list((words or {})[program.refs[a]])
            elif kind == gp.SIGMA:
                word = [values[a] * values[b]]
            elif kind == gp.BLOCK:
                offsets = range(values[c])
                if flags & gp.BLOCK_DESCENDING:
                    offsets = reversed(offsets)
                word = [values[d] * (values[b] + t + x)
                        for t in offsets for x in range(values[a])]
            elif kind == gp.HALF_TWIST:
                word = [x for length in range(values[a] - 1, 0, -1)
                        for x in range(values[b], values[b] + length)]
            else:
                assert kind == gp.SLIDER
                word = [("slider", values[a], values[b],
                         bool(flags & gp.SLIDER_MIRRORED),
                         bool(flags & gp.SLIDER_REVERSED),
                         bool(flags & gp.SLIDER_BARRED))]
            if flags & gp.REVERSE:
                word.reverse()
            transformed = []
            for token in word:
                if type(token) is int:
                    if flags & gp.NEGATE:
                        token = -token
                    if flags & gp.REFLECT:
                        token = (1 if token > 0 else -1) * (env["n"] - abs(token))
                elif flags & gp.REFLECT:
                    _, size, start, mirrored, rev, bar = token
                    token = ("slider", size, env["n"] - size - start + 2,
                             not mirrored, rev, bar)
                transformed.append(token)
            output.extend(transformed)
        for power, coefficient in term.coefficients:
            result[(values[power], tuple(output))] += coefficient
    return {key: value for key, value in result.items() if value}


def test_every_frozen_rule_and_derived_bm_matches_catalog():
    catalog = Catalog.bundled()
    bindings = {}
    for r in catalog.rules:
        base = r["base_id"]
        if base not in bindings:
            for i, j, k, l, eps in product(range(1, 6), range(1, 8),
                                         range(6), range(1, 8), (-1, 1)):
                env = dict(n=8, i=i, j=j, k=k, l=l, eps=eps)
                if condition(r["condition"], env):
                    bindings[base] = env
                    break
            else:
                pytest.fail("No binding found for " + base)
        env = bindings[base]
        compiled = gp.compile_rule(r)
        assert emit(compiled, env) == expand_polynomial(r["rhs"], env)
        if r["id"].startswith("zip.") and r["id"].endswith(".M"):
            bm = catalog.semilinear_rule(r["id"], env)
            assert emit(gp.compile_rule(bm), env) == expand_polynomial(bm["rhs"], env)
    assert len(bindings) == 32


def test_nested_wrappers_words_zero_factors_and_expression_operations():
    block = {"type": "block", "length": 2, "start": 2, "count": 2,
             "sign": -1, "order": "descending"}
    source = {"type": "concat", "items": [
        {"type": "identity"}, sigma(), block,
        {"type": "fat", "length": 2, "start": 1},
        {"type": "half_twist", "size": 3, "start": 2},
        {"type": "word_ref", "name": "U"},
        {"type": "word_ref", "name": "U"},
    ]}
    source = {"type": "reverse", "word": {
        "type": "reflect", "n": 8, "word": {
            "type": "bar", "word": {"type": "inverse", "word": source}}}}
    polynomial = rule(source, [
        {"integer": -1, "q_power": {"op": "mul", "args": [
            "i", {"op": "sub", "args": ["n", 1]}, 2]}},
        {"integer": 0, "q_power": {"op": "add", "args": [1, 2, 3]}},
        {"integer": 1, "q_power": "eps"},
    ])
    env = dict(n=8, i=3, eps=-1)
    words = {"U": (1, -3, 4)}
    program = gp.compile_rule(polynomial)
    assert emit(program, env, words) == expand_polynomial(polynomial["rhs"], env, words=words)
    assert program.refs == ("U",)
    assert all(all(type(v) is int for v in row) for row in program.expressions)
    assert len(program.expressions) == len(set(program.expressions))


def test_slider_flags_size_one_and_parameterized_large_sizes():
    word = {"type": "slider", "size": "i", "start": "j", "orientation": "mirrored"}
    for name in ["reverse", "bar", "reflect", "reflect", "bar", "reverse"]:
        word = {"type": name, "word": word}
    program = gp.compile_rule(rule(word))
    for i in [1, 3, 1000]:
        env = dict(n=i + 5, i=i, j=2)
        assert emit(program, env) == expand_polynomial(rule(word)["rhs"], env)
    assert gp.compile_rule(rule(word)) == program
    # All flags are symbolic: no finite slider-size table is compiled.
    a = program.terms[0].primitives[0][1]
    assert program.expressions[a][0] == gp.PARAM


def test_reflection_width_and_empty_polynomials():
    word = {"type": "reflect", "n": 3, "word": sigma()}
    program = gp.compile_rule(rule(word))
    assert emit(program, {"n": 3}) == {(0, (2,)): 1}
    with pytest.raises(CatalogError, match="ambient n"):
        gp.evaluate(program, {"n": 4})
    assert gp.compile_rule({"rhs": []}).terms == ()
    assert emit(gp.compile_rule(rule({"type": "identity"})), {"n": 1}) == {(0, ()): 1}


@pytest.mark.parametrize("word", [
    None, {}, {"type": "absent"}, {"type": "concat", "items": 1},
    {"type": "inverse", "word": None},
    {"type": "inverse", "word": {"type": "concat", "items": 1}},
    {"type": "identity", "extra": 1},
    {"type": "sigma", "index": 1, "extra": 0},
    {"type": "block", "length": 1, "start": 1, "count": 1, "order": "bad"},
    {"type": "slider", "size": 2, "start": 1, "orientation": "bad"},
    {"type": "word_ref", "name": 1},
    {"type": "inverse", "word": {"type": "slider", "size": 2, "start": 1}},
    {"type": "inverse", "word": {"type": "concat", "items": [
        sigma(), {"type": "bar", "word": {"type": "slider", "size": 2, "start": 1}}]}},
])
def test_malformed_and_unsupported_words_fail_closed(word):
    with pytest.raises(CatalogError):
        gp.compile_rule(rule(word))


@pytest.mark.parametrize("expr", [
    True, 1.0, "absent", [], {}, {"op": "div", "args": [1, 2]},
    {"op": "add", "args": []}, {"op": "sub", "args": [1, 2, 3]},
    {"op": "mul", "args": "bad"}, 1 << 63, -(1 << 63) - 1,
])
def test_malformed_and_overflowing_integer_syntax(expr):
    with pytest.raises(CatalogError):
        gp.compile_rule(rule(sigma(expr)))


@pytest.mark.parametrize("polynomial", [
    {}, {"rhs": None}, {"rhs": [None]}, {"rhs": [{"word": sigma(), "coefficient": None}]},
    rule(sigma(), [{"integer": True, "q_power": 0}]),
    rule(sigma(), [{"integer": 2, "q_power": 0}]),
    rule(sigma(), [{"integer": -1, "q_power": 0, "extra": 1}]),
])
def test_malformed_polynomials_and_nonbundled_factors(polynomial):
    with pytest.raises(CatalogError):
        gp.compile_rule(polynomial)


def test_evaluation_rejects_missing_parameters_and_each_intermediate_overflow():
    program = gp.compile_rule(rule(sigma("i")))
    with pytest.raises(CatalogError, match="Missing"):
        gp.evaluate(program, {})
    for value in [True, 1.0, 1 << 63]:
        with pytest.raises(CatalogError):
            gp.evaluate(program, {"i": value})
    # The final mathematical value could fit, but overflow in an earlier
    # device register is rejected rather than silently wrapped.
    expr = {"op": "sub", "args": [
        {"op": "add", "args": ["i", 1]}, 1]}
    program = gp.compile_rule(rule(sigma(expr)))
    with pytest.raises(CatalogError):
        gp.evaluate(program, {"i": gp.I64_MAX})
    program = gp.CompiledRule(((99, 0, 0, 0),), (), ())
    with pytest.raises(CatalogError, match="Unknown"):
        gp.evaluate(program, {})


def test_compilation_uses_value_semantics_and_never_expands_a_successor(monkeypatch):
    from zipper_algebra.zipping_algorithm import catalog
    r = rule({"type": "concat", "items": [
        {"type": "block", "length": "i", "start": "j", "count": "k", "order": "descending"},
        {"type": "half_twist", "size": "i", "start": "j"}]})
    before = deepcopy(r)
    monkeypatch.setattr(catalog, "expand_word", lambda *a, **k: pytest.fail("CPU expansion"))
    program = gp.compile_rule(r)
    assert r == before
    assert len(program.terms[0].primitives) == 2
    assert gp.evaluate(program, dict(n=10000, i=1000, j=1, k=1000))
