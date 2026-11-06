"""Host admission rejects invalid occurrences before GPU buffer allocation."""

from copy import deepcopy

import pytest

from zipper_algebra.zipping_algorithm import gpu_rewrite_host as host
from zipper_algebra.zipping_algorithm import gpu_rewrite_program as gp
from zipper_algebra.zipping_algorithm.catalog import CatalogError
from zipper_algebra.zipping_algorithm.interpreter import Branch, Limits, LimitExceeded, ReplayError


def term(word, coefficient=1):
    return {"word": word, "coefficient": [{"integer": coefficient, "q_power": 0}]}


def case(rhs=None):
    return {
        "root": 0, "id": 1, "batch": 2, "n": 8,
        "rule": {"id": "test", "execution": {"algebra_verified": True},
                 "condition": True,
                 "lhs": [term({"type": "sigma", "index": 1})],
                 "rhs": [term({"type": "sigma", "index": 2})] if rhs is None else rhs},
        "parameters": {"n": 8}, "offset": 1,
        "predecessor": Branch((-1, 1, -2), ((-1, 3), (2, -7)), '{"saved":true}').to_json(),
    }


def admit(value):
    return host.admit(value, gp.compile_rule(value["rule"]))


def test_common_contract_preserves_context_limits_and_control():
    value = case()
    before = deepcopy(value)
    accepted = admit(value)
    assert accepted.lengths == (3,)
    assert accepted.end == 2
    assert accepted.limits.max_steps == accepted.limits.max_terms == 100000
    assert accepted.limits.expansion_budget == 1000000
    assert accepted.branch.control == {"saved": True}
    assert host.expand_cpu([value])[0]["branches"] == [
        Branch((-1, 2, -2), ((-1, 3), (2, -7)), '{"saved":true}').to_json()]
    assert value == before
    accepted.raw["parameters"]["n"] = 9
    assert value == before


@pytest.mark.parametrize("change", [
    lambda c: c.pop("root"),
    lambda c: c.update(extra=1),
    lambda c: c.update(root=True),
    lambda c: c.update(id=-1),
    lambda c: c.update(batch=1.5),
    lambda c: c.update(rule=None),
    lambda c: c["rule"].update(id=3),
    lambda c: c.update(parameters=None),
    lambda c: c.update(limits=None),
])
def test_malformed_common_contracts_fail_before_compilation(change):
    value = case()
    change(value)
    with pytest.raises(ReplayError):
        host.inputs(value)
    with pytest.raises(ReplayError):
        host.expand_cpu([value])
    with pytest.raises(ReplayError):
        host.inputs(None)


@pytest.mark.parametrize("change,message", [
    (lambda c: c["parameters"].update(n=9), "ambient"),
    (lambda c: c["rule"]["execution"].update(algebra_verified=False), "verification"),
    (lambda c: c["rule"].update(condition=False), "guard"),
    (lambda c: c["rule"].update(lhs=[]), "monomial"),
    (lambda c: c["rule"]["lhs"].append(term({"type": "sigma", "index": 2})), "monomial"),
    (lambda c: c["rule"]["lhs"][0]["coefficient"][0].update(q_power=1), "nonempty"),
    (lambda c: c["rule"]["lhs"][0]["coefficient"][0].update(integer=-1), "nonempty"),
    (lambda c: c["rule"]["lhs"][0].update(word={"type": "identity"}), "nonempty"),
    (lambda c: c.update(offset=0), "literal"),
])
def test_guard_lhs_and_literal_occurrence_checks(change, message):
    value = case()
    change(value)
    with pytest.raises(ReplayError, match=message):
        admit(value)


def test_all_primitive_lengths_empty_blocks_and_unbounded_slider_sizes():
    value = case([
        term({"type": "block", "length": 2, "start": 2, "count": 3,
              "order": "descending", "sign": -1}),
        term({"type": "half_twist", "size": 4, "start": 1}),
        term({"type": "slider", "size": 8, "start": 1}),
        term({"type": "word_ref", "name": "U"}),
        term({"type": "fat", "length": 0, "start": 100}),
    ])
    value["words"] = {"U": [1, -2, 4]}
    assert admit(value).lengths == (8, 8, 3, 5, 2)
    value = case([term({"type": "slider", "size": 1000, "start": 1})])
    value["n"] = value["parameters"]["n"] = 1000
    assert admit(value).lengths == (3,)
    value["predecessor"] = Branch((("slider", 2, 1, False, False, False), 1)).to_json()
    value["offset"] = 1
    assert admit(value).lengths == (2,)


def test_word_and_coefficient_budgets_fail_closed():
    value = case([term({"type": "fat", "length": 6, "start": 1})])
    value["limits"] = vars(Limits(expansion_budget=5))
    with pytest.raises(LimitExceeded, match="expansion word"):
        admit(value)
    value = case([term({"type": "concat", "items": [
        {"type": "fat", "length": 4, "start": 1},
        {"type": "fat", "length": 4, "start": 1}]})])
    value["limits"] = vars(Limits(expansion_budget=5))
    with pytest.raises(LimitExceeded, match="concatenation"):
        admit(value)
    value = case([term({"type": "fat", "length": 4, "start": 1})])
    value["limits"] = vars(Limits(max_word_size=4))
    with pytest.raises(LimitExceeded, match="successor"):
        admit(value)
    value = case()
    value["limits"] = vars(Limits(expansion_budget=1))
    with pytest.raises(LimitExceeded, match="coefficient expansion"):
        admit(value)


def test_unknown_manual_program_and_oracle_identity_fail_closed():
    value = case()
    program = gp.compile_rule(value["rule"])
    term_ = gp.CompiledTerm(((99, 0, 0, 0, 0, 0),), ((0, 1),))
    bad = gp.CompiledRule(program.expressions, (term_,), ())
    with pytest.raises(CatalogError, match="Unknown"):
        host.admit(value, bad)
    with pytest.raises(ReplayError, match="identity"):
        host._RuleView(value["rule"]).get("different")


def test_admission_never_expands_rhs_and_cpu_oracle_never_compiles(monkeypatch):
    value = case()
    actual_expand = host.expand_polynomial
    seen = []

    def lhs_only(polynomial, *args, **kwargs):
        assert polynomial == value["rule"]["lhs"]
        seen.append(True)
        return actual_expand(polynomial, *args, **kwargs)

    monkeypatch.setattr(host, "expand_polynomial", lhs_only)
    admit(value)
    assert seen == [True]
    monkeypatch.setattr(gp, "compile_rule", lambda *_: pytest.fail("oracle compiled GPU program"))
    assert host.expand_cpu([value])[0]["branches"]
