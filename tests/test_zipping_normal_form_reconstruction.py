"""Compare current-output reconstruction against independent full executions."""
from copy import deepcopy

import pytest

from zipper_algebra.zipping_algorithm.collected import CollectedSession
from zipper_algebra.zipping_algorithm.input_normal_form import input_normal_form, leaf_requests
from zipper_algebra.zipping_algorithm.normal_form_reconstruction import (
    execution_cache_key, execution_identity, reconstruct_output,
)
from zipper_algebra.zipping_algorithm.output_normal_form import from_certificate


@pytest.fixture(scope="module")
def session():
    return CollectedSession()


def prepared(session, word, n):
    form = input_normal_form(word, n=n)
    requests = leaf_requests(form)
    results = {r["key"]: session.reduce(r["word"], n=r["n"]) for r in requests}
    target = session.reduce(word, n=n)
    return form, requests, results, target, execution_identity(target)


@pytest.mark.parametrize("word,n", [([], 1), ([], 4), ([-1], 2), ([1,-1,1], 2),
    ([1,3], 4), ([3,1,-3], 4), ([2,2,4], 6), ([1,2,1,5], 6),
    ([1,-3,1,3,1,-3,2,1], 4), ([1,3,5,1], 6)])
def test_reconstructs_complete_full_fields(session, word, n):
    form, _, results, target, identity = prepared(session, word, n)
    before = deepcopy((form, results))
    assert reconstruct_output(form, results, source_identity=identity) == from_certificate(target)
    assert (form, results) == before


def test_duplicate_requests_reuse_one_run(session):
    form, requests, results, target, identity = prepared(session, [1,3], 4)
    assert len(requests) == 1
    assert reconstruct_output(form, results, source_identity=identity) == from_certificate(target)
    key = execution_cache_key(requests[0], identity)
    other = deepcopy(identity)
    other["catalog_sha256"] = "0" * 64
    assert execution_cache_key(requests[0], other) != key
    request = dict(requests[0], word=[-1])
    with pytest.raises(ValueError, match="key disagrees"):
        execution_cache_key(request, identity)


def test_bindings_and_complete_requirement(session):
    form, requests, results, _, identity = prepared(session, [1,3], 4)
    key = requests[0]["key"]
    with pytest.raises(ValueError, match="keys must exactly"):
        reconstruct_output(form, {}, source_identity=identity)
    with pytest.raises(ValueError, match="keys must exactly"):
        reconstruct_output(form, dict(results, extra=results[key]), source_identity=identity)
    wrong = deepcopy(results)
    wrong[key]["input_word"] = [-1]
    with pytest.raises(ValueError, match="literal input"):
        reconstruct_output(form, wrong, source_identity=identity)
    wrong = deepcopy(results)
    wrong[key]["input_word"] = [True]
    with pytest.raises(ValueError, match="literal input"):
        reconstruct_output(form, wrong, source_identity=identity)
    wrong = deepcopy(results)
    wrong[key]["status"] = "incomplete"
    with pytest.raises(ValueError, match="completed"):
        reconstruct_output(form, wrong, source_identity=identity)
    wrong = deepcopy(results)
    wrong[key]["catalog_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="producer identity"):
        reconstruct_output(form, wrong, source_identity=identity)


@pytest.mark.parametrize("where", ["term", "component", "record", "parity"])
def test_unknown_embedding_semantics_rejected(session, where):
    form, requests, results, _, identity = prepared(session, [1,3], 4)
    results = deepcopy(results)
    term = results[requests[0]["key"]]["output"][0]
    target = {"term": term, "component": term["components"][0],
              "record": term["components"][0]["record"],
              "parity": term["components"][0]["record"]["canonical_parities"]}[where]
    target["external_geometry"] = [1,2,3]
    with pytest.raises(ValueError, match="extension"):
        reconstruct_output(form, results, source_identity=identity)


def test_identity_validation(session):
    _, requests, _, target, identity = prepared(session, [1], 2)
    bad = deepcopy(target)
    del bad["source_sha256"]
    with pytest.raises(ValueError, match="identity"):
        execution_identity(bad)
    for alteration in ({}, dict(identity, catalog_sha256="bad"),
                       dict(identity, execution={"unexact": 1.2})):
        with pytest.raises(ValueError):
            execution_cache_key(requests[0], alteration)


def test_configuration_boolean_is_not_integer(session):
    form = input_normal_form([1], n=2)
    request = leaf_requests(form)[0]
    result = session.reduce([1], n=2, batch_size=1)
    identity = execution_identity(result)
    identity["execution"]["batch_size"] = True
    with pytest.raises(ValueError, match="producer identity"):
        reconstruct_output(form, {request["key"]: result}, source_identity=identity)


@pytest.mark.parametrize("malformed_request", [None, {}, {"key": "x", "n": True, "word": [1]},
    {"key": "x", "n": 0, "word": []}, {"key": "x", "n": 2, "word": "1"},
    {"key": "x", "n": 2, "word": [True]}, {"key": "x", "n": 2, "word": [0]},
    {"key": "x", "n": 2, "word": [2]}])
def test_invalid_request_shape(session, malformed_request):
    identity = execution_identity(session.reduce([], n=1))
    with pytest.raises(ValueError, match="invalid primitive request"):
        execution_cache_key(malformed_request, identity)


def test_nested_execution_identity_is_exact_and_copied(session):
    form, requests, _, target, identity = prepared(session, [1], 2)
    target = deepcopy(target)
    target["execution"]["external_policy"] = [None, [], {"integers": [1, True], "flag": False}]
    result = execution_identity(target)
    assert result["execution"] == target["execution"]
    target["execution"]["external_policy"].append("changed")
    assert result["execution"] != target["execution"]
    assert execution_cache_key(requests[0], result) != execution_cache_key(requests[0], identity)
    result["execution"]["external_policy"] = {1: "non-string key"}
    with pytest.raises(ValueError, match="exact JSON"):
        execution_cache_key(requests[0], result)


def test_zero_terminal_sum_composes_without_inventing_terms(session):
    form, requests, results, _, identity = prepared(session, [1, 3], 4)
    # Structural adapters cannot certify this zero as a braid evaluation. They
    # must nevertheless preserve an explicitly supplied exact zero formal sum.
    results = deepcopy(results)
    results[requests[0]["key"]]["output"] = []
    result = reconstruct_output(form, results, source_identity=identity)
    assert result["completion"] == "complete" and result["terms"] == []
