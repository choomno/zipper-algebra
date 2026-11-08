"""Fresh protocol parsing keeps exact coefficients and independent mutable trees."""
from copy import deepcopy
import io
import json
from types import SimpleNamespace

import pytest

from zipper_algebra.zipping_algorithm import native
from zipper_algebra.zipping_algorithm.native_verifier import VerifierError, _validate_candidate


def decode(payload):
    return json.loads(payload, object_hook=native._decode_coefficient_object)


def test_nested_coefficients_extensions_and_arbitrary_width():
    digits = "9" * 10000
    value = {
        "output": [{"coefficient": [[-3, "-" + digits]], "word": [1, -1],
                    "extension": {"coefficient": [], "label": "123"}}],
        "transitions": [{"nested": {"coefficient": [[0, digits]]}}],
        "catalog": {"coefficient": {"integer": -1, "q_power": "i"}},
    }
    payload = json.dumps(value)
    result = decode(payload)
    assert result == native._decode_coefficients(value)
    assert result["output"][0]["coefficient"] == [[-3, 1 - 10**10000]]
    assert result["transitions"][0]["nested"]["coefficient"] == [[0, 10**10000 - 1]]
    assert result["catalog"] == value["catalog"]
    assert result["output"][0]["extension"] == {"coefficient": [], "label": "123"}
    assert value["output"][0]["coefficient"][0][1] == "-" + digits


@pytest.mark.parametrize("scalar", [True, 1.0, None, {}, [], "", "-", "+1", "1.0", "１２"])
def test_malformed_scalar_is_not_coerced(scalar):
    with pytest.raises(native.NativeError, match="Malformed native integer coefficient"):
        decode(json.dumps({"coefficient": [[0, scalar]]}))


@pytest.mark.parametrize("exponent", [True, 0.0])
def test_noninteger_exponent_remains_invalid_for_verification(exponent):
    result = decode(json.dumps({"coefficient": [[exponent, "2"]]}))
    assert type(result["coefficient"][0][0]) is type(exponent)
    assert result["coefficient"][0][1] == "2"
    with pytest.raises(VerifierError, match="integer scalars and exponents"):
        _validate_candidate(result)


def test_repeated_loads_and_generic_decoder_keep_independent_trees():
    value = {"coefficient": [[0, "12"]], "child": [{"coefficient": [[1, "-2"]]}]}
    saved = deepcopy(value)
    first = decode(json.dumps(value))
    second = decode(json.dumps(value))
    generic = native._decode_coefficients(value)
    first["child"][0]["coefficient"][0][1] = 99
    second["coefficient"][0][1] = 77
    assert generic == {"coefficient": [[0, 12]], "child": [{"coefficient": [[1, -2]]}]}
    assert first["coefficient"] == [[0, 12]]
    assert second["child"][0]["coefficient"] == [[1, -2]]
    assert value == saved


def test_discarded_duplicate_key_subtree_can_fail_closed():
    # Hooks visit the first child before the enclosing duplicate key discards it.
    # Authentic native and stored canonical JSON do not contain duplicate keys.
    payload = '{"extension":{"coefficient":[[0,"bad"]]},"extension":null}'
    assert native._decode_coefficients(json.loads(payload)) == {"extension": None}
    with pytest.raises(native.NativeError, match="Malformed"):
        decode(payload)


def test_protocol_decoding_is_opt_in_and_produces_fresh_results():
    session = object.__new__(native.NativeSession)
    line = '{"output":[{"coefficient":[[0,"12"]]}]}\n'
    session._process = SimpleNamespace(stdin=io.StringIO(), stdout=io.StringIO(line * 3))
    session._stderr = None
    try:
        raw = session._request({"op": "diagnostic"})
        first = session._request({"op": "reduce"}, decode_coefficients=True)
        second = session._request({"op": "reduce_reused"}, decode_coefficients=True)
        assert raw["output"][0]["coefficient"] == [[0, "12"]]
        assert first == second == {"output": [{"coefficient": [[0, 12]]}]}
        first["output"][0]["coefficient"][0][1] = 99
        assert second["output"][0]["coefficient"] == [[0, 12]]
    finally:
        session._process = None


def test_malformed_coefficient_closes_native_protocol():
    actions = []
    process = SimpleNamespace(
        stdin=io.StringIO(), stdout=io.StringIO('{"coefficient":[[0,true]]}\n'),
        kill=lambda: actions.append("kill"), wait=lambda timeout: actions.append("wait"),
    )
    session = object.__new__(native.NativeSession)
    session._process, session._stderr = process, None
    with pytest.raises(native.NativeError, match="Malformed"):
        session._request({"op": "reduce"}, decode_coefficients=True)
    assert session._process is None and actions == ["kill", "wait"]
    assert process.stdin.closed and process.stdout.closed
    with pytest.raises(native.NativeError, match="closed"):
        session._request({"op": "reduce"}, decode_coefficients=True)
