"""Genuine device word expansion and exact signed coefficient propagation."""

from copy import deepcopy

import pytest

from zipper_algebra.zipping_algorithm import gpu_rewrite as gpu
from zipper_algebra.zipping_algorithm.catalog import CatalogError
from zipper_algebra.zipping_algorithm.interpreter import Branch, ReplayError


def case(rhs, *, n=16, coefficient=None, controls=None, index=0):
    rule = {
        "id": "test.packed", "condition": True,
        "execution": {"algebra_verified": True},
        "lhs": [{"coefficient": [{"integer": 1, "q_power": 0}],
                 "word": {"type": "sigma", "index": 1}}],
        "rhs": rhs,
    }
    return {
        "root": 0, "id": index, "batch": 0, "n": n, "rule": rule,
        "parameters": {"n": n, "i": 10, "j": 2, "k": 0, "l": 7, "eps": -1},
        "offset": 1, "words": {"W": [-1, 2, -3]},
        "predecessor": Branch((-2, 1, 3), coefficient or ((-7, 2**130), (5, -19)),
                              '{"keep": [1, 2]}').to_json(),
        "successor_controls": controls,
    }


def term(word, pairs=((0, 1),)):
    return {"word": word,
            "coefficient": [{"integer": scalar, "q_power": exponent}
                            for exponent, scalar in pairs]}


def constructor_without_device():
    """Only initialize host-side state; never label this as GPU execution."""
    result = object.__new__(gpu.OpenClRewriteExpander)
    result._templates = []
    result._template_keys = {}
    result._template_arrays = {}
    result._template_buffers = {}
    result._program_term_indices = []
    result._term_factor_indices = []
    result._templates_dirty = False
    return result


def test_failed_template_compilation_does_not_poison_cache():
    expander = constructor_without_device()
    bad = case([term({"type": "unsupported"})])
    good = case([term({"type": "sigma", "index": 2})])
    with pytest.raises(CatalogError):
        expander._pack([bad])
    assert not expander._template_keys and not expander._templates
    expander._pack([good])
    with pytest.raises(CatalogError):
        expander._pack([bad])
    assert len(expander._template_keys) == len(expander._templates) == 1


def test_input_identities_are_required_and_not_merged():
    expander = constructor_without_device()
    good = case([term({"type": "sigma", "index": 2})])
    for invalid in [None, {}, {**good, "root": False}, {**good, "id": "0"}]:
        with pytest.raises(ReplayError, match="identities"):
            expander._pack([invalid])
    with pytest.raises(ReplayError, match="duplicate"):
        expander._pack([good, deepcopy(good)])
    second = deepcopy(good)
    second["root"] = 1
    packed = expander._pack([good, second])
    assert len(packed.admitted) == 2
    assert len(expander._templates) == 1
    huge_unused = deepcopy(good)
    huge_unused["parameters"]["i"] = 2**1000
    assert len(expander._pack([huge_unused]).admitted) == 1


def hardware_cases():
    primitive = {"type": "concat", "items": [
        {"type": "sigma", "index": 7, "sign": -1},
        {"type": "block", "length": 3, "start": 2, "count": 3,
         "order": "descending", "sign": -1},
        {"type": "half_twist", "size": 5, "start": 3},
        {"type": "slider", "size": "i", "start": 2, "orientation": "mirrored"},
        {"type": "word_ref", "name": "W"},
    ]}
    transformed = {"type": "reflect", "word": {"type": "reverse", "word": {
        "type": "bar", "word": primitive}}}
    inverse = {"type": "inverse", "word": {"type": "concat", "items": [
        {"type": "half_twist", "size": 10, "start": 2},
        {"type": "block", "length": 2, "start": 7, "count": 3, "order": "descending"},
        {"type": "fat", "length": 3, "start": 2},
        {"type": "word_ref", "name": "W"},
    ]}}
    values = [
        case([term(transformed, ((-10, -1), (0, 1))), term(inverse)],
             coefficient=((-(2**63) + 10, -(2**32769)), (5, 19))),
        case([
            term({"type": "sigma", "index": 1}),
            term({"type": "sigma", "index": 2}),
            term({"type": "sigma", "index": 1}, ((0, -1), (1, 1))),
        ], controls=[{"label": "B"}, {"label": "A"}], index=1),
        case([term({"type": "identity"})], index=2),
        case([], index=3),
        case([term({"type": "sigma", "index": 2}, ((0, 0),))], index=4),
        case([term({"type": "concat", "items": [
            {"type": "half_twist", "size": 1, "start": 1},
            {"type": "block", "length": 0, "count": 3, "start": 1, "order": "ascending"},
            {"type": "word_ref", "name": "empty"},
            {"type": "slider", "size": 1, "start": 1},
        ]})], index=5),
    ]
    values[-1]["words"]["empty"] = []
    ninth = case([term({"type": "slider", "size": "i", "start": 2})], index=6)
    ninth["parameters"]["i"] = 9
    values.append(ninth)
    values.append(case([
        term({"type": "sigma", "index": 1}, ((0, 0),)),
        term({"type": "sigma", "index": 2}),
        term({"type": "sigma", "index": 1}),
    ], controls=[{"label": "A"}, {"label": "B"}], index=7))
    values.append(case([term({"type": "sigma", "index": 2}, ((1, 0),))],
                       coefficient=(((2**63) - 1, 3),), index=8))
    return values


def test_actual_gpu_expands_words_controls_and_wide_coefficients(monkeypatch):
    try:
        expander = gpu.OpenClRewriteExpander()
    except gpu.GPUUnavailableError as error:
        pytest.skip(str(error))
    cases = hardware_cases()
    expected = gpu.expand_cpu(cases)
    cold = expander.expand(cases)
    warm = expander.expand(cases)
    assert cold.output == warm.output == expected
    assert cold.metadata["new_templates"] > 0
    assert warm.metadata["new_templates"] == 0
    assert cold.metadata["h2d_bytes"] > warm.metadata["h2d_bytes"]
    packed = expander._pack(cases)
    input_payload = (sum(max(len(data), 8) for data in packed.arrays.values())
                     + max(4 * len(cases), 8))
    assert warm.metadata["h2d_bytes"] == input_payload
    assert cold.metadata["h2d_bytes"] == input_payload + cold.metadata["template_upload_bytes"]
    assert cold.metadata["template_upload_happened"]
    assert not warm.metadata["template_upload_happened"]
    assert warm.metadata["template_upload_bytes"] == 0
    assert cold.metadata["emitted_word_tokens"] == cold.metadata["output_word_tokens"]
    assert cold.metadata["emitted_coefficient_records"] == cold.metadata["output_coefficient_monomials"]
    assert cold.metadata["emitted_coefficient_records"] > sum(
        len(branch["coefficient"]) for item in cold.output for branch in item["branches"])
    assert cold.metadata["max_input_coefficient_bits"] == 32770
    assert cold.metadata["verified"] and cold.metadata["backend"] == "opencl"
    assert all(value >= 0 for value in cold.timings.values())
    assert sum(v for k, v in cold.timings.items() if k != "total") <= cold.timings["total"]
    assert [b["control"]["label"] for b in cold.output[1]["branches"]] == ["B", "A"]
    assert cold.output[3]["branches"] == cold.output[4]["branches"] == []
    assert cold.output[5]["branches"][0]["word"][1][0] == "slider"
    assert cold.output[6]["branches"][0]["word"][1][1] == 9
    assert [b["control"]["label"] for b in cold.output[7]["branches"]] == ["A", "B"]
    assert cold.output[8]["branches"] == []
    assert expander.expand([]).output == ()
    with pytest.raises(ValueError, match="boolean"):
        expander.expand(cases, verify=1)
    monkeypatch.setattr(gpu, "expand_cpu", lambda _: ())
    with pytest.raises(ReplayError, match="independent contextual"):
        expander.expand(cases)
    assert expander.expand(cases, verify=False).output == expected
