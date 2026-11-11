"""Native C++ packed expansion remains exact and source-bound without a GPU."""

from copy import deepcopy
import json
from pathlib import Path
import shutil
import struct

import pytest

from experiments import zipper_native_rewrite as native
from zipper_algebra.zipping_algorithm.gpu_rewrite import expand_cpu
from zipper_algebra.zipping_algorithm.interpreter import Branch, ReplayError


def case(word, *, factors=((0, 1),), coefficient=None, index=0, controls=None):
    return {
        "root": 2, "id": index, "batch": 0, "n": 32,
        "rule": {"id": "native.exact", "condition": True,
                 "execution": {"algebra_verified": True},
                 "lhs": [{"coefficient": [{"integer": 1, "q_power": 0}],
                          "word": {"type": "sigma", "index": 1}}],
                 "rhs": [{"word": word, "coefficient": [
                     {"integer": c, "q_power": q} for q, c in factors]}]},
        "parameters": {"n": 32, "i": 10, "j": 2, "k": 0, "l": 4},
        "words": {"W": [-2, 3]}, "offset": 1,
        "predecessor": Branch((-1, 1, 2), coefficient or ((0, 1),),
                              '{"preserved": [3]}').to_json(),
        "successor_controls": controls,
    }


@pytest.fixture(scope="module")
def expander():
    if shutil.which("g++") is None:
        pytest.skip("C++ compiler is unavailable")
    return native.NativeRewriteExpander()


def test_exact_wide_integer_context_transformations_and_retained_templates(expander):
    word = {"type": "concat", "items": [
        {"type": "inverse", "word": {"type": "half_twist", "size": 10, "start": 2}},
        {"type": "block", "length": 3, "start": 2, "count": 4, "order": "descending"},
        {"type": "reflect", "word": {"type": "reverse", "word": {
            "type": "bar", "word": {"type": "concat", "items": [
                {"type": "slider", "size": "i", "start": 2},
                {"type": "word_ref", "name": "W"}]}}}},
    ]}
    original = case(word, factors=((-4, -1), (7, 1)),
                    coefficient=((-100, 1 << 32769), (4, -(1 << 32769) + 3)))
    snapshot = deepcopy(original)
    result = expander.expand([original])
    assert result.output == expand_cpu([original])
    assert original == snapshot
    assert result.metadata["backend"] == "cpp" and result.metadata["threads"] == 1
    assert result.metadata["max_input_coefficient_bits"] == 32770
    assert result.metadata["h2d_bytes"] == 0
    assert not any(k in result.timings for k in ("h2d_and_allocation", "d2h"))
    following = deepcopy(original)
    following["parameters"]["i"] = 24
    following["id"] += 1
    second = expander.expand([following])
    assert second.output == expand_cpu([following])
    assert second.metadata["new_templates"] == 0
    assert expander.expand([]).output == ()


def test_literal_zero_resurrection_and_control_order(expander):
    item = case({"type": "sigma", "index": 2}, controls=[{"order": "A"}, {"order": "B"}])
    item["rule"]["rhs"] = [
        {"word": {"type": "sigma", "index": 2}, "coefficient": [{"integer": 0, "q_power": 0}]},
        {"word": {"type": "sigma", "index": 3}, "coefficient": [{"integer": 1, "q_power": 0}]},
        {"word": {"type": "sigma", "index": 2}, "coefficient": [{"integer": -1, "q_power": 0}]},
    ]
    result = expander.expand([item])
    assert result.output == expand_cpu([item])
    assert [b["word"] for b in result.output[0]["branches"]] == [[-1, 2, 2], [-1, 3, 2]]
    assert [b["control"]["order"] for b in result.output[0]["branches"]] == ["A", "B"]


def test_device_contract_status_is_checked_for_native_overflow(expander):
    item = case({"type": "sigma", "index": 2}, factors=((1, 1),))
    packed = expander._pack([item])
    packed.arrays["old_q"] = struct.pack("<q", (1 << 63) - 1)
    with pytest.raises(ReplayError, match="checked arithmetic"):
        expander._execute(packed)


def test_compiler_cache_and_manifest_are_bound_to_binary(expander, tmp_path):
    path = Path(expander.metadata["library"])
    record = native._manifest(path)
    assert record["abi"] == native.ABI
    copy = tmp_path / path.name
    shutil.copyfile(path, copy)
    manifest = copy.with_suffix(".json")
    manifest.write_text(json.dumps(record))
    loaded = native.NativeRewriteExpander(library=copy)
    assert loaded.metadata["binary_sha256"] == record["binary_sha256"]
    record["binary_sha256"] = "wrong"
    manifest.write_text(json.dumps(record))
    with pytest.raises(native.NativeBuildError, match="provenance differs"):
        native.NativeRewriteExpander(library=copy)


def test_missing_compiler_and_manifest_fail_explicitly(tmp_path):
    with pytest.raises(native.NativeBuildError, match="compiler is unavailable"):
        native.build_library(compiler="zipper-no-such-compiler")
    with pytest.raises(native.NativeBuildError, match="valid matching provenance"):
        native.NativeRewriteExpander(library=tmp_path / "missing.so")


def test_cold_and_warm_native_template_metadata_never_claims_gpu_uploads(expander):
    # A fresh instance retains no program cache, even though its binary is cached.
    fresh = native.NativeRewriteExpander(library=Path(expander.metadata["library"]))
    item = case({"type": "sigma", "index": 2})
    cold = fresh.expand([item])
    warm = fresh.expand([item])
    assert cold.metadata["new_templates"] == 1
    assert warm.metadata["new_templates"] == 0
    for result in (cold, warm):
        assert result.metadata["template_upload_happened"] is False
        assert result.metadata["template_upload_bytes"] == 0
        assert result.metadata["h2d_bytes"] == result.metadata["d2h_bytes"] == 0
        assert result.metadata["native_template_buffer_bytes"] > 0
        assert result.metadata["native_output_bytes"] > 0
        assert "no host-to-device copies" in result.metadata["h2d_byte_measure"]
    assert cold.output == warm.output == expand_cpu([item])
