"""Acceptance and implementation-isolation checks for the independent C++ verifier."""

from copy import deepcopy
from pathlib import Path
import re

import pytest

from zipper_algebra.zipping_algorithm import collected, native, independent_replay
from zipper_algebra.zipping_algorithm.interpreter import Limits
from zipper_algebra.zipping_algorithm.native_verifier import CppVerifier, VerifierError

pytestmark = pytest.mark.timeout(300)


@pytest.fixture(scope="module")
def replay_material():
    with native.NativeSession() as producer:
        complete = producer.reduce(
            [1, -2, 1], n=3, trace=True, limits=Limits(max_steps=10000)
        )
        prefix = producer.reduce(
            [1, -2] * 5, n=3, trace=True, limits=Limits(max_steps=2)
        )
        assert complete["status"] == "complete"
        assert prefix["status"] == "incomplete"
        yield producer.metadata, complete, prefix


def test_verifier_executes_without_producer_or_python_math_callbacks(
    replay_material, monkeypatch
):
    metadata, complete, prefix = replay_material

    def forbidden(*args, **kwargs):
        raise AssertionError(
            "The independent C++ verifier called another implementation"
        )

    with CppVerifier(producer_metadata=metadata) as verifier:
        monkeypatch.setattr(native.NativeSession, "reduce", forbidden)
        monkeypatch.setattr(native.NativeSession, "check", forbidden)
        monkeypatch.setattr(native.NativeSession, "_request", forbidden)
        monkeypatch.setattr(collected.CollectedSession, "reduce", forbidden)
        monkeypatch.setattr(collected, "_transition", forbidden)
        monkeypatch.setattr(
            independent_replay.IndependentVerifier, "check_native", forbidden
        )
        assert verifier.check(complete) is True
        assert verifier.check(prefix) is True


@pytest.mark.parametrize(
    "field",
    [
        "output",
        "pending",
        "metrics",
        "transition",
        "scope",
        "extra",
        "n_float",
        "word_float",
        "trace_integer",
        "limit_float",
        "batch_float",
        "coefficient_string",
        "coefficient_bool",
    ],
)
def test_warm_verifier_rejects_tampering_without_poisoning_cache(
    replay_material, field
):
    metadata, complete, prefix = replay_material
    original = prefix if field == "pending" else complete
    bad = deepcopy(original)
    if field == "output":
        bad["output"][0]["coefficient"][0][1] += 1
    elif field == "pending":
        bad["pending"][0]["depth"] += 1
    elif field == "metrics":
        bad["metrics"]["steps"] += 1
    elif field == "transition":
        bad["transitions"][0]["predecessors"][0]["depth"] += 1
    elif field == "scope":
        bad["scope"] = "unproved_arbitrary_theorem"
    elif field == "extra":
        bad["unexpected"] = True
    elif field == "n_float":
        bad["n"] = float(bad["n"])
    elif field == "word_float":
        bad["input_word"][0] = float(bad["input_word"][0])
    elif field == "trace_integer":
        bad["trace_enabled"] = 1
    elif field == "limit_float":
        bad["limits"]["max_steps"] = float(bad["limits"]["max_steps"])
    elif field == "batch_float":
        bad["execution"]["batch_size"] = float(bad["execution"]["batch_size"])
    elif field == "coefficient_string":
        bad["output"][0]["coefficient"][0][1] = str(
            bad["output"][0]["coefficient"][0][1]
        )
    else:
        bad["output"][0]["coefficient"][0][1] = True
    with CppVerifier(producer_metadata=metadata) as verifier:
        assert verifier.check(original) is True
        with pytest.raises(VerifierError):
            verifier.check(bad)
        assert verifier.check(original) is True


def test_independent_sources_include_no_producer_math_headers():
    root = Path(__file__).resolve().parents[1] / "native" / "zipper_verifier"
    for path in sorted(root.glob("*")):
        if path.suffix not in (".cpp", ".hpp", ".h"):
            continue
        text = path.read_text()
        includes = re.findall(r'^\s*#\s*include\s*[<"]([^>"\n]+)', text, re.M)
        assert all("zipper_engine" not in item for item in includes), (path, includes)
        assert "zipper::" not in text, path
