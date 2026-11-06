"""The standalone command preserves scope and detects malformed saved data."""

from __future__ import annotations

import json
import subprocess
import sys

import pytest

from zipper_algebra.zipping_algorithm.__main__ import main
from zipper_algebra.zipping_algorithm.interpreter import Branch, Request, RunState


def invoke(capsys, args, expected=0):
    assert main(args) == expected
    captured = capsys.readouterr()
    if expected in (2, 3) and not captured.out:
        return json.loads(captured.err)
    assert not captured.err
    return json.loads(captured.out)


def save(tmp_path, data):
    path = tmp_path / "input.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return str(path)


def replay_fixture():
    state = RunState(4, (Branch((3, ("slider", 1, 1, False, False, False))),) * 2)
    requests = [
        Request(idx, 0, "single.far.positive.I", {"n": 4, "i": 1, "j": 1, "l": 3})
        for idx in range(2)
    ]
    return state, requests


def test_collected_run_replay_and_checkpoint(capsys, tmp_path):
    args = ["run-collected", "--n", "3", "--word", "[1,2]", "--batch-size", "8"]
    certificate = invoke(capsys, args + ["--trace"])
    assert certificate["schema"] == "zipper-collected-run-v1"
    assert certificate["status"] == "complete"
    assert invoke(capsys, ["check", save(tmp_path, certificate)])["pass"]
    checkpoint = invoke(capsys, args + ["--limits", '{"max_steps":1}'], 3)
    assert checkpoint["status"] == "incomplete"
    assert invoke(capsys, ["check", save(tmp_path, checkpoint)])["pass"]
    assert invoke(capsys, args + ["--workers", "0"], 2)


def test_instantiate(capsys):
    parameters = {"n": 5, "i": 3, "j": 1, "k": 2, "l": 1}
    value = invoke(
        capsys,
        [
            "instantiate",
            "zip.nearest.negative.i_ge_3.I",
            "--parameters",
            json.dumps(parameters),
        ],
    )
    assert value["algebra_verified"]
    assert value["parameters"] == parameters
    assert value["lhs"][0]["word"][-1] == ["slider", 3, 1, False, False, False]
    expanded = invoke(
        capsys,
        [
            "instantiate",
            "single.foe.right.positive.I",
            "--parameters",
            '{"n":2,"i":1,"j":1,"l":1}',
            "--expand-sliders",
        ],
    )
    assert all(isinstance(t, int) for term in expanded["rhs"] for t in term["word"])


def test_replay_complete_and_checkpoint(capsys, tmp_path):
    state, requests = replay_fixture()
    specification = {
        "state": state.to_json(),
        "requests": [x.to_json() for x in requests],
    }
    path = save(tmp_path, specification)
    certificate = invoke(capsys, ["replay", path])
    assert certificate["status"] == "schedule_exhausted"
    assert len(certificate["transitions"]) == 2
    assert invoke(capsys, ["check", save(tmp_path, certificate)]) == {
        "pass": True,
        "scope": "algebraic_catalog_replay",
    }
    specification["limits"] = {"max_steps": 1}
    checkpoint = invoke(capsys, ["replay", save(tmp_path, specification)], 3)
    assert checkpoint["status"] == "incomplete"
    assert invoke(capsys, ["check", save(tmp_path, checkpoint)])["pass"]


def test_upper_complete_and_checkpoint(capsys, tmp_path):
    args = ["upper", "--n", "5", "--size", "2", "--start", "1", "--friends", "0"]
    certificate = invoke(capsys, args + ["--left", "[]"])
    assert certificate["status"] == "upper_pass_complete"
    assert invoke(capsys, ["check", save(tmp_path, certificate)]) == {
        "pass": True,
        "scope": "normal_upward_pass_only",
    }
    checkpoint = invoke(capsys, args + ["--left", "[4,4]", "--max-steps", "1"], 3)
    assert checkpoint["status"] == "incomplete"


def test_two_strand_complete_and_checkpoint(capsys, tmp_path):
    certificate = invoke(capsys, ["two-strand", "--word", "[1,-1,1]"])
    assert certificate["status"] == "complete"
    assert invoke(capsys, ["check", save(tmp_path, certificate)])["pass"]
    checkpoint = invoke(
        capsys, ["two-strand", "--word", "[1,1]", "--max-steps", "1"], 3
    )
    assert checkpoint["status"] == "incomplete"
    assert invoke(capsys, ["check", save(tmp_path, checkpoint)])["pass"]


def test_transport_and_standalone_check(capsys, tmp_path):
    certificate = invoke(
        capsys,
        [
            "transport",
            "--n",
            "3",
            "--size",
            "2",
            "--start",
            "1",
            "--word",
            "[-2,-1]",
            "--mirrored",
        ],
    )
    assert certificate["witness"]["end"] == 2
    assert certificate["witness"]["word"] == [-2, -1]
    assert invoke(capsys, ["check", save(tmp_path, certificate)]) == {
        "pass": True,
        "scope": "whole_band_transport",
    }
    failure = invoke(
        capsys,
        ["transport", "--n", "4", "--size", "3", "--start", "1", "--word", "[1]"],
        3,
    )
    assert failure["status"] == "unsupported"
    assert "offset 0" in failure["diagnostic"]


@pytest.mark.parametrize(
    "payload",
    [
        None,
        {},
        {"schema": "unknown"},
        {"schema": "zipper-band-transport-v1", "witness": {}, "extra": 1},
        {"schema": "zipper-band-transport-v1", "witness": {}},
        {"schema": "zipper-band-transport-v1", "witness": {"word": [], "steps": None}},
        {
            "schema": "zipper-band-transport-v1",
            "witness": {"word": [], "steps": [None]},
        },
    ],
)
def test_invalid_certificate(capsys, tmp_path, payload):
    result = invoke(capsys, ["check", save(tmp_path, payload)], 2)
    assert result["status"] == "invalid"
    assert result["diagnostic"]


@pytest.mark.parametrize(
    "payload",
    [
        {},
        {"state": {}, "requests": [], "extra": True},
        {"state": {}, "requests": []},
        {"state": None, "requests": []},
    ],
)
def test_invalid_replay(capsys, tmp_path, payload):
    assert invoke(capsys, ["replay", save(tmp_path, payload)], 2)["status"] == "invalid"


def test_json_file_and_parameter_errors(capsys, tmp_path):
    path = tmp_path / "malformed.json"
    path.write_text("{", encoding="utf-8")
    for args in [
        ["check", str(path)],
        ["check", str(tmp_path / "missing.json")],
        ["instantiate", "unknown", "--parameters", "{}"],
        ["instantiate", "single.far.positive.I", "--parameters", "[]"],
        ["instantiate", "single.far.positive.I", "--parameters", "{"],
        ["two-strand", "--word", "[2]"],
    ]:
        assert invoke(capsys, args, 2)["status"] == "invalid"


def test_serialization_error_is_reported(capsys):
    params = '{"n":4,"i":1,"j":1,"l":3,"extra":NaN}'
    assert (
        invoke(
            capsys, ["instantiate", "single.far.positive.I", "--parameters", params], 2
        )["status"]
        == "invalid"
    )


def test_argument_parser_errors(capsys):
    for args in [[], ["unknown"], ["two-strand"], ["upper", "--n", "bad"]]:
        with pytest.raises(SystemExit) as error:
            main(args)
        assert error.value.code == 2
        assert "error:" in capsys.readouterr().err


def test_module_entrypoint():
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "zipper_algebra.zipping_algorithm",
            "two-strand",
            "--word",
            "[]",
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["status"] == "complete"


@pytest.mark.parametrize(
    "args",
    [
        ["two-strand", "--word", "[1]"],
        ["upper", "--n", "5", "--size", "2", "--start", "1", "--left", "[4]"],
    ],
)
def test_optional_trace(capsys, tmp_path, args):
    compact = invoke(capsys, args)
    detailed = invoke(capsys, args + ["--trace"])
    assert compact["trace_enabled"] is False
    assert compact["transitions"] == []
    assert detailed["trace_enabled"] is True
    assert detailed["transitions"]
    assert compact["status"] == detailed["status"]
    assert invoke(capsys, ["check", save(tmp_path, compact)])["pass"]
    assert invoke(capsys, ["check", save(tmp_path, detailed)])["pass"]


def test_verify_local_rule_and_reload(capsys, tmp_path):
    certificate = invoke(
        capsys,
        [
            "verify-rule",
            "single.foe.right.positive.I",
            "--parameters",
            '{"n":2,"i":1,"j":1,"l":1}',
        ],
    )
    assert certificate["schema"] == "zipper-local-original-relator-v1"
    assert invoke(capsys, ["check", save(tmp_path, certificate)])["pass"]
    failure = invoke(
        capsys,
        [
            "verify-rule",
            "single.foe.right.positive.I",
            "--parameters",
            '{"n":2,"i":1,"j":1,"l":1}',
            "--limits",
            '{"max_operations":1}',
        ],
        3,
    )
    assert failure["status"] == "incomplete"


def test_exhaust_literal_words(capsys, tmp_path):
    args = [
        "exhaust",
        str(tmp_path / "literal"),
        "--ranks",
        "[1,2]",
        "--seconds",
        "10",
        "--max-length",
        "1",
    ]
    summary = invoke(capsys, args)
    assert summary["stop_reason"] == "max_length"
    assert summary["new_or_resumed_words"] == 4
    resumed = invoke(capsys, args)
    assert resumed["audited_words"] == 4
    assert resumed["new_or_resumed_words"] == 0
    unsupported = invoke(
        capsys,
        ["exhaust", str(tmp_path / "unsupported"), "--ranks", "[3]", "--seconds", "1"],
        3,
    )
    assert unsupported["unsupported_ranks"] == [3]
    checkpoint = invoke(
        capsys,
        ["exhaust", str(tmp_path / "checkpoint"), "--ranks", "[2]", "--max-steps", "0"],
        3,
    )
    assert checkpoint["stop_reason"] == "word_limit"
    assert (
        invoke(capsys, ["exhaust", str(tmp_path / "bad"), "--ranks", "[]"], 2)["status"]
        == "invalid"
    )


def test_general_run_compact_trace_and_checkpoint(capsys, tmp_path):
    args = ["run", "--n", "3", "--word", "[1,2]"]
    compact = invoke(capsys, args)
    traced = invoke(capsys, args + ["--trace"])
    assert compact["status"] == traced["status"] == "complete"
    assert compact["input_word"] == [1, 2]
    assert compact["output"] == traced["output"]
    assert compact["transitions"] == []
    assert traced["transitions"]
    assert compact["source_sha256"]["algorithm.py"]
    assert compact["catalog_sha256"]
    for certificate in (compact, traced):
        assert invoke(capsys, ["check", save(tmp_path, certificate)]) == {
            "pass": True,
            "scope": "literal_word_spanning_procedure",
        }
    checkpoint = invoke(capsys, args + ["--limits", '{"max_steps":2}'], 3)
    assert checkpoint["status"] == "incomplete"
    assert checkpoint["pending"]
    assert invoke(capsys, ["check", save(tmp_path, checkpoint)])["pass"]
    assert (
        invoke(capsys, ["run", "--n", "3", "--word", "[3]"], 2)["status"] == "invalid"
    )


def test_general_run_four_strand_staircase(capsys, tmp_path):
    certificate = invoke(capsys, ["run", "--n", "4", "--word", "[1,2,3]"])
    assert certificate["status"] == "complete"
    assert certificate["input_word"] == [1, 2, 3]
    assert certificate["metrics"]["steps"] > 3
    assert invoke(capsys, ["check", save(tmp_path, certificate)])["pass"]


def test_general_literal_exhaustion(capsys, tmp_path):
    args = [
        "exhaust-general",
        str(tmp_path / "general"),
        "--ranks",
        "[3]",
        "--seconds",
        "10",
        "--max-length",
        "1",
        "--trace",
    ]
    summary = invoke(capsys, args)
    assert summary["stop_reason"] == "max_length"
    assert summary["new_or_resumed_words"] == 5
    assert summary["identity"]["trace"] is True
    assert summary["identity"]["algorithm_source_sha256"]
    resumed = invoke(capsys, args)
    assert resumed["audited_words"] == 5
    assert resumed["new_or_resumed_words"] == 0
    checkpoint = invoke(
        capsys,
        [
            "exhaust-general",
            str(tmp_path / "general-limit"),
            "--ranks",
            "[3]",
            "--seconds",
            "10",
            "--max-length",
            "1",
            "--limits",
            '{"max_steps":1}',
        ],
        3,
    )
    assert checkpoint["stop_reason"] == "word_limit"
