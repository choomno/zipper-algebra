"""Installed frontend refusals preserve existing certificates on malformed inputs."""
from __future__ import annotations

import json

from zipper_algebra.cli import main
from zipper_algebra.cofactor.artifacts import load_artifact, save_artifact
from zipper_algebra.cofactor.skeleton import Skeleton
from zipper_algebra.presentation import Presentation


def _arguments(tmp_path):
    return ["compose", "--n", "2", "--lift", str(tmp_path / "lift.json"),
            "--skeleton", str(tmp_path / "skeleton.json"),
            "--multipliers", str(tmp_path / "multipliers.json"),
            "--out", str(tmp_path / "certificate.txt"), "--json"]


def _other_inputs(tmp_path):
    for name in ("lift.json", "multipliers.json"):
        (tmp_path / name).write_text("{}")
    (tmp_path / "certificate.txt").write_text("previous verified certificate\n")


def test_frontend_compose_rejects_json_array_without_replacing_output(tmp_path, capsys):
    _other_inputs(tmp_path)
    (tmp_path / "skeleton.json").write_text("[]")
    assert main(_arguments(tmp_path)) == 2
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "invalid-input"
    assert report["error"] == "composition inputs must contain JSON objects"
    assert (tmp_path / "certificate.txt").read_text() == "previous verified certificate\n"


def test_frontend_compose_requires_target_for_valid_bound_artifact(tmp_path, capsys):
    _other_inputs(tmp_path)
    presentation = Presentation(2)
    skeleton = Skeleton(presentation.slug(), presentation.content_hash(), "X_2", (5,7,0), [])
    path = tmp_path / "skeleton.json"
    save_artifact(path, skeleton, presentation, "X_2", order="length-lex-index-v1")
    assert load_artifact(path, presentation, "X_2", kind="skeleton", order="length-lex-index-v1").to_json() == skeleton.to_json()
    assert main(_arguments(tmp_path)) == 2
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "invalid-input"
    assert report["error"] == "bound composition artifacts require an explicit --target"
    assert (tmp_path / "certificate.txt").read_text() == "previous verified certificate\n"
