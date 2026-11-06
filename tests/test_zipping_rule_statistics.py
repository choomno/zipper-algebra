"""Rule frequencies describe semantic work, including incomplete prefixes."""

import csv
import json
import sqlite3

import pytest

from experiments import zipper_rule_statistics as stats


def database(tmp_path, rows):
    path = tmp_path / "words.sqlite3"
    with sqlite3.connect(path) as db:
        db.execute("CREATE TABLE metadata(key TEXT,value TEXT)")
        db.execute("INSERT INTO metadata VALUES('identity','{\"source\":\"fixture\"}')")
        db.execute("CREATE TABLE ranks(n INTEGER,length INTEGER,ordinal INTEGER)")
        db.execute("INSERT INTO ranks VALUES(3,2,2)")
        db.execute("CREATE TABLE words(n INTEGER,length INTEGER,ordinal INTEGER,status TEXT,metrics TEXT,certificate_sha TEXT)")
        for ordinal, status, steps, applications in rows:
            db.execute("INSERT INTO words VALUES(3,2,?,?,?,?)", (
                ordinal, status, json.dumps({"steps": steps, "rule_applications": applications}),
                f"certificate-{ordinal}",
            ))
    return path


def test_exact_frequency_denominators_and_incomplete_prefix(tmp_path):
    path = database(tmp_path, [
        (0, "complete", 5, {"rule.a": 4, "rule.b": 1, "unused": 0}),
        (1, "incomplete", 3, {"rule.b": 3}),
    ])
    original = path.read_bytes()
    report = stats.statistics(path)
    assert path.read_bytes() == original
    assert report["inputs"] == 2
    assert report["statuses"] == {"complete": 1, "incomplete": 1}
    assert report["semantic_term_transitions"] == report["rule_application_sum"] == 8
    assert report["rules"] == [
        {"rule_id": "rule.a", "semantic_applications": 4,
         "share_of_semantic_applications": .5, "inputs_using_rule": 1,
         "applications_per_input": 2},
        {"rule_id": "rule.b", "semantic_applications": 4,
         "share_of_semantic_applications": .5, "inputs_using_rule": 2,
         "applications_per_input": 2},
    ]
    assert report["actual_evaluator_calls"] is report["cache_reuse"] is None
    assert report["by_rank_length"][0]["rules"] == report["rules"]
    assert report["by_rank_length"][0]["expected_words"] == 16
    assert report["by_rank_length"][0]["complete_words"] == 1
    assert report["by_rank_length"][0]["visited_fraction"] == .125
    assert report["by_rank_length"][0]["layer_complete"] is False
    csv_path = tmp_path / "statistics.csv"
    stats.write_csv(report, csv_path)
    with csv_path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 4
    assert [row["scope"] for row in rows] == ["all", "all", "rank_length", "rank_length"]
    assert rows[2]["n"] == "3" and rows[2]["length"] == "2"


@pytest.mark.parametrize("steps,applications", [
    (2, {"rule.a": 1}), (0, {"rule.a": -1}), (1, {"rule.a": True}),
    (False, {}), (0, []),
])
def test_malformed_or_inconsistent_metrics_are_rejected(tmp_path, steps, applications):
    path = database(tmp_path, [(0, "complete", steps, applications)])
    with pytest.raises(ValueError, match="semantic transition metrics"):
        stats.statistics(path)


def test_empty_campaign(tmp_path):
    report = stats.statistics(database(tmp_path, []))
    assert report["inputs"] == report["semantic_term_transitions"] == 0
    assert report["rules"] == report["by_rank_length"] == []


def test_catalog_only_distribution_excludes_control_and_transport(tmp_path):
    report = stats.statistics(database(tmp_path, [(0, "complete", 16, {
        "control.extract": 10, "layer.route.left": 2,
        "zip.internal.i2.negative.I": 3, "single.foe.left.positive.I": 1,
    })]))
    for group in (report, report["by_rank_length"][0]):
        categories = group["rule_categories"]
        assert categories["control_actions"]["semantic_applications"] == 10
        assert categories["transport_rewrites"]["semantic_applications"] == 2
        algebraic = categories["catalog_algebraic_rewrites"]
        assert algebraic["semantic_applications"] == 4
        assert algebraic["share_of_all_semantic_applications"] == .25
        assert [row["share_of_category_semantic_applications"] for row in algebraic["top10"]] == [.75, .25]
        assert categories["unclassified"]["rules"] == []


def test_mixed_source_versions_remain_visible_in_statistics(tmp_path):
    path = database(tmp_path, [(0, "complete", 1, {"control.extract": 1}),
                               (1, "complete", 1, {"control.extract": 1})])
    with sqlite3.connect(path) as db:
        db.execute("CREATE TABLE versions(hash TEXT,identity TEXT)")
        db.executemany("INSERT INTO versions VALUES(?,?)", [("old", '{"source":"old"}'), ("new", '{"source":"new"}')])
        db.execute("INSERT INTO metadata VALUES('active_version','new')")
        db.execute("CREATE TABLE configurations(hash TEXT,config TEXT)")
        db.executemany("INSERT INTO configurations VALUES(?,?)", [("a", '{"identity_sha256":"old"}'), ("b", '{"identity_sha256":"new"}')])
        db.execute("ALTER TABLE words ADD COLUMN config_sha TEXT")
        db.execute("UPDATE words SET config_sha=CASE WHEN ordinal=0 THEN 'a' ELSE 'b' END")
    report = stats.statistics(path)
    assert report["identity"] == {"source": "fixture"}
    assert report["active_version"] == "new"
    assert report["versions"] == {"old": {"source": "old"}, "new": {"source": "new"}}
    assert report["by_rank_length"][0]["version_counts"] == {"old": 1, "new": 1}
