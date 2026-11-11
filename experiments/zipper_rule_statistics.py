"""Count committed semantic zipper transitions from one SQLite snapshot.

Read the collected exhaustive campaign's saved metrics without replaying its
arithmetic. Counts are unweighted by Laurent coefficients and include the
committed prefix of incomplete calculations. Rule frequency is not a timing
measurement. Evaluator calls and cache reuse require separate instrumentation.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import time


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _histogram(counts, present, inputs):
    total = sum(counts.values())
    return [
        {
            "rule_id": rule,
            "semantic_applications": count,
            "share_of_semantic_applications": count / total if total else 0,
            "inputs_using_rule": present[rule],
            "applications_per_input": count / inputs if inputs else 0,
        }
        for rule, count in sorted(counts.items(), key=lambda pair: (-pair[1], pair[0]))
    ]


def rule_category(rule):
    """Separate control actions, named routing moves, and catalog rule IDs."""
    if rule.startswith("control."):
        return "control_actions"
    if rule.startswith("layer.route."):
        return "transport_rewrites"
    if rule.startswith(("single.", "zip.")):
        return "catalog_algebraic_rewrites"
    return "unclassified"


def rule_categories(group):
    """Recompute shares within each category, retaining the full denominator."""
    result = {}
    for category in ("control_actions", "transport_rewrites", "catalog_algebraic_rewrites", "unclassified"):
        rows = [row for row in group["rules"] if rule_category(row["rule_id"]) == category]
        total = sum(row["semantic_applications"] for row in rows)
        rules = [dict(row, share_of_category_semantic_applications=row["semantic_applications"] / total)
                 for row in rows]
        result[category] = {
            "semantic_applications": total,
            "share_of_all_semantic_applications": total / group["semantic_term_transitions"]
            if group["semantic_term_transitions"] else 0,
            "applications_per_input": total / group["inputs"] if group["inputs"] else 0,
            "rules": rules,
            "top10": rules[:10],
        }
    return result


def statistics(database):
    """Read a consistent snapshot; do not acquire the campaign writer lock."""
    database = Path(database).resolve()
    started = time.perf_counter()
    counts, present, statuses = Counter(), Counter(), Counter()
    cells = {}
    inputs = steps = 0
    index = hashlib.sha256()
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as db:
        db.execute("BEGIN")
        metadata = dict(db.execute("SELECT key,value FROM metadata"))
        identity = json.loads(metadata["identity"])
        version_counts = {}
        versions = {}
        if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='versions'").fetchone():
            versions = {digest: json.loads(encoded) for digest, encoded in db.execute("SELECT hash,identity FROM versions")}
            for n, length, version, count in db.execute(
                "SELECT w.n,w.length,json_extract(c.config,'$.identity_sha256'),count(*) "
                "FROM words w JOIN configurations c ON w.config_sha=c.hash "
                "GROUP BY w.n,w.length,json_extract(c.config,'$.identity_sha256')"
            ):
                version_counts.setdefault((n, length), {})[version] = count
        snapshot = datetime.now(timezone.utc).isoformat()
        for n, length, ordinal, status, encoded, digest in db.execute(
            "SELECT n,length,ordinal,status,metrics,certificate_sha "
            "FROM words ORDER BY n,length,ordinal"
        ):
            metrics = json.loads(encoded)
            applications = metrics["rule_applications"]
            row_steps = metrics["steps"]
            if (
                status not in ("complete", "incomplete", "control_obligation")
                or type(row_steps) is not int
                or row_steps < 0
                or not isinstance(applications, dict)
                or any(type(value) is not int or value < 0 for value in applications.values())
                or sum(applications.values()) != row_steps
            ):
                raise ValueError(f"invalid semantic transition metrics at {(n, length, ordinal)}")
            applications = {rule: value for rule, value in applications.items() if value}
            inputs += 1
            steps += row_steps
            statuses[status] += 1
            counts.update(applications)
            present.update(applications.keys())
            cell = cells.setdefault(
                (n, length),
                {"inputs": 0, "steps": 0, "statuses": Counter(),
                 "counts": Counter(), "present": Counter()},
            )
            cell["inputs"] += 1
            cell["steps"] += row_steps
            cell["statuses"][status] += 1
            cell["counts"].update(applications)
            cell["present"].update(applications.keys())
            index.update((_json([n, length, ordinal, digest]) + "\n").encode())
        cursors = [
            {"n": n, "length": length, "ordinal": ordinal}
            for n, length, ordinal in db.execute("SELECT n,length,ordinal FROM ranks ORDER BY n")
        ]
    report = {
        "schema": "zipper-rule-statistics-v1",
        "database": str(database),
        "snapshot_utc": snapshot,
        "identity": identity,
        "versions": versions,
        "active_version": metadata.get("active_version"),
        "analysis_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "inputs": inputs,
        "statuses": dict(statuses),
        "semantic_term_transitions": steps,
        "rule_application_sum": sum(counts.values()),
        "rules": _histogram(counts, present, inputs),
        "primary_comparison_scope": "rank_length; pooled rules summarize this workload snapshot only",
        "by_rank_length": [
            {"n": n, "length": length, "inputs": cell["inputs"],
             "expected_words": (2 * (n - 1)) ** length,
             "complete_words": cell["statuses"]["complete"],
             "visited_fraction": cell["inputs"] / (2 * (n - 1)) ** length,
             "layer_complete": cell["inputs"] == cell["statuses"]["complete"] == (2 * (n - 1)) ** length,
             "version_counts": version_counts.get((n, length), {}),
             "statuses": dict(cell["statuses"]), "semantic_term_transitions": cell["steps"],
             "rules": _histogram(cell["counts"], cell["present"], cell["inputs"])}
            for (n, length), cell in sorted(cells.items())
        ],
        "cursors": cursors,
        "certificate_index_sha256": index.hexdigest(),
        "index_encoding": "UTF-8 compact JSON [n,length,ordinal,certificate_sha] plus LF, ordered by n,length,ordinal",
        "actual_evaluator_calls": None,
        "cache_reuse": None,
        "interpretation": (
            "Committed semantic term transitions, unweighted by Laurent coefficients. "
            "Incomplete calculations contribute their committed prefix. Rule frequencies "
            "do not measure per-rule runtime or identify evaluator calls/cache reuse. "
            "Compare rank/length strata with their coverage; partial literal layers are lexicographic prefixes. "
            "One read-only consistent transaction; no arithmetic replay or certificate re-audit."
        ),
        "elapsed_seconds": time.perf_counter() - started,
    }
    for group in [report, *report["by_rank_length"]]:
        group["rule_categories"] = rule_categories(group)
    report["category_convention"] = (
        "Categories follow rule ID prefixes: control.*, layer.route.*, and catalog single.*/zip.*. "
        "Catalog-only shares exclude control actions and routing moves; a control action may itself use algebraic identities."
    )
    return report


def write_csv(report, path):
    """Write aggregate and rank/length rows with explicit denominators."""
    fields = ["scope", "n", "length", "inputs", "expected_words", "complete_words",
              "visited_fraction", "layer_complete", "semantic_term_transitions", "rule_id",
              "rule_category", "share_of_category_semantic_applications",
              "semantic_applications", "share_of_semantic_applications", "inputs_using_rule",
              "applications_per_input"]
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for group in [report, *report["by_rank_length"]]:
            for rule in group["rules"]:
                category = rule_category(rule["rule_id"])
                denominator = group["rule_categories"][category]["semantic_applications"]
                writer.writerow({
                    "scope": "rank_length" if "n" in group else "all",
                    "n": group.get("n", ""), "length": group.get("length", ""),
                    "inputs": group["inputs"],
                    **{key: group.get(key, "") for key in (
                        "expected_words", "complete_words", "visited_fraction", "layer_complete"
                    )},
                    "semantic_term_transitions": group["semantic_term_transitions"], **rule,
                    "rule_category": category,
                    "share_of_category_semantic_applications": rule["semantic_applications"] / denominator,
                })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database")
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--csv-output", type=Path)
    args = parser.parse_args()
    report = statistics(args.database)
    encoded = json.dumps(report, indent=2) + "\n"
    if args.json_output:
        args.json_output.write_text(encoded)
    else:
        print(encoded, end="")
    if args.csv_output:
        write_csv(report, args.csv_output)


if __name__ == "__main__":
    main()
