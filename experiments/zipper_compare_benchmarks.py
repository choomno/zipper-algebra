"""Compare exact rich output records and continuation states across executors."""

import argparse, json, hashlib, sqlite3, zlib
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("baseline")
p.add_argument("candidate")
p.add_argument("--output", required=True)
p.add_argument("--database")
a = p.parse_args()
b = Path(a.baseline)
c = Path(a.candidate)
rows = []
database = (
    Path(a.database).resolve()
    if a.database
    else c.resolve().parent / "benchmark-results.sqlite3"
)
root = database.parent


def raw(path):
    if path.exists():
        return path.read_bytes()
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as db:
        payload = db.execute(
            "SELECT payload FROM certificates WHERE alias=?",
            (str(path.resolve().relative_to(root)),),
        ).fetchone()
        if payload is None:
            raise FileNotFoundError(path)
        return zlib.decompress(payload[0])


paths = sorted(c.glob("pass*-case*.json"))
if not paths:
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as db:
        paths = [
            root / r[0]
            for r in db.execute(
                "SELECT alias FROM certificates WHERE run_alias=? ORDER BY pass,case_index",
                (str(c.resolve().relative_to(root)),),
            )
        ]


def canonical(records):
    return sorted(json.dumps(x, sort_keys=True, separators=(",", ":")) for x in records)


for candidate in paths:
    idx = candidate.name.split("-case")[1]
    baseline = b / ("pass0-case" + idx)
    x = json.loads(raw(baseline))
    y = json.loads(raw(candidate))
    checks = {k: x[k] == y[k] for k in ["n", "input_word", "status", "metrics"]}
    checks["output_records"] = canonical(x["output"]) == canonical(y["output"])
    checks["pending_records"] = canonical(x["pending"]) == canonical(y["pending"])
    rows.append(
        {
            "baseline": str(baseline),
            "candidate": str(candidate),
            "checks": checks,
            "baseline_sha256": hashlib.sha256(raw(baseline)).hexdigest(),
            "candidate_sha256": hashlib.sha256(raw(candidate)).hexdigest(),
        }
    )
result = {
    "scope": "Exact JSON rich-record comparison; independent baseline producer, not an algebraic quotient equivalence test or candidate replay.",
    "cases": rows,
    "all_pass": bool(rows) and all(all(r["checks"].values()) for r in rows),
}
Path(a.output).write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"comparisons": len(rows), "all_pass": result["all_pass"]}))
raise SystemExit(0 if result["all_pass"] else 1)
