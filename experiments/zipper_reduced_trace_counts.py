"""Exact candidate-workload counts after reduced-trace canonicalization.

This is a DIFFERENT procedure: first identify inverse cancellations and distant
commutations, then run the zipper algorithm on one chosen representative. These
are not established output invariances of the existing literal-input controller.
Adjacent Artin braid relations are not imposed.

Representatives are colex normal forms: read a word right to left and choose the
lexicographically least representative of its unique reduced trace in that order.
F forbids a new movable lex inversion. C forbids a cancellable inverse letter:
F(pa)=(F(p)&I(a)) | {b<a: b in I(a)}, C(pa)=(C(p)&I(a)) | {-a}.
The current position of original bottom strand n and an undercrossing flag allow
counting the representatives which cannot peel their first singleton. These are
residual computation jobs, not a count of equivalence classes after peeling.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import inspect
import json
from pathlib import Path
import sqlite3
from time import perf_counter

from zipper_algebra.zipping_algorithm.enumeration import alphabet
from zipper_algebra.zipping_algorithm.essential_enumeration import EssentialWords
from zipper_algebra.zipping_algorithm.peeling_enumeration import NonPeelableWords


def count_layers(n: int, maximum_length: int) -> list[tuple[int, int]]:
    """Return (full-support reduced traces, nonpeelable colex reps) by length."""
    if type(n) is not int or n < 2:
        raise ValueError("n must be an integer >=2")
    if type(maximum_length) is not int or maximum_length < 0:
        raise ValueError("maximum_length must be a nonnegative integer")
    letters = tuple(range(1 - n, 0)) + tuple(range(1, n))
    index = {a: j for j, a in enumerate(letters)}
    independent = tuple(sum(1 << j for j, b in enumerate(letters)
                            if abs(abs(a) - abs(b)) >= 2) for a in letters)
    smaller = tuple(mask & ((1 << j) - 1) for j, mask in enumerate(independent))
    inverse = tuple(1 << index[-a] for a in letters)
    support = tuple(1 << (abs(a) - 1) for a in letters)
    full = (1 << (n - 1)) - 1
    states: dict[tuple[int, int, int, int, bool], int] = {(0, 0, 0, n, False): 1}
    answer = []
    for length in range(maximum_length + 1):
        total = residual = 0
        for (_, _, used, _, under), count in states.items():
            if used == full:
                total += count
                residual += count if under else 0
        answer.append((total, residual))
        if length == maximum_length:
            break
        following: defaultdict[tuple[int, int, int, int, bool], int] = defaultdict(int)
        for (forbidden, cancellation, used, position, under), count in states.items():
            for j, a in enumerate(letters):
                if (forbidden | cancellation) & (1 << j):
                    continue
                adjacent = position in (abs(a), abs(a) + 1)
                active_over = (position == abs(a)) == (a > 0)
                next_position = 2 * abs(a) + 1 - position if adjacent else position
                key = ((forbidden & independent[j]) | smaller[j],
                       (cancellation & independent[j]) | inverse[j], used | support[j],
                       next_position, under or (adjacent and not active_over))
                following[key] += count
        states = dict(following)
    return answer


def _source_identity() -> dict[str, str]:
    files = {"script_sha256": Path(__file__).resolve(),
             "essential_enumeration_sha256": Path(inspect.getfile(EssentialWords)),
             "peeling_enumeration_sha256": Path(inspect.getfile(NonPeelableWords)),
             "enumeration_sha256": Path(inspect.getfile(alphabet))}
    return {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in files.items()}


def write_counts(database: Path, *, maximum_rank: int = 8, maximum_length: int = 10) -> str:
    """Append a source-bound count run to its own research SQLite database."""
    if type(maximum_rank) is not int or maximum_rank < 2:
        raise ValueError("maximum_rank must be an integer >=2")
    if type(maximum_length) is not int or maximum_length < 0:
        raise ValueError("maximum_length must be a nonnegative integer")
    identities = _source_identity()
    definition = {
        "schema": "reduced-trace-baseline-v2", "maximum_rank": maximum_rank,
        "maximum_length": maximum_length, "order": "colex-signed-integer",
        "relations": ["distant commutation", "inverse cancellation"],
        "current_controller_output_invariance": False,
        "primitive_filter": "original bottom n has an undercrossing",
        "count_semantics": "full-support representatives and residual jobs; not full braid classes",
        "count_storage": "arbitrary precision decimal strings",
        **identities,
    }
    payload = json.dumps(definition, sort_keys=True)
    run_id = hashlib.sha256(payload.encode()).hexdigest()
    database.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE IF NOT EXISTS runs (run_id TEXT PRIMARY KEY, definition TEXT NOT NULL, created_utc TEXT NOT NULL)")
        connection.execute("CREATE TABLE IF NOT EXISTS counts (run_id TEXT NOT NULL, n INTEGER NOT NULL, length INTEGER NOT NULL, literal_full_support TEXT NOT NULL, trace_full_support TEXT NOT NULL, reduced_trace_full_support TEXT NOT NULL, reduced_trace_nonpeelable TEXT NOT NULL, PRIMARY KEY(run_id,n,length))")
        connection.execute("CREATE TABLE IF NOT EXISTS current_comparison (run_id TEXT NOT NULL, n INTEGER NOT NULL, length INTEGER NOT NULL, current_nonpeelable TEXT NOT NULL, PRIMARY KEY(run_id,n,length))")
        connection.execute("INSERT OR IGNORE INTO runs VALUES (?,?,?)", (run_id, payload, datetime.now(timezone.utc).isoformat()))
        for n in range(2, maximum_rank + 1):
            rows = count_layers(n, maximum_length)
            for length, (reduced, primitive) in enumerate(rows):
                literal = EssentialWords(n, length, commute=False).count
                traces = EssentialWords(n, length, commute=True).count
                values = (run_id, n, length, str(literal), str(traces), str(reduced), str(primitive))
                old = connection.execute("SELECT * FROM counts WHERE run_id=? AND n=? AND length=?", values[:3]).fetchone()
                if old is not None and old != values:
                    raise ValueError("existing source-bound count disagrees")
                connection.execute("INSERT OR IGNORE INTO counts VALUES (?,?,?,?,?,?,?)", values)
                comparison = (run_id, n, length, str(NonPeelableWords(n, length).count))
                previous = connection.execute("SELECT * FROM current_comparison WHERE run_id=? AND n=? AND length=?", comparison[:3]).fetchone()
                if previous is not None and previous != comparison:
                    raise ValueError("existing source-bound comparison disagrees")
                connection.execute("INSERT OR IGNORE INTO current_comparison VALUES (?,?,?,?)", comparison)
        if _source_identity() != identities:
            raise ValueError("source changed during count generation")
        integrity = connection.execute("PRAGMA integrity_check").fetchone()
        if integrity != ("ok",):
            raise ValueError(f"SQLite integrity check failed: {integrity}")
    return run_id


def export_counts(database: Path, run_id: str, destination: Path) -> None:
    """Export only the chosen source-bound run, never sum across archived runs."""
    columns = ("literal_full_support", "trace_full_support", "reduced_trace_full_support",
               "reduced_trace_nonpeelable", "current_nonpeelable")
    with sqlite3.connect(database) as connection:
        # All metadata and rows must belong to one coherent SQLite snapshot.
        connection.execute("BEGIN")
        definition_row = connection.execute("SELECT definition FROM runs WHERE run_id=?", (run_id,)).fetchone()
        if definition_row is None:
            raise ValueError("unknown run identifier")
        if hashlib.sha256(definition_row[0].encode()).hexdigest() != run_id:
            raise ValueError("run definition does not match its identifier")
        definition = json.loads(definition_row[0])
        maximum_rank, maximum_length = definition["maximum_rank"], definition["maximum_length"]
        if (type(maximum_rank) is not int or maximum_rank < 2
                or type(maximum_length) is not int or maximum_length < 0):
            raise ValueError("invalid declared count grid")
        expected = {(n, length) for n in range(2, maximum_rank + 1)
                    for length in range(maximum_length + 1)}
        count_rows = connection.execute("SELECT n,length,literal_full_support,trace_full_support,reduced_trace_full_support,reduced_trace_nonpeelable FROM counts WHERE run_id=? ORDER BY n,length", (run_id,)).fetchall()
        comparison_rows = connection.execute("SELECT n,length,current_nonpeelable FROM current_comparison WHERE run_id=?", (run_id,)).fetchall()
        if {row[:2] for row in count_rows} != expected:
            raise ValueError("counts rows do not match the complete declared grid")
        if {row[:2] for row in comparison_rows} != expected:
            raise ValueError("comparison rows do not match the complete declared grid")
        comparison = {row[:2]: row[2] for row in comparison_rows}
        rows = [(*row, comparison[row[:2]]) for row in count_rows]
    records = [dict(n=row[0], length=row[1], **dict(zip(columns, map(int, row[2:])))) for row in rows]
    totals = {column: sum(row[column] for row in records) for column in columns}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps({"run_id": run_id, "definition": json.loads(definition_row[0]),
                                      "database": str(database.resolve()), "global_identity_extra": 1,
                                      "totals_excluding_identity": totals, "rows": records}, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--maximum-rank", type=int, default=8)
    parser.add_argument("--maximum-length", type=int, default=10)
    args = parser.parse_args()
    started = perf_counter()
    run_id = write_counts(args.database, maximum_rank=args.maximum_rank,
                          maximum_length=args.maximum_length)
    export_counts(args.database, run_id, args.json_output or args.database.with_suffix(".json"))
    print(json.dumps({"run_id": run_id, "database": str(args.database),
                      "seconds": perf_counter() - started}))


if __name__ == "__main__":
    main()
