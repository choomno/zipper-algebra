"""Count a revised terminating representative policy with local Artin reductions.

Start in colex reduced-trace normal form, replace a decreasing signed Artin
triple (optionally an Artin substitution exposing inverse cancellation),
and normalize again until no replacement applies. Length/colex strictly
decreases, so this is a valid representative-producing procedure. It is NOT
asserted confluent or canonical for braid equality. Counts are its fixed input
words, optionally with initial singleton peeling excluded, not braid classes.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
from itertools import product
import json
from pathlib import Path
import sqlite3
from time import perf_counter


def decreasing_triples(n: int) -> dict[tuple[int, ...], tuple[int, ...]]:
    """Orient the six valid signed RIII variants by signed integer lex order."""
    if type(n) is not int or n < 2:
        raise ValueError("n must be an integer >=2")
    rules: dict[tuple[int, ...], tuple[int, ...]] = {}
    for i in range(1, n - 1):
        for a, b in ((i, i + 1), (i + 1, i)):
            for e, f, g in product((-1, 1), repeat=3):
                if e == g and e != f:
                    continue
                left, right = (e * a, f * b, g * a), (g * b, f * a, e * b)
                if left > right:
                    rules[left] = right
    return rules


def decreasing_rules(n: int, *, length_four: bool = False) -> dict[tuple[int, ...], tuple[int, ...]]:
    """Optionally include Artin substitutions immediately exposing cancellation."""
    if type(length_four) is not bool:
        raise ValueError("length_four must be a boolean")
    rules = decreasing_triples(n)
    if length_four:
        for left, right in tuple(rules.items()):
            for lhs, rhs in ((left, right), (right, left)):
                for source, target in ((lhs + (-rhs[-1],), rhs[:-1]),
                                       ((-rhs[0],) + lhs, rhs[1:])):
                    rules[source] = min(rules.get(source, target), target)
    return rules


def _reduced_trace(word: tuple[int, ...]) -> tuple[int, ...]:
    remaining = list(word)
    while True:
        pair = next(((i, j) for i in range(len(remaining))
                     for j in range(i + 1, len(remaining))
                     if remaining[j] == -remaining[i]
                     and all(abs(abs(x) - abs(remaining[i])) >= 2
                             for x in remaining[i + 1:j])), None)
        if pair is None:
            break
        i, j = pair
        del remaining[j]
        del remaining[i]
    result = []
    while remaining:
        a, i = min((a, i) for i, a in enumerate(remaining)
                   if all(abs(abs(a) - abs(b)) >= 2 for b in remaining[:i]))
        result.append(a)
        del remaining[i]
    return tuple(result)


def representative(word: tuple[int, ...], *, n: int, length_four: bool = False) -> tuple[int, ...]:
    """Small-input reference policy; deliberately not a production normalizer."""
    rules = decreasing_rules(n, length_four=length_four)
    if any(type(a) is not int or not 0 < abs(a) < n for a in word):
        raise ValueError("invalid generator")
    current = _reduced_trace(tuple(reversed(word)))
    while True:
        found = next(((i, width) for width in (4, 3) for i in range(len(current) - width + 1)
                      if current[i:i + width] in rules), None)
        if found is None:
            return tuple(reversed(current))
        i, width = found
        replaced = current[:i] + rules[current[i:i + width]] + current[i + width:]
        following = _reduced_trace(replaced)
        assert (len(following), following) < (len(current), current)
        current = following


def count_layers(n: int, maximum_length: int, *, length_four: bool = False) -> list[tuple[int, int]]:
    """Exact (full-support fixed words, first-peel-excluded jobs), by length."""
    rules = decreasing_rules(n, length_four=length_four)
    if type(maximum_length) is not int or maximum_length < 0:
        raise ValueError("maximum_length must be a nonnegative integer")
    letters = tuple(range(1 - n, 0)) + tuple(range(1, n))
    index = {a: j for j, a in enumerate(letters)}
    independent = tuple(sum(1 << j for j, b in enumerate(letters)
                            if abs(abs(a) - abs(b)) >= 2) for a in letters)
    smaller = tuple(mask & ((1 << j) - 1) for j, mask in enumerate(independent))
    inverse = tuple(1 << index[-a] for a in letters)
    support = tuple(1 << (abs(a) - 1) for a in letters)
    prefixes = [()] + sorted({r[:k] for r in rules for k in range(1, len(r))})
    longest_prefix = max(map(len, prefixes))
    prefix_index = {p: i for i, p in enumerate(prefixes)}
    transitions = []
    for prefix in prefixes:
        row = []
        for a in letters:
            candidate = prefix + (a,)
            if any(candidate[-k:] in rules for k in (3, 4) if k <= len(candidate)):
                row.append(-1)
            else:
                matched_suffix = next((candidate[-k:] for k in range(min(longest_prefix, len(candidate)), 0, -1)
                                       if candidate[-k:] in prefix_index), ())
                row.append(prefix_index[matched_suffix])
        transitions.append(row)
    full, all_letters = (1 << (n - 1)) - 1, (1 << len(letters)) - 1
    # Position 0 means an undercrossing has occurred; its later position is immaterial.
    states = {(0, 0, 0, n, 0): 1}
    answer = []
    for length in range(maximum_length + 1):
        total = residual = 0
        for (_, _, used, position, _), count in states.items():
            if used == full:
                total += count
                residual += count if position == 0 else 0
        answer.append((total, residual))
        if length == maximum_length:
            break
        following: defaultdict[tuple[int, int, int, int, int], int] = defaultdict(int)
        for (forbidden, cancellation, used, position, suffix), count in states.items():
            allowed = all_letters & ~(forbidden | cancellation)
            while allowed:
                bit = allowed & -allowed
                allowed -= bit
                j = bit.bit_length() - 1
                next_suffix = transitions[suffix][j]
                if next_suffix < 0:
                    continue
                next_used = used | support[j]
                if (full ^ next_used).bit_count() > maximum_length - length - 1:
                    continue
                a = letters[j]
                next_position = position
                if position in (abs(a), abs(a) + 1):
                    next_position = 2 * abs(a) + 1 - position if (position == abs(a)) == (a > 0) else 0
                key = ((forbidden & independent[j]) | smaller[j],
                       (cancellation & independent[j]) | inverse[j], next_used,
                       next_position, next_suffix)
                following[key] += count
        states = dict(following)
    return answer


def write_counts(database: Path, *, maximum_rank: int = 8, maximum_length: int = 10,
                 length_four: bool = False) -> str:
    if type(maximum_rank) is not int or maximum_rank < 2:
        raise ValueError("maximum_rank must be an integer >=2")
    if type(maximum_length) is not int or maximum_length < 0:
        raise ValueError("maximum_length must be a nonnegative integer")
    if type(length_four) is not bool:
        raise ValueError("length_four must be a boolean")
    source = Path(__file__).read_bytes()
    definition = {"schema": "local-artin-baseline-v2", "length_four": length_four, "maximum_rank": maximum_rank,
                  "maximum_length": maximum_length, "order": "colex-signed-integer",
                  "policy": "reduced trace; optional length4 cancellation; decreasing signed Artin triple; repeat",
                  "braid_class_canonical": False, "current_output_invariance": False,
                  "count_semantics": "fixed words and first-peel-excluded recursive jobs",
                  "source_sha256": hashlib.sha256(source).hexdigest()}
    payload = json.dumps(definition, sort_keys=True)
    run_id = hashlib.sha256(payload.encode()).hexdigest()
    database.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(database) as c:
        c.execute("CREATE TABLE IF NOT EXISTS runs(run_id TEXT PRIMARY KEY, definition TEXT NOT NULL, created_utc TEXT NOT NULL)")
        c.execute("CREATE TABLE IF NOT EXISTS counts(run_id TEXT,n INTEGER,length INTEGER,normal_words TEXT,residual_jobs TEXT,PRIMARY KEY(run_id,n,length))")
        c.execute("INSERT OR IGNORE INTO runs VALUES(?,?,?)", (run_id, payload, datetime.now(timezone.utc).isoformat()))
        for n in range(2, maximum_rank + 1):
            start = perf_counter()
            rows = count_layers(n, maximum_length, length_four=length_four)
            for length, (normal, residual) in enumerate(rows):
                record = (run_id, n, length, str(normal), str(residual))
                old = c.execute("SELECT * FROM counts WHERE run_id=? AND n=? AND length=?", record[:3]).fetchone()
                if old is not None and old != record:
                    raise ValueError("existing source-bound count disagrees")
                c.execute("INSERT OR IGNORE INTO counts VALUES(?,?,?,?,?)", record)
            print(json.dumps({"n": n, "seconds": perf_counter() - start, "normal_total": sum(x for x, _ in rows), "residual_total": sum(y for _, y in rows), "length_max": rows[-1]}), flush=True)
        if Path(__file__).read_bytes() != source:
            raise ValueError("source changed during counting")
        if c.execute("PRAGMA integrity_check").fetchone() != ("ok",):
            raise ValueError("SQLite integrity check failed")
    return run_id


def export_counts(database: Path, run_id: str, destination: Path) -> None:
    """Export exactly one complete source-bound run."""
    with sqlite3.connect(database) as c:
        c.execute("BEGIN")
        record = c.execute("SELECT definition FROM runs WHERE run_id=?", (run_id,)).fetchone()
        if record is None:
            raise ValueError("unknown run identifier")
        definition = json.loads(record[0])
        if hashlib.sha256(json.dumps(definition, sort_keys=True).encode()).hexdigest() != run_id:
            raise ValueError("run definition hash mismatch")
        rows = c.execute("SELECT n,length,normal_words,residual_jobs FROM counts WHERE run_id=? ORDER BY n,length", (run_id,)).fetchall()
    expected = [(n, length) for n in range(2, definition["maximum_rank"] + 1)
                for length in range(definition["maximum_length"] + 1)]
    if [(r[0], r[1]) for r in rows] != expected:
        raise ValueError("incomplete count grid")
    result = {"run_id": run_id, "definition": definition,
              "global_identity_extra": 1,
              "normal_words": sum(int(r[2]) for r in rows),
              "residual_jobs": sum(int(r[3]) for r in rows),
              "rows": [dict(n=n, length=length, normal_words=int(normal), residual_jobs=int(residual))
                       for n, length, normal, residual in rows]}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--length-four", action="store_true")
    parser.add_argument("--maximum-rank", type=int, default=8)
    parser.add_argument("--maximum-length", type=int, default=10)
    args = parser.parse_args()
    run_id = write_counts(args.database, maximum_rank=args.maximum_rank, maximum_length=args.maximum_length, length_four=args.length_four)
    export_counts(args.database, run_id, args.json_output or args.database.with_suffix(".json"))
    print(run_id)


if __name__ == "__main__":
    main()
