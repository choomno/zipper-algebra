"""Capture actual admitted catalog rewrites and benchmark packed GPU expansion.

Only the local catalog expansion stage is offloaded. Scheduler admission,
geometric transport, collection and recursive control stay on the CPU.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
from copy import deepcopy
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import platform
import sqlite3
import sys
import time
from typing import Any
import zlib

from zipper_algebra.zipping_algorithm import algorithm, catalog, collected, framed, interpreter

SCHEMA = "zipper-gpu-rewrite-capture-v1"
DEFAULT_WORKLOADS = [
    {"name": "alternating-B3-L10", "n": 3, "word": [1, -2] * 5},
    {"name": "positive-B5-L6", "n": 5, "word": [1, 2, 3, 4, 1, 2]},
    {"name": "positive-B8-L7", "n": 8, "word": list(range(1, 8))},
]


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def digest(value: Any) -> str:
    return hashlib.sha256(canonical(value)).hexdigest()


def encode(value: Any) -> tuple[bytes, str]:
    raw = canonical(value)
    return zlib.compress(raw), hashlib.sha256(raw).hexdigest()


def decode(payload: bytes, expected: str) -> Any:
    raw = zlib.decompress(payload)
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("compressed record hash mismatch")
    return json.loads(raw)


class RuleView:
    """Exact captured catalog record, without derived-rule regeneration."""

    def __init__(self, rule):
        self.rule = rule

    def get(self, rule_id):
        if rule_id != self.rule["id"]:
            raise ValueError("captured rule ID mismatch")
        return self.rule


def reference_expand(cases: list[dict[str, Any]]) -> tuple[dict[str, Any], ...]:
    """Independent contextual replay with arbitrary-precision Python integers."""
    result = []
    for case in cases:
        predecessor = interpreter.Branch.from_json(case["predecessor"])
        request = interpreter.Request(
            0, case["offset"], case["rule"]["id"], case["parameters"],
            words=case["words"], successor_controls=case["successor_controls"],
        )
        following, _ = interpreter.apply_rule(
            interpreter.RunState(case["n"], (predecessor,)), request,
            catalog=RuleView(case["rule"]), limits=interpreter.Limits(**case["limits"]),
            trace=False,
        )
        result.append({"root": case["root"], "id": case["id"],
                       "branches": [branch.to_json() for branch in following.branches]})
    return tuple(result)


def single_projection(instance, saved):
    """Adapt the source's serialized single-rule ordering to contextual words.

    This is a conceptual full-exterior identity, not the single scheduler
    operation. Serialized catalog order and apply_rule order may differ.
    """
    rhs = {}
    for term in instance["rhs"]:
        rhs.setdefault(canonical(term["word"]), Counter())[term["q_power"]] += term["integer"]
    predecessor = saved["predecessor"]
    offset = saved["offset"]
    # Single-rule LHS is one crossing followed by its one-strand coupon.
    prefix, suffix = predecessor["word"][:offset], predecessor["word"][offset + 2:]
    return [
        {"word": prefix + json.loads(word_key) + suffix,
         "coefficient": [list(x) for x in interpreter._multiply(
             tuple(tuple(x) for x in saved["selected_work_coefficient"]),
             tuple((q, c) for q, c in sorted(coefficient.items()) if c))],
         "control": {}}
        for word_key, coefficient in rhs.items()
    ]


def word_control_permutation(actual, projected):
    """Return projected-to-source indices after exact branch-multiset equality."""
    locations = defaultdict(deque)
    for index, branch in enumerate(actual):
        locations[canonical(branch)].append(index)
    permutation = []
    for branch in projected:
        key = canonical(branch)
        if not locations[key]:
            raise ValueError("single projection differs as a word/coefficient/control multiset")
        permutation.append(locations[key].popleft())
    if any(locations.values()):
        raise ValueError("single projection has unmatched source branches")
    return permutation


class RewriteCaptureSession(collected.CollectedSession):
    """Observe local calls in this process without changing scheduler sources."""

    def __init__(self, root: int):
        super().__init__()
        self.root = root
        self.batches: list[dict[str, Any]] = []
        self.templates: dict[str, dict[str, Any]] = {}
        self.frontier: dict[str, Any] = {}
        self.output: dict[str, Any] = {}
        self.before: tuple[dict[str, Any], dict[str, Any]] = ({}, {})
        self.cases: list[dict[str, Any]] = []
        self.next_case = 0
        self.weight: list[list[int]] = [[0, 1]]
        self.limits = interpreter.Limits()
        self.term = 0
        self.selected_key = ""

    def _record(self, *, rule, parameters, words, predecessor, offset,
                controls, actual, kind, used, unit_limits):
        weighted = deepcopy(predecessor)
        weighted["coefficient"] = [list(pair) for pair in interpreter._multiply(
            tuple(tuple(x) for x in predecessor["coefficient"]),
            tuple(tuple(x) for x in self.weight))]
        monomials = max(1, sum(len(term["coefficient"]) for term in rule["rhs"]))
        pairs = max(1, len(weighted["coefficient"]) * monomials)
        projection_limits = asdict(unit_limits)
        projection_limits.update(
            max_terms=max(unit_limits.max_terms, len(rule["rhs"])),
            max_coefficient_support=max(unit_limits.max_coefficient_support, pairs),
            max_integer_bits=max(
                unit_limits.max_integer_bits,
                max(abs(c).bit_length() for _, c in weighted["coefficient"])
                + max(abs(c["integer"]).bit_length() for term in rule["rhs"] for c in term["coefficient"])
                + (pairs - 1).bit_length() + 2),
        )
        case = {"root": self.root, "id": self.next_case, "batch": len(self.batches),
                "n": parameters["n"], "rule": rule, "parameters": dict(parameters),
                "words": {name: list(word) for name, word in words.items()},
                "offset": offset, "predecessor": weighted,
                "successor_controls": controls, "limits": projection_limits}
        expected = reference_expand([case])[0]
        if actual is not None:
            scaled = []
            for branch in actual:
                row = branch.to_json()
                row["coefficient"] = [list(pair) for pair in interpreter._multiply(
                    branch.coefficient, tuple(tuple(x) for x in self.weight))]
                scaled.append(row)
            if scaled != expected["branches"]:
                raise ValueError("weighted captured expansion disagrees with actual local result")
        rule_sha = digest(rule)
        self.templates.setdefault(rule_sha, deepcopy(rule))
        saved = dict(case)
        del saved["rule"]
        saved.update(rule_sha256=rule_sha, expected=expected, capture_kind=kind,
                     used_by_scheduler=used, selected_term=self.term,
                     selected_work_coefficient=deepcopy(self.weight),
                     unit_predecessor_coefficient=deepcopy(predecessor["coefficient"]),
                     actual_unit_execution_limits=asdict(unit_limits),
                     projection_limits_policy="conservative bounds for weighted contextual projection; original atomic scheduler limits remain unchanged",
                     projection_semantics=(
                         "actual contextual call with original selected coefficient restored"
                         if kind == "contextual_apply_rule" else
                         "conceptual weighted full-exterior single identity; exact branch multiset with source-order permutation"),
                     selected_state_key_sha256=hashlib.sha256(self.selected_key.encode()).hexdigest())
        self.cases.append(saved)
        self.next_case += 1

    def _profile(self, frame, event, returned):
        if event != "return":
            return
        values = frame.f_locals
        if frame.f_code is interpreter.apply_rule.__code__:
            if not (isinstance(returned, tuple) and len(returned) == 2
                    and isinstance(returned[0], interpreter.RunState)):
                return
            request, state = values["request"], values["state"]
            count = len(returned[0].branches) - len(state.branches) + 1
            self._record(
                rule=deepcopy(values["rule"]), parameters=dict(request.parameters),
                words=dict(request.words), predecessor=values["predecessor"].to_json(),
                offset=request.offset, controls=request.to_json()["successor_controls"],
                actual=returned[0].branches[request.branch:request.branch + count],
                kind="contextual_apply_rule", used=True, unit_limits=values["limits"],
            )
        elif frame.f_code is catalog.instantiate_rule.__code__ and isinstance(returned, dict):
            caller = frame.f_back
            caller = caller.f_back if caller is not None else None
            if caller is None or caller.f_code is not algorithm._single_step.__code__:
                return
            work = caller.f_locals["work"]
            predecessor = interpreter.Branch(
                tuple(work["left"]) + (("slider", 1, work["position"], False, False, False),)
                + tuple(work["processed"]),
                tuple(tuple(x) for x in work["coefficient"]),
            )
            self._record(
                rule=deepcopy(values["rule"]), parameters=dict(values["env"]), words={},
                predecessor=predecessor.to_json(), offset=len(work["left"]) - 1,
                controls=None, actual=None, kind="single_catalog_instantiation",
                used=bool(caller.f_locals["offensive"]), unit_limits=caller.f_locals["limits"],
            )
            saved = self.cases[-1]
            saved["single_projected_to_source_order"] = word_control_permutation(
                single_projection(returned, saved), saved["expected"]["branches"])
            self.cases[-1]["actual_single_instance_sha256"] = digest(returned)

    def _execute_batch(self, selected, limits, trace, workers, batch_size, pool):
        if workers != 1:
            raise ValueError("rewrite capture requires one term worker")
        if sys.getprofile() is not None:
            raise ValueError("rewrite capture cannot replace an existing profiler")
        if not self.batches:
            self.frontier = {key: deepcopy(item) for key, item in selected}
        self.before = (dict(self.frontier), dict(self.output))
        self.cases = []
        self.limits = limits
        results = []
        for self.term, (self.selected_key, item) in enumerate(selected):
            self.weight = item["coefficient"]
            sys.setprofile(self._profile)
            try:
                results.append(collected._transition(item, limits, trace, self._catalog))
            finally:
                sys.setprofile(None)
        rules = Counter()
        for result in results:
            if "error" not in result:
                event = result["unit_evidence"]
                rules[str(event.get("rule_id", event.get("catalog_transition", {}).get(
                    "request", {}).get("rule_id")))] += 1
        self.batches.append({"batch": len(self.batches), "selected_steps": len(selected),
                             "cases": self.cases, "rule_applications": dict(rules)})
        if not any("error" in result for result in results):
            removed = {key for key, _ in selected}
            for pending, current, incoming, excluded in (
                (True, self.frontier, [x for r in results for x in r["children"]], removed),
                (False, self.output, [x for r in results for x in r["output"]], set()),
            ):
                updates, _, _ = collected._stage(current, incoming, pending=pending, removed=excluded)
                collected._commit(current, updates, excluded)
        return results, pool

    def reconcile(self, certificate):
        committed = certificate["metrics"]["batches"]
        if len(self.batches) == committed + 1:
            self.frontier, self.output = self.before
        elif len(self.batches) != committed:
            raise ValueError("unexpected captured/committed batch count")
        if list(self.frontier.values()) != certificate["pending"] or list(self.output.values()) != certificate["output"]:
            raise ValueError("final capture shadow differs from scheduler certificate")
        self.batches = self.batches[:committed]
        counted = Counter(self.templates[c["rule_sha256"]]["id"]
                          for batch in self.batches for c in batch["cases"])
        expected_counts = {name: count for name, count in
                           certificate["metrics"]["rule_applications"].items()
                           if name.startswith(("single.", "zip."))}
        if counted != expected_counts:
            raise ValueError("captured catalog calls differ from scheduler rule histogram")


def initialize_database(db):
    db.executescript("""
        CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
        CREATE TABLE templates(sha TEXT PRIMARY KEY,payload BLOB NOT NULL);
        CREATE TABLE workloads(root INTEGER PRIMARY KEY,name TEXT NOT NULL,n INTEGER NOT NULL,
            word TEXT NOT NULL,status TEXT NOT NULL,certificate BLOB NOT NULL,
            certificate_sha TEXT NOT NULL,metadata TEXT NOT NULL);
        CREATE TABLE batches(root INTEGER NOT NULL,batch INTEGER NOT NULL,payload BLOB NOT NULL,
            payload_sha TEXT NOT NULL,cases INTEGER NOT NULL,selected_steps INTEGER NOT NULL,
            PRIMARY KEY(root,batch),FOREIGN KEY(root) REFERENCES workloads(root));
        CREATE TABLE benchmarks(id INTEGER PRIMARY KEY,created_utc TEXT DEFAULT CURRENT_TIMESTAMP,
            backend TEXT NOT NULL,result TEXT NOT NULL);
    """)
    db.execute("INSERT INTO metadata VALUES('schema',?)", (json.dumps(SCHEMA),))


def capture(database: Path, workloads: list[dict[str, Any]], *,
            max_steps=20000, batch_size=1024, max_terms=10000):
    database = database.resolve()
    if database.exists() and database.stat().st_size:
        raise ValueError("capture requires a new database")
    if not workloads:
        raise ValueError("at least one workload is required")
    database.parent.mkdir(parents=True, exist_ok=True)
    limits = interpreter.Limits(max_steps=max_steps, max_terms=max_terms,
                                expansion_budget=1000000)
    source_sha = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    manifest = []
    with sqlite3.connect(database) as db:
        initialize_database(db)
        db.execute("INSERT INTO metadata VALUES('intended_workloads',?)",
                   (json.dumps(workloads),))
        db.execute("INSERT INTO metadata VALUES('complete', 'false')")
        db.commit()
        for root, work in enumerate(workloads):
            session = RewriteCaptureSession(root)
            started = time.perf_counter()
            cert = session.reduce(work["word"], n=work["n"], limits=limits,
                                  trace=False, workers=1, batch_size=batch_size)
            capture_wall = time.perf_counter() - started
            cert = json.loads(json.dumps(cert))
            session.reconcile(cert)
            started = time.perf_counter()
            collected.CollectedSession().check(cert)
            replay_wall = time.perf_counter() - started
            payload, sha = encode(cert)
            metadata = {"capture_source_sha256": source_sha,
                        "source_sha256": cert["source_sha256"],
                        "catalog_sha256": cert["catalog_sha256"], "limits": asdict(limits),
                        "batch_size": batch_size, "capture_wall": capture_wall,
                        "full_certificate_replay_wall": replay_wall,
                        "final_shadow_reconciled": True, "certificate_reloaded_and_replayed": True,
                        "coefficient_restoration": "actual selected work weight multiplied by unit-local predecessor coefficient",
                        "metrics": cert["metrics"]}
            db.execute("INSERT INTO workloads VALUES(?,?,?,?,?,?,?,?)",
                       (root, work.get("name", str(root)), work["n"], json.dumps(work["word"]),
                        cert["status"], payload, sha, json.dumps(metadata)))
            for rule_sha, rule in session.templates.items():
                compressed, encoded_sha = encode(rule)
                if rule_sha != encoded_sha:
                    raise ValueError("template hash mismatch")
                db.execute("INSERT OR IGNORE INTO templates VALUES(?,?)", (rule_sha, compressed))
            case_count = 0
            for batch in session.batches:
                compressed, hashed = encode(batch)
                case_count += len(batch["cases"])
                db.execute("INSERT INTO batches VALUES(?,?,?,?,?,?)", (
                    root, batch["batch"], compressed, hashed, len(batch["cases"]),
                    batch["selected_steps"]))
            row = {"root": root, "n": work["n"], "name": work.get("name", str(root)),
                   "word": work["word"], "status": cert["status"],
                   "committed_batches": len(session.batches), "rewrite_cases": case_count,
                   **metadata}
            manifest.append(row)
            db.commit()
            print(json.dumps({k: v for k, v in row.items()
                              if k not in {"metrics", "source_sha256"}}), flush=True)
        db.execute("UPDATE metadata SET value='true' WHERE key='complete'")
    return {"database": str(database), "workloads": manifest}


def capture_tables_identity(db: sqlite3.Connection) -> str:
    """Commit current capture tables within the caller's SQLite transaction.

    This is a prepared-snapshot identity, not an original capture-time signature.
    Benchmark append rows are intentionally excluded.
    """
    result = hashlib.sha256()

    def field(value):
        if value is None:
            encoded, tag = b"", b"null"
        elif isinstance(value, bytes):
            encoded, tag = value, b"bytes"
        elif isinstance(value, str):
            encoded, tag = value.encode(), b"text"
        elif type(value) is int:
            encoded, tag = str(value).encode(), b"integer"
        else:
            raise ValueError("unsupported SQLite value in capture identity")
        result.update(len(tag).to_bytes(8, "big"))
        result.update(tag)
        result.update(len(encoded).to_bytes(8, "big"))
        result.update(encoded)

    for table, ordering in (
        ("metadata", "key"), ("workloads", "root"),
        ("templates", "sha"), ("batches", "root,batch"),
    ):
        schema = db.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name=?",
                            (table,)).fetchone()
        if schema is None:
            raise ValueError("capture table is missing during identity validation")
        field(table)
        field(schema[0])
        for row in db.execute(f"SELECT * FROM {table} ORDER BY {ordering}"):
            field("row")
            field(len(row))
            for value in row:
                field(value)
        field("end")
    return result.hexdigest()


def load_cases(database: Path, roots: set[int] | None = None):
    chunks = []
    metadata_rows = []
    expected_rows = []
    bundled = catalog.Catalog.bundled()
    source = collected._source()
    checked_bindings = set()
    with sqlite3.connect(database.resolve().as_uri() + "?mode=ro", uri=True) as db:
        db.execute("BEGIN")
        meta = {k: json.loads(v) for k, v in db.execute("SELECT key,value FROM metadata")}
        if meta.get("schema") != SCHEMA or meta.get("complete") is not True:
            raise ValueError("unsupported or incomplete rewrite capture")
        intended = meta["intended_workloads"]
        workload_rows = db.execute(
            "SELECT root,n,word,status,certificate,certificate_sha,metadata FROM workloads ORDER BY root"
        ).fetchall()
        if [row[0] for row in workload_rows] != list(range(len(intended))):
            raise ValueError("missing or orphan workload root")
        available = {row[0] for row in workload_rows}
        if roots is not None and (not roots or not roots <= available):
            raise ValueError("unknown or empty requested roots")
        if {row[0] for row in db.execute("SELECT DISTINCT root FROM batches")} - available:
            raise ValueError("orphan batch root")
        templates = {sha: decode(payload, sha) for sha, payload in db.execute("SELECT sha,payload FROM templates")}
        for root, n, word, status, cert_payload, cert_sha, saved_meta in workload_rows:
            if roots is not None and root not in roots:
                continue
            certificate = decode(cert_payload, cert_sha)
            saved_metadata = json.loads(saved_meta)
            work = intended[root]
            if (n != work["n"] or json.loads(word) != work["word"]
                    or certificate["n"] != n or certificate["input_word"] != work["word"]
                    or certificate["status"] != status):
                raise ValueError("workload/certificate/manifest mismatch")
            if (certificate["source_sha256"] != source
                    or certificate["catalog_sha256"] != bundled.sha256
                    or saved_metadata["source_sha256"] != certificate["source_sha256"]
                    or saved_metadata["catalog_sha256"] != certificate["catalog_sha256"]
                    or saved_metadata["limits"] != certificate["limits"]
                    or saved_metadata["batch_size"] != certificate["execution"]["batch_size"]
                    or saved_metadata["final_shadow_reconciled"] is not True
                    or saved_metadata["certificate_reloaded_and_replayed"] is not True):
                raise ValueError("capture source/catalog/certificate metadata mismatch")
            batches = db.execute(
                "SELECT batch,payload,payload_sha,cases,selected_steps FROM batches WHERE root=? ORDER BY batch",
                (root,)).fetchall()
            if [row[0] for row in batches] != list(range(certificate["metrics"]["batches"])):
                raise ValueError("missing or noncontiguous committed batch")
            ids = []
            steps = 0
            case_count = 0
            counted_rules = Counter()
            for index, payload, sha, count, selected in batches:
                batch = decode(payload, sha)
                if (batch["batch"] != index or len(batch["cases"]) != count
                        or batch["selected_steps"] != selected):
                    raise ValueError("batch count/identity mismatch")
                cases, expected, single_orders = [], [], {}
                for saved in batch["cases"]:
                    if (saved["root"] != root or saved["batch"] != index
                            or not 1 <= saved["n"] <= n
                            or saved["parameters"].get("n") != saved["n"]
                            or saved["rule_sha256"] not in templates
                            or not 0 <= saved["selected_term"] < selected):
                        raise ValueError("case binding mismatch")
                    if saved["predecessor"]["coefficient"] != [
                        list(pair) for pair in interpreter._multiply(
                            tuple(tuple(x) for x in saved["unit_predecessor_coefficient"]),
                            tuple(tuple(x) for x in saved["selected_work_coefficient"]))]:
                        raise ValueError("selected coefficient restoration mismatch")
                    case = {key: saved[key] for key in (
                        "root", "id", "batch", "n", "parameters", "words", "offset",
                        "predecessor", "successor_controls", "limits")}
                    case["rule"] = templates[saved["rule_sha256"]]
                    binding = (saved["rule_sha256"], canonical(case["parameters"]))
                    if binding not in checked_bindings:
                        accepted = framed.rule_for_event(
                            bundled, {"rule_id": case["rule"]["id"],
                                      "parameters": case["parameters"]})
                        if canonical(accepted) != canonical(case["rule"]):
                            raise ValueError("captured template differs from bound accepted catalog rule")
                        checked_bindings.add(binding)
                    if saved["capture_kind"] == "single_catalog_instantiation":
                        instance = catalog.instantiate_rule(
                            case["rule"], case["parameters"], words=case["words"],
                            max_terms=case["limits"]["expansion_budget"], require_verified=True)
                        if digest(instance) != saved["actual_single_instance_sha256"]:
                            raise ValueError("single source instantiation hash mismatch")
                        permutation = word_control_permutation(
                            single_projection(instance, saved), saved["expected"]["branches"])
                        if ("single_projected_to_source_order" in saved
                                and saved["single_projected_to_source_order"] != permutation):
                            raise ValueError("single source-order adapter mismatch")
                        single_orders[case["id"]] = permutation
                    cases.append(case)
                    expected.append(saved["expected"])
                    ids.append(case["id"])
                    counted_rules[case["rule"]["id"]] += 1
                if reference_expand(cases) != tuple(expected):
                    raise ValueError("saved expansion differs from exact contextual replay")
                chunks.append({"root": root, "batch": index, "cases": cases,
                               "expected": expected, "selected_steps": selected,
                               "kinds": [c["capture_kind"] for c in batch["cases"]],
                               "used": [c["used_by_scheduler"] for c in batch["cases"]],
                               "single_projected_to_source_orders": single_orders,
                               "used_cases": sum(bool(c["used_by_scheduler"]) for c in batch["cases"])})
                expected_rows.extend(expected)
                steps += selected
                case_count += count
            if ids != list(range(len(ids))) or steps != certificate["metrics"]["steps"]:
                raise ValueError("case sequence or selected-step denominator mismatch")
            if counted_rules != {name: count for name, count in
                                 certificate["metrics"]["rule_applications"].items()
                                 if name.startswith(("single.", "zip."))}:
                raise ValueError("captured rule counts differ from certificate")
            metadata_rows.append({"root": root, "n": n, "word": work["word"], "status": status,
                                  "certificate_sha256": cert_sha, "rewrite_cases": case_count,
                                  "selected_steps": steps, "metadata": json.loads(saved_meta)})
        capture_identity = capture_tables_identity(db)
    if not any(chunk["cases"] for chunk in chunks):
        raise ValueError("capture has no catalog expansion cases")
    return chunks, {"workloads": metadata_rows, "chunks": len(chunks),
                    "prepared_capture_tables_sha256": capture_identity,
                    "capture_identity_binding": "capture-table content observed in the preparation read transaction; not a capture-time signature",
                    "rewrite_cases": sum(len(c["cases"]) for c in chunks),
                    "used_cases": sum(c["used_cases"] for c in chunks),
                    "capture_kind_counts": dict(Counter(
                        kind for chunk in chunks for kind in chunk["kinds"])),
                    "single_order_adapters_revalidated": sum(
                        len(chunk["single_projected_to_source_orders"]) for chunk in chunks),
                    "nonidentity_single_order_permutations": sum(
                        order != list(range(len(order)))
                        for chunk in chunks for order in chunk["single_projected_to_source_orders"].values()),
                    "single_projection_validation": "weighted word/coefficient/control multiset and explicit projected-to-source permutation; conceptual operation beyond literal single controller",
                    "selected_steps": sum(c["selected_steps"] for c in chunks),
                    "input_sha256": digest([c["cases"] for c in chunks]),
                    "expected_sha256": digest(expected_rows)}


def expansion_sources():
    directory = Path(collected.__file__).parent
    return {**collected._source(),
            **{name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
               for name in ("gpu_rewrite.py", "gpu_rewrite_program.py",
                            "gpu_rewrite_host.py", "gpu_collection.py")},
            "zipper_gpu_rewrite_benchmark.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            **{name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
               for name in ("zipper_native_rewrite.py",)},
            "native/zipper_rewrite_cpu.cpp": hashlib.sha256(
                (Path(__file__).parent / "native/zipper_rewrite_cpu.cpp").read_bytes()).hexdigest()}


@dataclass(frozen=True)
class PreparedSelection:
    database: Path
    roots: frozenset[int] | None
    mode: str
    kind: str
    chunks: list[list[dict[str, Any]]]
    workload: dict[str, Any]
    preparation_wall: float
    source_sha256: dict[str, str]
    payload_sha256: str


def prepare_benchmark(database: Path, *, roots=None, mode="sequential", kind="contextual"):
    """Verify a source-bound immutable input snapshot once for paired trials.

    Every timed trial still runs the full external exact expansion verifier.
    The prepared snapshot hash is rechecked before use; SQLite benchmark appends
    do not invalidate its explicitly recorded workload/input identities.
    """
    if mode not in {"sequential", "pooled"} or kind not in {"contextual", "used", "all"}:
        raise ValueError("invalid benchmark selection")
    loaded = time.perf_counter()
    chunks, workload = load_cases(database, roots)
    selected_chunks, selected_expected = [], []
    for chunk in chunks:
        selected = [(case, expected) for case, expected, captured_kind, used in zip(
            chunk["cases"], chunk["expected"], chunk["kinds"], chunk["used"], strict=True)
                    if kind == "all" or (kind == "used" and used)
                    or (kind == "contextual" and captured_kind == "contextual_apply_rule")]
        if selected:
            selected_chunks.append([case for case, _ in selected])
            selected_expected.extend(expected for _, expected in selected)
    if not selected_chunks:
        raise ValueError("no rewrite cases selected by kind")
    workload.update(captured_all_rewrite_cases=workload["rewrite_cases"],
                    selected_rewrite_cases=sum(map(len, selected_chunks)),
                    expected_sha256=digest(selected_expected),
                    selected_input_sha256=digest(selected_chunks), selection_kind=kind)
    workload["selected_rule_counts"] = dict(Counter(
        case["rule"]["id"] for chunk in selected_chunks for case in chunk))
    workload["selected_local_ranks"] = dict(Counter(
        case["n"] for chunk in selected_chunks for case in chunk))
    chunks = selected_chunks
    if mode == "pooled":
        chunks = [[case for chunk in chunks for case in chunk]]
    return PreparedSelection(
        database.resolve(), None if roots is None else frozenset(roots), mode, kind,
        chunks, workload, time.perf_counter() - loaded, expansion_sources(),
        digest({"chunks": chunks, "workload": workload}),
    )


def benchmark(database: Path, *, backend: str, repeats=3, roots=None,
              mode="sequential", kind="contextual", expander=None, prepared=None):
    if (backend not in {"cpu", "opencl", "cpp"} or mode not in {"sequential", "pooled"}
            or kind not in {"contextual", "used", "all"} or repeats < 1):
        raise ValueError("invalid benchmark configuration")
    from zipper_algebra.zipping_algorithm.gpu_rewrite import OpenClRewriteExpander, expand_cpu

    reused_preparation = prepared is not None
    if prepared is None:
        prepared = prepare_benchmark(database, roots=roots, mode=mode, kind=kind)
    revalidation_started = time.perf_counter()
    requested_roots = None if roots is None else frozenset(roots)
    if (not isinstance(prepared, PreparedSelection)
            or prepared.database != database.resolve()
            or prepared.roots != requested_roots or prepared.mode != mode or prepared.kind != kind):
        raise ValueError("prepared selection configuration differs")
    if (prepared.source_sha256 != expansion_sources()
            or digest({"chunks": prepared.chunks, "workload": prepared.workload}) != prepared.payload_sha256):
        raise ValueError("prepared selection source/input binding differs")
    revalidation_wall = time.perf_counter() - revalidation_started
    chunks, workload = prepared.chunks, prepared.workload
    if expander is not None and (backend == "cpu" or expander.metadata["backend"] != backend):
        raise ValueError("reused expander backend differs")
    retained = expander is not None
    setup = time.perf_counter()
    if backend == "opencl" and expander is None:
        expander = OpenClRewriteExpander()
    elif backend == "cpp" and expander is None:
        from experiments.zipper_native_rewrite import NativeRewriteExpander
        expander = NativeRewriteExpander()
    setup_wall = 0.0 if retained else time.perf_counter() - setup
    rows = []
    for iteration in range(repeats):
        started = time.perf_counter()
        timings = Counter()
        outputs = []
        per_chunk = []
        for chunk in chunks:
            chunk_started = time.perf_counter()
            if backend == "cpu":
                output = expand_cpu(chunk)
                phase = {"expansion": time.perf_counter() - chunk_started}
                device = {"backend": "cpu", "arithmetic": "Python arbitrary precision"}
            else:
                result = expander.expand(chunk, verify=False)
                output, phase, device = result.output, dict(result.timings), dict(result.metadata)
            before_verify = time.perf_counter()
            if output != reference_expand(chunk):
                raise ValueError("expansion differs from independent contextual replay")
            phase["verification"] = time.perf_counter() - before_verify
            device["verified"] = True
            device["external_verification"] = "same exact contextual CPU replay for both backends"
            phase["total"] = time.perf_counter() - chunk_started
            timings.update(phase)
            outputs.extend(output)
            per_chunk.append({"cases": len(chunk), "timings": phase, "metadata": device})
        verify = time.perf_counter()
        output_sha = digest(outputs)
        if output_sha != workload["expected_sha256"]:
            raise ValueError("expansion differs from captured scheduler output")
        timings["captured_output_verification"] = time.perf_counter() - verify
        timings["total"] = time.perf_counter() - started
        rows.append({"iteration": iteration, "timings": dict(timings),
                     "output_sha256": output_sha, "chunks": per_chunk})
    result = {"schema": "zipper-gpu-rewrite-benchmark-v1", "backend": backend,
              "mode": mode, "kind": kind, "workload": workload, "chunk_count": len(chunks),
              "largest_chunk_cases": max(map(len, chunks)),
              "loading_integrity_and_reference_replay_wall": 0.0 if reused_preparation else prepared.preparation_wall,
              "preparation_original_wall": prepared.preparation_wall,
              "preparation_reused": reused_preparation,
              "prepared_snapshot_sha256": prepared.payload_sha256,
              "prepared_snapshot_revalidation_wall": revalidation_wall,
              "setup_wall": setup_wall, "expander_reused": retained,
              "first_use_including_setup_wall": None if retained else rows[0]["timings"]["total"] + setup_wall,
              "timing_scope": "contextual catalog expansion only; CPU scheduling/geometric transport/collection excluded",
              "expected_projection_semantics": "selected scheduler scalar restored to contextual identities; optional single identities are conceptual full-exterior projections with verified source-order permutation",
              "single_identity_selection_scope": (
                  "excluded" if kind == "contextual" else
                  "conceptual weighted identities, including single-rule expansion beyond the single controller's literal operation"),
              "pooling_scope": "actual available batches" if mode == "sequential" else "offline pooling of successive scheduler times",
              "hardware": {"platform": platform.platform(), "python": platform.python_version()},
              "benchmark_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "expansion_source_sha256": prepared.source_sha256,
              "rows": rows}
    with sqlite3.connect(database.resolve().as_uri() + "?mode=rw", uri=True) as db:
        db.execute("BEGIN IMMEDIATE")
        identity_started = time.perf_counter()
        current_capture = capture_tables_identity(db)
        if current_capture != workload["prepared_capture_tables_sha256"]:
            raise ValueError("benchmark destination capture differs from prepared snapshot")
        result["capture_write_identity_validation_wall"] = time.perf_counter() - identity_started
        result["capture_write_binding"] = "prepared capture-table digest authenticated on the same BEGIN IMMEDIATE connection before benchmark insertion"
        db.execute("INSERT INTO benchmarks(backend,result) VALUES(?,?)", (backend, json.dumps(result)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    cap = commands.add_parser("capture")
    cap.add_argument("database", type=Path)
    cap.add_argument("--workloads", type=Path)
    cap.add_argument("--max-steps", type=int, default=20000)
    cap.add_argument("--max-terms", type=int, default=10000)
    cap.add_argument("--batch-size", type=int, default=1024)
    run = commands.add_parser("benchmark")
    run.add_argument("database", type=Path)
    run.add_argument("--backend", choices=["cpu", "opencl", "cpp"], required=True)
    run.add_argument("--repeats", type=int, default=3)
    run.add_argument("--roots", type=int, nargs="+")
    run.add_argument("--mode", choices=["sequential", "pooled"], default="sequential")
    run.add_argument("--kind", choices=["contextual", "used", "all"], default="contextual")
    stat = commands.add_parser("status")
    stat.add_argument("database", type=Path)
    args = parser.parse_args()
    if args.command == "capture":
        result = capture(args.database,
                         json.loads(args.workloads.read_text()) if args.workloads else DEFAULT_WORKLOADS,
                         max_steps=args.max_steps, max_terms=args.max_terms, batch_size=args.batch_size)
    elif args.command == "benchmark":
        result = benchmark(args.database, backend=args.backend, repeats=args.repeats,
                           roots=set(args.roots) if args.roots else None, mode=args.mode,
                           kind=args.kind)
    else:
        _, result = load_cases(args.database)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
