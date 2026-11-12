"""Reconstruct current-procedure output from certified input plans and leaf runs.

These routines validate bindings and structural witnesses, not rewrite proofs.
Replay leaf certificates independently before treating them as verified results.
Only completed runs from one explicitly supplied execution identity are accepted.
The result is a canonical terminal sum, never a fabricated run certificate.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import hashlib
import json
from typing import Any, Mapping

from .input_normal_form import POLICY, leaf_requests, reconstruction_plan
from .output_normal_form import canonical_output, engine_output, from_certificate

_CONTEXT_FIELDS = {"rank", "size", "component_word", "top_route", "bottom_route", "record"}
_RECORD_FIELDS = {"bottom_ports", "top_ports", "original_bottom_ports", "original_top_ports",
                  "twist_powers", "orientation", "top_route", "bottom_route",
                  "canonical_top_route", "canonical_bottom_route", "canonical_parities",
                  "inherited_top_route", "inherited_bottom_route", "input_hierarchy", "output_hierarchy"}


def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _sha(value: Any) -> str:
    return hashlib.sha256(_json(value).encode("ascii")).hexdigest()


def _hash(value: Any) -> bool:
    return type(value) is str and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _identity(identity: Any) -> dict[str, Any]:
    fields = {"schema", "source_sha256", "catalog_sha256", "execution"}
    if (type(identity) is not dict or set(identity) != fields
            or type(identity["schema"]) is not str
            or type(identity["source_sha256"]) is not dict or not identity["source_sha256"]
            or any(type(k) is not str or not _hash(v) for k, v in identity["source_sha256"].items())
            or not _hash(identity["catalog_sha256"]) or type(identity["execution"]) is not dict):
        raise ValueError("explicit source/catalog execution identity required")
    # Execution configuration is an exact JSON object, with no floats/custom objects.
    def exact(value: Any) -> None:
        if value is None or type(value) in (str, bool, int):
            return
        if type(value) is list:
            for x in value:
                exact(x)
            return
        if type(value) is dict and all(type(k) is str for k in value):
            for x in value.values():
                exact(x)
            return
        raise ValueError("execution identity requires exact JSON values")
    exact(identity["execution"])
    return deepcopy(identity)


def completed_output(result: Any) -> dict[str, Any]:
    """Structurally import a standalone certificate or native completed DAG result.

    DAG results keep their own schema and cannot be used as standalone traces.
    Neither adapter independently verifies mathematical correctness.
    """
    if type(result) is dict and result.get("schema") == "zipper-native-reuse-v1":
        if (result.get("scope") != "literal_word_native_reused_spanning_procedure"
                or result.get("status") != "complete"
                or result.get("coefficient_ring") != "Z[q,q^-1]"
                or result.get("trace_enabled") is not False
                or result.get("pending", []) != []
                or type(result.get("execution")) is not dict
                or result["execution"].get("schedule") != "intrinsic-completion-dag-v1"):
            raise ValueError("invalid completed native reuse result")
        return canonical_output(result.get("output"), n=result["n"], complete=True)
    return from_certificate(result)


def execution_identity(certificate: Any) -> dict[str, Any]:
    """Extract strict source/catalog/configuration binding, without replaying.

    Python and native producers have different identities. Configuration is
    retained conservatively, including batch size; this does not claim it changes
    completed mathematical outputs. It prevents implicit cross-policy reuse.
    """
    completed_output(certificate)
    return _identity({"schema": certificate["schema"],
                      "source_sha256": certificate.get("source_sha256"),
                      "catalog_sha256": certificate.get("catalog_sha256"),
                      "execution": certificate.get("execution", {})})


def _request(request: Any) -> dict[str, Any]:
    if (type(request) is not dict or set(request) != {"key", "n", "word"}
            or type(request["n"]) is not int or request["n"] < 1
            or type(request["word"]) is not list
            or any(type(x) is not int or not 0 < abs(x) < request["n"] for x in request["word"])):
        raise ValueError("invalid primitive request")
    expected = _sha({"policy": POLICY, "rank": request["n"], "word": request["word"]})
    if request["key"] != expected:
        raise ValueError("primitive request key disagrees with its literal input")
    return request


def execution_cache_key(request: Any, source_identity: Any) -> str:
    """Bind a mathematical primitive request to a particular producer policy."""
    return _sha({"schema": "zipper-current-execution-cache-v1", "request": _request(request),
                 "execution_identity": _identity(source_identity)})


def _supported(term: dict[str, Any]) -> None:
    if set(term) != {"word", "components", "coefficient"}:
        raise ValueError("cannot embed unknown terminal-term extension fields")
    for context in term["components"]:
        if set(context) != _CONTEXT_FIELDS or not set(context["record"]) <= _RECORD_FIELDS:
            raise ValueError("cannot embed unknown component extension fields")
        if set(context["record"]["canonical_parities"]) != {"upper", "lower_inverse"}:
            raise ValueError("cannot embed unknown parity extension fields")


def _singleton(bottom: int, top: int) -> dict[str, Any]:
    return {"size": 1, "component_word": [], "record": {
        "original_bottom_ports": [bottom], "original_top_ports": [top],
        "twist_powers": [], "orientation": "normal", "top_route": [], "bottom_route": [],
        "canonical_top_route": [], "canonical_bottom_route": [],
        "canonical_parities": {"upper": [], "lower_inverse": []}}}


def _over_route(n: int, selected: list[int]) -> list[int]:
    labels = list(range(1, n + 1))
    chronological = []
    for destination, label in enumerate(selected):
        position = labels.index(label)
        while position > destination:
            chronological.append(-position)
            labels[position - 1], labels[position] = labels[position], labels[position - 1]
            position -= 1
    return chronological[::-1]


def _shift(word: list[Any], amount: int) -> list[Any]:
    return [(1 if x > 0 else -1) * (abs(x) + amount) if type(x) is int
            else [x[0], x[1], x[2] + amount, *x[3:]] for x in word]


def _assemble(components: list[dict[str, Any]], *, n: int, scalar: dict[int, int]) -> dict[str, Any]:
    top, bottom = list(range(1, n + 1)), list(range(1, n + 1))
    contexts = deepcopy(components)
    for context in contexts:
        context["rank"] = len(bottom)
        for side, labels in (("top", top), ("bottom", bottom)):
            original = context["record"]["original_" + side + "_ports"]
            local = [labels.index(p) + 1 for p in original]
            context["record"][side + "_ports"] = local
            route = _over_route(len(labels), local)
            context[side + "_route"] = [-a for a in reversed(route)] if side == "top" else route
            labels[:] = [p for p in labels if p not in original]
    # The checked plan and leaf partitions cover every endpoint. The final
    # canonical_output call independently validates the assembled partition.
    word: list[Any] = []
    for context in reversed(contexts):
        word = (context["top_route"] + context["component_word"]
                + _shift(word, context["size"]) + context["bottom_route"])
    return {"coefficient": [[q, c] for q, c in sorted(scalar.items()) if c],
            "word": word, "components": contexts}


def reconstruct_output(form: Any, leaf_results: Mapping[str, Any], *,
                       source_identity: Any) -> dict[str, Any]:
    """Reconstruct the current output using completed, explicitly bound leaf runs.

    ``leaf_results`` maps each primitive request key to its original run
    certificate, whose rank, literal input and producer identity must match.
    Extra or missing keys, partial runs and unknown mathematical extensions are
    rejected. Certificates are structurally checked here, not independently
    replayed. Store their verification evidence separately. An empty recipe
    needs no leaf runs but still takes an explicit target execution identity.
    """
    identity = _identity(source_identity)
    requests = leaf_requests(form)
    plan = reconstruction_plan(form)
    if not isinstance(leaf_results, Mapping) or set(leaf_results) != {r["key"] for r in requests}:
        raise ValueError("leaf result keys must exactly match the primitive requests")
    decoded: dict[str, list[dict[str, Any]]] = {}
    for request in requests:
        certificate = leaf_results[request["key"]]
        normal = completed_output(certificate)
        if normal["completion"] != "complete":
            raise ValueError("reconstruction requires completed primitive runs")
        literal = certificate.get("input_word")
        if (certificate["n"] != request["n"] or type(literal) is not list
                or any(type(a) is not int for a in literal) or literal != request["word"]):
            raise ValueError("leaf result does not match the requested literal input and rank")
        if _json(execution_identity(certificate)) != _json(identity):
            raise ValueError("leaf result producer identity disagrees with the expected identity")
        output = engine_output(normal)
        for term in output:
            _supported(term)
        decoded[request["key"]] = output
    products: list[tuple[dict[int, int], list[dict[str, Any]]]] = [({0: 1}, [])]
    for atom in plan["atoms"]:
        if atom["kind"] == "singleton":
            piece = _singleton(atom["bottom"], atom["top"])
            products = [(scalar, components + [piece]) for scalar, components in products]
            continue
        following = []
        for scalar, components in products:
            for term in decoded[atom["key"]]:
                coefficient: defaultdict[int, int] = defaultdict(int)
                for p, a in scalar.items():
                    for q, b in term["coefficient"]:
                        coefficient[p + q] += a * b
                combined = {q: c for q, c in coefficient.items() if c}
                # Canonical leaf coefficients are nonzero; their product in
                # Z[q,q^-1] is nonzero, although individual powers may cancel.
                embedded = deepcopy(term["components"])
                for component in embedded:
                    for side in ("top", "bottom"):
                        field = "original_" + side + "_ports"
                        component["record"][field] = [atom[side + "_labels"][p - 1]
                                                       for p in component["record"][field]]
                following.append((combined, components + embedded))
        products = following
    terms = [_assemble(components, n=plan["n"], scalar=scalar) for scalar, components in products]
    return canonical_output(terms, n=plan["n"], complete=True)
