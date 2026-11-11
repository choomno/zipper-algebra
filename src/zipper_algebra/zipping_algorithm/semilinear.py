"""Derived mirrored rules with the exact semilinear symmetry q -> q^-1.

B preserves product order, sends sigma_l to sigma_l^-1 and q to q^-1.
Directly from the slider definition,
B(X_i)=(-1)^(i-1) q^(-i(i-1)/2) X_i, in the free braid-word algebra.
Every upward zipper term contains exactly one slider of the same rank, so
this common unit cancels. B(Z_i)=-q^(-(i-1))*B(X_i)/X_i*Z_i,
including the exceptional constant in Z_2. Thus B preserves the original
zipper ideal. B alone with q fixed is not the symmetry used here.

Applying B to a frozen M rule gives a BM rule with negative mirrored
friends. The frozen records remain unchanged. Half twists use the braid
identity B(Delta)=Delta^-1; their full retained factors are never omitted.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from typing import Any, Mapping

from .catalog import CatalogError, condition, integer

PROOF_VERSION = "zipper-semilinear-bar-v1"


def _negative(value: Any) -> Any:
    return -value if type(value) is int else {"op": "sub", "args": [0, deepcopy(value)]}


def _sliders(word: Mapping[str, Any], *, inverted: bool = False) -> list[Any]:
    kind = word.get("type")
    if kind == "slider":
        if inverted:
            raise CatalogError("semilinear zipper transformation cannot invert a slider")
        return [word["size"]]
    result: list[Any] = []
    if "word" in word:
        result.extend(_sliders(word["word"], inverted=inverted or kind == "inverse"))
    for child in word.get("items", []):
        result.extend(_sliders(child, inverted=inverted))
    return result


def semilinear_word(word: Mapping[str, Any]) -> dict[str, Any]:
    """Transform an explicit word AST, with a common slider unit factored out.

    Pure braid named words are barred explicitly. Reverse, reflection and bar
    wrappers commute with B. Inverse wrappers must not contain sliders.
    """
    node = deepcopy(dict(word))
    kind = node.get("type")
    if kind in ("identity", "slider"):
        return node
    if kind == "sigma":
        node["sign"] = _negative(node.get("sign", 1))
        return node
    if kind in ("fat", "block"):
        node["sign"] = _negative(node.get("sign", 1))
        return node
    if kind == "half_twist":
        return {"type": "inverse", "word": node}
    if kind == "concat":
        node["items"] = [semilinear_word(child) for child in node["items"]]
        return node
    if kind in ("inverse", "bar", "reverse", "reflect"):
        if kind == "inverse":
            _sliders(node)
        node["word"] = semilinear_word(node["word"])
        return node
    if kind == "word_ref":
        return {"type": "bar", "word": node}
    raise CatalogError("unknown word in semilinear zipper transformation")


def semilinear_zip_rule(rule: Mapping[str, Any], parameters: Mapping[str, int]) -> dict[str, Any]:
    """Derive a same-rank BM upward identity from an accepted frozen M record.

    Bindings keep their original upward meanings. Every term is checked to
    contain one rank-i slider before its common B(X_i) unit is canceled.
    The result records its exact source hash and the relator transport unit;
    no additional author review is inferred from this derivation.
    """
    source = deepcopy(dict(rule))
    env = dict(parameters)
    rid = source.get("id")
    if (not isinstance(rid, str) or not rid.startswith("zip.") or not rid.endswith(".M")
            or source.get("execution", {}).get("algebra_verified") is not True):
        raise CatalogError("an accepted frozen mirrored zipper rule is required")
    if (type(env.get("i")) is not int or env["i"] < 2
            or not condition(source["condition"], env)):
        raise CatalogError("valid original bindings and a true zipper guard are required")
    for side in ("lhs", "rhs"):
        if not source[side]:
            raise CatalogError("nonempty zipper polynomial required")
        for term in source[side]:
            ranks = _sliders(term["word"])
            if len(ranks) != 1 or integer(ranks[0], env) != env["i"]:
                raise CatalogError("every zipper term must contain one slider of the same rank")
    source_hash = hashlib.sha256(json.dumps(source, sort_keys=True, separators=(",", ":"),
                                           allow_nan=False).encode()).hexdigest()
    result = deepcopy(source)
    result["id"] = rid[:-2] + ".BM"
    result["source_record_id"] = rid
    result["source_rule_sha256"] = source_hash
    for side in ("lhs", "rhs"):
        for term in result[side]:
            term["word"] = semilinear_word(term["word"])
            for coefficient in term["coefficient"]:
                coefficient["q_power"] = _negative(coefficient["q_power"])
    i = env["i"]
    result["semilinear"] = {
        "proof_version": PROOF_VERSION, "source_record_id": rid,
        "source_rule_sha256": source_hash, "parameters": env,
        "coefficient_map": "q -> q^-1", "crossing_map": "sigma_l -> sigma_l^-1",
        "preserves_product_order": True,
        "common_slider_unit": {"integer": (-1) ** (i - 1), "q_power": -i * (i - 1) // 2},
        "normalized_relator_unit": {"integer": -1, "q_power": -(i - 1)},
        "half_twist_identity": "bar(Delta) = inverse(Delta)",
        "author_review": "inherits partial source review; semilinear derivation recorded separately",
    }
    result["symmetry"] = dict(result["symmetry"], name="BM", invert_q=True,
                              invert_crossing_sign=True, semilinear=True)
    result["execution"] = dict(result["execution"], derived_semilinear=True)
    # The M source carries base-notation cofactors, not directly the transformed
    # cofactors. Preserve that evidence under an explicit source label.
    for key in ("original_relator_cofactor", "cofactor_original_notation", "cofactor_transport"):
        if key in result:
            result["semilinear_source_" + key] = result.pop(key)
    return result
