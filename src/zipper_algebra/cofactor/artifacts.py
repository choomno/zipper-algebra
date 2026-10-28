"""Versioned mathematical data bound to an explicit presentation and target.

An artifact records data and their interpretation. Reading one does not prove its
identity, certify a basis, or validate a fitted coefficient function. Producers
must perform the corresponding mathematical checks before reporting success.
Historical unbound formats remain available through their existing explicit
readers; this module accepts only the bound envelope defined here.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import tempfile

from .atoms import AtomSet
from .cofactor_io import Cofactor, symbolic_target, validate_word
from .lift import Lift
from .skeleton import Skeleton


SCHEMA = "zipper-mathematical-artifact-v1"


def _canonical(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False).encode()


def _digest(value) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def presentation_identity(pres) -> dict:
    """Identify ordered defining data, independently of any monomial order."""
    canonical = pres.canonical_text()
    return {"schema": "zipper-presentation-v1", "coefficient_domain": "Q[q,t]",
            "generators": pres.var_names(), "relation_labels": pres.labels(),
            "canonical_relations": canonical,
            "sha256": hashlib.sha256(canonical.encode()).hexdigest(),
            "legacy_hash": pres.content_hash()}


def _numerator(poly) -> list:
    return [[list(degree), str(Fraction(coefficient))]
            for degree, coefficient in sorted(poly.items())]


def target_identity(pres, target) -> dict:
    """Identify the actual polynomial; equivalent target names have one identity."""
    polynomial = symbolic_target(pres, target)
    terms = [[list(word), _numerator(coefficient)]
             for word, coefficient in sorted(polynomial.items())]
    return {"coefficient_domain": "Q[q,t]", "terms": terms, "sha256": _digest(terms)}


def _kind(value) -> str:
    for cls, name in ((Skeleton, "skeleton"), (Lift, "lift"), (Cofactor, "cofactor")):
        if isinstance(value, cls):
            return name
    raise TypeError("an artifact must contain a Skeleton, Lift or Cofactor")


def _nonnegative(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")


def _validate(value, pres, target):
    kind = _kind(value)
    if kind == "cofactor":
        value.validate(pres)
    else:
        if (len(value.anchor) != 3 or any(type(x) is not int for x in value.anchor)
                or value.anchor[2] != 0):
            raise ValueError("an artifact anchor must have integer q,t and characteristic zero")
        if kind == "lift":
            # These are the same domain and reference checks used before composition.
            from .compose import _order, _validate_lift
            _validate_lift(value, pres)
            _order(value.entries)
            if any(not isinstance(key, str) or not key for key in value.entries):
                raise ValueError("artifact lift keys must be nonempty strings")
            for index in value.anchor_index.values():
                if index is not None:
                    _nonnegative(index, "anchor basis index")
            _nonnegative(value.samples_used, "sample count")
            _nonnegative(value.holdout, "held-out count")
        else:
            _nonnegative(value.remainder_terms, "remainder term count")
            if (any(not isinstance(key, str) or not key for key in value.basis_keys)
                    or len(set(value.basis_keys)) != len(value.basis_keys)):
                raise ValueError("skeleton basis keys must be distinct nonempty strings")
            for step in value.steps:
                if step["basis_key"] not in value.basis_keys:
                    raise ValueError("a skeleton step must name a key in its proof closure")
                validate_word(step["left"], pres.n)
                validate_word(step["right"], pres.n)
                if isinstance(step["mult"], bool) or not isinstance(step["mult"], (int, Fraction)):
                    raise ValueError("skeleton multipliers must be exact rational numbers")
    recorded = getattr(value, "presentation_hash", "")
    if recorded and recorded != pres.content_hash():
        raise ValueError("artifact data belong to a different presentation")
    named = getattr(value, "target", "")
    if named and target_identity(pres, named) != target_identity(pres, target):
        raise ValueError("artifact data name a different target")


def _payload(value) -> dict:
    if not isinstance(value, Cofactor):
        return value.to_json()
    return {"n": value.n, "atoms": value.atoms.to_json(),
            "presentation": value.presentation, "presentation_hash": value.presentation_hash,
            "target": value.target, "producer": value.producer,
            "anchor": list(value.anchor) if value.anchor is not None else None,
            "notes": list(value.notes),
            "terms": [{"left": list(left), "relation": relation, "right": list(right),
                       "num": _numerator(num), "denom": list(den)}
                      for left, relation, right, num, den in value.terms]}


def _decode(kind, payload):
    if kind == "skeleton":
        return Skeleton.from_json(payload)
    if kind == "lift":
        return Lift.from_json(payload)
    if kind != "cofactor":
        raise ValueError("unknown mathematical artifact kind")
    terms = [(tuple(term["left"]), term["relation"], tuple(term["right"]),
              {tuple(degree): Fraction(coef) for degree, coef in term["num"]},
              tuple(term["denom"])) for term in payload["terms"]]
    return Cofactor(terms, payload["n"], AtomSet.from_json(payload["atoms"]),
                    payload["presentation"], payload["presentation_hash"],
                    payload["target"], payload["producer"], payload["anchor"], payload["notes"])


def _order_identity(kind, order):
    if kind == "cofactor":
        if order is not None:
            raise ValueError("a cofactor identity does not depend on a monomial order")
    elif not isinstance(order, (str, dict)) or not order:
        raise ValueError("skeleton and lift artifacts require an explicit monomial order")
    # Copy through JSON so mutable caller data cannot alter the envelope later.
    return json.loads(_canonical(order))


def write_json_atomic(path, data) -> None:
    """Replace one JSON artifact atomically; failed writes preserve the old file."""
    content = json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + "\n"
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    stream = tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix="." + path.name + ".", delete=False)
    temporary = Path(stream.name)
    try:
        with stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def save_artifact(path, value, pres, target, *, order=None) -> None:
    """Write validated, bound data; no verification outcome is implied."""
    _validate(value, pres, target)
    kind = _kind(value)
    payload = _payload(value)
    data = {"schema": SCHEMA, "kind": kind, "presentation": presentation_identity(pres),
            "target": target_identity(pres, target), "order": _order_identity(kind, order),
            "payload": payload, "payload_sha256": _digest(payload)}
    write_json_atomic(path, data)


def load_artifact(path, pres, target, *, kind=None, order=None, atoms=None):
    """Read bound data with the requested interpretation and optional atom family.

    Digests detect content changes, not authorship or mathematical correctness.
    Unknown or historical unbound schemas are rejected; use their explicit readers.
    """
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    fields = {"schema", "kind", "presentation", "target", "order", "payload", "payload_sha256"}
    if not isinstance(data, dict) or set(data) != fields or data["schema"] != SCHEMA:
        raise ValueError("unknown or unbound mathematical artifact schema")
    actual_kind = data["kind"]
    if kind is not None and actual_kind != kind:
        raise ValueError("mathematical artifact kind mismatch")
    if (_canonical(data["presentation"]) != _canonical(presentation_identity(pres))
            or _canonical(data["target"]) != _canonical(target_identity(pres, target))):
        raise ValueError("artifact presentation or target does not match the requested data")
    if _canonical(data["order"]) != _canonical(_order_identity(actual_kind, order)):
        raise ValueError("artifact monomial order does not match")
    payload = data["payload"]
    if data["payload_sha256"] != _digest(payload):
        raise ValueError("artifact payload digest mismatch")
    value = _decode(actual_kind, payload)
    _validate(value, pres, target)
    if _canonical(_payload(value)) != _canonical(payload):
        raise ValueError("artifact payload is not in canonical exact form")
    if atoms is not None:
        atoms.validate()
        if not hasattr(value, "atoms") or _canonical(value.atoms.to_json()) != _canonical(atoms.to_json()):
            raise ValueError("artifact denominator atoms do not match")
    return value
