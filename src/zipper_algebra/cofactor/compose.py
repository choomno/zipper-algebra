#!/usr/bin/env python3
"""Compose a Q(q,t) cofactor of the target over the defining relations.

Phase 4.  Input: a Lift (each needed basis element's non-expanded proof, with every coefficient
a rational function) and a Skeleton (the target as a combination of basis elements).  Output: a
Cofactor -- a finite list of tuples (L, relation, R, coefficient) with

    sum_j  c_j(q,t) * L_j * rel_{i_j} * R_j  =  target        in Q(q,t)<x_1..y_{n-1}>.

Two substitutions happen here.  First each basis element's proof is expanded recursively over the
*input relations* (a proof entry may reference another basis element).  Then the skeleton's steps
are applied.  Nothing is specialised: the arithmetic is exact in Q(q,t), with denominators carried
as exponent vectors over the AtomSet the lift was fitted with.

Ported from compose_x3_cofactor_qt.py, which was X_3-at-n=4 only in three ways, all now parameters:
the basis indices 16/17/102 and the multipliers q^3, -q^2, q(q-1) were typed into the source; the
denominator exponent tuples were positional over six hard-coded atoms; and the output path was
fixed.  The rational arithmetic is deliberately unchanged -- LCM cross-multiplication with no GCD
reduction, zero numerators dropped only at the end of each bucket and at the end of the
composition -- because the shipped certificate's denominators are exactly what that produces.

**The multipliers are lifted, not assumed.**  Each is a ratio of leading coefficients, hence a
rational function like any other; at n=4/X_3 they happen to be polynomials in q alone, and
WORKFLOW.md records that this is algebra-dependent.  `lift_multipliers` fits them from the same
sample set by the same `fit_rational`, and the composition uses the fitted values.  Passing them
in is supported for reproducing a legacy run.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from collections import defaultdict
from collections.abc import Mapping
from fractions import Fraction

from .atoms import AtomSet
from .cofactor_io import (Cofactor, resolve_target,
                         validate_word, validate_qt_polynomial)
from .lift import NeedMoreSamples, fit_rational
from ..presentation import Presentation
from ..coefficients import NumericRing


# ---------------------------------------------------------------------------
# Q(q,t) rationals as (numerator {(dq,dt): Fraction}, denominator exponent tuple)
# ---------------------------------------------------------------------------

def poly_add(a: dict, b: dict) -> dict:
    out = dict(a)
    for m, c in b.items():
        v = out.get(m, Fraction(0)) + c
        if v:
            out[m] = v
        elif m in out:
            del out[m]
    return out


def poly_mul(a: dict, b: dict) -> dict:
    out: dict = defaultdict(Fraction)
    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            out[(i1 + i2, j1 + j2)] += c1 * c2
    return {m: c for m, c in out.items() if c}


def poly_mul_atoms(num: dict, exps, atoms: AtomSet) -> dict:
    """num * prod_k atom_k^{exps_k}, expanded."""
    out = num
    for atom, e in zip(atoms.atoms, exps):
        for _ in range(e):
            nxt: dict = defaultdict(Fraction)
            for (mi, mj), c0 in out.items():
                for (ai, aj), ac in atom.coefs.items():
                    nxt[(mi + ai, mj + aj)] += c0 * ac
            out = dict(nxt)
    return {m: c for m, c in out.items() if c}


class Rat:
    """The (num, den) pair plus its AtomSet, so the arithmetic needs no globals."""

    __slots__ = ("atoms",)

    def __init__(self, atoms: AtomSet):
        self.atoms = atoms

    def zero(self):
        return ({}, self.atoms.zero())

    def const(self, c):
        c = Fraction(c)
        return ({(0, 0): c}, self.atoms.zero()) if c else self.zero()

    def add(self, r1, r2):
        n1, d1 = r1
        n2, d2 = r2
        d = AtomSet.lcm(d1, d2)
        n1s = poly_mul_atoms(n1, tuple(a - b for a, b in zip(d, d1)), self.atoms)
        n2s = poly_mul_atoms(n2, tuple(a - b for a, b in zip(d, d2)), self.atoms)
        return (poly_add(n1s, n2s), d)

    def mul(self, r1, r2):
        return (poly_mul(r1[0], r2[0]), AtomSet.mul(r1[1], r2[1]))

    def eval(self, r, q, t):
        num, den = r
        if not num:
            return Fraction(0)
        n = sum(c * Fraction(q) ** i * Fraction(t) ** j for (i, j), c in num.items())
        return n / Fraction(self.atoms.eval_denom(den, q, t))


# ---------------------------------------------------------------------------
# Expansion of the lifted proofs over the input relations
# ---------------------------------------------------------------------------

def _order(entries: dict) -> list:
    """Keys in an order where every basis reference precedes its referrer.

    The legacy code relied on ascending integer index (proof entries referencing only smaller
    indices).  Keys are leading words, so the order has to be computed; a genuine cycle is an
    error, a forward reference is not.
    """
    state, out = {}, []

    def visit(k):
        st = state.get(k)
        if st == 2:
            return
        if st == 1:
            raise ValueError(f"cycle in the proof DAG at {k!r}")
        state[k] = 1
        for _den, _num, _l, kind, ref, _r in entries.get(k, ()):
            if kind == "basis" and ref in entries:
                visit(ref)
        state[k] = 2
        out.append(k)

    for k in entries:
        visit(k)
    return out


def _validate_context(word, pres=None):
    if pres is not None:
        validate_word(word, pres.n)
    elif (not isinstance(word, (tuple, list))
          or any(type(v) is not int or v < 1 for v in word)):
        raise ValueError("proof contexts must contain positive integer variable indices")


def _validate_lift(lift, pres=None):
    """Validate every proof entry before zero terms can hide malformed data."""
    lift.atoms.validate()
    if not isinstance(lift.entries, Mapping):
        raise ValueError("lift entries must be a mapping")
    limit = len(pres.labels()) if pres is not None else None
    for key, terms in lift.entries.items():
        if not isinstance(terms, (tuple, list)):
            raise ValueError(f"lift entry {key!r} must contain a list of proof terms")
        for term in terms:
            if not isinstance(term, (tuple, list)) or len(term) != 6:
                raise ValueError(f"lift entry {key!r} has a malformed proof term")
            den, num, left, kind, ref, right = term
            lift.atoms.validate_exponents(den)
            validate_qt_polynomial(num)
            _validate_context(left, pres)
            _validate_context(right, pres)
            if kind == "input":
                if type(ref) is not int or ref < 0 or (limit is not None and ref >= limit):
                    raise ValueError(f"lift entry {key!r} has invalid input relation index {ref!r}")
            elif kind == "basis":
                try:
                    exists = ref in lift.entries
                except TypeError:
                    raise ValueError(f"lift entry {key!r} has an invalid basis key") from None
                if not exists:
                    raise KeyError(f"{key!r} references {ref!r}, which is not in the lift")
            else:
                raise ValueError(f"lift entry {key!r} has unknown proof kind {kind!r}")


def expand_over_inputs(lift, keys=None, pres=None) -> dict:
    """{key: {(left_word, input_index, right_word): rat}} -- each needed basis element written
    over the input relations alone."""
    _validate_lift(lift, pres)
    rat = Rat(lift.atoms)
    entries = lift.entries
    keys = set(entries) if keys is None else set(keys)
    expanded: dict = {}
    for k in _order(entries):
        if k not in keys and not any(k == ref for _e in entries.values()
                                     for (_d, _n, _l, kd, ref, _r) in _e if kd == "basis"):
            continue
        bucket: dict = {}
        for den, num, left, kind, ref, right in entries[k]:
            coef = (dict(num), tuple(den))
            if kind == "input":
                key = (tuple(left), ref, tuple(right))
                bucket[key] = rat.add(bucket.get(key, rat.zero()), coef)
            else:
                for (s_l, s_i, s_r), s_rat in expanded[ref].items():
                    nk = (tuple(left) + s_l, s_i, s_r + tuple(right))
                    bucket[nk] = rat.add(bucket.get(nk, rat.zero()), rat.mul(coef, s_rat))
        # zeros are dropped here and not inside the loop: an entry that cancels has already
        # raised the LCM of its key, and the shipped denominators reflect that.
        expanded[k] = {kk: v for kk, v in bucket.items() if v[0]}
    return expanded


# ---------------------------------------------------------------------------
# The skeleton multipliers, as rational functions
# ---------------------------------------------------------------------------

def lift_multipliers(skeleton, samples, atoms: AtomSet, max_pq=5, max_pt=5) -> list:
    """Fit each skeleton step's multiplier over Q(q,t) from its values at the sample points.

    samples: {(q, t): [Fraction, ...]}, in skeleton-step order, as produced by
    re-reducing the target with the same contextual step sequence at each point.
    A basis-key mapping is accepted for legacy callers only when all step keys are unique:
    it cannot identify separate occurrences of the same basis element.
    """
    keys = [s["basis_key"] for s in skeleton.steps]
    ordered_samples = []
    for point, vals in sorted(samples.items()):
        if isinstance(vals, Mapping):
            if len(set(keys)) != len(keys):
                raise ValueError("basis-key multiplier samples are ambiguous for repeated "
                                 "basis elements; use values in skeleton-step order")
            if set(vals) != set(keys):
                raise ValueError(f"multiplier sample {point!r} does not carry every step key")
            vals = [vals[k] for k in keys]
        else:
            vals = list(vals)
        if len(vals) != len(keys):
            raise ValueError(f"multiplier sample {point!r} has {len(vals)} values for "
                             f"{len(keys)} skeleton steps")
        ordered_samples.append((point, vals))

    out = []
    for index, key in enumerate(keys):
        pts = [(point, vals[index]) for point, vals in ordered_samples]
        if len(pts) < 3:
            raise NeedMoreSamples(f"only {len(pts)} points carry a multiplier for step {index} "
                             f"({key!r})")
        res = fit_rational(pts, atoms, max_pq=max_pq, max_pt=max_pt)
        if res is None:
            raise NeedMoreSamples(f"no model in the current class for the multiplier of step "
                             f"{index} ({key!r})")
        den, num = res
        out.append((dict(num), tuple(den)))
    return out


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------

def compose(lift, skeleton, pres: Presentation, multipliers=None, samples=None,
            producer="zipper_algebra.cofactor.compose", verbose=True) -> Cofactor:
    """The target's cofactor over the input relations, as a Cofactor."""
    if getattr(skeleton, "remainder_terms", 0) != 0:
        raise ValueError("cannot compose a membership certificate from a nonzero remainder")
    for name, data in (("skeleton", skeleton), ("lift", lift)):
        recorded_hash = getattr(data, "presentation_hash", "")
        recorded_presentation = getattr(data, "presentation", "")
        if recorded_hash and recorded_hash != pres.content_hash():
            raise ValueError(f"{name} belongs to a different presentation (hash mismatch)")
        if recorded_presentation and recorded_presentation != pres.slug():
            raise ValueError(f"{name} belongs to presentation {recorded_presentation!r}, "
                             f"not {pres.slug()!r}")
    atoms = lift.atoms
    rat = Rat(atoms)
    for step in skeleton.steps:
        if not isinstance(step, Mapping) or any(k not in step for k in ("basis_key", "left", "right")):
            raise ValueError("a skeleton step needs a basis key and left/right contexts")
        _validate_context(step["left"], pres)
        _validate_context(step["right"], pres)
        try:
            hash(step["basis_key"])
        except TypeError:
            raise ValueError("skeleton basis keys must be hashable") from None
    keys = [s["basis_key"] for s in skeleton.steps]
    missing = [k for k in keys if k not in lift.entries]
    if missing:
        raise KeyError(f"the lift is missing skeleton elements {missing}")

    expanded = expand_over_inputs(lift, pres=pres)
    if verbose:
        for k in keys:
            print(f"  {k}: {len(expanded[k])} R-input terms")

    if multipliers is None:
        if samples is None:
            raise ValueError("compose needs either multipliers= or samples= to lift them from")
        multipliers = lift_multipliers(skeleton, samples, atoms)
    multipliers = list(multipliers)
    if len(multipliers) != len(skeleton.steps):
        raise ValueError(f"{len(multipliers)} multipliers for {len(skeleton.steps)} "
                         "skeleton steps")
    for mult in multipliers:
        if not isinstance(mult, (tuple, list)) or len(mult) != 2:
            raise ValueError("each multiplier must have a numerator and denominator vector")
        validate_qt_polynomial(mult[0])
        atoms.validate_exponents(mult[1])

    acc: dict = {}
    for mult, step in zip(multipliers, skeleton.steps):
        for (left, idx, right), r in expanded[step["basis_key"]].items():
            k2 = (tuple(step["left"]) + left, idx, right + tuple(step["right"]))
            acc[k2] = rat.add(acc[k2], rat.mul(mult, r)) if k2 in acc else rat.mul(mult, r)
    acc = {k: v for k, v in acc.items() if v[0]}

    labels = pres.labels()
    terms = [(left, labels[idx], right, num, den)
             for (left, idx, right), (num, den) in
             sorted(acc.items(), key=lambda x: (x[0][1], len(x[0][0]), x[0][0],
                                                len(x[0][2]), x[0][2]))]
    if verbose:
        print(f"{skeleton.target}: {len(terms)} (left, relation, right) terms")
    return Cofactor.from_presentation(pres, skeleton.target, atoms, terms, producer,
                                      anchor=skeleton.anchor)


def numeric_check(cof: Cofactor, pres: Presentation, points=((5, 7), (3, 11), (7, 13),
                                                             (11, 17), (10, 7)),
                  target=None, verbose=True) -> bool:
    """Substitute the input relations back at several (q,t) and compare with the target.

    A specialised check of a symbolic object: it can only refute, never certify.  The
    certificate is verify.check.
    """
    ok_all = True
    for q, t in points:
        ring = NumericRing(q, t, 0)
        tgt = resolve_target(pres, target if target is not None else cof.target, ring)
        acc = cof.evaluate(q, t, pres)
        diff = acc - tgt
        ok = not diff.terms
        ok_all &= ok
        if verbose:
            print(f"  ({q},{t}): {'OK' if ok else f'MISMATCH ({len(diff.terms)} terms)'}")
    return ok_all


MULTIPLIER_SCHEMA = "zipper-cofactor-multipliers-v1"


def _metadata_json(value):
    """Compare JSON metadata without equating booleans, integers and floats."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def multipliers_to_json(multipliers, skeleton, atoms, pres):
    """Encode exact rational multipliers in contextual skeleton-occurrence order.

    The full SHA-256 binds the current canonical presentation bytes without changing
    its legacy display hash. Ordered atom definitions and contexts bind coefficient
    interpretation and the particular reduction occurrences.
    """
    atoms.validate()
    records = []
    for (num, den), step in zip(multipliers, skeleton.steps, strict=True):
        validate_qt_polynomial(num)
        atoms.validate_exponents(den)
        validate_word(step["left"], pres.n)
        validate_word(step["right"], pres.n)
        records.append({"basis_key": step["basis_key"], "left": list(step["left"]),
                        "right": list(step["right"]), "denom": list(den),
                        "num": [[list(m), str(c)] for m, c in sorted(num.items())]})
    return {"schema": MULTIPLIER_SCHEMA,
            "presentation_sha256": hashlib.sha256(pres.canonical_text().encode()).hexdigest(),
            "target": skeleton.target, "atoms": atoms.to_json(), "steps": records}


def load_multipliers(path, skeleton, atoms, pres):
    """Read explicit symbolic multipliers, rejecting mismatched or ambiguous data."""
    with open(path) as source:
        data = json.load(source)
    return multipliers_from_json(data, skeleton, atoms, pres)


def multipliers_from_json(data, skeleton, atoms, pres):
    """Decode a multiplier snapshot with its presentation and occurrence checks."""
    if not isinstance(data, dict) or data.get("schema") != MULTIPLIER_SCHEMA:
        raise ValueError("unknown multiplier data schema")
    digest = hashlib.sha256(pres.canonical_text().encode()).hexdigest()
    if data.get("presentation_sha256") != digest or data.get("target") != skeleton.target:
        raise ValueError("multipliers belong to a different presentation or target")
    if _metadata_json(data.get("atoms")) != _metadata_json(atoms.to_json()):
        raise ValueError("multiplier denominator atoms do not match the lift")
    records = data["steps"]
    if not isinstance(records, list) or len(records) != len(skeleton.steps):
        raise ValueError("multipliers must match the number of contextual skeleton occurrences")
    out = []
    for record, step in zip(records, skeleton.steps):
        expected = {"basis_key": step["basis_key"], "left": list(step["left"]),
                    "right": list(step["right"])}
        if (not isinstance(record, dict)
                or _metadata_json({k: record.get(k) for k in expected}) != _metadata_json(expected)):
            raise ValueError("multiplier occurrence does not match the reduction skeleton")
        num = {}
        for monomial, coefficient in record["num"]:
            monomial = tuple(monomial)
            if not isinstance(coefficient, str) or monomial in num:
                raise ValueError("multiplier numerators need unique powers and exact rational strings")
            num[monomial] = Fraction(coefficient)
        validate_qt_polynomial(num)
        atoms.validate_exponents(record["denom"])
        out.append((num, tuple(record["denom"])))
    return out


def _save_atomic(cofactor, output, *, producer=None, inputs=()):
    from ..provenance import capture_provenance, write_output_provenance
    producer = capture_provenance() if producer is None else producer
    output = Path(output)
    descriptor, temporary = tempfile.mkstemp(prefix="." + output.name + "-", suffix=".tmp",
                                             dir=output.parent)
    os.close(descriptor)
    try:
        cofactor.save(temporary, provenance=False)
        os.replace(temporary, output)
        write_output_provenance(output, producer=producer, inputs=inputs)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def add_arguments(parser):
    """Shared installed/module command arguments for exact symbolic composition."""
    from ..cli import _presentation_arguments
    parser.add_argument("--lift", required=True)
    parser.add_argument("--skeleton", required=True)
    parser.add_argument("--multipliers", required=True,
                        help="versioned symbolic multipliers in contextual occurrence order")
    parser.add_argument("--out", required=True)
    parser.add_argument("--target", help="explicit target required for bound artifact files")
    parser.add_argument("--solver-order", choices=("degree-lex",), default="degree-lex")
    _presentation_arguments(parser, require_n=True)


def _composition_inputs(args, pres, documents=None):
    from .artifacts import artifact_from_json, target_identity
    from .lift import Lift
    from .skeleton import Skeleton
    if documents is None:
        documents = [json.loads(Path(path).read_bytes()) for path in (args.skeleton, args.lift)]
    if any(not isinstance(data, dict) for data in documents):
        raise ValueError("composition inputs must contain JSON objects")
    bound = any("schema" in data for data in documents)
    if bound and args.target is None:
        raise ValueError("bound composition artifacts require an explicit --target")
    order = {"solver": args.solver_order, "reduction": "length-lex-index-v1"}
    skeleton = (artifact_from_json(documents[0], pres, args.target, kind="skeleton", order=order)
                if "schema" in documents[0] else Skeleton.from_json(documents[0]))
    target = skeleton.target if args.target is None else args.target
    if target_identity(pres, target) != target_identity(pres, skeleton.target):
        raise ValueError("requested target differs from the reduction skeleton")
    fitted = (artifact_from_json(documents[1], pres, target, kind="lift", order=order)
              if "schema" in documents[1] else Lift.from_json(documents[1]))
    return fitted, skeleton


def run(args):
    """Compose and publish an exact identity using already parsed command arguments."""
    from ..cli import _presentation
    from .verify import check
    from ..provenance import capture_provenance
    try:
        producer = capture_provenance()
        paths = (args.skeleton, args.lift, args.multipliers)
        snapshots = [Path(path).read_bytes() for path in paths]
        inputs = [{"name": Path(path).name, "bytes": len(content),
                   "sha256": hashlib.sha256(content).hexdigest()}
                  for path, content in zip(paths, snapshots)]
        documents = [json.loads(content) for content in snapshots]
        pres = _presentation(args)
        lift, skeleton = _composition_inputs(args, pres, documents[:2])
        multipliers = multipliers_from_json(documents[2], skeleton, lift.atoms, pres)
        cofactor = compose(lift, skeleton, pres, multipliers=multipliers, verbose=not args.json)
        report = check(cofactor, pres, verbose=not args.json)
        if not report["pass"]:
            result = {"pass": False, "status": "false-identity", "exact": report,
                      "producer": producer, "inputs": inputs}
            print(json.dumps(result) if args.json else "Exact identity verification failed; output was not written.")
            return 1
        _save_atomic(cofactor, args.out, producer=producer, inputs=inputs)
        result = {"pass": True, "status": "verified", "exact": report, "output": args.out,
                  "producer": producer, "inputs": inputs}
        print(json.dumps(result) if args.json else "Verified exact identity; wrote " + args.out)
        return 0
    except (OSError, ValueError, TypeError, KeyError, ZeroDivisionError) as error:
        if args.json:
            print(json.dumps({"pass": False, "status": "invalid-input", "error": str(error)}))
        else:
            print(f"compose: {error}", file=sys.stderr)
        return 2


def main(argv=None):
    """Compatibility module command using the installed frontend's implementation."""
    import argparse
    parser = argparse.ArgumentParser(description="compose and exactly verify a Q(q,t) cofactor")
    add_arguments(parser)
    return run(parser.parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
