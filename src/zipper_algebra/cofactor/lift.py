#!/usr/bin/env python3
"""Lift the basis elements a skeleton needs from per-point proofs to Q(q,t).

Consolidates the fitting half of symbolic_lift_x3.py (SPEC.md, lift.py).  The pieces:

  fit_poly_qt / fit_rational   exact full-rank fitting: for each candidate
                               denominator in AtomSet.candidates() order (ascending exponent
                               sum), scale the sampled values, and take the first bidegree in
                               ascending total degree whose determined numerator agrees on
                               every fitting point. Same enumeration order
                               as before, so the first-found denominator is the same.
  fit_rational_known           the same scaling + numerator fit for ONE candidate, plus a
                               reducedness check; used to try a reconstructed denominator
                               before the full enumeration (the SPEC optimisation).
  line_denominator             heuristic denominator exponent proposal from univariate
                               reconstruction along grid lines. Specialization can merge
                               factors, so this proposal need not bound the generic exponents.
  stable_points                the stability filter by LM re-identification: a point is stable
                               for a key if its basis has an element with that leading word and
                               the element's proof has the anchor's term structure
                               (left, kind, key-or-input-index, right).  Indices are never
                               compared, so a renumbered basis is still stable.
  lift                         the held-out protocol of symbolic_lift_x3.main: HOLDOUT =
                               max(4, n // 10) points are never fitted and every fit must predict
                               them; rejected fits raise NeedMoreSamples, with WrongClass
                               retained as a compatibility subtype for held-out rejection.

Words are whatever the Proof carries (tuples of 1-based variable indices per SPEC; the legacy
parser's variable-name lists are converted by the caller); lift never interprets them, nor
the key text, which is owned by proof.Basis.

Measured on the 100-point X_3 grid on disk (tests_lift.py, 231 slots, atoms-n4-Q.json box
[6,4,2,1,2,1] = 1260 candidates): the enumeration route takes about 1.4 s in total, the
line-reconstruction route about 100 s, both giving byte-identical entries.  The optimisation
pays only when the candidate count is large (it grows as the product of the box widths), so
`lift` uses it only when asked (`known=`) -- see the report in notes/refactor/lift-report.md.
"""
from __future__ import annotations

import json
import time
from fractions import Fraction

from .atoms import (AtomSet, infer as infer_atoms, reconstruct_univariate,
                    exact_fraction, normalise_points)


Anchor = tuple  # (q, t, char); char 0 = Q


# ----------------------------------------------------------------------------
# result types (were exit codes 3 / 4 of symbolic_lift_x3.py)
# ----------------------------------------------------------------------------

class NeedMoreSamples(Exception):
    """Reconstruction is inconclusive for the current samples and search class.

    More independent samples or a revised class may help. A candidate fit still requires
    exact certificate verification; this exception makes no non-membership claim.
    """

    def __init__(self, msg, slots=()):
        super().__init__(msg)
        self.slots = list(slots)


class WrongClass(NeedMoreSamples):
    """Compatibility name for a rejected held-out fit, an inconclusive reconstruction.

    More independent data or a revised search class can help. This exception does not
    establish that the entire hypothesis class is wrong. New callers catch NeedMoreSamples.
    """

    def __init__(self, msg, slots=()):
        super().__init__(msg)
        self.slots = list(slots)


# ----------------------------------------------------------------------------
# the fitter (symbolic_lift_x3.py, semantics unchanged)
# ----------------------------------------------------------------------------

def poly_eval_qt(coefs, q, t):
    """Evaluate {(i,j): coef} as sum coef * q^i * t^j, exactly."""
    return sum(c * Fraction(q) ** i * Fraction(t) ** j for (i, j), c in coefs.items())


def fit_poly_qt(points, max_deg_q, max_deg_t):
    """Fit a polynomial p(q,t) of bidegree <= (max_deg_q, max_deg_t) to the
    samples. Returns a uniquely determined coefficient dictionary, or None.
    Duplicate locations count once; conflicting or inexact samples are rejected."""
    if any(type(d) is not int or d < 0 for d in (max_deg_q, max_deg_t)):
        raise ValueError("degree bounds must be nonnegative integers")
    points = normalise_points(points, 2)
    basis_mons = [(i, j) for i in range(max_deg_q + 1)
                  for j in range(max_deg_t + 1)]
    M = len(basis_mons)
    if len(points) < M:
        return None
    # Augmented matrix, exact Fractions
    aug = []
    for (q, t), val in points:
        row = [Fraction(q ** i * t ** j) for (i, j) in basis_mons]
        row.append(Fraction(val))
        aug.append(row)
    n = len(aug)
    rank = 0
    pivot_col = []
    for col in range(M):
        pv = None
        for r in range(rank, n):
            if aug[r][col] != 0:
                pv = r
                break
        if pv is None:
            continue
        aug[rank], aug[pv] = aug[pv], aug[rank]
        inv = Fraction(1, 1) / aug[rank][col]
        aug[rank] = [x * inv for x in aug[rank]]
        for r in range(n):
            if r != rank and aug[r][col] != 0:
                f = aug[r][col]
                aug[r] = [aug[r][i] - f * aug[rank][i] for i in range(M + 1)]
        pivot_col.append(col)
        rank += 1
    # Check consistency
    for i in range(n):
        if all(aug[i][c] == 0 for c in range(M)) and aug[i][M] != 0:
            return None
    if rank != M:
        return None
    coefs = {m: Fraction(0) for m in basis_mons}
    for ri, col in enumerate(pivot_col):
        coefs[basis_mons[col]] = aug[ri][M]
    return coefs


def _scale_by_denominator(points, atoms, denom_exps):
    """Scale exact rational samples by a candidate denominator, rejecting its poles.

    Polynomial numerators have coefficients in Q; scaled values need not be integers.
    """
    scaled = []
    for (q, t), val in points:
        d_val = atoms.eval_denom(denom_exps, q, t)
        if d_val == 0:
            return None
        v = val * Fraction(d_val)
        scaled.append(((q, t), v))
    return scaled


def _fit_numerator(scaled, max_pq, max_pt):
    """Try bidegrees by total degree, using every sample for rank and consistency.

    At least two samples beyond the coefficient count are required. A degenerate prefix
    must not hide independent rows later in the data. Returns a coefficient dict or None.
    """
    for total_deg in range(max_pq + max_pt + 1):
        for dq in range(min(total_deg, max_pq) + 1):
            dt = total_deg - dq
            if dt > max_pt:
                continue
            M = (dq + 1) * (dt + 1)
            if len(scaled) < M + 2:
                continue
            res = fit_poly_qt(scaled, dq, dt)
            if res is None:
                continue
            ok = True
            for (q, t), val in scaled:
                pred = sum(res[(i, j)] * Fraction(q ** i * t ** j)
                           for (i, j) in res)
                if pred != val:
                    ok = False
                    break
            if ok:
                return res
    return None


def fit_rational(points, atoms: AtomSet, max_pq=5, max_pt=5, candidates=None):
    """Try (denom, numerator-degrees) until a perfect fit is found that matches ALL
    sample points.  `candidates` defaults to atoms.candidates(); pass it in when fitting
    many slots against one atom set.  Returns (denom_exps, num_coefs) or None."""
    if any(type(d) is not int or d < 0 for d in (max_pq, max_pt)):
        raise ValueError("degree bounds must be nonnegative integers")
    atoms.validate()
    points = normalise_points(points, 2)
    if candidates is None:
        candidates = atoms.candidates()
    for denom_exps in candidates:
        scaled = _scale_by_denominator(points, atoms, denom_exps)
        if scaled is None:
            continue
        res = _fit_numerator(scaled, max_pq, max_pt)
        if res is not None:
            return denom_exps, res
    return None


def _divisible(num: dict, atom: dict) -> bool:
    """Whether the polynomial {(i,j): c} `atom` divides `num` in Q[q,t].  Division by one
    polynomial under lex (q before t): a single polynomial is a Groebner basis of the ideal
    it generates, so the remainder is zero iff it divides.  A leading term the divisor's
    leading term does not divide would go to the remainder, so we stop there."""
    r = {m: Fraction(c) for m, c in num.items() if c != 0}
    if not r:
        return True
    lm_a = max(atom)
    c_a = Fraction(atom[lm_a])
    while r:
        lm_r = max(r)
        if lm_r[0] < lm_a[0] or lm_r[1] < lm_a[1]:
            return False
        shift = (lm_r[0] - lm_a[0], lm_r[1] - lm_a[1])
        f = r[lm_r] / c_a
        for (i, j), c in atom.items():
            m = (i + shift[0], j + shift[1])
            v = r.get(m, Fraction(0)) - f * c
            if v == 0:
                r.pop(m, None)
            else:
                r[m] = v
    return True


def is_reduced(atoms: AtomSet, denom_exps, num: dict) -> bool:
    """Compatibility predicate: no whole used atom divides the numerator.

    This is not a general polynomial gcd check. Composite atoms can share a proper factor
    with the numerator. With irreducible pairwise nonassociate atoms, the predicate does
    imply coprimality of the represented numerator and denominator.
    """
    return not any(e > 0 and _divisible(num, a.coefs) for a, e in zip(atoms.atoms, denom_exps))


def fit_rational_known(points, atoms: AtomSet, denom_exps, max_pq=5, max_pt=5):
    """Try one bounded denominator, with no removable whole declared atom.

    Return (denom_exps, num_coefs), or None. As for fit_rational, the fitted numerator is
    determined by the exact samples in its selected degree box. This shortcut need not
    choose the same representation as enumeration. For composite atoms, a proper common
    factor can remain; is_reduced checks only cancellation of whole declared atoms.
    """
    if any(type(d) is not int or d < 0 for d in (max_pq, max_pt)):
        raise ValueError("degree bounds must be nonnegative integers")
    atoms.validate()
    points = normalise_points(points, 2)
    denom_exps = tuple(denom_exps)
    atoms.validate_exponents(denom_exps)
    if any(e > m for e, m in zip(denom_exps, atoms.max_exp)):
        return None                     # outside the box: the enumeration could never find it
    scaled = _scale_by_denominator(points, atoms, denom_exps)
    if scaled is None:
        return None
    res = _fit_numerator(scaled, max_pq, max_pt)
    if res is None or not is_reduced(atoms, denom_exps, res):
        return None
    return denom_exps, res


# ----------------------------------------------------------------------------
# per-slot denominator from line reconstruction (the SPEC optimisation's input)
# ----------------------------------------------------------------------------

def _specialise(coefs: dict, fixed: str, x0) -> list:
    """{(i,j): c} restricted to t = x0 ('t') or q = x0 ('q'), as a univariate coefficient
    list low->high in the other variable.  (atoms._specialise, kept private there.)"""
    out = {}
    for (i, j), c in coefs.items():
        if fixed == "t":
            out[i] = out.get(i, 0) + c * Fraction(x0) ** j
        else:
            out[j] = out.get(j, 0) + c * Fraction(x0) ** i
    deg = max(out) if out else 0
    poly = [out.get(k, Fraction(0)) for k in range(deg + 1)]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def _poly_divmod(num, den):
    """Exact univariate division over Fraction, coefficient lists low->high.  (atoms._poly_divmod.)"""
    num = [Fraction(c) for c in num]
    den = [Fraction(c) for c in den]
    while den and den[-1] == 0:
        den.pop()
    if not den:
        raise ZeroDivisionError
    q = [Fraction(0)] * max(0, len(num) - len(den) + 1)
    r = num[:]
    while len(r) >= len(den) and any(r):
        if r[-1] == 0:
            r.pop(); continue
        k = len(r) - len(den)
        c = r[-1] / den[-1]
        q[k] = c
        for i, d in enumerate(den):
            r[k + i] -= c * d
        r.pop()
    while r and r[-1] == 0:
        r.pop()
    return q, r


def line_denominator(points, atoms: AtomSet, min_line=6, max_lines=None):
    """Propose denominator exponents from reduced rational reconstructions on grid lines.

    On sufficiently long fixed-q or fixed-t lines, successively divide the reconstructed
    denominator by specialized atoms and retain maximal observed multiplicities. Factors
    can collide or cancel under specialization, so this is a heuristic proposal, not a
    componentwise lower bound on generic exponents. fit_rational_known validates a proposed
    representation; enumeration is the fallback. Return None when no line reconstructs.
    max_lines limits the number of longest lines considered in each direction.
    """
    points = normalise_points(points, 2)
    by_t, by_q = {}, {}
    for (q, t), v in points:
        by_t.setdefault(t, []).append((q, v))
        by_q.setdefault(q, []).append((t, v))
    exps = [0] * len(atoms)
    seen = False
    for fixed, groups in (("t", by_t), ("q", by_q)):
        lines = [(x0, line) for x0, line in groups.items() if len(line) >= min_line]
        lines.sort(key=lambda x: len(x[1]), reverse=True)
        if max_lines is not None:
            lines = lines[:max_lines]
        for x0, line in lines:
            rec = reconstruct_univariate(line)
            if rec is None:
                continue
            seen = True
            D = rec[1]
            for k, a in enumerate(atoms.atoms):
                a1 = _specialise(a.coefs, fixed, x0)
                if len(a1) < 2:
                    continue
                m = 0
                while len(D) >= len(a1):
                    qq, r = _poly_divmod(D, a1)
                    if r:
                        break
                    D, m = qq, m + 1
                exps[k] = max(exps[k], m)
    return tuple(exps) if seen else None


# ----------------------------------------------------------------------------
# stability by LM re-identification, and the closure of a skeleton's keys
# ----------------------------------------------------------------------------

def term_structure(terms, basis):
    """(left, kind, key-or-input-index, right) per term: a basis reference is named by the
    referenced element's key, so the structure is comparable across points."""
    return tuple((tuple(left), kind, basis.key[idx] if kind == "basis" else idx, tuple(right))
                 for _, left, kind, idx, right in terms)


def closure(seed_keys, proof, basis):
    """Keys of every basis element the seeds transitively reference in `proof`, in
    ascending anchor-index order (symbolic_lift_x3.transitive_deps, by key)."""
    needed = set()
    work = [basis.index_by_key[k] for k in seed_keys]
    while work:
        k = work.pop()
        if k in needed:
            continue
        needed.add(k)
        for _, _, kind, idx, _ in proof.terms[k]:
            if kind == "basis":
                work.append(idx)
    return [basis.key[k] for k in sorted(needed)]


def stable_points(keys, anchor_proof, anchor_basis, samples):
    """The points of `samples` (in samples.points order) whose proofs re-identify every key
    with the anchor's term structure.  `samples` provides points, proof(q,t), basis(q,t)."""
    want = {key: term_structure(anchor_proof.terms[anchor_basis.index_by_key[key]], anchor_basis)
            for key in keys}
    out = []
    for q, t in samples.points:
        basis = samples.basis(q, t)
        proof = samples.proof(q, t)
        for key in keys:
            j = basis.index_by_key.get(key)
            if j is None or j >= len(proof.terms) \
                    or term_structure(proof.terms[j], basis) != want[key]:
                break
        else:
            out.append((q, t))
    return out


def holdout_size(n_samples: int) -> int:
    return max(4, n_samples // 10)


# ----------------------------------------------------------------------------
# the result
# ----------------------------------------------------------------------------

class Lift:
    """entries[key] = [(den_exps, num, left, kind, key_or_idx, right), ...] in the anchor
    proof's term order; num is {(dq,dt): Fraction} without zero coefficients.  anchor_index
    maps each key to its index at the anchor, for traceability only."""

    def __init__(self, atoms: AtomSet, entries: dict, samples_used: int, holdout: int,
                 anchor, anchor_index: dict, atom_report=None, timing=None):
        self.atoms = atoms
        self.entries = entries
        self.samples_used = samples_used
        self.holdout = holdout
        self.anchor = tuple(anchor)
        self.anchor_index = dict(anchor_index)
        self.atom_report = atom_report
        self.timing = timing or {}

    def keys(self):
        return list(self.entries)

    def to_json(self):
        entries = {}
        for key, terms in self.entries.items():
            out = []
            for den, num, left, kind, ref, right in terms:
                item = {"denom": list(den),
                        "num": [[list(m), str(c)] for m, c in num.items() if c != 0],
                        "left": list(left), "kind": kind, "ref": ref, "right": list(right)}
                if kind == "basis":
                    item["ref_anchor_index"] = self.anchor_index.get(ref)
                out.append(item)
            entries[key] = {"anchor_index": self.anchor_index.get(key), "terms": out}
        return {"atoms": self.atoms.to_json(), "anchor": list(self.anchor),
                "samples_used": self.samples_used, "holdout": self.holdout,
                "atom_report": self.atom_report, "timing": self.timing, "entries": entries}

    @classmethod
    def from_json(cls, d):
        atoms = AtomSet.from_json(d["atoms"])
        entries, anchor_index = {}, {}
        for key, ent in d["entries"].items():
            anchor_index[key] = ent["anchor_index"]
            terms = []
            for it in ent["terms"]:
                num = {tuple(m): Fraction(c) for m, c in it["num"]}
                terms.append((tuple(it["denom"]), num, tuple(it["left"]), it["kind"],
                              it["ref"], tuple(it["right"])))
            entries[key] = terms
        return cls(atoms, entries, d["samples_used"], d["holdout"], tuple(d["anchor"]),
                   anchor_index, d.get("atom_report"), d.get("timing"))

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.to_json(), f, indent=1)

    @classmethod
    def load(cls, path):
        with open(path) as f:
            return cls.from_json(json.load(f))

    def __eq__(self, other):
        return isinstance(other, Lift) and self.entries == other.entries \
            and self.atoms.names() == other.atoms.names()


# ----------------------------------------------------------------------------
# the lift
# ----------------------------------------------------------------------------

def lift(skeleton, samples, atoms: AtomSet | None = None, holdout: int | None = None,
         known=None, max_pq=5, max_pt=5, verbose=True) -> Lift:
    """Lift every basis element `skeleton` needs to Q(q,t) from the stable points of `samples`.

    skeleton : has .anchor (q, t, char) and .steps [{basis_key, ...}], optionally
               .basis_keys (the closure; checked against the anchor proof if present).
    samples  : has .points [(q,t)], .proof(q,t), .basis(q,t); the anchor point must be one
               of them.  Order matters: the stable points keep samples.points order and the
               LAST `holdout` of them are held out (symbolic_lift_x3: glob order, last
               HOLDOUT).
    atoms    : the hypothesis class.  None = atoms.infer over the stable samples' values,
               seeded with AtomSet.DEFAULT, exactly as symbolic_lift_x3.main does; the
               report is kept in Lift.atom_report.
    known    : the SPEC optimisation.  "lines" computes each slot's denominator by
               line_denominator() and tries it before the enumeration; a dict
               {(key, i): exps} supplies denominators from elsewhere; None (default)
               enumerates only. The shortcut can select a different valid candidate; it
               must pass the same held-out and final identity checks.

    Raises NeedMoreSamples (including its legacy WrongClass subtype) with rejected
    slots. Held-out rejection is inconclusive; more data or a revised class may help.
    """
    if holdout is not None and (type(holdout) is not int or holdout < 1):
        raise ValueError("holdout must be a positive integer")
    if (len(skeleton.anchor) != 3 or type(skeleton.anchor[2]) is not int
            or skeleton.anchor[2] != 0):
        raise ValueError("rational reconstruction requires characteristic-zero samples")
    sample_char = getattr(samples, "char", 0)
    if type(sample_char) is not int or sample_char != 0:
        raise ValueError("rational reconstruction requires characteristic-zero samples")
    if any(type(d) is not int or d < 0 for d in (max_pq, max_pt)):
        raise ValueError("degree bounds must be nonnegative integers")
    normalise_points([(skeleton.anchor[:2], 0)], 2)
    t_total = time.time()
    q0, t0 = skeleton.anchor[0], skeleton.anchor[1]
    anchor_proof, anchor_basis = samples.proof(q0, t0), samples.basis(q0, t0)
    seeds = [s["basis_key"] for s in skeleton.steps]
    keys = closure(seeds, anchor_proof, anchor_basis)
    given = getattr(skeleton, "basis_keys", None)
    if given and set(given) != set(keys):
        raise ValueError("skeleton.basis_keys is not the closure of its steps under the anchor "
                         f"proof at ({q0},{t0}): {sorted(set(given) ^ set(keys))}")
    anchor_index = {k: anchor_basis.index_by_key[k] for k in keys}
    if verbose:
        print(f"  need to lift {len(keys)} basis elements")

    # Prefer the SampleSet's own filter: it compares each key against the MODAL shape across
    # points rather than against the anchor's, so one degenerate anchor coefficient costs a
    # point instead of the whole grid (samples.shape_census). Falls back to the anchor-based
    # criterion for a samples object that does not provide it.
    if hasattr(samples, "stable"):
        stable = samples.stable(skeleton)
    else:
        stable = stable_points(keys, anchor_proof, anchor_basis, samples)
    stable = list(dict.fromkeys(normalise_points([(point, 0)], 2)[0][0]
                                for point in stable))
    if verbose:
        print(f"  {len(stable)} stable samples of {len(samples.points)}")
    # values per slot at every stable point, then the proofs can go
    n_terms = {k: len(anchor_proof.terms[anchor_index[k]]) for k in keys}
    values = {(k, i): [] for k in keys for i in range(n_terms[k])}
    for q, t in stable:
        proof, basis = samples.proof(q, t), samples.basis(q, t)
        for k in keys:
            terms = proof.terms[basis.index_by_key[k]]
            for i in range(n_terms[k]):
                values[(k, i)].append(((q, t), exact_fraction(terms[i][0])))

    n = len(stable)
    h = holdout_size(n) if holdout is None else holdout
    if n <= h + 2:
        raise NeedMoreSamples(f"{n} samples cannot support a held-out check (need > {h + 2})")

    atom_report = None
    if atoms is None:
        if verbose:
            print("Inferring the atom set from the sampled denominators...")
        atoms, atom_report = infer_atoms(
            {slot: points[:n - h] for slot, points in values.items()},
            seed=AtomSet.DEFAULT, verbose=verbose)
        if verbose:
            print(f"  atoms: {atoms.names()}  box {atoms.max_exp}")
            if atom_report["unexplained"]:
                print(f"  WARNING: {atom_report['unexplained']} reconstructed denominators have "
                      f"a factor outside the search box {atom_report['search']}; the lift below "
                      f"may need more data or a revised search class in atoms.infer().")

    n_fit = n - h
    if verbose:
        print(f"  fitting on {n_fit}, holding out {h} for verification")

    candidates = atoms.candidates()
    entries, failed = {}, []
    routes = {"known": 0, "enum": 0}
    t_known = 0.0
    for k in keys:
        t_k = time.time()
        terms_out = []
        for i, (_c0, left, kind, idx, right) in enumerate(anchor_proof.terms[anchor_index[k]]):
            pts = values[(k, i)]
            pts_fit, pts_hold = pts[:n_fit], pts[n_fit:]
            res = None
            if known is not None:
                t_g = time.time()
                exps = line_denominator(pts_fit, atoms) if known == "lines" else known.get((k, i))
                if exps is not None:
                    res = fit_rational_known(pts_fit, atoms, exps, max_pq, max_pt)
                t_known += time.time() - t_g
                if res is not None:
                    routes["known"] += 1
            if res is None:
                res = fit_rational(pts_fit, atoms, max_pq, max_pt, candidates)
                if res is not None:
                    routes["enum"] += 1
            if res is None:
                terms_out.append(None)
                failed.append((k, i))
                continue
            den, num = res
            for (hq, ht), want in pts_hold:
                d_val = atoms.eval_denom(den, hq, ht)
                if d_val == 0 or poly_eval_qt(num, hq, ht) / Fraction(d_val) != want:
                    terms_out.append(None)
                    failed.append((k, i, "holdout"))
                    break
            else:
                ref = anchor_basis.key[idx] if kind == "basis" else idx
                terms_out.append((tuple(den), {m: c for m, c in num.items() if c != 0},
                                  tuple(left), kind, ref, tuple(right)))
        entries[k] = terms_out
        if verbose:
            n_ok = sum(1 for x in terms_out if x is not None)
            print(f"  basis[{anchor_index[k]}] {k}: {n_ok}/{n_terms[k]} terms fit "
                  f"({time.time() - t_k:.2f}s)", flush=True)

    timing = {"seconds": round(time.time() - t_total, 3), "routes": routes,
              "known_route_seconds": round(t_known, 3), "n_slots": len(values)}
    if failed:
        held = [f for f in failed if len(f) > 2]
        nofit = [f for f in failed if len(f) == 2]
        if held:
            raise WrongClass(f"{len(held)} coefficient(s) failed held-out validation; "
                             f"reconstruction is inconclusive and needs more independent "
                             f"data or a revised search class", held)
        raise NeedMoreSamples(f"{len(nofit)} coefficient(s) had no model in the current class "
                              f"at {n} points", nofit)
    return Lift(atoms, entries, n, h, tuple(skeleton.anchor), anchor_index, atom_report, timing)
