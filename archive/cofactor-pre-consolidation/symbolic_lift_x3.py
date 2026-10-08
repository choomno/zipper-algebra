#!/usr/bin/env python3
"""Symbolic lifting of the X_3 cofactor to Q(q,t).

Strategy (from user task):
  1. At anchor (q,t)=(5,7), X_3 reduces in 3 steps against the GB:
       X_3 = q^3 * basis[17] - q^2 * basis[16] + q*(q-1) * basis[102]
     (these multipliers verified across all 226 samples).
  2. basis[16], basis[17], basis[102] and their recursive non-expanded-proof
     dependencies form a closed set of 25 basis elements.
  3. The non-expanded proof skeleton (left, kind, idx, right) for these 25
     elements is identical across 217/226 (q,t) samples.
  4. For each (k, term i), we collect the rational coefficient at every
     stable sample (q,t). We fit it as N(q,t)/D(q,t) with denominator from
     a candidate list and numerator a polynomial in (q,t).
  5. Holding out random samples we cross-validate, and at the end we have a
     symbolic basis[k] = sum (rational(q,t)) * left * R_idx * right.
  6. Compose with X_3 = q^3 * basis[17] - q^2 * basis[16] + q*(q-1) * basis[102]
     to get the cofactor in Q(q,t) terms.
  7. Verify at fresh (q,t) samples.
"""

from __future__ import annotations
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py lives
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import sys, os, glob, re, json, time
from fractions import Fraction
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from extract_cofactor import parse_noex_proof


# cofactor/ -> tools/ -> zipper-algebra/ . These scripts were written when they
# lived at the repo root; a later reorganisation moved them without updating ROOT.
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PROOF_DIR = os.path.join(ROOT, 'logs_and_traces', 'f4ncgb')


def skel(terms):
    return tuple((tuple(left), kind, idx, tuple(right))
                 for _, left, kind, idx, right in terms)


def transitive_deps(anchor_main, seeds):
    needed = set()
    work = list(seeds)
    while work:
        k = work.pop()
        if k in needed:
            continue
        needed.add(k)
        for _, _, kind, idx, _ in anchor_main[k]:
            if kind == 'basis':
                work.append(idx)
    return needed


def load_stable_samples(anchor_main, anchor_skel, needed):
    qt_re = re.compile(r'proof-n4-Q-q(\d+)-t(\d+)-noex.txt')
    samples = []
    for pf in sorted(glob.glob(os.path.join(PROOF_DIR,
                                            'proof-n4-Q-q*-noex.txt'))):
        m = qt_re.search(pf)
        if not m:
            continue
        q, t = int(m.group(1)), int(m.group(2))
        pr = parse_noex_proof(pf)
        if len(pr) != len(anchor_skel):
            continue
        if any(skel(pr[k]) != anchor_skel[k] for k in needed):
            continue
        samples.append((q, t, pr))
    return samples


# The denominator atoms are DATA (tools/cofactor/atoms.py; prop:atoms). They are the
# hypothesis class of fit_rational, and a wrong class cannot be repaired by sampling --
# lift_until_verified.py failure mode (b). Loaded from the atom file written for this
# (presentation, target) if it exists, else the historical six as a prior; main() then
# INFERS from the sampled denominators and appends anything the prior lacks. Existing
# atoms are never removed or reordered, so every earlier lift/cofactor file stays valid.
from atoms import AtomSet, infer as infer_atoms  # noqa: E402

ATOMS_PATH = os.path.join(PROOF_DIR, 'atoms-n4-Q.json')
ATOMS = AtomSet.load(ATOMS_PATH) if os.path.exists(ATOMS_PATH) else AtomSet.DEFAULT

# Kept as names for the rest of this file; all delegate to ATOMS.
DENOM_FACTORS = [(a.name, a.eval) for a in ATOMS.atoms]


def poly_eval_qt(coefs, q, t):
    """Evaluate {(i,j): coef} as sum coef * q^i * t^j, exactly."""
    return sum(c * Fraction(q) ** i * Fraction(t) ** j for (i, j), c in coefs.items())


def denom_eval_from_exps(exps, q, t):
    return ATOMS.eval_denom(exps, q, t)


def gen_denoms_sorted():
    # Exponent tuples inside the atom set's box, ascending total degree. For the historical
    # six this is exactly the enumeration the pipeline always used (box (5,4,2,1,2,1)).
    return ATOMS.candidates()


def gen_denoms_for_value(val_q, val_t, factor_vals):
    """Given a Fraction value at one (q,t), enumerate denominator exponents
    that COULD work (denominator must divide a known multiple of factor product).
    For speed we just return DENOMS sorted by total degree."""
    return DENOMS


DENOMS = gen_denoms_sorted()


def fit_poly_qt(points, max_deg_q, max_deg_t):
    """Fit a polynomial p(q,t) of bidegree <= (max_deg_q, max_deg_t) to the
    samples. Returns the dict of coefficients (or None if infeasible)."""
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
    coefs = {m: Fraction(0) for m in basis_mons}
    for ri, col in enumerate(pivot_col):
        coefs[basis_mons[col]] = aug[ri][M]
    return coefs


def fit_rational(points, max_pq=5, max_pt=5):
    """Try (denom, numerator-degrees) until a perfect fit is found that
    matches ALL sample points (not just the training subset).

    Prune denominators whose `val * denom` is non-integer at ALL sample points
    quickly: we just check that scaled values clear to integers.
    """
    points = list(points)
    points.sort()
    for denom_exps in DENOMS:
        # Compute scaled at every sample; if at any sample the scaled value
        # has denominator != 1 (i.e. denominator doesn't actually clear),
        # this candidate cannot give a polynomial — skip immediately.
        scaled = []
        ok_scale = True
        for (q, t), val in points:
            d_val = denom_eval_from_exps(denom_exps, q, t)
            if d_val == 0:
                ok_scale = False; break
            v = val * Fraction(d_val)
            if v.denominator != 1:
                ok_scale = False; break
            scaled.append(((q, t), v))
        if not ok_scale:
            continue
        # Try bidegrees by ascending total degree
        for total_deg in range(max_pq + max_pt + 1):
            for dq in range(min(total_deg, max_pq) + 1):
                dt = total_deg - dq
                if dt > max_pt:
                    continue
                M = (dq + 1) * (dt + 1)
                if len(scaled) < M + 2:
                    continue
                fit_pts = scaled[:M + 2]
                res = fit_poly_qt(fit_pts, dq, dt)
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
                    return denom_exps, res
    return None


def poly_qt_str(coefs):
    """Pretty form of a poly given as {(i,j): coef}."""
    terms = []
    for (i, j), c in sorted(coefs.items()):
        if c == 0:
            continue
        s = ''
        if c > 0:
            s += '+'
        s += str(c.numerator) if c.denominator == 1 else f'{c.numerator}/{c.denominator}'
        if i:
            s += f'*q^{i}'
        if j:
            s += f'*t^{j}'
        terms.append(s)
    return ' '.join(terms) if terms else '0'


def denom_str(key):
    return ATOMS.denom_str(key)


def main():
    t_total = time.time()
    print('Loading anchor non-expanded proof...')
    anchor_main = parse_noex_proof(os.path.join(PROOF_DIR,
                                                 'proof-n4-Q-noex.txt'))
    anchor_skel = [skel(t) for t in anchor_main]
    # The skeleton is DATA, emitted by extract_cofactor.py (Phase 2), not constants.
    # These indices are a property of this target against this anchor and change
    # with either, so hard-coding them is the first thing that breaks for any
    # target other than X_3. Falls back to the historical values with a loud
    # warning if Phase 2 has not been run.
    _skel_path = os.path.join(PROOF_DIR, 'skeleton-n4-Q.json')
    if os.path.exists(_skel_path):
        with open(_skel_path) as _f:
            _skel = json.load(_f)
        _seeds = _skel['basis_indices']
        print(f"  skeleton from {os.path.basename(_skel_path)}: "
              f"indices {_seeds}, {_skel['n_steps']} steps, "
              f"remainder {_skel['remainder_terms']}")
        if _skel.get('remainder_terms'):
            raise SystemExit('REFUSING: skeleton has a nonzero remainder — the target '
                             'is not in the ideal at this anchor.')
    else:
        _seeds = [16, 17, 102]
        print('  WARNING: no skeleton-n4-Q.json; falling back to hard-coded '
              'indices [16, 17, 102]. Run extract_cofactor.py (Phase 2) first.')
    needed = transitive_deps(anchor_main, _seeds)
    print(f'  need to lift {len(needed)} basis elements')

    print('Loading stable samples...')
    samples = load_stable_samples(anchor_main, anchor_skel, needed)
    print(f'  {len(samples)} stable samples')

    # ATOM SET (planned feature, implemented 2026-09-17). Infer the denominator atoms from
    # the sampled coefficients themselves and extend the prior with anything it lacks.
    # Done once per (presentation, target) and persisted next to the skeleton; re-run with
    # INFER_ATOMS=1. This is search; the proof is Phase 5 on the lift built from it.
    global ATOMS, DENOMS, DENOM_FACTORS
    if not os.path.exists(ATOMS_PATH) or os.environ.get('INFER_ATOMS') == '1':
        print('Inferring the atom set from the sampled denominators...')
        slots = {(k, i): [((q, t), pr[k][i][0]) for q, t, pr in samples]
                 for k in sorted(needed) for i in range(len(anchor_main[k]))}
        ATOMS, atom_report = infer_atoms(slots, seed=ATOMS)
        ATOMS.save(ATOMS_PATH)
        with open(ATOMS_PATH[:-5] + '-report.json', 'w') as _f:
            json.dump(atom_report, _f, indent=1)
        print(f'  atoms: {ATOMS.names()}  box {ATOMS.max_exp}  -> {os.path.basename(ATOMS_PATH)}')
        if atom_report['unexplained']:
            print(f'  WARNING: {atom_report["unexplained"]} reconstructed denominators have a '
                  f'factor outside the search box {atom_report["search"]}; the lift below may '
                  f'fail with WRONG_CLASS. Widen the box in atoms.infer().')
    else:
        print(f'  atoms from {os.path.basename(ATOMS_PATH)}: {ATOMS.names()}  box {ATOMS.max_exp}')
    DENOM_FACTORS = [(a.name, a.eval) for a in ATOMS.atoms]
    DENOMS = gen_denoms_sorted()

    # DYNAMIC GRID SIZING (2026-09-16). A fixed grid is what produced the wrong
    # cofactor: the abandoned "unified" model was fit on run_lift_grid.py's
    # hard-coded 37 points and disagreed with the exact answer on all 695 slots.
    # fit_rational() requires agreement on ALL supplied points, but with fewer
    # points than unknowns many polynomials pass through all of them -- fitting is
    # not verification. So we hold points OUT, fit on the rest, and require the fit
    # to PREDICT the held-out values. If any coefficient fails, the grid is too
    # small and the answer is "sample more", not "ship it".
    HOLDOUT = max(4, len(samples) // 10)
    if len(samples) <= HOLDOUT + 2:
        print(f'NEED_MORE_SAMPLES: {len(samples)} samples cannot support a held-out check '
              f'(need > {HOLDOUT + 2}).')
        raise SystemExit(3)      # 3 == grid too small; see lift_until_verified.py
    samples_fit, samples_hold = samples[:-HOLDOUT], samples[-HOLDOUT:]
    print(f'  fitting on {len(samples_fit)}, holding out {len(samples_hold)} for verification')

    # For each (k, term i), fit a rational coef(q,t)
    lifted = {}
    failed = []
    for k in sorted(needed):
        terms = anchor_main[k]
        lifted_terms = []
        t0 = time.time()
        for i, (_coef0, left, kind, idx, right) in enumerate(terms):
            pts = [((q, t), pr[k][i][0]) for q, t, pr in samples_fit]
            res = fit_rational(pts)
            if res is None:
                lifted_terms.append(None)
                failed.append((k, i))
                continue
            denom_key, num_coefs = res
            # held-out verification: the fit must PREDICT points it never saw
            for (hq, ht, hpr) in samples_hold:
                want = hpr[k][i][0]
                d_val = denom_eval_from_exps(denom_key, hq, ht)
                if d_val == 0:
                    continue
                got = poly_eval_qt(num_coefs, hq, ht) / Fraction(d_val)
                if got != want:
                    lifted_terms.append(None)
                    failed.append((k, i, 'holdout'))
                    break
            else:
                lifted_terms.append((denom_key, num_coefs, left, kind, idx, right))
            continue
            lifted_terms.append((denom_key, num_coefs, left, kind, idx, right))
        lifted[k] = lifted_terms
        n_ok = sum(1 for x in lifted_terms if x is not None)
        print(f'  basis[{k}]: {n_ok}/{len(terms)} terms fit '
              f'({time.time()-t0:.1f}s)', flush=True)

    print(f'\nTotal fit failures (k, i): {len(failed)}')
    if failed:
        print('  ', failed[:20])
        held = [f for f in failed if len(f) > 2]      # mispredicted held-out points
        nofit = [f for f in failed if len(f) == 2]    # no model found at all
        if os.environ.get('FORCE_LIFT') != '1':
            # ANY unfitted coefficient means the lift is incomplete. Writing it
            # produces a cofactor missing a term, which composes and then fails the
            # symbolic check far downstream -- after clobbering verified artifacts.
            # (Observed 2026-09-16: basis[68] term 3 unfitted at 60 points, written
            # anyway, because the guard only checked held-out failures.)
            if held:
                print(f'WRONG_CLASS: {len(held)} coefficient(s) fit every point shown and still '
                      f'mispredict unseen ones. More points cannot fix this -- widen the '
                      f'hypothesis class (fit_rational max_pq/max_pt, or the denominator atoms).')
                raise SystemExit(4)  # 4 == hypothesis class too narrow
            print(f'NEED_MORE_SAMPLES: {len(nofit)} coefficient(s) had no model in the current '
                  f'class at {len(samples)} points. Refusing to write an incomplete lift.')
            raise SystemExit(3)      # 3 == grid too small; see lift_until_verified.py

    # Write symbolic basis representation
    out_path = os.path.join(PROOF_DIR, 'basis-symbolic-Qqt.txt')
    with open(out_path, 'w') as f:
        f.write('# Symbolic Q(q,t) lift of basis[k] for k in transitive deps of {16,17,102}.\n')
        f.write(f'# Anchor: proof-n4-Q-noex.txt; stable samples used: {len(samples)}\n')
        f.write(f'# Atoms (denominator factors, positional): {ATOMS.names()}\n')
        f.write(f'# Format per line: k = ..., term i: (numerator) / (denom) * left * [kind:idx] * right\n#\n')
        for k in sorted(needed):
            f.write(f'\n## basis[{k}]\n')
            for i, item in enumerate(lifted[k]):
                if item is None:
                    f.write(f'  term {i}: FAILED\n')
                    continue
                denom_key, num, left, kind, idx, right = item
                num_s = poly_qt_str(num)
                den_s = denom_str(denom_key)
                ls = '*'.join(left) if left else '()'
                rs = '*'.join(right) if right else '()'
                f.write(f'  term {i}: ({num_s}) / ({den_s}) * {ls} * [{kind}:{idx}] * {rs}\n')
    print(f'Wrote symbolic basis: {out_path}')

    # Save machine-readable too
    out_json = os.path.join(PROOF_DIR, 'basis-symbolic-Qqt.json')
    data = {
        'needed': sorted(needed),
        'atoms': ATOMS.to_json(),
        'lifted': {}
    }
    for k in sorted(needed):
        l = []
        for item in lifted[k]:
            if item is None:
                l.append(None); continue
            denom_key, num, left, kind, idx, right = item
            l.append({
                'denom': list(denom_key),
                'num': [[list(m), str(c)] for m, c in num.items() if c != 0],
                'left': left, 'kind': kind, 'idx': idx, 'right': right
            })
        data['lifted'][str(k)] = l
    with open(out_json, 'w') as f:
        json.dump(data, f, indent=1)
    print(f'Wrote JSON: {out_json}')

    print(f'\nTotal time: {time.time()-t_total:.1f}s')
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
