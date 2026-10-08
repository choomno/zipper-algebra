#!/usr/bin/env python3
"""Compose X_3 cofactor in Q(q,t) using the symbolically-lifted basis.

X_3 = q^3 * basis[17] - q^2 * basis[16] + q*(q-1) * basis[102]

basis[k] for k in transitive deps = sum c_{k,i}(q,t) * left * R_idx * right (or basis[j])

Recursively substitute basis references to get an expression purely over the
31 R-inputs, then collect by (left, idx, right) keys, producing a cofactor

   X_3 = sum_t  C_t(q,t) * left_t * R_{i_t} * right_t

where C_t(q,t) lives in Q(q,t). Verify at fresh (q,t) samples by evaluating
each C_t and comparing against X_3 specialized at the same (q,t).
"""

from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)  # tools/, where zipper_in_ms.py lives
for _p in (HERE, PARENT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import json
import time
from fractions import Fraction
from collections import defaultdict
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from extract_cofactor import parse_noex_proof, build_originals_and_x3
from zipper_in_ms import Poly, Coef
from relations_naming import label
# cofactor/ -> tools/ -> zipper-algebra/ . These scripts were written when they
# lived at the repo root; a later reorganisation moved them without updating ROOT.
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
JSON_PATH = os.path.join(ROOT, 'logs_and_traces', 'f4ncgb',
                         'basis-symbolic-Qqt.json')

VAR_NAMES = ['x1', 'x2', 'x3', 'y1', 'y2', 'y3']
VAR_INDEX = {v: i + 1 for i, v in enumerate(VAR_NAMES)}


def names_to_mon(names):
    return tuple(VAR_INDEX[n] for n in names)


# The denominator atoms travel WITH the lift (tools/cofactor/atoms.py): the lift JSON
# carries the atom list its exponent tuples are positional over. Older lift files, written
# before the atom set was data, mean the historical six.
from atoms import AtomSet  # noqa: E402
ATOMS = AtomSet.DEFAULT


def load_lifted():
    global ATOMS
    with open(JSON_PATH) as f:
        data = json.load(f)
    ATOMS = AtomSet.from_json(data['atoms']) if 'atoms' in data else AtomSet.DEFAULT
    print(f'  atoms: {ATOMS.names()}')
    lifted = {}
    for k_s, terms in data['lifted'].items():
        k = int(k_s)
        out = []
        for term in terms:
            if term is None:
                out.append(None); continue
            denom = tuple(term['denom'])
            num = {tuple(m): Fraction(c) for m, c in term['num']}
            out.append({
                'denom': denom,
                'num': num,
                'left': term['left'],
                'kind': term['kind'],
                'idx': term['idx'],
                'right': term['right']
            })
        lifted[k] = out
    return lifted, set(data['needed'])


# Represent a Q(q,t) rational as (numerator_poly_dict, denom_tuple).
# To add two rationals with different denoms, we cross-multiply.
# Denoms are products of (q,t,q-1,t-q,t+q^2,q+1) — we use exponent tuples.

def denom_eval(denom, q, t):
    return ATOMS.eval_denom(denom, q, t)


def denom_lcm(d1, d2):
    return AtomSet.lcm(d1, d2)


def denom_div_exps(big, small):
    """Return exponents of (big / small) — assumes big >= small componentwise."""
    return tuple(b - s for b, s in zip(big, small))


def poly_mul_factor(num, factor_exps):
    """Multiply a polynomial dict {(i,j): coef} by prod_k atom_k^{e_k}, expanding."""
    out = num
    for atom, e in zip(ATOMS.atoms, factor_exps):
        for _ in range(e):
            nxt = defaultdict(Fraction)
            for (mi, mj), c0 in out.items():
                for (ai, aj), ac in atom.coefs.items():
                    nxt[(mi + ai, mj + aj)] += c0 * ac
            out = dict(nxt)
    return {m: cc for m, cc in out.items() if cc != 0}


def poly_add(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, Fraction(0)) + v
        if out[k] == 0:
            del out[k]
    return out


def poly_mul(a, b):
    out = defaultdict(Fraction)
    for (i1, j1), c1 in a.items():
        for (i2, j2), c2 in b.items():
            out[(i1 + i2, j1 + j2)] += c1 * c2
    return {k: v for k, v in out.items() if v != 0}


def poly_scalar(a, c):
    if c == 0: return {}
    return {k: v * c for k, v in a.items()}


def rat_add(r1, r2):
    """Add two (num,denom) rationals; result still in form (num',denom_lcm).
    No GCD simplification — that's OK for correctness."""
    n1, d1 = r1
    n2, d2 = r2
    d_lcm = denom_lcm(d1, d2)
    # Multiply n1 by d_lcm / d1, n2 by d_lcm / d2
    mul1 = denom_div_exps(d_lcm, d1)
    mul2 = denom_div_exps(d_lcm, d2)
    n1s = poly_mul_factor(n1, mul1)
    n2s = poly_mul_factor(n2, mul2)
    return (poly_add(n1s, n2s), d_lcm)


def rat_mul(r1, r2):
    n1, d1 = r1
    n2, d2 = r2
    return (poly_mul(n1, n2), AtomSet.mul(d1, d2))


def rat_zero():
    return ({}, ATOMS.zero())


def rat_const(c):
    """Return Q(q,t) rational of constant Fraction c."""
    if c == 0: return rat_zero()
    return ({(0, 0): Fraction(c)}, ATOMS.zero())


def rat_eval(rat, q, t):
    num, denom = rat
    if not num: return Fraction(0)
    n_val = sum(c * Fraction(q ** mi * t ** mj) for (mi, mj), c in num.items())
    d_val = Fraction(denom_eval(denom, q, t))
    return n_val / d_val


def main():
    print('Loading lifted basis...')
    lifted, needed = load_lifted()
    print(f'  {len(lifted)} basis elements lifted, needed = {len(needed)}')

    # Compute the symbolic Q(q,t) cofactor of each lifted basis[k] over the
    # 31 R-inputs, iteratively (basis indices reference smaller indices).
    print('Recursive symbolic expansion of basis[k] over R-inputs...')
    # expanded[k] = dict keyed by (left_mon_tuple, R_idx, right_mon_tuple)
    # -> (num_poly, denom_exps) rational
    expanded = {}
    failed_terms = 0
    for k in sorted(needed):
        bucket = {}
        for term in lifted[k]:
            if term is None:
                failed_terms += 1
                continue
            denom_key = term['denom']
            num = term['num']
            l_mon = names_to_mon(term['left'])
            r_mon = names_to_mon(term['right'])
            coef_rat = (dict(num), denom_key)
            if term['kind'] == 'input':
                key = (l_mon, term['idx'], r_mon)
                bucket[key] = rat_add(bucket.get(key, rat_zero()), coef_rat)
            else:
                # Substitute basis[idx]'s already-expanded cofactor
                sub_idx = term['idx']
                if sub_idx not in expanded:
                    print(f'  ERROR: basis[{k}] references basis[{sub_idx}] '
                          f'but {sub_idx} not yet expanded')
                    failed_terms += 1
                    continue
                for skey, srat in expanded[sub_idx].items():
                    s_l, s_i, s_r = skey
                    new_left = l_mon + s_l
                    new_right = s_r + r_mon
                    new_rat = rat_mul(coef_rat, srat)
                    nk = (new_left, s_i, new_right)
                    bucket[nk] = rat_add(bucket.get(nk, rat_zero()), new_rat)
        # Drop zero entries
        bucket = {k_: v for k_, v in bucket.items() if v[0]}
        expanded[k] = bucket
        print(f'  basis[{k}]: {len(bucket)} R-input terms')
    print(f'Total failed lifted terms: {failed_terms}')

    if 16 not in expanded or 17 not in expanded or 102 not in expanded:
        print('Cannot compose X_3 cofactor: missing key basis elements.')
        return 1

    # Compose: X_3 = q^3 * basis[17] - q^2 * basis[16] + q*(q-1) * basis[102]
    print('\nComposing X_3 cofactor...')
    cofactor = {}
    for mult_rat, k in [
        (({(3, 0): Fraction(1)}, ATOMS.zero()), 17),     # +q^3
        (({(2, 0): Fraction(-1)}, ATOMS.zero()), 16),    # -q^2
        # +q*(q-1) = q^2 - q
        (({(2, 0): Fraction(1), (1, 0): Fraction(-1)}, ATOMS.zero()), 102),
    ]:
        for key, rat in expanded[k].items():
            new_rat = rat_mul(mult_rat, rat)
            if key in cofactor:
                cofactor[key] = rat_add(cofactor[key], new_rat)
            else:
                cofactor[key] = new_rat
    cofactor = {k: v for k, v in cofactor.items() if v[0]}
    print(f'X_3 cofactor: {len(cofactor)} (left, R_idx, right) terms')

    # Verify at several (q,t)
    print('\nVerifying X_3 cofactor at fresh (q,t)...')
    test_qts = [(5, 7), (3, 11), (7, 13), (11, 17), (10, 7)]
    for q, t in test_qts:
        try:
            originals, x3 = build_originals_and_x3(0, q, t)
        except Exception as e:
            print(f'  ({q},{t}): build error {e}'); continue
        acc = Poly(0)
        for (lmon, idx, rmon), rat in cofactor.items():
            cval = rat_eval(rat, q, t)
            cc = Coef(cval, 0)
            lp = Poly(0); lp.add_term(lmon, Coef.one(0))
            rp = Poly(0); rp.add_term(rmon, Coef.one(0))
            prod = (lp * originals[idx] * rp).scalar_mul(cc)
            acc = acc + prod
        diff = acc - x3
        ok = (not diff.terms)
        print(f'  ({q},{t}): {"OK" if ok else "MISMATCH"}'
              + (f' (diff terms = {len(diff.terms)})' if not ok else ''))

    # Write output
    out_path = os.path.join(ROOT, 'logs_and_traces', 'f4ncgb',
                            'cofactor-X3-n4-Qqt-symbolic.txt')
    with open(out_path, 'w') as f:
        f.write('# Symbolic Q(q,t) cofactor of X_3 over the n=4 Bigelow zipper algebra.\n')
        f.write('# Built by composing 3-step GB reduction (X_3 = q^3*B17 - q^2*B16 + q(q-1)*B102)\n')
        f.write(f'# with symbolic lift of {len(needed)} basis elements (transitive deps).\n')
        f.write(f'# Total cofactor terms: {len(cofactor)}\n#\n')
        # Sort
        items = sorted(cofactor.items(), key=lambda x: (x[0][1], len(x[0][0]), x[0][0], len(x[0][2]), x[0][2]))
        for (lmon, idx, rmon), (num, denom) in items:
            num_s = ' + '.join(f'{c}*q^{mi}*t^{mj}' for (mi, mj), c in sorted(num.items())) if num else '0'
            d_s = ATOMS.denom_str(denom)
            ls = '*'.join(VAR_NAMES[i - 1] for i in lmon) if lmon else '()'
            rs = '*'.join(VAR_NAMES[i - 1] for i in rmon) if rmon else '()'
            f.write(f'{ls} * {label(idx)} * {rs}  ::  ({num_s}) / ({d_s})\n')
    print(f'\nWrote: {out_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
