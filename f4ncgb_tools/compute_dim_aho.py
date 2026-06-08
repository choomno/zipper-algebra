#!/usr/bin/env python3
"""
Efficient Hilbert series / dimension computation for finitely-presented free
associative algebras, from a noncommutative Gröbner basis.

The set of standard monomials (free monomials not containing any leading
monomial as a contiguous substring) forms a regular language. We build an
Aho--Corasick automaton for the leading monomial set, then iteratively count
length-d paths through the non-rejecting subautomaton. This replaces
compute_dim.py's BFS-over-standard-monomials with a matrix iteration: each
degree costs O(|states| * |alphabet|) regardless of dim.

Usage:
    python3 compute_dim_aho.py path/to/basis.txt --max-deg 50
"""
from __future__ import annotations
import argparse
import re
from collections import deque


def parse_leading_monomials(path):
    """Parse f4ncgb basis output; return list of leading monomials, each as
    a tuple of variable names."""
    lms = []
    with open(path) as f:
        for raw in f:
            line = raw.strip().rstrip(",")
            if not line:
                continue
            # Leading term is everything before the first ' + ' or ' - ' separator.
            best_idx = -1
            for sep in (" + ", " - "):
                idx = line.find(sep, 1)
                if idx > 0 and (best_idx == -1 or idx < best_idx):
                    best_idx = idx
            if best_idx > 0:
                line = line[:best_idx]
            m = re.match(r"^([0-9]+/[0-9]+|[0-9]+)\*?", line)
            if m and (line[m.end():].lstrip("*") or "").strip() != "":
                mon = line[m.end():].lstrip("*")
            else:
                mon = line
            mon = mon.strip()
            if mon == "" or mon.isdigit():
                return [()]
            vs = tuple(v.strip() for v in mon.split("*") if v.strip())
            if vs:
                lms.append(vs)
    return lms


def collect_variables(lms):
    """Variables that actually appear in some leading monomial."""
    s = set()
    for lm in lms:
        s.update(lm)
    return sorted(s)


def build_aho_corasick(lms, variables):
    """Build the Aho-Corasick automaton for the pattern set lms over the
    alphabet `variables`. Returns:
        delta : dict[(state, var)] -> next_state
        is_match : list[bool] indexed by state
        nstates : int
        root : int
    `is_match[s]` is True iff arriving at state s means a leading monomial
    has just been matched as a substring. We will exclude these states (and
    their descendants) from the path-counting.
    """
    # Trie nodes are integers; root = 0.
    children = [dict()]   # children[s][v] = next_state (only direct trie edges)
    is_match = [False]    # is_match[s] = True if some leading monomial ends here
    # Build trie
    for lm in lms:
        s = 0
        for v in lm:
            if v not in children[s]:
                children.append(dict())
                is_match.append(False)
                new = len(children) - 1
                children[s][v] = new
            s = children[s][v]
        is_match[s] = True
    nstates = len(children)
    # Build failure function & goto function via BFS (standard AC construction).
    fail = [0] * nstates
    delta = {}
    # Root's children: failure = root; delta from root on v: child if exists else root.
    q = deque()
    for v in variables:
        if v in children[0]:
            c = children[0][v]
            fail[c] = 0
            delta[(0, v)] = c
            q.append(c)
        else:
            delta[(0, v)] = 0
    while q:
        s = q.popleft()
        for v in variables:
            if v in children[s]:
                c = children[s][v]
                # failure of c: follow fail[s] on v
                fs = fail[s]
                while fs != 0 and v not in children[fs]:
                    fs = fail[fs]
                if v in children[fs] and children[fs][v] != c:
                    fail[c] = children[fs][v]
                elif v in children[fs] and children[fs][v] == c:
                    # Same as parent: fail goes back to root (shouldn't happen here)
                    fail[c] = 0
                else:
                    fail[c] = 0
                # Propagate match: if any state on the failure chain is a match,
                # arriving at c also implies we have matched a leading monomial.
                if is_match[fail[c]]:
                    is_match[c] = True
                delta[(s, v)] = c
                q.append(c)
            else:
                # No direct child; follow failure
                delta[(s, v)] = delta.get((fail[s], v), 0)
    return delta, is_match, nstates


def count_hilbert(delta, is_match, nstates, variables, max_deg):
    """Iteratively count standard monomials of each degree via the AC automaton.
    counts[d] = number of length-d paths from root (state 0) that never visit
    a match state.
    """
    counts = [0] * (max_deg + 1)
    counts[0] = 1
    # f[s] = number of length-d paths to state s avoiding match states
    f = [0] * nstates
    f[0] = 1
    for d in range(1, max_deg + 1):
        nxt = [0] * nstates
        for s in range(nstates):
            if f[s] == 0 or is_match[s]:
                continue
            for v in variables:
                t = delta[(s, v)]
                if is_match[t]:
                    continue
                nxt[t] += f[s]
        f = nxt
        counts[d] = sum(f)
        if counts[d] == 0 and all(x == 0 for x in f):
            break
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("basis", help="path to f4ncgb basis output")
    ap.add_argument("--max-deg", type=int, default=40,
                    help="upper bound on monomial degree to enumerate")
    ap.add_argument("--vars", help="comma-separated variable names override")
    args = ap.parse_args()

    lms = parse_leading_monomials(args.basis)
    print(f"Read {len(lms)} leading monomials.")
    if args.vars:
        variables = args.vars.split(",")
    else:
        variables = collect_variables(lms)
    print(f"Variables (count {len(variables)}): {variables}")

    delta, is_match, nstates = build_aho_corasick(lms, variables)
    print(f"AC automaton: {nstates} states, "
          f"{sum(is_match)} match states excluded.")

    counts = count_hilbert(delta, is_match, nstates, variables, args.max_deg)
    total = sum(counts)
    # Find last nonzero
    last = max((d for d, c in enumerate(counts) if c != 0), default=0)
    print(f"Quotient dimension (up to deg {args.max_deg}): {total}")
    print(f"Hilbert series (deg 0..{args.max_deg}): {counts}")
    print(f"  (last nonzero: deg {last}, count {counts[last]})")
    if last == args.max_deg and counts[last] > 0:
        print("  WARNING: last nonzero is at the truncation bound;")
        print("           algebra may be infinite-dim. Try higher --max-deg.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
