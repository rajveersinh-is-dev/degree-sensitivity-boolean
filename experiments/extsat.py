"""extsat.py -- SAT-based fixed-g extension search (pysat + Cadical).

Given a (k,d)-function g (y-json), search r (free bits, levels 1..d-1) such that
the derived r(S) in {0,1} for all d <= |S| <= k, where
    r(S) = -c_g(S) - sum_{T subsetneq S, 1<=|T|<=d-1} (-1)^{|S|-|T|} r(T).
Equivalently (k+1,d) exists iff such r exists.

SAT => builds and verifies the (k+1,d) function, saves to experiments.

Usage: python extsat.py <k> <d> <gfile> [--time TL] [--tag t]
"""
from itertools import combinations
import json
import os
import sys
import time

from pysat.card import CardEnc, EncType
from pysat.formula import CNF
from pysat.solvers import Cadical153
import extchain
import ksearch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))



def solve(k, d, gfile, time_limit=300, tag=""):
    """Solve.
    
    Args:
        k:
        d:
        gfile:
        time_limit (int):
        tag (str):
    
    Returns:
        The computed result
    
    """
    with open(gfile) as fh:
        yraw = json.load(fh)
    g = extchain.g_values_from_y(k, d, yraw)

    # c_g(S) for |S| >= d (0 for |S| > d since deg g <= d)
    c_g = {}
    for m in range(d, k + 1):
        for S in combinations(range(k), m):
            cg = 0
            for j in range(0, m + 1):
                for T in combinations(S, j):
                    cg += (-1) ** (m - j) * g[T]
            c_g[S] = cg

    # free bits: T with 1 <= |T| <= d-1
    free = []
    free_idx = {}
    for m in range(1, d):
        for T in combinations(range(k), m):
            free_idx[T] = len(free)
            free.append(T)

    # derived S: d <= |S| <= k
    derived = []
    for m in range(d, k + 1):
        for S in combinations(range(k), m):
            derived.append(S)

    # precompute for each S: lits (signed literals) and lo/hi bounds
    # R(S) = -c_g(S) - sum_{T subsetneq S} coeff*r(T) in {0,1}
    # coeff = (-1)^{|S|-|T|}
    # L = sum_{coeff +1} r(T) + sum_{coeff -1} (1 - r(T)) in {target-1+C, target+C}
    # where target = -c_g(S), C = #negative coeff free bits in S
    cnf = CNF()
    top_id = len(free)
    for S in derived:
        target = -c_g[S]
        lits = []
        C = 0
        for j in range(1, min(d, len(S)) + 1):
            if j >= d:
                break
            for T in combinations(S, j):
                if T not in free_idx:
                    continue
                coeff = 1 if (len(S) - j) % 2 == 0 else -1
                if coeff == 1:
                    lits.append(free_idx[T] + 1)
                else:
                    lits.append(-(free_idx[T] + 1))
                    C += 1
        lo = target - 1 + C
        hi = target + C
        if lo > 0:
            cnf.extend(CardEnc.atleast(lits, lo, top_id=top_id,
                                       encoding=EncType.seqcounter))
            top_id = cnf.nv
        if hi < len(lits):
            cnf.extend(CardEnc.atmost(lits, hi, top_id=top_id,
                                      encoding=EncType.seqcounter))
            top_id = cnf.nv

    # solve
    solver = Cadical153()
    t0 = time.time()
    solver.append_formula(cnf)
    print(f"CNF: {len(cnf.clauses)} clauses, {cnf.nv} vars; solving...", flush=True)
    sat = solver.solve()
    elapsed = time.time() - t0
    print("SAT=%s elapsed=%.2f clauses=%d free=%d derived=%d"
          % (sat, elapsed, len(cnf.clauses), len(free), len(derived)), flush=True)
    if not sat:
        solver.delete()
        return None
    model = solver.get_model()
    solver.delete()

    # r values for free bits
    r_free = {free[i]: (1 if model[i] > 0 else 0) for i in range(len(free))}
    r_full = {}
    for m in range(0, 0 + 1):
        for T in combinations(range(k), m):
            r_full[T] = 0
    for m in range(1, d):
        for T in combinations(range(k), m):
            r_full[T] = r_free[T]
    # derived r(S) for |S| >= d, in increasing size order
    for m in range(d, k + 1):
        for S in combinations(range(k), m):
            val = -c_g[S]
            for j in range(1, m):
                for T in combinations(S, j):
                    coeff = (-1) ** (m - j)
                    val -= coeff * r_full[T]
            r_full[S] = val

    bad = [S for S in derived if r_full[S] not in (0, 1)]
    if bad:
        print(f"WARNING: {len(bad)} derived r not in {{0,1}}: {...}", flush=True)
        return None

    # u = g + r (full cube values)
    uvals = {}
    for m in range(0, k + 1):
        for T in combinations(range(k), m):
            uvals[T] = g[T] + r_full[T]

    # y-values of f on [k+1]
    ysol = {}
    for m in range(2, d + 1):
        for B in combinations(range(k + 1), m):
            if k not in B:
                ysol[B] = g[B]
            else:
                T = tuple(x for x in B if x != k)
                ysol[B] = g[T] + 1 - uvals[T]

    ok, err = ksearch.verify(k + 1, d, ysol)
    print("verify:", "OK" if ok else err, flush=True)
    if not ok:
        return None
    fn = f"experiments/ksol_sat_{k + 1}_{d}{tag}.json"
    with open(fn, "w") as fh:
        json.dump({",".join(map(str, sorted(B))): v for B, v in ysol.items()}, fh)
    print("saved", fn, flush=True)
    return ysol


if __name__ == "__main__":
    k = int(sys.argv[1])
    d = int(sys.argv[2])
    gfile = sys.argv[3]
    tl = 300
    tag = ""
    i = 4
    while i < len(sys.argv):
        if sys.argv[i] == "--time":
            tl = float(sys.argv[i + 1]); i += 2; continue
        if sys.argv[i] == "--tag":
            tag = sys.argv[i + 1]; i += 2; continue
        i += 1
    # (time limit not enforced by Cadical in this wrapper; rely on small instances)
    ysol = solve(k, d, gfile, time_limit=tl, tag=tag)
    print("done")
