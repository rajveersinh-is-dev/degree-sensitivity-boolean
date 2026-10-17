"""localsearch.py -- local search (SA) for the fixed-g extension problem.

Given a (k,d)-function g (y-json), search r (free bits, levels 1..d-1, r(0)=0)
minimizing the penalty  P = sum_{S: d<=|S|<=k} dist(r_derived(S), {0,1})
where r_derived(S) = -c_g(S) - sum_{T subsetneq S} (-1)^{|S|-|T|} r(T)
(r(T) = r_free[T] if |T|<=d-1 else r_derived[T]).

P == 0  <=>  r is a valid co-factor  <=>  (k+1,d) exists (with this g).

On success: builds f, verifies, saves.

Usage: python localsearch.py <k> <d> <gfile> [--time TL] [--restarts R] [--seed S]
"""
import json
import os
import random
import sys
import time
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extchain
import ksearch

import numpy as np


def setup(k, d, g):
    # c_g(S) for |S| >= d
    c_g = {}
    for m in range(d, k + 1):
        for S in combinations(range(k), m):
            cg = 0
            for j in range(0, m + 1):
                for T in combinations(S, j):
                    cg += (-1) ** (m - j) * g[T]
            c_g[S] = cg
    # target t_S = -c_g(S)
    # free bits: T with 1 <= |T| <= d-1
    free = []
    free_idx = {}
    for m in range(1, d):
        for T in combinations(range(k), m):
            free_idx[T] = len(free)
            free.append(T)
    nfree = len(free)
    # derived S: d <= |S| <= k
    derived = []
    for m in range(d, k + 1):
        for S in combinations(range(k), m):
            derived.append(S)
    # for each derived S: list of (free_bit_index, coeff) for T subset S, 1<=|T|<=d-1
    # coeff in r(S) contribution: (−1)^{|S|-|T|}; R(S) changes by -coeff when r(T) flips
    # precompute affected: for each free bit, list of derived indices containing it
    # and per-derived: list of (free_idx, coeff)
    S_free = []   # list of list of (fi, coeff) for each derived
    for S in derived:
        lst = []
        for j in range(1, min(d, len(S)) + 1):
            if j >= d:
                break
            for T in combinations(S, j):
                if T not in free_idx:
                    continue
                coeff = 1 if (len(S) - j) % 2 == 0 else -1
                lst.append((free_idx[T], coeff))
        S_free.append(lst)
    affected = [[] for _ in range(nfree)]
    for si, S in enumerate(derived):
        for (fi, coeff) in S_free[si]:
            affected[fi].append((si, coeff))
    # targets
    targets = [-c_g[S] for S in derived]
    return free, free_idx, derived, S_free, affected, targets, nfree, len(derived)


def cost_of(R):
    # R: int; dist to {0,1}
    if R == 0 or R == 1:
        return 0
    if R < 0:
        return -R
    return R - 1


def init_state(nfree, nderv, S_free, targets, affected):
    r = [random.randrange(2) for _ in range(nfree)]
    # r[0] for T=empty is 0 (not in free list)
    R = list(targets)  # R(S) when all free bits = 0
    # apply r flips
    for fi in range(nfree):
        if r[fi] == 1:
            for (si, coeff) in affected[fi]:
                # flipping free bit 0->1: R[si] -= coeff
                R[si] -= coeff
    P = sum(cost_of(R[si]) for si in range(nderv))
    return r, R, P


def main():
    k = int(sys.argv[1])
    d = int(sys.argv[2])
    gfile = sys.argv[3]
    tl = 120.0
    restarts = 5
    seed = 0
    i = 4
    while i < len(sys.argv):
        if sys.argv[i] == "--time":
            tl = float(sys.argv[i + 1]); i += 2; continue
        if sys.argv[i] == "--restarts":
            restarts = int(sys.argv[i + 1]); i += 2; continue
        if sys.argv[i] == "--seed":
            seed = int(sys.argv[i + 1]); i += 2; continue
        i += 1
    random.seed(seed)

    with open(gfile) as fh:
        yraw = json.load(fh)
    g = extchain.g_values_from_y(k, d, yraw)
    # verify g fully sensitive
    if g[()] != 0:
        print("g(0) != 0", flush=True); return
    for i in range(k):
        if g[(i,)] != 1:
            print("g(e_%d) != 1" % i, flush=True); return

    free, free_idx, derived, S_free, affected, targets, nfree, nderv = setup(k, d, g)
    print("free bits=%d derived=%d penalty0-target" % (nfree, nderv), flush=True)

    best_overall = None
    best_P_overall = 10**9
    t_start = time.time()
    ri = 0
    while time.time() - t_start < tl:
        ri += 1
        if ri == 1:
            r = [0] * nfree
            R = list(targets)
            for fi in range(nfree):
                if r[fi] == 1:
                    for (si, coeff) in affected[fi]:
                        R[si] -= coeff
            P = sum(cost_of(R[si]) for si in range(nderv))
        else:
            r, R, P = init_state(nfree, nderv, S_free, targets, affected)
        best_r = list(r)
        best_P = P
        no_improve = 0
        temp = 2.0
        steps = 0
        while time.time() - t_start < tl and no_improve < 2000 and P > 0:
            steps += 1
            # pick a random free bit
            fi = random.randrange(nfree)
            # compute delta P if flip
            # affected derived indices for fi
            aff = affected[fi]
            # old contribution
            old_cost = 0
            new_cost = 0
            # compute old R and new R for affected
            deltas = []
            for (si, coeff) in aff:
                # flipping r[fi] 0->1: R[si] -= coeff; 1->0: R[si] += coeff
                if r[fi] == 0:
                    dR = -coeff
                else:
                    dR = coeff
                old = R[si]
                new = old + dR
                old_cost += cost_of(old)
                new_cost += cost_of(new)
                deltas.append((si, dR))
            dP = new_cost - old_cost
            accept = (dP <= 0) or (random.random() < (2.71828 ** (-dP / max(temp, 1e-6))))
            if accept:
                # apply flip
                r[fi] ^= 1
                for (si, dR) in deltas:
                    R[si] += dR
                P += dP
                if P < best_P:
                    best_r = list(r)
                    best_P = P
                    no_improve = 0
                else:
                    no_improve += 1
            else:
                no_improve += 1
            temp *= 0.9999
            if temp < 0.05:
                temp = 0.05
        if P == 0:
            # found solution
            print("FOUND solution in restart %d (steps=%d)" % (ri, steps), flush=True)
            # build u and f
            r_free_map = {free[i]: best_r[i] for i in range(nfree)}
            r_full = {}
            for m in range(0, 1):
                for T in combinations(range(k), m):
                    r_full[T] = 0
            for m in range(1, d):
                for T in combinations(range(k), m):
                    r_full[T] = r_free_map[T]
            for m in range(d, k + 1):
                for S in combinations(range(k), m):
                    val = -c_g_  # placeholder
            # compute derived r via recurrence
            r_full2 = dict(r_full)
            for m in range(d, k + 1):
                for S in combinations(range(k), m):
                    val = -c_g[S if False else S]
                    for j in range(0, m):
                        for T in combinations(S, j):
                            if j == 0:
                                val -= (-1) ** (m - j) * 0
                            else:
                                val -= (-1) ** (m - j) * r_full2[T]
                    r_full2[S] = val
            # verify derived in {0,1}
            bad = [S for S in derived if r_full2[S] not in (0, 1)]
            if bad:
                print("post-hoc derived not in {0,1}: %d" % len(bad), flush=True)
                continue
            uvals = {}
            for m in range(0, k + 1):
                for T in combinations(range(k), m):
                    uvals[T] = g[T] + r_full2[T]
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
            if ok:
                fn = "experiments/ksol_ls_%d_%d.json" % (k + 1, d)
                with open(fn, "w") as fh:
                    json.dump({",".join(map(str, sorted(B))): v for B, v in ysol.items()}, fh)
                print("saved", fn, flush=True)
                return
        if best_P < best_P_overall:
            best_P_overall = best_P
            best_overall = list(best_r)
        print("restart %d: best_P=%d (overall best=%d) steps=%d"
              % (ri, best_P, best_P_overall, steps), flush=True)
    print("TIMEOUT. overall best penalty=%d" % best_P_overall, flush=True)
    # save best r for inspection
    if best_overall is not None:
        with open("experiments/lsbest_%d_%d.json" % (k, d), "w") as fh:
            json.dump({"free": [list(T) for T in free],
                       "r": best_overall, "penalty": best_P_overall}, fh)
        print("saved best r to lsbest_%d_%d.json" % (k, d), flush=True)


if __name__ == "__main__":
    main()
