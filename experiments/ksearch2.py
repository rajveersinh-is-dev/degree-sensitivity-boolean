"""ksearch2.py -- exact (n,d) search with exact level-count bounds.

Identical model to ksearch.py (y_B = f(B) for 2 <= |B| <= d; derived f(A) for
|A| > d; deg <= d  <=>  every derived value in {0,1}), plus implied constraints
extracted from the exact level polytope:

  * exact integer bounds on every level count  b_m = #{x: |x|=m, f(x)=1},
    computed by exact vertex enumeration of the level system
        w_m = m + sum_{k=2..d} a_k C(m,k) in [0,1]   (m = 2..n),
        a_2 in [-2,-1], a_3 in [0,4], a_4 in [-14,9], a_5 in [-80,86].
    Every one of these rows is a separately valid necessary condition for a
    fully-sensitive-at-0 degree-<=d Boolean function, hence the resulting
    b_m-bounds are implied by the model (pure strengthening).

  * optional WLOG pinning y_{01} = 1 when the exact bound gives lo_2 >= 1
    (valid because we may relabel coordinates).

Usage: python ksearch2.py <n> <d> <time_limit> <tag> <workers> [--fix01] [--log]
"""
import json
import os
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations
from math import ceil, comb, floor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ksearch
import amanotype

from ortools.sat.python import cp_model


def level_bounds(n, d, verbose=True):
    """Exact integer bounds on b_m from the level polytope (vertex enum)."""
    rows = amanotype.base_rows(n, d, ranges=True)
    nv = d - 1
    t0 = time.time()
    for i in range(nv):
        for sgn in (1, -1):
            c = [Fr(0)] * nv
            c[i] = Fr(sgn)
            st, val, _ = amanotype.linprog_dir(rows, nv, c)
            if st != 0:
                raise RuntimeError("level polytope unbounded: a_%d dir %+d" % (i + 2, sgn))
    verts = amanotype.enumerate_vertices(rows, nv)
    if not verts:
        raise RuntimeError("no vertices found (system infeasible?)")
    bounds = {}
    for m in range(2, n + 1):
        vals = [amanotype.val_of_func(v, amanotype.wvec(n, d, m)) for v in verts]
        lo, hi = min(vals), max(vals)
        Cm = comb(n, m)
        bounds[m] = (ceil(lo * Cm), floor(hi * Cm))
    if verbose:
        print("level polytope: %d vertices (%.1fs); exact integer bounds:"
              % (len(verts), time.time() - t0))
        for m in range(2, n + 1):
            lo, hi = bounds[m]
            print("  b_%-2d in [%d, %d] of %d%s"
                  % (m, lo, hi, comb(n, m), "  FORCED" if lo == hi else ""))
    return bounds, verts


def build2(n, d, bounds, fix01=False, slices=None):
    model, y, fvars, big_As = ksearch.build(n, d)
    ys = {}
    for B, v in y.items():
        ys.setdefault(len(B), []).append(v)
    for m in range(2, d + 1):
        lo, hi = bounds[m]
        model.Add(sum(ys[m]) >= lo)
        model.Add(sum(ys[m]) <= hi)
    fs = {}
    for A, v in fvars.items():
        fs.setdefault(len(A), []).append(v)
    for m in range(d + 1, n + 1):
        lo, hi = bounds[m]
        model.Add(sum(fs[m]) >= lo)
        model.Add(sum(fs[m]) <= hi)
    if fix01:
        if bounds[2][0] < 1:
            print("WARNING: lo_2 = 0; cannot pin y_01 = 1")
        else:
            model.Add(y[(0, 1)] == 1)
    if slices:
        # per-slice level-count bounds: every m0-subset S carries a restriction
        # g_S which is again fully sensitive at 0 of degree <= d, so its level
        # counts are inside the exact (m0,d) polytope range.
        lo_m0, hi_m0 = slices
        fv_map = dict(fvars)
        for m0 in range(max(lo_m0, 2), hi_m0 + 1):
            if m0 < d:
                continue
            sb, _v = level_bounds(m0, d, verbose=False)
            for S in combinations(range(n), m0):
                for m in range(2, m0 + 1):
                    l, h = sb[m]
                    if m <= d:
                        terms = [y[T] for T in combinations(S, m)]
                    else:
                        terms = [fv_map[T] for T in combinations(S, m)]
                    model.Add(sum(terms) >= l)
                    model.Add(sum(terms) <= h)
    return model, y, fvars


if __name__ == "__main__":
    args = list(sys.argv[1:])
    slices = None
    if "--slices" in args:
        i = args.index("--slices")
        slices = (int(args[i + 1]), int(args[i + 2]))
        del args[i:i + 3]
    flags = set(a for a in args if a.startswith("--"))
    pos = [a for a in args if not a.startswith("--")]
    n = int(pos[0])
    d = int(pos[1])
    tl = float(pos[2])
    tag = pos[3]
    wk = int(pos[4])
    bounds, verts = level_bounds(n, d)
    os.makedirs("experiments", exist_ok=True)
    with open("experiments/bounds_%d_%d.json" % (n, d), "w") as fh:
        json.dump({"n": n, "d": d, "n_vertices": len(verts),
                   "bounds": {str(m): list(bounds[m]) for m in bounds}}, fh,
                  indent=1)
    model, y, fvars = build2(n, d, bounds, fix01=("--fix01" in flags),
                             slices=slices)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = tl
    solver.parameters.num_search_workers = wk
    if "--log" in flags:
        solver.parameters.log_search_progress = True
    st = solver.Solve(model)
    name = solver.StatusName(st)
    print("n=%d d=%d status=%s wall=%.1f (ksearch2%s%s)"
          % (n, d, name, solver.WallTime(),
             " fix01" if "--fix01" in flags else "",
             " slices=%d-%d" % slices if slices else ""))
    if name in ("OPTIMAL", "FEASIBLE"):
        sol = {B: solver.Value(v) for B, v in y.items()}
        ok, err = ksearch.verify(n, d, sol)
        print("independent verification:", "OK" if ok else err)
        with open("experiments/ksol2_%d_%d%s.json" % (n, d, tag), "w") as fh:
            json.dump({",".join(map(str, sorted(B))): v for B, v in sol.items()},
                      fh)
        print("saved experiments/ksol2_%d_%d%s.json" % (n, d, tag))
    sys.stdout.flush()
