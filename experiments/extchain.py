"""extchain.py -- exact extension-chain search, equivalent to (k+1,d) existence.

Theory (Mobius decomposition in the last coordinate; proved in research notes):

  Any Boolean f on {0,1}^{k+1} splits as  f(x,y) = g(x) + y*h(x)  with
  g = f(.,0),  h = f(.,1) - f(.,0).  If deg f <= d then deg h <= d-1
  (the top part of y*h would otherwise have degree > d).

  Putting u := 1 - h  (integer-valued, u(x) in {0,1,2}, u(0) = 0) yields the
  exact equivalence

     (k+1,d) exists  <=>  exists g: (k,d)-function (fully sensitive at 0,
     deg <= d)  and  u: {0,1}^k -> {0,1,2}  with
        g <= u <= g + 1  pointwise,   u(0) = 0,   deg u <= d-1.

  Reconstruction:  f(x,y) = g(x) + y*(1 - u(x))  is Boolean, deg <= d and
  fully sensitive at 0 (f(0)=0, f(e_i)=1, f(e_new)=1).

  Equivalently (r := u - g, a Boolean function): exists g and Boolean r with
  r(0) = 0 and deg(g + r) <= d-1.

Usage:
  python extchain.py <k> <d> <tl> <workers> [--tag t] [--gfile <json>] [--log]
    no --gfile    : full model, g searched too (with exact level bounds)
    --gfile f.json: g fixed (y-values of a found (k,d) solution, 2..d-sets)

SAT => writes the (k+1,d) y-values to experiments/ksol_e<k+1>_<d><tag>.json
after independent verification (ksearch.verify).
"""
import json
import os
import sys
import time
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ksearch
import ksearch2

from ortools.sat.python import cp_model


def g_values_from_y(k, d, yvals):
    """Full cube values g(F), F subset [k], from y-json (2..d-sets)."""
    y = {}
    for key, v in yvals.items():
        y[tuple(sorted(int(x) for x in str(key).split(",")))] = v
    c = {}
    for m in range(0, d + 1):
        for T in combinations(range(k), m):
            if m == 0:
                c[T] = 0
            elif m == 1:
                c[T] = 1
            else:
                val = 0
                for j in range(0, m + 1):
                    for U in combinations(T, j):
                        base = 0 if j == 0 else (1 if j == 1 else y[U])
                        val += (-1) ** (m - j) * base
                c[T] = val
    g = {}
    for m in range(0, k + 1):
        for F in combinations(range(k), m):
            g[F] = sum(c[T] for j in range(0, min(m, d) + 1)
                       for T in combinations(F, j))
    return g


def u_degree_constraints(model, u, k, d):
    """deg u <= d-1: vanishing Mobius coefficients at all |S| >= d."""
    for m in range(d, k + 1):
        for S in combinations(range(k), m):
            terms = []
            for j in range(1, m + 1):
                cf = 1 if (m - j) % 2 == 0 else -1
                for T in combinations(S, j):
                    terms.append(cf * u[T])
            model.Add(sum(terms) == 0)


def build_fixed(k, d, g, log=False):
    model = cp_model.CpModel()
    u = {}
    for m in range(0, k + 1):
        for T in combinations(range(k), m):
            gv = int(g[T])
            if m == 0:
                u[T] = model.NewIntVar(0, 0, "u_")
            else:
                u[T] = model.NewIntVar(gv, gv + 1, "u_" + "_".join(map(str, T)))
    u_degree_constraints(model, u, k, d)
    return model, u


def build_full(k, d, log=False):
    bounds, verts = ksearch2.level_bounds(k, d, verbose=log)
    model, y, fvars = ksearch2.build2(k, d, bounds)
    u = {}
    gexpr = {}
    for m in range(0, k + 1):
        for T in combinations(range(k), m):
            u[T] = model.NewIntVar(0, 2, "u_" + "_".join(map(str, T)))
            if m == 0:
                gexpr[T] = 0
            elif m == 1:
                gexpr[T] = 1
            elif m <= d:
                gexpr[T] = y[T]
            else:
                gexpr[T] = fvars[T]
            model.Add(u[T] >= gexpr[T])
            model.Add(u[T] <= gexpr[T] + 1)
    u_degree_constraints(model, u, k, d)
    model.Add(u[()] == 0)
    return model, u, (y, fvars, None)


def extract(k, d, g, uvals):
    """y-values of f on [k+1] from g on [k] and u on [k]."""
    ysol = {}
    for m in range(2, d + 1):
        for B in combinations(range(k + 1), m):
            if k not in B:
                ysol[B] = g[B]
            else:
                T = tuple(B[:-1]) if B[-1] == k else tuple(x for x in B if x != k)
                ysol[B] = g[T] + 1 - uvals[T]
    return ysol


if __name__ == "__main__":
    args = list(sys.argv[1:])
    gfile = None
    if "--gfile" in args:
        i = args.index("--gfile")
        gfile = args[i + 1]
        del args[i:i + 2]
    if "--hintg" in args:
        i = args.index("--hintg")
        hint_path = args[i + 1]
        del args[i:i + 2]
    else:
        hint_path = None
    flags = set(a for a in args if a.startswith("--"))
    pos = [a for a in args if not a.startswith("--")]
    k = int(pos[0])
    d = int(pos[1])
    tl = float(pos[2])
    wk = int(pos[3])
    tag = pos[4] if len(pos) > 4 else ""
    log = "--log" in flags

    t0 = time.time()
    if gfile:
        with open(gfile) as fh:
            yraw = json.load(fh)
        g = g_values_from_y(k, d, yraw)
        model, u = build_fixed(k, d, g)
        print("FIXED-g model k=%d d=%d: %d u-vars (%.1fs)"
              % (k, d, len(u), time.time() - t0), flush=True)
    else:
        model, u, gvars = build_full(k, d, log)
        print("FULL model k=%d d=%d: %d u-vars (%.1fs)"
              % (k, d, len(u), time.time() - t0), flush=True)

    solver = cp_model.CpSolver()
    if "--branch" in flags:
        free_u = [u[T] for m in range(1, d) for T in combinations(range(k), m)]
        rest_u = [u[T] for m in range(d, k + 1) for T in combinations(range(k), m)]
        model.add_decision_strategy(free_u + rest_u,
                                   cp_model.CHOOSE_FIRST,
                                   cp_model.SELECT_MIN_VALUE)
        solver.parameters.search_branching = cp_model.FIXED_SEARCH
    solver.parameters.max_time_in_seconds = tl
    solver.parameters.num_search_workers = wk
    if log:
        solver.parameters.log_search_progress = True
    if "--lns" in flags:
        solver.parameters.use_lns_only = True
        solver.parameters.randomize_search = True
    if hint_path:
        raw = json.load(open(hint_path))
        d_old = max(len(str(key).split(",")) for key in raw)
        gvals = g_values_from_y(k, d_old, raw)
        if gfile:
            for T, v in u.items():
                model.add_hint(v, gvals[T])
        else:
            y_t, f_t = gvars[0], gvars[1]
            for B, v in y_t.items():
                model.add_hint(v, gvals[B])
            for A, v in f_t.items():
                model.add_hint(v, gvals[A])
            for T, v in u.items():
                model.add_hint(v, gvals[T])
        print("hint applied from %s (d_old=%d)" % (hint_path, d_old), flush=True)
    st = solver.Solve(model)
    name = solver.StatusName(st)
    print("k=%d d=%d status=%s wall=%.1f (extchain%s)"
          % (k, d, name, solver.WallTime(), " fixed" if gfile else ""),
          flush=True)
    if name in ("OPTIMAL", "FEASIBLE"):
        uvals = {T: solver.Value(v) for T, v in u.items()}
        if gfile:
            ysol = extract(k, d, g, uvals)
        else:
            # rebuild g from the model solution (y-vars on 2..d-sets)
            y_t = gvars[0]
            ysol_tmp = {B: solver.Value(v) for B, v in y_t.items()}
            g = g_values_from_y(k, d,
                                {",".join(map(str, sorted(B))): v
                                 for B, v in ysol_tmp.items()})
            ysol = extract(k, d, g, uvals)
        ok, err = ksearch.verify(k + 1, d, ysol)
        print("independent verification:", "OK" if ok else err, flush=True)
        fn = "experiments/ksol_e%d_%d%s.json" % (k + 1, d, tag)
        with open(fn, "w") as fh:
            json.dump({",".join(map(str, sorted(B))): v
                       for B, v in ysol.items()}, fh)
        print("saved", fn, flush=True)
