"""Exact search for (n,d): does there exist a Boolean f on {0,1}^n of degree <= d
with f(0)=0 and f(e_i)=1 for all i (all n directions sensitive at the origin)?

Clean formulation (derived in research notes):
Let y_B = f(B) for subsets B with 2 <= |B| <= d (Booleans; y_empty = 0, y_{e_i} = 1
are fixed).  Then for every A with |A| > d:

    f(A) = sum_{B subset A, |B| <= d} (Mobius-transform of y at B)
         = sum_{B subset A, 2<=|B|<=d} K(A,B) y_B  +  sum_{i in A} K(A,{i}),
    K(A,B) = (-1)^{d-|B|} C(|A|-|B|-1, d-|B|),

and the requirement "deg f <= d" is equivalent to f(A) in {0,1} for all |A| > d,
which (given y is {0,1} on small sets) is exactly: the above expression in {0,1}.

We add f(A) as integer variables in [0,1] and post the equality.
Independent verification recomputes everything from scratch (see verify()).
"""
import json
import sys
from itertools import combinations
from math import comb

from ortools.sat.python import cp_model


def build(n, d):
    model = cp_model.CpModel()
    y = {}
    for k in range(2, d + 1):
        for B in combinations(range(n), k):
            y[B] = model.NewBoolVar("y_" + "_".join(map(str, B)))
    for k in range(1, d + 1):
        pass  # singletons fixed =1, empty=0

    big_As = []
    for k in range(d + 1, n + 1):
        for A in combinations(range(n), k):
            big_As.append(A)

    fvars = {}
    for A in big_As:
        fa = model.NewIntVar(0, 1, "f_" + "_".join(map(str, A[:3])) + ("_" + str(len(A)) if len(A) > 3 else ""))
        fvars[A] = fa
        const = 0
        terms = []
        for k in range(2, d + 1):
            for B in combinations(A, k):
                K = comb(len(A) - k - 1, d - k) * (1 if (d - k) % 2 == 0 else -1)
                if K:
                    terms.append(K * y[B])
        # singletons contribution
        K1 = comb(len(A) - 2, d - 1) * (1 if (d - 1) % 2 == 0 else -1)
        const = len(A) * K1
        model.Add(sum(terms) + const == fa)
    return model, y, fvars, big_As


def solve(n, d, time_limit=300.0, workers=8, hint=None, log=False):
    model, y, fvars, big_As = build(n, d)
    if hint:
        for B, v in hint.items():
            if B in y:
                model.AddHint(y[B], v)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    solver.parameters.log_search_progress = log
    st = solver.Solve(model)
    name = solver.StatusName(st)
    if name in ("OPTIMAL", "FEASIBLE"):
        sol = {B: solver.Value(v) for B, v in y.items()}
        return name, sol, solver
    return name, None, solver


def verify(n, d, sol):
    """Fully independent verification, recomputing from the raw y data only.

    1. rebuild f on the whole cube from the Mobius transform of the y-values,
    2. check f is Boolean, f(0)=0, f({i})=1 for all i,
    3. check degree <= d by direct Mobius inversion of f (independent route),
    4. check the FrontierMath inequality  n > d^(log6/log3)  exactly (with Fractions-free
       high precision) and n >= 14 etc. as required.
    """
    import math
    from fractions import Fraction as Fr
    yfull = {}
    for B, v in sol.items():
        yfull[frozenset(B)] = v
    # rebuild on all subsets: use Mobius/up-set sums. Define c(B) = mobius(y)(B) for |B|<=d,
    # and y(A) = sum_{B subset A} c(B) for all A.
    c = {}
    for k in range(0, d + 1):
        for B in combinations(range(n), k):
            Bs = frozenset(B)
            if k == 0:
                c[Bs] = 0
            elif k == 1:
                c[Bs] = 1
            else:
                val = 0
                for j in range(0, k + 1):
                    for C in combinations(B, j):
                        Cs = frozenset(C)
                        base = 0 if j == 0 else (1 if j == 1 else yfull.get(Cs))
                        if base is None:
                            return False, "missing y for %s" % (Cs,)
                        val += (-1) ** (k - j) * base
                c[Bs] = val
    f = {}
    for k in range(0, n + 1):
        for A in combinations(range(n), k):
            As = frozenset(A)
            v = sum(c[Bs] for j in range(0, min(k, d) + 1) for Bs in (frozenset(C) for C in combinations(A, j)))
            f[As] = v
            if v not in (0, 1):
                return False, "f not Boolean at %s: %d" % (A, v)
    if f[frozenset()] != 0:
        return False, "f(0) != 0"
    for i in range(n):
        if f[frozenset([i])] != 1:
            return False, "f(e_%d) != 1" % i
    # independent degree check: Mobius inversion of the *rebuilt* f
    for m in range(d + 1, n + 1):
        for A in combinations(range(n), m):
            s = 0
            for j in range(0, m + 1):
                for B in combinations(A, j):
                    s += (-1) ** (m - j) * f[frozenset(B)]
            if s != 0:
                return False, "degree > d: c[%s] = %d" % (A, s)
    # FrontierMath quality gate
    from exponent_math import pow_gt
    if not pow_gt(n, d):
        return True, "OK (but n <= d^c)"
    return True, None


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
    tag = sys.argv[4] if len(sys.argv) > 4 else ""
    wk = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    name, sol, solver = solve(n, d, time_limit=tl, workers=wk)
    print("n=%d d=%d status=%s wall=%.1f" % (n, d, name, solver.WallTime()))
    if sol:
        ok, err = verify(n, d, sol)
        if ok and err is not None:
            print("independent verification:", err)
        elif ok:
            print("independent verification: OK")
        else:
            print("independent verification:", err)
        with open("experiments/ksol_%d_%d%s.json" % (n, d, tag), "w") as fh:
            json.dump({",".join(map(str, sorted(B))): v for B, v in sol.items()}, fh)
        print("saved experiments/ksol_%d_%d%s.json" % (n, d, tag))
