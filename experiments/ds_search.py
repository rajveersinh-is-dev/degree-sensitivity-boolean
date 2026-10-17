"""Search for Boolean multilinear polynomials of degree d with many sensitive coordinates at 0.

For a given (n, d) we ask: does there exist an integer multilinear polynomial P on n variables
with P Boolean on {0,1}^n, P(0)=0, P(e_i)=1 for all i, and deg(P) <= d?

Equivalently: does there exist a family F subset 2^[n], with empty set not in F, all singletons
in F, whose Moebius transform vanishes on all sets of size > d?

We print solutions (as coefficient dictionaries) to inspect their structure.
"""
from __future__ import annotations

from itertools import combinations
import json
import sys

from ortools.sat.python import cp_model



BOUND = 10 ** 6


def build(n, d) -> tuple:
    """Create the object.
    
    Args:
        n:
        d:
    
    Returns:
        tuple: Result of type tuple
    
    """
    model = cp_model.CpModel()
    coeffs = {}
    for k in range(2, d + 1):
        for A in combinations(range(n), k):
            coeffs[frozenset(A)] = model.NewIntVar(-BOUND, BOUND, "c" + "_".join(map(str, A)))
    for size in range(2, n + 1):
        for S in combinations(range(n), size):
            expr = size + sum(coeffs[frozenset(A)] for k in range(2, min(d, size) + 1)
                              for A in combinations(S, k))
            model.Add(expr >= 0)
            model.Add(expr <= 1)
    return model, coeffs


def solve(n, d, time_limit, hint=None, workers=8) -> tuple:
    """Solve.
    
    Args:
        n:
        d:
        time_limit:
        hint:
        workers (int):
    
    Returns:
        tuple: Result of type tuple
    
    """
    model, coeffs = build(n, d)
    if hint:
        for A, v in hint.items():
            if A in coeffs:
                model.AddHint(coeffs[A], v)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    st = solver.Solve(model)
    name = solver.StatusName(st)
    if name in ("OPTIMAL", "FEASIBLE"):
        return name, {A: solver.Value(c) for A, c in coeffs.items()}, solver
    return name, None, solver


def verify(n, d, coeffs):
    """Independently verify Booleanness and degree."""
    bad = []
    for size in range(2, n + 1):
        for S in combinations(range(n), size):
            val = size + sum(v for A, v in coeffs.items() if A <= frozenset(S))
            if val not in (0, 1):
                bad.append((S, val))
    return bad


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    tl = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
    name, sol, solver = solve(n, d, tl)
    print("n=%d d=%d status=%s wall=%.1f" % (n, d, name, solver.WallTime()))
    if sol:
        bad = verify(n, d, sol)
        print("independent verification:", "OK" if not bad else bad[:5])
        with open(f"experiments/sol_{n}_{d}.json", "w") as f:
            json.dump({",".join(map(str, sorted(A))): v for A, v in sol.items()}, f, indent=0)
        print(f"saved experiments/sol_{n}_{d}.json")
