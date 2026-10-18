"""
Degree vs point-sensitivity feasibility study.

We look for a multilinear polynomial P with integer coefficients on n variables such that
  * P is Boolean on the whole cube {0,1}^n, i.e. P(S) in {0,1} for every S subset [n],
  * P(empty) = 0,
  * P({i}) = 1 for all i in [n]  (so the sensitivity at 0 is exactly n),
  * deg(P) <= d.

Because P(S) = sum_{A subset S} c_A is linear in the unknown coefficients c_A, feasibility is
an integer linear program.  We solve it with OR-Tools CP-SAT.

Reporting: for each (n, d) whether a solution exists (and one solution if so).
"""
from __future__ import annotations

from itertools import combinations

from ortools.sat.python import cp_model


BOUND = 10**6


def solve(n, d, time_limit=60.0, extra_sensitive=None) -> tuple:
    """Return (status_name, coeffs) where coeffs is a dict keyed by frozenset.

    extra_sensitive: if not None, the set of indices forced sensitive; default = all of [n].
    """
    if extra_sensitive is None:
        sensitive = set(range(n))
    else:
        sensitive = set(extra_sensitive)

    model = cp_model.CpModel()
    coeffs = {}
    for k in range(1, d + 1):
        for A in combinations(range(n), k):
            key = frozenset(A)
            if k == 1:
                val = 1 if A[0] in sensitive else 0
                coeffs[key] = val
            else:
                v = model.NewIntVar(-BOUND, BOUND, f"c_{'_'.join(map(str, A))}")
                coeffs[key] = v

    # constraints: for every subset S of [n], 0 <= P(S) <= 1
    for size in range(0, n + 1):
        for S in combinations(range(n), size):
            if size == 0:
                continue  # P(empty)=0 automatically since we omitted the constant term
            terms = []
            for k in range(1, min(d, size) + 1):
                for A in combinations(S, k):
                    c = coeffs[frozenset(A)]
                    if isinstance(c, int):
                        if c:
                            terms.append(c)
                    else:
                        terms.append(c)
            expr = sum(terms) if terms else 0
            if isinstance(expr, int):
                if expr not in (0, 1):
                    return ("INTEGER_INFEASIBLE", None)
            else:
                model.Add(expr >= 0)
                model.Add(expr <= 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = 8
    status = solver.Solve(model)
    name = solver.StatusName(status)
    if name in ("OPTIMAL", "FEASIBLE"):
        sol = {A: (c if isinstance(c, int) else solver.Value(c)) for A, c in coeffs.items()}
        return ("SAT", sol)
    return (name, None)


def check(n, d, sol) -> tuple:
    """Verify the returned solution independently."""
    for size in range(0, n + 1):
        for S in combinations(range(n), size):
            val = 0
            for A, c in sol.items():
                if A <= frozenset(S):
                    val += c
            if val not in (0, 1):
                return False, (S, val)
    return True, None


if __name__ == "__main__":
    print("n d status  (all singletons sensitive)")
    for d in range(2, 7):
        for n in range(d, 3 * d + 3):
            status, sol = solve(n, d, time_limit=20.0)
            ok = ""
            if status == "SAT":
                good, bad = check(n, d, sol)
                ok = "verified" if good else f"VERIFY-FAIL {bad}"
            print(f"{n} {d} {status} {ok}", flush=True)
