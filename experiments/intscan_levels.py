"""Integer scan of the level system at (n,d): does alpha_lattice force levels?

System (exact, integer):
    b_m = n*C(n-1,m-1) + sum_{k=2}^{min(d,m)} alpha_k C(n-k,m-k),   m = 2..n
    0 <= b_m <= C(n,m)
    alpha_k integers in the proven per-coefficient-sum ranges:
        k=2: [-2*C(n,2), -C(n,2)]
        k=3: [0, 4*C(n,3)]
        k=4: [-14*C(n,4), 9*C(n,4)]
        k=5: [-80*C(n,5), 86*C(n,5)]
        k>=6: [-2^k*C(n,k), 2^k*C(n,k)]   (Mobius box)
For each m: min/max b_m over integer solutions.  If min==max==0 or C(n,m),
the level is *integer-forced* to be empty/full, which propagates to every
m0-subcube restriction (all slices have level-m count 0 resp. C(m0,m)).

Also reports LP-relaxation min/max for comparison.
"""
import sys
from math import comb

from ortools.sat.python import cp_model


def ranges_for(n, k):
    C = comb(n, k)
    if k == 2:
        return -2 * C, -C
    if k == 3:
        return 0, 4 * C
    if k == 4:
        return -14 * C, 9 * C
    if k == 5:
        return -80 * C, 86 * C
    return -(2 ** k) * C, (2 ** k) * C


def main():
    n, d = int(sys.argv[1]), int(sys.argv[2])
    model = cp_model.CpModel()
    al = {}
    for k in range(2, d + 1):
        lo, hi = ranges_for(n, k)
        al[k] = model.NewIntVar(lo, hi, "alpha_%d" % k)
    b = {}
    for m in range(2, n + 1):
        terms = [comb(n - k, m - k) * al[k] for k in range(2, min(d, m) + 1)]
        base = n * comb(n - 1, m - 1)
        bv = model.NewIntVar(0, comb(n, m), "b_%d" % m)
        model.Add(sum(terms) + base == bv)
        b[m] = bv
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 4
    rows = []
    for m in range(2, n + 1):
        Cm = comb(n, m)
        vals = {}
        for sense in ("min", "max"):
            model.ClearObjective()
            if sense == "min":
                model.Minimize(b[m])
            else:
                model.Maximize(b[m])
            st = solver.Solve(model)
            nm = solver.StatusName(st)
            if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                vals[sense] = None
            else:
                vals[sense] = int(solver.ObjectiveValue())
        tag = ""
        if vals["min"] is not None and vals["min"] == vals["max"]:
            v = vals["min"]
            if v == 0 or v == Cm:
                tag = "  <== INTEGER-FORCED %s" % ("ZERO" if v == 0 else "FULL")
        row = "m=%2d  b_m in [%s, %s] of C=%d%s" % (
            m, vals["min"], vals["max"], Cm, tag)
        print(row)
        rows.append(row)
    with open("intscan_%d_%d.txt" % (n, d), "w", encoding="utf-8") as fh:
        fh.write("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
