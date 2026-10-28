"""Scan level-LP feasibility over ALL sizes n' in [d, n_min(d)] for each d.

If the level system is infeasible at some (n', d) with d <= n' <= n_min(d),
then no fully sensitive degree-d function exists on n_min(d) variables either
(restriction to n' coords containing a degree-d monomial), and hence the whole
degree-d case of the FrontierMath problem is dead.

Prints a matrix; for infeasible entries tries an exact Farkas certificate.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amanotype as AT
from exponent_math import nmin


def main():
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    ranges = "ranges" in sys.argv
    lines = []
    for d in range(2, dmax + 1):
        nm, _ = nmin(d)
        bad = []
        row = []
        for n in range(d, nm + 1):
            rows = AT.base_rows(n, d, ranges)
            feas, pt = AT.linprog_feasible(rows, d - 1)
            if feas is False:
                bad.append(n)
            row.append("F" if feas else (". " if feas is None else "I"))
        # compress
        msg = "d=%2d n_min=%3d  [n=d..n_min] %s   INFEASIBLE at: %s" % (
            d, nm, "".join(row), bad)
        print(msg)
        lines.append(msg)
        for n in bad:
            rows = AT.base_rows(n, d, ranges)
            y, err = AT.farkas_numeric_support(rows, d - 1)
            if y is None:
                m2 = "      (n=%d) no exact cert: %s" % (n, err)
            else:
                val = sum(rows[i][1] * y[i] for i in range(len(y)))
                nz = sum(1 for z in y if z != 0)
                m2 = "      (n=%d) EXACT FARKAS CERT: %d nonzeros, b.y=%s" % (
                    n, nz, val)
            print(m2)
            lines.append(m2)
    with open("scan2_out.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
