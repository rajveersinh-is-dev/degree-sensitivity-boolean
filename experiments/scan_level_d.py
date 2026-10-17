"""Scan the level-LP feasibility at the minimal candidate (n_min(d), d).

n_min(d) = floor(d^c) + 1, c = log_6(log 3)=ln6/ln3  (the least n with n > d^c).

If the level system is infeasible at (n_min(d), d), then by the restriction
argument no fully sensitive degree-<=d function exists on ANY n >= n_min(d),
so that degree d is dead for the FrontierMath problem.

Rows are built WITHOUT any coefficient-range extras (pure count constraints
0 <= w_m <= 1 plus the (star) relations), so infeasibility certificates only
use the most elementary necessary conditions.
"""
import sys
sys.path.insert(0, r"C:\Users\davea\Downloads\newp\experiments")
import amanotype as AT
from decimal import Decimal, getcontext

getcontext().prec = 80
C = Decimal(6).ln() / Decimal(3).ln()


def nmin(d):
    x = Decimal(d) ** C
    n = int(x) + 1
    while Decimal(n) <= x:
        n += 1
    while Decimal(n - 1) > x:
        n -= 1
    return n, x


def main():
    out = []
    lo, hi = 2, int(sys.argv[1]) if len(sys.argv) > 1 else 14
    ranges = "ranges" in sys.argv
    for d in range(lo, hi + 1):
        n, x = nmin(d)
        rows = AT.base_rows(n, d, ranges)
        feas, pt = AT.linprog_feasible(rows, d - 1)
        msg = "d=%2d n_min=%3d  d^c=%s  feasible=%s" % (
            d, n, str(x)[:22], feas)
        print(msg)
        out.append(msg)
        if feas is False:
            y, err = AT.farkas_numeric_support(rows, d - 1)
            if y is None:
                m2 = "     numeric infeasible, no exact cert: %s" % err
            else:
                nz = [(i, y[i]) for i in range(len(y)) if y[i] != 0]
                val = sum(rows[i][1] * y[i] for i in range(len(y)))
                m2 = ("     EXACT FARKAS CERT ok: %d nonzeros, b.y=%s"
                      % (len(nz), val))
            print(m2)
            out.append(m2)
    with open("scan_level_d_out.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
