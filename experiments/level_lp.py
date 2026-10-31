"""Level-density LP relaxation for the sensitivity-vs-degree problem.

Derivation (see research notes):
Let f: {0,1}^n -> {0,1}, f(0)=0, f(e_i)=1 for all i (all coords sensitive at 0),
deg f <= d.  Let F = f^{-1}(1), b_k = #(F intersect level k), w_k = b_k / C(n,k).
The Mobius coefficients c_A = sum_{B subset A} (-1)^{|A|-|B|} f(B) vanish for |A| >= d+1.
Summing that condition over all m-subsets A (m = |A|) gives the necessary system

  (*)  sum_{k=2}^{m} (-1)^k w_k C(m,k) = m      for m = d+1,...,n.

Let T(m) := sum_{k=0}^m (-1)^k C(m,k) w_k (with w_0 = 0, w_1 = 1); (*) says T(m)=0
for m in [d+1,n].  Since the binomial transform is an involution, this is equivalent to
the existence of a real polynomial

  H(x) = C(x,1) + sum_{m=2}^{d} b_m C(x,m)     (b_m real, free)

satisfying the LINEAR INEQUALITIES

  0 <= H(k) <= 1   for all k = 2,...,n       (indeed H(k) = w_k).

So LP(n,d) asks: does there exist a degree-<=d polynomial with binomial-basis
coefficient of C(x,1) equal to 1, taking values in [0,1] at the integers 2..n?

Necessity: every genuine (n,d) example implies LP(n,d) feasible.  Hence the true
max sensitivity N(d) satisfies N(d) <= n_LP(d) := max{n : LP(n,d) feasible},
and LP feasibility is monotone in n (a solution for n' restricts to one for n<=n').
We probe n_LP(d) for small d by LP.
"""
from math import comb, log
import numpy as np
from scipy.optimize import linprog
from exponent_math import C_CRIT, Decimal


def lp_feasible(n: int, d: int):
    """LP(n,d) feasibility. Returns (status, x) where status 0 = optimal (feasible)."""
    assert d >= 2 and n >= 2
    if n <= d:
        # No equalities beyond degree; the all-ones level structure w_k=1 is fine?
        # Still must satisfy inequalities H(k) in [0,1]; try the explicit interpolant
        # H(x) = x (b_m = 0) which gives H(k) = k > 1 for k >= 2.  So still need LP.
        pass
    num = d - 1  # variables b_2..b_d
    A = np.array([[comb(k, m) for m in range(2, d + 1)] for k in range(2, n + 1)], dtype=float)
    lb = np.array([-float(k) for k in range(2, n + 1)])
    ub = np.array([1.0 - float(k) for k in range(2, n + 1)])
    A_ub = np.vstack([A, -A])
    b_ub = np.concatenate([ub, -lb])
    row_norms = np.abs(A_ub).max(axis=1)
    row_norms[row_norms == 0] = 1.0
    res = linprog(c=np.zeros(num), A_ub=A_ub / row_norms[:, np.newaxis], b_ub=b_ub / row_norms,
                  bounds=[(None, None)] * num, method="highs")
    return res.status, (res.x if res.status == 0 else None)


def n_lp(d: int, nmax: int = 300, tol: float = 1e-7):
    """max n with LP(n,d) feasible, scanning up from 2 (monotone in n)."""
    n = 2
    last_ok = None
    while n <= nmax:
        st, _ = lp_feasible(n, d)
        if st == 0:
            last_ok = n
        else:
            return last_ok, n   # last feasible, first infeasible
        n += 1
    return last_ok, None


if __name__ == "__main__":
    out = []
    out.append(f"c = log6/log3 = {C_CRIT!r}")
    out.append(f"{'d':>3} {'n_LP(d)':>8} {'d^c':>9} {'need n>=':>9} {'verdict':>8}")
    for d in range(2, 33):
        mx, firstbad = n_lp(d, nmax=400)
        dc = Decimal(d) ** C_CRIT
        need = int(dc) + 1
        verdict = "KILLED" if (mx is not None and mx < need) else "OPEN"
        out.append(f"{d:>3} {str(mx):>8} {float(dc):>9.3f} {need:>9} {verdict:>8}")
    text = "\n".join(out)
    print(text)
    with open("level_lp_out.txt", "w", encoding="utf-8") as fh:
        fh.write(text + "\n")
