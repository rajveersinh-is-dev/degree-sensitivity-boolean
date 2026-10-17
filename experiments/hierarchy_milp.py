"""Integrality-augmented hierarchy: all level counts are integers.

All variables:
  b_k  := #{B in F, |B|=k}, k=2..n            (integers)
  c^{(i)}_k := #{B in F, i not in B, |B|=k}   (integers), i=1..n, k=2..n-1
with the master/level-1 systems in density form parametrised by free real
x = (w_2..w_d) and x^{(i)} = (w^{(i)}_2..w^{(i)}_d), respectively:

  b_k  = C(n,k)  * w_k(x)                (affine in x; w_k for k>d determined)
  c^{(i)}_k = C(n-1,k) * w^{(i)}_k(x^{(i)})

Constraints: 0 <= b_k <= C(n,k); 0 <= c^{(i)}_k <= C(n-1,k);
  sum_i c^{(i)}_k = (n-k) b_k            (k=2..n-1)
  0 <= b_k - c^{(i)}_k <= C(n-1,k-1)     (each i,k)

Integrality of all b_k, c^{(i)}_k is REQUIRED for any genuine example.
This MILP is a pure relaxation (all constraints are necessary), so
INFEASIBLE => no degree-<=d Boolean function has all n directions sensitive.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import comb
import sys

from scipy.optimize import milp, LinearConstraint, Bounds
import numpy as np



def build_system(N, d) -> tuple:
    """Create system.
    
    Args:
        N:
        d:
    
    Returns:
        tuple: Result of type tuple
    
    """
    nv = d - 1
    aff = {}
    for m in range(2, d + 1):
        c = [Fr(0)] * nv
        c[m - 2] = Fr(1)
        aff[m] = (Fr(0), tuple(c))
    for m in range(d + 1, N + 1):
        const = Fr(m)
        coef = [Fr(0)] * nv
        for k in range(2, m):
            ck = Fr(comb(m, k) * (1 if k % 2 == 0 else -1))
            cc, vv = aff[k]
            const -= ck * cc
            for i in range(nv):
                coef[i] -= ck * vv[i]
        sm = 1 if m % 2 == 0 else -1
        aff[m] = (sm * const, tuple(sm * z for z in coef))
    return aff, nv


def milp_feasible(n, d, master_only=False, time_limit=None) -> tuple:
    """Return (status, detail). status in {0 feasible, 2 infeasible, other}."""
    nv = d - 1
    N = n - 1
    aff_m, _ = build_system(n, d)
    aff_i, _ = build_system(N, d)

    # variable layout: x (nv), then b: indices k=2..n (n-1 vars),
    # then c: for each i (n) and k=2..N (N-1 vars)
    nb = n - 1
    nc = 0 if master_only else n * (N - 1)
    nx = nv if not master_only else nv  # x always needed (both systems share free x? no!)
    # NOTE: master x and each c-system x^{(i)} are separate free vectors.
    if master_only:
        nx_tot = nv
        layout = dict(x0=0, b=nv, c=None)
    else:
        nx_tot = nv * (1 + n)
        layout = dict(x0=0, x_i=lambda i: nv * (1 + i), b=nv * (1 + n), c=nv * (1 + n) + nb)
    ntot = (nx_tot if not master_only else nv) + nb + nc

    A_eq = []
    b_eq = []
    lb = np.full(ntot, -np.inf)
    ub = np.full(ntot, np.inf)

    # b_k = C(n,k) w_k(x0)
    for k in range(2, n + 1):
        const, coef = aff_m[k]
        row = np.zeros(ntot)
        for j, z in enumerate(coef):
            row[layout["x0"] + j] = -float(z) * comb(n, k)
        row[layout["b"] + (k - 2)] = 1.0
        A_eq.append(row)
        b_eq.append(float(const) * comb(n, k))
        lb[layout["b"] + (k - 2)] = 0.0
        ub[layout["b"] + (k - 2)] = comb(n, k)

    if not master_only:
        for i in range(n):
            xo = layout["x_i"](i)
            for k in range(2, N + 1):
                const, coef = aff_i[k]
                row = np.zeros(ntot)
                for j, z in enumerate(coef):
                    row[xo + j] = -float(z) * comb(N, k)
                row[layout["c"] + i * (N - 1) + (k - 2)] = 1.0
                A_eq.append(row)
                b_eq.append(float(const) * comb(N, k))
                idx = layout["c"] + i * (N - 1) + (k - 2)
                lb[idx] = 0.0
                ub[idx] = comb(N, k)

    # ---- pure t-space constraints ----
    A_rows = []
    b_rows = []

    def bcol(k):
        """Bcol.
        
        Args:
            k:
        
        Returns:
            The computed result
        
        """
        return layout["b"] + (k - 2)

    if not master_only:
        def ccol(i, k):
            """Ccol.
            
            Args:
                i:
                k:
            
            Returns:
                The computed result
            
            """
            return layout["c"] + i * (N - 1) + (k - 2)

        for k in range(2, N + 1):
            # sum_i c^{(i)}_k = (n-k) b_k
            row = np.zeros(ntot)
            for i in range(n):
                row[ccol(i, k)] += 1.0
            row[bcol(k)] -= (n - k)
            A_rows.append(row)
            b_rows.append(0.0)
            A_rows.append(-row)
            b_rows.append(0.0)
            # 0 <= b_k - c^{(i)}_k <= C(N,k-1)
            for i in range(n):
                row2 = np.zeros(ntot)
                row2[bcol(k)] += 1.0
                row2[ccol(i, k)] -= 1.0
                A_rows.append(row2)
                b_rows.append(comb(N, k - 1))
                A_rows.append(-row2)
                b_rows.append(0.0)

    constraints = [LinearConstraint(np.array(A_eq), np.array(b_eq), np.array(b_eq))]
    if A_rows:
        constraints.append(LinearConstraint(np.array(A_rows), -np.inf, np.array(b_rows)))

    integrality = np.zeros(ntot)
    integrality[layout["b"]:layout["b"] + nb] = 1
    if not master_only:
        integrality[layout["c"]:] = 1

    options = {}
    if time_limit:
        options["time_limit"] = time_limit
    res = milp(c=np.zeros(ntot), constraints=constraints, bounds=Bounds(lb, ub),
               integrality=integrality, options=options)
    return res.status, res


if __name__ == "__main__":
    cases = [(9, 4), (10, 4), (11, 4), (14, 5), (15, 5), (16, 5)]
    print("master_only / full hierarchy integrality MILP")
    for (n, d) in cases:
        for mo in (False,):
            st, res = milp_feasible(n, d, master_only=mo, time_limit=120)
            print(f"n={n} d={d} master_only={mo}: status={st}")
            if st == 0:
                nv = d - 1
                # recover b values
                nb = n - 1
                off_b = nv * (1 + n) if not mo else nv
                bs = res.x[off_b:off_b + nb]
                print("   b_k =", ", ".join(f"{comb(n,k)}*w_{k}" for k in range(2, n + 1)))
                print("   approx b_k =", ", ".join(f"{v:.4f}" for v in bs))
        sys.stdout.flush()
