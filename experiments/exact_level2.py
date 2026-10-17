"""Exact feasibility of LP(n,d) (level-density necessary conditions).

Strategy: build the exact rational system, solve numerically with HiGHS, then
(a) if feasible: rationalize an interior point (max-slack LP) and verify it
    exactly with Fractions -> certificate of FEASIBILITY;
(b) if infeasible: search for an exact Farkas certificate y >= 0 with
    A^T y = 0 and b.y < 0, solving a small exact linear system on a chosen
    support -> certificate of INFEASIBILITY (complete, since exact).

Output to UTF-8 file (Windows console is cp1252).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import comb
import itertools

from scipy.optimize import linprog
import numpy as np


OUT = []


def pr(*a):
    """Pr.
    
    """
    OUT.append(" ".join(str(x) for x in a))


def build_system(n, d) -> tuple:
    """Create system.
    
    Args:
        n:
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
    for m in range(d + 1, n + 1):
        const = Fr(m)
        coef = [Fr(0)] * nv
        for k in range(2, m):
            ck = Fr(comb(m, k) * (1 if k % 2 == 0 else -1))  # (-1)^k C(m,k)
            cc, vv = aff[k]
            const -= ck * cc
            for i in range(nv):
                coef[i] -= ck * vv[i]
        sign_m = 1 if m % 2 == 0 else -1
        aff[m] = (sign_m * const, tuple(sign_m * z for z in coef))
    rows = []
    for m in range(2, n + 1):
        c, v = aff[m]
        # w_m >= 0  <=>  -(c + v.x) <= 0  <=>  (-v).x <= c
        rows.append((tuple(-z for z in v), c))
        # w_m <= 1  <=>  v.x <= 1 - c
        rows.append((tuple(v), Fr(1) - c))
    return rows, aff


def check_rows(rows, x) -> bool:
    """Check whether rows.
    
    Args:
        rows:
        x:
    
    Returns:
        bool: Result of type bool
    
    """
    for (c, b) in rows:
        s = sum(cc * xx for cc, xx in zip(c, x))
        if s > b:
            return False
    return True


def rat(v, max_den=10 ** 8) -> tuple:
    """Rat.
    
    Args:
        v:
        max_den:
    
    Returns:
        tuple: Result of type tuple
    
    """
    return tuple(Fr(float(x)).limit_denominator(max_den) for x in v)


def exact_solve(mat, rhs) -> tuple:
    """Solve square rational system A x = b by Gaussian elimination.
    Returns solution tuple or None if singular."""
    m = len(mat)
    A = [[Fr(x) for x in row] + [Fr(rhs[i])] for i, row in enumerate(mat)]
    for col in range(m):
        piv = None
        for r in range(col, m):
            if A[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        for r in range(m):
            if r != col and A[r][col] != 0:
                f = A[r][col] / pv
                for k in range(col, m + 1):
                    A[r][k] -= f * A[col][k]
    return tuple(A[i][m] / A[i][i] for i in range(m))


def feas_point(rows, nv) -> tuple:
    """Feas point.
    
    Args:
        rows:
        nv:
    
    Returns:
        tuple: Result of type tuple
    
    """
    A = np.array([[float(z) for z in c] for (c, b) in rows])
    b = np.array([float(x) for (c, x) in rows])
    # maximin slack:  A x + t <= b,  maximize t
    A2 = np.hstack([A, np.ones((A.shape[0], 1))])
    res = linprog(c=np.array([0.0] * nv + [-1.0]), A_ub=A2, b_ub=b,
                  bounds=[(None, None)] * nv + [(None, None)], method="highs")
    if res.status == 0 and res.x[nv] > 0:
        pt = rat(res.x[:nv])
        if check_rows(rows, pt):
            return pt, res.x[nv]
        # fallback: round more aggressively a few times
        for md in (10 ** 6, 10 ** 4, 10 ** 2):
            pt = rat(res.x[:nv], md)
            if check_rows(rows, pt):
                return pt, res.x[nv]
    # plain feasibility
    res = linprog(c=np.zeros(nv), A_ub=A, b_ub=b, bounds=[(None, None)] * nv,
                  method="highs")
    if res.status == 0:
        pt = rat(res.x[:nv])
        if check_rows(rows, pt):
            return pt, 0.0
    return None, None


def farkas(rows, nv) -> tuple:
    """Find exact y >= 0 with A^T y = 0, b.y = -1 (Farkas certificate of
    infeasibility of {Ax <= b})."""
    m = len(rows)
    A = np.array([[float(z) for z in c] for (c, b) in rows])
    b = np.array([float(x) for (c, x) in rows])
    # numeric: minimize b.y  s.t. A^T y = 0, y >= 0, sum y = 1
    Aeq = np.vstack([A.T, np.ones((1, m))])
    beq = np.zeros(nv + 1)
    beq[-1] = 1.0
    res = linprog(c=b, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * m, method="highs")
    if res.status != 0 or b @ res.x > -1e-9:
        return None, "no numeric ray (status %s, val %.3g)" % (res.status, b @ res.x)
    dbl = np.where(np.abs(res.x) > 1e-7)[0]
    # exact: choose nv+1 indices from support, solve A^T y = 0 (nv rows), b.y = -1
    best = None
    sup = list(dbl)
    for sub in itertools.combinations(sup, nv + 1):
        mat = [[Fr(A[j, i]) for j in sub] for i in range(nv)]
        mat.append([Fr(b[j]) for j in sub])
        rhs = [Fr(0)] * nv + [Fr(-1)]
        sol = exact_solve(mat, rhs)
        if sol is None:
            continue
        if all(s >= 0 for s in sol):
            best = {j: s for j, s in zip(sub, sol)}
            break
    if best is None:
        return None, f"no exact rational ray on support {list(dbl)}"
    # verify exactly: A^T y = 0, y >= 0, b.y < 0
    y = [best.get(j, Fr(0)) for j in range(m)]
    for i in range(nv):
        s = sum(Fr(A[j, i]) * y[j] for j in range(m))
        if s != 0:
            return None, f"certificate failed row {i}"
    s = sum(Fr(b[j]) * y[j] for j in range(m))
    if s >= 0 or any(t < 0 for t in y):
        return None, "certificate failed sign/sum"
    return (y, s), None


def run_case(n, d) -> bool:
    """Worker function for case.
    
    Args:
        n:
        d:
    
    Returns:
        bool: Result of type bool
    
    """
    pr("=" * 72)
    pr(f"case (n={n}, d={d})")
    rows, aff = build_system(n, d)
    nv = d - 1
    pt, slack = feas_point(rows, nv)
    if pt is not None:
        pr("  FEASIBLE (exact rational point verified), max-slack ~ %.3g" % (slack if slack else 0))
        pr("  x = (w_2..w_d) =", ", ".join(str(z) for z in pt))
        w = {}
        for m in range(2, n + 1):
            c0, v = aff[m]
            w[m] = c0 + sum(a * b for a, b in zip(v, pt))
        pr("  w_k:", ", ".join(f"w_{k}={...}" for k in range(2, n + 1)))
        pr("  b_k:", ", ".join(f"b_{k}={... * comb(n, k)}" for k in range(2, n + 1)))
        return True
    cert, err = farkas(rows, nv)
    if cert is not None:
        y, s = cert
        pr(f"  INFEASIBLE (exact Farkas certificate, b.y = {s})")
        nonzero = [(i, y[i]) for i in range(len(y)) if y[i] != 0]
        pr(f"  certificate support: {len(nonzero)} rows of {len(rows)}")
        return False
    pr(f"  UNKNOWN: numerical infeasibility but no exact certificate found ({err})")
    return None


if __name__ == "__main__":
    cases = [(6, 3), (7, 3), (8, 3), (9, 4), (10, 4), (11, 4),
             (13, 5), (14, 5), (15, 5), (16, 5), (18, 6), (19, 6)]
    for (n, d) in cases:
        try:
            run_case(n, d)
        except Exception as e:
            pr("  ERROR:", repr(e))
    with open("exact_level2_out.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("\n".join(OUT))
