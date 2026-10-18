"""Exact (rational) analysis of the level-density necessary conditions.

Derivation recap (all PROVED, see research notes):
Let f: {0,1}^n -> {0,1} be Boolean with f(0)=0 and f(e_i)=1 for all i (all n
coordinate directions sensitive at the origin), and deg f <= d.  Write
F = f^{-1}(1) subset 2^[n], b_k = #{B in F : |B| = k}, w_k = b_k / C(n,k).

Mobius coefficients: c_A = sum_{B subset A} (-1)^{|A|-|B|} f(B) = 0 whenever
|A| >= d+1.   For |A| = m >= d+1 this reads
      sum_{k=2..m} (-1)^{m-k} N_k(A) = (-1)^m m,      N_k(A)=#(F cap C(A,k)).
Summing over all m-subsets A gives the *aggregate* system
  (*)  sum_{k=2..m} (-1)^k w_k C(m,k) = m          for all m = d+1,...,n.

For fixed n, the map (w_2,...,w_d) -> (w_{d+1},...,w_n) is affine and uniquely
determined by (*), so (*) plus the constraints 0 <= w_k <= 1 for all k is a
linear feasibility problem LP(n,d) in d-1 variables.  It is a *necessary*
condition (a relaxation), and feasible sets are nested in n (monotone).

Additionally b_k = w_k C(n,k) must be an *integer* for every k.

This script:
 1. builds LP(n,d) exactly with Fraction arithmetic;
 2. decides feasibility by exact Fourier-Motzkin elimination;
 3. if feasible, recovers a rational interior point / vertices;
 4. optionally tests integrality of b_k with an exact MILP-style case check
    using scipy's MILP on the same system (floating point, cross-checked).
Output is written to a UTF-8 file (Windows console is cp1252).
"""
from __future__ import annotations

from fractions import Fraction as Fr
from math import comb
import sys

OUT = []


def pr(*a):
    """Pr.
    
    """
    OUT.append(" ".join(str(x) for x in a))


# ---------- exact LP construction ----------


def build_system(n, d) -> tuple:
    """Variables x_i = w_{2+i}, i=0..d-2 (i.e. w_2..w_d).
    Rows: (coeffs tuple, rhs, sense) meaning  sum coeffs*x <= rhs  (sense '<=')
    representing 0 <= w_m <= 1 for m = 2..n, where w_m for m > d is the affine
    function of the variables given by (*).
    Returns (rows, affines) with affines[m] = (const, coeffs) for w_m.
    """
    nv = d - 1
    # affine representation of w_m: dict m -> (const Fraction, tuple of Fraction)
    aff = {}
    for m in range(2, d + 1):
        c = [Fr(0)] * nv
        c[m - 2] = Fr(1)
        aff[m] = (Fr(0), tuple(c))
    for m in range(d + 1, n + 1):
        # w_m = (-1)^m [ m - sum_{k=2}^{m-1} (-1)^k w_k C(m,k) ]
        const = Fr(m)
        coef = [Fr(0)] * nv
        for k in range(2, m):
            sgn = Fr((-1) ** k) if k % 2 == 0 else Fr(-1)
            ck = Fr(comb(m, k) * (1 if k % 2 == 0 else -1))  # (-1)^k*C(m,k)
            cc, vv = aff[k]
            const -= ck * cc
            for i in range(nv):
                coef[i] -= ck * vv[i]
        sign_m = 1 if m % 2 == 0 else -1
        const *= sign_m
        coef = [sign_m * z for z in coef]
        aff[m] = (const, tuple(coef))
    rows = []
    for m in range(2, n + 1):
        c, v = aff[m]
        # w_m >= 0  <=>  -w_m <= 0
        rows.append((tuple(-z for z in v), -c))
        # w_m <= 1
        rows.append((tuple(v), Fr(1) - c))
    return rows, aff


def fourier_motzkin(rows) -> tuple:
    """Exact FM elimination. rows: list of (coeff tuple, rhs) for coeff*x <= rhs.
    Returns (feasible: bool, new_rows).  Variables eliminated left to right.
    Returns also a witness point if feasible (back-substitution is not done here)."""
    if not rows:
        return True, rows
    nv = len(rows[0][0])
    for j in range(nv):
        pos, neg, zer = [], [], []
        for (c, b) in rows:
            a = c[j]
            if a > 0:
                pos.append((c, b))
            elif a < 0:
                neg.append((c, b))
            else:
                zer.append((c, b))
        if not pos and not neg:
            rows = zer
            continue
        new = list(zer)
        for (cp, bp) in pos:
            for (cn, bn) in neg:
                ap, an = cp[j], cn[j]
                # (p.x <= P, p_j>0), (n.x <= N, n_j<0):
                # lower <= upper  <=>  sum_i (p_j n_i - n_j p_i) x_i <= p_j N - n_j P
                c = tuple(ap * cn[i] - an * cp[i] for i in range(nv))
                b = ap * bn - an * bp
                # normalize
                g = 0
                for z in c:
                    g = gcd(g, abs(z.numerator) if z else 0)
                if g == 0:
                    if b < 0:
                        return False, None  # 0 <= b violated
                    continue
                c = tuple(z for z in c)
                new.append((c, b))
        # dedupe / drop trivial 0<=b rows with b>=0
        seen = set()
        rows2 = []
        for (c, b) in new:
            if all(z == 0 for z in c):
                if b < 0:
                    return False, None
                continue
            key = (c, b)
            if key in seen:
                continue
            seen.add(key)
            rows2.append((c, b))
        rows = rows2
        if len(rows) > 40000:
            raise RuntimeError(f"FM blowup at variable {j}: {len(rows)} rows")
    for (c, b) in rows:
        if b < 0:
            return False, None
    return True, rows


def gcd(a, b):
    """Gcd.
    
    Args:
        a:
        b:
    
    Returns:
        The computed result
    
    """
    while b:
        a, b = b, a % b
    return a


def vertices_and_point(n, d) -> tuple:
    """Numerically find vertices of LP(n,d), then verify exactly."""
    import numpy as np
    from scipy.optimize import linprog
    rows, aff = build_system(n, d)
    nv = d - 1
    A = np.array([[float(z) for z in c] for (c, b) in rows])
    b = np.array([float(x) for (_, x) in rows])
    pts = []
    # vertex enumeration: for each subset of nv constraints, solve equality system
    from itertools import combinations as C
    idx = range(len(rows))
    checked = 0
    for sub in C(idx, nv):
        M = A[list(sub)]
        try:
            det = np.linalg.det(M)
        except Exception:
            continue
        if abs(det) < 1e-9:
            continue
        sol = np.linalg.solve(M, b[list(sub)])
        checked += 1
        if np.all(A @ sol <= b + 1e-7):
            # verify exactly
            ok = exact_check(rows, sol)
            if ok:
                fr = tuple(Fr(float(x)).limit_denominator(10 ** 9) for x in sol)
                if exact_check_rows(rows, fr):
                    pts.append(fr)
    return rows, aff, pts


def exact_check(rows, sol) -> bool:
    """Exact check.
    
    Args:
        rows:
        sol:
    
    Returns:
        bool: Result of type bool
    
    """
    import numpy as np
    return bool(np.all(np.array([[float(z) for z in c] for (c, _) in rows]) @ sol
                       <= np.array([float(x) for (_, x) in rows]) + 1e-7))


def exact_check_rows(rows, fr) -> bool:
    """Exact check rows.
    
    Args:
        rows:
        fr:
    
    Returns:
        bool: Result of type bool
    
    """
    for (c, b) in rows:
        s = sum(cc * x for cc, x in zip(c, fr))
        if s > b:
            return False
    return True


def run_case(n, d, do_milp=True):
    """Worker function for case.
    
    Args:
        n:
        d:
        do_milp (bool):
    
    """
    pr("=" * 70)
    pr(f"case (n={n}, d={d})")
    rows, aff = build_system(n, d)
    feas, res = fourier_motzkin(rows)
    pr("  LP feasibility (exact FM):", "FEASIBLE" if feas else "INFEASIBLE")
    if feas:
        # find a feasible point numerically for reporting
        import numpy as np
        from scipy.optimize import linprog
        nv = d - 1
        A = np.array([[float(z) for z in c] for (c, b) in rows])
        b = np.array([float(x) for (_, x) in rows])
        r = linprog(c=np.zeros(nv), A_ub=A, b_ub=b, bounds=[(None, None)] * nv,
                    method="highs")
        pt = None
        if r.status == 0:
            pt = tuple(Fr(float(x)).limit_denominator(10 ** 7) for x in r.x)
            while not exact_check_rows(rows, pt):
                pt = tuple(x + Fr(1, 10 ** 6) for x in pt)  # nudge
        pr("  numeric LP status:", r.status)
        if pt:
            pr("  rational point w_2..w_d:", ", ".join(str(x) for x in pt))
            w = {}
            for m in range(2, n + 1):
                c0, v = aff[m]
                w[m] = c0 + sum(a * b for a, b in zip(v, pt))
            pr("  implied w_k =", ", ".join(f"w_{k}={...}" for k in range(2, n + 1)))
            pr("  implied b_k =", ", ".join(f"b_{k}={... * comb(n, k)}"
                                            for k in range(2, n + 1)))
        if do_milp:
            milp_integrality(n, d, rows, aff)


def milp_integrality(n, d, rows, aff):
    """Ask an exact question with floating MILP: does there exist a rational
    feasible point x with b_k = w_k*C(n,k) integer for EVERY k = d+1..n
    (k<=d b's are free, only need 0<=b<=C(n,k); integrality at free levels is
    automatic if we choose rational points? -> we require integrality only for
    k > d where b_k is affine in x; for k <= d b_k=C(n,k)*x_k can be continuous).
    If INFEASIBLE -> no level counts with all determined b_k integral exist."""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    nv = d - 1
    # variables: x (nv continuous) then integer t_6..t_n (n-d variables)
    nin = n - d
    m = nv + nin
    A_rows = [[float(z) for z in c] + [0.0] * nin for (c, b) in rows]
    b_rows = [float(x) for (_, x) in rows]
    # equalities: t_k - C(n,k)*w_k(x) = 0  ->  -C(n,k)*v.x + t_k = C(n,k)*c0
    A_eq = []
    b_eq = []
    for k in range(d + 1, n + 1):
        c0, v = aff[k]
        row = [-float(z) * comb(n, k) for z in v] + [0.0] * nin
        row[nv + (k - d - 1)] = 1.0
        A_eq.append(row)
        b_eq.append(float(c0) * comb(n, k))
    constraints = [LinearConstraint(np.array(A_rows), -np.inf, np.array(b_rows)),
                   LinearConstraint(np.array(A_eq), np.array(b_eq), np.array(b_eq))]
    lb = [-np.inf] * nv + [0.0] * nin
    ub = [np.inf] * nv + [float(comb(n, k)) for k in range(d + 1, n + 1)]
    integrality = [0] * nv + [1] * nin
    res = milp(c=np.zeros(m), constraints=constraints, bounds=Bounds(lb, ub),
               integrality=integrality)
    pr("  MILP (integrality of b_k for k>d):", res.status,
       "(0=feasible, 2=infeasible, other=unknown)")
    if res.status == 0:
        pr("    integer b_k:", ", ".join(
            f"b_{k}={round(...)}" for k in range(d + 1, n + 1)))


if __name__ == "__main__":
    cases = [(6, 3), (7, 3), (9, 4), (10, 4), (14, 5), (15, 5), (16, 5), (19, 6)]
    if len(sys.argv) > 2:
        cases = [(int(sys.argv[1]), int(sys.argv[2]))]
    for (n, d) in cases:
        try:
            run_case(n, d)
        except RuntimeError as e:
            pr("  FM blowup:", e)
    with open("exact_level_out.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("\n".join(OUT))
