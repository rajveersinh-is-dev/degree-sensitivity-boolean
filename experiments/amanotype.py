"""Amano-type analysis of the level-density (symmetrization) system.

Mathematical background (see research/notes; also Amano EATCS 101 (2010)).

Let f: {0,1}^n -> {0,1} be Boolean, fully sensitive at the origin:
f(0)=0, f(e_i)=1 for all i.  Write p for its multilinear polynomial, with
integer coefficients a_T.  Define

    alpha_k := sum_{|T|=k} a_T          (integer),
    b_m     := #{x : |x| = m, f(x) = 1} (integer, 0 <= b_m <= C(n,m)),
    w_m     := b_m / C(n,m).

Then the following identity holds for every m = 2..n:

    b_m = n*C(n-1,m-1) + sum_{k=2}^{min(d,m)} alpha_k * C(n-k, m-k),        (*)

equivalently  w_m = m + sum_{k=2}^{min(d,m)} a_k C(m,k),  a_k := alpha_k/C(n,k).

(Proof: sum_{|x|=m} p(x) = sum_T a_T * #{|x|=m : T subset supp x}
 = sum_k a_k C(n-k,m-k); on Boolean points p(x)=f(x) so the sum is b_m; and
 alpha_0 = 0 since p(0)=0, alpha_1 = n since all linear coefficients are 1
 by full sensitivity at the origin.)

Necessary conditions used (each proven above in comments):
  (i)   0 <= w_m <= 1 for m = 2..n.            [counts]
  (ii)  alpha_k integer.                        [Mobius]
  (iii) per-coefficient ranges:
           pair coefficients in {-2,-1}   => a_2 in [-2,-1],
           triple coefficient = f(T) + #{pairs in T with coeff -2} in {0..4}
                                          => a_3 in [0,4],
           a_4 in [-14, 9],  a_5 in [-80, 86]   (derived below),
           general Mobius box |alpha_k| <= 2^k C(n,k).
  Proofs of the degree-4/5 ranges (uses only weight<=5 Boolean-ness):
    triple: f(T) = 3 + pairsum + a_T, pairsum = -(3) + (-1)*r  (r = # -2 pairs)
            => a_T = f(T) + r in {0..4}.
    quad S: f(S) = 4 + pairsum + triplesum + a_S, pairsum = -6 - r2, r2<=6,
            triplesum = sum of the 4 triple coeffs in [0,16]
            => a_S = f(S) + 2 + r2 - triplesum in [-14, 9].
    quint Q: f(Q) = 5 + pairsum + triplesum + quadsum + a_Q,
            pairsum >= -20, triplesum in [0,40], quadsum in [-70,45]
            => a_Q in [5-20-40-45+0, 5-10+0+70+1] = [-80, 86].

Restriction fact: fixing all coordinates outside a subset S (|S|=m0) to 0
gives g fully sensitive at the origin on m0 vars with deg g <= d and
alpha'_k = sum_{T subset S, |T|=k} a_T integer.  The same system (*) holds
with n -> m0.  Levels m with b_m = C(n,m) (all-ones) or b_m = 0 (all-zero)
propagate: for m <= m0, b'_m = C(m0,m) resp. 0.

Commands:
  python amanotype.py geom  <n> <d> [ranges|noranges]
  python amanotype.py test  <n> <d> <m0> [ranges|noranges] [time_limit_s]

Output is ASCII-only and written both to stdout and out file (utf-8).
"""
from __future__ import annotations

import sys
import itertools
from fractions import Fraction as Fr
from math import comb

# ----------------------------------------------------------------------------
# exact linear algebra helpers
# ----------------------------------------------------------------------------


def exact_solve(mat, rhs) -> tuple:
    """Solve square rational system A x = b (Gauss-Jordan). None if singular."""
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


def check_rows(rows, x) -> bool:
    """Check whether rows.
    
    Args:
        rows:
        x:
    
    Returns:
        bool: Result of type bool
    
    """
    for (c, b) in rows:
        if sum(cc * xx for cc, xx in zip(c, x)) > b:
            return False
    return True


def linprog_feasible(rows, nv) -> tuple:
    """numeric feasibility check; returns True/False (None if undecided)."""
    import numpy as np
    from scipy.optimize import linprog
    A = np.array([[float(z) for z in c] for (c, b) in rows])
    b = np.array([float(x) for (c, x) in rows])
    res = linprog(c=[0.0] * nv, A_ub=A, b_ub=b,
                  bounds=[(None, None)] * nv, method="highs")
    if res.status == 0:
        return True, res.x
    if res.status == 2:
        return False, None
    return None, None


def linprog_dir(rows, nv, cvec) -> tuple:
    """numeric min of cvec.x over {rows}; returns (status, value)."""
    import numpy as np
    from scipy.optimize import linprog
    A = np.array([[float(z) for z in c] for (c, b) in rows])
    b = np.array([float(x) for (c, x) in rows])
    res = linprog(c=[float(z) for z in cvec], A_ub=A, b_ub=b,
                  bounds=[(None, None)] * nv, method="highs")
    return res.status, res.fun, res.x


def farkas_numeric_support(rows, nv) -> tuple:
    """Find y >= 0 with A^T y = 0, b.y = -1 exactly (Farkas certificate of
    infeasibility).  Returns (y, b.y) or (None, reason)."""
    import numpy as np
    from scipy.optimize import linprog
    m = len(rows)
    A = np.array([[float(z) for z in c] for (c, b) in rows])
    b = np.array([float(x) for (c, x) in rows])
    Aeq = np.vstack([A.T, np.ones((1, m))])
    beq = np.zeros(nv + 1)
    beq[-1] = 1.0
    res = linprog(c=b, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * m,
                  method="highs")
    if res.status != 0 or b @ res.x > -1e-9:
        return None, "no numeric ray (status %s val %.3g)" % (res.status,
                                                              b @ res.x)
    dbl = np.where(np.abs(res.x) > 1e-7)[0]
    sup = list(dbl)
    if len(sup) < nv + 1:
        return None, "support too small"
    for sub in itertools.combinations(sup, nv + 1):
        mat = [[Fr(A[j, i]) for j in sub] for i in range(nv)]
        mat.append([Fr(b[j]) for j in sub])
        rhs = [Fr(0)] * nv + [Fr(-1)]
        sol = exact_solve(mat, rhs)
        if sol is None:
            continue
        if all(s >= 0 for s in sol):
            y = [Fr(0)] * m
            for j, s in zip(sub, sol):
                y[j] = s
            # exact verification
            ok = True
            for i in range(nv):
                if sum(Fr(A[j, i]) * y[j] for j in range(m)) != 0:
                    ok = False
                    break
            if ok and sum(Fr(b[j]) * y[j] for j in range(m)) < 0 and \
                    all(t >= 0 for t in y):
                return y, None
    return None, f"no exact rational ray on support {list(dbl)}"


# ----------------------------------------------------------------------------
# level system in a-space
# ----------------------------------------------------------------------------


def base_rows(n, d, ranges=True):
    """rows (coeffs, rhs) over unknowns a_2..a_d meaning coeffs . a <= rhs."""
    nv = d - 1
    rows = []
    for m in range(2, n + 1):
        row = [Fr(comb(m, k)) for k in range(2, d + 1)]
        rows.append((tuple(-z for z in row), Fr(m)))        # p~(m) >= 0
        rows.append((tuple(row), Fr(1 - m)))                 # p~(m) <= 1
    if ranges:
        def box(idx, lo, hi):
            """Box.
            
            Args:
                idx:
                lo:
                hi:
            
            """
            r1 = [Fr(0)] * nv
            r1[idx] = Fr(1)
            rows.append((tuple(r1), Fr(hi)))
            r2 = [Fr(-1) if i == idx else Fr(0) for i in range(nv)]
            rows.append((tuple(r2), Fr(-lo)))
        if d >= 2:
            box(0, -2, -1)      # a_2
        if d >= 3:
            box(1, 0, 4)        # a_3
        if d >= 4:
            box(2, -14, 9)      # a_4
        if d >= 5:
            box(3, -80, 86)     # a_5
    return rows


def wvec(n, d, m) -> tuple:
    """w_m = m + sum a_k C(m,k) as (constant, coeffs)."""
    return Fr(m), tuple(Fr(comb(m, k)) for k in range(2, d + 1))


def enumerate_vertices(rows, nv) -> list:
    """All vertices of {A x <= b} (assumed bounded): exact, complete."""
    verts = {}
    R = len(rows)
    for sub in itertools.combinations(range(R), nv):
        mat = [rows[i][0] for i in sub]
        rhs = [rows[i][1] for i in sub]
        x = exact_solve(mat, rhs)
        if x is None:
            continue
        if check_rows(rows, x):
            verts[x] = True
    return list(verts.keys())


def val_of_func(vt, cfun):
    """Val of func.
    
    Args:
        vt:
        cfun:
    
    Returns:
        The computed result
    
    """
    c0, cv = cfun
    return c0 + sum(a * b for a, b in zip(cv, vt))


def geom(n, d, ranges=True, out=None) -> tuple:
    """Geom.
    
    Args:
        n:
        d:
        ranges (bool):
        out:
    
    Returns:
        tuple: Result of type tuple
    
    """
    def pr(*a):
        """Pr.
        
        """
        s = " ".join(str(x) for x in a)
        print(s)
        if out is not None:
            out.append(s)

    nv = d - 1
    rows = base_rows(n, d, ranges)
    pr("=" * 72)
    pr(f"GEOM n={n} d={d} ranges={ranges}  (#rows={len(rows)})")

    # boundedness diagnostics + numeric feasibility
    bounded = True
    for i in range(nv):
        for sgn in (1, -1):
            c = [Fr(0)] * nv
            c[i] = Fr(sgn)
            st, val, _ = linprog_dir(rows, nv, c)
            if st != 0:
                bounded = False
                pr("  numeric LP status %s for dir %+d*a_%d" % (st, sgn, i + 2))
    feas, xd = linprog_feasible(rows, nv)
    pr(f"  numeric feasibility: {feas} ; bounded(dir-scan): {bounded}")

    if feas is False:
        y, err = farkas_numeric_support(rows, nv)
        if y is None:
            pr("  NUMERIC INFEASIBLE, no exact certificate found: " + str(err))
            return None
        nz = [(i, y[i]) for i in range(len(y)) if y[i] != 0]
        pr("  EXACT FARKAS CERTIFICATE of INFEASIBILITY: %d nonzero rows" %
           len(nz))
        pr("  b.y =", sum(rows[i][1] * y[i] for i in range(len(y))))
        return ("infeasible", y)

    if not bounded:
        pr("  WARNING: unboundedness detected; vertex enumeration unsafe.")
        return None

    verts = enumerate_vertices(rows, nv)
    pr(f"  vertices: {len(verts)}")
    # dimension check: affine rank of vertices
    rank = 0
    if verts:
        # affine rank via Gaussian elimination on differences
        base = verts[0]
        pts = [[v[i] - base[i] for i in range(nv)] for v in verts[1:]]
        mat = [row for row in pts]
        rk = 0
        for col in range(nv):
            piv = None
            for r in range(rk, len(mat)):
                if mat[r][col] != 0:
                    piv = r
                    break
            if piv is not None:
                mat[rk], mat[piv] = mat[piv], mat[rk]
                for r, item in enumerate(mat):
                    if r != rk and item[col] != 0:
                        f = item[col] / mat[rk][col]
                        for k in range(col, nv):
                            item[k] -= f * mat[rk][k]
                rk += 1
        rank = rk
    pr("  affine dimension (rank of vertex hull): %d (full = %d)" %
       (rank, nv))
    if len(verts) <= 8:
        for v in verts:
            pr(f'    vertex a = ({", ".join(...)})')

    # exact min/max of each a_k and w_m over the polytope (= over vertices)
    pr("  exact coordinate ranges over polytope:")
    for i in range(nv):
        vals = [v[i] for v in verts]
        pr("    a_%d in [%s, %s]   (width %.6g)" %
           (i + 2, min(vals), max(vals), float(max(vals) - min(vals))))
    forced = {}
    for m in range(2, n + 1):
        vals = [val_of_func(v, wvec(n, d, m)) for v in verts]
        lo, hi = min(vals), max(vals)
        tag = ""
        if lo == hi:
            tag = "   FORCED EXACT"
            forced[m] = lo
        pr("    w_%02d in [%s, %s] = [%.6f, %.6f]%s" %
           (m, lo, hi, float(lo), float(hi), tag))

    # exact exclusion certificates for forced 0/1 levels
    pr("  exact certificates that w_m can be neither below nor above each "
       "forced value:")
    for m, v in sorted(forced.items()):
        if v == 0:
            # exclude b_m >= 1  <=>  w_m >= 1/C(n,m)
            _, coef = wvec(n, d, m)
            row_rhs = Fr(m) - Fr(1, comb(n, m))
            rows2 = rows + [(tuple(-z for z in coef), row_rhs)]
            y, err = farkas_numeric_support(rows2, nv)
            pr("    w_%d == 0 certified: %s" %
               (m, "YES (Farkas)" if y is not None else "NO: " + str(err)))
        elif v == 1:
            # exclude b_m <= C(n,m)-1  <=>  w_m <= 1 - 1/C(n,m)
            _, coef = wvec(n, d, m)
            row_rhs = Fr(1 - m) - Fr(1, comb(n, m))
            rows2 = rows + [(tuple(coef), row_rhs)]
            y, err = farkas_numeric_support(rows2, nv)
            pr("    w_%d == 1 certified: %s" %
               (m, "YES (Farkas)" if y is not None else "NO: " + str(err)))
    ones = sorted(m for m, v in forced.items() if v == 1)
    zeros = sorted(m for m, v in forced.items() if v == 0)
    pr(f"  FORCED ONES: {ones}")
    pr(f"  FORCED ZEROS: {zeros}")
    return ("ok", verts, forced)


# ----------------------------------------------------------------------------
# restricted system integrality test (CP-SAT, exact integers)
# ----------------------------------------------------------------------------


def restrict_test(n, d, m0, forced, time_limit=600, out=None) -> tuple:
    """Restrict test.
    
    Args:
        n:
        d:
        m0:
        forced:
        time_limit (int):
        out:
    
    Returns:
        tuple: Result of type tuple
    
    """
    def pr(*a):
        """Pr.
        
        """
        s = " ".join(str(x) for x in a)
        print(s)
        if out is not None:
            out.append(s)

    from ortools.sat.python import cp_model
    pr("=" * 72)
    pr(f"RESTRICT TEST n={n} d={d} m0={m0} time={time_limit}s")

    # forced equalities for the m0-subset: (level m, value v), 2 <= m <= m0
    fq = []
    for m, v in sorted(forced.items()):
        if 2 <= m <= m0 and v in (0, 1):
            fq.append((m, int(v)))
    pr(f"  forced values on restriction: {fq}")

    # valid beta ranges
    def C(k):
        """C.
        
        Args:
            k:
        
        Returns:
            The computed result
        
        """
        return comb(m0, k)
    rng = {2: (-2 * C(2), -1 * C(2)), 3: (0, 4 * C(3)), 4: (-14 * C(4), 9 * C(4)),
           5: (-80 * C(5), 86 * C(5))}

    # scale by D = lcm of C(m0,k), k=2..d
    def lcm(a, b):
        """Lcm.
        
        Args:
            a:
            b:
        
        Returns:
            The computed result
        
        """
        from math import gcd
        return a * b // gcd(a, b)
    D = 1
    for k in range(2, d + 1):
        D = lcm(D, C(k) if k <= m0 else 1)
    pr(f"  scale D = {D}")

    model = cp_model.CpModel()
    beta = {}
    for k in range(2, d + 1):
        lo, hi = rng.get(k, (-10 ** 12, 10 ** 12))
        beta[k] = model.NewIntVar(int(lo), int(hi), f"beta_{k}")
    for m in range(2, m0 + 1):
        terms = []
        for k in range(2, min(d, m) + 1):
            coef = comb(m, k) * (D // C(k))
            terms.append(coef * beta[k])
        expr = m * D + sum(terms)
        model.Add(expr >= 0)
        model.Add(expr <= D)
    for (m, v) in fq:
        terms = []
        for k in range(2, min(d, m) + 1):
            coef = comb(m, k) * (D // C(k))
            terms.append(coef * beta[k])
        model.Add(m * D + sum(terms) == v * D)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_workers = 4
    st = solver.Solve(model)
    pr(f"  CP-SAT status: {solver.StatusName(st)}")
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        sol = [solver.Value(beta[k]) for k in range(2, d + 1)]
        pr(f"  solution beta = alpha' = {sol}")
        # cross-check rational feasibility exactly
        ok = True
        for m in range(2, m0 + 1):
            val = Fr(m) + sum(Fr(sol[k - 2]) * Fr(comb(m, k), C(k))
                              for k in range(2, min(d, m) + 1))
            if val < 0 or val > 1:
                ok = False
        pr(f"  exact check of found solution: {ok}")
        return ("SAT", sol)
    if st == cp_model.INFEASIBLE:
        pr("  => restriction integrality INFEASIBLE (computer proof; "
           "to be re-derived exactly)")
        return ("UNSAT", None)
    return ("UNKNOWN", None)


if __name__ == "__main__":
    args = sys.argv[1:]
    out = []
    try:
        cmd = args[0]
        if cmd == "geom":
            n, d = int(args[1]), int(args[2])
            ranges = ("noranges" not in args)
            geom(n, d, ranges, out)
        elif cmd == "test":
            n, d, m0 = int(args[1]), int(args[2]), int(args[3])
            ranges = ("noranges" not in args)
            tl = 600
            for a in args[4:
                ]:
                if a.isdigit():
                    tl = int(a)
            res = geom(n, d, ranges, out)
            if res is not None and res[0] == "ok":
                _, verts, forced = res
                restrict_test(n, d, m0, forced, tl, out)
            else:
                print("geom failed; no test")
    finally:
        fn = "amanotype_out.txt"
        with open(fn, "w", encoding="utf-8") as fh:
            fh.write("\n".join(out) + "\n")
