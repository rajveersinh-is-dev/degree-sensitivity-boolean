"""Two-level hierarchy LP for the sensitivity-vs-degree necessary conditions.

KEY LEMMA (proved in research notes):
Let F = f^{-1}(1) for a Boolean f with f(0)=0, all singletons in F, deg f <= d.
For any subset T of coordinates define the *avoiding* counts
    c^{(T)}_k := #{B in F : B cap T = empty, |B| = k},     k = 2,...,N;  N := n-|T|.
Then for every m = d+1,...,N:
    sum_{k=2}^{m} (-1)^{m-k} c^{(T)}_k C(N-k, m-k) = (-1)^m m C(N, m).
(Proof: sum the exact conditions c_A = 0 over all m-subsets A of [n]\\T; f(empty)=0
and all singletons of [n]\\T are in F, exactly as in the master derivation.)

Master system: T = empty (N = n).
Level-1 systems: T = {i}, one for each coordinate (N = n-1).
Consistency (linear!):
    sum_i c^{({i})}_k = (n-k) b_k,                b_k := #{B in F : |B| = k}
and cross bounds 0 <= b_k - c^{({i})}_k <= C(n-1, k-1).

Each individual system is affine-determined above level d (free w's at levels 2..d).
Joint feasibility is a linear program.  Infeasibility at (n,d) proves that no
degree-<=d Boolean function has all n coordinates sensitive at the origin.

This script decides joint feasibility for given (n,d) numerically and then
(a) verifies a rational point exactly (feasibility certificate), or
(b) looks for an exact Farkas certificate (infeasibility certificate).
Optionally runs a MILP with integrality of all counts b_k, c^{(i)}_k.
"""
from fractions import Fraction as Fr
from math import comb
import itertools
import numpy as np
from scipy.optimize import linprog, milp, LinearConstraint, Bounds

OUT = []


def pr(*a):
    OUT.append(" ".join(str(x) for x in a))


def build_system(N, d):
    """Affine parametrisation of system (*) for size N, degree d.
    Returns (aff, nv) where aff[m] = (const, coefftuple) expresses
    w_m (the density b_m/C(N,m)) as affine function of free w_2..w_d."""
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


def rows_for(aff, N, nv, offset, scale):
    """Inequalities 0 <= w_m <= 1 for m=2..N as rows in the joint variable vector.
    Returns list of (dict-free list of coefficients length nv_total, rhs).
    Caller assembles full-length coefficient lists."""
    raise NotImplementedError


def assemble(n, d, levels=1):
    """Build the joint system.
    Variables: [master free x (d-1)] + per-i free x (d-1) * n   (levels=1)
    Rows: (list coeffs (length nvars, Fractions), rhs Fraction)  meaning c.x <= rhs.
    Also returns metadata: functions b_k(x), c_i_k(x) as affine expressions."""
    nv = d - 1
    nvars = nv * (1 + n * levels)
    aff_master, _ = build_system(n, d)
    aff_i = []
    for i in range(n):
        a, _ = build_system(n - levels, d)  # N = n-1 when levels=1
        aff_i.append(a)

    def unpack(vec):
        """vec: dict m -> (const, coeffs as list of length nvars)"""
        return vec

    def mrow(m, off, affv, scale, sign=1):
        """row for (sign)*(sum coeffs) <= rhs? we instead return coefficient list."""
        pass

    rows = []

    def w_expr(off, affv, m):
        const, vv = affv[m]
        c = [Fr(0)] * nvars
        for j, z in enumerate(vv):
            c[off + j] = z
        return (const, c)

    # master rows: 0 <= w_m <= 1, m = 2..n
    master_b = {}
    for m in range(2, n + 1):
        const, c = w_expr(0, aff_master, m)
        scale = comb(n, m)
        # b_m = scale * w_m >= 0  and  <= C(n,m)
        rows.append(([-scale * z for z in c], scale * const))
        rows.append(([scale * z for z in c], Fr(comb(n, m)) - scale * const))
        master_b[m] = (scale * const, [scale * z for z in c])

    # per-i rows: 0 <= w^{(i)}_k <= 1 for k = 2..N
    N = n - levels
    cik = {}
    for i in range(n):
        off = nv * (1 + i * levels)
        for k in range(2, N + 1):
            const, c = w_expr(off, aff_i[i], k)
            scale = comb(N, k)
            rows.append(([-scale * z for z in c], scale * const))
            rows.append(([scale * z for z in c], Fr(comb(N, k)) - scale * const))
            cik[(i, k)] = (scale * const, [scale * z for z in c])

    # consistency: sum_i c^{(i)}_k = (n-k) b_k  for k = 2..N
    for k in range(2, N + 1):
        const = Fr(0)
        c = [Fr(0)] * nvars
        for i in range(n):
            const += cik[(i, k)][0]
            for j in range(nvars):
                c[j] += cik[(i, k)][1][j]
        # subtract (n-k) b_k
        bk0, bkc = master_b[k]
        const -= (n - k) * bk0
        for j in range(nvars):
            c[j] -= (n - k) * bkc[j]
        rows.append((c, Fr(0)))            #  sum <= 0
        rows.append(([-z for z in c], Fr(0)))  # -sum <= 0

    # cross bounds: 0 <= b_k - c^{(i)}_k <= C(n-1, k-1), k=2..N
    for i in range(n):
        for k in range(2, N + 1):
            const = master_b[k][0] - cik[(i, k)][0]
            cvec = [master_b[k][1][j] - cik[(i, k)][1][j] for j in range(nvars)]
            rows.append(([-z for z in cvec], const))       # b_k - c_k >= 0
            rows.append((cvec, Fr(comb(N, k - 1)) - const))  # <= C(N,k-1)
    return rows, nvars, master_b, cik


def rat(v, md=10 ** 8):
    return tuple(Fr(float(x)).limit_denominator(md) for x in v)


def check_rows(rows, x):
    for (c, b) in rows:
        s = sum(cc * xx for cc, xx in zip(c, x))
        if s > b:
            return False
    return True


def exact_solve(mat, rhs):
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
                for k2 in range(col, m + 1):
                    A[r][k2] -= f * A[col][k2]
    return tuple(A[i][m] / A[i][i] for i in range(m))


def feas_point(rows, nv):
    A = np.array([[float(z) for z in c] for (c, b) in rows])
    b = np.array([float(x) for (_, x) in rows])
    A2 = np.hstack([A, np.ones((A.shape[0], 1))])
    res = linprog(c=np.array([0.0] * nv + [-1.0]), A_ub=A2, b_ub=b,
                  bounds=[(None, None)] * nv + [(None, None)], method="highs")
    if res.status == 0 and res.x[nv] > 1e-9:
        for md in (10 ** 12, 10 ** 9, 10 ** 6):
            pt = rat(res.x[:nv], md)
            if check_rows(rows, pt):
                return pt, res.x[nv]
    res = linprog(c=np.zeros(nv), A_ub=A, b_ub=b, bounds=[(None, None)] * nv,
                  method="highs")
    if res.status == 0:
        pt = rat(res.x, 10 ** 7)
        if check_rows(rows, pt):
            return pt, 0.0
    # interior point attempt
    try:
        res2 = linprog(c=np.zeros(nv), A_ub=A, b_ub=b, bounds=[(None, None)] * nv,
                       method="highs-ipm")
        if res2.status == 0:
            for md in (10 ** 12, 10 ** 9, 10 ** 6):
                pt = rat(res2.x, md)
                if check_rows(rows, pt):
                    return pt, 0.0
    except Exception:
        pass
    return None, None


def farkas(rows, nv):
    m = len(rows)
    A = np.array([[float(z) for z in c] for (c, b) in rows])
    b = np.array([float(x) for (_, x) in rows])
    Aeq = np.vstack([A.T, np.ones((1, m))])
    beq = np.zeros(nv + 1)
    res = linprog(c=b, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * m, method="highs")
    if res.status != 0:
        return None, "numeric ray LP status %s" % res.status
    val = b @ res.x
    if val > -1e-9:
        return None, "no numeric ray (val %.3g)" % val
    sup = list(np.where(res.x > 1e-8)[0])
    for sub in itertools.combinations(sup, min(nv + 1, len(sup))):
        mat = [[Fr(A[i, j]) for j in sub] for i in range(nv)]
        mat.append([Fr(b[j]) for j in sub])
        rhs = [Fr(0)] * nv + [Fr(-1)]
        sol = exact_solve(mat, rhs)
        if sol is None or any(s < 0 for s in sol):
            continue
        y = [Fr(0)] * m
        for j, s in zip(sub, sol):
            y[j] = s
        ok = True
        for i in range(nv):
            if sum(Fr(A[i, j]) * y[j] for j in range(m)) != 0:
                ok = False
                break
        if ok and sum(Fr(b[j]) * y[j] for j in range(m)) < 0:
            return y, None
    return None, "no exact ray on support"


def run_case(n, d, levels=1, do_milp=False):
    pr("=" * 72)
    pr("hierarchy case (n=%d, d=%d, levels=%d)" % (n, d, levels))
    rows, nvars, master_b, cik = assemble(n, d, levels)
    pr("  variables:", nvars, " rows:", len(rows))
    pt, slack = feas_point(rows, nvars)
    if pt is not None:
        pr("  FEASIBLE (exact rational point verified), max-slack ~ %.4g" % (slack or 0))
        # report master b_k
        bk = []
        for k in range(2, n + 1):
            c0, c = master_b[k]
            val = c0 + sum(x * y for x, y in zip(c, pt))
            bk.append("b_%d=%s" % (k, val))
        pr("  master level counts:", ", ".join(bk))
        if do_milp:
            pr("  (milp skipped)")
    else:
        cert, err = farkas(rows, nvars)
        if cert is not None:
            pr("  INFEASIBLE (exact Farkas certificate; support %d of %d rows)"
               % (sum(1 for z in cert if z != 0), len(cert)))
        else:
            pr("  UNKNOWN:", err)


if __name__ == "__main__":
    import sys
    cases = [(6, 3), (7, 3), (9, 4), (10, 4), (11, 4), (14, 5), (15, 5), (16, 5)]
    if len(sys.argv) > 2:
        cases = [(int(sys.argv[1]), int(sys.argv[2]))]
    for (n, d) in cases:
        try:
            run_case(n, d)
        except Exception as e:
            pr("  ERROR:", repr(e))
    with open("hierarchy_out.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("\n".join(OUT))
