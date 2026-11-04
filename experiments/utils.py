import numpy as np
from scipy.optimize import linprog

def linprog_scaled(c, A_ub=None, b_ub=None, A_eq=None, b_eq=None, bounds=None, method="highs", **kwargs):
    """
    Wrapper for scipy.optimize.linprog that applies infinity-norm row scaling
    to the constraint matrices to prevent 'Model error' and `res.x = None`
    on large configurations. Result objects are guaranteed to have `.x`
    defensively set to `None` if it is missing or optimization failed.
    """

    scaled_A_ub = None
    scaled_b_ub = None
    if A_ub is not None and b_ub is not None:
        A_ub = np.asarray(A_ub)
        b_ub = np.asarray(b_ub)
        row_norms_ub = np.abs(A_ub).max(axis=1)
        row_norms_ub[row_norms_ub == 0] = 1.0
        scaled_A_ub = A_ub / row_norms_ub[:, np.newaxis]
        scaled_b_ub = b_ub / row_norms_ub

    scaled_A_eq = None
    scaled_b_eq = None
    if A_eq is not None and b_eq is not None:
        A_eq = np.asarray(A_eq)
        b_eq = np.asarray(b_eq)
        row_norms_eq = np.abs(A_eq).max(axis=1)
        row_norms_eq[row_norms_eq == 0] = 1.0
        scaled_A_eq = A_eq / row_norms_eq[:, np.newaxis]
        scaled_b_eq = b_eq / row_norms_eq

    res = linprog(c=c, A_ub=scaled_A_ub, b_ub=scaled_b_ub, A_eq=scaled_A_eq, b_eq=scaled_b_eq, bounds=bounds, method=method, **kwargs)

    # Defensive x handling
    if not hasattr(res, 'x'):
        res.x = None

    return res
