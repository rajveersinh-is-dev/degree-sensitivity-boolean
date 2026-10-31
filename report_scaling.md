# FINAL REPORT

## Repository and baseline
- **Repository:** `degree-sensitivity-boolean`
- **Branch:** `jules-6954957818189371045-615a5c46`
- **Starting commit:** Initial state
- **Scope examined:** The mathematical modeling scripts in the `experiments/` directory, specifically focusing on linear programming checks which used `scipy.optimize.linprog` with the `highs` solver. The task objective explicitly noted that scaling issues with large matrices (e.g., d >= 15) caused `linprog` to fail with `Model error` and `res.x = None`.

## Improvements implemented
- **Problem discovered:** For large configuration sizes, `scipy.optimize.linprog` (using the `highs` method) encounters severe matrix scaling issues (with matrix values spanning many orders of magnitude up to `1e16`), which leads to HiGHS failing with status `2` (Model error) and `res.x` being `None`.
- **Root cause:** The coefficients in the linear programs (e.g., combinations and factorials in exact level constraint matrices) grow extremely fast, causing severe conditioning issues in the `A_ub` and `A_eq` matrices.
- **Solution implemented:** Normalization steps were added to all usages of `linprog` in `experiments/amanotype.py`, `experiments/exact_level.py`, `experiments/exact_level2.py`, `experiments/hierarchy_lp.py`, and `experiments/level_lp.py`. We compute the maximum absolute value for each row of the constraint matrix and scale both the row and the right-hand side (`b` or `beq`) before passing them to `linprog`. We also use defensive `.max(axis=1)` handling to prevent division by zero.
- **Files or components affected:**
  - `experiments/amanotype.py`
  - `experiments/exact_level.py`
  - `experiments/exact_level2.py`
  - `experiments/hierarchy_lp.py`
  - `experiments/level_lp.py`
- **Expected and demonstrated benefit:** The numerical LPs for larger configurations now solve optimally without throwing `Model error` (e.g., verified on `(83, 15)` and `(93, 16)` configurations directly and in tests like `scan_level_d.py`).

## Verification results
| Check or command | Actual outcome | Relevant output or evidence | Explanation |
| :--- | :--- | :--- | :--- |
| `python tests/test_small.py` | PASS | `ALL PASS` | All regression checks and tests on small configurations are preserved. |
| `python experiments/scan_level_d.py 16` | PASS | For `d=15` and `d=16` returned `feasible=True`, instead of crashing or returning `None`. | Confirmed that `linprog` effectively solves LP instances beyond `d=14` successfully. |

*Note: All verification was done locally.*

## Measurable impact
- The `linprog` queries for `d=15` (n=83) and `d=16` (n=93) now return valid feasibility outcomes rather than stalling at `Model error`. The script `experiments/scan_level_d.py` successfully completes up to `d=16` where it previously crashed or stalled.

## Remaining risks
- Since we normalized the matrices by row max absolute values, we rely on `scipy`'s interior tolerance. The system involves exact `Fraction` math derived from extreme magnitude combinations. If exact fractions demand extreme tolerance, numerical feasibility bounds might theoretically differ; however, this is a heuristic scan step before exact verification via certificates.

## Next-run priorities
1. Further investigate exact integrality and CP-SAT solvers handling large coefficient bounds.
2. Complete full independent verification and benchmarking of `hierarchy_milp.py` for `(14,5)`.
3. Check and robustify parsing of results in `ksearch2.py` against potential timeouts or unbounded variables.

## Overall verdict
**VERIFIED IMPROVEMENT**
