FINAL REPORT

A. Repository and baseline
- Repository: degree-sensitivity-boolean
- Branch: bugfix/linprog-res-x
- Starting commit: Current head
- Resulting commit: To be created
- Scope examined: `experiments/*.py` scripts calling `scipy.optimize.linprog` and `milp`.

B. Improvements implemented
- Problem discovered: Using `scipy.optimize.linprog` (especially the HIGHS solver) and `scipy.optimize.milp` could crash the scripts with a `TypeError` when checking `res.status != 0` because it implicitly assumed `res.x` was populated even on solver failures (e.g. `Model error`), but the result could have no `x` vector.
- Root cause: Missing defensive checks in multiple scripts when evaluating `res.x`.
- Solution implemented: Added `getattr(res, 'x', None) is not None` checks to conditional branching when utilizing `res.x` across five critical files: `exact_level2.py`, `exact_level.py`, `hierarchy_lp.py`, `hierarchy_milp.py`, and `level_lp.py`. Also safely retrieved `b @ res.x` by defaulting to `float('nan')` if `x` was `None`.
- Files or components affected:
  - `experiments/exact_level.py`
  - `experiments/exact_level2.py`
  - `experiments/hierarchy_lp.py`
  - `experiments/hierarchy_milp.py`
  - `experiments/level_lp.py`
- Expected and demonstrated benefit: Prevents scripts from crashing due to unexpected `OptimizeResult` states during mathematical experiments on large inputs.

C. Verification results
| Check or command | Actual outcome | Relevant output | Explanation |
|---|---|---|---|
| `python tests/test_small.py` | PASS | ALL PASS | Tests executed successfully, validating behavior remains sound. |
| `python experiments/level_lp.py` | PASS | (Valid Output) | The LP experiments run effectively. |
| `python experiments/hierarchy_milp.py` | PASS | (Valid Output) | The MILP experiments run effectively without failing in previously crashing edges. |
| `python experiments/exact_level.py 7 3` | PASS | (Valid Output) | Manual execution passes cleanly. |
| `python experiments/exact_level2.py` | PASS | (Valid Output) | Evaluated up to n=19, terminating execution without errors. |

D. Measurable impact
- Scripts no longer raise arbitrary `TypeError` unhandled exceptions under solver failures, enabling broader, fault-tolerant scans over parametric spaces to establish Boolean function degree sensitivities.

E. Remaining risks
- The repository executes complex continuous and integral relaxations natively on large problem configurations. Numerical instability in edge cases can still arise due to floating-point truncation (although the solver issues exact rational warnings where necessary).

F. Next-run priorities
- Consolidate common linear programming wrapper calls in a shared `utils.py` to prevent logic duplication.
- Scale `scan_level_d.py` into larger spaces (`d > 6`) with a robust timeout mechanism via `pytest`.
- Increase `hierarchy_milp.py` search depths using multiprocessing or distributed computing if further validation demands it.

G. Overall verdict
VERIFIED IMPROVEMENT
