FINAL REPORT

A. Repository and baseline
- Repository: degree-sensitivity-boolean
- Branch: jules-2608941107094228634-cab4211d
- Starting commit: Initial state at run start.
- Scope examined: Explored `experiments/amanotype.py` robustness when running the mathematical verification via `scipy.optimize.linprog` for higher dimensions, finding it throws exceptions instead of failing gracefully.

B. Improvements implemented
- Problem discovered: The `experiments/scan_level_d.py` script crashed for `d=15` due to `scipy.optimize.linprog` returning a HiGHS Model Error (`status = 2`). Because the mathematical formulation for higher levels results in combination sizes like `math.comb(83, 40)` that can create massively unscaled coefficients (e.g. `1e16`), `linprog` exits correctly without assigning `res.x`. Consequently, calculations expecting `res.x` to be a vector (such as `b @ res.x`) throw `ValueError`. In addition, `linprog_feasible` erroneously treated any `status == 2` result as definitively "infeasible" instead of "failed to solve numerically", preventing the script from correctly surfacing mathematically indeterminate states due to numerical constraints.
- Root cause: Missing safeguards around `scipy.optimize.linprog` return states, assuming that any status failure still contains a valid `.x` matrix, or assuming that `status == 2` universally maps to mathematical infeasibility when using the Highs solver.
- Solution implemented: Replaced direct accesses of `res.x` and `res.fun` in `linprog_feasible`, `linprog_dir`, and `farkas_numeric_support` in `experiments/amanotype.py` with defensive usage like `getattr(res, 'x', None)`. Implemented checks in `linprog_feasible` to inspect `res.message` to accurately distinguish a "Model error" (undecided) from true feasibility logic (returning `None, None`).
- Files or components affected: `experiments/amanotype.py`
- Expected and demonstrated benefit: Ensures absolute structural stability for large configuration `d` searches. Prevents script crashes when matrices become poorly scaled, properly surfacing `feasible=None` when the model errors out natively.

C. Verification results

| Check or command | Actual outcome | Relevant output or evidence | Explanation of limitation |
| :--- | :--- | :--- | :--- |
| `python experiments/scan_level_d.py 15` | PASS | Script evaluates through `d=15` returning `feasible=None` without exception. | Local verification. |
| `python tests/test_small.py` | PASS | `ALL PASS` for core tests. | Local verification. |
| `git status` | PASS | Removed testing artifacts. | Local verification. |

D. Measurable impact
Before: `python experiments/scan_level_d.py 15` crashed with `ValueError: matmul: Input operand 1 does not have enough dimensions`.
After: It successfully continues processing and correctly prints `d=15 n_min= 83  d^c=82.8175211537761530757  feasible=None`.

E. Remaining risks
The solver currently yields "Model error" due to scaling and cannot evaluate large degrees like `d=15`. It might be possible to rescale the matrix by powers of 10 prior to linear programming constraints without altering mathematical invariants, restoring explicit Farkas verification.

F. Next-run priorities
1. Research strategies for better `linprog` matrix scaling within `amanotype.py` to overcome HiGHS `Model error` limitations for values beyond `d=14`.
2. Determine if the Fraction-based exhaustive Farkas verification loops can replace numerical LP where feasible for specific constraints.

G. Overall verdict
VERIFIED IMPROVEMENT
