FINAL REPORT

A. Repository and baseline
- Repository: degree-sensitivity-boolean
- Branch: jules-6307820464392629802-b18bb73f
- Starting commit: Initial state at run start.
- Scope examined: Evaluated execution precision of FrontierMath bounds specifically `log(6)/log(3)` and the exact calculations of `$d^c$` in `experiments/level_lp.py`.

B. Improvements implemented
- Problem discovered: Using `math.log` to define `C_CRIT = log(6.0)/log(3.0)` caused floating-point underflow errors such as `$3^{\log_3 6}$` evaluating to `5.999999999999999`. A Cyrillic `с` homoglyph was also used for `dс` making it susceptible to typos. Because of the precision issue, `need n>=` was calculated as `6` instead of `7` when $d=3$.
- Root cause: Standard `float64` imprecision causing edge-case integer comparisons to fail.
- Solution implemented: Replaced `float64` math with arbitrary precision using Python's `decimal` module, configuring 60 digits of precision. Computed `C_CRIT = Decimal(6).ln() / Decimal(3).ln()` and extracted bounds safely. Cleaned up the homoglyph.
- Files or components affected: `experiments/level_lp.py`
- Expected and demonstrated benefit: Ensures absolute integrity for the Kushilevitz exponent calculation. Running it establishes that `need` correctly identifies the exact needed `n` without underestimating due to precision issues (e.g., correctly showing `KILLED` for $d=3$ with `need n>=7`).

C. Verification results

| Check or command | Actual outcome | Relevant output or evidence | Explanation of limitation |
| :--- | :--- | :--- | :--- |
| `python experiments/level_lp.py` | PASS | See `level_lp_out.txt` showing `7` for $d=3$. | Local verification. |
| `python tests/test_small.py` | PASS | `ALL PASS` for core tests. | Local verification. |
| `git status` | PASS | Only `experiments/level_lp.py` is modified and staged. | Local verification. |

D. Measurable impact
Before: `python experiments/level_lp.py` printed `3 6 6.000 6 OPEN` for $d=3$.
After: It prints `3 6 6.000 7 KILLED` for $d=3$, properly determining that $6 < 7$ kills the exponent viability for $n=6, d=3$.

E. Remaining risks
The other files inside `experiments` also calculate `C = Decimal(6).ln() / Decimal(3).ln()`, which might need 60 precisions explicitly set across the board (like in `scan_level_d.py` where it sets 80). Further inspection could synchronize the precisions identically if necessary.

F. Next-run priorities
1. Synchronize the `decimal.getcontext().prec` across all scripts defining `C_CRIT` to an agreed upon value (e.g. 60 or 80) in a centralized utility to prevent mismatches.
2. Extend `test_small.py` with rigorous mathematical tests testing the exponent explicitly with varying `n` and `d` variables to guarantee floating-point regressions do not occur in CI.
3. Review `amanotype.py` for similar precision bounds in `Fraction` usage.

G. Overall verdict
VERIFIED IMPROVEMENT
