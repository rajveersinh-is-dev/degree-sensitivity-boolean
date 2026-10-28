"""Reproducibility checks for the small (n,d) cases.

Run:  python tests/test_small.py
Expected:
  (6,5)  SAT   (verified)
  (7,5)  SAT   (verified)
  (8,5)  SAT   (verified)
  (10,4) INFEASIBLE (Amano; our machine proof)
"""
import json
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable


def run(n, d, tl=300, workers=2):
    cmd = [PY, "-u", os.path.join(ROOT, "experiments", "ksearch2.py"),
           str(n), str(d), str(tl), "test", str(workers)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=tl + 120,
                       cwd=ROOT)
    out = r.stdout + r.stderr
    status = None
    for line in out.splitlines():
        if line.startswith("n=%d d=%d" % (n, d)):
            status = line.split("status=")[1].split()[0]
    return status, out


def main():
    cases = [(6, 5, "OPTIMAL"), (7, 5, "OPTIMAL"), (8, 5, "OPTIMAL"),
             (10, 4, "INFEASIBLE")]
    ok = True
    for n, d, expected in cases:
        st, out = run(n, d)
        good = (st == expected)

        # If we expect a solution, ensure it passes independent verification (Boolean/degree checks)
        if expected == "OPTIMAL" and good:
            verif_ok = False
            for line in out.splitlines():
                if line.startswith("independent verification:"):
                    if "OK" in line:
                        verif_ok = True
                    break
            if not verif_ok:
                good = False
                print("FAIL verification for (n=%d,d=%d). Output:" % (n, d))
                print(out[-500:])

        ok = ok and good
        print("%s (n=%d,d=%d): got %s, expected %s %s"
              % ("PASS" if good else "FAIL", n, d, st, expected,
                 "" if good else "\n" + out[-500:]))

    # Mathematical tests for exponent bounds (issue robustness)
    sys.path.insert(0, os.path.join(ROOT, "experiments"))
    from exponent_math import pow_gt

    # 3^c is exactly 6.000 (when using precise bounds), so pow_gt(6, 3) must be False, and pow_gt(7, 3) True
    math_cases = [
        (6, 3, False),
        (7, 3, True),
        (4, 2, True),
        (3, 2, False),
        (10, 4, True),
        (9, 4, False)
    ]
    math_ok = True
    for n, d, expected in math_cases:
        res = pow_gt(n, d)
        if res != expected:
            print("FAIL pow_gt(n=%d, d=%d): got %s, expected %s" % (n, d, res, expected))
            math_ok = False
        else:
            print("PASS pow_gt(n=%d, d=%d)" % (n, d))

    ok = ok and math_ok

    print("ALL PASS" if ok else "SOME FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
