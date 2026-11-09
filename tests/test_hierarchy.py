import sys
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from hierarchy_milp import milp_feasible

def main():
    cases = [(9, 4, 0), (10, 4, 0), (11, 4, 2), (14, 5, 0), (15, 5, 0), (16, 5, 2)]
    ok = True
    for n, d, expected in cases:
        st, _ = milp_feasible(n, d, master_only=False, time_limit=120)
        good = (st == expected)
        ok = ok and good
        print("%s (n=%d,d=%d): got %s, expected %s" % ("PASS" if good else "FAIL", n, d, st, expected))

    print("ALL PASS" if ok else "SOME FAILED")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
