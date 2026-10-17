"""Slice-level bounds: for each m0 <= 14, compute the exact (m0,5) level-polytope
bounds on the slice level counts.  Any forced slice value is a theorem about
restrictions of (possibly) fully sensitive degree-5 functions on 14 vars."""
import os
import sys
from math import comb

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "experiments"))
import ksearch2

print("(m0,5) level-polytope: b_m ranges for m <= 5 and forced levels")
for m0 in range(5, 15):
    b, v = ksearch2.level_bounds(m0, 5, verbose=False)
    row = []
    for m in range(2, min(6, m0) + 1):
        lo, hi = b[m]
        tag = "*" if lo == hi else ""
        row.append("b_%d[%d,%d]%s" % (m, lo, hi, tag))
    forced = [m for m in b if b[m][0] == b[m][1]]
    print("m0=%2d verts=%4d | %s | forced: %s"
          % (m0, len(v), " ".join(row), forced))
