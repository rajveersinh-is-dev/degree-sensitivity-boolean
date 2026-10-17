"""orinv.py -- OR-extension: from a found (k,d) fully sensitive function g,
construct f = g OR x_{k+1}, which is a (k+1, d+1) fully sensitive function
(as a polynomial: f = g + y(1-g), degree = d+1).

Verifies with ksearch.verify and saves.

Usage: python orinv.py <input.json> <k> <d>
   e.g.: python orinv.py experiments/ksol_9_4.json 9 4   ->  (10,5)
"""
import json
import os
import sys
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extchain
import ksearch


def or_extend(k, d, yvals):
    g = extchain.g_values_from_y(k, d, yvals)
    dp = d + 1
    ysol = {}
    for m in range(2, dp + 1):
        for B in combinations(range(k + 1), m):
            if k in B:
                ysol[B] = 1
            else:
                ysol[B] = g[B]
    return ysol


if __name__ == "__main__":
    src = sys.argv[1]
    k = int(sys.argv[2])
    d = int(sys.argv[3])
    with open(src) as fh:
        yraw = json.load(fh)
    ysol = or_extend(k, d, yraw)
    ok, err = ksearch.verify(k + 1, d + 1, ysol)
    print("OR-extension (%d,%d) -> (%d,%d): verify %s"
          % (k, d, k + 1, d + 1, "OK" if ok else err))
    fn = "experiments/ksol_or_%d_%d.json" % (k + 1, d + 1)
    with open(fn, "w") as fh:
        json.dump({",".join(map(str, sorted(B))): v for B, v in ysol.items()}, fh)
    print("saved", fn)
