"""restrict.py -- restrict a (k,d) y-json to its first k_new coordinates.

Usage: python restrict.py <src.json> <k_new> <dst.json>
"""
import json
import sys

src, knew, dst = sys.argv[1], int(sys.argv[2]), sys.argv[3]
d = json.load(open(src))
out = {}
for key, v in d.items():
    idx = [int(x) for x in str(key).split(",")]
    if all(i < knew for i in idx):
        out[key] = v
json.dump(out, open(dst, "w"))
print("kept", len(out), "of", len(d), "keys")
