"""Analyse the Kushilevitz function h on 6 variables: sensitivity 6, degree 3.

h(z1..z6) = sum z_i - sum_{i<j} z_i z_j + (z1z3z4 + z1z2z5 + z1z4z5 + z2z3z4 + z2z3z5
                                          + z1z2z6 + z1z3z6 + z2z4z6 + z3z5z6 + z4z5z6)
"""
from itertools import combinations, permutations

TRIPLES = [(1, 3, 4), (1, 2, 5), (1, 4, 5), (2, 3, 4), (2, 3, 5),
           (1, 2, 6), (1, 3, 6), (2, 4, 6), (3, 5, 6), (4, 5, 6)]


def h(S):
    S = set(S)
    val = sum(1 for i in S)
    val -= sum(1 for i, j in combinations(sorted(S), 2))
    val += sum(1 for t in TRIPLES if set(t) <= S)
    return val


def main():
    F = []
    for k in range(0, 7):
        for S in combinations(range(1, 7), k):
            if h(S) == 1:
                F.append(S)
    print("|F| =", len(F), " of 64")
    from collections import Counter
    print("sizes:", sorted(Counter(len(s) for s in F).items()))
    print("F =", F)
    # check degree: compute Moebius transform
    Fset = set(F)
    def moebius(A):
        A = set(A)
        tot = 0
        for k in range(0, len(A) + 1):
            for B in combinations(sorted(A), k):
                if set(B) in Fset:
                    tot += (-1) ** (len(A) - k) * 1
        return tot
    deg = max(len(A) for k in range(0, 7) for A in combinations(range(1, 7), k)
              if (moebius(A) != 0))
    print("degree =", deg)
    print("Moebius on sets of size 2:", set(moebius(A) for A in combinations(range(1, 7), 2)))
    print("Moebius on sets of size 3:", set(moebius(A) for A in combinations(range(1, 7), 3)))
    # automorphisms
    autos = []
    for p in permutations(range(1, 7)):
        m = {i: p[i - 1] for i in range(1, 7)}
        if all(tuple(sorted(m[x] for x in s)) in Fset for s in F):
            autos.append(p)
    print("|Aut| =", len(autos))
    print("autos:", autos[:24])
    # complements?
    print("complement family size:", sum(1 for s in F if tuple(sorted(set(range(1, 7)) - set(s))) in Fset))
    # F as upward/downward closed?
    up = all(tuple(sorted(set(s) | {x})) in Fset for s in F for x in range(1, 7))
    print("upward closed:", up)


if __name__ == "__main__":
    main()
