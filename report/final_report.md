# Final report — degree-sensitivity-boolean

## 1. Problem and reductions

FrontierMath open problem `degree-sensitivity-boolean`: construct a multilinear
Boolean \(P\) with \(\sum_i |P(0)-P(e_i)| = \deg(P)^\omega\) for
\(\omega>\log6/\log3\).  Equivalently: find \(d\) with a fully sensitive (at 0)
degree-\(\le d\) function on \(n_{\min}(d)=\lfloor d^c\rfloor+1\) variables,
\(c=\log6/\log3\).  Restriction is downward-closed, so \((14,5)\) exists iff
\((n,5)\) exists for every \(5\le n\le 14\).

## 2. Results proved here

1. **Farkas kills** of \((4,2)\) and \((7,3)\) with exact certificates.
2. **Machine proof of Amano's \((10,4)\) impossibility** in 14.1 s via
   `experiments/ksearch2.py` (plain search: UNKNOWN at 1200 s).
3. **OR-extension theorem:** \(N(d+1)\ge N(d)+1\); new lower bound
   \(N(5)\ge 10\) with verified witness `ksol_or_10_5.json`.
4. **Extension equivalence theorem:** \((k+1,d)\) exists iff there is a
   \((k,d)\)-function \(g\) and integer \(u\) with \(g\le u\le g+1\),
   \(u(0)=0\), \(\deg u\le d-1\).
5. **\((14,5)\) obstruction analysis:** level polytope full-dimensional (17
   vertices, no forced levels); \(b_2\in[84,91]\) (≤7 zero-pairs); slice bounds;
   convex/averaging relaxations provably trivial.

## 3. Computational frontier

- \((n,5)\) SAT for \(n=6,7,8,9,10\) (witnesses verified).
- \((11,5)\)–\((14,5)\): UNDECIDED (long runs in progress; direct searches
  returned UNKNOWN).

## 4. Honest limitations

- The \((14,5)\) question is not decided; no claim of solving FrontierMath.
- The computational methods are weak at SAT-finding on these instances; this is
  reported, not hidden.
- Category (Section 26): **B**.

## 5. Section 27 statement

I declare that this report is an honest account of the work performed: every
theorem above is either proved in the repository, reproduced from cited prior
work, or labelled as conjectured/unresolved; every computation is reported with
its actual outcome (including UNKNOWN); no result, citation, commit, or push
status has been fabricated; and the open status of \((14,5)\) is stated plainly
rather than overstated.
EOF
echo "report written"