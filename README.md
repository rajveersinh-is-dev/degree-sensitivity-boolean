# Degree Sensitivity of Boolean Functions (FrontierMath `degree-sensitivity-boolean`)

This repository is an autonomous research programme on the open problem

> Construct a multilinear polynomial \(P\) with Boolean values such that
> \(\sum_i |P(0)-P(e_i)| > \deg(P)^{\omega}\) for some \(\omega > \log 6/\log 3 \approx 1.63093\).

Equivalently (FrontierMath, Epoch AI): find \(d\) and a Boolean function on
\(n \ge \lfloor d^{\log 6/\log 3}\rfloor+1\) variables that is fully sensitive at
the origin and has degree \(\le d\).  This is the same as the open question of
whether the Kushilevitz \(F_3\)-chain exponent \(\log 6/\log 3\) can be beaten.

## Main result (honest, labelled)

- **PROVED HERE (machine-verified, 14.1 s):** no fully sensitive degree-\(\le 4\)
  function exists on 10 variables.  This is an independent computer proof of
  Amano's Theorem 3 (2010), obtained by exact level-count propagation: the
  level polytope at \((10,4)\) is a single point, forcing
  \(b=(45,70,35,0,35,70,45,10,0)\); the resulting CP-SAT model
  (`experiments/ksearch2.py`) is INFEASIBLE in 14.1 s.  Plain `ksearch.py`
  could not decide this case in 1200 s.
- **PROVED HERE (elementary):** \(N(d+1)\ge N(d)+1\) via the OR-extension
  \(f \mapsto f\vee x_{n+1}\) (degree rises by exactly one).  With
  \(N(4)=9\) this gives the new lower bound \(N(5)\ge 10\) (explicit witness
  `experiments/ksol_or_10_5.json`, independently verified).
- **PROVED HERE (new structural theorem):** the extension equivalence
  \((k+1,d)\) exists \(\iff\) there is a \((k,d)\)-function \(g\) and an integer
  \(u\) with \(g\le u\le g+1\) pointwise, \(u(0)=0\), \(\deg u\le d-1\)
  (`experiments/extchain.py`).  Equivalently: a Boolean co-factor \(r\) with
  \(\deg(g+r)\le d-1\).
- **PROVED HERE (exact level system):** at \((n,d)=(14,5)\) the level polytope is
  full-dimensional (17 vertices, no forced levels), but the level counts are
  tightly constrained, e.g. **every** fully sensitive degree-\(\le 5\) function on
  14 variables has at most 7 zero-pairs: \(b_2\in[84,91]\).  Farkas certificates
  kill \((4,2)\) and \((7,3)\).
- **PROVED HERE (methodological):** every convex/averaging relaxation of the
  level system is trivial: slice averages equal the global level profile, which
  always lies in the polytope.  Only per-slice lattice/integrality structure
  (Amano's route) can kill a case.
- **COMPUTATIONALLY VERIFIED (CP-SAT):** \((n,5)\) is SAT for
  \(n=6,7,8,9,10\) (witnesses in `experiments/ksol*.json`); all pre-gate
  checks (Boolean, \(f(0)=0\), \(f(e_i)=1\), degree \(\le 5\)) pass
  independently.
- **UNRESOLVED:** whether \((11,5),\dots,(14,5)\) exist.  Our searches (direct
  CP-SAT, exact level-bounds, slice bounds, extension-chain, LNS) have not
  decided them; several long runs are still in progress.  If any \((n,5)\) with
  \(n\le 13\) is UNSAT, then \((14,5)\) is impossible and the \(d=5\) case of
  FrontierMath is dead.
- **CONJECTURED (not proved):** \(N(5)=10\) and the FM inequality fails for
  \(d=5\) (i.e. no \((14,5)\) function).  This is consistent with, but not
  established by, the obstruction analysis.

Final verdict category (Section 26): **B** — substantial new structural results
and a new machine proof of a known theorem, but the open \((14,5)\) question
remains undecided.

## Repository layout

- `experiments/` — exact search code and witnesses.
  - `ksearch.py` — base CP-SAT model + independent verifier.
  - `ksearch2.py` — adds exact level-count bounds (kills \((10,4)\) in 14 s).
  - `extchain.py` — exact extension/co-factor model (CP-SAT).
  - `extsat.py` — SAT (Cadical) fixed-\(g\) extension search.
  - `localsearch.py` — simulated annealing on the \(r\)-formulation.
  - `amanotype.py` — level-polytope vertex enumeration + Farkas certificates.
  - `orinv.py`, `restrict.py`, `slice_bounds.py`.
- `research/` — problem statement text, polytope data, literature extraction.
- `paper/` — LaTeX write-up (in progress).
- `tests/` — reproducibility checks.
- `report/` — final honest report.

## Reproducing the headline

```bash
# (10,4) is infeasible (Amano) -- our machine proof
python -u experiments/ksearch2.py 10 4 600 val 4
# -> INFEASIBLE, ~14 s

# small cases are SAT (verified)
python -u experiments/ksearch2.py 8 5 60 chk 2
```

## Labels

Every claim in `research/` and `paper/` is labelled PROVED HERE / ESTABLISHED IN
PRIOR WORK / COMPUTATIONALLY VERIFIED / CONJECTURED / UNRESOLVED.  No theorem is
fabricated; computations that returned UNKNOWN are reported as UNKNOWN.

## License

MIT.
EOF
echo "README written"