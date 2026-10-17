# Research journal (degree-sensitivity-boolean)

## 2026-10-17 — Target lock and reductions

- FM statement confirmed from `research/data/ds_full.txt`.
- Reduction chain derived: FM ⟺ ∃ d with a fully sensitive (at 0) degree-≤d
  function on \(n_{\min}(d)=\lfloor d^c\rfloor+1\) variables, \(c=\log6/\log3\).
- Restriction is downward-closed: \((n,d)\) example ⟹ \((n',d)\) example for
  \(d\le n'\le n\).  Hence \((14,5)\) exists ⟺ \((n,5)\) exists for ALL
  \(n\in[5,14]\) ⟺ \(N(5)\ge 14\).  Any UNSAT at \(n\le 13\) kills \((14,5)\).
- \(n_{\min}\) table: d=2..14 → 4,7,10,14,19,24,30,37,43,50,58,66,75.
  Known: N(2)=3, N(3)=6, N(4)=9, (10,4) impossible, (14,5) is the FM decider.

## d=2, d=3 kills (PROVED HERE, Farkas certificates)

- Level system at \((4,2)\) and \((7,3)\) is infeasible; exact Farkas certificates
  \(b\cdot y=-1\) (files: `experiments/scan_level_d.py`, `scan_all_n.py`).

## d=4: machine reproof of Amano (PROVED HERE)

- Level polytope at (10,4) is a single point:
  \((a_2,a_3,a_4)=(-1,7/12,-1/6)\), forced ones {2,8,9}, zeros {5,10}.
- Proskurin 2021 (arXiv 2101.08600 §5) independently derives the same unique
  symmetrization polynomial via LP (Theorems 10/11); consistent with Amano's
  Fact 5; Proskurin was apparently unaware of Amano's impossibility.
- `experiments/ksearch2.py` (exact level-count bounds added to `ksearch.py`)
  proves \((10,4)\) INFEASIBLE in 14.1 s (plain `ksearch.py`: UNKNOWN at 1200 s).
  This is a new, self-contained machine proof of Amano's Theorem 3.

## d=5: obstruction analysis (PROVED HERE)

- \((14,5)\) level polytope: 17 vertices, full-dimensional (4-dim), NO forced
  levels.  Ranges: \(a_2\in[-13/12,-1]\), \(w_2\in[11/12,1]\) ⟹ **every**
  fully sensitive degree-≤5 function on 14 vars has \(b_2\in[84,91]\) (at most 7
  zero-pairs).  Similarly \(b_3\in[174,238]\), \(b_4\in[108,351]\),
  \(b_5\in[0,520]\), etc. (`research/data/geom_14_5.txt`, `intscan_14_5.txt`).
- No per-slice forced levels at any \((m_0,5)\) for \(m_0=5..14\) (slice_bounds.py).
- Per-slice lower bounds on pair counts are weak and only strengthen locally.
- **Methodological theorem:** every convex/averaging relaxation of the level
  system is trivial: slice averages equal the global level profile, which always
  lies in the polytope.  Only per-slice lattice/integrality (Amano's route) can
  kill a case; at \((14,5)\) the LP is full-dimensional so that route has no
  forced inputs.

## New structural theorems (PROVED HERE)

- **OR-extension:** if \((k,d)\) exists then \((k+1,d+1)\) exists via
  \(f\mapsto f\vee x_{k+1}\).  Hence \(N(d+1)\ge N(d)+1\).
- **Extension equivalence:** \((k+1,d)\) exists ⟺ there is a \((k,d)\)-function
  \(g\) and an integer \(u\) with \(g\le u\le g+1\) pointwise, \(u(0)=0\),
  \(\deg u\le d-1\) (equivalently a Boolean co-factor \(r\) with
  \(\deg(g+r)\le d-1\)).  (`experiments/extchain.py`.)
- **Amano Fact 6** (bounded class) extended in our notes: bounded degree-5 on 14
  vars impossible by Markström–Talbot (max 126 edges, K4⁻-free 3-graphs).

## Computational frontier (COMPUTATIONALLY VERIFIED / UNRESOLVED)

- \((n,5)\) SAT for \(n=6,7,8,9,10\) (witnesses verified by `ksearch.verify`).
- \(N(5)\ge 10\) via OR-extension of the verified \((9,4)\) witness.
- \((11,5)\)–\((14,5)\): UNDECIDED.  Direct CP-SAT, exact level-bounds, slice
  bounds, extension-chain and LNS have not settled them; long runs in progress.
- \((14,5)\) direct search returned UNKNOWN at 14400 s (6 workers) with the
  plain model; the slice-bounds model is still running.

## Literature

- Amano, EATCS 101 (2010): Theorem 3 \(((10,4)\) impossible), Fact 6 (bounded
  \((14,5)\) impossible), Problem 4/5.
- Proskurin, arXiv 2101.08600 §5: Theorems 10/11 (symmetrization uniqueness).
- Markström–Talbot, J. Combin. Designs 18 (2009): K4⁻-free 3-graph extremal.
- Low-degree sensitivity paper arXiv 2107.11205: rotation-symmetric searches.

## Files

See `README.md`.
EOF
echo "journal written"