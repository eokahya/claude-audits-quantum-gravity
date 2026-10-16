# 09 — Schwinger–Keldysh mode equation and late-time mode function (Paper II)

Scripts for **Task 18** (notes chapter "Schwinger-Keldysh mode equation and the late-time mode function (Paper II)").
They insert the corrected Table 1 of Paper II (= corrected Table X of Paper I) into the Schwinger–Keldysh (SK)
mode equation (II.29), re-derive Tables 2–5 of Paper II, and redo the late-time analysis without the
simplification (II.44).

All scripts use units H = 1 and represent ln a by the symbol `L_a` (`extract.py` prints it in LaTeX form as `\ln a`,
and k stands for k/H). Unless stated otherwise they work in the
**Table-X (paper) convention**, i.e. with the sign of the printed tables; the physical self-mass-squared in the
**definition sign** of Paper I, Eq. (I.168), is −1 times this (the notes quote physical results in the definition sign).

## Requirements and running

Python ≥ 3.10 with `sympy` and `mpmath` (see `../requirements.txt`). No network access is needed.

```sh
cd code/09_schwinger_keldysh
python3 tables234.py          # or any other script below
sh run_all.sh                 # runs everything and rewrites expected_output/ (about 2 minutes)
```

## Files

| File | What it does | Runtime |
|---|---|---|
| `engine.py` | Helper module. Exact closed-form evaluation of the SK integrals: after a' = a s every integral becomes ∫ s^α (1−s)^β {1, ln s, ln(1−s)} ds, done with harmonic numbers and −π²/6; result is an exact expansion in 1/a. Also the operator 𝒟 = −(a²∂²_a + 4a∂_a) (= □/H² at k = 0) and truncation in powers of a. | import |
| `rows.py` | Helper module. The 13 rows of the corrected Table 1 (constant of row 7: −54 → +18; ζ row: 10/3 → 34/3, −24 → −40) and the printed values. | import |
| `tables234.py` | Tables 2–4 of Paper II: the 12 integrals S_K[f], K = 3, 4, 5, f = 1/x, ln x/x, ln x, ln²x; compares with the printed entries and prints the next order. | ~1 s |
| `numcheck.py` | Engine versus direct mpmath quadrature (30 digits) at a = 50 and 400 for seven integrals. | ~4 s |
| `table5.py` | Table 5 of Paper II: the 19 entries (rows without ∇², simplification (II.44)) and their total. | ~1 s |
| `kexp.py N` | Beyond (II.44): exact u(η',k), sin(kr)/(kr), the k²/a² part of □ and all ∇² rows, expanded exactly to order k^N; prints the non-falling part row by row and in total. | N=2: 2 s, 4: 6 s, 6: 10 s, 10: 27 s |
| `numcheck_k.py` | The k⁸-truncated series of `kexp.py` versus quadrature with the exact k dependence (a = 200, k/H = 0.3, 1). | ~5 s |
| `analysis2.py` | Printed versus corrected table through k⁶; contribution of the ζ correction; corrected total through O(1/a). | ~10 s |
| `analysis3.py` | Coefficients k^n a^−m of the full right-hand side for the corrected and the printed table. | ~10 s |
| `extract.py` | The numbers quoted in the notes (definition sign): Table 5″ per row, late-time form of R(η,k) through k¹⁰, 1/a tails of the printed table. | ~35 s |
| `run_all.sh` | Runs all scripts and writes `expected_output/`. | ~2 min |

## Expected results (see `expected_output/`)

- `tables234.txt`: `paper - derived = 0` for all 12 entries of Tables 2–4.
- `numcheck.txt`: |exact expansion − quadrature| ≤ 1.02e-19 in all 14 cases (largest for q=3, ln²x, ν=2, m=1, a=50).
- `table5.txt`: `paper-derived: 0` for all 19 entries; `TOTAL (non-falling): 0`.
- `kexp_N2/N4/N6/N10.txt`: non-falling total 0 at every order k⁰ … k^N. Rows 10 and 11 separately contain
  ±(4/15) k² a (they cancel); row 12 (corrected ζ row) gives +1/70 at k⁴ (Table-X convention).
- `numcheck_k.txt`: differences 1e-10 – 2e-8 at k/H = 0.3 and 5e-6 – 1e-3 at k/H = 1, the size of the k⁹ truncation error.
- `analysis2.txt`: printed table non-falling part `[0, 0, 0, 0, 0, 0, 0]`; ζ correction
  k⁴·4(a²−15)/(45a³) + k⁶·2(6a²−35)/(525a³) (only falling terms); corrected total has no 1/a term.
- `analysis3.txt`: corrected total starts at a⁻² (k⁰: −12, k²: (24 L_a + 103)/3, …, Table-X convention); the printed
  table additionally has k⁴ a⁻¹: −4/45 and k⁶ a⁻¹: −4/175.
- `extract.txt`: definition-sign coefficients; non-falling totals `[0, …, 0]` (k⁰ … k¹⁰) for both tables;
  a⁻² coefficients 12, 0, −8 ln a − 103/3, −8i, 8/3, 32i/45, −2/15, …; a⁻³ only −12i at k³.

## Agreement with the original run

All scripts were re-run for the notes (Python 3.13, sympy 1.14, mpmath 1.3). Every output is **identical,
line by line**, to the output recorded in the original Task 18 command log (for `kexp.py 6` and `kexp.py 10`,
whose full output the log stored only in files, to those original files).

## Changes with respect to the original scripts

- Absolute scratch paths removed; an unused helper (`Bz`) removed from `engine.py`.
- The original scripts passed intermediate results through pickle files (`S_tables.pkl`, `analysis2.pkl`,
  `kexp_rows_N10.pkl`). Each script now recomputes what it needs, so every script runs on its own.
- `kexp.py` was wrapped in a class `KExpansion(nmax)` so that `numcheck_k.py`, `analysis2.py` and `extract.py`
  can choose the expansion order without manipulating `sys.argv`. The mathematics is unchanged.
- A first attempt (`sk_core.py`, `check_tables234.py`), which used sympy's generic `integrate` followed by a
  series expansion in 1/a, failed (no asymptotic expansion of the resulting polylogarithms) and was replaced by
  `engine.py`; it is not included.
