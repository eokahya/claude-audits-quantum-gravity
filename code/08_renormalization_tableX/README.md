# 08 — Renormalization and the final self-mass (Paper I, Eqs. (I.146)–(I.167), Tables VI–X)

Scripts used in Tasks 16 and 17 of the re-derivation of

- **Paper I**: E. O. Kahya and R. P. Woodard, *Quantum gravity corrections to the one loop scalar self-mass during inflation*, Phys. Rev. D 76, 124005 (2007).

The full derivation is in the notes PDF (`../../notes/notes.pdf`), chapter "Renormalization and the final self-mass: Tables VI-X, Eqs. (I.146)-(I.167)". A summary page is in `../../docs/08-renormalized-self-mass-table-X.md`.

## Requirements

Python 3.10+ and SymPy (tested with SymPy 1.14.0, mpmath 1.3.0, Python 3.13). No network access, no other files.

## How to run

```bash
cd code/08_renormalization_tableX
python3 check_146_IX.py         # ~3 s
python3 check_160_167.py        # ~2 s
python3 check_tableX_global.py  # ~17 s (imports check_global.py from this folder)
python3 check_global.py         # ~14 s (optional: the helper's own Task 13 self-test)
```

The output of each script is stored in `expected_output/<script>.txt`.

## Scripts

| Script | Task | What it verifies | Runtime |
|---|---|---|---|
| `check_146_IX.py` | 16 | (A) pointwise parts and residues of (148)–(150) and the delta coefficient of (147); (B) the "add zero" step (151) and the D→4 expansions (152), (153); (C) flat Euclidean analogue of (151), exact in D; (D) Table VI from Table IV for general D; (E) Table VII; (F) the kernel identities (154), (155) and six simpler a′→a conversions on a test function; (G) Tables VIII and IX, rows and totals; (H) (156)–(159); (I) flat-space loop numerator V·P·V = p²k²/2 | ~3 s |
| `check_160_167.py` | 17 | (A) Table X assembled from Tables III, V, VII and (159); (B) pointwise identities used in the global test; (C) β₁, β₂, β₃ of (161)–(165) for general D and at D = 4, with the corrected and definition-sign values | ~2 s |
| `check_tableX_global.py` | 17 | Global pointwise test at D = 4: direct sum of the normal contributions P₁…P₇ plus the nonlocal term of (67) versus the 13 operators of Table X, at three points with 30 digits; printed table, with the ζ-row fix, and with ζ-row fix plus 72γ₂ | ~17 s |
| `check_global.py` | 13 (helper) | Helper module providing `run(...)`, which builds P₁…P₇ by brute-force differentiation and the ten external operators acting on the appendix totals. Run on its own, it repeats the Task 13 generic-power self-test (diff = 0.0) | ~14 s |

## What the outputs mean

**`check_146_IX.txt`**: every pass/fail check prints `True` or `0`; the other lines print the derived expressions (residues, the rows and totals of Tables VI–IX, the (157) coefficient). In section (D) the line "with Table IV's printed gamma2 sign ..." shows the γ₂ row you would get from the printed sign of Table IV. That row is the negative of the one in the paper's Table VI, so the paper used the correct sign in its calculation and the Table IV sign is a typo. In section (G), "Table IX total = ... value at D=4: -16/3" is the forbidden H⁴ a^D δ^D divergence that the printed γ₂ sign would produce. With the correct sign the total carries a factor (D−4).

**`check_160_167.txt`**:
- (A): all 13 lines read `assembled - paper = 0`. The assembled ζ row deliberately leaves out the Table III term, and the last line prints that term separately: `-16*x*log(x) + 8*log(x)`. It is the contribution the printed Table X omits.
- (B): three zeros. At D = 4 the nonlocal parts of (152) and (153) reproduce (4/y)² and (4/y)³ pointwise, and the two terms of (159) cancel pointwise. This is why (159) does not enter the global test.
- (C): reproduces −61/3 and −4/(D−4) + 58/3 + 2γ (`e` stands for D−4). It also gives the corrected values −49/3 and −4/(D−4) + 52/3 + 2γ (paper convention) and the definition-sign values, which are −1 times these.

**`check_tableX_global.txt`**: one line per point. The differences "lhs − TableX(printed)" (−11.5, 64.1, −3.53) and "+ zeta fix" (1.18, 65.3, 0.016) are nonzero. "+ zeta + 72 gamma2" gives about 1e-31, i.e. zero at 30-digit precision. So Table X is reproduced only when **both** corrections are applied.

**`check_global.txt`**: four lines with `diff = 0.0`. This is the Task 13 self-test of the helper.

## Agreement with the original runs

All four scripts were re-run in this folder after cleaning (docstrings added, unused variables and unused imports, including `mpmath`, removed). The outputs are **identical character for character** to the outputs recorded in the original command logs of the tasks (Task 16 and Task 17; for `check_global.py`, Task 13). These logs are working files and are not part of this repository.

During Task 16, the first two runs of `check_146_IX.py` failed because of three script bugs. All three were fixed before the final run. The fixes changed only the script, not the results:
1. An integer had no `.subs` method.
2. `H` was not declared positive, so SymPy could not simplify (−1/(Hη))^D. This produced spurious `False` results for (154), (155) and □ ln a = −(D−1)H².
3. The coefficient of the formal symbol `Box_H2` was read from a factored expression instead of an expanded one. This dropped the β□/H² and γ₁□/H² rows and gave spurious mismatches of the Table VIII/IX totals.

The first version of check (C) also used numerical quadrature. It lost accuracy as D → 4 because the integrand r^(3−D) approaches non-integrability: the difference was 2e-16 at D = 3.5, 1.3e-5 at D = 3.8 and 7.5 at D = 3.95. It was replaced by the exact Gamma-function evaluation that is in the script now.
