# Code for Chapter 6: Tables IV and V, Eqs. (I.94)–(I.125) (Task 13)

These scripts are the independent symbolic and numerical checks used in Task 13. They check, but do not replace, the analytic derivation in the notes (chapter "Tables IV and V: Eqs. (I.94)-(I.125)" of `../../notes/notes.pdf`).

Requirements: Python 3.10+ and sympy (tested with Python 3.13.1, sympy 1.14.0). No network access is needed.

Run each script from this folder:

```
cd code/06_tables_IV_V
python3 check_94_125.py
python3 check_global.py
python3 check_tables45.py
python3 check_tables45b.py
python3 check_table4.py
python3 diag_leading_pair.py
python3 superseded_test_E.py
```

`check_tables45.py` and `check_tables45b.py` import `run()` from `check_global.py` in this folder.

## Scripts

| Script | Label in notes | What it verifies | Runtime |
|---|---|---|---|
| `check_94_125.py` | (A)–(D) | (A) Eqs. (I.94)–(I.98), (I.100), (I.101), (I.102)–(I.103), pointwise in explicit coordinates at D = 4 and D = 5 with test functions A = y^β, B = y^γ, f = y^α, 30-digit arithmetic. (B) (I.99), (I.108), (I.112). (C) Exact-D coefficients of (4/y)^D and (4/y)^(D-1) in the eight products (I.115), and the D = 4 closed forms (I.116)–(I.123). (D) (I.124)–(I.125). | 93 s |
| `check_global.py` | (F4a,b) | The ten external operators acting on the totals of the appendix Tables XI–XX (as printed, with five typo corrections) reproduce the direct sum of the parts P1–P7. Uses generic powers A = y^(-7/3), B = y^(-5/4), ΔC = y^(2/7) at D = 5 and D = 4. | 14 s |
| `check_tables45.py` | (F4c), (T5) | The same global check with the actual D = 4 propagators. Then, for each operator, the derived D = 4 function minus the printed [Table IV(D = 4) + Table V], and the printed tables evaluated against the direct sum. | 21 s |
| `check_tables45b.py` | (T5) | Hypothesis tests: the printed tables; the γ2 sign of Table IV corrected; plus −12x in the δ row; plus +72 in the γ2 constant; only −12x. Also prints the derived D = 4 operator functions. | 31 s |
| `check_table4.py` | (T4) | Table IV in general D from the appendix totals and the exact-D leading coefficients. | 2 s |
| `diag_leading_pair.py` | diagnostic | A first version of part (C) extracted coefficients from the full series by matching exponents, and printed wrong coefficients for the Δ-type products. This run prints the leading-pair coefficient −8 of A'ΔC' at (4/y)^(D-1) and the six individual exponent-(D-1) terms of A'(C−A)' from the full series. The six terms (their sum is not printed) add up to −4(D−3)/(D−2), i.e. to −8 after the factor 2(D−2)/(D−3), the value of (I.122). | 2 s |
| `superseded_test_E.py` | superseded (E) | The first global test, from the first version of `check_94_125.py`. It used the printed tables and the direct sum with the derived signs of Task 7. Neither "direct + tables" nor "direct − tables" vanishes. It was replaced by the structured checks above and is kept for the record. | 16 s |

Runtimes were measured on an Apple-silicon Mac. No script takes longer than about 1.5 minutes, so all of them were re-run for this repository.

## What the outputs mean

- `check_94_125.py`: every entry in parts (A) and (B) prints `True`. Part (C) prints the coefficients listed in the table of leading coefficients in the notes. Coefficients with a 1/(D−4) pole (AA'', A''ΔB, A''ΔC at (4/y)^(D-1)) combine with the constant −π cot(πD/2) term into the D = 4 logarithms. Part (D) prints the constant (24 ln 2 − 36) ln 2, which α = (aa')^D □² annihilates, and `[True, True, True, True]`.
- `check_global.py`, and the first block of `check_tables45.py`: `diff = 0.0` at every point. This means zero to the 30-digit working precision.
- `check_tables45.py`, per-operator block: the differences are constants, convention-dependent ln 2 terms, `4/(3*x**2)` for γ2 (the Table IV sign), and `-12*x` in δ. The printed tables miss the direct sum by 0.0356, 17.05 and 0.607 at the three points.
- `check_tables45b.py`: only the two combined corrections (γ2 sign plus −12x in δ, or equivalently γ2 sign plus a γ2 constant of +54) give `[0.0, 0.0, 0.0]`.
- `check_table4.py`: every line ends in `True True` except γ2, which ends in `True False`. Our value is +(D³−3D²−4D+8)/(4(D−1)); the paper prints −(D³−3D²−4D+8)/(4(D−1)).

All tables are in the paper's sign convention, the sign of P1–P7, i.e. of the printed (I.47)–(I.49).

## Expected output and reproducibility

`expected_output/<script>.txt` holds the output of each script, re-run for this repository. For `check_94_125.py`, `check_global.py`, `check_tables45.py`, `check_tables45b.py` and `check_table4.py` it is byte-for-byte identical to the output recorded in the original Task 13 log. The outputs of `diag_leading_pair.py` and `superseded_test_E.py` match the corresponding inline runs in that log. The random test points come from fixed seeds.
