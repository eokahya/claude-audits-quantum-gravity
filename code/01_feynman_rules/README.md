# 01 — Lagrangian expansion and Feynman rules (Tasks 1–2)

Computer-algebra checks for the first chapter of the notes
("Lagrangian expansion and Feynman rules: Eqs. (I.6)–(I.8), Tables I–II").
Paper I is E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

The scripts are **checks**, not the derivation itself. The derivation is done
by hand in the notes (`../../notes/notes.pdf`, chapter 1). Both scripts were run
in the original working session; they were reconstructed from the command logs
of Tasks 1 and 2 (for `check_eq6.py`, the final version, after the in-session
edit that renamed one print label and removed an unused line). The versions
here have been cleaned up: the original header comment was replaced by a
docstring and the code was reformatted (spacing, one import per line). The
original scripts contained no absolute paths and read or write no files. The
computations and the printed output are unchanged.

## Requirements

- Python 3.10 or later
- sympy (tested with 1.14.0)

No network access and no other packages are needed.

## Scripts

| Script | Checks | Method | Runtime |
|---|---|---|---|
| `check_eq6.py` | Eq. (I.6). It also checks (I.4) and (I.5) implicitly, and confirms that the truncated inverse metric is correct up to terms of order κ³. | Computes the exact inverse metric (adjugate/determinant) and the exact `sqrt(-det g)` for `g = eta + kappa h`, Taylor-expands both to κ², and compares the result with the derived bracket of (I.6). Cases: symbolic `h` for D = 2, 3; random rational `h` for D = 4, 5, 6; the conformal case `h = c eta` for symbolic D; and the residual `g g^{-1} - 1 = kappa^3 h h h` for D = 4. | ~22 s |
| `check_vertices.py` | Eqs. (I.7), (I.8), Table I and Table II | (A) Builds the full vertices from (I.6) by polarization with distinct legs, for D = 3, 4. (B) Compares the tables with the **exact** metric, without using (I.6), for D = 2, 3 (symbolic) and D = 4, 5 (random rational). (C) Fixes the table coefficients uniquely from the exact metric alone. The sign of `U_2` can be switched to compare the derived value `+i/4` with the printed `-i/4`. | ~10 s |

## How to run

```bash
cd code/01_feynman_rules
python3 check_eq6.py
python3 check_vertices.py
```

The outputs recorded in `expected_output/` were produced by these exact files
(Python 3.13.1, sympy 1.14.0, macOS):

- `expected_output/check_eq6.txt`
- `expected_output/check_vertices.txt`

Both files agree line for line with the outputs in the original command logs
of Tasks 1 and 2.

## What the output means

**`check_eq6.py`**: every check line should end in `True`. The line
`order-2 coefficient of exact: c**2*(D - 4)*(D - 2)/8` is the κ² coefficient
of `(1 + kappa c)^{(D-2)/2}` (conformal case). The derived bracket of (I.6)
reproduces it for general D.

**`check_vertices.py`**:

- **(A)** `True` / `0` means the derived Table I and Table II give the full
  vertex exactly. For the full 4-point vertex this is `U[3,4] + U[4,3]`.
  With the printed `U_2 = -i/4` the script reports a non-zero residual
  `I*(p0*q0 - p1*q1 - ...)*(a00*b00 - 2*a01*b01 + ...)`. In covariant form
  this is `-i (p1.p2)(h3^{rs} h4_{rs})`, which is exactly the sign flip of
  `U_2`.
- **(B)** All `True` means the derived tables rebuild the exact
  O(κ) and O(κ²) Lagrangian. The residual `paper - exact` is
  `-(1/4)(dphi.dphi)(h^{rs}h_{rs})`, without the factor κ² a^{D-2}. For
  D = 4, 5 it appears as a rational number times `p.p`, because `h` there is
  a random rational matrix.
- **(C)** The unique solution is `x = 1, y = -1/2, u = (-1/8, 1/4, 1/2, -1)`,
  i.e. `V_1 = i kappa a^{D-2} d1^a d2^b`, `V_2 = -(i/2) kappa a^{D-2} eta^{ab} d1.d2`,
  `U_I = i kappa^2 a^{D-2} u_I x (structure_I)`. The printed Table II has
  `u_2 = -1/4` and fails this test.

## Findings these scripts support

- Eqs. (I.4), (I.5) and (I.6) are reproduced exactly, for general D and with
  no gauge condition.
- Table I, Eq. (I.8) and rows 1, 3 and 4 of Table II are reproduced exactly.
- The left-hand side of (I.7) is printed with `-1/2 kappa`; the correct factor
  is `+1/2 kappa`. This is a misprint: `V_1` itself is correct.
- Table II, row 2 is printed as `-i/4`; the correct value is `+i/4`. Task 3
  later showed that the paper's own (I.33) uses `+i/4`, so this too is a
  misprint only.
