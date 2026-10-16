# 05 — Normal contributions and Table III: Eqs. (I.68)–(I.93)

These are the code checks for the notes chapter *"Normal contributions and Table III: Eqs. (I.68)-(I.93)"*
(Tasks 11 and 12 of the re-derivation of Kahya & Woodard, Phys. Rev. D 76, 124005 (2007), "Paper I").
They check identities that were first derived by hand. They do not replace the analytic derivation in the notes.

| Script | Task | What it checks | Runtime |
|---|---|---|---|
| `check_68_80.py` | 11 | Derivative-extraction identities (I.70)–(I.80); the delta-function content of $`I^2[fA'']`$ (its $`y^{1-D/2}`$ component equals $`c_f A_{\rm lead}`$); the power counting after (I.80) | about 18 s |
| `check_81_93.py` | 12 | Time-derivative identities (I.81)–(I.88), the d'Alembertian (I.90), the reductions (I.91)–(I.92), the identity (I.93), and (I.89) plus Table III with the true $`D=4`$ function $`A(y)`$ of (I.62) | about 20 s |

Requirements: Python 3.10 or later, `sympy` (tested with 1.14.0) and `mpmath`. Both scripts run offline.

## How to run

```bash
cd code/05_normal_terms_tableIII
python3 check_68_80.py
python3 check_81_93.py
```

The outputs we got are in `expected_output/check_68_80.txt` and `expected_output/check_81_93.txt`.

## Method

The general identities are tested with power-law functions, $`f=y^\alpha`$ and $`A=y^\beta`$, with generic exponents, so no integration constants come up. Each test uses explicit conformal coordinates and the exact
$`y(x;x')=aa'H^2\big(\|\vec x-\vec x'\|^2-(\eta-\eta')^2\big)`$ with $`a=-1/(H\eta)`$.
Both sides are compared at random rational points with 25 significant digits. The tests run at $`D=4`$, and at $`D=5`$ where stated. Some identities are checked exactly in symbolic form; the docstring of each script says which.
The check of Table III uses the true $`D=4`$ propagator function $`A(y)=\frac{H^2}{16\pi^2}\{4/y-2\ln(y/4)-1\}`$.

## What the output means

**`check_68_80.py`**

- (A)–(E): every line should print `True`. These lines cover (I.70) (all 16 components), (I.71)–(I.72), (I.74), (I.75)–(I.78) with naive derivatives, and (I.79)–(I.80).
- (F): the $`y^{1-D/2}`$ part of $`I^2[fA'']`$ is `2**(D - 2)*N*Y**(1 - D/2)*c*c0`. This equals $`c\,\mathcal N c_0 4^{D/2-1}y^{1-D/2}=c_f A_{\rm lead}`$, so the comparison prints `True`.
- (G): at $`D=4`$, only $`y^{2-D}`$ in $`I^2[AA'']`$ is flagged `DIVERGENT`, and only $`y^{1-D}`$ and $`y^{2-D}`$ in $`I[A'^2]`$.

**`check_81_93.py`**

- Every identity line should print `True`.
- In the (I.89)/Table III block, "derived(67)_nonlocal + Table III(paper)" is exactly 0.0 at all three points.
  "derived(67)" is the nonlocal term of (I.67) as re-derived in Task 10. It uses the *definition sign* fixed by (I.168), while the printed Table III uses the paper's convention.
  So the derived Table III is $`-1`$ times the printed one, which is the overall sign inherited from (I.47).
- "(89) with H^2 restored == -derived" prints `True` at all three points.
  "(89) as printed (no H^2) == -derived" prints `False` at all three points: the second term of the printed (I.89) is missing a factor $`H^2`$.
- The last three lines expand the coefficient functions in $`x=y/4`$ (paper convention). Constants vanish under $`\nabla^2`$; in the third line, which multiplies $`\nabla^4`$, the term linear in $`y`$ (here $`8x`$) also vanishes, because $`\nabla^4y=0`$. Dropping these gives the three printed rows of Table III: $`-12/x+24\ln x`$, $`4/x`$ and $`8\ln x-16x\ln x`$.

## Comparison with the original runs (task logs)

- `check_81_93.py`: the output is identical to the final run in the Task 12 log.
  In the *first* run in that log, $`H`$ was drawn as `randint(5,20)/10` and came out as $`H=1`$ at two of the three points. At those points the printed (I.89), which is missing $`H^2`$, matched by accident.
  The final version, which is the one included here, draws $`H`$ from $`\{7/10,13/10,3/2,17/10\}`$.
  The derived term and Table III do not depend on $`H`$ at fixed $`(\eta,\eta',\vec x,\vec x')`$, so the printed scales are the same in both runs.
- `check_68_80.py`: the output is identical to the final run in the Task 11 log, with one cosmetic difference. In line (G), two entries both have exponent $`0`$ at $`D=4`$ ($`y^{4-D}`$ and $`y^{2-D/2}`$). The original sort left their order to Python's set ordering, which can change from one run to the next.
  We added a secondary sort key, `str`, so the order is now deterministic. The list of powers and their classification are unchanged.
- No absolute paths. The only edits to the scripts are the docstrings, two inline comments, the sort key above, and the removal of an unused placeholder line in `check_81_93.py`. None of these changes the logic.
