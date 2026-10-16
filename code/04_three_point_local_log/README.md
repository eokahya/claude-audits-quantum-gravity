# 04 — 3-point local and logarithm terms, Eqs. (I.54)–(I.67)

Code checks for the chapter "3-point local and logarithm terms: Eqs. (I.54)–(I.67)" of the
derivation notes (`../../notes/notes.pdf`). Here "Paper I" is E. O. Kahya and R. P. Woodard,
Phys. Rev. D 76, 124005 (2007).

The scripts support the analytic derivation in the notes. They do not replace it. Each script
prints `True` or `False` for every identity it tests, and prints explicit coefficients where
the notes compare numbers with the paper.

| Script | Task | Verifies | Runtime |
|---|---|---|---|
| `check_54_59.py` | 9 | Local 3-point terms (I.54)–(I.59): (A) the local polynomial taken straight from (I.45) and Table I for D = 4, 5, 6; (B) the outer-derivative rules on $`a^{D-2}\delta^D`$ and (I.55)–(I.56); (C) assembly of (I.59): finite parts $`-2`$ (for $`a^D\square`$) and $`+5/2`$ (for $`a^{D-2}\nabla^2`$) in the definition sign. Also checks that the paper's own (I.54)+(I.57)+(I.58) give the finite part $`-5/2`$ for the $`a^{D-2}\nabla^2`$ coefficient, while (I.59) prints $`-7/2`$, and that the paper's terms are $`-1`$ times the derived ones. | ~1 s |
| `check_60_67.py` | 10 | Logarithm terms (I.60)–(I.67): (a) (I.62) from (I.14) at $`D=4+\epsilon`$; (b) $`\square A=(D-1)kH^2`$ at $`D=4`$; (c) the delta-function lemma (D5) and helper identities for symbolic D; (d) a pointwise 30-digit check of the non-local parts of the four structures $`\mathcal V_1,\mathcal W,\mathcal V_3,\mathcal V_4`$ at three random points, using the actual $`A(y)`$; (e) the general-D assembly and the cancellation against the log-type local terms of (I.44), which leaves a factor $`(D-4)(D+1)/[4(D-3)]`$; (f) the paper's (I.63)–(I.66) are $`-1`$ times the derived terms and add up to the printed (I.67). | ~8 s |

## How to run

```bash
cd code/04_three_point_local_log
python3 check_54_59.py
python3 check_60_67.py
```

Requirements: Python 3.10+ and SymPy (see `../requirements.txt`; tested with Python 3.13.1, SymPy 1.14.0,
mpmath 1.3.0).
No network access is needed and no files are read or written.

## Expected output

`expected_output/check_54_59.txt` and `expected_output/check_60_67.txt` hold the outputs of
the runs made while preparing this folder. They are identical, character for character, to the
outputs recorded in the original working logs of Tasks 9 and 10. The random points of
`check_60_67.py`, section (d), come from a fixed seed, so the output is reproducible.

How to read the output:

- `check_54_59.py`, section (C): `derived coeff` gives the exact-D coefficients of
  $`a^D\square\,\delta^D`$ and $`a^{D-2}\nabla^2\delta^D`$ in units of $`i\kappa^2 k`$, in the
  definition sign of Paper I's (I.168). `pihat` stands for
  $`\hat\pi=-\pi\cot(\pi D/2)=-2/(D-4)+O(D-4)`$. The line `finite parts at D=4` of section (C)
  gives the finite parts at D = 4. The paper's printed (I.59) is in the paper's own convention, which is $`-1`$ times the
  definition sign.
- `check_60_67.py`, section (d): the numbers in the first list are the absolute differences
  between the brute-force derivatives and the derived non-local expressions. They should be
  0.0. The "scales" list gives the absolute value of each derived expression at that point, for
  comparison.
- `check_60_67.py`, section (e): `LB`, `LL`, `M`, `Q4`, `T`, `Tp` stand for
  $`\ln(aa')\,a^D\square\,i\delta^D`$, $`\ln(aa')\,a^{D-2}\nabla^2 i\delta^D`$, $`H^2a^D i\delta^D`$,
  $`(aa')^{D-2}\ln(aa')\nabla^4A`$, $`H(aa')^{D-2}a\partial_0\nabla^2A`$ and
  $`H(aa')^{D-2}a'\partial_0'\nabla^2A`$, in units of $`\kappa^2 k`$. Section (f) uses the same
  symbols in units of $`\kappa^2H^2/(8\pi^2)`$ at $`D=4`$.
