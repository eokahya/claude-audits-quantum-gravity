# 8. Renormalization and the final self-mass: Eqs. (I.146)–(I.167), Tables VI–X

<sub>[Overview (README)](../README.md) · Previous: [de Sitter symmetry and counterterms](07-de-sitter-symmetry-and-counterterms.md) · Next: [Schwinger–Keldysh and late times](09-schwinger-keldysh-late-time.md)</sub>

*Tasks 16 and 17 of the re-derivation of Paper I (Kahya and Woodard, Phys. Rev. D 76, 124005 (2007)).*

## What this step computes

This is the last stage of the calculation in Paper I. The "normal" 3-point contributions of Table IV contain the powers $`(4/y)^{D-1}`$ and $`(4/y)^{D-2}`$, which are not integrable in four dimensions. Pulling d'Alembertians out of them splits each one into a local divergent piece, proportional to $`\delta^D(x-x')`$, and a finite nonlocal piece built from $`\ln x/x`$, with $`x=y/4`$ (Task 16: Eqs. (I.146)–(I.159), Tables VI–IX).

All local pieces are then collected into three coefficients $`\beta_1,\beta_2,\beta_3`$, which the BPHZ counterterms absorb. All finite nonlocal pieces are collected into Table X. Together they give the renormalized self-mass-squared, Eq. (I.167) (Task 17).

We also tested Table X globally. At three pairs of points $`(x,x')`$ we compared it with a direct, brute-force sum of the normal contributions, with 30-digit precision.

Signs: all tables are in the paper's sign convention. The printed self-mass is $`-1`$ times the one defined by Paper I's own (I.168) (established in Tasks 3–13; both the 4-point and the 3-point parts are reversed, so their relative sign is consistent), so in that "definition sign" every entry is multiplied by $`-1`$.

IDs such as S2 or C17 refer to the [complete list of findings](10-complete-list-of-findings.md).

## Key steps

**1. Extracting d'Alembertians.** For functions of $`y`$, $`\square`$ acts as a differential operator in $`y`$ plus a delta function. The delta function comes only from the component $`y^{1-D/2}`$ (Eq. (I.147)). Eq. (I.149) writes $`(4/y)^{D-2}`$ as $`\square/H^2`$ acting on $`(4/y)^{D-3}`$. "Adding zero" in the form of the identity (I.150) for $`(4/y)^{D/2-1}`$, which carries the delta function, then moves the whole pole onto an explicit delta function (Eqs. (I.151)–(I.152)):

```math
\Big(\frac4y\Big)^{D-2}=\frac{i(4\pi)^{D/2}H^{-D}\,(D-2)}{(D-3)(D-4)\,\Gamma(\frac D2)}\,\frac{\delta^D(x-x')}{a^D}
-\frac{\square}{H^2}\Big\{\frac4y\ln\frac y4\Big\}+2\,\frac4y\ln\frac y4-\frac4y+O(D-4).
```

The power $`(4/y)^{D-1}`$ follows from (I.148) in the same way (Eq. (I.153)). The same mechanism was checked exactly in $`D`$ on the flat Euclidean analogue, using a Gaussian test function.

**2. Converting $`a'`$ into $`a`$.** An external operator acting on $`a^{-D}\delta^D(x-x')`$ turns $`a'^q`$ into $`a^q`$, but the factor ends up inside the $`\square`$. For example (Eq. (I.155)):

```math
(aa')^{D-1}(a^2+a'^2)H^2\square\big[a^{-D}\delta^D(x-x')\big]=2a^D\big(H^2\square-H^4\big)\delta^D(x-x').
```

**3. Two consistency checks (Tables VIII and IX).** The total coefficient of $`a^D\square^2\delta^D`$ vanishes. We confirmed independently that it must: the flat-space loop numerator reduces to $`\tfrac12p^2k^2`$, so the flat-space divergence is a scaleless integral and is zero. The total coefficient of $`a^DH^4\delta^D`$ carries a factor $`(D-4)`$:

```math
\text{Table IX total}=\frac{(D-4)(-D^4+5D^3-16D+16)}{4(D-2)}\qquad\text{(finite, as required by }\phi\to\phi+\text{const)}.
```

With the sign of the $`\gamma_2`$ entry printed in Table IV, the four $`D=4`$ values would add up to $`-16/3`$ instead of $`0`$, giving a divergence that shift symmetry forbids. This confirms the corrected ($`+`$) sign of Task 13 independently. The paper's own Tables VI–IX use that sign, so the Table IV sign is only a typo.

**4. The local coefficients.** In the paper's convention, in units of $`H^{D-4}/(4\pi)^{D/2}`$:

```math
\beta_1=0,\qquad \beta_2=-\frac{61}{3}+O(D-4),\qquad \beta_3=-\frac{4}{D-4}+\frac{58}{3}+2\gamma+O(D-4).
```

These reproduce Eqs. (I.161)–(I.165) exactly from the printed inputs. The two earlier finite errors, in Eq. (I.40) (B1) and in the arithmetic of Eq. (I.59) (B2), shift only the finite parts: $`\beta_2\to-\tfrac{49}{3}`$ and $`\beta_3\to-\tfrac{4}{D-4}+\tfrac{52}{3}+2\gamma`$ (paper convention; B3, B4). The pole $`-4/(D-4)`$ is unchanged. The counterterms absorb these shifts, and the only pole sits in the de Sitter-breaking $`\beta_3`$.

**5. Table X, row by row.** Assembling Tables III, V, VII and Eq. (I.159) reproduces 12 of the 13 operator rows. In the $`(aa')^2\nabla^4`$ row (the $`\zeta`$ row) the Table III term $`8\ln x-16x\ln x`$ was left out (S2):

```math
(aa')^2\nabla^4:\qquad \frac{34}{3}\ln x-40x\ln x+24x^2\ln x-36x^2
\qquad\Big(\text{printed: }\frac{10}{3}\ln x-24x\ln x+24x^2\ln x-36x^2\Big).
```

**6. Global test.** We compared the direct sum of the normal parts $`P_1,\dots,P_7`$, plus the nonlocal term of (I.67), with the 13 operators of Table X acting pointwise. The direct sum does not use the appendix's external-operator reduction, so the test checks the nonlocal entries of Table X (at $`x\neq x'`$) independently of that reduction. They agree to 30 digits only when two corrections are applied together: the $`(aa')^2\nabla^4`$ correction (S2), and the $`72\gamma_2`$ correction found in Task 13 (S3; in Table X, the constant of the $`(aa')^3(a^2+a'^2)H^4`$ row changes from $`-54`$ to $`+18`$). Neither correction alone is enough. The test is numerical, at three sample point pairs; it is not a symbolic identity, and it checks only the sum of the rows.

## Result

> **Renormalized one-loop self-mass-squared (definition sign of Eq. (I.168), $`D=4`$, $`x=y/4`$):**
>
> ```math
> -iM^2_{\rm ren}(x;x')=i\kappa^2a^4\Big(\Delta\alpha_1\square^2+\Delta\alpha_2H^2\square+\Delta\alpha_3H^2\frac{\nabla^2}{a^2}\Big)\delta^4(x-x')\;-\;\big[\text{corrected Table X}\big],
> ```
>
> Here the $`\Delta\alpha_i`$ are arbitrary finite constants. The corrected Table X equals the printed one (paper's convention) with three changes:
> - $`(aa')^2\nabla^4`$ becomes $`\tfrac{34}{3}\ln x-40x\ln x+24x^2\ln x-36x^2`$;
> - $`(aa')^3(a^2+a'^2)H^4`$ becomes $`\tfrac{4\ln x}{3x}-\tfrac{32}{3x}+18`$;
> - the duplicated $`(aa')^2(a^2+a'^2)\nabla^2\square`$ row is removed.
>
> The physically meaningful nonlocal part therefore has the opposite overall sign to the paper's (S1), and its $`(aa')^2\nabla^4`$ row contains an extra non-analytic term (S2). The $`72\gamma_2`$ change (S3) is analytic and drops out of every retarded quantity.
>
> **Late times.** The overall sign reverses the sign of the one-loop correction $`\Phi_1`$ to the mode function. At late times, however, the non-falling total of the source is zero ([Schwinger–Keldysh page](09-schwinger-keldysh-late-time.md)), so only decaying terms change sign and Paper II's late-time conclusion is unchanged. The $`\zeta`$-row correction adds only falling terms.

## Comparison with the paper

| Item | Paper I | This work | Type |
|---|---|---|---|
| Eqs. (I.146)–(I.159), Tables VI–IX | as printed | reproduced exactly (general $`D`$, and $`D\to4`$ where the paper takes the limit); Table VI with the corrected $`\gamma_2`$ sign | agreement |
| Table IV, $`\gamma_2`$ entry | sign $`-`$ | sign $`+`$; Tables VI–IX use $`+`$ (shift symmetry confirms it) | typo (C15) |
| Table VII header | coefficient of $`i\kappa^2H^4/(4\pi)^4`$ | coefficient of $`\kappa^2H^4/(4\pi)^4`$ (no $`i`$) | typo (C16) |
| Eqs. (I.160), (I.167) | $`\beta_2\square+\beta_3\nabla^2/a^2`$ | $`\beta_2H^2\square+\beta_3H^2\nabla^2/a^2`$ (same for $`\Delta\alpha_{2,3}`$) | typo (C17, C18) |
| $`\beta_1`$, $`\beta_2`$, $`\beta_3`$ (Eqs. (I.161)–(I.165)), paper convention, units $`H^{D-4}/(4\pi)^{D/2}`$ | $`0,\ -\tfrac{61}{3},\ -\tfrac{4}{D-4}+\tfrac{58}{3}+2\gamma`$ | reproduced; after fixing (I.40), (I.59): $`0,\ -\tfrac{49}{3},\ -\tfrac{4}{D-4}+\tfrac{52}{3}+2\gamma`$; pole unchanged | absorbable (B3, B4; finite parts only) |
| Table X, $`(aa')^2\nabla^4`$ row | $`\tfrac{10}{3}\ln x-24x\ln x+24x^2\ln x-36x^2`$ | $`\tfrac{34}{3}\ln x-40x\ln x+24x^2\ln x-36x^2`$ | substantive, non-analytic (S2) |
| Table X, $`(aa')^3(a^2+a'^2)H^4`$ row | constant $`-54`$ | constant $`+18`$ (from Table V, Task 13) | substantive, analytic (S3); drops out of the retarded kernel |
| Table X, $`(aa')^2(a^2+a'^2)\nabla^2\square`$ | printed twice | once | typo (C19) |
| Overall sign of $`-iM^2`$ | paper's convention | $`-1`$ times, by definition (I.168) | substantive (S1, inherited) |
| Eq. (I.170) (effect of the local terms on the mode equation) | $`\Delta\alpha_2\square+\Delta\alpha_3\nabla^2/a^2`$ | $`\Delta\alpha_2H^2\square+\Delta\alpha_3H^2\nabla^2/a^2`$; result $`-\kappa^2a^4\Delta\alpha_3H^2k^2/a^2\,\Phi_0`$ | typo (C22; noticed while refereeing the write-up, confirmed on the printed page) |

## Code

All scripts are in [`../code/08_renormalization_tableX/`](../code/08_renormalization_tableX/). They need Python 3.10+ and SymPy. Their outputs are stored in [`expected_output/`](../code/08_renormalization_tableX/expected_output/).

| Script | What it does | Runtime |
|---|---|---|
| [`check_146_IX.py`](../code/08_renormalization_tableX/check_146_IX.py) | Symbolic checks of Eqs. (I.147)–(I.159) and Tables VI–IX for general $`D`$, the flat-space analogue, and the flat-space loop numerator | ~3 s |
| [`check_160_167.py`](../code/08_renormalization_tableX/check_160_167.py) | Assembles Table X from Tables III, V, VII and (I.159); checks $`\beta_1,\beta_2,\beta_3`$ for general $`D`$ and at $`D=4`$, with the corrected and definition-sign values | ~2 s |
| [`check_tableX_global.py`](../code/08_renormalization_tableX/check_tableX_global.py) | Global pointwise test of Table X against the direct sum of the normal contributions (30 digits, 3 point pairs) | ~17 s |
| [`check_global.py`](../code/08_renormalization_tableX/check_global.py) | Helper from Task 13, imported by the global test; run on its own, it repeats the Task 13 self-test | ~14 s on its own |

```bash
cd code/08_renormalization_tableX
python3 check_146_IX.py
python3 check_160_167.py
python3 check_tableX_global.py
python3 check_global.py   # optional
```

The re-run outputs are identical to the original runs.

## Full derivation

The full step-by-step derivation is in the notes: [`../notes/notes.pdf`](../notes/notes.pdf), chapter "Renormalization and the final self-mass: Tables VI-X, Eqs. (I.146)-(I.167)". It includes every intermediate table, the index and coefficient bookkeeping, and the calculation records.

Partially open: the possible double counting flagged in Tasks 11–13 has exactly the counterterm form $`c(D)\,(a^D\square-a^{D-2}\nabla^2)\delta^D`$, so it is absorbed into $`\alpha_2,\alpha_3`$ and does not change the renormalized result. We have not shown, however, that the external-operator reduction of the appendix produces no other local terms at coincidence. That needs a full distributional analysis. The vanishing $`a^D\square^2`$ total and the finiteness of Table IX support it only partially.
