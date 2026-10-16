# 3. The 3-point contribution: setup and Eqs. (I.45)–(I.53)

<sub>[Overview (README)](../README.md) · Previous: [4-point contribution](02-four-point.md) · Next: [3-point local and log terms](04-three-point-local-and-log.md)</sub>

*Tasks 6, 7 and 8 of the re-derivation of Paper I (Kahya and Woodard, Phys. Rev. D 76, 124005 (2007)).*

## What this step computes

The one-loop scalar self-mass-squared gets a contribution from two cubic scalar–scalar–graviton vertices joined by a scalar
and a graviton propagator (Fig. 2 of Paper I). This step builds that contribution from the effective action, fixes its
normalization and the rules for which derivatives act on what, inserts the graviton propagator of Paper I, and sorts the
result into nine terms of three types ($`AA`$, $`AB`$, $`AC`$). It then prepares everything that follows: the identity
(I.52) for the mixed second derivative of the massless minimally coupled scalar propagator, which splits each term into a
*local*, a *logarithm* and a *normal* part. The printed (I.45), (I.46) and (I.50)–(I.53) are confirmed. The printed (I.47)–(I.49)
have the opposite overall sign, in all nine terms, to what follows from the paper's own (I.46). Separately, the first coefficient
of the propagator function (I.14) contains a typographical error.

## Key steps

**1. Effective action and (I.45).** With $`\phi=\Phi+\varphi`$, only the cross term $`\Phi\varphi h`$ of the cubic vertices gives a
1PI diagram. Collecting the factors $`\tfrac12\cdot 4\cdot(1/2i)^2\cdot i^2\cdot 2=1`$ gives

```math
-iM^2_{\rm 3pt}(x;x')=\sum_{I,J=1}^{2}V_I^{\alpha\beta}(x)\,V_J^{\rho\sigma}(x')\;i\big[{}_{\alpha\beta}\Delta_{\rho\sigma}\big](x;x')\;i\Delta_A(x;x'),
```

where a derivative on an *outer* leg becomes $`-\partial`$ acting on everything at its vertex (one partial integration per
vertex), and a derivative on an *inner* leg acts only on $`i\Delta_A`$. The two outer signs multiply to $`(-1)(-1)=+1`$.

**2. Expanding the vertices, (I.46).** With $`V_1^{\alpha\beta}=i\kappa a^{D-2}\partial_1^\alpha\partial_2^\beta`$ and
$`V_2^{\alpha\beta}=-\tfrac i2\kappa a^{D-2}\eta^{\alpha\beta}\partial_1\cdot\partial_2`$:

```math
\begin{aligned}
-iM^2_{\rm 3pt}={}&-\kappa^2\partial_\alpha\partial'_\rho\big\{(aa')^{D-2}\,i[{}^{\alpha\beta}\Delta^{\rho\sigma}]\,\partial_\beta\partial'_\sigma i\Delta_A\big\}
+\tfrac{\kappa^2}{2}\partial^\mu\partial'_\rho\big\{(aa')^{D-2}\,i[{}^{\alpha}_{\alpha}\Delta^{\rho\sigma}]\,\partial_\mu\partial'_\sigma i\Delta_A\big\}\\
&+\tfrac{\kappa^2}{2}\partial_\alpha\partial'^\nu\big\{(aa')^{D-2}\,i[{}^{\alpha\beta}\Delta^{\rho}{}_{\rho}]\,\partial_\beta\partial'_\nu i\Delta_A\big\}
-\tfrac{\kappa^2}{4}\partial^\mu\partial'^\nu\big\{(aa')^{D-2}\,i[{}^{\alpha}_{\alpha}\Delta^{\rho}{}_{\rho}]\,\partial_\mu\partial'_\nu i\Delta_A\big\}.
\end{aligned}
```

This agrees term by term with the paper. It also confirms that the $`-\tfrac12\kappa`$ on the left-hand side of (I.7) is only a typo:
with that sign the 2nd and 3rd terms would flip.

**3. Inserting the graviton propagator, (I.47)–(I.49).** Inserting $`i[{}_{\alpha\beta}\Delta_{\rho\sigma}]=\sum_{I}T^I_{\alpha\beta\rho\sigma}\,i\Delta_I`$
($`I=A,B,C`$) and contracting gives nine terms. For example, the $`T^C`$ parts of the four terms of (I.46) combine into one
product proportional to $`(D-2)^2\,P_0Q_0\,P'_0Q'_0`$ ($`P,P'`$ outer and $`Q,Q'`$ inner derivatives), which gives

```math
AC:\quad -2\kappa^2\,\frac{D-2}{D-3}\;\partial_0\partial'_0\Big[(aa')^{D-2}\,i\Delta_C\,\partial_0\partial'_0\,i\Delta_A\Big],
```

while the first $`AA`$ term is $`-\kappa^2\,\nabla\cdot\nabla'\big[(aa')^{D-2}\,i\Delta_A\,\nabla\cdot\nabla' i\Delta_A\big]`$. An extra $`AA`$ pair
$`\partial_i\partial'_j\{\ldots(\partial_i\partial'_j-\partial_j\partial'_i)\,i\Delta_A\}`$ vanishes because $`i\Delta_A`$ is spatially homogeneous.
In all nine terms the paper prints the opposite overall sign.

**4. The mixed derivative of the scalar propagator, (I.52).** With the $`i\delta`$ prescription ($`\delta>0`$), $`y`$ depends on time through
$`|\Delta\eta|`$ with $`\Delta\eta=\eta-\eta'`$. Its kink gives $`\partial_0\partial'_0(|\Delta\eta|-i\delta)^2=-2+4i\delta\,\delta(\Delta\eta)`$. Together with the
nascent delta function $`\delta\,(r^2+\delta^2)^{-D/2}\to\pi^{D/2}\Gamma(D/2)^{-1}\,\delta^{D-1}(\vec x-\vec x')`$ this gives

```math
\partial_\rho\partial'_\sigma\, i\Delta_A(x;x')=\delta^0_\rho\delta^0_\sigma\,\frac{i}{a^{D-2}}\,\delta^D(x-x')
+A''(y)\,\frac{\partial y}{\partial x^\rho}\frac{\partial y}{\partial x'^\sigma}+A'(y)\,\frac{\partial^2 y}{\partial x^\rho\,\partial x'^\sigma}.
```

The term $`k\ln(aa')`$ drops out because $`\partial_\rho\partial'_\sigma[\ln a+\ln a']=0`$. Two further, independent routes give the same local
coefficient: the propagator equation (I.12), which needs no knowledge of the normalization of $`A(y)`$, and the canonical
equal-time commutator, i.e. the contact term between the $`T`$ and $`T^*`$ products.

**5. The first coefficient of (I.14).** In the direct route the local coefficient is $`i\,c_0(\tfrac D2-1)\,a^{2-D}/\Gamma(\tfrac D2)`$,
where $`c_0`$ is the coefficient of the leading power $`(4/y)^{D/2-1}`$ inside the braces of (I.14). The two other routes fix the local coefficient to
$`i\,a^{2-D}`$, so $`c_0=\Gamma(\tfrac D2-1)`$. The same value follows from (I.26) and from the recursion that the printed
second coefficient $`c_1=\Gamma(\tfrac D2+1)/(\tfrac D2-2)`$ imposes. The paper prints $`c_0=\Gamma(\tfrac D2-1)/(\tfrac D2-1)`$.

## Result

> **(I.45)–(I.46)** are derived exactly as printed (normalization 1; outer legs $`\to-\partial`$ on everything, inner legs only on $`i\Delta_A`$).
>
> **(I.47)–(I.49)** (derived, definition sign):
> $`AA`$: $`-\kappa^2\,\nabla\cdot\nabla'[\ldots\nabla\cdot\nabla']`$, $`+\kappa^2\frac{D-1}{D-3}\,\partial_0\partial'_0[\ldots\partial_0\partial'_0]`$,
> $`-\kappa^2\,\partial_i\partial'_0[\ldots\partial_i\partial'_0]`$, $`-\kappa^2\,\partial_0\partial'_i[\ldots\partial_0\partial'_i]`$;
> $`AB`$: all four terms with $`+\kappa^2`$; $`AC`$: $`-2\kappa^2\frac{D-2}{D-3}\,\partial_0\partial'_0[\ldots\partial_0\partial'_0]`$.
> That is exactly $`-1`$ times the printed equations.
>
> **(I.52)** holds as printed. The local term $`\delta^0_\rho\delta^0_\sigma\, i\,a^{2-D}\,\delta^D(x-x')`$ is confirmed along three independent routes.
>
> **(I.50), (I.51), (I.53)** are correct, and the split into local, logarithm and normal parts is complete and free of overlap.
> The local parts come from exactly three terms: the 2nd term of (I.47), the 4th term of (I.48), and (I.49).

**Status of the overall sign.** The derivative rule of (I.45) is correct. The sign error enters at the step from (I.46)
to (I.47)–(I.49). In Task 9 ([next chapter](04-three-point-local-and-log.md)) we find that the paper carried the printed signs of (I.47)–(I.49) into the
rest of the 3-point computation. The 4-point contribution of Paper I is reversed as well ([Chapter 2](02-four-point.md)). The two contributions are
therefore consistent with each other, and the printed bare $`-iM^2`$ is $`-1`$ times the one defined by Paper I's own (I.168).
At the end of Task 6, before this was known, we concluded that the 4-point contribution entered with the opposite sign
relative to the 3-point one. That conclusion was later revised.

## Comparison with the paper

| Item | Our result | Paper I | Verdict |
|---|---|---|---|
| (I.45) normalization and derivative rules | 1; outer legs $`-\partial`$, inner legs on $`i\Delta_A`$ only | same | agrees |
| (I.46), four terms | coefficients $`-\kappa^2,\ +\tfrac{\kappa^2}{2},\ +\tfrac{\kappa^2}{2},\ -\tfrac{\kappa^2}{4}`$ | same | agrees |
| (I.7), left-hand side | $`+\tfrac12\kappa`$ (consistent with (I.46)) | $`-\tfrac12\kappa`$ | typo (C3) |
| (I.47)–(I.49), nine terms | structures, $`\tfrac{D-1}{D-3}`$, $`\tfrac{2(D-2)}{D-3}`$ as printed | same magnitudes | **overall sign reversed** (S1) |
| antisymmetric $`AA`$ pair | vanishes by homogeneity | not written | agrees |
| (I.50), (I.51), (I.53), split (i)–(iii) | correct, complete, no overlap | — | agrees |
| (I.52), local term | $`\delta^0_\rho\delta^0_\sigma\,i\,a^{2-D}\,\delta^D`$ (3 routes) | same | agrees |
| (I.52), regular part | chain rule | same | agrees |
| (I.14), first coefficient | $`\Gamma(\tfrac D2-1)`$ | $`\Gamma(\tfrac D2-1)/(\tfrac D2-1)`$ | typo (C4; printing only, (I.52) is right) |

The IDs (S1, C3, C4) refer to the [complete list of findings](10-complete-list-of-findings.md).

Had the printed (I.14) coefficient been used, the local term would carry an extra factor $`2/(D-2)`$. This factor equals 1 at
$`D=4`$, but multiplied by the pole of $`i\Delta_A(x;x)`$ it would leave a finite error $`2k/(D-2)`$.

## Code

All scripts are in [`../code/03_three_point_setup/`](../code/03_three_point_setup/) (Python 3, sympy, mpmath). Run them from
that folder. Each output is stored in `expected_output/` and is identical, character for character, to the output recorded in the
original task logs (for `check_52.py`, the final run of Task 8).

| Script | What it checks | Runtime |
|---|---|---|
| [`check_45_46.py`](../code/03_three_point_setup/check_45_46.py) | index algebra (I.45) to (I.46); normalization; derivative rule against a Gateaux derivative in a 1D toy model | about 3 min |
| [`check_47_49.py`](../code/03_three_point_setup/check_47_49.py) | (I.46) route equals (I.45) route; derived nine terms; result is minus the printed (I.47)–(I.49), for D = 4, 5, 6 | about 2 s |
| [`check_52.py`](../code/03_three_point_setup/check_52.py) | kink lemma, nascent delta, local coefficient, recursion for the coefficients of (I.14), box operator on f(y), derivatives of y | about 2 s |

```bash
cd code/03_three_point_setup
python3 check_45_46.py
python3 check_47_49.py
python3 check_52.py
```

## Full derivation

Every step (index tracking, contraction tables, the three derivations of the local coefficient, and the calculation records)
is in the notes: [`../notes/notes.pdf`](../notes/notes.pdf), chapter "The 3-point contribution: setup and Eqs. (I.45)-(I.53)".
