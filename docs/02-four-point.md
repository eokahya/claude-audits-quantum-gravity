# 2. The 4-point contribution: Eqs. (I.31)–(I.44)

<sub>[Overview (README)](../README.md) · Previous: [Feynman rules](01-feynman-rules.md) · Next: [3-point setup](03-three-point-setup.md)</sub>

*Paper I: E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007). This page covers Tasks 3, 4 and 5 of the re-derivation.*

## What this step computes

The simplest one-loop diagram for the scalar self-mass-squared has one 4-point
("seagull") vertex at $`x`$ with a graviton loop closing on itself. Its value is the
vertex operators $`U_I`$ of Table II, contracted with the graviton propagator at
coincident points, times $`\delta^D(x-x')`$. We re-derived three things. First, how the
derivatives of Table II act on the delta function (Task 3, Eqs. (I.31)–(I.33)). Second, the
coincidence limits of the scalar propagators and of the four contractions of the graviton
propagator in dimensional regularization (Task 4, Eqs. (I.34)–(I.40)). Third, the final
local expression after these are substituted (Task 5, Eqs. (I.41)–(I.44)). The result is a
divergent local term proportional to $`\pi\cot(\pi D/2)`$ (a pole at $`D=4`$), finite $`\ln(aa')`$ terms and a
$`9H^2a^D`$ term. Later chapters deal with these: the pole is removed by counterterms, and the
finite $`\ln(aa')`$ and $`9H^2a^D`$ terms cancel exactly against local terms of the 3-point
contribution ([Chapter 4](04-three-point-local-and-log.md)).

## Key steps

**1. One derivative per leg: the kernel of a quadratic form.** Averaging the 4-point
Lagrangian over the graviton loop gives an effective kinetic term
$`\int \tfrac{1}{2i}\,\mathcal K^{\mu\nu}\partial_\mu\Phi\,\partial_\nu\Phi`$. Its kernel follows from
a general lemma (Lemma A) for any quadratic form with one derivative on each field:

```math
\frac{\delta^2 F}{\delta\phi(x)\,\delta\phi(x')} = -\partial_\mu\Big[(K+K^T)^{\mu\nu}(x)\,\partial_\nu\delta^D(x-x')\Big],
\qquad F=\int K^{\mu\nu}\partial_\mu\phi\,\partial_\nu\phi .
```

One derivative ends up acting on everything as $`-\partial`$. The other acts on the delta function as
$`+\partial`$ **with respect to** $`x`$, which equals $`-\partial'`$. The paper's rule for (I.32) uses
$`+\partial'`$ instead. That one sign flips the whole 4-point contribution. Applied to the
symmetric $`\mathcal K`$, the lemma gives $`-iM^2_{\rm 4pt}=-\partial_\mu\big[\mathcal K^{\mu\nu}(x)\,\partial_\nu\delta^D(x-x')\big]`$.

**2. Our version of (I.33)** (sign fixed by the paper's own definition (I.168)):

```math
\begin{aligned}
-iM^2_{\rm 4pt} ={}& +\tfrac{i}{8}\kappa^2\partial^\mu\{a^{D-2}\,i[{}^\alpha_\alpha\Delta^\rho_\rho](x;x)\,\partial_\mu\delta^D\}
-\tfrac{i}{4}\kappa^2\partial^\mu\{a^{D-2}\,i[{}^{\alpha\beta}\Delta_{\alpha\beta}](x;x)\,\partial_\mu\delta^D\}\\
&-\tfrac{i}{2}\kappa^2\partial_\rho\{a^{D-2}\,i[{}^\alpha_\alpha\Delta^{\rho\sigma}](x;x)\,\partial_\sigma\delta^D\}
+i\kappa^2\partial^\alpha\{a^{D-2}\,i[{}_{\alpha\rho}\Delta^{\rho\sigma}](x;x)\,\partial_\sigma\delta^D\}.
\end{aligned}
```

The magnitudes, tensor contractions and relative signs agree with the paper; the overall
sign is opposite. The paper's (I.33) has $`+\tfrac i4`$ in the second term. Under the paper's
own rule a vertex $`U_2=c\,(\dots)`$ gives $`+c`$ there, so the paper used $`U_2=+\tfrac i4`$, the
value we derive. The $`-\tfrac i4`$ printed for $`U_2`$ in Table II is therefore a typo.

**3. Coincidence limits.** In dimensional regularization every power $`y^{p(D)}`$ with
$`p\not\equiv 0`$ vanishes at $`y=0`$. Only the constant terms survive:

```math
i\Delta_A(x;x)=k\Big\{-\pi\cot\big(\tfrac{\pi D}{2}\big)+2\ln a\Big\},\qquad
i\Delta_B(x;x)=-\frac{k}{D-2},\qquad i\Delta_C(x;x)=\frac{k}{(D-2)(D-3)},
```

with $`k=\frac{H^{D-2}}{(4\pi)^{D/2}}\frac{\Gamma(D-1)}{\Gamma(D/2)}`$. Here
$`-\pi\cot(\pi D/2)=-\tfrac{2}{D-4}+O(D-4)`$ is the pole.

**4. The mixed contraction (I.40)**, exactly in $`D`$ (all $`i\Delta_I`$ at $`(x;x)`$):

```math
i[{}_{\alpha\rho}\Delta^{\rho\sigma}](x;x)=\tfrac{D^2-3D-2}{D-3}\,\bar\delta^\sigma_\alpha\, i\Delta_A
+\big[(D-1)\delta^0_\alpha\delta^\sigma_0+\bar\delta^\sigma_\alpha\big]\,i\Delta_B
+\tfrac{2}{(D-2)(D-3)}\big[(D-3)^2\delta^0_\alpha\delta^\sigma_0+\bar\delta^\sigma_\alpha\big]\,i\Delta_C .
```

At $`D=4`$ its finite part is $`-2\,\delta^0_\alpha\delta^\sigma_0\,H^2/16\pi^2`$. The paper prints
$`+2`$. Its trace must equal (I.38) identically. The printed (I.40) gives $`+2H^2/16\pi^2`$
for that trace, while the printed (I.38) has $`-2H^2/16\pi^2`$.

**5. Assembling (I.41).** The $`1/(D-3)`$ factors of the graviton propagator cancel exactly:

```math
\mathcal K^{\mu\nu}_A = i\kappa^2 a^{D-2}\, i\Delta_A(x;x)\Big[\tfrac14 D(D-1)\,\eta^{\mu\nu}-D\,\bar\eta^{\mu\nu}\Big],
\qquad
\mathcal K^{\mu\nu}_{B,C} = i\kappa^2 a^{D-2}\,i\Delta_B\,\frac{D-4}{2(D-3)}\Big[(D-2)\bar\eta^{\mu\nu}+D\,u^\mu u^\nu\Big].
```

Here $`\mathcal K=\mathcal K_A+\mathcal K_{B,C}`$, $`\bar\eta^{\mu\nu}`$ is the spatial part of $`\eta^{\mu\nu}`$ and
$`u^\mu u^\nu=\delta^\mu_0\delta^\nu_0`$. With the correct (I.40), the finite $`B`$- and $`C`$-type parts are
exactly $`O(D-4)`$ and drop out.

**6. Extracting the logarithm, (I.43)**, also in general $`D`$:

```math
\partial^\mu\big(a^{D-2}\ln a\,\partial_\mu\delta^D(x-x')\big)=\tfrac12\ln(aa')\,a^D\square\,\delta^D(x-x')+\tfrac12(D-1)H^2a^D\,\delta^D(x-x'),
```

which follows from $`\square\ln a=-(D-1)H^2`$. It is valid only with the symmetric $`\ln(aa')`$,
not with $`\ln a^2`$.

## Result

> **Final 4-point contribution (sign of the definition (I.168)):**
>
> ```math
> -iM^2_{\rm 4pt}(x;x') = -i\kappa^2 k\Big\{\Big[-\tfrac14 D(D-1)\pi\cot\big(\tfrac{\pi}{2}D\big)+3\ln(aa')\Big]a^D\square
> +\Big[D\pi\cot\big(\tfrac{\pi}{2}D\big)-4\ln(aa')\Big]a^{D-2}\nabla^2+9H^2a^D+O(D-4)\Big\}\delta^D(x-x')
> ```
>
> In the paper's sign convention (multiply by $`-1`$), this is the printed (I.44) **without**
> the constants $`-2`$ and $`+2`$ inside the two brackets.

The algebra (I.41) → (I.44) in the paper is correct. Starting from the printed (I.33) and the
printed (I.40), we reproduce the printed (I.41) and (I.44) exactly. Both differences come
from earlier steps.

What this means for the final result (from later tasks): the 3-point contributions of the paper
turn out to be reversed as well (Tasks 6, 7, 9). So the relative sign of the 4-point and 3-point
parts is consistent, and the paper's whole bare self-mass is $`-1`$ times the one defined by
(I.168). The spurious $`\mp2`$ terms have exactly the form of the counterterms (I.131) and (I.145),
with $`\Delta\alpha_2=-\tfrac{1}{4\pi^2}`$ and $`\Delta\alpha_3=+\tfrac{1}{4\pi^2}`$ (paper convention). They only shift finite
local coefficients: $`\beta_2`$ by $`+4`$ and $`\beta_3`$ by $`-4`$ (Task 17).

## Comparison with the paper

| Equation | Our result | Paper I | Status |
|---|---|---|---|
| (I.31) | structure, normalization 1 | same | agrees |
| (I.32) | $`-\tfrac i8\kappa^2\partial^\mu\{\dots\partial'_\mu\delta^D\}`$ | $`+\tfrac i8`$ | overall sign (S1; derivative rule for $`\partial_2`$) |
| (I.33) | coefficients $`(+\tfrac i8,-\tfrac i4,-\tfrac i2,+i)`$ | $`(-\tfrac i8,+\tfrac i4,+\tfrac i2,-i)`$ | overall sign (S1) |
| Table II, $`U_2`$ | $`+\tfrac i4`$ | $`-\tfrac i4`$ | typo (C5; the paper's (I.33) uses $`+\tfrac i4`$) |
| (I.34)–(I.38) | exact in $`D`$; finite parts at $`D=4`$ | same | agree |
| (I.39) | finite part $`[2\delta^\rho_0\delta^\sigma_0+2\bar\eta^{\rho\sigma}]h`$ | $`\bar\eta_{\rho\sigma}`$ | index placement typo (C6; same matrix) |
| (I.40) | finite part $`-2\,\delta^0_\alpha\delta^\sigma_0\,h`$ | $`+2\,\delta^0_\alpha\delta^\sigma_0\,h`$ | sign error; finite, absorbable in $`\alpha_2,\alpha_3`$ (B1) |
| (I.41) | no $`\{-\partial^\mu(a^2\partial_\mu\delta^4)+a^2\nabla^2\delta^4\}`$ bracket; overall $`-`$ | bracket present; overall $`+`$ | inherited from (I.40) (B1) and (I.32) (S1) |
| (I.42), (I.43) | reproduced (plus general-$`D`$ (I.43)) | same | agree |
| (I.44) | no $`-2`$, $`+2`$; overall $`-`$ | $`-2`$, $`+2`$; overall $`+`$ | inherited from (I.40) (B1) and (I.32) (S1) |

Here $`h=H^2/16\pi^2`$. The IDs (S1, B1, C5, C6) refer to the [complete list of findings](10-complete-list-of-findings.md).

## Code

All scripts are in [`../code/02_four_point/`](../code/02_four_point/). They are supplementary checks
of the analytic derivation. Each one runs standalone; the expected output is in `expected_output/`.

| Script | What it checks | Runtime |
|---|---|---|
| [`check_4pt.py`](../code/02_four_point/check_4pt.py) | Lemma A with explicit test functions; lattice Hessian and positivity; operator-level (I.33) in $`D=3`$ with a generic coincident propagator; calibration with the paper's own counterterm vertex (I.131) | ~25 s |
| [`check_coincidence.py`](../code/02_four_point/check_coincidence.py) | (I.35)–(I.36) from Gauss's theorem; numerical coincidence limits at $`D=0.6`$; explicit-component contractions (I.37)–(I.40) for $`D=4,\dots,7`$ and trace identities; $`D=4`$ finite parts against the printed ones | ~2 s |
| [`check_41_44.py`](../code/02_four_point/check_41_44.py) | symbolic-$`D`$ assembly of (I.41); explicit $`D=4`$ kernel; (I.43) as a kernel identity; (I.41) → (I.44) on a test function | ~1 s |

```bash
cd code/02_four_point
python3 check_4pt.py
python3 check_coincidence.py
python3 check_41_44.py
```

Requirements: Python 3.10+, `sympy`, `mpmath`.

## Full derivation

The complete step-by-step derivation, with every intermediate result, index contraction,
comparison table and calculation record, is in the notes:
[`../notes/notes.pdf`](../notes/notes.pdf), chapter **"The 4-point contribution: Eqs. (I.31)-(I.44)"**.
