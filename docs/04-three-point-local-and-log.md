# 4. 3-point local and logarithm terms: Eqs. (I.54)–(I.67)

<sub>[Overview (README)](../README.md) · Previous: [3-point setup](03-three-point-setup.md) · Next: [Normal terms and Table III](05-normal-terms-table-III.md)</sub>

*Paper I*: E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007). This page covers Tasks 9 and 10 of the re-derivation.

## What this step computes

The 3-point diagram has two graviton–scalar vertices joined by a graviton propagator and a scalar propagator $`i\Delta_A`$. Paper I splits its contribution to the scalar self-mass into three parts:

- **local** terms, from the contact term in $`\partial_0\partial_0' i\Delta_A`$;
- **logarithm** terms, from the de Sitter-breaking piece $`k\ln(aa')`$ of the outer propagator $`i\Delta_A = A(y) + k\ln(aa')`$;
- **normal** terms, which contain everything else.

This step re-derives the first two parts at general dimension $`D`$. Task 9 covers the local part, Eqs. (I.54)–(I.59). Task 10 covers the logarithm part, Eqs. (I.60)–(I.67).

In short: Eqs. (I.55)–(I.56) and (I.60)–(I.62) agree with the paper as printed. Every contribution to the self-mass, (I.54), (I.57)–(I.59) and (I.63)–(I.67), has the same structure as printed but the opposite overall sign. The source is a sign slip in the step (I.46) → (I.47)–(I.49), which the paper carried into all later 3-point algebra. Apart from the sign, the magnitudes agree except for one arithmetic slip in (I.59), which a counterterm absorbs. There is also one typo in the text. Finally, we confirm that the $`\ln(aa')`$ and $`H^2a^D`$ local terms of the 4-point result (I.44) cancel against those of (I.67).

Notation: $`k = \frac{H^{D-2}}{(4\pi)^{D/2}}\frac{\Gamma(D-1)}{\Gamma(D/2)} \to \frac{H^2}{8\pi^2}`$ and $`\hat\pi = -\pi\cot(\pi D/2) = -\frac{2}{D-4}+O(D-4)`$. All derived formulas below use the **definition sign**, fixed by Paper I's own Eq. (I.168). Paper I prints the opposite overall sign.

## Key steps

**1. Local terms (Task 9).** Only three of the nine terms of (I.47)–(I.49) have the inner pair $`\partial_0\partial_0'`$. In these, the contact term $`\partial_0\partial_0' i\Delta_A \supset i\,a^{2-D}\delta^D(x-x')`$ sets the outer propagator equal to its coincidence limit: $`A\to k\hat\pi`$, $`i\Delta_B\to -\tfrac{k}{D-2}`$ and $`i\Delta_C\to\tfrac{k}{(D-2)(D-3)}`$. Setting $`A(y)\,\delta^D = k\hat\pi\,\delta^D`$ is an assumption of dimensional regularization. Following the paper's split, only $`A(y)`$ is used here; the $`k\ln(aa')`$ piece of $`i\Delta_A`$ goes to the logarithm terms. The outer derivatives then follow from

```math
\partial_0\partial_0'\big[a^{D-2}\delta^D\big] = \big(a^D\square - a^{D-2}\nabla^2\big)\delta^D,
\qquad
\nabla\!\cdot\!\nabla'\big[a^{D-2}\delta^D\big] = -a^{D-2}\nabla^2\delta^D .
```

Summing the three terms at exact $`D`$ gives

```math
-iM^2_{\rm 3pt,loc} = i\kappa^2 k\Big\{\Big[\tfrac{D-1}{D-3}\hat\pi-\tfrac{2}{(D-3)^2}\Big]a^D\square
+\Big[-\tfrac{D-1}{D-3}\hat\pi+\tfrac{1}{D-2}+\tfrac{2}{(D-3)^2}\Big]a^{D-2}\nabla^2\Big\}\delta^D(x-x').
```

**2. Logarithm terms (Task 10).** Here the propagator equation (I.12) does most of the work:

```math
\partial_0\big(a^{D-2}\partial_0\, i\Delta_A\big) = -i\delta^D(x-x') + a^{D-2}\nabla^2 A(y).
```

The time-dependent factor $`\ln(aa')`$ in front of a delta function is moved with a distributional lemma:

```math
-2i\,\partial_0\big[a^{D-2}\ln a\;\partial_0\delta^D\big] = i\ln(aa')\big(a^D\square - a^{D-2}\nabla^2\big)\delta^D + i(D-1)H^2a^D\delta^D .
```

The four terms of (I.47) then add up, at general $`D`$, to

```math
\begin{aligned}
-iM^2_{\rm 3pt,log} = \kappa^2 k\Big\{&\ln(aa')\Big[\tfrac{D-1}{D-3}a^D\square-\tfrac{4}{D-3}a^{D-2}\nabla^2\Big]i\delta^D+\tfrac{(D-1)^2}{D-3}H^2a^D\, i\delta^D\\
&-\tfrac{2(D-4)}{D-3}(aa')^{D-2}\ln(aa')\nabla^4A+\tfrac{2}{D-3}H(aa')^{D-2}\big(a\partial_0+a'\partial_0'\big)\nabla^2A\Big\}.
\end{aligned}
```

The $`\nabla^4 A`$ term has a coefficient of order $`D-4`$ and drops out at $`D=4`$. No pole can compensate: the term involves only one propagator, and $`A(y)\sim\Delta x^{2-D}`$ is locally integrable, so its derivatives are well-defined distributions with no pole at $`D=4`$. For the same reason, setting $`D=4`$ in the other coefficients changes them only at order $`D-4`$.

**3. Cancellation against the 4-point terms.** Add the $`\ln(aa')`$ and $`H^2a^D`$ local terms of the correct (I.44) to the local terms above. The sum is

```math
-i\kappa^2 k\,\frac{(D-4)(D+1)}{4(D-3)}\Big\{(D-1)\ln(aa')\,a^D\square-4\ln(aa')\,a^{D-2}\nabla^2+(D-1)^2H^2a^D\Big\}\delta^D \;\xrightarrow{D\to4}\; 0 .
```

This remainder is not multiplied by a pole, so it vanishes at $`D=4`$: the terms $`3\ln(aa')a^D\square`$, $`-4\ln(aa')a^{D-2}\nabla^2`$ and $`9H^2a^D`$, which are not of counterterm form, cancel exactly. The cancellation also confirms that the relative sign of the 4-point and 3-point contributions is the one fixed by the definitions; with the other relative sign these terms would add. The only local 4-point remainder left over is the $`\mp2`$ term from the slip in (I.40) (Chapter 2), which has counterterm form.

## Result

> **Local part (definition sign):**
> ```math
> -iM^2_{\rm 3pt,loc}(x;x') = i\kappa^2k\Big\{\Big[-\tfrac{D-1}{D-3}\pi\cot\big(\tfrac{\pi}{2}D\big)-2\Big]a^D\square+\Big[\tfrac{D-1}{D-3}\pi\cot\big(\tfrac{\pi}{2}D\big)+\tfrac52\Big]a^{D-2}\nabla^2+O(D-4)\Big\}\delta^D(x-x')
> ```
> **Logarithm part (definition sign):**
> ```math
> -iM^2_{\rm 3pt,log}(x;x') = \frac{\kappa^2H^2}{8\pi^2}\Big\{\ln(aa')\big[3a^4\square-4a^2\nabla^2\big]i\delta^4(x-x')+9H^2a^4\, i\delta^4(x-x')+2H(aa')^2\big(a\partial_0+a'\partial_0'\big)\nabla^2A(y)+O(D-4)\Big\}
> ```
> where $`A(y)=\frac{H^2}{16\pi^2}\big\{\frac4y-2\ln\frac y4-1\big\}`$ at $`D=4`$, which agrees with (I.62). We obtain it with the corrected first coefficient of (I.14); the printed coefficient changes it only at order $`D-4`$.

## Comparison with the paper

Derived results use the definition sign. Paper I's expressions are as printed and use the paper's convention, which is $`-1`$ times the definition sign.

| Equation | Derived | Paper I (printed) | Verdict |
|---|---|---|---|
| (I.54) | $`\ldots\{+a^D\square-a^{D-2}\nabla^2\}`$ | $`\ldots\{-a^D\square+a^{D-2}\nabla^2\}`$ | overall sign (inherited, S1) |
| (I.55)–(I.56) | same | same | agree |
| (I.57) | $`+\frac{i\kappa^2H^2}{16\pi^2}a^2\nabla^2\delta^4`$ | $`-\frac{i\kappa^2H^2}{16\pi^2}a^2\nabla^2\delta^4`$ | overall sign (inherited) |
| (I.58) | $`-\frac{i\kappa^2H^2}{4\pi^2}\{a^4\square-a^2\nabla^2\}\delta^4`$ | $`+\frac{i\kappa^2H^2}{4\pi^2}\{\ldots\}`$ | overall sign (inherited) |
| (I.59), $`a^D\square`$ | $`-\frac{D-1}{D-3}\pi\cot-2`$ | $`+\frac{D-1}{D-3}\pi\cot+2`$ | overall sign (inherited) |
| (I.59), $`a^{D-2}\nabla^2`$ | $`+\frac{D-1}{D-3}\pi\cot+\frac52`$ | $`-\frac{D-1}{D-3}\pi\cot-\frac72`$ | sign, and $`\frac72`$ should be $`\frac52`$: an **arithmetic slip**, absorbable (B2; $`\Delta\alpha_3=-\frac{1}{8\pi^2}`$, shifts $`\beta_3`$ by $`+2`$) |
| text before (I.54) | the 2nd term of (I.47) | "the 4th term of (47)" | typo (C7) |
| (I.60)–(I.62) | same | same | agree |
| (I.63)–(I.66) | four terms, derived at general $`D`$ | the same four terms, each with opposite sign | overall sign (inherited) |
| (I.67) | $`\{3,\,-4,\,9,\,2\}`$ structure | $`\{-3,\,4,\,-9,\,-2\}`$ | overall sign (inherited); magnitudes agree |
| local cancellation (I.44) + (I.67) | holds at general $`D`$, remainder of order $`D-4`$ | stated in the text | agree (holds in both conventions) |

The IDs (S1, B2, C7) refer to the [complete list of findings](10-complete-list-of-findings.md).

Taken together, the 4-point and 3-point contributions in Paper I are consistent with each other: both are reversed relative to the definition (I.168). Paper I's bare $`M^2`$ is therefore $`-1`$ times the defined one. The sign slip in (I.47)–(I.49) is not a new finding of this step; it was located in Task 7. What Task 9 adds is that the paper carried the printed signs of (I.47)–(I.49) into the rest of the 3-point calculation, which answers the question Task 7 left open. This also revises the Task 6 conclusion that the 4-point contribution entered with the opposite sign relative to the 3-point one. Task 10 found no new error.

Still open after this step: the normal contributions, (I.68) and later, and the treatment of the term $`2H(aa')^2(a\partial_0+a'\partial_0')\nabla^2A`$ of (I.67).

## Code

Scripts are in [`code/04_three_point_local_log/`](../code/04_three_point_local_log/). Both reproduce the original outputs of Tasks 9 and 10 exactly.

| Script | What it checks | Runtime |
|---|---|---|
| [`check_54_59.py`](../code/04_three_point_local_log/check_54_59.py) | Takes the local polynomial straight from (I.45) and Table I for D = 4, 5, 6. Checks the outer-derivative rules and (I.55)–(I.56). Assembles (I.59) and audits the paper's own arithmetic. | ~1 s |
| [`check_60_67.py`](../code/04_three_point_local_log/check_60_67.py) | Checks (I.62) from (I.14) and the delta-function lemma. Runs a 30-digit pointwise check of the non-local parts of (I.63)–(I.66). Checks the general-D sum and the cancellation against (I.44), and that the paper's terms are minus the derived ones. | ~8 s |

```bash
cd code/04_three_point_local_log
python3 check_54_59.py
python3 check_60_67.py
```

The outputs we got are in `expected_output/`; the [folder README](../code/04_three_point_local_log/README.md) explains how to read them.

## Full derivation

See the chapter "3-point local and logarithm terms: Eqs. (I.54)–(I.67)" in the [derivation notes (PDF)](../notes/notes.pdf). It contains every step, the index bookkeeping, the general-D intermediate results, the comparison tables and the calculation records.
