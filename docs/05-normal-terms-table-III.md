# 5. Normal contributions and Table III: Eqs. (I.68)–(I.93)

<sub>[Overview (README)](../README.md) · Previous: [3-point local and log terms](04-three-point-local-and-log.md) · Next: [Tables IV and V](06-tables-IV-V.md)</sub>

*Tasks 11 and 12. Paper I = Kahya & Woodard, Phys. Rev. D 76, 124005 (2007).*

## What this step computes

In the "normal" 3-point contributions of Paper I, two derivatives act on one propagator function $`A(y)`$ and the result is multiplied by a second propagator function $`f(y)`$.
Paper I does not evaluate these terms directly. It pulls the derivatives outside first, using repeated indefinite integrals $`I[f](y)=\int^y dy'\,f(y')`$, so that only "external operators" act on functions of the de Sitter length $`y(x;x')`$.
This step re-derives that machinery:

- the derivative-extraction identities (I.68)–(I.80);
- the time-derivative identities (I.81)–(I.88);
- their first application, which turns the nonlocal logarithm term of (I.67) into **Table III**;
- the tools (I.90)–(I.93) that convert double time derivatives into covariant d'Alembertians.

All identities were derived in general $`D`$ and then checked with SymPy.

## Key steps

**1. Extracting two derivatives** (I.70). For any $`f(y)`$ and $`A(y)`$:

```math
f\Big\{A''\,\frac{\partial y}{\partial x^\rho}\frac{\partial y}{\partial x'^\sigma}+A'\,\frac{\partial^2y}{\partial x^\rho\partial x'^\sigma}\Big\}
=\partial_\rho\partial'_\sigma\, I^2[fA'']+\frac{\partial^2y}{\partial x^\rho\partial x'^\sigma}\,I[f'A'] .
```

This holds as an identity between *naive* (pointwise) derivatives, as Paper I intends; see its footnote to (I.75).
Taken as distributions, both sides also pick up the same $`\delta^D(x-x')`$ term. That term comes from the $`y^{1-D/2}`$ component $`c_f A_{\rm lead}`$ of $`I^2[fA'']`$, where $`c_f`$ is the constant part of $`f`$.
It is exactly the local contribution (I.54)–(I.58), which is counted separately. So the inner derivatives in the normal terms must be taken naively, or that local term is counted twice.

**2. Mixed derivative of** $`y`$ (I.71), with the lower-index convention $`\Delta x_\sigma=\eta_{\sigma\nu}(x-x')^\nu`$ (the paper does not state it explicitly):

```math
\frac{\partial^2y}{\partial x^\rho\partial x'^\sigma}=H^2aa'\big\{y\,\delta^0_\rho\delta^0_\sigma-2a\,\delta^0_\rho H\Delta x_\sigma+2a'H\Delta x_\rho\delta^0_\sigma-2\eta_{\rho\sigma}\big\}.
```

**3. Trading a time derivative for a Laplacian** (I.88), valid in general $`D`$:

```math
(a\partial_0+a'\partial'_0)f(y)=H(a+a')^2\,y\,f'(y)+(D-1)Haa'f(y)-\frac{\nabla^2}{2H}\,I[f](y).
```

**4. Applying it to the nonlocal log term of (I.67).** With the re-derived (I.67), in the definition sign of (I.168), at $`D=4`$:

```math
\frac{\kappa^2H^2}{8\pi^2}\,2H(aa')^2(a\partial_0+a'\partial'_0)\nabla^2A
=\frac{\kappa^2H^2}{16\pi^2}\Big\{12(aa')^3H^2\nabla^2A+4H^2(aa')^2(a+a')^2\nabla^2(yA')-2(aa')^2\nabla^4I[A]\Big\}.
```

Then insert $`A=\tfrac{H^2}{16\pi^2}\{4/y-2\ln(y/4)-1\}`$ from (I.62). Constants are dropped because $`\nabla^2`$ kills them. In the $`\nabla^4`$ row, terms linear in $`y`$ are dropped as well, because $`\nabla^2y`$ is constant in space and so $`\nabla^4y=0`$.

**5. The covariant d'Alembertian on functions of** $`y`$ (I.93), used in the next chapter:

```math
\square f(y)=H^2\big[(4y-y^2)f''(y)+D(2-y)f'(y)\big]=\square' f(y).
```

This is a pointwise identity. In (I.92) the $`\square`$ is a distributional operator. If it is evaluated pointwise with (I.93), a legitimate local term $`-ik\hat\pi\,a^{D-2}\nabla^2\delta^D`$ is lost. This point was followed up in later tasks (see below).

## Result

> **Table III** (coefficients of $`\kappa^2H^4/(4\pi)^4`$, with $`x=y/4`$), as derived:
>
> | External operator | Definition sign | Paper convention |
> |---|---|---|
> | $`(aa')^3H^2\nabla^2`$ | $`+\frac{12}{x}-24\ln x`$ | $`-\frac{12}{x}+24\ln x`$ |
> | $`(aa')^2(a+a')^2H^2\nabla^2`$ | $`-\frac{4}{x}`$ | $`+\frac{4}{x}`$ |
> | $`(aa')^2\nabla^4`$ | $`-8\ln x+16x\ln x`$ | $`+8\ln x-16x\ln x`$ |
>
> In the paper convention this is exactly the printed Table III. The overall sign difference is the one inherited from (I.47): the printed $`-iM^2`$ is $`-1`$ times the one defined by (I.168).
> The printed (I.89) has, besides this inherited sign, a typo: its second term lacks a factor $`H^2`$. Table III itself has the $`H^2`$.
> All other identities, (I.68)–(I.88) and (I.90)–(I.93), agree with the paper.

## Comparison with the paper

| Item | Our result | Paper | Verdict |
|---|---|---|---|
| (I.68)–(I.80) | derived in general $`D`$ | same | agree |
| (I.69): "never a delta function" | true for naive derivatives; distributional derivatives add $`c_f\,\delta^D`$ terms that equal the separately counted local part | same statement | agree (naive reading) |
| Divergence count after (I.80) | not locally integrable at $`D=4`$: only $`y^{2-D}`$ in $`I^2[AA'']`$, and $`y^{1-D}`$, $`y^{2-D}`$ in $`I[A'^2]`$ | same | agree |
| (I.81)–(I.88), (I.90)–(I.93) | derived in general $`D`$ | same | agree |
| (I.89), second term (paper convention) | $`-4H^2(aa')^2(a+a')^2\nabla^2(yA')`$ | $`-4(aa')^2(a+a')^2\nabla^2(yA')`$ | **typo** (C8): $`H^2`$ missing (dimensions; numerical check); Table III is correct |
| (I.89) and Table III, overall sign | definition sign | paper convention | sign inherited from (I.47) (S1), not a new error |

The IDs (S1, S2, C8) refer to the [complete list of findings](10-complete-list-of-findings.md).

**Later findings (Tasks 13, 16 and 17), recorded here for orientation.**

- *Local terms (steps 1 and 5).* Task 13 takes the inner derivatives of the normal contributions naively, following the paper's appendix text and its footnote to (I.75).
  Task 16 finds that the possible double counting has exactly the counterterm form $`c\,(a^D\square-a^{D-2}\nabla^2)\delta^D`$, with $`c`$ a constant and no $`\ln a`$. It is therefore absorbed into $`\alpha_2`$ and $`\alpha_3`$ and does not change the renormalized result.
  Whether the external-operator reduction produces no *other* local differences is recorded there as only partly settled.
- *Table III in Table X (Task 17).* Paper I builds Table X from Tables III, V and VII and Eq. (I.159).
  It adds the first two rows of Table III, but leaves the third row, $`8\ln x-16x\ln x`$ (paper convention), out of the $`(aa')^2\nabla^4`$ row of Table X (S2).
  The missing term is non-analytic, so it changes the physical nonlocal self-mass (it enters the retarded kernel, although at late times it adds only falling terms); see [Renormalization and Table X](08-renormalized-self-mass-table-X.md).

## Code

Folder: [`../code/05_normal_terms_tableIII/`](../code/05_normal_terms_tableIII/). Python 3.10 or later, with sympy and mpmath.

| Script | What it checks | Runtime |
|---|---|---|
| [`check_68_80.py`](../code/05_normal_terms_tableIII/check_68_80.py) | (I.70)–(I.80) with $`f=y^\alpha`$, $`A=y^\beta`$ in explicit coordinates ((I.70) at $`D=4`$, (I.71)–(I.72) exactly, (I.74)–(I.80) at $`D=4,5`$); the $`c_f A_{\rm lead}`$ delta content; the power counting after (I.80) | about 18 s |
| [`check_81_93.py`](../code/05_normal_terms_tableIII/check_81_93.py) | (I.81)–(I.93); (I.89) and Table III at three points with the true $`A(y)`$ (printed table $`=-`$derived; (I.89) holds only with $`H^2`$ restored) | about 20 s |

```bash
cd code/05_normal_terms_tableIII
python3 check_68_80.py
python3 check_81_93.py
```

The outputs we got are in `expected_output/`; the [folder README](../code/05_normal_terms_tableIII/README.md) explains each line.

## Full derivation

Every step, with index tracking, the distributional analysis and the calculation records, is in [the notes PDF](../notes/notes.pdf), chapter "Normal contributions and Table III: Eqs. (I.68)-(I.93)".
