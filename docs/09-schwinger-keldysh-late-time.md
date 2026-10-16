# 9. Schwinger–Keldysh mode equation and the late-time mode function (Paper II)

<sub>[Overview (README)](../README.md) · Previous: [Renormalized self-mass and Table X](08-renormalized-self-mass-table-X.md) · Next: [Complete list of findings](10-complete-list-of-findings.md)</sub>

*Task 18 of the re-derivation. Paper II: E. O. Kahya and R. P. Woodard, "Scalar Field Equations from Quantum Gravity
during Inflation", arXiv:0710.5282.*

## What this step computes

Paper II puts the one-loop self-mass-squared of Paper I into the Schwinger–Keldysh (in-in) effective field
equation for the plane-wave mode function of the massless, minimally coupled scalar. It concludes that the mode
function receives no significant one-loop correction at late times. We redo this analysis with **our corrected
self-mass-squared**. Three things changed in the self-mass table with respect to the printed one: the overall sign
(the printed self-mass is $`-1`$ times the one defined by Paper I's own (I.168); S1), the ζ row of Table X (S2), and the
constant 72γ₂ (S3). In addition, the finite parts of the counterterm coefficients $`\beta_{2,3}`$ shifted (B3–B6). The
IDs refer to the [complete list of findings](10-complete-list-of-findings.md). The analysis has four parts:

- build the self-mass-squared for all four SK polarities;
- re-derive Tables 2–5 of Paper II without using them as input;
- go **beyond Paper II's simplification (II.44)** by keeping the full momentum dependence;
- ask whether our corrections change the late-time mode function.

The answer to the last question is no. The paper's conclusion holds with the correct input. However, the
paper's justification of (II.44) does not hold term by term (S4).

Notation: $`x=y/4`$, $`\Delta\eta=\eta-\eta'`$, $`r=\Vert\vec x'-\vec x\Vert`$, $`\hat k=k/H`$, and $`H=1`$ in intermediate
steps. The *Table-X convention* is the sign of the printed tables. The *definition sign* is that of (I.168) and equals
−1 times the Table-X convention. Physical statements below use the definition sign.

## Key steps

**1. Schwinger–Keldysh self-masses.** Call $`\mathcal T(y)`$ the sum of the corrected Table 1 of Paper II, in the
Table-X convention. The $`-`$ vertex is the complex conjugate of the $`+`$ vertex and gives a factor −1. Mixed-polarity
diagrams have no local terms. In the definition sign this gives

```math
-iM^2_{++}=i\kappa^2a^4\Big(\Delta\alpha_1\square^2+\Delta\alpha_2H^2\square+\Delta\alpha_3H^2\frac{\nabla^2}{a^2}\Big)\delta^4(x-x')-\mathcal T(y_{++}),\qquad
-iM^2_{+-}=+\,\mathcal T(y_{+-}) .
```

For $`\eta'<\eta`$ the nonlocal part of $`M^2_{++}+M^2_{+-}`$ is $`2\,\mathrm{Im}\,\mathcal T(y_{++})`$, which is real. It
vanishes for $`\eta'>\eta`$ and outside the light cone, so the kernel is real and causal.

**2. The mode equation** (II.29). Write $`\Phi=u+\kappa^2\Phi_1`$ (Paper II, (II.23)), where $`u`$ is the tree-order mode
function, and let the subscript 1 denote the one-loop self-mass-squared with its factor $`\kappa^2`$ stripped off. Then

```math
\Big[\partial_t^2+3H\partial_t+\frac{k^2}{a^2}\Big]\Phi_1=\mathcal R_1(\eta,k)\equiv-\frac{1}{a^4}\int_{\eta_i}^{0}\!d\eta'\,u(\eta',k)\int\! d^3x'\,\big\{M^2_{++}+M^2_{+-}\big\}_1\,e^{i\vec k\cdot(\vec x'-\vec x)} .
```

The source $`\mathcal R_1`$ and $`\Phi_1`$ carry no factor $`\kappa^2`$: $`\kappa^2`$ appears once, in $`\Phi=u+\kappa^2\Phi_1`$.
Table 5 entries and $`\mathcal R_1`$ are quoted in units of $`H^4u(0,k)/16\pi^2`$.

**3. SK differences** (II.39)–(II.41). Write $`\Delta f=f(y_{++}/4)-f(y_{+-}/4)`$. Inside the light cone:

```math
\Delta\ln x=2\pi i\,\theta(\Delta\eta-r),\qquad \Delta\ln^2x=4\pi i\,\theta(\Delta\eta-r)\ln|x|,\qquad
\Delta\frac1x=-8\pi i\,\delta(y),\qquad \Delta(\text{const},\,x,\,x^2)=0 .
```

The last relation means that the corrected 72γ₂ constant drops out of the SK kernel identically.

**4. Exact integrals.** Substitute $`a'=as`$. Then $`X_0=\tfrac14aa'\Delta\eta^2=(1-s)^2/(4s)`$ no longer depends on $`a`$.
Every integral reduces to $`\int s^p\{1,\ln s,\ln(1-s)\}ds`$, for example

```math
S_K[\ln x]\equiv-\frac{i}{16\pi^2}\int d^4x'\,a'^K\,\Delta\ln x=\frac16\int_1^a da'\,a'^{K-2}\Big(\frac1{a'}-\frac1a\Big)^3 .
```

The only source of $`\pi^2`$ in Paper II's tables is $`\int_0^1 \ln(1-s)\,ds/s=-\pi^2/6`$.

**5. Beyond (II.44).** We keep the exact $`u(\eta',k)`$, the phase $`e^{i\vec k\cdot\Delta\vec x}`$, the $`k^2/a^2`$ part
of $`\square`$, and all $`\nabla^2`$ rows, and we expand exactly in $`k`$. Some single rows then **grow** with $`a`$
(Table-X convention):

```math
\text{row }(aa')^2(a^2{+}a'^2)H^2\nabla^2:\ +\tfrac{4}{15}\hat k^2a,\qquad
\text{row }(aa')^2(a^2{+}a'^2)\nabla^2\square:\ -\tfrac{4}{15}\hat k^2a .
```

These two rows cancel each other. In the full sum, every non-falling term cancels separately at each order
$`k^0,\dots,k^{10}`$. At order $`k^4`$ this cancellation **requires** the $`24x^2\ln x`$ term of the ζ row, which
contributes $`-\hat k^4/70`$ in the definition sign ($`+\hat k^4/70`$ in the Table-X convention). This term is present
in both the printed and the corrected row, so the non-falling part also cancels with the printed table; the ζ correction
itself changes only falling terms (step 6).

**6. Late-time right-hand side** (definition sign, corrected table, through $`k^{10}`$). There is no $`1/a`$ term:

```math
\mathcal R_1=\frac{H^4u(0,k)}{16\pi^2}\Big\{\frac{1}{a^2}\Big[12-\Big(8\ln a+\tfrac{103}{3}\Big)\hat k^2-8i\hat k^3+\tfrac83\hat k^4+\dots\Big]-\frac{12i\hat k^3}{a^3}+O\Big(\frac{\ln a}{a^4}\Big)\Big\}.
```

With the printed ζ row the right-hand side would also contain $`\tfrac{4}{45}\hat k^4/a+\tfrac{4}{175}\hat k^6/a+\dots`$
(inside the same braces). The ζ correction removes these $`1/a`$ tails, so with the corrected table the right-hand side
falls like $`\ln a/a^2`$. The $`a^{-2}`$ coefficients depend on the initial time and on $`\hat k`$, i.e. on the initial
state. (The original Turkish working report wrote an extra explicit $`\kappa^2`$ in the source and in the particular
solution, counting $`\kappa^2`$ twice. The notes count it once, as here, and record the slip in a footnote.)

**7. Solving for Φ₁.** With $`\theta=a\partial_a`$ (so $`\partial_t=H\theta`$) the left-hand side is
$`H^2(\theta^2+3\theta)`$ plus the $`k^2/a^2`$ term, which adds only $`O(a^{-4})`$ to the $`a^{-2}`$ solution. For
$`(\theta^2+3\theta)F=S`$, a source $`S=a^{-2}(\alpha\ln a+\beta)`$ gives
$`F=a^{-2}\big(-\tfrac{\alpha}{2}\ln a+\tfrac{\alpha}{4}-\tfrac{\beta}{2}\big)`$. A source $`\propto a^{-3}`$ is resonant
and gives $`a^{-3}\ln a`$. Only falling terms are generated, plus a constant multiple $`C\,u(\eta,k)`$ of the homogeneous
solution; the initial condition removes the other homogeneous solution, $`u^*`$. The result is

```math
\Phi_1=C\,u(\eta,k)+\frac{H^2u(0,k)}{16\pi^2}\Big\{\frac1{a^2}\Big[4\hat k^2\ln a-6+\tfrac{91}{6}\hat k^2+4i\hat k^3+\dots\Big]+\frac{4i\hat k^3\ln a}{a^3}+\dots\Big\},
\qquad C\equiv\frac{H^2}{16\pi^2}\,\tilde C .
```

## Result

> **Late-time mode function (definition sign, corrected self-mass, no use of (II.44)):**
>
> ```math
> \Phi(\eta,k)=u(\eta,k)+\kappa^2\Phi_1(\eta,k)\ \xrightarrow{\ a\to\infty\ }\ u(0,k)\Big\{1+\frac{\kappa^2H^2}{16\pi^2}\Big[\tilde C+O\Big(\frac{\ln a}{a^2}\Big)\Big]\Big\}
> ```
>
> $`\tilde C`$ is a constant (a field-strength renormalization, from the homogeneous part $`C\,u`$ of $`\Phi_1`$, with
> $`\tilde C=16\pi^2C/H^2`$). Nothing grows like $`\ln a`$, $`\ln^2 a`$ or $`a`$, so **there is no significant one-loop
> late-time correction**, as Paper II states.
> None of our corrections to Paper I changes this. The overall sign (S1) flips every Table 5 entry and the particular
> solution $`\Phi_1`$; since the non-falling total is zero, only decaying terms change sign. The ζ correction (S2)
> enters only falling terms (bounded by a constant times $`\hat k^4/a`$ for every $`k`$) and removes the $`1/a`$ tails. The
> 72γ₂ constant (S3) drops out identically. The finite shifts of $`\beta_{2,3}`$ (B3–B6) are absorbed into
> $`\Delta\alpha_{2,3}`$.

## Comparison with the paper

| Item | Paper II | This work | Status |
|---|---|---|---|
| SK rules and effective equation (II.15)–(II.21) | stated | derived; written in the definition sign | ✓ |
| Finite parts of $`\alpha_2`$, $`\alpha_3`$, (II.8)–(II.9) | $`\{\tfrac{61}{3}\}`$, $`\{\tfrac{4}{D-4}-\tfrac{58}{3}-2\gamma\}`$ | $`\{\tfrac{49}{3}\}`$, $`\{\tfrac{4}{D-4}-\tfrac{52}{3}-2\gamma\}`$ (Table-X convention); pole unchanged | absorbable (B5, B6) |
| Local terms of (II.10) | $`\Delta\alpha_2\square+\Delta\alpha_3\nabla^2/a^2`$ | $`\Delta\alpha_2H^2\square+\Delta\alpha_3H^2\nabla^2/a^2`$ | typo (C20) |
| Local term (II.34) | only $`\Delta\alpha_3`$ contributes, $`\propto k^2u/a^2`$ (printed without $`H^2`$) | same, with $`H^2`$: $`\mathcal R_{1,\rm local}=-\Delta\alpha_3H^2k^2u/a^2\to0`$ | ✓ (typo C21) |
| Identities (II.35)–(II.41) | stated | re-derived; $`1/x`$ also directly from $`\delta(y)`$ | ✓ |
| Tables 2–4 (12 integrals) | printed | identical (exact closed form; quadrature agrees to $`\sim10^{-19}`$) | ✓ |
| Table 5 (19 entries), total 0 | printed in the Table-X convention | identical in the Table-X convention; each entry −1× in the definition sign; total 0 | ✓ (sign convention, S1) |
| Simplification (II.44): ∇² rows negligible | argued from extra factors of 1/a′ or 1/a | false row by row (two rows grow like $`\hat k^2 a`$); the total cancels order by order through $`k^{10}`$ | justification ✗ (S4), result ✓ |
| Fall-off of the right-hand side | at least $`\ln^n a/a`$ | $`\ln a/a^2`$ with the corrected table; $`\hat k^4/a`$ with the printed one | ✓ |
| Late-time mode function | no significant correction | no significant correction | ✓ |

**Limitations.**

- The cancellation of the non-falling terms is exact through $`k^{10}`$. Beyond that, only a convergent-series argument
  is available (entire integrand, compact region); there is no all-orders proof.
- The $`a^{-2}`$ coefficients of $`\mathcal R_1`$ and $`\Phi_1`$ depend on the initial state.
- The explanation of the $`k=0`$ cancellation ($`M^2=\partial\partial'K`$) is a heuristic, labelled as an assumption in
  the notes. The operator ordering (scale factors on the left) is assumed; it was verified in the Task 17 pointwise test.
- The direct light-cone evaluation of $`S_K[1/x]`$ was done by hand only.

The complete list of findings and limitations is on the [next page](10-complete-list-of-findings.md).

## Code

Folder: [`../code/09_schwinger_keldysh/`](../code/09_schwinger_keldysh/). See the README there. All outputs agree
line by line with the original run.

| Script | What it checks |
|---|---|
| [`tables234.py`](../code/09_schwinger_keldysh/tables234.py) | Tables 2–4 of Paper II (12 integrals), including the next order |
| [`numcheck.py`](../code/09_schwinger_keldysh/numcheck.py) | exact closed-form integrals versus direct quadrature (largest difference 1.0e-19) |
| [`table5.py`](../code/09_schwinger_keldysh/table5.py) | Table 5 of Paper II (19 entries, total 0) |
| [`kexp.py`](../code/09_schwinger_keldysh/kexp.py) | full right-hand side beyond (II.44), expanded exactly in $`k`$ through order $`k^N`$ |
| [`numcheck_k.py`](../code/09_schwinger_keldysh/numcheck_k.py) | $`k`$-series versus quadrature with the exact $`k`$ dependence |
| [`analysis2.py`](../code/09_schwinger_keldysh/analysis2.py), [`analysis3.py`](../code/09_schwinger_keldysh/analysis3.py) | printed versus corrected table; contribution of the ζ correction |
| [`extract.py`](../code/09_schwinger_keldysh/extract.py) | the coefficients quoted above (definition sign) |

Helper modules: [`engine.py`](../code/09_schwinger_keldysh/engine.py) (exact integrals) and
[`rows.py`](../code/09_schwinger_keldysh/rows.py) (corrected Table 1).

```sh
cd code/09_schwinger_keldysh
python3 tables234.py && python3 table5.py && python3 kexp.py 10
# compare the whole folder with the stored outputs (about 2 minutes):
python3 ../run_all.py --folder 09
# note: sh run_all.sh in this folder regenerates (overwrites) expected_output/
```

## Full derivation

For every step, the index and sign bookkeeping, all tables, and the calculation record, see
[the derivation notes (PDF)](../notes/notes.pdf), chapter "Schwinger-Keldysh mode equation and the late-time mode function (Paper II)".
