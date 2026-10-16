# 10. Complete list of findings

<sub>[Overview (README)](../README.md) · Previous: [Schwinger–Keldysh and late times](09-schwinger-keldysh-late-time.md)</sub>

This page lists every result of the re-derivation of two papers:

- **Paper I**: E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007). It computes the one-loop quantum-gravity
  self-mass-squared of a massless, minimally coupled scalar on de Sitter space, in the gauge (I.15), which is not de Sitter
  invariant.
- **Paper II**: E. O. Kahya and R. P. Woodard, arXiv:0710.5282. It uses that self-mass in the Schwinger–Keldysh mode equation
  and studies the late-time behaviour of the scalar mode function.

The work was done step by step in 18 tasks, with Claude Opus 5.5 (Anthropic) as a checking tool. Each item below points to
the chapter of the [derivation notes](../notes/notes.pdf) where it is derived and to the scripts in [`code/`](../code/) that
check it. The full derivations, comparison tables and calculation records are in the notes.

**Contents:**
[At a glance](#at-a-glance) ·
[Summary for non-experts](#summary-for-non-experts) ·
[What changed and what did not](#what-changed-and-what-did-not) ·
[How to read the tables](#how-to-read-the-tables) ·
[Conventions](#conventions) ·
[Table A](#table-a-substantive-findings) ·
[Table B](#table-b-finite-errors-absorbed-into-counterterms) ·
[Table C](#table-c-typographical-errors) ·
[Table D](#table-d-confirmed-as-printed) ·
[Table E](#table-e-late-time-results-paper-ii) ·
[Limitations](#methodological-notes-and-limitations)

---

## At a glance

| Category | Number | IDs | Physical effect |
|---|---|---|---|
| [Substantive findings](#table-a-substantive-findings) | 4 | S1–S4 | None changes the late-time conclusion of Paper II |
| [Finite errors absorbed into counterterms](#table-b-finite-errors-absorbed-into-counterterms) | 6 | B1–B6 | None (finite parts are arbitrary) |
| [Typographical errors](#table-c-typographical-errors) | 23 | C1–C23 | None |
| [Confirmed as printed](#table-d-confirmed-as-printed) | 31 items | D1–D31 | — |
| [Late-time (Schwinger–Keldysh) results](#table-e-late-time-results-paper-ii) | 10 items | E1–E10 | Main conclusion of Paper II confirmed |

> **Main conclusion.** With the corrected self-mass, and without Paper II's simplification (II.44), the scalar mode function
> tends to its tree-level value times a constant:
>
> ```math
> \Phi=u+\kappa^2\Phi_1\ \xrightarrow{\ a\to\infty\ }\ u(0,k)\Big\{1+\frac{\kappa^2H^2}{16\pi^2}\Big[\tilde C+O\Big(\frac{\ln a}{a^2}\Big)\Big]\Big\},
> ```
>
> where $`\tilde C`$ is a field-strength renormalization constant. One-loop quantum gravity gives no significant correction
> to the scalar mode function at late times.

## Summary for non-experts

| Question | Paper I / II | This work | Does it matter physically? |
|---|---|---|---|
| **Does one-loop quantum gravity noticeably change the inflationary scalar's mode function at late times?** | **Paper II: no. The mode function tends to its tree-level value times a constant.** | **Confirmed**, with the corrected self-mass and without Paper II's shortcut (II.44): $`\Phi\to u(0,k)\{1+\frac{\kappa^2H^2}{16\pi^2}[\tilde C+O(\ln a/a^2)]\}`$. | **This is the main conclusion, and it stands.** |
| Overall sign of the computed self-mass | Printed as "local terms + Table X" | The printed result is $`-1`$ times the one defined by the paper's own (I.168). The physical nonlocal part is $`-`$(corrected Table X). | It flips only terms that die away at late times. No late-time effect. |
| Is the final table of nonlocal terms (Table X) right? | As printed | Two entries need fixing: a missing logarithmic term in one row, and a constant ($`-54\to+18`$) in another. | The log term changes the self-mass but only adds decaying terms (the source falls faster). The constant has no effect at all. |
| Finite counterterm numbers | $`\beta_2=-61/3`$, $`\beta_3=-4/(D-4)+58/3+2\gamma`$ | $`-49/3`$ and $`-4/(D-4)+52/3+2\gamma`$; the divergent part is unchanged. | No: these finite parts are arbitrary and are absorbed into the counterterms. |
| Which divergences appear, and which counterterms remove them? | Two invariant counterterms plus one de Sitter-breaking counterterm (I.145); $`\beta_1=0`$; only $`\beta_3`$ diverges | Confirmed. | Confirmed as printed. |
| Paper II's integral tables (Tables 2–5) | As printed | Tables 2–4 exact. Table 5 exact in the paper's sign; in the definition sign every entry flips, and the total is still 0. | Confirms the late-time cancellation. |
| Paper II's simplification (II.44) | Justified by a fall-off argument | The argument fails term by term, but the full calculation (through $`k^{10}`$) gives the same result. | No change to the conclusion. |
| Misprints | — | 23 typographical errors (signs, missing $`H^2`$, indices, appendix rows). | None. |

## What changed and what did not

We re-derived Paper I from the Lagrangian to the renormalized self-mass-squared, and the late-time analysis of Paper II.
The target equations were not used as input (the exceptions are listed under
[limitations](#methodological-notes-and-limitations), item 1).

**Confirmed.** Most of the calculation checks out:

- the Feynman rules (apart from misprints) and the coincidence limits;
- the 3-point split; Table III; Table IV (apart from one sign typo that the paper did not use in its calculation); all
  non-analytic terms of Table V; Tables VI–IX;
- the counterterm analysis: two invariant counterterms plus the single noninvariant one, (I.145);
- $`\beta_1=0`$ and the pole $`-4/(D-4)`$ of $`\beta_3`$;
- Tables 2–5 of Paper II.

**Changed.** Three things change in the self-mass.

1. **Overall sign (S1).** The printed self-mass is $`-1`$ times the one defined by (I.168). Both the 4-point part (the
   derivative rule in (I.32)) and the 3-point part (a sign slip in the step (I.46) → (I.47)–(I.49)) are reversed. Their
   relative sign, and the local cancellation stated after (I.67), are therefore consistent. The physical nonlocal part is
   $`-`$(corrected Table X).
2. **Table X (S2, S3).** Two corrections are needed. The $`(aa')^2\nabla^4`$ row lacks the Table III term
   $`8\ln x-16x\ln x`$ (non-analytic). The constant of the $`(aa')^3(a^2+a'^2)H^4`$ row is $`+18`$, not $`-54`$ (analytic).
   Only with both corrections does Table X match the direct sum of the normal contributions pointwise to 30 digits.
3. **Finite local errors (B1–B6).** Errors in (I.40) and (I.59) shift $`\beta_2:-\tfrac{61}{3}\to-\tfrac{49}{3}`$ and
   $`\beta_3:+\tfrac{58}{3}\to+\tfrac{52}{3}`$ (paper convention). The pole is unchanged, and the arbitrary finite
   counterterms absorb these shifts.

There are also 23 typographical errors with no effect (C1–C23).

**Paper II.** None of these corrections changes the late-time physics. We used the corrected self-mass in the definition
sign and did not use the simplification (II.44). The justification of (II.44) fails term by term (S4), although its
conclusion holds through order $`k^{10}`$. The mode function tends to
$`u(0,k)\{1+\frac{\kappa^2H^2}{16\pi^2}[\tilde C+O(\ln a/a^2)]\}`$. **The main conclusion of Paper II stands: one-loop
quantum gravity gives no significant correction to the scalar mode function at late times.**

## How to read the tables

- **Ch. N** is chapter N of the [derivation notes](../notes/notes.pdf). Each chapter has a summary page in `docs/` and a
  code folder in `code/`:

  | Ch. | Tasks | Summary page | Code |
  |---|---|---|---|
  | 1 | 1–2 | [Feynman rules](01-feynman-rules.md) | [`01_feynman_rules`](../code/01_feynman_rules/) |
  | 2 | 3–5 | [4-point contribution](02-four-point.md) | [`02_four_point`](../code/02_four_point/) |
  | 3 | 6–8 | [3-point setup](03-three-point-setup.md) | [`03_three_point_setup`](../code/03_three_point_setup/) |
  | 4 | 9–10 | [3-point local and log terms](04-three-point-local-and-log.md) | [`04_three_point_local_log`](../code/04_three_point_local_log/) |
  | 5 | 11–12 | [Normal terms and Table III](05-normal-terms-table-III.md) | [`05_normal_terms_tableIII`](../code/05_normal_terms_tableIII/) |
  | 6 | 13 | [Tables IV and V](06-tables-IV-V.md) | [`06_tables_IV_V`](../code/06_tables_IV_V/) |
  | 7 | 14–15 | [de Sitter symmetry and counterterms](07-de-sitter-symmetry-and-counterterms.md) | [`07_desitter_counterterms`](../code/07_desitter_counterterms/) |
  | 8 | 16–17 | [Renormalized self-mass and Table X](08-renormalized-self-mass-table-X.md) | [`08_renormalization_tableX`](../code/08_renormalization_tableX/) |
  | 9 | 18 | [Schwinger–Keldysh and late times](09-schwinger-keldysh-late-time.md) | [`09_schwinger_keldysh`](../code/09_schwinger_keldysh/) |

- A script name such as `check_4pt.py` refers to the code folder of the chapter cited next to it. Letters in parentheses,
  for example `check_41_44.py` (A),(D), are the labelled parts of that script's output.
- Equations of the papers are cited as (I.*n*) and (II.*n*). A prime, as in (I.33′), marks our corrected version of a
  printed equation.

## Conventions

Every number on this page is quoted together with its sign convention and units, as defined here.

- **Paper convention** (also called the **Table-X convention**): the sign of the self-mass contributions as Paper I prints
  them.
  - 4-point part: (I.32), (I.33), (I.41), (I.44);
  - 3-point part: (I.47)–(I.67), Tables III–X, (I.160)–(I.167);
  - Paper II: Table 1, (II.10), Table 5.

  The printed (I.45)–(I.46) agree with the definition sign. Identities such as (I.55)–(I.56) and (I.60)–(I.62) do not
  depend on the convention.
- **Definition sign**: the sign fixed by Paper I's own definition (I.168),

  ```math
  \partial_\mu\big(\sqrt{-g}\,g^{\mu\nu}\partial_\nu\Phi\big)-\int d^4x'\,M^2(x;x')\,\Phi(x')=0,\qquad -iM^2=\frac{\delta^2(i\Gamma)}{\delta\Phi\,\delta\Phi}.
  ```

  **Definition sign $`=-1\times{}`$ paper convention.** Physical statements are made in the definition sign.
- **Units.**
  - Table X and Paper II Table 1 entries are coefficients of $`\kappa^2H^4/(4\pi)^4`$, with $`x\equiv y/4`$.
  - $`\beta_i`$ and $`\alpha_i`$ are in units of $`H^{D-4}/(4\pi)^{D/2}`$.
  - $`\hat k\equiv k/H`$ and $`h\equiv H^2/16\pi^2`$.
- **One $`\kappa^2`$ convention** (Paper II, (II.23)). We write $`\Phi=u+\kappa^2\Phi_1`$, where $`u`$ is the tree-order
  mode function and $`\Phi_1`$ is sourced by the one-loop self-mass $`M^2_1`$ with $`\kappa^2`$ stripped off:

  ```math
  \Big[\partial_t^2+3H\partial_t+\frac{k^2}{a^2}\Big]\Phi_1=\mathcal R_1\equiv-\frac{1}{a^4}\int d\eta'\,u(\eta',k)\int d^3x'\,\big\{M^2_{++}+M^2_{+-}\big\}_1\,e^{i\vec k\cdot(\vec x\,'-\vec x)} .
  ```

  Table 5 entries and $`\mathcal R_1`$ are quoted in units of $`H^4u(0,k)/16\pi^2`$. The constant in the late-time result
  is $`\tilde C\equiv16\pi^2C/H^2`$, where $`C\,u(\eta,k)`$ is the homogeneous part of $`\Phi_1`$ (E2). The original
  Turkish Task 18 report wrote an extra explicit $`\kappa^2`$ inside $`\mathcal R`$ and inside the particular solution;
  the notes count $`\kappa^2`$ once, as here, and record the slip in a footnote of Ch. 9.

---

## Table A. Substantive findings

| ID | Where | Printed → corrected | Type | Effect at late times |
|---|---|---|---|---|
| **S1** | Overall sign of $`-iM^2`$: the 4-point rule (I.32) and the step (I.46) → (I.47)–(I.49) | the printed $`-iM^2`$ is $`-1`$ times the one defined by (I.168); the physical nonlocal part is $`-`$(corrected Table X) | substantive (overall sign) | none: the non-falling total is zero, so only falling terms change sign |
| **S2** | Table X, row $`(aa')^2\nabla^4`$ (the $`\zeta`$ row) | $`\tfrac{10}{3}\ln x-24x\ln x+\dots`$ → $`\tfrac{34}{3}\ln x-40x\ln x+\dots`$ | substantive (non-analytic) | only falling terms; the source then falls like $`\ln a/a^2`$ |
| **S3** | Table X, row $`(aa')^3(a^2+a'^2)H^4`$ ($`72\gamma_2`$) | constant $`-54`$ → $`+18`$ | substantive (analytic) | identically none |
| **S4** | Justification of (II.44) in Paper II | the fall-off argument fails term by term | argument only | none; the conclusion is verified through $`k^{10}`$ |

### S1. Overall sign of the self-mass-squared

- **Where.** 4-point origin: the rule (I.32)–(I.33). The paper reads $`\partial_2`$ as $`+\partial'`$ acting on
  $`\delta^D`$; the correct reading is $`+\partial`$ with respect to $`x`$, i.e. $`-\partial'`$. 3-point origin: the step
  (I.46) → (I.47)–(I.49), an overall sign slip in all nine terms. The 3-point rule itself ($`-\partial`$ on each outer leg)
  is correct. The sign is carried into (I.41), (I.44), (I.54)–(I.59), (I.63)–(I.67), (I.89), Tables III–X,
  (I.160)–(I.167), Paper II Table 1, (II.10) and Table 5.
- **Printed.** (I.33) coefficients $`(-\tfrac i8,+\tfrac i4,+\tfrac i2,-i)`$. (I.47)–(I.49): AA
  $`(+\kappa^2,-\kappa^2\tfrac{D-1}{D-3},+\kappa^2,+\kappa^2)`$, AB all $`-\kappa^2`$, AC $`+2\kappa^2\tfrac{D-2}{D-3}`$.
  (I.167): $`-iM^2_{\rm ren}={}`$local $`+`$ Table X.
- **Corrected (definition sign).** (I.33′) $`(+\tfrac i8,-\tfrac i4,-\tfrac i2,+i)`$; (I.47′)–(I.49′) $`={}-`$printed in all
  nine terms; and

  ```math
  -iM^2_{\rm ren}=i\kappa^2a^4\Big(\Delta\alpha_1\square^2+\Delta\alpha_2H^2\square+\Delta\alpha_3H^2\frac{\nabla^2}{a^2}\Big)\delta^4(x-x')-\big[\text{corrected Table X}\big].
  ```

- **Effect on physics.** Both parts are reversed, so the relative 4-point/3-point sign is consistent and the local
  cancellation of (I.44) against (I.67) holds. The sign of the local part is unphysical (the $`\Delta\alpha_i`$ are
  arbitrary). The nonlocal part has the opposite sign. In Paper II every Table 5 entry flips (the total is still 0) and the
  particular solution $`\Phi_1`$ flips sign. Since the non-falling total is zero, only falling terms are affected:
  **no late-time effect.**
- **Evidence.** Ch. 2, Task 3: Lemma A (two proofs), the field-equation route, a cross-test with the paper's counterterm
  vertex (I.131), test functions, the lattice Hessian and positivity;
  [`check_4pt.py`](../code/02_four_point/check_4pt.py) (a)–(d). Ch. 3, Task 7: two routes;
  [`check_47_49.py`](../code/03_three_point_setup/check_47_49.py). Ch. 4, Tasks 9–10: the paper's right-hand sides evaluate
  left-hand sides that carry the printed (I.47) signs; [`check_54_59.py`](../code/04_three_point_local_log/check_54_59.py)
  (A),(C), [`check_60_67.py`](../code/04_three_point_local_log/check_60_67.py) (f). Ch. 5:
  [`check_81_93.py`](../code/05_normal_terms_tableIII/check_81_93.py) (f). Ch. 8, Task 17: the definition-sign result.
  Ch. 9, Table 5′: [`table5.py`](../code/09_schwinger_keldysh/table5.py).

### S2. Table X, the $`(aa')^2\nabla^4`$ ($`\zeta`$) row

- **Where.** Paper I Table X, row $`(aa')^2\nabla^4`$ (= Paper II Table 1, row 12).
- **Printed.** $`\tfrac{10}{3}\ln x-24x\ln x+24x^2\ln x-36x^2`$.
- **Corrected (paper convention).** $`\tfrac{34}{3}\ln x-40x\ln x+24x^2\ln x-36x^2`$. The Table III $`\zeta`$-row
  contribution $`8\ln x-16x\ln x`$ was left out. Table III rows 1–2 were added to the $`(aa')^3H^2\nabla^2`$ and
  $`(aa')^2(a^2+a'^2)H^2\nabla^2`$ rows, and the text says that Table X includes Table III.
- **Effect on physics.** Non-analytic, so it enters the retarded kernel and changes the physical nonlocal self-mass. In
  Paper II it drops out under (II.44), since it is a $`\nabla^4`$ row. In the full Schwinger–Keldysh calculation it
  contributes only falling terms (paper convention):

  ```math
  \delta\mathcal R_1=\frac{H^4u(0,k)}{16\pi^2}\Big[\hat k^4\Big(\frac{4}{45a}-\frac{4}{3a^3}\Big)+\hat k^6\Big(\frac{4}{175a}-\frac{2}{15a^3}\Big)+\dots\Big],
  ```

  bounded by $`C\hat k^4/a`$ for every $`k`$. With the correction, the right-hand side of (II.29) falls like
  $`\ln a/a^2`$ instead of $`\hat k^4/a`$: the correction removes the $`\hat k^4/a`$, $`\hat k^6/a`$, … tails that the
  printed row leaves. **No effect on the late-time mode function.** The order-$`k^4`$ cancellation of the non-falling terms
  (E9) uses the $`24x^2\ln x`$ term of this row, which is the same in the printed and in the corrected row; the correction
  itself changes only falling terms.
- **Evidence.** Ch. 8, Task 17, the assembly of Table X:
  [`check_160_167.py`](../code/08_renormalization_tableX/check_160_167.py) (A) prints the missing
  $`-16x\ln x+8\ln x`$. Global test: [`check_tableX_global.py`](../code/08_renormalization_tableX/check_tableX_global.py)
  (see the [note below](#the-global-pointwise-test)). Ch. 5, Task 12, Table III:
  [`check_81_93.py`](../code/05_normal_terms_tableIII/check_81_93.py) (f). Ch. 9, the contribution of the $`\zeta`$
  correction: [`analysis2.py`](../code/09_schwinger_keldysh/analysis2.py),
  [`analysis3.py`](../code/09_schwinger_keldysh/analysis3.py), [`extract.py`](../code/09_schwinger_keldysh/extract.py),
  plus an analytic bound.

### S3. Table X, the constant of the $`(aa')^3(a^2+a'^2)H^4`$ row ($`72\gamma_2`$)

- **Where.** Found globally in Paper I Table V (Task 13). It appears in the Table X row $`(aa')^3(a^2+a'^2)H^4`$
  (= Paper II Table 1, row 7).
- **Printed.** Table X: $`\tfrac{4\ln x}{3x}-\tfrac{32}{3x}-54`$. Table V: $`\gamma_2`$ row $`-\tfrac5x-18`$; $`\delta`$ row
  without $`-12x`$.
- **Corrected (paper convention).** Table X: $`\tfrac{4\ln x}{3x}-\tfrac{32}{3x}+18`$. Equivalent Table V forms: add
  $`-12x`$ to the $`\delta`$ row ($`\delta[-12x]=72\gamma_2[1]`$, because $`\nabla^2\square x=-6aa'H^4`$), or make the
  $`\gamma_2`$ row $`-\tfrac5x+54`$. Where the term belongs is an interpretation: the pointwise test fixes only the sum.
- **Effect on physics.** Analytic ($`y`$-independent). It is identical in the $`++`$ and $`+-`$ kernels, so it drops out of
  the retarded kernel identically. **No effect on (II.21), Table 5 or the mode function.**
- **Evidence.** Ch. 6, Table V and hypothesis tests: [`check_tables45.py`](../code/06_tables_IV_V/check_tables45.py),
  [`check_tables45b.py`](../code/06_tables_IV_V/check_tables45b.py). The printed Tables IV + V miss the direct sum by
  0.0356, 17.05 and 0.607; the residual is zero to 30 digits only with the $`\gamma_2`$ sign fix (C15) plus
  $`72\gamma_2`$. Ch. 8, the global test. Ch. 9, the Schwinger–Keldysh differences.

### S4. The justification of (II.44)

- **Where.** Paper II, the text after (II.44).
- **Printed.** The deviations $`u(\eta')-u(0)`$ and $`e^{i\vec k\cdot\Delta\vec x}-1`$ bring at least two factors of
  $`1/a`$ or $`1/a'`$, so the $`\nabla^2`$ rows of Table 1 may be dropped.
- **Corrected.** The argument fails term by term. The rows $`(aa')^2(a^2+a'^2)H^2\nabla^2`$ and
  $`(aa')^2(a^2+a'^2)\nabla^2\square`$ each grow like $`\mp\tfrac{4}{15}\hat k^2a`$ (definition sign). The $`k^2`$
  corrections of the rows without derivatives give $`O(1)`$ constants. We kept the exact $`u(\eta',k)`$, the phase, the
  $`k^2/a^2`$ part of $`\square`$ and all $`\nabla^2`$ rows. The non-falling part then cancels separately at each order
  $`k^0,\dots,k^{10}`$, for the corrected and for the printed table.
- **Effect on physics.** The conclusion of (II.44) holds (verified through $`k^{10}`$); only its justification is invalid.
  No change to the physics.
- **Evidence.** Ch. 9, "Without the approximation (II.44)": the hand computation
  $`S_2[x\ln x]=-\tfrac{a}{180}+\tfrac1{24}+\dots`$; Table 5″. Scripts:
  [`kexp.py`](../code/09_schwinger_keldysh/kexp.py) ($`N=2,4,6,10`$),
  [`numcheck_k.py`](../code/09_schwinger_keldysh/numcheck_k.py), [`extract.py`](../code/09_schwinger_keldysh/extract.py).

### The global pointwise test

This test is the evidence for S2 and S3 (Ch. 8,
[`check_tableX_global.py`](../code/08_renormalization_tableX/check_tableX_global.py)).

- **Left-hand side:** the direct sum $`\sum_kP_k`$ of the normal contributions plus the nonlocal term of (I.67), with
  naive inner derivatives and the actual $`D=4`$ propagators. It does not depend on the appendix's outer-operator
  reduction.
- **Right-hand side:** the 13 operators of Table X acting on their coefficient functions.
- **Setup:** three random rational spacelike pairs (seed 29), 30-digit arithmetic. The size of the left-hand side is
  7.10, 0.107 and 8.60 at the three points.

| Table X used | Residual (point 1) | Residual (point 2) | Residual (point 3) |
|---|---|---|---|
| as printed | $`-11.52`$ | $`64.14`$ | $`-3.526`$ |
| with the $`\zeta`$ correction (S2) only | $`1.182`$ | $`65.28`$ | $`0.0161`$ |
| with the $`\zeta`$ and $`72\gamma_2`$ corrections (S2 and S3) | $`4\times10^{-31}`$ | $`0`$ | $`-3\times10^{-31}`$ |

Neither correction alone is enough. The test is numerical, at sample points; it is not a symbolic identity.

---

## Table B. Finite errors absorbed into counterterms

| ID | Where | Printed | Corrected | Effect | Evidence |
|---|---|---|---|---|---|
| **B1** | (I.40), finite part, and its consequences in (I.41), (I.44) | $`+2\,\delta^0_\alpha\delta^\sigma_0\,H^2/16\pi^2`$ | $`-2\,\delta^0_\alpha\delta^\sigma_0\,H^2/16\pi^2`$ (exact-$`D`$ form (I.40′)). The printed (I.40) contradicts the printed (I.38) through the trace identity ($`+2h`$ vs $`-2h`$). | With the correct (I.40) the $`B`$/$`C`$ finite parts are $`O(D-4)`$. Two printed pieces are then absent: the extra bracket $`+\tfrac{i\kappa^2H^2}{4\pi^2}\{-\partial^\mu(a^2\partial_\mu\delta^4)+a^2\nabla^2\delta^4\}`$ of (I.41), and the constants $`-2`$, $`+2`$ in the brackets of (I.44). The extra terms $`i\kappa^2k(-2a^D\square+2a^{D-2}\nabla^2)\delta^D`$ have counterterm form ((I.131), (I.145)), with $`\Delta\alpha_2=-\tfrac1{4\pi^2}`$ and $`\Delta\alpha_3=+\tfrac1{4\pi^2}`$ (paper convention). The correction shifts $`\beta_2`$ by $`+4`$ and $`\beta_3`$ by $`-4`$. | Ch. 2, Tasks 4–5: exact-$`D`$ contraction, components for $`D=4,\dots,7`$. [`check_coincidence.py`](../code/02_four_point/check_coincidence.py), [`check_41_44.py`](../code/02_four_point/check_41_44.py) (A),(D): the paper's (I.33) plus the printed (I.40) give exactly the paper's (I.41) and (I.44). |
| **B2** | (I.59), finite constant of the $`a^{D-2}\nabla^2`$ coefficient | $`-\tfrac72`$ (paper convention) | $`-\tfrac52`$ in the paper convention ($`+\tfrac52`$ in the definition sign). This is the sum of the paper's own (I.57) and (I.58): $`-\tfrac12-2`$. | Printed minus own sum $`=-\tfrac{i\kappa^2H^2}{8\pi^2}a^2\nabla^2\delta^4`$, which has the form of (I.145) with $`\Delta\alpha_3=-\tfrac1{8\pi^2}`$. The correction shifts $`\beta_3`$ by $`+2`$. | Ch. 4, Task 9: [`check_54_59.py`](../code/04_three_point_local_log/check_54_59.py) (C). |
| **B3** | $`\beta_2`$, (I.163) | $`-\tfrac{61}{3}`$ | $`-\tfrac{49}{3}`$ in the paper convention; $`+\tfrac{49}{3}`$ in the definition sign. | Finite; absorbed into $`\Delta\alpha_2`$. No pole either way. | Ch. 8, Task 17, "Effect of the earlier corrections on the β's": [`check_160_167.py`](../code/08_renormalization_tableX/check_160_167.py) (C). |
| **B4** | $`\beta_3`$, (I.165) | $`-\tfrac{4}{D-4}+\tfrac{58}{3}+2\gamma`$ | $`-\tfrac{4}{D-4}+\tfrac{52}{3}+2\gamma`$ in the paper convention; $`+\tfrac{4}{D-4}-\tfrac{52}{3}-2\gamma`$ in the definition sign. | **Pole unchanged**; the finite shift is absorbed into $`\Delta\alpha_3`$. | Same as B3. |
| **B5** | Paper II (II.8), $`\alpha_2`$ | $`\{\tfrac{61}{3}\}+\Delta\alpha_2`$ | $`\{\tfrac{49}{3}\}+\Delta\alpha_2`$ in the Table-X convention; $`\{-\tfrac{49}{3}\}+\Delta\alpha_2`$ in the definition sign. | Absorbed; enters no result. | Ch. 9, "Corrected counterterms". |
| **B6** | Paper II (II.9), $`\alpha_3`$ | $`\{\tfrac{4}{D-4}-\tfrac{58}{3}-2\gamma\}+\Delta\alpha_3`$ | $`\{\tfrac{4}{D-4}-\tfrac{52}{3}-2\gamma\}+\Delta\alpha_3`$ in the Table-X convention; $`\{-\tfrac{4}{D-4}+\tfrac{52}{3}+2\gamma\}+\Delta\alpha_3`$ in the definition sign. | Absorbed. The only local effect in the mode equation is $`\mathcal R_{1,\rm local}=-\Delta\alpha_3H^2k^2u/a^2\to0`$. | Same as B5. |

Units: $`\beta_i`$ and $`\alpha_i`$ in units of $`H^{D-4}/(4\pi)^{D/2}`$; $`h=H^2/16\pi^2`$.

**Notes on Table B.**

- (a) The printed (I.161)–(I.165) are reproduced **exactly** from the printed (I.44) and (I.59), so the paper did use those
  values. The shifts come only from correcting these two inputs.
- (b) $`\alpha_1=0+\Delta\alpha_1`$ ((II.7)) and $`\beta_1=0`$ are unchanged.
- (c) The outer $`\square`$ might act as a distribution, which could count a local term twice (see
  [limitations](#methodological-notes-and-limitations), item 3). Any such double counting would have exact counterterm
  form, $`c(D)(a^D\square-a^{D-2}\nabla^2)\delta^D`$, and would be absorbed into $`\alpha_2,\alpha_3`$. It is a potential
  issue, not an established error.

---

## Table C. Typographical errors

None of these affects any result. C1–C19, C22 and C23 are in Paper I; C20 and C21 are in Paper II.

| ID | Where | Printed | Correct | Note / evidence |
|---|---|---|---|---|
| C1 | (I.1), first term | $`-\tfrac12\partial_\mu\phi\,\partial_\mu\phi\,g^{\mu\nu}\sqrt{-g}`$ | $`-\tfrac12\partial_\mu\phi\,\partial_\nu\phi\,g^{\mu\nu}\sqrt{-g}`$ | Free index. Ch. 1, Task 1. |
| C2 | (I.6), left-hand side | $`-\tfrac12\partial_\mu\phi\,\partial_\nu g^{\mu\nu}\sqrt{-g}`$ | $`-\tfrac12\partial_\mu\phi\,\partial_\nu\phi\,g^{\mu\nu}\sqrt{-g}`$ | The right-hand side is correct. Noticed while preparing the notes; the Task 1 report retyped it as if correct. Ch. 1. |
| C3 | (I.7), left-hand side | $`-\tfrac12\kappa a^{D-2}\partial_\alpha\phi\,\partial_\beta\phi\,h^{\alpha\beta}`$ | $`+\tfrac12\kappa a^{D-2}\partial_\alpha\phi\,\partial_\beta\phi\,h^{\alpha\beta}`$ | $`V_1=i\kappa a^{D-2}\partial_1^\alpha\partial_2^\beta`$ is correct and was used in (I.46). Ch. 1, Task 2; Ch. 3, Task 6 ([`check_vertices.py`](../code/01_feynman_rules/check_vertices.py), [`check_45_46.py`](../code/03_three_point_setup/check_45_46.py)). |
| C4 | (I.14), first coefficient of $`A(y)`$ | $`\Gamma(\tfrac D2-1)/(\tfrac D2-1)`$ | $`\Gamma(\tfrac D2)/(\tfrac D2-1)=\Gamma(\tfrac D2-1)`$ | Required by (I.52), by (I.26), and by the recursion with the printed second coefficient. (I.52) carries the correct value. Ch. 3, Task 8 ([`check_52.py`](../code/03_three_point_setup/check_52.py) (e),(f)). |
| C5 | Table II, $`U_2`$ | $`-\tfrac i4\kappa^2a^{D-2}\eta^{\alpha\rho}\eta^{\beta\sigma}\partial_1\cdot\partial_2`$ | $`+\tfrac i4\kappa^2a^{D-2}\eta^{\alpha\rho}\eta^{\beta\sigma}\partial_1\cdot\partial_2`$ | The printed (I.33) uses $`+\tfrac i4`$. Ch. 1, Task 2; Ch. 2, Task 3 ([`check_vertices.py`](../code/01_feynman_rules/check_vertices.py), [`check_4pt.py`](../code/02_four_point/check_4pt.py) (c)). |
| C6 | (I.39), finite part | $`[2\delta^\rho_0\delta^\sigma_0+2\bar\eta_{\rho\sigma}]H^2/16\pi^2`$ | $`[2\delta^\rho_0\delta^\sigma_0+2\bar\eta^{\rho\sigma}]H^2/16\pi^2`$ | Index placement only (same matrix). Ch. 2, Task 4. |
| C7 | Text before (I.54) | "the 4th term of (47)" | "the 2nd term of (47)" | Ch. 4, Task 9. |
| C8 | (I.89), second term | $`-4(aa')^2(a+a')^2\nabla^2(yA')`$ | $`-4H^2(aa')^2(a+a')^2\nabla^2(yA')`$ (paper convention) | Dimensions; Table III has the $`H^2`$. Ch. 5, Task 12 ([`check_81_93.py`](../code/05_normal_terms_tableIII/check_81_93.py) (f)). |
| C9 | Appendix, Table XII header | $`\beta=(aa')^{D-2}(a^2+a'^2)H^2\square`$ | $`(aa')^{D-1}(a^2+a'^2)H^2\square`$, as in (I.105) | Weight counting; global check. Ch. 6. |
| C10 | Appendix, (A18) | $`-2(D-1)aa'H^2\,I[f'A]`$ | $`-2(D-1)aa'H^2\,I[f'A']`$ | Trace of (I.78). Ch. 6. |
| C11 | Appendix, Table XIII, row $`P_{7c}`$ | repeats part of the total | $`-\tfrac12(D-1)(D-2)\,y\,I[A'\Delta C']-\tfrac12(D-1)^2(D-2)\,I^2[A'\Delta C']`$ | The printed total is correct. Ch. 6. |
| C12 | Appendix, Table XV, row $`P_{4b}`$ | $`(D-1)\,y'(A'B')'`$ | $`(D-1)\,y^2(A'B')'`$ | From (A36). Ch. 6. |
| C13 | Appendix, Table XVIII, total | $`\tfrac14\,y\,I^2[A'^2\Delta C']`$ | $`\tfrac14\,y\,I^2[A'\Delta C']`$ | Row $`P_{7d}`$. Ch. 6. |
| C14 | (I.128), list of four-derivative invariants | five terms | the reduced list also contains $`\phi^2R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}`$ and $`\phi^2\square R`$ | Omission. The shift symmetry removes both, so (I.129) is unaffected. Ch. 7, Task 14. |
| C15 | Table IV, $`\gamma_2`$ row, second column | $`-\tfrac{D^3-3D^2-4D+8}{4(D-1)}`$ | $`+\tfrac{D^3-3D^2-4D+8}{4(D-1)}`$ | The paper's Tables VI–IX use $`+`$, and shift symmetry requires $`+`$: the printed sign would leave a forbidden $`-\tfrac{16}{3}`$ in Table IX at $`D=4`$. Ch. 6, Ch. 8 ([`check_table4.py`](../code/06_tables_IV_V/check_table4.py), [`check_146_IX.py`](../code/08_renormalization_tableX/check_146_IX.py) (D),(G)). |
| C16 | Table VII, header | "Coefficient of $`i\kappa^2H^4/(4\pi)^4`$" | "Coefficient of $`\kappa^2H^4/(4\pi)^4`$" | Table IV units carry no $`i`$, and the Table X header has none. Ch. 8. |
| C17 | (I.160) | $`\beta_2\square+\beta_3\nabla^2/a^2`$ | $`\beta_2H^2\square+\beta_3H^2\nabla^2/a^2`$ | Dimensions. (I.131) and (I.145) carry $`H^2`$, and (I.166) needs $`\alpha_i=-\beta_i`$. Ch. 8. |
| C18 | (I.167) | $`\Delta\alpha_2\square+\Delta\alpha_3\nabla^2/a^2`$ | $`\Delta\alpha_2H^2\square+\Delta\alpha_3H^2\nabla^2/a^2`$ | Same as C17. Ch. 8. |
| C19 | Table X | row $`(aa')^2(a^2+a'^2)\nabla^2\square`$ printed twice | printed once | Paper II Table 1 lists it once. Ch. 8. |
| C20 | Paper II, (II.10) | $`\Delta\alpha_2\square+\Delta\alpha_3\nabla^2/a^2`$ | $`\Delta\alpha_2H^2\square+\Delta\alpha_3H^2\nabla^2/a^2`$ | Same omission as C17–C18. Flagged in the write-up; the Task 18 report restored $`H^2`$ without comment. Ch. 9. |
| C21 | Paper II, (II.34) | $`\Delta\alpha`$ term without $`H^2`$ | with $`H^2`$ | Same as C20. Ch. 9. |
| C22 | (I.170) | $`\Delta\alpha_2\square+\Delta\alpha_3\nabla^2/a^2`$ (and the result $`-\kappa^2a^4\Delta\alpha_3\,k^2/a^2\,\Phi_0`$) | with $`H^2`$: $`\Delta\alpha_2H^2\square+\Delta\alpha_3H^2\nabla^2/a^2`$, and the result $`-\kappa^2a^4\Delta\alpha_3H^2k^2/a^2\,\Phi_0`$ | Same omission as C17–C18. Noticed while refereeing the write-up; confirmed on the printed page 124005-15. |
| C23 | Appendix, text before (A31) | "Applying the second identity (A18) to $`P_1`$ and $`P_6`$" | "… to $`P_1`$ and $`P_4`$" (the subparts that follow are $`P_{1a},P_{1b},P_{4a},\dots`$) | Wording only. Noticed while refereeing the write-up; confirmed on the printed page 124005-17. |

**Notational remarks (not counted as errors).**

- (I.6): $`h^\nu_\rho`$ is unstaggered. This is unambiguous, since $`h_{\mu\nu}`$ is symmetric.
- (I.71): $`\Delta x_\sigma`$ must be read as $`\eta_{\sigma\nu}\Delta x^\nu`$, so $`\Delta x_0=-\Delta\eta`$.
- (I.69): "never a delta function" holds for naive (pointwise) derivatives. The distributional derivative adds exactly the
  local term that is counted separately in (I.54)–(I.58). So inner derivatives must be naive, as the footnote to (I.75)
  says.
- (I.92)–(I.93): the $`\square`$ in (I.92) is distributional; (I.93) is pointwise.
- (I.125): the finite bracket is defined only up to the constant $`(24\ln2-36)\ln2`$, which $`\alpha=(aa')^D\square^2`$
  annihilates.
- Table V: the polynomial and $`\ln2`$ terms depend on integration constants and cancel globally, except $`72\gamma_2`$ (S3).
- (I.134): the dilatation parameter needs $`k>0`$.
- (I.135): $`x\cdot x`$ is the Minkowski product $`\Vert\vec x\Vert^2-\eta^2`$.

---

## Table D. Confirmed as printed

| ID | Item | Status | Evidence |
|---|---|---|---|
| D1 | (I.4), (I.5) and the bracket of (I.6) | Exact, general $`D`$, no gauge or trace condition | Ch. 1 ([`check_eq6.py`](../code/01_feynman_rules/check_eq6.py)) |
| D2 | Table I ($`V_1`$, $`V_2`$), (I.8), Table II rows $`U_1,U_3,U_4`$ | Exact | Ch. 1 ([`check_vertices.py`](../code/01_feynman_rules/check_vertices.py)) |
| D3 | (I.31): structure and normalization 1 | Exact | Ch. 2, Task 3 |
| D4 | (I.34)–(I.38); (I.39) up to C6; $`i\Delta_A`$ coefficient of (I.40) | Exact in $`D`$; $`D=4`$ finite parts $`+4h`$, $`-2h`$, $`[2,2]h`$ | Ch. 2, Task 4 ([`check_coincidence.py`](../code/02_four_point/check_coincidence.py)) |
| D5 | (I.42), (I.43) (plus a general-$`D`$ form); the algebra (I.41) → (I.44) | Correct, given the printed inputs | Ch. 2, Task 5 ([`check_41_44.py`](../code/02_four_point/check_41_44.py)) |
| D6 | (I.45), (I.46) | Identical (printed = definition sign) | Ch. 3, Task 6 ([`check_45_46.py`](../code/03_three_point_setup/check_45_46.py)) |
| D7 | (I.47)–(I.49): structures, magnitudes and $`D`$-dependence of all nine terms | Agree (overall sign: S1) | Ch. 3, Task 7 |
| D8 | (I.50)–(I.53) and the local/log/normal split; (I.52) along three independent routes | Correct; the split is complete and has no overlap | Ch. 3, Task 8 ([`check_52.py`](../code/03_three_point_setup/check_52.py)) |
| D9 | (I.55)–(I.56), (I.60)–(I.62) | Exact (convention-independent) | Ch. 4 |
| D10 | (I.63)–(I.67): magnitudes (3, 4, 9, 2) and cancellation of $`\nabla^4A`$. Local cancellation of $`3\ln(aa')a^D\square`$, $`-4\ln(aa')a^{D-2}\nabla^2`$ and $`9H^2a^D`$ between (I.44) and (I.67) | Confirmed at general $`D`$ (remainder $`\propto(D-4)(D+1)/(4(D-3))`$); holds in both conventions | Ch. 4, Task 10 ([`check_60_67.py`](../code/04_three_point_local_log/check_60_67.py) (e)) |
| D11 | (I.68)–(I.88), (I.90)–(I.93); divergence count after (I.80) | Exact | Ch. 5 ([`check_68_80.py`](../code/05_normal_terms_tableIII/check_68_80.py), [`check_81_93.py`](../code/05_normal_terms_tableIII/check_81_93.py)) |
| D12 | Table III | Equals the printed table in the paper convention | Ch. 5, Task 12 |
| D13 | (I.94)–(I.124); (I.125) up to a constant | Correct | Ch. 6 ([`check_94_125.py`](../code/06_tables_IV_V/check_94_125.py)) |
| D14 | Appendix Tables XI–XX, totals | Exact after C9–C13, via the global pointwise test (an audit; see [limitations](#methodological-notes-and-limitations), item 2) | Ch. 6 ([`check_global.py`](../code/06_tables_IV_V/check_global.py)) |
| D15 | Table IV | 19 of 20 entries in general $`D`$ (C15 is a typo) | Ch. 6 ([`check_table4.py`](../code/06_tables_IV_V/check_table4.py)) |
| D16 | Table V | All logarithmic and inverse-power terms ($`1/x`$, $`\ln x/x`$, $`\ln x`$, $`\ln^2x`$, $`x\ln x`$, $`x^2\ln x`$) | Ch. 6 ([`check_tables45.py`](../code/06_tables_IV_V/check_tables45.py)) |
| D17 | (I.126)–(I.145): the two invariant counterterms (I.129) and their vertices (I.130)–(I.131); the de Sitter group with $`\tfrac12D(D+1)`$ parameters in the four classes (I.132)–(I.135); the gauge (I.15) breaks only the spatial special conformal transformations; exactly three noninvariant counterterms (I.144), of which only the last survives the flat limit; the vertex (I.145); the statements after (I.145) | All reproduced (remarks: C14 and the notational remarks under Table C) | Ch. 7 ([`check_126_135.py`](../code/07_desitter_counterterms/check_126_135.py), [`check_136_145.py`](../code/07_desitter_counterterms/check_136_145.py)) |
| D18 | (I.146)–(I.159), Tables VI–IX | Exact (with the $`+`$ sign of $`\gamma_2`$, which the paper used) | Ch. 8, Task 16 ([`check_146_IX.py`](../code/08_renormalization_tableX/check_146_IX.py)) |
| D19 | Vanishing $`a^D\square^2`$ total of Table VIII, hence $`\beta_1=0`$ | Confirmed. Independent argument: the flat-space numerator $`V\cdot P\cdot V=\tfrac12p^2k^2`$ gives a scaleless integral | Ch. 8 (check (I)) |
| D20 | Table IX is finite; (I.157)–(I.158) $`\to\frac{\kappa^2H^4}{8\pi^2}ia^4\delta^4`$ | Confirmed | Ch. 8 |
| D21 | (I.161)–(I.165) from the printed inputs: $`\beta_1=0`$, $`\beta_2=-\tfrac{61}{3}`$ (no pole), $`\beta_3=-\tfrac{4}{D-4}+\tfrac{58}{3}+2\gamma`$ | Reproduced exactly. The pole $`-4/(D-4)`$ (paper convention) survives the corrections. The only divergence is in the de Sitter-breaking $`\beta_3`$ | Ch. 8, Task 17 ([`check_160_167.py`](../code/08_renormalization_tableX/check_160_167.py) (C)) |
| D22 | (I.166), $`\alpha_i=-\beta_i+\Delta\alpha_i`$ | Correct (convention-independent) | Ch. 8 |
| D23 | Table X | 12 of 13 rows follow exactly from Tables III, V, VII and (I.159); 11 rows are unchanged after the corrections | Ch. 8 ([`check_160_167.py`](../code/08_renormalization_tableX/check_160_167.py) (A)) |
| D24 | Paper II (II.15)–(II.21): SK rules and effective equation; reality and causality of the retarded kernel | Derived (written in the definition sign) | Ch. 9 |
| D25 | (II.29) and the particular solutions (II.30)–(II.33) (leading terms) | Reproduced | Ch. 9 |
| D26 | (II.34): only $`\Delta\alpha_3`$ contributes locally | Confirmed (C21 aside) | Ch. 9 |
| D27 | (II.35)–(II.41) | Confirmed; $`1/x`$ also computed directly from the light-cone $`\delta(y)`$ | Ch. 9 |
| D28 | Tables 2–4 (12 entries) | Identical through the order kept; next order obtained; quadrature agreement $`\lesssim10^{-19}`$ | Ch. 9 ([`tables234.py`](../code/09_schwinger_keldysh/tables234.py), [`numcheck.py`](../code/09_schwinger_keldysh/numcheck.py)) |
| D29 | Table 5 (19 entries), total 0 | Identical in the paper (Table-X) convention (sign: S1) | Ch. 9 ([`table5.py`](../code/09_schwinger_keldysh/table5.py)) |
| D30 | Statement after Table 5: the right-hand side falls at least like $`1/a`$ times powers of $`\ln a`$ | Consistent (corrected table: $`\ln a/a^2`$; printed table: $`\hat k^4/a`$) | Ch. 9 |
| D31 | Conclusion of Paper II: no significant one-loop late-time correction | Confirmed (E1) | Ch. 9 |

---

## Table E. Late-time results (Paper II)

All results in this table are in the definition sign, with $`\kappa^2`$ counted once (see [Conventions](#conventions)).

**E1. Late-time mode function (the main physics conclusion).**

```math
\Phi=u+\kappa^2\Phi_1\ \xrightarrow{\ a\to\infty\ }\ u(0,k)\Big\{1+\frac{\kappa^2H^2}{16\pi^2}\Big[\tilde C+O\Big(\frac{\ln a}{a^2}\Big)\Big]\Big\},
```

where $`\tilde C`$ is a field-strength renormalization constant. There is no growing ($`a`$, $`\ln a`$, $`\ln^2a`$) or
drifting term. **There is no significant one-loop correction at late times.** This confirms Paper II with the corrected
input and without the approximation (II.44). (Ch. 9, "Late-time effect on the mode function".)

**E2. The one-loop correction** ($`\kappa^2`$ stripped off):

```math
\Phi_1=C\,u(\eta,k)+\frac{H^2u(0,k)}{16\pi^2}\Big\{\frac1{a^2}\Big[4\hat k^2\ln a-6+\tfrac{91}{6}\hat k^2+4i\hat k^3+\dots\Big]+\frac{4i\hat k^3\ln a}{a^3}+\dots\Big\}.
```

$`C\,u`$ is a homogeneous solution that the late-time analysis does not fix, with $`C=\tfrac{H^2}{16\pi^2}\tilde C`$. The
initial condition removes the other homogeneous solution, $`u^*`$. (Ch. 9, "Solving for $`\Phi_1`$".)

**E3. The late-time source** with the corrected table, through $`k^{10}`$:

```math
\mathcal R_1=\frac{H^4u(0,k)}{16\pi^2}\Big\{\frac1{a^2}\Big[12-\Big(8\ln a+\tfrac{103}{3}\Big)\hat k^2-8i\hat k^3+\tfrac83\hat k^4+\tfrac{32i}{45}\hat k^5-\tfrac2{15}\hat k^6+\dots\Big]-\frac{12i\hat k^3}{a^3}+O\Big(\frac{\ln a}{a^4}\Big)\Big\}.
```

**There is no $`a^{-1}`$ term.** With the printed $`\zeta`$ row there would be the extra terms

```math
\frac{H^4u(0,k)}{16\pi^2}\Big[\tfrac{4}{45}\hat k^4+\tfrac{4}{175}\hat k^6+\tfrac{2i}{135}\hat k^7-\tfrac{8}{1323}\hat k^8+\dots\Big]\frac1a .
```

The $`a^{-2}`$ coefficients depend on the initial time and on $`\hat k`$, i.e. they are sensitive to the initial state.
(Ch. 9, [`extract.py`](../code/09_schwinger_keldysh/extract.py).)

| ID | Result | Evidence |
|---|---|---|
| E4 | Local part: $`\mathcal R_{1,\rm local}=-\Delta\alpha_3H^2k^2u(\eta,k)/a^2\to0`$. The corrected finite parts of $`\beta_{2,3}`$ only redefine $`\Delta\alpha_{2,3}`$. | Ch. 9 |
| E5 | SK polarities built from the $`++`$ result: $`-iM^2_{++}={}`$local $`-\,\mathcal T(y_{++})`$, $`-iM^2_{+-}=+\mathcal T(y_{+-})`$, $`-iM^2_{--}=(-iM^2_{++})^*`$, $`-iM^2_{-+}=(-iM^2_{+-})^*`$. Here $`\mathcal T`$ is the corrected Table X (= Paper II Table 1) in the paper convention. The retarded kernel is $`M^2_{++}+M^2_{+-}={}`$local $`+\,2\,{\rm Im}\,\mathcal T(y_{++})`$ for $`\eta'<\eta`$; its nonlocal part is $`0`$ for $`\eta'>\eta`$. It is real and causal (support on and inside the past light cone). This is Paper II's rule below (II.19), written in the definition sign. | Ch. 9, "Schwinger-Keldysh self-masses" |
| E6 | SK differences: $`\Delta\ln x=2\pi i\theta`$, $`\Delta\ln^2x=4\pi i\theta\ln\vert x\vert`$, $`\Delta(x^j\ln x)=2\pi i\,x^j\theta`$, $`\Delta(1/x)=-8\pi i\,\delta(y_r)`$, $`\Delta(\text{analytic})=0`$. Hence the $`72\gamma_2`$ correction and the terms $`-120`$, $`96x`$, $`12x`$, $`-36x^2`$ drop out identically. | Ch. 9, (II.39)–(II.41) |
| E7 | Tables 2–4: all 12 entries in exact closed form, identical to Paper II; next order obtained. The largest difference from quadrature is $`1.0\times10^{-19}`$, at $`a=50`$ and $`a=400`$. | [`tables234.py`](../code/09_schwinger_keldysh/tables234.py), [`numcheck.py`](../code/09_schwinger_keldysh/numcheck.py) |
| E8 | Table 5′ (within (II.44)): each of the 19 entries is $`-1`$ times the printed one. For example, the printed $`-18,\ 2,\ 3,\ -9,\ 33-18\ln a`$ become $`+18,\ -2,\ -3,\ +9,\ -33+18\ln a`$. The total is $`0`$, with the $`\ln^2a`$, $`\pi^2`$, $`\ln a`$ and constant parts cancelling separately. The $`\zeta`$ and $`72\gamma_2`$ corrections do not enter Table 5′. | [`table5.py`](../code/09_schwinger_keldysh/table5.py) |
| E9 | Beyond (II.44): the non-falling part vanishes at each order $`k^0,\dots,k^{10}`$, for the corrected and for the printed table. Rows 10 and 11 each grow like $`\mp\tfrac{4}{15}\hat k^2a`$ and cancel each other. The $`\zeta`$ row contributes $`-\hat k^4/70`$ at order $`k^4`$ and is needed for the $`k^4`$ cancellation. This contribution comes from its $`24x^2\ln x`$ term, which is the same in the printed and in the corrected row. The $`k^8`$ series agrees with exact-$`k`$ quadrature at $`a=200`$ to within $`10^{-10}`$ to $`2\times10^{-8}`$ at $`\hat k=0.3`$, and $`5\times10^{-6}`$ to $`9.3\times10^{-4}`$ at $`\hat k=1`$. These differences are the expected $`k^9`$ truncation error. | [`kexp.py`](../code/09_schwinger_keldysh/kexp.py), [`numcheck_k.py`](../code/09_schwinger_keldysh/numcheck_k.py), [`extract.py`](../code/09_schwinger_keldysh/extract.py) |
| E10 | Effect of each correction on the late-time result. Overall sign (S1): every Table 5 entry flips and $`\Phi_1\to-\Phi_1`$; not observable at late times, since the non-falling total is zero and only falling terms change sign. $`\zeta`$ row (S2): only falling terms ($`\le C\hat k^4/a`$ for every $`k`$); it removes the $`1/a`$ tails. $`72\gamma_2`$ (S3): identically none. $`\beta_{2,3}`$ shifts (B3–B6): none. Duplicated row (C19): none, since it is already single in Paper II Table 1. | Ch. 9, the "Answer to the question" table |

---

## Methodological notes and limitations

1. **How the work was done.** Each step was re-derived from its inputs, with Claude (an AI model) as a checking tool.
   - The analytic derivations are kept separate from the code checks.
   - $`D`$ is kept general wherever the working reports kept it general.
   - Exceptions to "the target is not used as input": Tables IV–V were computed from the printed appendix totals after
     those totals were checked. Table V enters the Table X assembly as printed; its non-analytic terms were confirmed in
     Task 13.
   - Tables 2–5 of Paper II were never used as input.
   - All scripts in [`code/`](../code/) were re-run and reproduce the original logged outputs. The one cosmetic difference
     is the order of two tied entries in `check_68_80.py`, now fixed by a sort key.
2. **The appendix "Extracting derivatives" was audited, not redone.**
   - Worked examples (the $`P_{3c}`$/$`\gamma_3`$ row, the $`\beta`$/$`\gamma_1`$ totals) were checked by hand.
   - Tables XII and XIV and (A1)–(A9) were read from the printed pages.
   - The other tables were used through their totals, which the global pointwise test verified (generic powers at
     $`D=4,5`$; actual $`D=4`$ propagators).
   - The Task 13 report's statement that every row was reproduced by hand was an overstatement.
3. **A possible extra local term from the appendix operator reduction is only partially resolved.**
   - If the outer $`\square`$ in $`\alpha`$ and $`\delta`$ acts as a distribution, the local term (I.54)+(I.58) could be
     counted twice.
   - Task 16 showed that any such double counting has exactly counterterm form,
     $`c(D)(a^D\square-a^{D-2}\nabla^2)\delta^D`$ (no $`\ln a`$), and is absorbed into $`\alpha_2,\alpha_3`$.
   - It was not shown that (A35)–(A38) produce no other local difference at coincidence (for example a finite
     $`H^4a^D\delta^D`$). That needs a full distributional analysis.
   - Partial support: the vanishing $`a^D\square^2`$ total and the finiteness of Table IX.
4. **Global checks are numerical, at sample points.**
   - The Table V/appendix test and the Table X test use 2–3 random spacelike points at 30-digit precision. They are not
     symbolic identities, and they check only sums.
   - Some other checks are dimension-specific:
     - the gauge-breaking Euler–Lagrange check is done in $`D=3`$;
     - counterterm completeness was checked for all 15 candidates only at $`D=4,5,6`$. General-$`D`$ formulas were derived
       for 11 of them, contrary to the Task 15 report's wording;
     - the pointwise check of the nonlocal part of (I.60)–(I.67) is done at $`D=4`$.
5. **Where $`72\gamma_2`$ belongs** (a dropped $`-12x`$ in the $`\delta`$ row of Table V) is an interpretation. The tests
   fix only the sum.
6. **The Schwinger–Keldysh cancellation is exact through $`k^{10}`$, but not proven to all orders.**
   - Beyond $`k^{10}`$ only the convergent-series argument is available (entire integrand, compact region).
   - The $`a^{-2}`$ coefficients of $`\mathcal R_1`$ and $`\Phi_1`$ are sensitive to the initial state.
   - The explanation of the $`k=0`$ cancellation ($`M^2=\partial\partial'K`$) is a heuristic, labelled as an assumption.
   - The operator ordering (scale factors on the left) is assumed. It was verified in the Task 17 pointwise test.
   - The direct light-cone evaluation of $`S_K[1/x]`$ was done by hand only.
7. **Assumptions shared with the paper:**
   - dimensional regularization preserves the shift symmetry and the residual symmetries of the gauge-fixed action;
   - ultraviolet divergences are polynomial in $`H`$;
   - coincidence limits are taken in dimensional regularization;
   - $`D=4`$ is set only in terms not multiplied by a pole.
8. **The overall sign rests on the paper's own definition** (I.168), with
   $`-iM^2=\delta^2(i\Gamma)/\delta\Phi\,\delta\Phi`$. Several independent lines support it: Lemma A, the field-equation
   route, the counterterm vertex (I.131), test functions, the lattice Hessian, and the consistency of the 3-point steps.
   Once (I.168) is adopted the sign of the nonlocal part is physical, because the nonlocal part has no counterterm freedom.
9. **Corrections made to the Turkish working reports in the English write-up.** None changes a result.
   - The Task 1 report retyped the left-hand side of (I.6) as if it were correct (C2).
   - The Task 6 report cited (I.174) for the flat 3-point expression; the correct equation is (I.175).
   - Tasks 7 and 9 later revised Task 6's relative-sign conclusion. This is a revision between reports, not an error
     inside one.
   - The Task 12 remark that $`\square y^{1-D/2}`$ vanishes pointwise is inexact:
     $`\square y^{1-D/2}=\tfrac{D(D-2)}{4}H^2y^{1-D/2}`$. The conclusion is unaffected.
   - The Task 13 and Task 15 overstatements (items 2 and 4 above).
   - In Task 18: $`\kappa^2`$ written twice and the $`C`$ vs $`\tilde C`$ notation (see [Conventions](#conventions)). The
     quadrature agreement was quoted as $`10^{-20}`$, $`\lesssim10^{-19}`$ and $`\le1.2\times10^{-19}`$ in different
     places; the logged maximum is $`1.02\times10^{-19}`$.
10. **Findings first noticed during the write-up**, marked as such in the notes: C2 (the left-hand side of (I.6)) and the
    missing $`H^2`$ in (II.10)/(II.34) (C20–C21). C22 and C23 were noticed while refereeing the write-up and confirmed on
    the printed pages.

---

<sub>See also: the [overview (README)](../README.md), the per-stage pages [01](01-feynman-rules.md)–[09](09-schwinger-keldysh-late-time.md),
the [derivation notes (PDF)](../notes/notes.pdf), whose last chapter, "Summary of findings", gives a pointer to the section
where each item is derived, and the [code README](../code/README.md).</sub>
