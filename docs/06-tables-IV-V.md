# 6. Tables IV and V: Eqs. (I.94)–(I.125)

<sub>[Overview (README)](../README.md) · Previous: [Normal terms and Table III](05-normal-terms-table-III.md) · Next: [de Sitter symmetry and counterterms](07-de-sitter-symmetry-and-counterterms.md)</sub>

*Task 13 of the re-derivation of Paper I (Kahya and Woodard, Phys. Rev. D 76, 124005 (2007)).*

## What this step computes

The "normal" part of the 3-point contribution consists of terms in which outer derivatives act on one propagator function times inner derivatives of another, for example $`\partial_0\partial'_0\big[(aa')^{D-2}A\,\partial_0\partial'_0A\big]`$. Paper I rewrites every such term so that all derivatives act *outside*, as one of ten **external operators** acting on a function of the de Sitter length $`y`$ alone:
$`\alpha,\ \beta,\ \gamma_1,\ \gamma_2,\ \gamma_3,\ \delta,\ \epsilon_1,\ \epsilon_2,\ \epsilon_3,\ \zeta`$.
The bookkeeping is done in the appendix of the paper ("Extracting derivatives"): the normal contribution is split into seven parts $`P_1,\dots,P_7`$, and Tables XI–XX list what each operator acts on. Substituting the propagator expansions then gives **Table IV** (the divergent coefficients, in general $`D`$) and **Table V** (the finite remainder at $`D=4`$).

We re-derived the identities (I.94)–(I.125) in general $`D`$, rebuilt Table IV and Table V, and tested the appendix decomposition (the ten table totals) pointwise against the direct sum of $`P_1,\dots,P_7`$. Almost everything agrees. Two things in Tables IV–V do not: one sign in Table IV (C15), and one finite term missing from Table V (S3). The appendix also has five typographical errors in its tables and equations (C9–C13). The IDs refer to the [complete list of findings](10-complete-list-of-findings.md).

## Key steps

**1. The ubiquitous reduction (I.101).** Any single mixed time derivative is split with $`\partial_0\partial'_0=\tfrac12(\partial_0+\partial'_0)^2-\tfrac12\partial_0^2-\tfrac12\partial_0'^2`$, and the identities of the previous chapter reduce each piece to $`\square`$, $`\nabla^2`$ and ordinary $`y`$-derivatives:

```math
\begin{aligned}
\partial_0\partial'_0\big[(aa')^{D-2}f(y)\big]={}&\tfrac12(aa')^{D-2}(a^2+a'^2)\big[\square f+H^2yf'\big]-(aa')^{D-2}\nabla^2f
+\tfrac14(D-2)(aa')^{D-2}\nabla^2I[f]\\
&+\tfrac12(D-2)(D-3)H^2(aa')^{D-1}f+\tfrac12H^2(a+a')^2(aa')^{D-2}\big[(D-1)yf'+y^2f''\big].
\end{aligned}
```

Here $`I[f]`$ is the indefinite integral of $`f`$ with respect to $`y`$. Similar manipulations with the same identities give (I.94)–(I.96) for the mixed $`A`$ terms and (I.102)–(I.103) for the two $`B`$ terms.

**2. The ten external operators (I.104)–(I.113)**, for example

```math
\alpha=(aa')^D\square^2,\qquad \beta=(aa')^{D-1}(a^2+a'^2)H^2\square,\qquad
\gamma_2=(aa')^{D-1}(a^2+a'^2)H^4,\qquad \delta=(aa')^{D-2}(a^2+a'^2)\nabla^2\square .
```

**3. The appendix decomposition.** With $`\mathcal T_{\rm op}`$ the total of the corresponding appendix table,

```math
\sum_{k=1}^{7}P_k=\alpha\,\mathcal T_\alpha+\beta\,\mathcal T_\beta+\sum_{i=1}^{3}\gamma_i\,\mathcal T_{\gamma_i}+\delta\,\mathcal T_\delta+\sum_{i=1}^{3}\epsilon_i\,\mathcal T_{\epsilon_i}+\zeta\,\mathcal T_\zeta .
```

All totals are built from eight products of propagator functions, $`A'^2,\ AA'',\ A'B',\ A''B,\ A'\Delta B',\ A''\Delta B,\ A'\Delta C',\ A''\Delta C`$, with $`\Delta B=B-A`$ and $`\Delta C=\tfrac{2(D-2)}{D-3}(C-A)`$.

**4. Leading powers in exact $`D`$.** Only the two leading terms of each propagator reach the powers $`z^D`$ and $`z^{D-1}`$ of the products, which are the ones that can give divergences. With $`z=4/y`$, for example

```math
A'B'=\frac{\Gamma^2(\frac D2)}{16}\frac{H^{2D-4}}{(4\pi)^D}\Big\{z^D+(D-2)\,z^{D-1}+\dots\Big\}.
```

**5. Table IV from simple rules.** The coefficients of $`z^{D-1}`$ and $`z^{D-2}`$ in each table total follow from these leading coefficients by the rules $`I[z^q]=-\tfrac{4z^{q-1}}{q-1}`$, $`y\,z^q=4z^{q-1}`$ and $`\tfrac{d}{dy}z^q=-\tfrac q4z^{q+1}`$. For the second column of the $`\gamma_2`$ row:

```math
\frac1{16}\Big[4D(D-2)+\frac{8D}{D-1}-32\Big]=+\frac{D^3-3D^2-4D+8}{4(D-1)} .
```

The paper prints this entry with a minus sign.

**6. Table V and the missing finite term.** At $`D=4`$ we computed every operator function exactly and compared it with the printed Tables IV + V. Apart from the $`\gamma_2`$ entry of Table IV (step 5), all logarithmic and inverse-power terms agree. The other differences are polynomial; they depend on the choice of integration constants, and each is either annihilated by its operator or cancels in the sum over the ten operators, except one. That one is a term $`-12x`$ which our $`\delta`$ function contains and the printed $`\delta`$ row lacks. It survives because at $`D=4`$

```math
\nabla^2\square x=-6aa'H^4\neq0\qquad\Longrightarrow\qquad \delta[-12x]=72\,(aa')^3(a^2+a'^2)H^4=72\,\gamma_2[1],\qquad x=\tfrac y4 .
```

## Result

> **Table IV**, row $`\gamma_2`$, second column: $`+\dfrac{D^3-3D^2-4D+8}{4(D-1)}`$ (printed with $`-`$). A later step ([Task 16](08-renormalized-self-mass-table-X.md)) showed that the paper's Table VI already uses the $`+`$ sign, and that shift symmetry requires it, so the printed sign is a typo.
>
> **Table V**: a finite term equivalent to $`72\,\gamma_2`$ acting on a constant is missing. It is most naturally a term $`-12x`$ dropped from the $`\delta`$ row; equivalently, the $`\gamma_2`$ row becomes $`-\tfrac5x+54`$ instead of $`-\tfrac5x-18`$. This term does not depend on $`y`$, so it cancels in the retarded (Schwinger–Keldysh) kernel. It reappears in Task 17 as the constant of the $`(aa')^3(a^2+a'^2)H^4`$ row of Table X ($`-54\to+18`$).
>
> With both corrections, Tables IV + V reproduce the direct sum $`\sum_k P_k`$ at three random points to 30 digits. All other entries stand as printed.
>
> **Sign convention:** the tables carry the paper's sign (that of the printed (I.47)–(I.49)). In the definition-sign convention of (I.168) every entry is multiplied by $`-1`$. This sign difference is inherited from the earlier chapters.
>
> **Point to follow (flagged in Task 13):** if the $`\square`$ in $`\alpha`$ and $`\delta`$ is applied outside the integral as a distribution, the local contribution (I.54)+(I.58) may be counted a second time. This was revisited in Task 16 and is only partially resolved.

## Comparison with the paper

| Item | Our result | Status |
|---|---|---|
| (I.94)–(I.103), (A8)–(A9) | derived; identical | agrees |
| (I.104)–(I.115) | definitions; (I.108), (I.112) correct; weights consistent | agrees |
| (I.116)–(I.123) | exact-$`D`$ leading coefficients and $`D=4`$ forms correct | agrees |
| (I.124)–(I.125) | correct up to a constant annihilated by $`\alpha`$ | agrees |
| Appendix Tables XI–XX | totals exact after 5 corrections: Table XII header ($`(aa')^{D-2}\to(aa')^{D-1}`$), (A18) ($`I[f'A]\to I[f'A']`$), Table XIII row $`P_{7c}`$, Table XV row $`P_{4b}`$ ($`y'\to y^2`$), Table XVIII total ($`A'^2\Delta C'\to A'\Delta C'`$) | typos (C9–C13) |
| Appendix, text before (A31) | the paper writes "Applying the second identity (A18) to $`P_1`$ and $`P_6`$"; it should read "... to $`P_1`$ and $`P_4`$" (the subparts that follow are $`P_{1a},P_{1b},P_{4a},\dots`$) | typo, wording only (C23; noticed while refereeing the write-up, confirmed on the printed page) |
| Table IV | 19 of 20 entries agree; $`\gamma_2`$ second column has the wrong sign | typo (C15; resolved in Task 16) |
| Table V | all non-analytic terms agree; finite term $`72\gamma_2`$ missing | substantive, analytic (S3); drops out of the retarded kernel |
| Overall sign | derived $`=-`$ paper | inherited (S1) |

A note on scope: the appendix was audited, not redone. Worked example rows were checked by hand, and the table totals were verified by the global pointwise test (numerical, at sample points; it checks only sums). It was not re-derived row by row. Where the $`72\gamma_2`$ term belongs (a $`-12x`$ dropped from the $`\delta`$ row, or a changed $`\gamma_2`$ constant) is an interpretation: the test fixes only the sum.

## Code

All scripts are in [`../code/06_tables_IV_V/`](../code/06_tables_IV_V/). They need Python 3.10+ and sympy. Run each one from that folder, for example `cd code/06_tables_IV_V && python3 check_table4.py`.

| Script | What it checks | Runtime |
|---|---|---|
| [`check_94_125.py`](../code/06_tables_IV_V/check_94_125.py) | (I.94)–(I.103) pointwise at $`D=4,5`$; (I.99), (I.108), (I.112); leading coefficients and $`D=4`$ forms (I.116)–(I.123); (I.124)–(I.125) | ~1.5 min |
| [`check_global.py`](../code/06_tables_IV_V/check_global.py) | appendix totals under the ten operators equal $`\sum_k P_k`$ (generic powers, $`D=4,5`$) | ~15 s |
| [`check_tables45.py`](../code/06_tables_IV_V/check_tables45.py) | same with the actual $`D=4`$ propagators; per-operator differences from the printed Tables IV + V | ~20 s |
| [`check_tables45b.py`](../code/06_tables_IV_V/check_tables45b.py) | hypothesis tests: which corrections make Tables IV + V consistent | ~30 s |
| [`check_table4.py`](../code/06_tables_IV_V/check_table4.py) | Table IV in general $`D`$ (19/20, $`\gamma_2`$ sign) | ~2 s |
| [`diag_leading_pair.py`](../code/06_tables_IV_V/diag_leading_pair.py) | diagnostic run for the $`A'\Delta C'`$ coefficient | ~2 s |
| [`superseded_test_E.py`](../code/06_tables_IV_V/superseded_test_E.py) | the first global test, later replaced (kept for the record) | ~15 s |

Expected outputs are in [`expected_output/`](../code/06_tables_IV_V/expected_output/). They match the outputs recorded during the original calculation line by line.

## Full derivation

The full step-by-step derivation, with every intermediate result, is in the chapter **"Tables IV and V: Eqs. (I.94)-(I.125)"** of the notes: [`../notes/notes.pdf`](../notes/notes.pdf).
