# 1. Lagrangian expansion and Feynman rules: Eqs. (I.6)–(I.8), Tables I–II

<sub>[Overview (README)](../README.md) · Next: [4-point contribution](02-four-point.md)</sub>

## What this step computes

Every later step of the one-loop calculation rests on the scalar–graviton
interaction vertices. This step derives them from scratch. We start from the
massless, minimally coupled scalar Lagrangian
$`\mathcal L = -\tfrac12\partial_\mu\phi\,\partial_\nu\phi\,g^{\mu\nu}\sqrt{-g}`$
and from the perturbed de Sitter metric
$`g_{\mu\nu}=a^2(\eta_{\mu\nu}+\kappa h_{\mu\nu})`$.

- **Task 1** expands the kinetic term to second order in $`\kappa`$, which
  gives Eq. (I.6) of Paper I. Along the way it gives the inverse metric (I.4)
  and the volume element (I.5).
- **Task 2** turns the order-$`\kappa`$ and order-$`\kappa^2`$ terms into the
  vertex operators of Paper I: Eq. (I.7) and Table I (three-point), and
  Eq. (I.8) and Table II (four-point).

The dimension $`D`$ is kept general throughout. No gauge or trace condition is
imposed, and the printed results are never used as input.

## Key steps

**Inverse metric and volume element.** Write the metric as
$`g = a^2\,\eta\,(1+\kappa\mathsf h)`$ with $`(\mathsf h)^\mu{}_\nu = h^\mu{}_\nu`$
(indices moved with $`\eta`$). The Neumann series and
$`\det = \exp\operatorname{Tr}\ln`$ then give

```math
g^{\mu\nu}=\frac{1}{a^2}\Big(\eta^{\mu\nu}-\kappa h^{\mu\nu}+\kappa^2 h^{\mu}{}_{\rho}h^{\rho\nu}+O(\kappa^3)\Big),
\qquad
\sqrt{-g}=a^D\Big(1+\tfrac12\kappa h+\tfrac18\kappa^2h^2-\tfrac14\kappa^2h^{\rho\sigma}h_{\rho\sigma}+O(\kappa^3)\Big).
```

We derived each of these a second way: by coefficient matching for the inverse
metric, and by the Leibniz expansion of the determinant.

**The product.** At second order the only cross term is $`-\tfrac12 h\,h^{\mu\nu}`$:

```math
\sqrt{-g}\,g^{\mu\nu}=a^{D-2}\Big\{\eta^{\mu\nu}-\kappa h^{\mu\nu}+\tfrac12\kappa\eta^{\mu\nu}h
+\tfrac18\kappa^2\eta^{\mu\nu}h^2-\tfrac14\kappa^2\eta^{\mu\nu}h^{\rho\sigma}h_{\rho\sigma}
-\tfrac12\kappa^2 h h^{\mu\nu}+\kappa^2 h^{\mu\rho}h_{\rho}{}^{\nu}+O(\kappa^3)\Big\}.
```

**Vertex normalization.** Paper I extracts a factor $`\tfrac12`$ for the two
identical scalars and none for the two gravitons. Each vertex also carries a
factor $`i`$ from $`e^{iS}`$ (the standard convention, consistent with the
explicit factors of $`i`$ in the paper's tables):

```math
\kappa\mathcal L^{(1)}=\frac{1}{2i}\sum_I\big[V_I^{\alpha\beta}\phi_1\phi_2\big]h_{\alpha\beta},
\qquad
\kappa^2\mathcal L^{(2)}=\frac{1}{2i}\sum_I\big[U_I^{\alpha\beta\rho\sigma}\phi_1\phi_2\big]h_{\alpha\beta}h_{\rho\sigma}.
```

This leads to a compact rule: every vertex equals $`-i`$ times the
corresponding coefficient $`b_I`$ in the curly bracket of (I.6), times
$`\kappa^n a^{D-2}`$ and the contracted tensor structure.

**Independent check from the exact metric.** Two special cases fix all six
vertex coefficients without using (I.6) at all:
$`h_{\mu\nu}=c\,\eta_{\mu\nu}`$, where
$`\sqrt{-g}\,g^{\mu\nu}=a^{D-2}(1+\kappa c)^{(D-2)/2}\eta^{\mu\nu}`$ exactly, and
the case where only $`h_{00}=b`$ is non-zero.

## Result

> **Task 1.** Eq. (I.6) is reproduced exactly, together with (I.4) and (I.5):
> structure by structure, coefficient by coefficient and sign by sign, for
> general $`D`$.
>
> **Task 2.** The derived vertex operators are
>
> ```math
> V_1^{\alpha\beta}=i\kappa a^{D-2}\partial_1^\alpha\partial_2^\beta,\qquad
> V_2^{\alpha\beta}=-\tfrac{i}{2}\kappa a^{D-2}\eta^{\alpha\beta}\partial_1\!\cdot\partial_2,
> ```
>
> ```math
> \begin{aligned}
> U_1&=-\tfrac{i}{8}\kappa^2a^{D-2}\eta^{\alpha\beta}\eta^{\rho\sigma}\partial_1\!\cdot\partial_2, &
> U_2&=+\tfrac{i}{4}\kappa^2a^{D-2}\eta^{\alpha\rho}\eta^{\beta\sigma}\partial_1\!\cdot\partial_2,\\
> U_3&=+\tfrac{i}{2}\kappa^2a^{D-2}\eta^{\alpha\beta}\partial_1^\rho\partial_2^\sigma, &
> U_4&=-i\kappa^2a^{D-2}\partial_1^\alpha\eta^{\beta\rho}\partial_2^\sigma .
> \end{aligned}
> ```
>
> Table I and (I.8) agree with the paper, and so do rows 1, 3 and 4 of
> Table II. There are two sign misprints. The left-hand side of (I.7) should
> carry $`+\tfrac12\kappa`$; the vertex $`V_1`$ itself is correct. Table II,
> row 2 should read $`+\tfrac{i}{4}`$. Neither misprint affects the vertices
> the paper goes on to use: $`V_1`$ is printed correctly, and the next step
> (Task 3) shows that the paper's own Eq. (I.33) uses $`+\tfrac{i}{4}`$.

## Comparison with the paper

| Item | Derived here | Printed in Paper I | Status |
|---|---|---|---|
| (I.4) inverse metric | as above | same | agrees |
| (I.5) volume element | as above | same | agrees |
| (I.6) bracket, 7 coefficients | as above | same | agrees |
| (I.6) left-hand side | $`-\tfrac12\partial_\mu\phi\,\partial_\nu\phi\,g^{\mu\nu}\sqrt{-g}`$ | second $`\phi`$ missing | typo (C2; noticed while preparing these notes) |
| (I.1) first term | $`\partial_\mu\phi\,\partial_\nu\phi\,g^{\mu\nu}`$ | $`\partial_\mu\phi\,\partial_\mu\phi\,g^{\mu\nu}`$ | typo (C1) |
| (I.7) left-hand side | $`+\tfrac12\kappa a^{D-2}\partial_\alpha\phi\,\partial_\beta\phi\,h^{\alpha\beta}`$ | $`-\tfrac12\kappa\ldots`$ | typo, sign (C3) |
| (I.7) right-hand side, $`V_1`$ | $`i\kappa a^{D-2}\partial_1^\alpha\partial_2^\beta`$ | same | agrees |
| Table I, $`V_2`$ | $`-\tfrac{i}{2}\kappa a^{D-2}\eta^{\alpha\beta}\partial_1\!\cdot\partial_2`$ | same | agrees |
| (I.8) and Table II, $`U_1`$ | $`-\tfrac{i}{8}\kappa^2a^{D-2}\eta^{\alpha\beta}\eta^{\rho\sigma}\partial_1\!\cdot\partial_2`$ | same | agrees |
| Table II, $`U_2`$ | $`+\tfrac{i}{4}\kappa^2a^{D-2}\eta^{\alpha\rho}\eta^{\beta\sigma}\partial_1\!\cdot\partial_2`$ | $`-\tfrac{i}{4}\ldots`$ | typo, sign (C5); (I.33) uses $`+\tfrac{i}{4}`$ |
| Table II, $`U_3`$, $`U_4`$ | as above | same | agrees |

The IDs C1–C5 refer to the [complete list of findings](10-complete-list-of-findings.md).

With the printed $`U_2`$, the rebuilt Lagrangian would differ from the exact one
by $`-\tfrac14\kappa^2a^{D-2}\,h^{\rho\sigma}h_{\rho\sigma}\,\partial^\mu\phi\,\partial_\mu\phi`$.

## Code

All scripts are in [`../code/01_feynman_rules/`](../code/01_feynman_rules/).
They need Python 3.10+ and sympy.

| Script | What it checks | Runtime |
|---|---|---|
| [`check_eq6.py`](../code/01_feynman_rules/check_eq6.py) | Compares the exact inverse metric and determinant, expanded to $`O(\kappa^2)`$, with the bracket of (I.6). Cases: symbolic $`h`$ for $`D=2,3`$, random rational $`h`$ for $`D=4,5,6`$, the conformal case for symbolic $`D`$, and the residual of the truncated inverse. | ~22 s |
| [`check_vertices.py`](../code/01_feynman_rules/check_vertices.py) | (A) Full vertices from (I.6) by polarization with distinct legs. (B) Derived and printed tables against the exact metric. (C) Solves for the table coefficients from the exact metric alone. | ~10 s |

```bash
cd code/01_feynman_rules
python3 check_eq6.py
python3 check_vertices.py
```

The expected outputs are in `expected_output/`. They match the outputs of the
original working session line for line.

## Full derivation

The complete step-by-step derivation is in the notes PDF
([`../notes/notes.pdf`](../notes/notes.pdf)), chapter "Lagrangian expansion and
Feynman rules: Eqs. (I.6)-(I.8), Tables I-II". It includes all index
manipulations, the two independent routes for (I.4) and (I.5), the analytic
special-case checks, and the normalization argument for the vertices.
