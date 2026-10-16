# 7. de Sitter symmetries, the gauge, and noninvariant counterterms: Eqs. (I.126)–(I.145)

<sub>[Overview (README)](../README.md) · Previous: [Tables IV and V](06-tables-IV-V.md) · Next: [Renormalized self-mass and Table X](08-renormalized-self-mass-table-X.md)</sub>

**Paper I, Eqs. (I.126)–(I.145)** · Tasks 14 and 15 · [Full derivation (notes PDF)](../notes/notes.pdf) · [Code](../code/07_desitter_counterterms/)

## What this step computes

Before the divergences of the self-mass can be subtracted, one has to know which counterterms are allowed.
If the gauge were de Sitter invariant, two invariant counterterms would be enough.
The graviton gauge of Paper I, Eq. (I.15), is not invariant, so we first work out exactly which part of the de Sitter group it keeps.
To do this we solve the Killing equation in conformal coordinates exactly.
The solution splits the $`\tfrac12D(D+1)`$ isometries into four classes:
spatial translations, rotations, dilatation, and spatial special conformal transformations.
We then prove that the gauge-fixing term is invariant under the first three classes and not under the fourth.
The three surviving symmetries restrict the possible noninvariant counterterms to three structures, Eq. (I.144).
In the flat-space limit the gauge becomes Poincaré invariant, so the noninvariant counterterms must vanish there; this keeps only one of them, $`H^2a^{D-2}\nabla\varphi\cdot\nabla\varphi`$, whose vertex is Eq. (I.145).

## Key steps

**1. Invariant counterterms (I.126)–(I.131).**
Dimension counting at one loop (in $`D=4`$) leaves four mass units besides $`\kappa^2\phi\phi`$.
The shift symmetry $`\phi\to\phi+\text{const}`$, masslessness and the background equation $`\hat R_{\mu\nu}=(D-1)H^2\hat g_{\mu\nu}`$ remove everything except

```math
\tfrac12\alpha_1\kappa^2\,\square\phi\,\square\phi\,a^D
\qquad\text{and}\qquad
-\tfrac12\alpha_2\kappa^2H^2\,\partial_\mu\phi\,\partial^\mu\phi\,a^{D-2}.
```

Their vertices, $`i\alpha_1\kappa^2a^D\square^2\delta^D(x-x')`$ and $`i\alpha_2\kappa^2H^2a^D\square\,\delta^D(x-x')`$, match the paper.

**2. All Killing vectors of de Sitter in conformal coordinates.**
With $`\hat g_{\mu\nu}=a^2\eta_{\mu\nu}`$, $`a=-1/(H\eta)`$, the Killing equation becomes $`\partial_\mu\xi_\nu+\partial_\nu\xi_\mu=(2\xi^0/\eta)\,\eta_{\mu\nu}`$, with $`\xi_\mu\equiv\eta_{\mu\nu}\xi^\nu`$ (flat lowering). Its general solution is

```math
\xi^0=\eta\,(\lambda+2\vec b\cdot\vec x),\qquad
\xi^i=c^i+\omega^i{}_j x^j+\lambda x^i+2(\vec b\cdot\vec x)\,x^i-(x\cdot x)\,b^i ,
```

with $`x\cdot x=\Vert\vec x\Vert^2-\eta^2`$. The parameters $`c`$, $`\omega`$, $`\lambda`$, $`\vec b`$ are the four classes, and

```math
(D-1)+\tfrac12(D-1)(D-2)+1+(D-1)=\tfrac12D(D+1).
```

Their finite forms are exactly Eqs. (I.132)–(I.135).

**3. Where the gauge breaks the symmetry.**
The gauge-fixing term can be written covariantly as $`\mathcal L_{\rm GF}=-\tfrac12\sqrt{-\hat g}\,\hat g^{\mu\nu}F_\mu F_\nu`$, and its vector splits as

```math
F_\mu=G_\mu+E_\mu,\qquad
E_\mu=H\big(2\psi_{\mu\nu}u^\nu-u_\mu\psi\big),\qquad
u^\mu=\delta^\mu_0/a,\quad \psi_{\mu\nu}=a^2h_{\mu\nu}.
```

Here $`G_\mu`$ is the covariant de Donder vector of the background, $`\psi=\hat g^{\rho\sigma}\psi_{\rho\sigma}=h`$, and $`u^\mu`$ is the unit normal of the equal-time slices.
So the gauge term is invariant exactly under the isometries that preserve $`u`$.

**4. Only the spatial special conformal transformations move the slicing.**

```math
[\xi,u]=\big(0,\;2H\eta^2\,\vec b\big),
\qquad
\delta F_\mu=\mathcal{L}_\xi F_\mu+\Delta_\mu,\quad
\Delta_\mu=-\big(4h_{\mu i}b^i-2\delta^i_\mu b_i h\big).
```

Here $`\mathcal{L}_\xi`$ is the Lie derivative and $`\delta`$ the induced change of $`F_\mu`$ when $`\psi_{\mu\nu}=a^2h_{\mu\nu}`$ transforms as a tensor.
The breaking is real: a spatial special conformal transformation maps a configuration with $`F_\mu=0`$ to one with $`\delta F_\mu=\Delta_\mu\neq0`$ (for generic $`h_{\mu\nu}`$), and the resulting change of the gauge-fixing action is not a boundary term (its Euler–Lagrange derivatives with respect to $`h_{00}`$ and $`h_{12}`$ are nonzero; this part was checked symbolically in $`D=3`$).

**5. The allowed noninvariant counterterms (I.136)–(I.144).**
Following the paper, the scalar is written $`\varphi`$ from here on.
Translations, rotations and dilatation force every counterterm into the form $`\kappa^2a^D\,Q(a^{-1}\partial_0,\,a^{-1}\partial_i,\,H)[\varphi,\varphi]`$.
Using $`a^{-2}\partial_0^2=-\square-(D-2)Ha^{-1}\partial_0+a^{-2}\nabla^2`$ and partial integrations such as

```math
Ha^{D-1}\partial_0\varphi\,\square\varphi\;\to\;\tfrac12(D-1)H^2\,\partial_\mu\varphi\,\partial_\nu\varphi\,g^{\mu\nu}\sqrt{-g}-H^2a^{D-2}\nabla\varphi\cdot\nabla\varphi ,
```

every candidate reduces, up to total derivatives, to five independent terms: the two invariant ones plus

```math
\tfrac12\kappa^2a^{D-2}\,\square\varphi\,\nabla^2\varphi,\qquad
\tfrac12\kappa^2a^{D-4}\,\nabla^2\varphi\,\nabla^2\varphi,\qquad
-\tfrac12\kappa^2H^2a^{D-2}\,\nabla\varphi\cdot\nabla\varphi .
```

Extent of the check: general-$`D`$ reduction formulas were derived for eleven candidate densities; the identities (I.136)–(I.143) were verified with symbolic $`D`$; the reduction of all fifteen candidates, and the independence of the five terms, were verified by code at $`D=4,5,6`$ only.

**6. Flat-space limit and the vertex (I.145).**
As $`H\to0`$ at fixed comoving time the gauge becomes Poincaré invariant.
Assuming that the ultraviolet divergences depend polynomially on $`H`$, the coefficients of the noninvariant counterterms must therefore vanish as $`H\to0`$.
The first two terms carry four derivatives, so their coefficients are pure numbers that cannot contain $`H`$; they must be zero.
The third term is proportional to $`H^2`$ and survives. Its vertex is

```math
-\tfrac12\alpha_3\kappa^2H^2a^{D}\,\frac{\nabla}{a}\varphi\cdot\frac{\nabla}{a}\varphi
\;\longrightarrow\;
i\alpha_3\kappa^2H^2a^{D-2}\,\nabla^2\delta^D(x-x').
```

## Result

> **Two invariant counterterms and one noninvariant counterterm suffice:**
> $`\tfrac12\alpha_1\kappa^2(\square\phi)^2a^D`$, $`-\tfrac12\alpha_2\kappa^2H^2\partial_\mu\phi\,\partial^\mu\phi\,a^{D-2}`$, and $`-\tfrac12\alpha_3\kappa^2H^2a^{D-2}\nabla\varphi\cdot\nabla\varphi`$.
> The gauge (I.15) preserves spatial translations, rotations and dilatation, and breaks only the spatial special conformal transformations (I.135), because these bend the equal-time slices.
> All of Eqs. (I.126)–(I.145) are reproduced; **no discrepancy with the paper** affects any result. There is one omission in the list (I.128) (typo C14 in the [complete list of findings](10-complete-list-of-findings.md); no effect) and two notational remarks, all listed in the table below.
>
> Assumptions used (implicit in the paper's argument): dimensional regularization preserves the symmetries of the gauge-fixed action, and the ultraviolet divergences are polynomial in $`H`$.
> The vertex rule is "vertex $`=i\,\delta^2S/\delta\phi\,\delta\phi`$" (the definition sign of Eq. (I.168)). The overall sign of the printed self-mass (see the earlier chapters) changes only the signs of the fitted $`\alpha_i`$, not the structure of the counterterms.

## Comparison with the paper

| Item | Paper I | Our result | Status |
|---|---|---|---|
| Invariants (I.126)–(I.127) | listed | derived by dimension counting | agrees |
| Four-derivative list (I.128) | 5 terms | also $`\phi^2R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}`$ and $`\phi^2\square R`$ belong to the list; both are removed by the shift symmetry | omission (C14), no effect on (I.129) |
| (I.129)–(I.131) | two counterterms and their vertices | derived | agrees |
| Four classes, $`\tfrac12D(D+1)`$ | stated | exact solution of the Killing equation | agrees |
| Dilatation (I.134) | $`\eta'=k\eta`$ | needs $`k>0`$ (otherwise $`\eta<0`$ is mapped to $`\eta>0`$) | remark |
| Special conformal (I.135) | stated | isometry proved; $`x\cdot x`$ must be read as the Minkowski product $`\Vert\vec x\Vert^2-\eta^2`$ | agrees (remark) |
| "Gauge breaks only (I.135)" | stated with a reference | proved via $`F=G+E`$ and $`[\xi,u]`$ | agrees |
| (I.136)–(I.143) | reductions | derived, every discarded total derivative explicit | agrees |
| (I.144) "only three" | stated | completeness: 5-dimensional space modulo total derivatives (code check at $`D=4,5,6`$) | agrees |
| Flat-space argument | stated | derived under the stated assumption | agrees |
| Vertex (I.145) | printed | derived | agrees |

## Code

All scripts are in [`code/07_desitter_counterterms/`](../code/07_desitter_counterterms/) (Python 3.10+, sympy). Run them from that folder:

```bash
cd code/07_desitter_counterterms
python3 check_126_135.py   # about 10 s
python3 check_136_145.py   # about 5 s
```

- [`check_126_135.py`](../code/07_desitter_counterterms/check_126_135.py) — background curvature; the commutator identity; isometry of (I.132)–(I.135); the general Killing vector and its parameter count ($`D=3,4,5`$); $`F=G+E`$; $`[\xi,u]`$; invariance of the gauge-fixing term under translations, rotations and dilatation, and a genuine (non-boundary) anomaly under spatial special conformal transformations ($`D=3`$); the Lie algebra.
- [`check_136_145.py`](../code/07_desitter_counterterms/check_136_145.py) — (I.136)–(I.143) with every discarded current (symbolic $`D`$); completeness of the basis of five terms via Euler–Lagrange operators (15 candidate densities, at $`D=4,5,6`$); the vertex (I.145); the remark that no $`Ha^{D-3}\nabla^2\partial_0\delta^D`$ term can arise.

Expected outputs are in [`expected_output/`](../code/07_desitter_counterterms/expected_output/); they are identical to the outputs recorded in the original tasks.

## Full derivation

The complete step-by-step derivation, with index tracking, all intermediate results and the code-check records, is in the notes:
[notes.pdf](../notes/notes.pdf), chapter **"de Sitter symmetries, the gauge, and noninvariant counterterms: Eqs. (I.126)–(I.145)"**.
