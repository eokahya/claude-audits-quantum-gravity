# 07 — de Sitter symmetries, the gauge, and noninvariant counterterms (Paper I, Eqs. (I.126)–(I.145))

Code checks for Tasks 14 and 15 of the re-derivation of
E. O. Kahya and R. P. Woodard, *Quantum gravity corrections to the one loop scalar self-mass during inflation*,
Phys. Rev. D **76**, 124005 (2007) ("Paper I").

These scripts back up the analytic derivation in the notes
(chapter "de Sitter symmetries, the gauge, and noninvariant counterterms: Eqs. (I.126)–(I.145)").
They do not replace it: every result was first derived analytically, and the scripts test the derived statements independently.

## Requirements

Python 3.10+ and `sympy` (tested with Python 3.13.1 and sympy 1.14.0). No network access is needed.

## Scripts

| Script | Task | What it verifies | Runtime |
|---|---|---|---|
| `check_126_135.py` | 14 | (A) Christoffel symbols and $`R_{\mu\nu}=(D-1)H^2\hat g_{\mu\nu}`$ for $`D=4,5`$; (B) the commutator identity used to reduce $`\phi_{;\mu\nu}\phi^{;\mu\nu}`$; (C) the finite maps (I.132)–(I.135) are isometries and time translation is not; (D) the general Killing vector solves the Killing equation for $`D=3,4,5`$ with $`D(D+1)/2`$ parameters, and a polynomial ansatz of degree $`\le3`$ gives a solution space of dimension exactly 6 ($`D=3`$) and 10 ($`D=4`$); (E) infinitesimal and $`\eta\to0`$ limits of (I.135); (F) $`F_\mu=G_\mu+E_\mu`$ (covariant de Donder vector plus a piece built on the unit normal $`u^\mu`$); (G) $`[\xi,u]=(0,2H\eta^2\vec b)`$; (H) the gauge-fixing Lagrangian is invariant (up to a total derivative) under translations, rotations and dilatation, and not under spatial special conformal transformations, whose anomaly is not a boundary term; (I) the Lie algebra of the Killing vectors. | ~10 s |
| `check_136_145.py` | 15 | (I.136), (I.137) and the partial integrations (I.138)–(I.143), each discarded current checked exactly (symbolic $`D`$); completeness: the Euler–Lagrange operators of the basis $`B_1,\dots,B_5`$ are independent and those of 15 admissible candidate densities lie in their span ($`D=4,5,6`$); the vertex (I.145); the remark after (I.145) that $`c(\eta)\,\partial_0\varphi\,\nabla^2\varphi`$ has no $`\partial_0`$ in its vertex. | ~5 s |

## How to run

```bash
cd code/07_desitter_counterterms
python3 check_126_135.py
python3 check_136_145.py
```

## Expected output

`expected_output/check_126_135.txt` and `expected_output/check_136_145.txt` hold the output of the runs made for these notes.
Both are identical, line by line, to the outputs recorded when the original Tasks 14 and 15 were carried out.

How to read the output:

- In `check_126_135.py` every check line prints `True` (or the expected expression), except the contrast test
  `for contrast, time translation t -> t + c : False` (time translation is not an isometry in conformal coordinates).
  Section (A) also prints $`R=12H^2`$ and $`20H^2`$; in (D), `fail: [True, True, True]` means that time translation, the boost and the temporal special conformal transformation are indeed not Killing vectors;
  section (G) prints the commutator as the list `[0, 2*H*b1*t**2, 2*H*b2*t**2, 2*H*b3*t**2]` (`t` is the conformal time $`\eta`$);
  section (H) prints `delta L_GF - d_l(xi^l L_GF) = 0` for translation, rotation and dilatation and, for the spatial special conformal transformation, that the anomaly is nonzero, that it equals
  $`-a^{D-2}\eta^{\mu\nu}F_\mu\Delta_\nu`$ with $`\Delta_\mu=-(4h_{\mu i}b^i-2b_\mu h)`$, and that its Euler–Lagrange derivatives with respect to $`h_{00}`$ and $`h_{12}`$ are nonzero.
- In `check_136_145.py` each candidate density is followed by its coordinate vector `c(B1..B5)`,
  meaning that the density is equal, up to a total derivative, to $`c_1B_1+\dots+c_5B_5`$, where
  $`B_1=a^D(\square\varphi)^2`$, $`B_2=H^2a^{D-2}\eta^{\mu\nu}\partial_\mu\varphi\partial_\nu\varphi`$,
  $`B_3=a^{D-2}\square\varphi\,\nabla^2\varphi`$, $`B_4=a^{D-4}(\nabla^2\varphi)^2`$, $`B_5=H^2a^{D-2}\nabla\varphi\cdot\nabla\varphi`$.
  No line may read `NOT in span`.

## Scope notes

- In `check_136_145.py` the field depends on conformal time and three spatial coordinates; the spacetime dimension $`D`$ enters through the powers of $`a`$ and through $`\square=a^{-2}[-\partial_0^2-(D-2)Ha\,\partial_0+\nabla^2]`$, which are the only places where $`D`$ appears in these reductions. The identities (I.136)–(I.143) are checked with symbolic $`D`$; the completeness test is run at $`D=4,5,6`$ only.
- In `check_126_135.py`, test (A) is done in $`D=4,5`$, test (D) in $`D=3,4,5`$ (polynomial ansatz in $`D=3,4`$), test (H) in $`D=3`$ with a generic $`h_{\mu\nu}(\eta,x_1,x_2)`$, and tests (B), (C), (E), (F), (G), (I) in $`D=4`$. The general-$`D`$ statements are proved analytically in the notes.
