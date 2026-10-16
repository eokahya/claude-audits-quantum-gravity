# 03 – The 3-point contribution: setup and Eqs. (I.45)–(I.53)

Code checks for Chapter 3 of the derivation notes (Tasks 6, 7 and 8).
Paper I = E. O. Kahya and R. P. Woodard, *Quantum gravity corrections to the one loop scalar self-mass during inflation*,
Phys. Rev. D **76**, 124005 (2007).

These scripts are **supplementary checks** of the analytic derivations in the notes
(`../../notes/notes.pdf`, chapter "The 3-point contribution: setup and Eqs. (I.45)-(I.53)").
They do not replace the derivations.

## Requirements and how to run

Python 3.10+ with `sympy` and `mpmath` (see `../requirements.txt`). No network access is needed.
Each script is standalone; run it from this folder:

```bash
cd code/03_three_point_setup   # from the repository root
python3 check_47_49.py     # ~2 s
python3 check_52.py        # ~1-3 s
python3 check_45_46.py     # ~3 min
```

The output of each script is stored in `expected_output/<script>.txt`.

## Scripts

| Script | Task | What it verifies | Runtime |
|---|---|---|---|
| `check_45_46.py` | 6 | (a) index algebra (I.45) → (I.46) with a generic within-pair-symmetric tensor in D = 3 and formal derivative vectors; the printed sign of the left-hand side of (I.7) would flip the 2nd and 3rd terms; (b) overall normalization = 1; (c) the derivative rule of (I.45) (both outer legs → −∂, on everything) against a Gateaux second derivative in a 1D toy model with explicit smooth kernels; a rule with one outer leg of opposite sign (the type used in (I.32)) gives the opposite value | ~170 s |
| `check_47_49.py` | 7 | for D = 4, 5, 6: (I.46) route = (I.45) route; result = derived nine terms (I.47′)–(I.49′) + an antisymmetric AA pair that vanishes by spatial homogeneity of iΔ_A; result = −(printed (I.47)–(I.49)); homogeneity of A(y) + k ln(aa′) | ~1.5 s |
| `check_52.py` | 8 | (a) ∂₀∂′₀ ln(aa′) = 0; (b) kink lemma ∂₀∂′₀(\|Δη\| − iδ)² = −2 + 4iδ δ(Δη); (c) flat-space check of the regular part; (d) nascent delta δ(r²+δ²)^(−D/2) → π^(D/2)/Γ(D/2) δ^(D−1) (closed form, general D, numerics); (e) local coefficient i a^(2−D) with c₀ = Γ(D/2−1) and 2i a^(2−D)/(D−2) with the printed c₀ of (I.14); (f) the printed second coefficient of (I.14) forces c₀ = Γ(D/2−1); (g) □f(y) = H²[(4y−y²)f″ + D(2−y)f′] for D = 4, 5; (h) six derivatives of y used later | ~1.3 s |

## What the output means

* `check_45_46.txt`: `True` in (a) means the four terms of (I.46) follow from (I.45) with the derived Table I;
  `False ... (difference = -2*(t2+t3): True)` shows that the printed left-hand side of (I.7) (−½κ) would reverse
  the 2nd and 3rd terms, so (I.46) was built with the correct vertex V₁ = +iκa^(D−2)∂₁^α∂₂^β.
  In (c) the value 0.779492709611 appears twice (Gateaux derivative and the (I.45) rule agree) and −0.779492709611 for the wrong rule.
* `check_47_49.txt`: for each D the last line, `(46)-route == -(paper's (47)-(49)) ... : True`, is the finding that all nine printed terms
  of (I.47)–(I.49) carry the opposite overall sign to what follows from the paper's own (I.46).
* `check_52.txt`: all identities `True`; the numerical nascent-delta integrals 9.6498, 9.8474, 9.8674 approach π² = 9.8696;
  (e) and (f) establish the typographical error in the first coefficient of (I.14) (printed Γ(D/2−1)/(D/2−1), correct Γ(D/2−1)).

## Provenance

The three scripts are the ones actually run in Tasks 6–8 (15 October 2026); in each, the original one-line header comment
was replaced by a docstring and the code is unchanged. `check_52.py` is the final version of Task 8 (an earlier run printed
check (c) in an unsimplified form).
Re-running them (Python 3.13.1, SymPy 1.14.0, mpmath 1.3.0, macOS) reproduces the outputs recorded in the original task logs
character for character.
