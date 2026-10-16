# 02 — The 4-point contribution, Eqs. (I.31)–(I.44)

These scripts are supplementary computer-algebra checks for Tasks 3, 4 and 5 of the
re-derivation of Paper I (E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007)).
They do **not** replace the analytic derivation, which is in the notes
(`../../notes/notes.pdf`, chapter "The 4-point contribution: Eqs. (I.31)-(I.44)").

All three scripts are reconstructed verbatim from the original command logs of the
tasks, including the small edits that were already made within the logs (one in
`check_4pt.py`: the eigenvalue print statement; two clean-up edits in `check_41_44.py`;
none in `check_coincidence.py`). Each edit is noted in the script's docstring. The only
additions are the module docstrings. The scripts import only `sympy`, `mpmath` and the
standard library, use no file paths, and write nothing to disk.

## Requirements

Python 3.10+ with `sympy` and `mpmath`. No network access. Tested with Python 3.13.1,
SymPy 1.14.0, mpmath 1.3.0.

## How to run

```bash
cd code/02_four_point
python3 check_4pt.py          # Task 3, ~25 s
python3 check_coincidence.py  # Task 4, ~2 s
python3 check_41_44.py        # Task 5, ~1 s
```

The outputs of these runs are stored in `expected_output/<script>.txt`. They are
identical (byte for byte) to the outputs recorded in the original task logs.

## Scripts

| Script | Task | Equations | What it checks | Runtime |
|---|---|---|---|---|
| `check_4pt.py` | 3 | (I.31)–(I.33), (I.131) | (a) Lemma A for the second functional derivative of a quadratic form with derivatives, with explicit test functions; (b) lattice Hessian and positivity; (c) operator-level check of (I.33) in D=3 with a generic pair-symmetric coincident propagator; (d) calibration with the paper's counterterm vertex (I.131) | ~25 s |
| `check_coincidence.py` | 4 | (I.34)–(I.40) | (1) Gauss's theorem on the hypergeometric forms (I.29)–(I.30); (2) numerical coincidence limits at D=0.6; (3) explicit-component contractions of the graviton tensor factors for D=4,5,6,7 and trace identities; (4) D=4 finite parts compared with the printed ones | ~2 s |
| `check_41_44.py` | 5 | (I.41)–(I.44) | (A) symbolic-D bookkeeping of (I.33) in the basis P[g], S[g]; (B) explicit D=4 kernel on a test function; (C) Eq. (I.43) as a kernel identity, D=4 and general D; (D) (I.41) → (I.44) on a test function | ~1 s |

## What the output means

### `check_4pt.py`

* `(a)`: `A = 97*sqrt(6)*pi*exp(-1/2)/54 ~ 8.384` is the bilinear term of the quadratic
  form; the kernel `-∂_μ[(K+Kᵀ)^{μν} ∂_ν δ]` gives the same value (`True`), the
  opposite-sign kernel gives `-A` (`True`).
* `(b)`: the Hessian equals the lattice version of `-∂(2c ∂ δ)` (`True`), not of
  `+∂(2c ∂ δ)` (`False`); the eigenvalues for `c=1` are all `>= 0`.
* `(c)`: `derived kernel == derived (33)` and `== -(paper's 33)` are `True`, so the
  derived (I.33′) is exactly −1 times the printed (I.33). The paper's reading rule
  combined with the derived `U2 = +i/4` reproduces the printed (I.33) (`True`); with the
  printed `U2 = -i/4` it does not (`False`). So the Table II sign of `U2` is a typo
  and the overall sign comes from the paper's rule for `∂_2`.
* `(d)`: Lemma A reproduces the paper's own counterterm vertex (I.131) (`True`); the
  paper's 4-point rule applied to the same vertex does not (`False`).

### `check_coincidence.py`

* `(1)`: both `True`: (I.35) and (I.36) follow from the hypergeometric forms.
* `(2)`: series and hypergeometric forms agree to ~1e-31 at y=0.5; the differences
  from (I.34)–(I.36) go to zero like (y/4)^0.7.
* `(3)`: all `True` for D = 4, 5, 6, 7.
* `(4)`: finite parts at D=4 (h = H²/16π²): (37) `4h`, (38) `-2h`, (39) `[2,2,0]h`, (40) `[-2,0,0]h`.
  The paper prints `+2h` for the `(α,σ)=(0,0)` entry of (I.40); its trace would then be
  `+2h`, inconsistent with the paper's own `-2h` in (I.38).

### `check_41_44.py`

* `(A)`: the iΔ_A-part coefficients `I*D*(D-1)/4` (of P) and `-I*D*(D-5)/4` (of S)
  give `-iκ²{¼D(D-1)∂^μ(a^{D-2}Λ∂_μ) - D a^{D-2}Λ∇²}`; the finite B/C parts are
  proportional to `(D-4)`. The paper's chain ((I.33) as printed plus the printed (I.40))
  gives the finite part `4*I*h` times `P[a²]`, which is the extra bracket
  `{-∂^μ(a²∂_μδ⁴)+a²∇²δ⁴}` of the printed (I.41); with the derived (I.40) it is `0`.
* `(B)`: `True`, `True`.
* `(C)`: (I.43) holds in D=4 and in general D (`True`), fails with `ln(a²)` in place
  of `ln(aa')` (`False`); `Box ln a = -(D-1)H²` is `True`.
* `(D)`: derived (I.41′) equals derived (I.44′) and printed (I.41) equals printed
  (I.44) (both `True`): the algebra (I.41) → (I.44) in the paper is correct.
