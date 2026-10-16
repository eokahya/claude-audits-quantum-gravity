"""check_eq6.py -- supplementary computer-algebra check of Eq. (I.6) (Task 1).

Purpose
    Verify the expansion of the scalar kinetic-term density sqrt(-g) g^{mu nu} to
    second order in kappa, i.e. the curly bracket of Eq. (I.6) of
    E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007):

        a^{-(D-2)} sqrt(-g) g^{mu nu} = eta^{mu nu} - k h^{mu nu} + (1/2) k eta^{mu nu} h
            + (1/8) k^2 eta^{mu nu} h^2 - (1/4) k^2 eta^{mu nu} h^{rs} h_{rs}
            - (1/2) k^2 h h^{mu nu} + k^2 h^{mu rho} h_rho^nu + O(k^3),

    with g_{mu nu} = a^2 (eta_{mu nu} + k h_{mu nu}) [Eq. (I.3)] and all indices of
    h moved with eta.  As a by-product it also checks Eqs. (I.4) and (I.5)
    implicitly (the exact inverse and the exact determinant are used).

Method
    The exact inverse metric (adjugate / determinant) and the exact sqrt(-det g)
    are computed by sympy -- no Neumann series and no "det = exp Tr ln" input --
    then Taylor-expanded in kappa to O(kappa^2) and compared entry by entry with
    the derived bracket.  The a-dependence is exactly the common factor a^{D-2},
    so a = 1 is used.
      (1) D = 2, 3: fully symbolic general symmetric h_{mu nu} (3 and 6 components).
      (2) D = 4, 5, 6: random rational symmetric h_{mu nu}, exact arithmetic
          (random.seed(1)).
      (3) conformal special case h_{mu nu} = c eta_{mu nu} for symbolic D:
          (1 + k c)^{(D-2)/2} versus the derived bracket evaluated on h = c eta.
      (4) D = 4, symbolic: the truncated inverse satisfies
          g_{mu rho} g^{rho nu} - delta_mu^nu = k^3 h_{mu rho} h^rho_sigma h^{sigma nu}.

Expected result
    Every line prints True; the kappa^2 coefficient in the conformal case is
    c**2*(D - 4)*(D - 2)/8.

Runtime
    About 20-25 s (sympy 1.14, Apple silicon).

Run
    python3 check_eq6.py
"""
import random

import sympy as sp


def eta(D):
    return sp.diag(*([-1] + [1]*(D-1)))


def formula(hl, D, k):
    """Derived bracket: a^{-(D-2)} sqrt(-g) g^{mu nu} to O(k^2), indices raised with eta."""
    E = eta(D); Ei = E.inv()
    hu = Ei*hl*Ei                  # h^{mu nu}
    hm = Ei*hl                     # h^mu_nu  (row up, column down)
    tr = hm.trace()                # h
    hh = sum(hl[i, j]*hu[i, j] for i in range(D) for j in range(D))   # h^{rs} h_{rs}
    hmrhrn = hu*E*hu               # h^{mu rho} eta_{rho sigma} h^{sigma nu} = h^{mu rho} h_rho^nu
    return Ei + k*(-hu + sp.Rational(1, 2)*Ei*tr) + k**2*(hmrhrn - sp.Rational(1, 2)*tr*hu
            + sp.Rational(1, 8)*Ei*tr**2 - sp.Rational(1, 4)*Ei*hh)


def exact_series(hl, D, k):
    """Exact a^{-(D-2)} sqrt(-g) g^{mu nu} with g = eta + k h (a-factors stripped), series to k^2."""
    g = eta(D) + k*hl
    detg = g.det()
    ginv = g.adjugate()/detg
    M = sp.sqrt(-detg)*ginv
    return M.applyfunc(lambda e: sp.series(e, k, 0, 3).removeO())


k = sp.symbols('kappa')

# (1) fully symbolic, general symmetric h_{mu nu}, D = 2 and D = 3
for D in (2, 3):
    hl = sp.zeros(D, D)
    for i in range(D):
        for j in range(i, D):
            s = sp.Symbol(f'h{i}{j}'); hl[i, j] = s; hl[j, i] = s
    diff = (exact_series(hl, D, k) - formula(hl, D, k)).applyfunc(sp.expand)
    print(f"symbolic D={D}: all entries zero ->", all(e == 0 for e in diff))

# (2) random rational symmetric h_{mu nu}, exact arithmetic, D = 4,5,6
random.seed(1)
for D in (4, 5, 6):
    hl = sp.zeros(D, D)
    for i in range(D):
        for j in range(i, D):
            v = sp.Rational(random.randint(-9, 9), random.randint(1, 7)); hl[i, j] = v; hl[j, i] = v
    diff = (exact_series(hl, D, k) - formula(hl, D, k)).applyfunc(sp.expand)
    print(f"random rational D={D}: all entries zero ->", all(e == 0 for e in diff))

# (3) conformal special case h = c eta, symbolic D not possible as matrix size; check closed form in D
c, Dsym = sp.symbols('c D')
exact = sp.series((1 + k*c)**(Dsym/2 - 1), k, 0, 3).removeO()
gen = 1 + k*(-c + sp.Rational(1, 2)*c*Dsym) + k**2*(c**2 - sp.Rational(1, 2)*(c*Dsym)*c
        + sp.Rational(1, 8)*(c*Dsym)**2 - sp.Rational(1, 4)*c**2*Dsym)
print("conformal case, general D:", sp.simplify(sp.expand(exact - gen)) == 0)
print("  order-2 coefficient of exact:", sp.factor(sp.expand(exact).coeff(k, 2)))

# (4) the separate g g^{-1} truncation check: residual of order kappa^3
D = 4
hl = sp.zeros(D, D)
for i in range(D):
    for j in range(i, D):
        s = sp.Symbol(f'h{i}{j}'); hl[i, j] = s; hl[j, i] = s
E = eta(D); Ei = E.inv(); hu = Ei*hl*Ei
ginv2 = Ei - k*hu + k**2*(Ei*hl*Ei*hl*Ei)    # h^mu_sigma h^{sigma nu}
res = ((E + k*hl)*ginv2 - sp.eye(D)).applyfunc(sp.expand)
expected = (k**3 * hl*Ei*hl*Ei*hl*Ei).applyfunc(sp.expand)   # h_{mu rho} h^rho_sigma h^{sigma nu}
print("g_{mu rho} g^{rho nu} - delta == kappa^3 h_{mu rho} h^rho_sigma h^{sigma nu} (D=4 symbolic):",
      all(sp.expand(r - e) == 0 for r, e in zip(res, expected)))
