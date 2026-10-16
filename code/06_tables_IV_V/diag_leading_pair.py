"""diag_leading_pair.py -- diagnostic for the (4/y)^(D-1) coefficient of A' Delta C'  [Task 13, debugging run]

Purpose
  A first version of part (C) of check_94_125.py expanded the full series (I.14), (I.27), (I.28)
  and extracted "exact-power" coefficients by matching exponents; it printed wrong-looking
  coefficients for the Delta-type products (e.g. a quartic rational function for A' Delta C').
  This diagnostic (run inline during Task 13) computes
    1. the coefficient of (4/y)^(D-1) in A' Delta C' from the two leading terms only
       (A_lead + A_2, C_lead + C_2), in units Gamma^2(D/2) N^2/16;
    2. every individual term of exponent D-1 in A'(C-A)' from the full series (two orders).
  The six terms of item 2 sum to -4(D-3)/(D-2); times 2(D-2)/(D-3) this is -8, the value of
  item 1 and of Eq. (I.122). Part (C) was then rewritten to use the two leading terms only.

Expected result
  "analytic leading-pair coefficient A'dC' at (4/y)^(D-1): -8" followed by six terms.
Runtime: a few seconds.
Run:     python3 diag_leading_pair.py
"""
import sympy as sp
D, Ns = sp.symbols('D N'); z = sp.Symbol('z', positive=True)
Gm = sp.gamma
cf = Ns*Gm(D/2-1)*z**(D/2-1)
A2 = Ns*Gm(D/2+1)/(D/2-2)*z**(D/2-2)
C2 = -Ns*(3-D/2)*Gm(D/2-1)*z**(D/2-2)
dy = lambda e: sp.diff(e, z)*(-z**2/4)
unit = Ns**2*Gm(D/2)**2/16
val = dy(cf)*(dy(C2) - dy(A2))*2*(D-2)/(D-3)/unit/z**(D-1)
print("analytic leading-pair coefficient A'dC' at (4/y)^(D-1):", sp.factor(sp.simplify(sp.gammasimp(sp.powsimp(val, force=True)))))
# full series, inspect all terms with exponent D-1 in A'dC'
def ser(Dv, nmax=3):
    A = Ns*(Gm(Dv/2 - 1)*z**(Dv/2 - 1) + Gm(Dv/2 + 1)/(Dv/2 - 2)*z**(Dv/2 - 2) - sp.pi*sp.cot(sp.pi*Dv/2)*Gm(Dv - 1)/Gm(Dv/2)
            + sum(Gm(n + Dv - 1)/(n*Gm(n + Dv/2))*z**(-n) - Gm(n + Dv/2 + 1)/((n - Dv/2 + 2)*Gm(n + 2))*z**(-(n - Dv/2 + 2)) for n in range(1, nmax)))
    cf = Ns*Gm(Dv/2 - 1)*z**(Dv/2 - 1)
    C = cf + Ns*sum((n + 1)*Gm(n + Dv - 3)/Gm(n + Dv/2)*z**(-n) - (n - Dv/2 + 3)*Gm(n + Dv/2 - 1)/Gm(n + 2)*z**(-(n - Dv/2 + 2)) for n in range(0, nmax))
    return A, C
A_, C_ = ser(D)
P = sp.expand(dy(A_)*dy(C_ - A_))
for tm in sp.Add.make_args(P):
    c, e = tm.as_independent(z)
    ex = sp.simplify(sp.powsimp(sp.expand_power_base(sp.simplify(z*sp.diff(e, z)/e), force=True), force=True))
    if sp.simplify(ex - (D-1)) == 0:
        print("  term:", sp.simplify(c/unit), " exponent", ex, " e=", e)
