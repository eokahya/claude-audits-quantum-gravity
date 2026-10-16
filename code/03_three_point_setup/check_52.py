"""Task 8 -- checks of Eq. (I.52):  d_rho d'_sigma i Delta_A(x;x')  (local + regular parts), and of the
first coefficient of (I.14).

Paper I = E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

Purpose
    Supplementary code checks of the analytic derivation in Chapter 3, Section "Task 8" of the notes.
    They do not replace the derivation.

What is checked
    (a) d_0 d'_0 ln(a a') = 0 (inner mixed derivatives remove the de Sitter breaking term).
    (b) Kink lemma: d_0 d'_0 (|deta| - i delta)^2 = -2 + 4 i delta DiracDelta(deta).
    (c) Flat-space check: for F = X^(1-D/2), X = r^2 - (|deta| - i delta)^2, on both sides of deta = 0,
        eta^{rs} d_r X d'_s X = -4X and (d_r d'_s X)_reg = -2 eta_rs (explicit, 4 coordinates); with the
        trace eta^{rs}(-2 eta_rs) = -2D the contracted chain-rule part vanishes identically for general D.
    (d) Nascent delta: Omega_{D-2} int_0^oo s^(D-2) (1+s^2)^(-D/2) ds = pi^(D/2)/Gamma(D/2), closed form for
        D = 3,4,5,6 and for general D via the Beta function; numerically (D = 4) for delta = 1e-2, 1e-3, 1e-4.
    (e) Coefficient of the local term: i a^(2-D) with c0 = Gamma(D/2-1) (leading term = i Delta_cf, (I.26));
        2i a^(2-D)/(D-2) with the c0 printed in (I.14).
    (f) ODE recursion: the printed second coefficient c1 = Gamma(D/2+1)/(D/2-2) of (I.14) requires
        c0 = Gamma(D/2-1), not the printed Gamma(D/2-1)/(D/2-1).
    (g) Box f(y) = H^2[(4y - y^2) f'' + D(2 - y) f'] in explicit coordinates, D = 4, 5.
    (h) The six first and mixed second derivatives of y (delta -> 0) used later (derivative table of the
        regular part in the Task 8 section of the notes).

Expected result
    All identities True; numeric nascent-delta integrals 9.6498..., 9.8473..., 9.8673... -> pi^2 = 9.8696...

Runtime
    About 1-3 seconds.

Origin
    Final version of check_52.py of Task 8 (15 Oct 2026).  The original one-line header comment was replaced by
    this docstring; the code is unchanged.  (An earlier run in Task 8 printed check (c) in an unsimplified form;
    the expected output is that of the final version.)
"""
import sympy as sp, mpmath as mp
I = sp.I
D, H, dl = sp.symbols('D H delta', positive=True)

print("== (a) ln(aa') drops out of mixed derivatives ==")
t, tp = sp.symbols('t tp', negative=True)
a, ap = -1/(H*t), -1/(H*tp)
print("  d_0 d'_0 ln(a a') =", sp.simplify(sp.diff(sp.log(a*ap), t, tp)))

print("== (b) jump lemma: d_0 d'_0 (|deta| - i delta)^2 = -2 + 4 i delta DiracDelta(deta) ==")
u, v = sp.symbols('u v', real=True)
e = sp.diff((sp.Abs(u - v) - I*dl)**2, u, v)
e = sp.expand(e).subs(sp.sign(u - v)**2, 1)
e = e.replace(lambda z: isinstance(z, sp.Mul) and z.has(sp.DiracDelta) and z.has(sp.Abs), lambda z: z.subs(sp.Abs(u - v), 0))
print("  result:", sp.simplify(e))

print("== (c) flat check: eta^{rs}[F'' d_r X d'_s X + F' (d_r d'_s X)_reg] = 0 for F = X^(1-D/2), delta>0, on both sides of deta=0 ==")
x1, x2, x3, y1, y2, y3 = sp.symbols('x1 x2 x3 y1 y2 y3', real=True)
for s in (+1, -1):                                    # |deta| = s*deta on each side
    deta = u - v
    X = (x1 - y1)**2 + (x2 - y2)**2 + (x3 - y3)**2 - (s*deta - I*dl)**2
    xs, ys = (u, x1, x2, x3), (v, y1, y2, y3)
    eta = [-1, 1, 1, 1]
    Fp = (1 - D/2)*X**(-D/2); Fpp = (1 - D/2)*(-D/2)*X**(-D/2 - 1)
    grad = sp.simplify(sum(eta[r]*sp.diff(X, xs[r])*sp.diff(X, ys[r]) for r in range(4)))   # = -4 X (same in any D)
    regmix = [sp.simplify(sp.diff(X, xs[r], ys[r])) for r in range(4)]                       # = -2 eta_rr (regular part)
    tot_D = sp.simplify(Fpp*(-4*X) + Fp*(-2)*D)                                             # eta^{rs}(-2 eta_rs) = -2D
    print(f"  side s={s:+d}: eta^rs d_r X d'_s X + 4X = {sp.simplify(grad + 4*X)};  (d_r d'_r X)_reg = {regmix};"
          f"  contracted chain part (general D) = {tot_D}")

print("== (d) nascent delta: delta (r^2+delta^2)^(-D/2) -> pi^(D/2)/Gamma(D/2) * delta^{D-1}(x) ==")
s_ = sp.symbols('s', positive=True)
for Dv in (3, 4, 5, 6):
    Om = 2*sp.pi**(sp.Rational(Dv - 1, 2))/sp.gamma(sp.Rational(Dv - 1, 2))
    val = sp.simplify(Om*sp.integrate(s_**(Dv - 2)*(1 + s_**2)**sp.Rational(-Dv, 2), (s_, 0, sp.oo)))
    print(f"  D={Dv}: Omega_(D-2) * int s^(D-2)(1+s^2)^(-D/2) ds = {val}   vs  pi^(D/2)/Gamma(D/2) = {sp.simplify(sp.pi**sp.Rational(Dv, 2)/sp.gamma(sp.Rational(Dv, 2)))}")
gen = 2*sp.pi**((D - 1)/2)/sp.gamma((D - 1)/2)*sp.Rational(1, 2)*sp.gamma((D - 1)/2)*sp.gamma(sp.Rational(1, 2))/sp.gamma(D/2)
print("  general D (Beta function):", sp.simplify(sp.gammasimp(gen - sp.pi**(D/2)/sp.gamma(D/2))) == 0)
mp.mp.dps = 20
for d_ in ('1e-2', '1e-3', '1e-4'):
    dd = mp.mpf(d_)
    val = mp.quad(lambda r: 4*mp.pi*r**2*mp.exp(-r**2)*dd/(r**2 + dd**2)**2, [0, dd, 10*dd, 1, mp.inf])
    print(f"  numeric D=4, delta={d_}: int d^3x e^(-r^2) delta/(r^2+delta^2)^2 = {mp.nstr(val, 10)}   (pi^2 = {mp.nstr(mp.pi**2, 10)})")

print("== (e) coefficient of the local term ==")
aa, r = sp.symbols('a r', positive=True)
yy = sp.Symbol('y', positive=True)
def coeff(c0):     # c0 = coefficient of (4/y)^(D/2-1) inside the braces of (14)
    A_lead = H**(D - 2)/(4*sp.pi)**(D/2)*c0*(4/yy)**(D/2 - 1)
    Ap = sp.diff(A_lead, yy).subs(yy, aa**2*H**2*(r**2 + dl**2))        # A'(y) at deta=0, a'=a
    loc = -4*I*dl*aa**2*H**2*Ap                                          # -4 i delta a a' H^2 A'(y) delta(deta)
    # loc = K * delta*(r^2+delta^2)^(-D/2); replace delta*(r^2+delta^2)^(-D/2) by pi^(D/2)/Gamma(D/2) delta^D
    K = sp.simplify(loc/(dl*(r**2 + dl**2)**(-D/2)))
    return sp.simplify(sp.gammasimp(K*sp.pi**(D/2)/sp.gamma(D/2)))
c_correct = sp.gamma(D/2)/(D/2 - 1)          # = Gamma(D/2-1): leading term = i Delta_cf, Eq. (26)
c_printed = sp.gamma(D/2 - 1)/(D/2 - 1)      # as printed in (14)
print("  leading term = i Delta_cf (26):  coefficient =", sp.simplify(coeff(c_correct)), "  [expected i a^(2-D)]",
      sp.simplify(coeff(c_correct) - I*aa**(2 - D)) == 0)
print("  leading term as printed in (14): coefficient =", sp.simplify(coeff(c_printed)), " = i a^(2-D) * 2/(D-2)")

print("== (f) ODE recursion fixes the leading coefficient of A(y) from its second term ==")
c0, c1 = sp.symbols('c0 c1')
# A ~ c0 (4/y)^(D/2-1) + c1 (4/y)^(D/2-2) + ...; L[y^s] = A_s y^(s-1) + B_s y^s with L = (4y-y^2) d^2 + D(2-y) d
Ls = lambda s: (s*(4*(s - 1) + 2*D), -s*((s - 1) + D))           # L[y^s] = A_s y^{s-1} + B_s y^s
A0, B0 = Ls(1 - D/2); A1, B1 = Ls(2 - D/2)
print("  indicial (coefficient of y^(-D/2) from y^(1-D/2)) :", sp.simplify(A0))
cond = sp.simplify(c0*4**(D/2 - 1)*B0 + c1*4**(D/2 - 2)*A1)        # total coefficient of y^(1-D/2)
sol = sp.solve(sp.Eq(cond, 0), c0)[0]
c1_paper = sp.gamma(D/2 + 1)/(D/2 - 2)
print("  c0 required by c1 = Gamma(D/2+1)/(D/2-2):", sp.simplify(sp.gammasimp(sol.subs(c1, c1_paper))),
      "  == Gamma(D/2-1)?", sp.simplify(sp.gammasimp(sol.subs(c1, c1_paper) - sp.gamma(D/2 - 1))) == 0,
      "  == printed Gamma(D/2-1)/(D/2-1)?", sp.simplify(sp.gammasimp(sol.subs(c1, c1_paper) - c_printed)) == 0)

print("== (g) Box f(y) = H^2[(4y-y^2) f'' + D(2-y) f'] in explicit coordinates (D=4,5) ==")
fF = sp.Function('f')
for Dv in (4, 5):
    X = sp.symbols(f'X0:{Dv}', real=True); Y = sp.symbols(f'Y0:{Dv}', real=True)
    t0 = sp.Symbol('t0', negative=True); s0 = sp.Symbol('s0', negative=True)
    aX, aY = -1/(H*t0), -1/(H*s0)
    yv = aX*aY*H**2*(sum((X[i] - Y[i])**2 for i in range(1, Dv)) - (t0 - s0)**2)
    coords = (t0,) + X[1:]
    etad = [-1] + [1]*(Dv - 1)
    F = fF(yv)
    box = sum(sp.diff(aX**(Dv - 2)*etad[m]*sp.diff(F, coords[m]), coords[m]) for m in range(Dv))/aX**Dv
    yS = sp.Symbol('yS')
    target = H**2*((4*yv - yv**2)*sp.diff(fF(yS), yS, 2).subs(yS, yv) + Dv*(2 - yv)*sp.diff(fF(yS), yS).subs(yS, yv))
    print(f"  D={Dv}:", sp.simplify(box - target) == 0)

print("== (h) table of first and mixed second derivatives of y (delta -> 0), D=4 ==")
X = sp.symbols('X1:4', real=True); Y = sp.symbols('Y1:4', real=True)
yv = a*ap*H**2*(sum((X[i] - Y[i])**2 for i in range(3)) - (t - tp)**2)
de = t - tp
chk = {
 'd0 y   = H a y - 2 a a\' H^2 deta':        sp.diff(yv, t) - (H*a*yv - 2*a*ap*H**2*de),
 'd\'0 y  = H a\' y + 2 a a\' H^2 deta':     sp.diff(yv, tp) - (H*ap*yv + 2*a*ap*H**2*de),
 'd0 d\'0 y = a a\' H^2 [y + 2 + 2 H^2 a a\' deta^2]': sp.diff(yv, t, tp) - a*ap*H**2*(yv + 2 + 2*H**2*a*ap*de**2),
 'd0 d\'j y = -2 H^3 a^2 a\' dx_j':          sp.diff(yv, t, Y[0]) - (-2*H**3*a**2*ap*(X[0] - Y[0])),
 'di d\'0 y = 2 H^3 a a\'^2 dx_i':           sp.diff(yv, X[0], tp) - (2*H**3*a*ap**2*(X[0] - Y[0])),
 'di d\'j y = -2 a a\' H^2 delta_ij':        sp.diff(yv, X[0], Y[0]) + 2*a*ap*H**2,
}
for kname, val in chk.items():
    print(f"  {kname:48s}:", sp.simplify(val) == 0)
