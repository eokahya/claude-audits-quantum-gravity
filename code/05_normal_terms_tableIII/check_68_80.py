"""
check_68_80.py -- independent checks of the derivative-extraction identities
Eqs. (I.68)-(I.80) of Kahya & Woodard, Phys. Rev. D 76, 124005 (2007)  [Paper I].

Purpose
    Task 11 of the re-derivation.  The identities (I.68)-(I.80) were derived by
    hand (see notes, chapter "Normal contributions and Table III: Eqs. (I.68)-(I.93)",
    section "Task 11").  This script checks them with SymPy; it does NOT replace the
    analytic derivation.

What is checked (labels as in the notes / printed output)
    (A) Eq. (I.70): f{A'' dy d'y + A' dd'y} = dd' I^2[fA''] + dd'y I[f'A'],
        all 16 (rho, sigma) components, D = 4, plus the integration-by-parts
        consistency I[fA''] = fA' - I[f'A'].
    (B) Eqs. (I.71)-(I.72): the mixed second derivative of y(x;x') with the
        lower-index convention Dx_sigma = eta_{sigma nu}(x - x')^nu (exact, symbolic).
    (C) Eq. (I.74), D = 4 and 5.
    (D) Eqs. (I.75)-(I.78) (naive derivatives), D = 4 and 5; (I.78) for i = j and i != j.
    (E) Eqs. (I.79)-(I.80), D = 4 and 5, left side differentiated by brute force.
    (F) Delta-function content (notes, Task 11, subsection "Delta content of (I.70)";
        Eq. (4.1) of the original Task 11 report): the y^{1-D/2} component of
        I^2[f A''] equals c_f * A_lead, where c_f is the constant part of f.
    (G) Power counting after (I.80): which powers y^s of I^2[AA''] and I[A'^2]
        are not locally integrable at D = 4 (criterion s > -D/2).

Method
    f = y^alpha, A = y^beta with generic exponents (so no integration constants
    arise); explicit conformal coordinates and the true y(x;x') =
    a a' H^2 (|x - x'|^2 - (eta - eta')^2), a = -1/(H eta).  Both sides are
    evaluated at random rational points (alpha = 3/7, beta = -5/3) with 25 digits;
    (B) and the integration-by-parts line of (A) are exact symbolic checks.
    (F) and (G) use symbolic D with model power series for f and A.

Expected result
    Every check in (A)-(E) prints True; (F) prints the y^{1-D/2} part
    2**(D-2)*N*Y**(1-D/2)*c*c0 = c * N * c0 * 4^(D/2-1) * y^(1-D/2) and True;
    (G) flags only y^{2-D} in I^2[AA''] and y^{1-D}, y^{2-D} in I[A'^2] as DIVERGENT.

Runtime
    about 15-20 s (SymPy 1.14, Python 3.13).

Run
    python3 check_68_80.py
"""
import sympy as sp, random
I_, R = sp.I, sp.Rational
H = sp.Symbol('H', positive=True)
al, be = sp.symbols('alpha beta')          # f = y^alpha, A = y^beta (generic powers; no integration constants arise)
Yv = sp.Symbol('Y', positive=True)

def Iop(g, n=1):                            # indefinite integral in y, n times (pure powers -> no constants)
    for _ in range(n):
        g = sp.integrate(g, Yv, conds='none')
    return sp.simplify(g)

def setup(Dv):
    t, tp = sp.symbols('t tp', negative=True)
    xs = sp.symbols(f'x1:{Dv}', real=True); ys = sp.symbols(f'y1:{Dv}', real=True)
    a, ap = -1/(H*t), -1/(H*tp)
    r2 = sum((xs[i] - ys[i])**2 for i in range(Dv - 1))
    y = a*ap*H**2*(r2 - (t - tp)**2)
    X = (t,) + xs; Xp = (tp,) + ys
    return t, tp, xs, ys, a, ap, r2, y, X, Xp

def F_of(expr_in_Y, y):                     # substitute Y -> y(x;x')
    return expr_in_Y.subs(Yv, y)

random.seed(11)
def rnd_point(t, tp, xs, ys):
    pt = {H: R(random.randint(5, 20), 10), t: -R(random.randint(3, 9), 4), tp: -R(random.randint(10, 19), 4)}
    for v in xs + ys: pt[v] = R(random.randint(-9, 9), 5)
    return pt
num = lambda e, pt: sp.N(e.subs(pt).subs({al: R(3, 7), be: -R(5, 3)}), 25)

print("== (A) Eq. (70): f{A'' dy d'y + A' dd'y} = dd' I^2[fA''] + dd'y I[f'A'],  f=y^alpha, A=y^beta, D=4 ==")
t, tp, xs, ys, a, ap, r2, y, X, Xp = setup(4)
f = Yv**al; A = Yv**be
I2fA2 = Iop(f*sp.diff(A, Yv, 2), 2); I1fpAp = Iop(sp.diff(f, Yv)*sp.diff(A, Yv), 1)
print("  integration by parts consistency  I[fA''] == fA' - I[f'A'] :", sp.simplify(Iop(f*sp.diff(A, Yv, 2)) - (f*sp.diff(A, Yv) - I1fpAp)) == 0)
ok = True
for rho in range(4):
    for sig in range(4):
        lhs = F_of(f, y)*(F_of(sp.diff(A, Yv, 2), y)*sp.diff(y, X[rho])*sp.diff(y, Xp[sig]) + F_of(sp.diff(A, Yv), y)*sp.diff(y, X[rho], Xp[sig]))
        rhs = sp.diff(F_of(I2fA2, y), X[rho], Xp[sig]) + sp.diff(y, X[rho], Xp[sig])*F_of(I1fpAp, y)
        pt = rnd_point(t, tp, xs, ys)
        ok &= abs(num(lhs - rhs, pt)) < 1e-18*max(1, abs(num(lhs, pt)))
print("  all 16 components at random points :", ok)

print("== (B) Eqs. (71)-(73): mixed second derivative of y, lower-index Dx_sigma = eta_{sigma nu}(x-x')^nu ==")
eta = [-1, 1, 1, 1]
Dxl = [eta[m]*(X[m] - Xp[m]) for m in range(4)]
d0 = lambda m: 1 if m == 0 else 0
ok = all(sp.simplify(sp.diff(y, X[r], Xp[s]) - H**2*a*ap*(y*d0(r)*d0(s) - 2*a*d0(r)*H*Dxl[s] + 2*ap*H*Dxl[r]*d0(s) - 2*(eta[r] if r == s else 0))) == 0
         for r in range(4) for s in range(4))
print("  (71) all components :", ok)
print("  (72a) d0 d0' y = H^2 aa'[2 - y + 2 aa' H^2 |Dx|^2] :", sp.simplify(sp.diff(y, t, tp) - H**2*a*ap*(2 - y + 2*a*ap*H**2*r2)) == 0)

print("== (C) Eq. (74): aa'H^2|Dx|^2 f = -(D-1)/2 I[f] - grad.grad' I^2[f]/(4aa'H^2), D=4,5 ==")
for Dv in (4, 5):
    t, tp, xs, ys, a, ap, r2, y, X, Xp = setup(Dv)
    g2 = F_of(Iop(f, 2), y)
    lhs = a*ap*H**2*r2*F_of(f, y)
    rhs = -R(Dv - 1, 2)*F_of(Iop(f, 1), y) - sum(sp.diff(g2, xs[i], ys[i]) for i in range(Dv - 1))/(4*a*ap*H**2)
    pt = rnd_point(t, tp, xs, ys)
    print(f"  D={Dv}:", abs(num(lhs - rhs, pt)) < 1e-18)

print("== (D) Eqs. (75)-(78) (naive derivatives), f=y^alpha, A=y^beta, D=4,5 ==")
for Dv in (4, 5):
    t, tp, xs, ys, a, ap, r2, y, X, Xp = setup(Dv)
    Fy = lambda e: F_of(e, y)
    naive = lambda r, s: Fy(sp.diff(A, Yv, 2))*sp.diff(y, X[r])*sp.diff(y, Xp[s]) + Fy(sp.diff(A, Yv))*sp.diff(y, X[r], Xp[s])
    J1 = Iop(sp.diff(f, Yv)*sp.diff(A, Yv), 1); J2 = Iop(J1, 1); J3 = Iop(J2, 1)
    gg = lambda e: sum(sp.diff(Fy(e), xs[i], ys[i]) for i in range(Dv - 1))
    e75 = Fy(f)*naive(0, 0) - (sp.diff(Fy(I2fA2), t, tp) - R(1, 2)*gg(J3) + H**2*a*ap*((2 - y)*Fy(J1) - (Dv - 1)*Fy(J2)))
    e76 = Fy(f)*naive(0, 1) - (sp.diff(Fy(I2fA2), t, ys[0]) + H*a*sp.diff(Fy(J2), ys[0]))
    e77 = Fy(f)*naive(1, 0) - (sp.diff(Fy(I2fA2), xs[0], tp) + H*ap*sp.diff(Fy(J2), xs[0]))
    e78a = Fy(f)*naive(1, 1) - (sp.diff(Fy(I2fA2), xs[0], ys[0]) - 2*H**2*a*ap*Fy(J1))
    e78b = Fy(f)*naive(1, 2) - sp.diff(Fy(I2fA2), xs[0], ys[1])
    pt = rnd_point(t, tp, xs, ys)
    print(f"  D={Dv}: (75)", abs(num(e75, pt)) < 1e-15, " (76)", abs(num(e76, pt)) < 1e-15, " (77)", abs(num(e77, pt)) < 1e-15,
          " (78) i=j", abs(num(e78a, pt)) < 1e-15, " i!=j", abs(num(e78b, pt)) < 1e-15)

print("== (E) Eqs. (79)-(80): grad.grad'[(aa')^{D-2} A grad.grad' A], A=y^beta, D=4,5 ==")
for Dv in (4, 5):
    t, tp, xs, ys, a, ap, r2, y, X, Xp = setup(Dv)
    Fy = lambda e: F_of(e, y)
    Ay = Fy(A)
    gg = lambda e: sum(sp.diff(e, xs[i], ys[i]) for i in range(Dv - 1))
    lap = lambda e: sum(sp.diff(e, xs[i], 2) for i in range(Dv - 1))
    lhs = gg((a*ap)**(Dv - 2)*Ay*gg(Ay))
    I2AA = Fy(Iop(A*sp.diff(A, Yv, 2), 2)); IAp2 = Fy(Iop(sp.diff(A, Yv)**2, 1))
    r79 = (a*ap)**(Dv - 2)*gg(gg(I2AA)) - 2*(Dv - 1)*H**2*(a*ap)**(Dv - 1)*gg(IAp2)
    r80 = (a*ap)**(Dv - 2)*lap(lap(I2AA)) + 2*(Dv - 1)*H**2*(a*ap)**(Dv - 1)*lap(IAp2)
    pt = rnd_point(t, tp, xs, ys)
    sc = abs(num(lhs, pt))
    print(f"  D={Dv}: (79)", abs(num(lhs - r79, pt)) < 1e-15*max(1, sc), "  (80)", abs(num(lhs - r80, pt)) < 1e-15*max(1, sc))

print("== (F) delta content: the y^{1-D/2} part of I^2[f A''] is c_f * A_lead (c_f = constant part of f) ==")
D = sp.Symbol('D'); Nn, c, d, e, c0, c1 = sp.symbols('N c d e c0 c1')
fM = c + d*(4/Yv)**(D/2 - 1) + e*(4/Yv)**(D/2 - 2)                 # model outer function: constant + propagator-type powers
AM = Nn*(c0*(4/Yv)**(D/2 - 1) + c1*(4/Yv)**(D/2 - 2))
prod = sp.expand(sp.powsimp(sp.expand(fM*sp.diff(AM, Yv, 2)), force=True))
terms = sp.Add.make_args(prod)
I2 = sum(Iop(tm, 2) for tm in terms)
I2 = sp.expand(sp.powsimp(sp.expand(I2), force=True))
lead = Nn*c0*4**(D/2 - 1)*Yv**(1 - D/2)
coef = sp.simplify(sum(tm for tm in sp.Add.make_args(I2) if sp.simplify(sp.powsimp(tm/Yv**(1 - D/2), force=True)).free_symbols.isdisjoint({Yv})))
print("  y^{1-D/2} part of I^2[fA''] :", sp.factor(coef), "   expected c*N*c0*4^(D/2-1)*Y^(1-D/2):", sp.simplify(coef - c*lead) == 0)
expo = lambda tm: sp.simplify(sp.powsimp(sp.expand_power_base(sp.simplify(Yv*sp.diff(tm, Yv)/tm), force=True), force=True))
print("  powers present in I^2[fA''] (exponent of Y):", sorted({expo(tm) for tm in sp.Add.make_args(I2)}, key=str))

print("== (G) power counting after (80), D=4: y^s is locally integrable in D dims iff s > -D/2 ==")
A3 = Nn*(c0*(4/Yv)**(D/2 - 1) + c1*(4/Yv)**(D/2 - 2)) + c + d*Yv
def powers(expr):
    ex = sp.expand(sp.powsimp(sp.expand(expr), force=True))
    return sorted({expo(tm) for tm in sp.Add.make_args(ex) if tm != 0}, key=lambda s: (float(s.subs(D, 4)), str(s)))   # secondary key: deterministic order of ties
pAA = powers(sum(Iop(tm, 2) for tm in sp.Add.make_args(sp.expand(sp.powsimp(sp.expand(A3*sp.diff(A3, Yv, 2)), force=True)))))
pAp = powers(sum(Iop(tm, 1) for tm in sp.Add.make_args(sp.expand(sp.powsimp(sp.expand(sp.diff(A3, Yv)**2), force=True)))))
print("  I^2[AA'']:", [(s, 'DIVERGENT' if s.subs(D, 4) <= -2 else 'finite') for s in pAA])
print("  I[A'^2]  :", [(s, 'DIVERGENT' if s.subs(D, 4) <= -2 else 'finite') for s in pAp])
