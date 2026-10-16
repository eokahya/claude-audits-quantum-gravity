"""
check_81_93.py -- independent checks of Eqs. (I.81)-(I.93) and Table III of
Kahya & Woodard, Phys. Rev. D 76, 124005 (2007)  [Paper I].

Purpose
    Task 12 of the re-derivation.  The time-derivative identities (I.81)-(I.88),
    the reduction (I.89) of the nonlocal logarithm term of (I.67), Table III, the
    d'Alembertian (I.90), the reduction (I.91)-(I.92) of a mixed A term and the
    identity (I.93) were derived by hand (see notes, chapter "Normal contributions
    and Table III: Eqs. (I.68)-(I.93)", section "Task 12").  This script checks them
    with SymPy; it does NOT replace the analytic derivation.

What is checked (labels as in the notes / printed output)
    (a) (I.81)-(I.83): exact, symbolic.
    (b) (I.84), (I.85), (I.87): symbolic with a generic function f(y);
        (I.86), (I.88): f = y^alpha at a random point, 25 digits, D = 4.
    (c) (I.90): general D, from the definition a^D Box = d_mu(a^{D-2} eta^{mu nu} d_nu)
        of (I.42), acting on a test function g(eta, x).
    (d) (I.91)-(I.92): A = y^beta, D = 4 and 5; left side differentiated by brute
        force, Box applied in explicit coordinates via (I.42).
    (e) (I.93): Box and Box' on f = y^alpha, D = 4 and 5.
    (f) (I.89) and Table III at D = 4 with the true A(y) of (I.62), at three random
        points with H != 1:
          - derived nonlocal term of (I.67) + printed Table III  = 0  (so derived = -paper),
          - (I.89) with the H^2 of its second term restored  = -derived,
          - (I.89) as printed (no H^2)                        != -derived,
        and the expansion of the three Table III coefficients
        (-12 A/h, -4 y A'/h, 2 I[A]/h in terms of x = y/4).

Sign conventions
    "derived(67)" is the nonlocal term of (I.67) as re-derived in Task 10, i.e. in the
    *definition sign* fixed by (I.168).  "Table III(paper)" and "(89)" are the printed
    expressions (paper convention).  They differ by an overall factor -1 inherited from
    (I.47); see the notes.

Notes on the random points
    In the first run recorded in the task log, H was drawn as randint(5,20)/10 and took
    the value 1 at two of the three points, where the printed (I.89) (missing H^2)
    accidentally agreed.  The final version (this file) draws H from {7/10, 13/10, 3/2,
    17/10}, so H != 1 always.  The derived term and Table III are independent of H at
    fixed (eta, eta', x, x'), which is why the printed scales are unchanged.

Expected result
    Every check prints True except "(89) as printed (no H^2) == -derived", which prints
    False at all three points; "derived(67)_nonlocal + Table III(paper)" is 0.0 at all
    three points.  Coefficients:
        -12*A/h   = 24*log(x) + 12 - 12/x
        -4*y A'/h = 8 + 4/x
        2*I[A]/h  = -16*x*log(x) + 8*x + 8*log(x) + 16*log(2)

Runtime
    about 20 s (SymPy 1.14, Python 3.13).

Run
    python3 check_81_93.py
"""
import sympy as sp, random
R = sp.Rational
H = sp.Symbol('H', positive=True)
al, be = sp.symbols('alpha beta')
Yv = sp.Symbol('Y', positive=True)
Iop = lambda g, n=1: [g := sp.integrate(g, Yv, conds='none') for _ in range(n)][-1]   # n-fold indefinite integral in y

def setup(Dv):
    t, tp = sp.symbols('t tp', negative=True)
    xs = sp.symbols(f'x1:{Dv}', real=True); ys = sp.symbols(f'y1:{Dv}', real=True)
    a, ap = -1/(H*t), -1/(H*tp)
    r2 = sum((xs[i] - ys[i])**2 for i in range(Dv - 1))
    y = a*ap*H**2*(r2 - (t - tp)**2)
    return t, tp, xs, ys, a, ap, r2, y
random.seed(5)
def rnd(t, tp, xs, ys):
    pt = {H: random.choice([R(7, 10), R(13, 10), R(3, 2), R(17, 10)]), t: -R(random.randint(3, 9), 4), tp: -R(random.randint(10, 19), 4),
          al: R(3, 7), be: -R(5, 3)}
    for v in xs + ys: pt[v] = R(random.randint(-9, 9), 5)
    return pt
num = lambda e, pt: sp.N(e.subs(pt), 25)
F = lambda e, y: e.subs(Yv, y)

t, tp, xs, ys, a, ap, r2, y = setup(4)
de = t - tp
f = sp.Function('f')
fy = f(y)
print("== (81)-(83) ==")
print("  (81) d'_i y = -d_i y                          :", all(sp.simplify(sp.diff(y, ys[i]) + sp.diff(y, xs[i])) == 0 for i in range(3)))
print("  (82) d_0 y = Ha(y-2a'H deta) = Ha(y-2+2a'/a)  :", sp.simplify(sp.diff(y, t) - H*a*(y - 2*ap*H*de)) == 0, sp.simplify(sp.diff(y, t) - H*a*(y - 2 + 2*ap/a)) == 0)
print("  (83) d'_0 y = Ha'(y+2aH deta) = Ha'(y-2+2a/a'):", sp.simplify(sp.diff(y, tp) - H*ap*(y + 2*a*H*de)) == 0, sp.simplify(sp.diff(y, tp) - H*ap*(y - 2 + 2*a/ap)) == 0)

print("== (84)-(88), generic f(y) and f = y^alpha for I[f] ==")
Y = sp.Symbol('Y')
d0 = lambda e: sp.diff(e, t); d0p = lambda e: sp.diff(e, tp)
fpY = sp.Derivative(f(Y), Y)
print("  (84) (d0+d0')f = H(a+a') y f'        :", sp.simplify(d0(fy) + d0p(fy) - (H*(a + ap)*y*fpY.subs(Y, y)).doit()) == 0)
print("  (85) (a'd0+a d0')f = 2Haa' aa'H^2 r^2 f':", sp.simplify(ap*d0(fy) + a*d0p(fy) - (2*H*a*ap*a*ap*H**2*r2*fpY.subs(Y, y)).doit()) == 0)
fpow = Yv**al; Ifp = Iop(fpow)
lap = lambda e: sum(sp.diff(e, xi, 2) for xi in xs)
e86 = ap*d0(F(fpow, y)) + a*d0p(F(fpow, y)) - (-(4 - 1)*H*a*ap*F(fpow, y) + lap(F(Ifp, y))/(2*H))
e88 = a*d0(F(fpow, y)) + ap*d0p(F(fpow, y)) - (H*(a + ap)**2*y*F(sp.diff(fpow, Yv), y) + (4 - 1)*H*a*ap*F(fpow, y) - lap(F(Ifp, y))/(2*H))
pt = rnd(t, tp, xs, ys)
print("  (86) D=4, f=y^alpha                  :", abs(num(e86, pt)) < 1e-18)
print("  (87) operator identity               :", sp.simplify(a*d0(fy) + ap*d0p(fy) - ((a + ap)*(d0(fy) + d0p(fy)) - (ap*d0(fy) + a*d0p(fy)))) == 0)
print("  (88) D=4, f=y^alpha                  :", abs(num(e88, pt)) < 1e-18)

print("== (90): Box from (42) == -a^-2 d0^2 - (D-2)(H/a) d0 + a^-2 Lap, general D ==")
D = sp.Symbol('D'); g = sp.Function('g')(t, *xs)
boxg = (-sp.diff(a**(D - 2)*sp.diff(g, t), t) + a**(D - 2)*lap(g))/a**D
print("  ", sp.simplify(boxg - (-sp.diff(g, t, 2)/a**2 - (D - 2)*H/a*sp.diff(g, t) + lap(g)/a**2)) == 0)

print("== (91)-(92), A = y^beta (pointwise), D=4,5; Box applied via (42) in explicit coordinates ==")
for Dv in (4, 5):
    t, tp, xs, ys, a, ap, r2, y = setup(Dv)
    A = Yv**be
    I2AA = Iop(A*sp.diff(A, Yv, 2), 2); I2Ap2 = Iop(sp.diff(A, Yv)**2, 2)
    Ay = F(A, y)
    lapx = lambda e: sum(sp.diff(e, xi, 2) for xi in xs)
    lapxp = lambda e: sum(sp.diff(e, yi, 2) for yi in ys)
    box = lambda e: (-sp.diff(a**(Dv - 2)*sp.diff(e, t), t) + a**(Dv - 2)*lapx(e))/a**Dv
    lhs = sum(sp.diff((a*ap)**(Dv - 2)*Ay*sp.diff(Ay, t, ys[i]), t, ys[i]) for i in range(Dv - 1))
    G, K = F(I2AA, y), F(I2Ap2, y)
    r91 = lapxp(sp.diff((a*ap)**(Dv - 2)*sp.diff(G, t), t)) + H*lapxp(sp.diff(a**(Dv - 1)*ap**(Dv - 2)*K, t))
    r92 = (-a**Dv*ap**(Dv - 2)*lapx(box(G)) + (a*ap)**(Dv - 2)*lapx(lapx(G))
           + H*a**(Dv - 1)*ap**(Dv - 2)*lapx(sp.diff(K, t)) + (Dv - 1)*H**2*a**Dv*ap**(Dv - 2)*lapx(K))
    pt = rnd(t, tp, xs, ys); sc = abs(num(lhs, pt))
    print(f"  D={Dv}: (91)", abs(num(lhs - r91, pt)) < 1e-15*max(1, sc), "  (92)", abs(num(lhs - r92, pt)) < 1e-15*max(1, sc))

print("== (93): Box f(y) = H^2[(4y-y^2)f'' + D(2-y)f'] = Box' f(y), f = y^alpha, D=4,5 ==")
for Dv in (4, 5):
    t, tp, xs, ys, a, ap, r2, y = setup(Dv)
    fy_ = F(Yv**al, y)
    lapx = lambda e: sum(sp.diff(e, xi, 2) for xi in xs); lapxp = lambda e: sum(sp.diff(e, yi, 2) for yi in ys)
    box = (-sp.diff(a**(Dv - 2)*sp.diff(fy_, t), t) + a**(Dv - 2)*lapx(fy_))/a**Dv
    boxp = (-sp.diff(ap**(Dv - 2)*sp.diff(fy_, tp), tp) + ap**(Dv - 2)*lapxp(fy_))/ap**Dv
    tgt = H**2*F((4*Yv - Yv**2)*sp.diff(Yv**al, Yv, 2) + Dv*(2 - Yv)*sp.diff(Yv**al, Yv), y)
    pt = rnd(t, tp, xs, ys)
    print(f"  D={Dv}: Box", abs(num(box - tgt, pt)) < 1e-15, "  Box'", abs(num(boxp - tgt, pt)) < 1e-15)

print("== (89) and Table III: pointwise, D=4, actual A(y) of (62) ==")
t, tp, xs, ys, a, ap, r2, y = setup(4)
pi = sp.pi
h = H**2/(16*pi**2)
A62 = h*(4/Yv - 2*sp.log(Yv/4) - 1)
Ay = F(A62, y); x = y/4
lapx = lambda e: sum(sp.diff(e, xi, 2) for xi in xs)
L2A = lapx(Ay)
derived67 = H**2/(8*pi**2)*2*H*(a*ap)**2*(a*sp.diff(L2A, t) + ap*sp.diff(L2A, tp))       # derived (6.1) of Task 10, non-local term / kappa^2
pref = H**4/(4*pi)**4                                                                       # kappa^2 H^4/(4 pi)^4 without kappa^2
tab_paper = pref*((a*ap)**3*H**2*lapx(-12/x + 24*sp.log(x)) + (a*ap)**2*(a + ap)**2*H**2*lapx(4/x)
                  + (a*ap)**2*lapx(lapx(8*sp.log(x) - 16*x*sp.log(x))))
yA = F(Yv*sp.diff(A62, Yv), y); IA = F(Iop(A62), y)
p89_printed = H**2/(16*pi**2)*(-12*(a*ap)**3*H**2*lapx(Ay) - 4*(a*ap)**2*(a + ap)**2*lapx(yA) + 2*(a*ap)**2*lapx(lapx(IA)))
p89_withH2 = H**2/(16*pi**2)*(-12*(a*ap)**3*H**2*lapx(Ay) - 4*H**2*(a*ap)**2*(a + ap)**2*lapx(yA) + 2*(a*ap)**2*lapx(lapx(IA)))
for k in range(3):
    pt = rnd(t, tp, xs, ys); sc = abs(num(derived67, pt))
    print(f"  point {k}: derived(67)_nonlocal + Table III(paper) = {sp.Float(num(derived67 + tab_paper, pt), 3)}"
          f"  | derived - Table III(paper) = {sp.Float(num(derived67 - tab_paper, pt), 3)}  (scale {sp.Float(sc, 3)})")
    print(f"           (89) with H^2 restored == -derived : {abs(num(p89_withH2 + derived67, pt)) < 1e-15*sc}"
          f"   (89) as printed (no H^2) == -derived : {abs(num(p89_printed + derived67, pt)) < 1e-15*sc}")
print("  Table III coefficients from A(y) (dropping terms killed by Lap / Lap^2):")
Ah = A62/h
print("    -12*A/h        =", sp.expand(sp.expand_log((-12*Ah).subs(Yv, 4*sp.Symbol('x', positive=True)), force=True)))
print("    -4*y A'/h      =", sp.expand((-4*Yv*sp.diff(Ah, Yv)).subs(Yv, 4*sp.Symbol('x', positive=True))))
print("    2*I[A]/h       =", sp.expand(sp.expand_log((2*Iop(Ah)).subs(Yv, 4*sp.Symbol('x', positive=True)), force=True)))
