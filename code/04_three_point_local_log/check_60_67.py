#!/usr/bin/env python3
"""
check_60_67.py -- code checks for the logarithm contributions of the 3-point
diagram, Eqs. (I.60)-(I.67) of Kahya & Woodard, Phys. Rev. D 76, 124005 (2007)
("Paper I").

Purpose
    Supplementary (not a substitute for the analytic derivation in the notes,
    chapter "3-point local and logarithm terms: Eqs. (I.54)-(I.67)", Task 10).

What it checks
    (a) (I.62): the first three terms of the propagator function A(y) of (I.14),
        expanded at D = 4 + eps with the corrected first coefficient
        c0 = Gamma(D/2 - 1), give H^2/(16 pi^2) {4/y - 2 ln(y/4) - 1};
        the n-th series term of (I.14) vanishes identically at D = 4;
        the printed c0 = Gamma(D/2-1)/(D/2-1) changes A only at O(D-4).
    (b) A(y) of (I.62) solves (4y - y^2) A'' + 4(2 - y) A' = 3k with k = 2h,
        h = H^2/(16 pi^2), i.e. Box A = (D-1) k H^2 at D = 4.
    (c) Distributional identities on a test function g, symbolic D
        (l = ln a, p = D - 2):
            -2i d0[a^p l d0 delta] = i ln(aa') (a^D Box - a^{D-2} Lap) delta
                                     + i (D-1) H^2 a^D delta          (lemma D5)
            ln(aa') Lap delta = 2 l Lap delta                         (D3)
            d0(a^p d0 ln a) = (D-1) H^2 a^D
    (d) Pointwise check at D = 4 of the NON-LOCAL parts of the derived expressions
        for V1, W, V3, V4 (notes, Task 10, "The four terms at general D"):
        the four structures V1, W, V3, V4 are differentiated by brute force in
        explicit coordinates with the actual A(y) of (I.62) and
        i Delta_A = A + k ln(aa'), and compared at three random non-coincident
        points with 30-digit precision.
    (e) General-D assembly -V1 + (D-1)/(D-3) W - V3 - V4 (units kappa^2 k) and its
        sum with the ln(aa')-type local terms of the derived (I.44) (Task 5):
        the sum is proportional to (D-4)(D+1)/(4(D-3)).
    (f) Paper I's printed (I.63)-(I.66) are -1 times the derived terms, and their
        sum reproduces the printed (I.67).

Method
    SymPy exact algebra; delta functions represented by their action on a test
    function; in (d) explicit coordinates and 30-digit numerical evaluation.

Expected result
    Every True/False check prints True; differences in (d) are 0.0; coefficients in (e):
    (D-1)/(D-3), -4/(D-3), (D-1)^2/(D-3), -2(D-4)/(D-3), 2/(D-3), 2/(D-3).

Runtime
    About 8 s (SymPy 1.14, Python 3.10+).

Run
    cd code/04_three_point_local_log && python3 check_60_67.py
"""
# Supplementary checks for Eqs. (60)-(67): logarithm contributions of the 3-point diagram.
import sympy as sp, random
I, R = sp.I, sp.Rational
D, H, eps = sp.symbols('D H epsilon', positive=True)
y = sp.Symbol('y', positive=True)

print("== (a) Eq. (62): A(y) at D = 4 + eps from (14) with c0 = Gamma(D/2-1) ==")
Dv = 4 + eps
N = H**(Dv - 2)/(4*sp.pi)**(Dv/2)
c0 = sp.gamma(Dv/2 - 1); c1 = sp.gamma(Dv/2 + 1)/(Dv/2 - 2)
first3 = N*(c0*(4/y)**(Dv/2 - 1) + c1*(4/y)**(Dv/2 - 2) - sp.pi*sp.cot(sp.pi*Dv/2)*sp.gamma(Dv - 1)/sp.gamma(Dv/2))
lim = sp.simplify(sp.series(first3, eps, 0, 1).removeO())
target = H**2/(16*sp.pi**2)*(4/y - 2*sp.log(y/4) - 1)
print("  first three terms -> H^2/(16 pi^2){4/y - 2 ln(y/4) - 1} :", sp.simplify(sp.expand_log(lim - target, force=True)) == 0)
n = sp.Symbol('n', positive=True, integer=True)
tn = (1/n)*sp.gamma(n + D - 1)/sp.gamma(n + D/2)*(y/4)**n - 1/(n - D/2 + 2)*sp.gamma(n + D/2 + 1)/sp.gamma(n + 2)*(y/4)**(n - D/2 + 2)
print("  n-th series term at D=4 vanishes identically             :", sp.simplify(sp.gammasimp(tn.subs(D, 4))) == 0)
c0p = sp.gamma(Dv/2 - 1)/(Dv/2 - 1)
print("  printed c0 changes A only at O(D-4)                      :", sp.limit(c0p - c0, eps, 0) == 0)

print("== (b) A(y) of (62) solves (4y-y^2)A'' + 4(2-y)A' = 3k with k = 2h  (Box A = (D-1) k H^2) ==")
h = sp.Symbol('h', positive=True)
A4 = h*(4/y - 2*sp.log(y/4) - 1)
print("  ", sp.simplify((4*y - y**2)*sp.diff(A4, y, 2) + 4*(2 - y)*sp.diff(A4, y) - 3*2*h) == 0)

print("== (c) delta identities on a test function g, general D (l = ln a, p = D-2) ==")
t = sp.Symbol('t', negative=True); x1, x2, x3 = xs = sp.symbols('x1 x2 x3', real=True)
a = -1/(H*t); l = sp.log(a); p = D - 2
g = sp.Function('g')(t, *xs)
lap = lambda f: sum(sp.diff(f, xi, 2) for xi in xs)
Op = lambda f: -sp.diff(a**p*sp.diff(f, t), t)          # (a^D Box - a^{D-2} Lap) acting
lhs = -2*I*sp.diff(a**p*l*sp.diff(g, t), t)              # -2i d0[a^p l d0 delta] acting on g
sym = l*Op(g) + Op(l*g)                                  # ln(aa') (a^D Box - a^{D-2}Lap) delta acting on g (symmetric convention)
rhs = I*sym + I*(D - 1)*H**2*a**D*g
print("  -2i d0[a^p l d0 delta] == i ln(aa')(a^D Box - a^{D-2}Lap) delta + i(D-1)H^2 a^D delta :", sp.simplify(lhs - rhs) == 0)
print("  ln(aa') Lap delta acting on g == 2 l Lap g (time-only factor)        :", sp.simplify((l*lap(g) + lap(l*g)) - 2*l*lap(g)) == 0)
print("  d0(a^p d0 ln a) = (D-1) H^2 a^D                                      :", sp.simplify(sp.diff(a**p*sp.diff(l, t), t) - (D - 1)*H**2*a**D) == 0)

print("== (d) pointwise check of the non-local parts, D = 4, actual A(y) of (62), k = 2h ==")
tp = sp.Symbol('tp', negative=True); ys = sp.symbols('y1 y2 y3', real=True)
ap = -1/(H*tp)
yy = a*ap*H**2*(sum((xs[i] - ys[i])**2 for i in range(3)) - (t - tp)**2)
k = 2*h
A = h*(4/yy - 2*sp.log(yy/4) - 1)
iDA = A + k*sp.log(a*ap)
Lg = sp.log(a*ap); pre = (a*ap)**2*k*Lg
lapx = lambda f: sum(sp.diff(f, xi, 2) for xi in xs)
V1 = sum(sp.diff(pre*sum(sp.diff(iDA, xs[j], ys[j]) for j in range(3)), xs[i], ys[i]) for i in range(3))
W = sp.diff(pre*sp.diff(iDA, t, tp), t, tp)
V3 = sum(sp.diff(pre*sp.diff(iDA, xs[i], tp), xs[i], tp) for i in range(3))
V4 = sum(sp.diff(pre*sp.diff(iDA, t, ys[i]), t, ys[i]) for i in range(3))
L2A = lapx(A); L4A = lapx(L2A)
R1 = k*(a*ap)**2*Lg*L4A
RW = k*((a*ap)**2*Lg*L4A + H*(a*ap)**2*(a*sp.diff(L2A, t) + ap*sp.diff(L2A, tp)))
R3 = k*((a*ap)**2*Lg*L4A + H*(a*ap)**2*ap*sp.diff(L2A, tp))
R4 = k*((a*ap)**2*Lg*L4A + H*(a*ap)**2*a*sp.diff(L2A, t))
random.seed(3)
for trial in range(3):
    pt = {H: sp.Rational(random.randint(5, 20), 10), h: sp.Rational(random.randint(1, 9), 7),
          t: -sp.Rational(random.randint(3, 9), 4), tp: -sp.Rational(random.randint(10, 19), 4)}
    for v in xs + ys: pt[v] = sp.Rational(random.randint(-9, 9), 3)
    res = [sp.N((e1 - e2).subs(pt), 30) for e1, e2 in ((V1, R1), (W, RW), (V3, R3), (V4, R4))]
    scale = [sp.N(abs(e2.subs(pt)), 30) for e2 in (R1, RW, R3, R4)]
    print(f"  point {trial}: |V1-R1|,|W-RW|,|V3-R3|,|V4-R4| =", [sp.Float(abs(r), 3) for r in res], " (scales", [sp.Float(s, 3) for s in scale], ")")

print("== (e) general-D assembly, units kappa^2 k: -V1 + (D-1)/(D-3) W - V3 - V4 ==")
LB, LL, M, Q4, T = sp.symbols('LB LL M Q4 T')   # L a^D Box i delta, L a^{D-2} Lap i delta, H^2 a^D i delta, (aa')^p L Lap^2 A, H(aa')^p a d0 Lap A (+ primed)
Tp = sp.Symbol('Tp')
V1e = Q4
We = (LB - 2*LL) + (D - 1)*M + Q4 + (T + Tp)
V3e = -LL + Q4 + Tp
V4e = -LL + Q4 + T
tot = sp.expand(-V1e + (D - 1)/(D - 3)*We - V3e - V4e)
print("  coefficients: LB:", sp.factor(tot.coeff(LB)), " LL:", sp.factor(tot.coeff(LL)), " M:", sp.factor(tot.coeff(M)),
      " Q4:", sp.factor(tot.coeff(Q4)), " T:", sp.factor(tot.coeff(T)), " Tp:", sp.factor(tot.coeff(Tp)))
four = -(D*(D - 1)/4*LB - D*LL + D*(D - 1)**2/4*M)            # log-type local terms of the derived (44), units kappa^2 k
loc3 = tot.coeff(LB)*LB + tot.coeff(LL)*LL + tot.coeff(M)*M
print("  derived (44)_log-local + derived (67)_local =", sp.factor(sp.expand(four + loc3)))
print("  at D=4: LB, LL, M coefficients of (67) =", [tot.coeff(s).subs(D, 4) for s in (LB, LL, M)], " T,Tp:", [tot.coeff(s).subs(D, 4) for s in (T, Tp)])

print("== (f) the paper's (63)-(66) and (67), units kappa^2 H^2/(8 pi^2) at D=4 ==")
p63 = Q4; p64 = 3*(-(LB - 2*LL) - 3*M - Q4 - (T + Tp)); p65 = -LL + Q4 + Tp; p66 = -LL + Q4 + T
print("  each equals minus the derived term:", sp.expand(p63 + (-V1e)) == 0, sp.expand(p64 + 3*We.subs(D, 4)) == 0,
      sp.expand(p65 + (-V3e)) == 0, sp.expand(p66 + (-V4e)) == 0)
psum = sp.expand(p63 + p64 + p65 + p66)
print("  paper sum:", psum, "   (paper's (67): -3LB + 4LL - 9M - 2T - 2Tp)")
