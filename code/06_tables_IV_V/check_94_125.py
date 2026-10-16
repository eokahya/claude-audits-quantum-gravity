"""check_94_125.py -- identities and expansions of Paper I, Eqs. (I.94)-(I.125)  [Task 13, checks (A)-(D)]

Purpose
  Independent symbolic/numerical checks of the analytic derivations of Chapter 6
  ("Tables IV and V: Eqs. (I.94)-(I.125)") of the derivation notes.

What is checked
  (A) Eqs. (I.94), (I.95), (I.96), (I.97), (I.98), (I.100), (I.101) and (I.102)-(I.103):
      pointwise in explicit conformal coordinates for D = 4 and D = 5, with the
      generic test functions A = y^beta, B = y^gamma, f = y^alpha (beta=-5/3,
      gamma=2/9, alpha=3/7), at one random rational point per D, 30-digit arithmetic.
      Box is applied from its definition (I.42); I^n are sympy antiderivatives.
  (B) Eq. (I.99) as an operator identity; Eqs. (I.108), (I.112) (gamma3 = 2 gamma1 + gamma2,
      eps3 = 2 eps1 + eps2) for general D.
  (C) Eqs. (I.116)-(I.123): exact-D coefficients of (4/y)^D and (4/y)^(D-1) of the eight
      products (I.115), computed from the two leading terms of A, B, C (the only terms
      that reach these powers), in units Gamma^2(D/2) N^2/16, N = H^(D-2)/(4 pi)^(D/2).
      Then the D = 4 closed forms: every series term of (I.27), (I.28) vanishes at D = 4
      (so B = C = 4h/y), and the exact D = 4 products equal (h^2/16) x the paper's brackets.
  (D) Eqs. (I.124)-(I.125): the D = 4 bracket of (A + Delta C) A''; the double integral of its
      finite part reproduces the finite bracket of (I.125) up to a constant; the exact-D
      leading term gives D/((D-1)(D-2)^2) (4/y)^(D-2) (checked at D = 5, 6, 7, 9/2).

Expected result
  All (A), (B) entries True; (C) the coefficients of the table of leading coefficients
  in the notes (e.g. A'dC': 0 and -8; poles 1/(D-4) for AA'', A''dB, A''dC at (4/y)^(D-1)),
  and every D = 4 closed-form line True; (D) constant (24 ln 2 - 36) ln 2 and
  [True, True, True, True].

Runtime: about 1.5 minutes (dominated by part (A)).
Run:     python3 check_94_125.py
"""
# Supplementary checks for Eqs. (94)-(125) and Tables IV-V.
import sympy as sp, random, itertools
R = sp.Rational
H = sp.Symbol('H', positive=True)
al, be, ga = sp.symbols('alpha beta gamma')
Yv = sp.Symbol('Y', positive=True)
def Iop(g, n=1):
    for _ in range(n): g = sp.integrate(sp.expand(g), Yv, conds='none')
    return g
def setup(Dv):
    t, tp = sp.symbols('t tp', negative=True)
    xs = sp.symbols(f'x1:{Dv}', real=True); ys = sp.symbols(f'y1:{Dv}', real=True)
    a, ap = -1/(H*t), -1/(H*tp)
    y = a*ap*H**2*(sum((xs[i] - ys[i])**2 for i in range(Dv - 1)) - (t - tp)**2)
    return t, tp, xs, ys, a, ap, y
random.seed(17)
def rnd(t, tp, xs, ys, spacelike=False):
    while True:
        pt = {H: random.choice([R(7, 10), R(13, 10), R(3, 2)]), t: -R(random.randint(3, 9), 4), tp: -R(random.randint(10, 19), 4),
              al: R(3, 7), be: -R(5, 3), ga: R(2, 9)}
        for v in xs + ys: pt[v] = R(random.randint(-9, 9), 3)
        r2 = sum((pt[xs[i]] - pt[ys[i]])**2 for i in range(len(xs))); de2 = (pt[t] - pt[tp])**2
        if not spacelike or r2 > de2 + 1: return pt
num = lambda e, pt: sp.N(e.subs(pt), 30)
F = lambda e, y: e.subs(Yv, y)

print("== (A) Eqs. (94)-(96), (97)-(98), (100)-(101), (102)-(103): pointwise, D=4,5 ==")
for Dv in (4, 5):
    t, tp, xs, ys, a, ap, y = setup(Dv); p = Dv - 2
    lap = lambda e: sum(sp.diff(e, xi, 2) for xi in xs); lapp = lambda e: sum(sp.diff(e, yi, 2) for yi in ys)
    box = lambda e: (-sp.diff(a**p*sp.diff(e, t), t) + a**p*lap(e))/a**Dv
    boxp = lambda e: (-sp.diff(ap**p*sp.diff(e, tp), tp) + ap**p*lapp(e))/ap**Dv
    A = Yv**be; B = Yv**ga; f = Yv**al
    Ay, By, fy = F(A, y), F(B, y), F(f, y)
    G = F(Iop(A*sp.diff(A, Yv, 2), 2), y); K = F(Iop(sp.diff(A, Yv)**2, 2), y); K3 = F(Iop(sp.diff(A, Yv)**2, 3), y)
    lhs94 = sum(sp.diff((a*ap)**p*Ay*sp.diff(Ay, xs[i], tp), xs[i], tp) for i in range(Dv - 1))
    r94 = lap(sp.diff((a*ap)**p*sp.diff(G, tp), tp)) + H*lap(sp.diff(a**p*ap**(p + 1)*K, tp))
    r95 = (-a**p*ap**Dv*lap(boxp(G)) + (a*ap)**p*lap(lap(G)) + H*a**p*ap**(p + 1)*lap(sp.diff(K, tp)) + (Dv - 1)*H**2*a**p*ap**Dv*lap(K))
    r96 = (-a**p*ap**Dv*lap(box(G)) + (a*ap)**p*lap(lap(G)) - H*a**(p - 1)*ap**Dv*lap(sp.diff(K, t)) + R(1, 2)*a**(p - 1)*ap**(p + 1)*lap(lap(K3)))
    e97 = sp.diff((a*ap)**p*sp.diff(fy, t), t) - (-a**Dv*ap**p*box(fy) + (a*ap)**p*lap(fy))
    e98 = sp.diff((a*ap)**p*sp.diff(fy, tp), tp) - (-a**p*ap**Dv*box(fy) + (a*ap)**p*lap(fy))
    lhs100 = sp.diff((a*ap)**p*fy, t, tp)
    r100 = (a*ap)**p*sp.diff(fy, t, tp) + p*H*(a*ap)**p*(ap*sp.diff(fy, t) + a*sp.diff(fy, tp)) + p**2*H**2*(a*ap)**(p + 1)*fy
    fp, fpp = F(sp.diff(f, Yv), y), F(sp.diff(f, Yv, 2), y); If = F(Iop(f), y)
    r101 = (R(1, 2)*(a*ap)**p*(a**2 + ap**2)*(box(fy) + H**2*y*fp) - (a*ap)**p*lap(fy) + R(1, 4)*p*(a*ap)**p*lap(If)
            + R(1, 2)*p*(Dv - 3)*H**2*(a*ap)**(p + 1)*fy + R(1, 2)*H**2*(a + ap)**2*(a*ap)**p*((Dv - 1)*y*fp + y**2*fpp))
    lhs102 = sum(sp.diff((a*ap)**p*By*sp.diff(Ay, t, ys[i]), xs[i], tp) + sp.diff((a*ap)**p*By*sp.diff(Ay, xs[i], tp), t, ys[i]) for i in range(Dv - 1))
    I2AB = F(Iop(sp.diff(A, Yv, 2)*B, 2), y); IAByAB = F(Iop(sp.diff(A, Yv)*B + Yv*sp.diff(A, Yv, 2)*B, 1), y)
    I3ApBp = F(Iop(sp.diff(A, Yv)*sp.diff(B, Yv), 3), y); I1ApBp = F(Iop(sp.diff(A, Yv)*sp.diff(B, Yv), 1), y); I2ApBp = F(Iop(sp.diff(A, Yv)*sp.diff(B, Yv), 2), y)
    Apy, Appy = F(sp.diff(A, Yv), y), F(sp.diff(A, Yv, 2), y)
    r103 = (-(a**2 + ap**2)*(a*ap)**p*lap(box(I2AB) + H**2*IAByAB) + (a*ap)**p*lap(lap(2*I2AB - R(1, 2)*I3ApBp))
            - H**2*(a + ap)**2*(a*ap)**p*lap((Dv - 2)*IAByAB + y*Apy*By + y**2*Appy*By - y*I1ApBp)
            + (Dv - 1)*H**2*(a**2 + a*ap + ap**2)*(a*ap)**p*lap(I2ApBp))
    pt = rnd(t, tp, xs, ys); s94 = max(1, abs(num(lhs94, pt))); s100 = max(1, abs(num(lhs100, pt))); s102 = max(1, abs(num(lhs102, pt)))
    print(f"  D={Dv}: (94) {abs(num(lhs94 - r94, pt)) < 1e-20*s94}  (95) {abs(num(lhs94 - r95, pt)) < 1e-20*s94}  (96) {abs(num(lhs94 - r96, pt)) < 1e-20*s94}"
          f"  (97) {abs(num(e97, pt)) < 1e-20}  (98) {abs(num(e98, pt)) < 1e-20}  (100) {abs(num(lhs100 - r100, pt)) < 1e-20*s100}"
          f"  (101) {abs(num(lhs100 - r101, pt)) < 1e-20*s100}  (102)-(103) {abs(num(lhs102 - r103, pt)) < 1e-20*s102}")

print("== (B) (99), (108), (112) ==")
u, v = sp.symbols('u v'); g = sp.Function('g')(u, v)
print("  (99):", sp.simplify(sp.diff(g, u, v) - (R(1, 2)*(sp.diff(g, u, 2) + 2*sp.diff(g, u, v) + sp.diff(g, v, 2)) - R(1, 2)*sp.diff(g, u, 2) - R(1, 2)*sp.diff(g, v, 2))) == 0)
aa, bb, Dd = sp.symbols('a ap D', positive=True)
print("  (108) gamma3 = 2 gamma1 + gamma2 :", sp.simplify((aa*bb)**(Dd - 1)*(aa + bb)**2 - (2*(aa*bb)**Dd + (aa*bb)**(Dd - 1)*(aa**2 + bb**2))) == 0)
print("  (112) eps3 = 2 eps1 + eps2       :", sp.simplify((aa*bb)**(Dd - 2)*(aa + bb)**2 - (2*(aa*bb)**(Dd - 1) + (aa*bb)**(Dd - 2)*(aa**2 + bb**2))) == 0)

print("== (C) Eqs. (116)-(123): exact-D leading coefficients (two leading terms of A, B, C) and D=4 closed forms ==")
D = sp.Symbol('D'); Ns = sp.Symbol('N')
z = sp.Symbol('z', positive=True)              # z = 4/y
pp = D/2 - 1                                    # exponent of the leading term
c0 = sp.gamma(D/2 - 1); a1 = sp.gamma(D/2 + 1)/(D/2 - 2); b1 = sp.gamma(D/2); cc1 = -(3 - D/2)*sp.gamma(D/2 - 1)
A_ = Ns*(c0*z**pp + a1*z**(pp - 1)); B_ = Ns*(c0*z**pp + b1*z**(pp - 1)); C_ = Ns*(c0*z**pp + cc1*z**(pp - 1))
dB_, dC_ = B_ - A_, 2*(D - 2)/(D - 3)*(C_ - A_)
dy = lambda e: sp.diff(e, z)*(-z**2/4)          # d/dy = -(z^2/4) d/dz
unit = Ns**2*sp.gamma(D/2)**2/16
def cof(expr, q):
    e = sp.expand(sp.powsimp(sp.expand(expr/unit), force=True))
    tot = sp.Integer(0)
    for t in sp.Add.make_args(e):
        r = sp.simplify(sp.powsimp(t/z**q, force=True))
        if r.free_symbols.isdisjoint({z}): tot += r
    return sp.factor(sp.simplify(sp.gammasimp(tot)))
prods = {"A'^2": dy(A_)**2, "AA''": None, "A'B'": dy(A_)*dy(B_), "A''B": dy(dy(A_))*B_,
         "A'dB'": dy(A_)*dy(dB_), "A''dB": dy(dy(A_))*dB_, "A'dC'": dy(A_)*dy(dC_), "A''dC": dy(dy(A_))*dC_}
prods["AA''"] = A_*dy(dy(A_))
for name, Pr in prods.items():
    print(f"  {name:6s}: coeff (4/y)^D = {cof(Pr, D)},   coeff (4/y)^(D-1) = {cof(Pr, D - 1)}")
print("  (A-type products' (4/y)^(D-1) coefficient is not separately meaningful when it carries a 1/(D-4) pole: it combines with the")
print("   constant -pi cot(pi D/2) term into the D=4 logarithms; only the (4/y)^D coefficient of A''-type products enters Table IV.)")
# D=4 closed forms
h = sp.Symbol('h', positive=True)
A4 = h*(z - 2*sp.log(1/z) - 1); B4 = h*z; C4 = h*z
dB4, dC4 = B4 - A4, 4*(C4 - A4)
L = sp.log(1/z)                                   # ln(y/4)
paper = {"A'^2": z**4 + 4*z**3 + 4*z**2, "AA''": 2*z**4 - 4*z**3*L - 4*z**2*L - 2*z**2, "A'B'": z**4 + 2*z**3, "A''B": 2*z**4 + 2*z**3,
         "A'dB'": -2*z**3 - 4*z**2, "A''dB": 4*z**3*L + 2*z**3 + 4*z**2*L + 2*z**2, "A'dC'": -8*z**3 - 16*z**2, "A''dC": 16*z**3*L + 8*z**3 + 16*z**2*L + 8*z**2}
P4 = {"A'^2": dy(A4)**2, "AA''": A4*dy(dy(A4)), "A'B'": dy(A4)*dy(B4), "A''B": dy(dy(A4))*B4,
      "A'dB'": dy(A4)*dy(dB4), "A''dB": dy(dy(A4))*dB4, "A'dC'": dy(A4)*dy(dC4), "A''dC": dy(dy(A4))*dC4}
nn = sp.Symbol('n', integer=True, nonnegative=True); yy4 = sp.Symbol('u', positive=True)   # u = y/4
tB = sp.gamma(nn + D - 2)/sp.gamma(nn + D/2)*yy4**nn - sp.gamma(nn + D/2)/sp.gamma(nn + 2)*yy4**(nn - D/2 + 2)
tC = (nn + 1)*sp.gamma(nn + D - 3)/sp.gamma(nn + D/2)*yy4**nn - (nn - D/2 + 3)*sp.gamma(nn + D/2 - 1)/sp.gamma(nn + 2)*yy4**(nn - D/2 + 2)
print("  D=4: every series term of (27) and (28) vanishes, so B = C = i Delta_cf = 4h/y :",
      sp.simplify(sp.gammasimp(tB.subs(D, 4))) == 0, sp.simplify(sp.gammasimp(tC.subs(D, 4))) == 0)
for name in paper:
    print(f"  {name:6s}: exact D=4 product == (h^2/16) x paper's bracket at D=4 :", sp.simplify(sp.expand(P4[name] - h**2/16*paper[name])) == 0)

print("== (D) Eq. (124)-(125): alpha coefficient I^2[(A+dC)A''] ==")
br = sp.expand(paper["AA''"] + paper["A''dC"]); print("  (A+dC)A'' bracket at D=4:", br)
fin = br - 2*z**4                                               # remove the D-exact leading term
yy = sp.Symbol('y', positive=True)
I2fin = Iop(sp.expand(fin.subs(z, 4/Yv)), 2).subs(Yv, yy)
target = 6*(4/yy)*sp.log(yy/4) + 13*(4/yy) - 6*sp.log(yy/4)**2 - 18*sp.log(yy/4)
diff = sp.expand(sp.expand_log(I2fin/16 - target, force=True))
print("  I^2[finite part]/16 - (paper's finite bracket) =", sp.nsimplify(sp.expand_log(sp.simplify(diff), force=True)), "  (a constant: killed by alpha = (aa')^D Box^2)")
okl = []
for Dv in (5, 6, 7, sp.Rational(9, 2)):
    I2 = sp.integrate(sp.integrate(sp.Rational(1)*Dv/(Dv - 2)*(4/yy)**Dv, yy), yy)/16
    okl.append(sp.simplify(I2 - Dv/((Dv - 1)*(Dv - 2)**2)*(4/yy)**(Dv - 2)) == 0)
print("  I^2[D/(D-2)(4/y)^D]/16 == D/((D-1)(D-2)^2) (4/y)^(D-2) for D = 5, 6, 7, 9/2 :", okl)

