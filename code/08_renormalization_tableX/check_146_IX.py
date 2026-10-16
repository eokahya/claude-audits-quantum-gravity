"""check_146_IX.py -- supplementary symbolic checks for Paper I, Eqs. (I.146)-(I.159) and Tables VI-IX
(Task 16 of the re-derivation; chapter "Renormalization and the final self-mass: Tables VI-X").

Paper I: E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

What is checked (section letters as in the chapter):
  (A) pointwise parts of (I.148)-(I.150) for symbolic D; residues of y^(D/2-2) f; the delta coefficient
      of (I.147) is consistent with the propagator equation a^D Box iDelta_A = i delta^D.
  (B) (I.151) = (I.149) - 2/((D-3)(D-4)) x [(I.150) = 0]; exact delta coefficient; D -> 4 expansions
      leading to (I.152) and (I.153).
  (C) flat Euclidean analogue of the "add zero" step, exact in D with the test function exp(-r^2)
      (Gamma functions), and its simple pole at D = 4.
  (D) Table VI from Table IV (gamma3, eps3 eliminated), general D; the printed Table IV gamma2 sign
      would give the opposite gamma2 row.
  (E) Table VII at D = 4 from the nonlocal parts of (I.152), (I.153).
  (F) kernel identities (I.154), (I.155) and the six simple a' -> a conversions, general D, acting on
      an arbitrary test function phi(eta, x).
  (G) Tables VIII and IX (rows and totals); with the printed Table IV gamma2 sign Table IX would have
      a forbidden pole (value -16/3 at D = 4).
  (H) (I.156)-(I.159).
  (I) flat-space limit: the 3-point loop numerator V.P.V equals p^2 k^2/2 (D = 4, 5, 6, random rational
      momenta), so the flat-space divergence (and the a^D Box^2 total of Table VIII) must vanish.

Method: SymPy (exact rational/symbolic arithmetic).  All tables are in the paper's sign convention
(printed (I.47)-(I.49)); every check is linear and homogeneous, so it is independent of that sign.
Expected result: every line prints True / 0, the matches with the paper are True, and the
"printed gamma2 sign" diagnostics print the mismatching expressions (-16/3 at D = 4).
Runtime: about 3 s (Python 3.10+, SymPy 1.14).
"""
import random
import sympy as sp

D = sp.Symbol('D'); H = sp.Symbol('H', positive=True)
Y = sp.Symbol('y', positive=True); z = sp.Symbol('z', positive=True)   # z = 4/y
Lop = lambda f: sp.expand((4*Y - Y**2)*sp.diff(f, Y, 2) + D*(2 - Y)*sp.diff(f, Y))   # pointwise Box/H^2 on f(y)
Zp = lambda q: (4/Y)**q

print("== (A) pointwise parts of (148), (149), (150) and residues ==")
for name, q, rhs in (("(148)", D - 1, 2/(D - 2)**2*Lop(Zp(D - 2)) - 2/(D - 2)*Zp(D - 2)),
                     ("(149)", D - 2, 2/((D - 3)*(D - 4))*Lop(Zp(D - 3)) - 4/(D - 4)*Zp(D - 3))):
    print(f"  {name} pointwise:", sp.simplify(sp.powsimp(sp.expand(Zp(q) - rhs), force=True)) == 0)
print("  (150) pointwise: Box/H^2 (4/y)^(D/2-1) = D/2(D/2-1)(4/y)^(D/2-1):",
      sp.simplify(sp.powsimp(Lop(Zp(D/2 - 1)) - D/2*(D/2 - 1)*Zp(D/2 - 1), force=True)) == 0)
res = lambda f: sp.simplify(sp.powsimp(sp.expand(Y**(D/2 - 2)*f*Y), force=True))   # y * y^(D/2-2) f : constant iff residue
print("  y^(D/2-2) f for f=(4/y)^(D/2-1):", res(Zp(D/2 - 1)), "  -> residue 4^(D/2-1)")
print("  y^(D/2-2) f for f=(4/y)^(D-2):", res(Zp(D - 2)), "  for f=(4/y)^(D-3):", res(Zp(D - 3)), "  (y-dependent for generic D -> no residue)")
C = (4*sp.pi)**(D/2)*H**(-D)/sp.gamma(D/2 - 1)
print("  (147) delta coefficient for (4/y)^(D/2-1), divided by H^2:",
      sp.simplify(sp.gammasimp(4**(D/2 - 1)*4*sp.pi**(D/2)*H**(2 - D)/sp.gamma(D/2 - 1)/H**2 - C)) == 0)
N = H**(D - 2)/(4*sp.pi)**(D/2)
print("  consistency with a^D Box A_lead = i delta, A_lead = N Gamma(D/2-1)(4/y)^(D/2-1):",
      sp.simplify(sp.gammasimp(N*sp.gamma(D/2 - 1)*C*H**2 - 1)) == 0)

print("== (B) (151) from (149) and (150), then (152) and (153) as D -> 4 ==")
# (151) - (149) must equal -2/((D-3)(D-4)) * [(150) which is zero]: verify coefficient bookkeeping symbolically with formal objects
BoxA, BoxB, delta = sp.symbols('BoxA BoxB delta')   # BoxA = Box/H^2 (4/y)^(D-3), BoxB = Box/H^2 (4/y)^(D/2-1)
fA, fB = sp.symbols('fA fB')                         # fA = (4/y)^(D-3), fB = (4/y)^(D/2-1)
r149 = 2/((D - 3)*(D - 4))*BoxA - 4/(D - 4)*fA
zero150 = BoxB - D/2*(D/2 - 1)*fB - C*sp.I*delta
r151 = 2/((D - 3)*(D - 4))*(C*sp.I*delta + BoxA - BoxB) - 4/(D - 4)*(fA - D*(D - 2)/(8*(D - 3))*fB)
print("  (151) = (149) - 2/((D-3)(D-4)) x [(150)=0]:", sp.simplify(sp.expand(r151 - (r149 - 2/((D - 3)*(D - 4))*zero150))) == 0)
print("  delta coefficient 2C/((D-3)(D-4)) == (4pi)^(D/2) H^-D (D-2)/((D-3)(D-4)Gamma(D/2)):",
      sp.simplify(sp.gammasimp(2*C - (4*sp.pi)**(D/2)*H**(-D)*(D - 2)/sp.gamma(D/2))) == 0)
e = sp.Symbol('e'); Dd = 4 + e
b1 = sp.series(2/((Dd - 3)*e)*(z**(Dd - 3) - z**(Dd/2 - 1)), e, 0, 1).removeO()
b2 = sp.series(-4/e*(z**(Dd - 3) - Dd*(Dd - 2)/(8*(Dd - 3))*z**(Dd/2 - 1)), e, 0, 1).removeO()
print("  2/((D-3)(D-4)) [(4/y)^(D-3)-(4/y)^(D/2-1)] -> ", sp.simplify(sp.expand_log(b1, force=True)), "  (= -(4/y)ln(y/4))")
print("  -4/(D-4){...} ->", sp.simplify(sp.expand_log(b2, force=True)), "  (= 2(4/y)ln(y/4) - 4/y with ln z = -ln(y/4))")
# (153) from (148) with (152)
Lloc = sp.Symbol('Lloc'); NL = sp.Symbol('NL')     # (4/y)^(D-2) = Lloc + NL
BoxNL = sp.Symbol('BoxNL')
r153 = 2/(D - 2)**2*(sp.Symbol('BoxLloc') + BoxNL) - 2/(D - 2)*(Lloc + NL)
print("  (153): 2/(D-2)^2 Box/H^2 [(152)] - 2/(D-2)[(152)];  at D=4 coefficients of Box(NL), NL:", r153.subs(D, 4).coeff(BoxNL), r153.subs(D, 4).coeff(NL))
print("  local: 2/(D-2)^2 * (D-2) = 2/(D-2) (Box/H^2 term) and -2/(D-2)*(D-2) = -2  ->  {2/(D-2) Box/H^2 - 2}  OK")

print("== (C) flat Euclidean analogue of the 'add zero' step (151), exact in D (radial test function exp(-r^2)) ==")
# int_0^inf r^s e^{-r^2} dr = Gamma((s+1)/2)/2 ;  Lap e^{-r^2} = (4r^2 - 2D) e^{-r^2} ;  Lap r^{2-D} = -(D-2) Omega_D delta^D
Om = 2*sp.pi**(D/2)/sp.gamma(D/2)
I_ = lambda s_: sp.gamma((s_ + 1)/2)/2
lhsC = Om*I_(3 - D)                                     # <r^{-2(D-2)}, phi>
lapint = lambda pw: 4*I_(D - 1 + pw + 2) - 2*D*I_(D - 1 + pw)   # int r^{D-1} r^pw Lap(phi)
rhsC = (Om*(lapint(6 - 2*D) - lapint(2 - D)) - (D - 2)*Om*1)/(2*(D - 3)*(D - 4))
print("  <r^-2(D-2),phi> == [<r^-2(D-3) - r^(2-D), Lap phi> - (D-2) Omega phi(0)] / (2(D-3)(D-4))  for all D:",
      sp.simplify(sp.gammasimp(lhsC - rhsC)) == 0)
print("  pole structure: <r^-2(D-2),phi> =", sp.series(sp.gammasimp(lhsC), D, 4, 1).removeO().simplify(), " (simple pole at D=4)")

print("== (D) Table VI from Table IV via (152),(153), with gamma3 = 2 gamma1 + gamma2, eps3 = 2 eps1 + eps2 eliminated ==")
q = D**3 - 3*D**2 - 4*D + 8
T4 = {'alpha': (0, D/((D - 1)*(D - 2)**2)), 'beta': (-D/(4*(D - 1)), -q/(4*(D - 1)*(D - 2))),
      'gamma1': (-D*(D - 2)/4, -q/4), 'gamma2': (D/4, +q/(4*(D - 1))), 'gamma3': (-D/4, 0),
      'delta': (0, 0), 'eps1': (0, (D**2 - 6*D + 4)/(2*(D - 1)*(D - 2))), 'eps2': (0, (1 - 2*D)/((D - 1)*(D - 2))), 'eps3': (0, 0), 'zeta': (0, 0)}
eff = dict(T4)
eff['gamma1'] = tuple(sp.simplify(T4['gamma1'][i] + 2*T4['gamma3'][i]) for i in range(2))
eff['gamma2'] = tuple(sp.simplify(T4['gamma2'][i] + T4['gamma3'][i]) for i in range(2))
eff['eps1'] = tuple(sp.simplify(T4['eps1'][i] + 2*T4['eps3'][i]) for i in range(2))
eff['eps2'] = tuple(sp.simplify(T4['eps2'][i] + T4['eps3'][i]) for i in range(2))
B = sp.Symbol('Box_H2')   # Box/H^2 acting on delta/a^D
T6 = {k: sp.factor(sp.expand(c1*(2/(D - 2)*B - 2) + c2*(D - 2))) for k, (c1, c2) in eff.items() if k not in ('gamma3', 'eps3')}
paper6 = {'alpha': D/((D - 1)*(D - 2)), 'beta': -D/(2*(D - 1)*(D - 2))*B - (D + 2)*(D - 4)/4,
          'gamma1': -D**2/(2*(D - 2))*B - (D**4 - 5*D**3 + 16*D - 16)/4, 'gamma2': (D - 2)*q/(4*(D - 1)),
          'delta': 0, 'eps1': (D**2 - 6*D + 4)/(2*(D - 1)), 'eps2': (1 - 2*D)/(D - 1), 'zeta': 0}
for k in T6:
    print(f"  {k:7s}: derived = {T6[k]!s:60s} matches paper: {sp.simplify(T6[k] - paper6[k]) == 0}")
g2_printed = sp.factor(-q/(4*(D - 1))*(D - 2))
print("  (with Table IV's printed gamma2 sign the gamma2 row would be", g2_printed, " -> paper's Table VI uses the corrected sign)")

print("== (E) Table VII at D=4 from the nonlocal parts of (152),(153) ==")
x = sp.Symbol('x', positive=True); lx = sp.log(x)
B1, B2 = sp.symbols('B1 B2')   # B1 = Box/H^2 [..], B2 = Box^2/H^4 [..]  acting on the bracket that follows
# nonlocal (152): -B1[lnx/x] + 2 lnx/x - 1/x ; (153): -1/2 B2[lnx/x] + B1[2 lnx/x - 1/(2x)] - 2 lnx/x + 1/x
def nl(c1, c2):
    c1, c2 = sp.nsimplify(c1), sp.nsimplify(c2)
    return {'B2': sp.simplify(-c1/2*lx/x), 'B1': sp.simplify(c1*(2*lx/x - 1/(2*x)) - c2*lx/x),
            '1': sp.simplify(c1*(-2*lx/x + 1/x) + c2*(2*lx/x - 1/x))}
paper7 = {'alpha': {'B2': 0, 'B1': -lx/(3*x), '1': 2*lx/(3*x) - 1/(3*x)},
          'beta': {'B2': lx/(6*x), 'B1': -lx/(3*x) + 1/(6*x), '1': 0},
          'gamma1': {'B2': 2*lx/x, 'B1': -6*lx/x + 2/x, '1': 4*lx/x - 2/x},
          'gamma2': {'B2': 0, 'B1': -2*lx/(3*x), '1': 4*lx/(3*x) - 2/(3*x)},
          'delta': {'B2': 0, 'B1': 0, '1': 0}, 'eps1': {'B2': 0, 'B1': lx/(3*x), '1': -2*lx/(3*x) + 1/(3*x)},
          'eps2': {'B2': 0, 'B1': 7*lx/(6*x), '1': -7*lx/(3*x) + 7/(6*x)}, 'zeta': {'B2': 0, 'B1': 0, '1': 0}}
for k in paper7:
    c1, c2 = [sp.simplify(sp.sympify(v).subs(D, 4)) for v in eff[k]]
    d = nl(c1, c2)
    ok = all(sp.simplify(d[s] - paper7[k][s]) == 0 for s in d)
    print(f"  {k:7s}: (c1,c2)|_D=4 = ({c1},{c2})   matches paper Table VII: {ok}")

print("== (F) kernel identities (154), (155) and the other conversions, general D ==")
t = sp.Symbol('t', negative=True); xs = sp.symbols('x1 x2 x3', real=True)
a = -1/(H*t); phi = sp.Function('phi')(t, *xs)
lap = lambda f: sum(sp.diff(f, xi, 2) for xi in xs)
box = lambda f: (-sp.diff(a**(D - 2)*sp.diff(f, t), t) + a**(D - 2)*lap(f))/a**D
# external operator a^p a'^q [Box^n] acting on a^{-D} delta(x-x'), applied to a test function phi(x'):
#   int dx' a(x)^p a(x')^q Box_x^n [a(x)^{-D} delta(x-x')] phi(x') = a^p Box^n [a^{-D} a^q phi]
act = lambda p, qq, n: a**p*(lambda f: [f := box(f) for _ in range(n)][-1] if n else f)(a**(-D)*a**qq*phi)
lhs154 = act(D - 1 + 2, D - 1, 2) + act(D - 1, D + 1, 2)
rhs154 = 2*a**D*box(box(phi)) - 12*H**2*a**D*box(phi) + 8*a**(D - 2)*H**2*lap(phi) + 2*(D**2 - 2*D + 2)*H**4*a**D*phi
print("  (154):", sp.simplify(sp.expand(lhs154 - rhs154)) == 0)
lhs155 = H**2*(act(D - 1 + 2, D - 1, 1) + act(D - 1, D + 1, 1))
print("  (155):", sp.simplify(sp.expand(lhs155 - 2*a**D*(H**2*box(phi) - H**4*phi))) == 0)
print("  alpha = (aa')^D Box^2      -> a^D Box^2 delta     :", sp.simplify(act(D, D, 2) - a**D*box(box(phi))) == 0)
print("  gamma1 Box/H^2 = (aa')^D H^2 Box -> a^D H^2 Box   :", sp.simplify(H**2*act(D, D, 1) - H**2*a**D*box(phi)) == 0)
print("  gamma1 = (aa')^D H^4        -> a^D H^4             :", sp.simplify(H**4*act(D, D, 0) - H**4*a**D*phi) == 0)
print("  gamma2 = (aa')^{D-1}(a^2+a'^2)H^4 -> 2 a^D H^4     :", sp.simplify(H**4*(act(D + 1, D - 1, 0) + act(D - 1, D + 1, 0)) - 2*H**4*a**D*phi) == 0)
lapact = lambda p, qq: a**p*lap(a**(-D)*a**qq*phi)
print("  eps1 = (aa')^{D-1}H^2 Lap   -> a^{D-2} H^2 Lap      :", sp.simplify(H**2*lapact(D - 1, D - 1) - H**2*a**(D - 2)*lap(phi)) == 0)
print("  eps2 = (aa')^{D-2}(a^2+a'^2)H^2 Lap -> 2 a^{D-2}H^2 Lap:", sp.simplify(H**2*(lapact(D, D - 2) + lapact(D - 2, D)) - 2*H**2*a**(D - 2)*lap(phi)) == 0)

print("== (G) Tables VIII and IX from Table VI; totals; consistency checks ==")
# columns: a^D Box^2, a^D H^2 Box, a^{D-2} H^2 Lap, a^D H^4
conv = {   # (operator in Table VI) -> contribution to the four columns per unit coefficient
 'alpha': [1, 0, 0, 0], 'betaB': [2, -12, 8, 2*(D**2 - 2*D + 2)], 'beta': [0, 2, 0, -2], 'gamma1B': [0, 1, 0, 0],
 'gamma1': [0, 0, 0, 1], 'gamma2': [0, 0, 0, 2], 'eps1': [0, 0, 1, 0], 'eps2': [0, 0, 2, 0]}
E6 = {k: sp.expand(v) for k, v in T6.items()}
coef = {'alpha': T6['alpha'], 'betaB': E6['beta'].coeff(B), 'beta': E6['beta'].subs(B, 0), 'gamma1B': E6['gamma1'].coeff(B),
        'gamma1': E6['gamma1'].subs(B, 0), 'gamma2': T6['gamma2'], 'eps1': T6['eps1'], 'eps2': T6['eps2']}
tot = [sp.factor(sum(coef[k]*conv[k][j] for k in conv)) for j in range(4)]
for k in conv:
    print(f"  {k:8s}:", [sp.factor(coef[k]*conv[k][j]) for j in range(4)])
print("  totals:", tot)
paperVIII = [0, (D - 4)*(-D**3 + D - 4)/(2*(D - 1)*(D - 2)), (D**3 - 16*D**2 + 28*D - 16)/(2*(D - 1)*(D - 2))]
paperIX = (D - 4)*(-D**4 + 5*D**3 - 16*D + 16)/(4*(D - 2))
print("  Table VIII totals match paper:", [sp.simplify(tot[j] - paperVIII[j]) == 0 for j in range(3)])
print("  Table IX total matches paper:", sp.simplify(tot[3] - paperIX) == 0, "  (has factor (D-4): finite)")
g2bad = sp.factor(tot[3] - 4*coef['gamma2'])   # flip the sign of the gamma2 contribution (2*coef -> -2*coef)
print("  with Table IV's printed gamma2 sign, Table IX total =", g2bad, " -> value at D=4:", sp.simplify(g2bad.subs(D, 4)), "(nonzero: forbidden H^4 divergence)")

print("== (H) (156)-(159) ==")
print("  Box ln a = -(D-1)H^2  =>  a^D Box[k ln(aa')] = -(D-1) k H^2 a^D :", sp.simplify(box(sp.log(a)) + (D - 1)*H**2) == 0)
v157 = sp.factor(sp.simplify(paperIX/((D - 3)*(D - 4))))
print("  (157) coefficient (x Gamma(D/2) i kappa^2 H^{D-4}/(4pi)^{D/2}):", v157)
lim = sp.limit(v157*sp.gamma(D/2)*H**(D - 4)/(4*sp.pi)**(D/2), D, 4)
print("  (158) D->4:", sp.simplify(lim), " = 1/(8 pi^2)")
kk = H**2/(8*sp.pi**2); h = H**2/(16*sp.pi**2)
yy = sp.Symbol('yy', positive=True)
A4 = h*(4/yy - 2*sp.log(yy/4) - 1)
# coefficient of (aa')^4 Box in units kappa^2: H^4/(8 pi^2) * A(y); compare with H^6/(4pi)^4 [2*4/y - 4 ln(y/4)]
print("  (159) Box part: H^4/(8pi^2) A(y) - H^6/(4pi)^4 [2*4/y - 4 ln(y/4)] is y-independent:",
      sp.simplify(sp.diff(H**4/(8*sp.pi**2)*A4 - H**6/(4*sp.pi)**4*(2*4/yy - 4*sp.log(yy/4)), yy)) == 0)
print("  (159) H^4 part: -H^4/(8pi^2)*3kH^2 == -kappa^2H^4/(4pi)^4 * 12 H^4 :", sp.simplify(-H**4/(8*sp.pi**2)*3*kk*H**2 + H**4/(4*sp.pi)**4*12*H**4) == 0)

print("== (I) flat-space limit: numerator of the 3-point loop, V^{ab}(p,k) P_{abrs} V^{rs}(k,p) = (1/2) p^2 k^2 ==")
random.seed(3)
for Dv in (4, 5, 6):
    eta = sp.diag(*([-1] + [1]*(Dv - 1)))
    p = sp.Matrix([sp.Rational(random.randint(-9, 9), 7) for _ in range(Dv)]); kv = sp.Matrix([sp.Rational(random.randint(-9, 9), 5) for _ in range(Dv)])
    dot = lambda u, v: (u.T*eta*v)[0]
    V = sp.Matrix(Dv, Dv, lambda a_, b_: sp.Rational(1, 2)*(p[a_]*kv[b_] + p[b_]*kv[a_]) - sp.Rational(1, 2)*eta.inv()[a_, b_]*dot(p, kv))   # upper indices
    Vl = eta*V*eta                   # lower indices
    VV = sum(V[a_, b_]*Vl[a_, b_] for a_ in range(Dv) for b_ in range(Dv))
    trV = sum(V[a_, b_]*eta[a_, b_] for a_ in range(Dv) for b_ in range(Dv))
    Nn = VV - trV**2/(Dv - 2)        # P = 1/2(eta eta + eta eta) - 1/(D-2) eta eta
    print(f"  D={Dv}: N - p^2 k^2/2 =", sp.simplify(Nn - dot(p, p)*dot(kv, kv)/2))
print("  => Sigma_flat ~ (p^2/2) int d^Dk 1/(p-k)^2 = 0 in dimensional regularization (scaleless): the a^D Box^2 divergence must vanish")
