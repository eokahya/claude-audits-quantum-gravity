"""check_160_167.py -- symbolic checks for Paper I, Eqs. (I.160)-(I.167) and Table X
(Task 17 of the re-derivation; chapter "Renormalization and the final self-mass: Tables VI-X").

Paper I: E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

What is checked:
  (A) Table X assembled row by row from Table III (Task 12), Table V (as printed), Table VII (Task 16)
      and Eq. (I.159) (Task 16), after eliminating gamma3 = 2 gamma1 + gamma2 and eps3 = 2 eps1 + eps2.
      Prints "assembled - paper" for each of the 13 operators, and the Table III contribution to the
      (aa')^2 Lap^2 (zeta) row that the printed table does not contain (8 ln x - 16 x ln x).
      NOTE: the assembled zeta row deliberately leaves the Table III term out, so that all 13
      differences are 0; the omitted term is printed separately on the last line of (A).
  (B) the nonlocal parts of (I.152), (I.153) reproduce (4/y)^2, (4/y)^3 pointwise at D = 4, and the
      two terms of (I.159) cancel pointwise (L[8/y - 4 ln(y/4)] = 12).
  (C) beta_1 = 0; beta_2, beta_3 assembled from (I.44), (I.59) and Table VIII agree with (I.162),
      (I.164) for general D; D -> 4 expansions (I.163) = -61/3 and (I.165) = -4/(D-4) + 58/3 + 2 gamma;
      the values after the corrections of (I.40) and (I.59) (-49/3 and -4/(D-4) + 52/3 + 2 gamma) and the
      definition-sign values (x -1).

Conventions: paper's sign convention unless stated; x = y/4; operators named by their D = 4 form;
beta_i in units H^(D-4)/(4 pi)^(D/2) (the symbol e below is D - 4).
Expected result: 13 zero differences, the Table III zeta term -16 x log(x) + 8 log(x), True/0 in (B),
True for (162), (164), and the numbers -61/3, 2 EulerGamma + 58/3 - 4/e, -49/3, 2 EulerGamma + 52/3 - 4/e.
Runtime: about 2 s (Python 3.10+, SymPy 1.14).
"""
import sympy as sp

R = sp.Rational
x = sp.Symbol('x', positive=True); lx = sp.log(x)

print("== (A) Table X assembled from Tables III, V, VII and Eq. (159) (paper's convention, gamma3/eps3 eliminated) ==")
# operators named by their D=4 form
ops = ['a4B3', 'a4B2', 'a4H2B', 'a4H4', 'a3sB3', 'a3sB2', 'a3sH2B', 'a3sH4', 'a3H2L', 'a3LB', 'a2sH2L', 'a2sLB', 'a2L2']
names = {'a4B3': "(aa')^4 Box^3/H^2", 'a4B2': "(aa')^4 Box^2", 'a4H2B': "(aa')^4 H^2 Box", 'a4H4': "(aa')^4 H^4",
         'a3sB3': "(aa')^3(a^2+a'^2) Box^3/H^2", 'a3sB2': "(aa')^3(a^2+a'^2) Box^2", 'a3sH2B': "(aa')^3(a^2+a'^2) H^2 Box",
         'a3sH4': "(aa')^3(a^2+a'^2) H^4", 'a3H2L': "(aa')^3 H^2 Lap", 'a3LB': "(aa')^3 Lap Box",
         'a2sH2L': "(aa')^2(a^2+a'^2) H^2 Lap", 'a2sLB': "(aa')^2(a^2+a'^2) Lap Box", 'a2L2': "(aa')^2 Lap^2"}
X = {o: sp.Integer(0) for o in ops}
def add(o, f, src): X[o] += f
# Table VII (operator x [Box/H^2]^n): alpha=(aa')^4 Box^2, beta=(aa')^3(a^2+a'^2)H^2 Box, gamma1=(aa')^4H^4, gamma2=(aa')^3(a^2+a'^2)H^4,
#                                       eps1=(aa')^3H^2 Lap, eps2=(aa')^2(a^2+a'^2)H^2 Lap
add('a4B3', -lx/(3*x), 'VII alpha'); add('a4B2', 2*lx/(3*x) - 1/(3*x), 'VII alpha')
add('a3sB3', lx/(6*x), 'VII beta'); add('a3sB2', -lx/(3*x) + 1/(6*x), 'VII beta')
add('a4B2', 2*lx/x, 'VII gamma1'); add('a4H2B', -6*lx/x + 2/x, 'VII gamma1'); add('a4H4', 4*lx/x - 2/x, 'VII gamma1')
add('a3sH2B', -2*lx/(3*x), 'VII gamma2'); add('a3sH4', 4*lx/(3*x) - 2/(3*x), 'VII gamma2')
add('a3LB', lx/(3*x), 'VII eps1'); add('a3H2L', -2*lx/(3*x) + 1/(3*x), 'VII eps1')
add('a2sLB', 7*lx/(6*x), 'VII eps2'); add('a2sH2L', -7*lx/(3*x) + 7/(6*x), 'VII eps2')
# Table V
add('a4B2', 6*lx/x + 13/x - 6*lx**2 - 18*lx, 'V alpha'); add('a3sH2B', 5/x - 18*lx, 'V beta')
add('a4H4', 30/x - 108*lx - 36, 'V gamma1'); add('a3sH4', -5/x - 18, 'V gamma2')
g3 = -5/x - 36; add('a4H4', 2*g3, 'V gamma3->gamma1'); add('a3sH4', g3, 'V gamma3->gamma2')
add('a2sLB', -4*lx/x - R(49, 6)/x + 4*lx**2 + 10*lx + 12*x*lx, 'V delta')
add('a3H2L', -R(26, 3)/x + 60*lx - 120*x*lx + 72*x, 'V eps1'); add('a2sH2L', R(13, 6)/x - 12*lx + 12*x*lx, 'V eps2')
e3 = -R(11, 6)/x + 36*x*lx + 12*x; add('a3H2L', 2*e3, 'V eps3->eps1'); add('a2sH2L', e3, 'V eps3->eps2')
add('a2L2', R(10, 3)*lx - 24*x*lx + 24*x**2*lx - 36*x**2, 'V zeta')
# Table III: (aa')^3H^2Lap, (aa')^2(a+a')^2 H^2 Lap = eps3 -> 2 eps1 + eps2, (aa')^2 Lap^2
add('a3H2L', -12/x + 24*lx, 'III'); add('a3H2L', 2*(4/x), 'III eps3->eps1'); add('a2sH2L', 4/x, 'III eps3->eps2')
TIII_zeta = 8*lx - 16*x*lx
# (159)
add('a4H2B', 2/x - 4*lx, '(159)'); add('a4H4', -12, '(159)')
paperX = {'a4B3': -lx/(3*x), 'a4B2': 26*lx/(3*x) + 38/(3*x) - 6*lx**2 - 18*lx, 'a4H2B': -6*lx/x + 4/x - 4*lx,
          'a4H4': 4*lx/x + 18/x - 120 - 108*lx, 'a3sB3': lx/(6*x), 'a3sB2': -lx/(3*x) + 1/(6*x),
          'a3sH2B': -2*lx/(3*x) + 5/x - 18*lx, 'a3sH4': 4*lx/(3*x) - 32/(3*x) - 54,
          'a3H2L': -2*lx/(3*x) - 16/x + 84*lx - 48*x*lx + 96*x, 'a3LB': lx/(3*x),
          'a2sH2L': -7*lx/(3*x) + 11/(2*x) - 12*lx + 48*x*lx + 12*x,
          'a2sLB': -17*lx/(6*x) - 49/(6*x) + 4*lx**2 + 10*lx + 12*x*lx, 'a2L2': R(10, 3)*lx - 24*x*lx + 24*x**2*lx - 36*x**2}
for o in ops:
    d = sp.simplify(sp.expand(X[o] - paperX[o]))
    print(f"  {names[o]:30s}: assembled - paper = {d}")
print(f"  Table III zeta contribution not contained in paper's (aa')^2 Lap^2 row: {TIII_zeta}")

print("== (B) (152), (153) nonlocal parts reproduce (4/y)^2, (4/y)^3 pointwise at D=4 ==")
Y = sp.Symbol('y', positive=True)
L4 = lambda f: sp.expand((4*Y - Y**2)*sp.diff(f, Y, 2) + 4*(2 - Y)*sp.diff(f, Y))
f = (4/Y)*sp.log(Y/4)
N152 = -L4(f) + 2*f - 4/Y
N153 = -R(1, 2)*L4(L4(f)) + L4(2*f - R(1, 2)*4/Y) - 2*f + 4/Y
print("  N152 - (4/y)^2 =", sp.simplify(N152 - (4/Y)**2), ";  N153 - (4/y)^3 =", sp.simplify(N153 - (4/Y)**3))
print("  (159) pointwise: H^2 Box[8/y - 4 ln(y/4)]/H^4 - 12 =", sp.simplify(L4(8/Y - 4*sp.log(Y/4)) - 12))

print("== (C) beta coefficients (160)-(165) ==")
D = sp.Symbol('D'); e = sp.Symbol('e')
k_ = sp.gamma(D - 1)/sp.gamma(D/2)                      # k / (H^{D-2}/(4pi)^{D/2})
pih = -sp.pi*sp.cot(sp.pi*D/2)
U8 = sp.gamma(D/2)/((D - 3)*(D - 4))                    # Table VIII unit / (H^{D-4}/(4pi)^{D/2})
VIII = [0, (D - 4)*(-D**3 + D - 4)/(2*(D - 1)*(D - 2)), (D**3 - 16*D**2 + 28*D - 16)/(2*(D - 1)*(D - 2))]
# paper's (44) and (59), coefficient of i kappa^2 k a^D Box and a^{D-2} Lap (log and 9H^2 parts cancel against (67))
p44 = (D*(D - 1)/4*pih - 2, -D*pih + 2); p59 = (-(D - 1)/(D - 3)*pih + 2, (D - 1)/(D - 3)*pih - R(7, 2))
b2 = sp.simplify(k_*(p44[0] + p59[0]) + U8*VIII[1])
b3 = sp.simplify(k_*(p44[1] + p59[1]) + U8*VIII[2])
b2p = (-D**3 + D - 4)*sp.gamma(D/2 - 1)/(4*(D - 1)*(D - 3)) - (D + 1)*(D - 4)*sp.gamma(D)*sp.pi*sp.cot(sp.pi*D/2)/(4*(D - 3)*sp.gamma(D/2))
b3p_exact = (D**3 - 16*D**2 + 28*D - 16)*sp.gamma(D/2)/(2*(D - 1)*(D - 2)*(D - 3)*(D - 4)) \
            + (D**2 - 4*D + 1)*sp.gamma(D - 1)*sp.pi*sp.cot(sp.pi*D/2)/((D - 3)*sp.gamma(D/2)) - R(3, 2)*k_
print("  beta1 = Table VIII col1 total =", VIII[0])
print("  beta2 assembled == (162):", sp.simplify(sp.gammasimp(b2 - b2p)) == 0)
print("  beta3 assembled == (164) [with -3 written exactly as -(3/2)Gamma(D-1)/Gamma(D/2)]:", sp.simplify(sp.gammasimp(b3 - b3p_exact)) == 0)
ser = lambda expr: sp.simplify(sp.series(expr.subs(D, 4 + e), e, 0, 1).removeO())
print("  (163) beta2 at D=4:", sp.nsimplify(ser(b2p)))
s3 = ser(b3p_exact)
print("  (165) beta3 expansion:", sp.expand(sp.simplify(s3)), "   [paper: -4/(D-4) + 58/3 + 2 EulerGamma]")
# corrections: (40) -> no -2/+2 in (44); (59) Lap finite -5/2 instead of -7/2
c44 = (D*(D - 1)/4*pih, -D*pih); c59 = (p59[0], (D - 1)/(D - 3)*pih - R(5, 2))
b2c = k_*(c44[0] + c59[0]) + U8*VIII[1]; b3c = k_*(c44[1] + c59[1]) + U8*VIII[2]
print("  corrected (paper convention): beta2(D=4) =", sp.nsimplify(ser(b2c)), ";  beta3 =", sp.expand(sp.simplify(ser(b3c))))
print("  definition convention (x -1): beta2 =", -sp.nsimplify(ser(b2c)), ";  beta3 =", sp.expand(-sp.simplify(ser(b3c))))
print("  divergent parts: beta2 has no pole:", sp.limit(b2p, D, 4) == sp.Rational(-61, 3), "; beta3 pole residue:", sp.limit((D - 4)*b3p_exact, D, 4))
