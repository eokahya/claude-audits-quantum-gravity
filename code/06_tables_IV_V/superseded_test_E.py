"""superseded_test_E.py -- the first (superseded) global test (E) of Task 13

Purpose
  Reproduces part (E) of the FIRST version of check_94_125.py (Task 13 command log). It compared
  the printed Tables IV(D=4)+V, acted on by the ten external operators, with the direct normal
  contribution written with the DERIVED signs of Task 7 (i.e. -1 times sum_k P_k of the appendix),
  using the actual D = 4 propagators A = 4/y - 2 ln(y/4) - 1, B = C = 4/y (kappa^2 h^2 stripped).
  Neither "direct + tables" nor "direct - tables" vanished. The test was removed from
  check_94_125.py and replaced by the structured checks check_global.py / check_tables45*.py,
  which locate the two defects of the printed tables (Table IV gamma2 sign; a finite term
  equivalent to 72 gamma2 in Table V).

  The random generator is advanced exactly as in the original script (seed 17, two points drawn
  by part (A) for D = 4 and D = 5 before part (E)), so the three test points are the original ones.

Expected result (from the original log)
  point 0: direct = -1.26640166984,  Tables IV+V = 0.893758671558, direct + tables = -0.373, direct - tables = -2.16
  point 1: direct = 0.0189424187218, Tables IV+V = -0.893296577647, direct + tables = -0.874, direct - tables = 0.912
  point 2: direct = 0.000241338543333, Tables IV+V = -0.00541935945134, direct + tables = -0.00518, direct - tables = 0.00566
Runtime: about 15 s.
Run:     python3 superseded_test_E.py
"""
import sympy as sp, random
R = sp.Rational
H = sp.Symbol('H', positive=True)
al, be, ga = sp.symbols('alpha beta gamma')
Yv = sp.Symbol('Y', positive=True)
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
# advance the random state exactly as part (A) of the original script did (one point for D=4, one for D=5)
for Dv in (4, 5):
    t_, tp_, xs_, ys_, *_ = setup(Dv)
    rnd(t_, tp_, xs_, ys_)

print("== (E) global pointwise check of Tables IV(D=4)+V against the direct normal contribution, D=4 ==")
t, tp, xs, ys, a, ap, y = setup(4)
x = y/4
A0 = 4/y - 2*sp.log(y/4) - 1                    # A(y)/h at D=4 (outer and inner)
Bf = 4/y; Cf = 4/y
X = (t,) + xs; Xp = (tp,) + ys
inner = lambda m, n: sp.diff(A0, X[m], Xp[n])
T = lambda f, mn, rs: sp.diff((a*ap)**2*f*inner(*rs), X[mn[0]], Xp[mn[1]])
Sp = range(1, 4)
direct = (-sum(T(A0, (i, i), (j, j)) for i in Sp for j in Sp) + 3*T(A0, (0, 0), (0, 0))
          - sum(T(A0, (i, 0), (i, 0)) for i in Sp) - sum(T(A0, (0, i), (0, i)) for i in Sp)
          + sum(T(Bf, (0, 0), (j, j)) for j in Sp) + sum(T(Bf, (i, 0), (0, i)) for i in Sp)
          + sum(T(Bf, (0, i), (i, 0)) for i in Sp) + sum(T(Bf, (i, i), (0, 0)) for i in Sp)
          - 4*T(Cf, (0, 0), (0, 0)))                                    # derived signs (Task 7), kappa^2 h^2 stripped
lnx = sp.log(Yv/4); xv = Yv/4
c1 = {'al': 0, 'be': -R(1, 3), 'g1': -2, 'g2': 1, 'g3': -1, 'de': 0, 'e1': 0, 'e2': 0, 'e3': 0, 'ze': 0}
c2 = {'al': R(1, 3), 'be': -R(1, 3), 'g1': -2, 'g2': -R(2, 3), 'g3': 0, 'de': 0, 'e1': -R(1, 3), 'e2': -R(7, 6), 'e3': 0, 'ze': 0}
V = {'al': 6*lnx/xv + 13/xv - 6*lnx**2 - 18*lnx, 'be': 5/xv - 18*lnx, 'g1': 30/xv - 108*lnx - 36, 'g2': -5/xv - 18, 'g3': -5/xv - 36,
     'de': -4*lnx/xv - R(49, 6)/xv + 4*lnx**2 + 10*lnx + 12*xv*lnx, 'e1': -R(26, 3)/xv + 60*lnx - 120*xv*lnx + 72*xv,
     'e2': R(13, 6)/xv - 12*lnx + 12*xv*lnx, 'e3': -R(11, 6)/xv + 36*xv*lnx + 12*xv, 'ze': R(10, 3)*lnx - 24*xv*lnx + 24*xv**2*lnx - 36*xv**2}
Fk = {k: c1[k]/xv**3 + c2[k]/xv**2 + V[k] for k in V}                    # in units of kappa^2 H^4/(4pi)^4 = kappa^2 h^2
boxY = lambda f: H**2*((4*Yv - Yv**2)*sp.diff(f, Yv, 2) + 4*(2 - Yv)*sp.diff(f, Yv))   # Eq. (93), D=4 (pointwise)
lap = lambda e: sum(sp.diff(e, xi, 2) for xi in xs)
s2 = a**2 + ap**2
ops = {'al': (a*ap)**4*F(boxY(boxY(Fk['al'])), y), 'be': (a*ap)**3*s2*H**2*F(boxY(Fk['be']), y),
       'g1': (a*ap)**4*H**4*F(Fk['g1'], y), 'g2': (a*ap)**3*s2*H**4*F(Fk['g2'], y), 'g3': (a*ap)**3*(a + ap)**2*H**4*F(Fk['g3'], y),
       'de': (a*ap)**2*s2*lap(F(boxY(Fk['de']), y)), 'e1': (a*ap)**3*H**2*lap(F(Fk['e1'], y)),
       'e2': (a*ap)**2*s2*H**2*lap(F(Fk['e2'], y)), 'e3': (a*ap)**2*(a + ap)**2*H**2*lap(F(Fk['e3'], y)), 'ze': (a*ap)**2*lap(lap(F(Fk['ze'], y)))}
tables = sum(ops.values())
for k in range(3):
    pt = rnd(t, tp, xs, ys, spacelike=True)
    dv, tv = num(direct, pt), num(tables, pt)
    print(f"  point {k}: direct = {sp.Float(dv, 12)},  Tables IV+V = {sp.Float(tv, 12)},  direct + tables = {sp.Float(dv + tv, 3)},  direct - tables = {sp.Float(dv - tv, 3)}")
