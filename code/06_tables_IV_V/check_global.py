"""check_global.py -- global pointwise check of the appendix decomposition (Tables XI-XX)  [Task 13, check (F4a,b)]

Purpose
  The normal part of the 3-point contribution is split in the appendix of Paper I into seven
  parts P1..P7, Eqs. (A1)-(A7) (naive inner derivatives, paper's sign = sign of the printed
  (I.47)-(I.49)). The appendix reduces them to ten external operators alpha, beta, gamma1..3,
  delta, eps1..3, zeta [(I.104)-(I.113)] acting on the totals of Tables XI-XX.
  This script codes the ten totals (with the five typographical corrections described in the
  notes) and checks pointwise that

        sum_k P_k  ==  sum_ops  op[ total_op ]

  in explicit conformal coordinates.

Method
  run(D, A, B, dC, ...) builds both sides for given functions A(y), B(y), DeltaC(y).
  First integrals are linked by partial integration, I[f A''] = f A' - I[f' A'], so that the
  integration constants are consistent; higher integrals are sympy antiderivatives.
  Box acting on functions of y is evaluated with (I.93) (pointwise).
  When run as a script: generic powers A = y^(-7/3), B = y^(-5/4), DeltaC = y^(2/7)
  at D = 5 and D = 4, two random space-like points each, 30-digit arithmetic.
  The function run() is imported by check_tables45.py and check_tables45b.py.

Expected result: diff = 0.0 at every point (zero to the 30-digit working precision).
Runtime: about 15 s.
Run:     python3 check_global.py
"""
# Global check: sum of the paper's seven parts P1..P7 (appendix (A1)-(A7)) versus the ten external
# operators acting on the appendix totals (Tables XI-XX), and versus Tables IV+V.
import sympy as sp, random
R = sp.Rational
H = sp.Symbol('H', positive=True)
Y = sp.Symbol('Y', positive=True)

def setup(Dv):
    t, tp = sp.symbols('t tp', negative=True)
    xs = sp.symbols(f'x1:{Dv}', real=True); ys = sp.symbols(f'y1:{Dv}', real=True)
    a, ap = -1/(H*t), -1/(H*tp)
    y = a*ap*H**2*(sum((xs[i] - ys[i])**2 for i in range(Dv - 1)) - (t - tp)**2)
    return t, tp, xs, ys, a, ap, y

def run(Dv, A, B, dC, label, pts=3, Bfun_for_P=None):
    D = sp.Integer(Dv)
    t, tp, xs, ys, a, ap, y = setup(Dv); p = Dv - 2
    Ap, Bp, dCp = sp.diff(A, Y), sp.diff(B, Y), sp.diff(dC, Y)
    App = sp.diff(A, Y, 2)
    intg = lambda g: sp.integrate(sp.expand(g), Y, conds='none')
    base = {}
    base['A2'] = [intg(Ap**2)]; base['AApp'] = [A*Ap - base['A2'][0]]
    base['ApBp'] = [intg(Ap*Bp)]; base['BApp'] = [B*Ap - base['ApBp'][0]]
    base['ApdCp'] = [intg(Ap*dCp)]; base['dCApp'] = [dC*Ap - base['ApdCp'][0]]
    def In(name, n):
        if name == 'ApdBp': return In('ApBp', n) - In('A2', n)
        if name == 'dBApp': return In('BApp', n) - In('AApp', n)
        L = base[name]
        while len(L) < n: L.append(intg(L[-1]))
        return L[n - 1]
    prod = {'A2': Ap**2, 'ApBp': Ap*Bp, 'ApdBp': Ap*(Bp - Ap), 'ApdCp': Ap*dCp, 'BApp': B*App}
    dY = lambda e: sp.diff(e, Y)
    T = {}
    T['al'] = In('AApp', 2) + In('dCApp', 2)
    T['be'] = D*In('ApBp', 1) - In('ApdBp', 1) - R(1, 2)*Y*In('A2', 1) - (D - 1)/2*In('A2', 2) + In('ApdCp', 1) - R(1, 2)*Y*In('ApdCp', 1) - (D - 1)/2*In('ApdCp', 2)
    T['g1'] = ((D)*(D - 1)*(D - 2)*In('ApBp', 1) - (D - 1)*(D - 2)*In('ApdBp', 1) - R(1, 2)*(D - 1)*(D - 2)*Y*In('A2', 1) - R(1, 2)*(D - 1)**2*(D - 2)*In('A2', 2)
               + (D - 1)*(D - 2)*In('ApdCp', 1) - R(1, 2)*(D - 1)*(D - 2)*Y*In('ApdCp', 1) - R(1, 2)*(D - 1)**2*(D - 2)*In('ApdCp', 2))
    T['g2'] = (D*Y*prod['ApBp'] - Y*prod['ApdBp'] - R(1, 2)*Y**2*prod['A2'] - D/2*Y*In('A2', 1) + Y*prod['ApdCp'] - D/2*Y*In('ApdCp', 1) - R(1, 2)*Y**2*prod['ApdCp'])
    T['g3'] = ((D - 1)*D*Y*prod['ApBp'] + D*Y**2*dY(prod['ApBp']) - (D - 1)*Y*prod['ApdBp'] - Y**2*dY(prod['ApdBp']) - D*Y**2*prod['A2'] - R(1, 2)*Y**3*dY(prod['A2'])
               - R(1, 2)*(D - 1)*D*Y*In('A2', 1) + (D - 1)*Y*prod['ApdCp'] + Y**2*dY(prod['ApdCp']) - D*Y**2*prod['ApdCp'] - R(1, 2)*Y**3*dY(prod['ApdCp'])
               - R(1, 2)*(D - 1)*D*Y*In('ApdCp', 1))
    T['de'] = 2*In('dBApp', 2) + R(1, 4)*In('A2', 3) - In('dCApp', 2) + R(1, 4)*In('ApdCp', 3)
    T['e1'] = (-2*(D - 2)*In('ApdBp', 1) - Y*In('ApdBp', 1) - (5*D - 4)/2*In('ApdBp', 2) + D**2/2*In('ApBp', 2) - D/4*Y*In('A2', 2) + (D - 2)**2*In('BApp', 2)
               - R(3, 4)*(D - 2)*In('A2', 3) - 2*In('ApdCp', 1) + Y*In('ApdCp', 1) + (3*D - 2)/2*In('ApdCp', 2) - D/4*Y*In('ApdCp', 2) - R(3, 4)*(D - 2)*In('ApdCp', 3))
    T['e2'] = 2*Y*In('BApp', 1) + In('ApBp', 2) - (D - 1)*In('ApdBp', 2) + R(1, 4)*Y*In('A2', 2) + R(1, 4)*Y*In('ApdCp', 2)
    T['e3'] = (2*Y**2*prod['BApp'] + 2*(D - 1)*Y*In('BApp', 1) + (D - 2)*In('ApBp', 2) + Y*In('A2', 1) + (D - 1)/4*Y*In('A2', 2) + R(1, 4)*Y**2*In('A2', 1)
               + (D - 1)/4*Y*In('ApdCp', 2) + R(1, 4)*Y**2*In('ApdCp', 1))
    T['ze'] = -4*In('dBApp', 2) + In('ApdBp', 3) + (D - 2)/8*In('A2', 4) + In('dCApp', 2) - R(1, 2)*In('ApdCp', 3) + (D - 2)/8*In('ApdCp', 4)
    F = lambda e: e.subs(Y, y)
    lap = lambda e: sum(sp.diff(e, xi, 2) for xi in xs)
    boxY = lambda f: H**2*((4*Y - Y**2)*sp.diff(f, Y, 2) + D*(2 - Y)*sp.diff(f, Y))
    s2 = a**2 + ap**2
    ops = {'al': lambda f: (a*ap)**Dv*F(boxY(boxY(f))), 'be': lambda f: (a*ap)**(Dv - 1)*s2*H**2*F(boxY(f)),
           'g1': lambda f: (a*ap)**Dv*H**4*F(f), 'g2': lambda f: (a*ap)**(Dv - 1)*s2*H**4*F(f), 'g3': lambda f: (a*ap)**(Dv - 1)*(a + ap)**2*H**4*F(f),
           'de': lambda f: (a*ap)**p*s2*lap(F(boxY(f))), 'e1': lambda f: (a*ap)**(Dv - 1)*H**2*lap(F(f)),
           'e2': lambda f: (a*ap)**p*s2*H**2*lap(F(f)), 'e3': lambda f: (a*ap)**p*(a + ap)**2*H**2*lap(F(f)), 'ze': lambda f: (a*ap)**p*lap(lap(F(f)))}
    # direct: P1..P7 with naive inner derivatives (paper's signs)
    X = (t,) + xs; Xp = (tp,) + ys
    Ay, By, dCy = F(A), F(B), F(dC)
    inner = lambda m, n: sp.diff(Ay, X[m], Xp[n])
    Tm = lambda f, mn, rs: sp.diff((a*ap)**p*f*inner(*rs), X[mn[0]], Xp[mn[1]])
    Sp = range(1, Dv)
    P1 = sum(Tm(Ay, (i, i), (j, j)) for i in Sp for j in Sp)
    P2 = sum(Tm(Ay, (i, 0), (i, 0)) + Tm(Ay, (0, i), (0, i)) for i in Sp)
    P3 = Tm(Ay, (0, 0), (0, 0))
    P4 = -sum(Tm(By, (0, 0), (j, j)) for j in Sp)
    P5 = -sum(Tm(By, (i, 0), (0, i)) + Tm(By, (0, i), (i, 0)) for i in Sp)
    P6 = -sum(Tm(By, (i, i), (0, 0)) for i in Sp)
    P7 = Tm(dCy, (0, 0), (0, 0))
    direct = P1 + P2 + P3 + P4 + P5 + P6 + P7
    opsum = sum(ops[k](T[k]) for k in T)
    random.seed(23)
    out = []
    for kk in range(pts):
        while True:
            pt = {H: random.choice([R(7, 10), R(13, 10), R(3, 2)]), t: -R(random.randint(3, 9), 4), tp: -R(random.randint(10, 19), 4)}
            for v in xs + ys: pt[v] = R(random.randint(-9, 9), 3)
            if sum((pt[xs[i]] - pt[ys[i]])**2 for i in range(len(xs))) > (pt[t] - pt[tp])**2 + 1: break
        dv, ov = sp.N(direct.subs(pt), 30), sp.N(opsum.subs(pt), 30)
        out.append((dv, ov))
        print(f"  [{label}] point {kk}: sum P_k = {sp.Float(dv, 14)}   sum ops x appendix totals = {sp.Float(ov, 14)}   diff = {sp.Float(dv - ov, 3)}")
    return T, ops, direct, opsum, (t, tp, xs, ys, a, ap, y)

if __name__ == '__main__':
    print("== (F4a) generic powers A=y^b, B=y^g, dC=y^e, D=5 ==")
    run(5, Y**R(-7, 3), Y**R(-5, 4), Y**R(2, 7), 'powers D=5', pts=2)
    print("== (F4b) generic powers, D=4 ==")
    run(4, Y**R(-7, 3), Y**R(-5, 4), Y**R(2, 7), 'powers D=4', pts=2)
