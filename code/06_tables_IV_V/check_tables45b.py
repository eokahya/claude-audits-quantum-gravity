"""check_tables45b.py -- hypothesis tests for the Table IV+V inconsistency  [Task 13, check (T5)]

Purpose
  Evaluate several versions of Tables IV(D=4)+V under the ten external operators and compare
  with the direct sum of P1..P7 (actual D = 4 propagators) at three random space-like points:
    1. as printed;
    2. Table IV gamma2 second column with the sign corrected (-2/3 -> +2/3 at D = 4);
    3. gamma2 sign corrected + a term -12x added to the delta row of Table V;
    4. gamma2 sign corrected + gamma2 constant of Table V changed from -18 to +54 (i.e. +72 gamma2);
    5. only -12x in the delta row (gamma2 sign as printed).
  Also prints the derived D = 4 operator functions with ln y = ln x + ln 4 made explicit.

Expected result
  1: [0.03560, 17.05, 0.6071];  2: [0.03509, 16.97, 0.6024];  3 and 4: [0.0, 0.0, 0.0];
  5: [0.0005088, 0.07447, 0.004686].
Runtime: about 30 s.
Run:     python3 check_tables45b.py   (imports check_global.py from this folder)
"""
import sympy as sp, random
from check_global import run, Y, H, R
A = 4/Y - 2*sp.log(Y/4) - 1; B = 4/Y; C = 4/Y; dC = 4*(C - A)
T, ops, direct, opsum, geo = run(4, A, B, dC, 'D=4 actual', pts=0)
x = sp.Symbol('x', positive=True); lnx = sp.log(x)
c1 = {'al': 0, 'be': -R(1, 3), 'g1': -2, 'g2': 1, 'g3': -1, 'de': 0, 'e1': 0, 'e2': 0, 'e3': 0, 'ze': 0}
c2 = {'al': R(1, 3), 'be': -R(1, 3), 'g1': -2, 'g2': -R(2, 3), 'g3': 0, 'de': 0, 'e1': -R(1, 3), 'e2': -R(7, 6), 'e3': 0, 'ze': 0}
V = {'al': 6*lnx/x + 13/x - 6*lnx**2 - 18*lnx, 'be': 5/x - 18*lnx, 'g1': 30/x - 108*lnx - 36, 'g2': -5/x - 18, 'g3': -5/x - 36,
     'de': -4*lnx/x - R(49, 6)/x + 4*lnx**2 + 10*lnx + 12*x*lnx, 'e1': -R(26, 3)/x + 60*lnx - 120*x*lnx + 72*x,
     'e2': R(13, 6)/x - 12*lnx + 12*x*lnx, 'e3': -R(11, 6)/x + 36*x*lnx + 12*x, 'ze': R(10, 3)*lnx - 24*x*lnx + 24*x**2*lnx - 36*x**2}
t, tp, xs, ys, a, ap, y = geo
def glob(c2mod, Vmod, label):
    paper = {k: (c1[k]/x**3 + c2mod[k]/x**2 + Vmod[k]).subs(x, Y/4) for k in V}
    s = sum(ops[k](paper[k]) for k in paper)
    random.seed(23); res = []
    for kk in range(3):
        while True:
            pt = {H: random.choice([R(7, 10), R(13, 10), R(3, 2)]), t: -R(random.randint(3, 9), 4), tp: -R(random.randint(10, 19), 4)}
            for v in xs + ys: pt[v] = R(random.randint(-9, 9), 3)
            if sum((pt[xs[i]] - pt[ys[i]])**2 for i in range(len(xs))) > (pt[t] - pt[tp])**2 + 1: break
        dv, pv = sp.N(direct.subs(pt), 30), sp.N(s.subs(pt), 30)
        res.append(dv - pv)
    print(f"  {label}: sum P_k - tables at 3 points =", [sp.Float(r, 4) for r in res])
glob(c2, V, 'as printed')
c2f = dict(c2); c2f['g2'] = +R(2, 3)
glob(c2f, V, 'Table IV gamma2 col2 sign corrected')
# search: which single additional finite-part change would close the gap? print my totals in a normalized convention
print("== my D=4 totals with ln Y -> ln x + ln 4 made explicit (x = y/4) ==")
for k in T:
    e = sp.expand(sp.expand_log(T[k].subs(Y, 4*x), force=True))
    print(f"  {k}: {sp.collect(sp.nsimplify(e), [lnx])}")
print("== hypotheses for the remaining gap ==")
Vd = dict(V); Vd['de'] = V['de'] - 12*x
glob(c2f, Vd, 'gamma2 sign fixed + delta row gets -12x')
Vg = dict(V); Vg['g2'] = V['g2'] + 72
glob(c2f, Vg, 'gamma2 sign fixed + gamma2 constant -18 -> +54')
glob(c2, Vd, 'only delta -12x (gamma2 as printed)')
