"""check_tables45.py -- actual D = 4 propagators: appendix totals vs direct sum, and vs Tables IV+V  [Task 13, checks (F4c), (T5)]

Purpose
  (F4c) Repeat the global check of check_global.py with the actual D = 4 propagator functions
        (units H^2/16 pi^2):  A = 4/y - 2 ln(y/4) - 1  (I.62),  B = C = 4/y,  DeltaC = 4 (C - A).
  (T5)  For each of the ten external operators, print
            [derived D = 4 total]  -  [printed Table IV at D = 4  +  printed Table V],
        written in x = y/4, and evaluate the printed Tables IV+V under the ten operators at three
        random space-like points against the direct sum of P1..P7.

Expected result
  (F4c): diff = 0.0 at three points.
  (T5):  per-operator differences: alpha, beta: constants; gamma2: 4/(3 x^2) (the Table IV gamma2
         sign); delta: -12 x + (ln 2 terms) ...; the printed tables miss the direct sum by
         0.0356, 17.05, 0.607 at the three points.
Runtime: about 20 s.
Run:     python3 check_tables45.py   (imports check_global.py from this folder)
"""
import sympy as sp, random
from check_global import run, Y, H, R
A = 4/Y - 2*sp.log(Y/4) - 1; B = 4/Y; C = 4/Y; dC = 4*(C - A)
print("== (F4c) actual D=4 propagators (units H^2/16pi^2), linked antiderivatives ==")
T, ops, direct, opsum, geo = run(4, A, B, dC, 'D=4 actual', pts=3)
x = sp.Symbol('x', positive=True)
lnx = sp.log(x)
c1 = {'al': 0, 'be': -R(1, 3), 'g1': -2, 'g2': 1, 'g3': -1, 'de': 0, 'e1': 0, 'e2': 0, 'e3': 0, 'ze': 0}
c2 = {'al': R(1, 3), 'be': -R(1, 3), 'g1': -2, 'g2': -R(2, 3), 'g3': 0, 'de': 0, 'e1': -R(1, 3), 'e2': -R(7, 6), 'e3': 0, 'ze': 0}
V = {'al': 6*lnx/x + 13/x - 6*lnx**2 - 18*lnx, 'be': 5/x - 18*lnx, 'g1': 30/x - 108*lnx - 36, 'g2': -5/x - 18, 'g3': -5/x - 36,
     'de': -4*lnx/x - R(49, 6)/x + 4*lnx**2 + 10*lnx + 12*x*lnx, 'e1': -R(26, 3)/x + 60*lnx - 120*x*lnx + 72*x,
     'e2': R(13, 6)/x - 12*lnx + 12*x*lnx, 'e3': -R(11, 6)/x + 36*x*lnx + 12*x, 'ze': R(10, 3)*lnx - 24*x*lnx + 24*x**2*lnx - 36*x**2}
paper = {k: (c1[k]/x**3 + c2[k]/x**2 + V[k]).subs(x, Y/4) for k in V}
print("== per-operator: my D=4 total minus paper's Table IV(D=4)+V, written in x = y/4 ==")
for k in T:
    d = sp.expand(sp.expand_log(sp.expand(T[k] - paper[k]).subs(Y, 4*x), force=True))
    print(f"  {k}: {sp.nsimplify(sp.collect(d, [lnx]))}")
t, tp, xs, ys, a, ap, y = geo
paper_sum = sum(ops[k](paper[k]) for k in paper)
random.seed(23)
for kk in range(3):
    while True:
        pt = {H: random.choice([R(7, 10), R(13, 10), R(3, 2)]), t: -R(random.randint(3, 9), 4), tp: -R(random.randint(10, 19), 4)}
        for v in xs + ys: pt[v] = R(random.randint(-9, 9), 3)
        if sum((pt[xs[i]] - pt[ys[i]])**2 for i in range(len(xs))) > (pt[t] - pt[tp])**2 + 1: break
    dv, pv = sp.N(direct.subs(pt), 30), sp.N(paper_sum.subs(pt), 30)
    print(f"  point {kk}: sum P_k = {sp.Float(dv, 12)}  Tables IV+V = {sp.Float(pv, 12)}  diff = {sp.Float(dv - pv, 4)}")
