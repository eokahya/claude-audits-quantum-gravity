"""check_tableX_global.py -- global pointwise test of Table X of Paper I (Task 17).

Paper I: E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

Compares, at D = 4 and at three random spacelike-separated point pairs (x != x'), with 30-digit arithmetic:
  LHS = direct sum of the paper's normal contributions P1..P7 (naive inner derivatives, actual D = 4
        propagator functions A = h{4/y - 2 ln(y/4) - 1}, B = C = 4h/y, DeltaC = 4(C - A), in units
        h = H^2/16 pi^2), plus the paper's nonlocal term of (I.67),
        -(kappa^2 H^2/8 pi^2) 2H (aa')^2 (a d_0 + a' d_0') Lap A;
  RHS = the 13 operators of Table X acting pointwise on their coefficient functions.
Both sides are in the paper's sign convention and in units kappa^2 h^2 = kappa^2 H^4/(4 pi)^4; neither
depends on the external-operator reduction of the appendix.  Eq. (I.159) is not included because its
two terms cancel pointwise.

Printed for each point: LHS - Table X (as printed); the same after adding 8 ln x - 16 x ln x to the
(aa')^2 Lap^2 row ("zeta fix"); and after also adding the constant 72 to the (aa')^3(a^2+a'^2)H^4 row
("72 gamma2", the Table V inconsistency found in Task 13).
Expected result: nonzero differences for the printed table and with the zeta fix alone; differences of
order 1e-31 (zero to working precision) when both corrections are applied.
Requires check_global.py in the same folder.  Runtime: about 16 s (Python 3.10+, SymPy 1.14).
"""
import sympy as sp, random
from check_global import run, Y, H, R
A = 4/Y - 2*sp.log(Y/4) - 1; B = 4/Y; C = 4/Y; dC = 4*(C - A)          # units h = H^2/(16 pi^2)
T, ops_unused, direct, opsum, geo = run(4, A, B, dC, 'D=4', pts=0)
t, tp, xs, ys, a, ap, y = geo
F = lambda e: e.subs(Y, y)
lap = lambda e: sum(sp.diff(e, xi, 2) for xi in xs)
boxY = lambda f: H**2*((4*Y - Y**2)*sp.diff(f, Y, 2) + 4*(2 - Y)*sp.diff(f, Y))
# paper's (67) nonlocal: -(kappa^2 H^2/8pi^2) 2H (aa')^2 (a d0 + a' d0') Lap A ; with A = h Ahat, H^2/8pi^2 = 2h -> units kappa^2 h^2:
LA = lap(F(A))
log67 = -4*H*(a*ap)**2*(a*sp.diff(LA, t) + ap*sp.diff(LA, tp))
lhs = direct + log67
x = Y/4; lx = sp.log(x); s2 = a**2 + ap**2
def B_(f, n): 
    for _ in range(n): f = boxY(f)
    return f
rows = {
 'a4B3': ((a*ap)**4/H**2, 3, 0, -lx/(3*x)),
 'a4B2': ((a*ap)**4, 2, 0, 26*lx/(3*x) + R(38, 3)/x - 6*lx**2 - 18*lx),
 'a4H2B': ((a*ap)**4*H**2, 1, 0, -6*lx/x + 4/x - 4*lx),
 'a4H4': ((a*ap)**4*H**4, 0, 0, 4*lx/x + 18/x - 120 - 108*lx),
 'a3sB3': ((a*ap)**3*s2/H**2, 3, 0, lx/(6*x)),
 'a3sB2': ((a*ap)**3*s2, 2, 0, -lx/(3*x) + 1/(6*x)),
 'a3sH2B': ((a*ap)**3*s2*H**2, 1, 0, -2*lx/(3*x) + 5/x - 18*lx),
 'a3sH4': ((a*ap)**3*s2*H**4, 0, 0, 4*lx/(3*x) - R(32, 3)/x - 54),
 'a3H2L': ((a*ap)**3*H**2, 0, 1, -2*lx/(3*x) - 16/x + 84*lx - 48*x*lx + 96*x),
 'a3LB': ((a*ap)**3, 1, 1, lx/(3*x)),
 'a2sH2L': ((a*ap)**2*s2*H**2, 0, 1, -7*lx/(3*x) + R(11, 2)/x - 12*lx + 48*x*lx + 12*x),
 'a2sLB': ((a*ap)**2*s2, 1, 1, -17*lx/(6*x) - R(49, 6)/x + 4*lx**2 + 10*lx + 12*x*lx),
 'a2L2': ((a*ap)**2, 0, 2, R(10, 3)*lx - 24*x*lx + 24*x**2*lx - 36*x**2)}
def apply(pref, nB, nL, f):
    g = F(B_(f, nB))
    for _ in range(nL): g = lap(g)
    return pref*g
paper_sum = sum(apply(*v) for v in rows.values())
fix_zeta = apply((a*ap)**2, 0, 2, 8*lx - 16*x*lx)
fix_72 = apply((a*ap)**3*s2*H**4, 0, 0, sp.Integer(72))
random.seed(29)
pts = []
while len(pts) < 3:
    pt = {H: random.choice([R(7, 10), R(13, 10), R(3, 2)]), t: -R(random.randint(3, 9), 4), tp: -R(random.randint(10, 19), 4)}
    for v in xs + ys: pt[v] = R(random.randint(-9, 9), 3)
    if sum((pt[xs[i]] - pt[ys[i]])**2 for i in range(3)) > (pt[t] - pt[tp])**2 + 1: pts.append(pt)
N = lambda e, pt: sp.N(e.subs(pt), 30)
for i, pt in enumerate(pts):
    L, P, Z, S = N(lhs, pt), N(paper_sum, pt), N(fix_zeta, pt), N(fix_72, pt)
    print(f"  point {i}: lhs - TableX(printed) = {sp.Float(L - P, 6)};  + zeta fix: {sp.Float(L - P - Z, 6)};"
          f"  + zeta + 72 gamma2: {sp.Float(L - P - Z - S, 6)};  (scale |lhs| = {sp.Float(abs(L), 4)})")
