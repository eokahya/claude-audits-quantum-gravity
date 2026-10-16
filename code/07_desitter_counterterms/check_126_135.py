"""
check_126_135.py -- code checks for Task 14 (Paper I, Eqs. (I.126)-(I.135)).

Paper I: E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

Purpose
    Independent symbolic (SymPy) checks of the analytic derivation of Task 14:
    the invariant counterterms, the decomposition of the D(D+1)/2-dimensional
    de Sitter group in conformal coordinates into four classes (spatial
    translations, rotations, dilatation, spatial special conformal
    transformations), and the residual symmetry of the gauge-fixing functional
    (I.15).  The script supplements, and does not replace, the derivation in
    the notes (chapter "de Sitter symmetries, the gauge, and noninvariant
    counterterms").

What is checked (section labels as printed)
    (A) D=4,5 de Sitter background a = -1/(H eta): Christoffel symbols
        Gamma^l_mn = H a (d^l_m d^0_n + d^l_n d^0_m + eta_mn d^l_0),
        R_mn = (D-1) H^2 g_mn, R = D(D-1) H^2.
    (B) Commutator identity used in the list of invariants (D=4 de Sitter,
        generic phi): phi;mn phi;^mn - (Box phi)^2 + R_mn phi^;m phi^;n is the
        total divergence nabla_m(phi^;n phi;_n^;m - phi^;m Box phi).
    (C) D=4: the finite maps (I.132)-(I.135) preserve ds^2; time translation
        does not.
    (D) The general Killing vector
            xi^0 = eta (lam + 2 b.x),
            xi^i = c^i + w^i_j x^j + lam x^i + 2 (b.x) x^i - (x.x) b^i,
        x.x = |x|^2 - eta^2, solves the Killing equation for D=3,4,5 and has
        D(D+1)/2 parameters; time translation, boost and temporal special
        conformal map fail.  A polynomial ansatz of degree
        <= 3 gives a solution space of dimension 6 (D=3) and 10 (D=4).
    (E) The infinitesimal limit of (I.135) is the b-part of the general
        Killing vector in (D); for eta -> 0
        (I.135) reduces to the Euclidean special conformal map of R^{D-1}.
    (F) D=4, generic h_mn: F_mu = G_mu + E_mu (covariant de Donder vector of
        psi = a^2 h plus E_mu = H(2 psi_mn u^n - u_mu psi), u^mu = d^mu_0 / a);
        a^{D-2} eta^mn = sqrt(-g) g^mn.
    (G) D=4: [xi, u] = (0, 2 H eta^2 b^i) for the general Killing vector.
    (H) D=3, generic h_mn(eta, x1, x2): delta L_GF - d_l(xi^l L_GF) = 0 for
        translations, rotations and dilatation; for the spatial special
        conformal map it equals -a^{D-2} eta^mn F_m Delta_n with
        Delta_m = -(4 h_mi b^i - 2 b_m h), and its Euler-Lagrange derivatives
        with respect to h_00 and h_12 are nonzero (not a boundary term).
    (I) D=4 Lie algebra: [Dil,P]=-P, [Dil,K]=K, [K_i,P_j] = -2(d_ij Dil + M_ij),
        [K,K]=[P,P]=0, [M_12,K_1]=-K_2, [M_12,P_1]=-P_2.

Method
    Exact symbolic algebra with sympy (generic functions phi(x), h_mn(x)).
    Coordinates: t stands for conformal time eta < 0.

Expected result
    Every test prints True (and the expected structures); "False" appears only
    for the contrast test "time translation", which is not an isometry.
    See expected_output/check_126_135.txt.

Runtime
    About 10 s (Python 3.13, sympy 1.14).

Run
    python3 check_126_135.py
"""
import itertools

import sympy as sp

R = sp.Rational
H = sp.Symbol('H', positive=True)

def coords(Dv):
    t = sp.Symbol('t', negative=True)            # conformal time eta < 0
    xs = sp.symbols(f'x1:{Dv}', real=True)
    return (t,) + xs

def bg(Dv):
    X = coords(Dv); t = X[0]
    a = -1/(H*t)
    eta = sp.diag(*([-1] + [1]*(Dv - 1)))
    g = a**2*eta
    return X, a, eta, g

def christoffel(g, X):
    n = len(X); gi = g.inv()
    return [[[sp.simplify(sum(gi[l, k]*(sp.diff(g[k, m], X[nn]) + sp.diff(g[k, nn], X[m]) - sp.diff(g[m, nn], X[k])) for k in range(n))/2)
              for nn in range(n)] for m in range(n)] for l in range(n)]

print("== (A) background: Christoffels, R_mn = (D-1)H^2 g_mn, R = D(D-1)H^2 ==")
for Dv in (4, 5):
    X, a, eta, g = bg(Dv); n = Dv
    G = christoffel(g, X)
    d0 = lambda m: 1 if m == 0 else 0
    okG = all(sp.simplify(G[l][m][nn] - H*a*(int(l == m)*d0(nn) + int(l == nn)*d0(m) + eta[m, nn]*d0(l))) == 0
              for l in range(n) for m in range(n) for nn in range(n))
    Ric = sp.zeros(n, n)
    for m in range(n):
        for nn in range(n):
            Ric[m, nn] = sp.simplify(sum(sp.diff(G[l][m][nn], X[l]) - sp.diff(G[l][m][l], X[nn])
                                         + sum(G[l][l][k]*G[k][m][nn] - G[l][nn][k]*G[k][m][l] for k in range(n)) for l in range(n)))
    okR = sp.simplify(Ric - (Dv - 1)*H**2*g) == sp.zeros(n, n)
    Rs = sp.simplify(sum(g.inv()[m, nn]*Ric[m, nn] for m in range(n) for nn in range(n)))
    print(f"  D={Dv}: Gamma = Ha(d^l_m d^0_n + d^l_n d^0_m + eta_mn d^l_0): {okG};  R_mn = (D-1)H^2 g_mn: {okR};  R = {sp.factor(Rs)}")

print("== (B) phi;mn phi;^mn - (Box phi)^2 + R_mn phi^;m phi^;n is a total divergence (D=4 dS, generic phi) ==")
Dv = 4; X, a, eta, g = bg(Dv); n = Dv; gi = g.inv(); G = christoffel(g, X)
phi = sp.Function('phi')(*X); sqg = a**Dv
dphi = [sp.diff(phi, X[m]) for m in range(n)]
hess = [[sp.diff(phi, X[m], X[nn]) - sum(G[l][m][nn]*dphi[l] for l in range(n)) for nn in range(n)] for m in range(n)]
up = lambda v: [sum(gi[m, k]*v[k] for k in range(n)) for m in range(n)]
dU = up(dphi)
boxphi = sum(gi[m, nn]*hess[m][nn] for m in range(n) for nn in range(n))
hess_ud = [[sum(gi[m, k]*hess[k][nn] for k in range(n)) for nn in range(n)] for m in range(n)]   # phi;^m_n
lhs = sum(hess_ud[m][nn]*hess_ud[nn][m] for m in range(n) for nn in range(n)) - boxphi**2 + (Dv - 1)*H**2*sum(g[m, nn]*dU[m]*dU[nn] for m in range(n) for nn in range(n))
J = [sum(dU[nn]*hess_ud[m][nn] for nn in range(n)) - dU[m]*boxphi for m in range(n)]     # J^m = phi^;n phi;_n^m - phi^;m Box phi
div = sum(sp.diff(sqg*J[m], X[m]) for m in range(n))/sqg
print("  identity holds:", sp.simplify(sp.expand(lhs - div)) == 0)

print("== (C) finite transformations (132)-(135) are isometries of ds^2 = (-d eta^2 + dx^2)/(H^2 eta^2), D=4 ==")
Dv = 4; X, a, eta, g = bg(Dv); t, x1, x2, x3 = X
def pullback_ok(newX):
    Jm = sp.Matrix(Dv, Dv, lambda m, nn: sp.diff(newX[m], X[nn]))
    tnew = newX[0]
    gnew = (1/(H*tnew))**2*eta
    return sp.simplify(Jm.T*gnew*Jm - g) == sp.zeros(Dv, Dv)
e1, e2, e3, k, th1, th2, th3 = sp.symbols('e1 e2 e3 k theta1 theta2 theta3', real=True)
ang = sp.Symbol('phi_r', real=True)
Rot = sp.Matrix([[sp.cos(ang), -sp.sin(ang), 0], [sp.sin(ang), sp.cos(ang), 0], [0, 0, 1]])
xv = sp.Matrix([x1, x2, x3])
print("  (132) translation:", pullback_ok([t, x1 + e1, x2 + e2, x3 + e3]))
print("  (133) rotation   :", pullback_ok([t] + list(Rot*xv)))
print("  (134) dilatation :", pullback_ok([k*t, k*x1, k*x2, k*x3]))
xx = -t**2 + x1**2 + x2**2 + x3**2; th = sp.Matrix([th1, th2, th3]); thx = th.dot(xv); th2n = th.dot(th)
Om = 1 - 2*thx + th2n*xx
print("  (135) spatial SCT:", pullback_ok([t/Om] + [(xv[i] - th[i]*xx)/Om for i in range(3)]))
print("  for contrast, time translation t -> t + c :", pullback_ok([t + e1, x1, x2, x3]))

print("== (D) Killing vectors: general solution, failures, parameter count ==")
def lie_g(xi, X, g):
    n = len(X)
    return sp.Matrix(n, n, lambda m, nn: sp.simplify(sum(xi[l]*sp.diff(g[m, nn], X[l]) + g[l, nn]*sp.diff(xi[l], X[m]) + g[m, l]*sp.diff(xi[l], X[nn]) for l in range(n))))
for Dv in (3, 4, 5):
    X, a, eta, g = bg(Dv); t = X[0]; xs = X[1:]; m = Dv - 1
    c = sp.symbols(f'c1:{Dv}'); bb = sp.symbols(f'b1:{Dv}'); lam = sp.Symbol('lam')
    om = sp.zeros(m, m); oms = []
    for i in range(m):
        for j in range(i + 1, m):
            s = sp.Symbol(f'w{i}{j}'); om[i, j] = s; om[j, i] = -s; oms.append(s)
    r2 = sum(v**2 for v in xs); bx = sum(bb[i]*xs[i] for i in range(m))
    xi = [t*(lam + 2*bx)] + [c[i] + sum(om[i, j]*xs[j] for j in range(m)) + lam*xs[i] + 2*bx*xs[i] - (r2 - t**2)*bb[i] for i in range(m)]
    ok = lie_g(xi, X, g) == sp.zeros(Dv, Dv)
    npar = len(c) + len(oms) + 1 + len(bb)
    # failures among flat conformal Killing vectors
    a0, w0, b0 = sp.symbols('a0 w0 b0')
    tt = [a0] + [0]*m                                  # time translation
    boost = [w0*xs[0]] + [w0*t] + [0]*(m - 1)          # boost in (t, x1)
    xxm = -t**2 + r2
    tsct = [2*(-b0*t)*t - xxm*b0] + [2*(-b0*t)*xs[i] for i in range(m)]   # flat SCT with b^mu = (b0, 0): b.x = -b0 t
    fails = [lie_g(v, X, g) != sp.zeros(Dv, Dv) for v in (tt, boost, tsct)]
    print(f"  D={Dv}: general xi solves Killing eq: {ok};  parameters = {npar} = D(D+1)/2 = {Dv*(Dv + 1)//2};"
          f"  time translation / boost / temporal SCT fail: {fails}")
# completeness by polynomial ansatz (degree <= 3) for D=3,4
for Dv in (3, 4):
    X, a, eta, g = bg(Dv)
    monos = [mm for mm in itertools.product(range(4), repeat=Dv) if sum(mm) <= 3]
    coeffs = []; xi = []
    for mu in range(Dv):
        comp = 0
        for mm in monos:
            s = sp.Symbol(f'q{mu}_' + '_'.join(map(str, mm))); coeffs.append(s)
            comp += s*sp.prod([X[i]**mm[i] for i in range(Dv)])
        xi.append(comp)
    L = lie_g(xi, X, g)
    eqs = []
    for e in L:
        num = sp.numer(sp.together(sp.expand(e)))
        eqs += sp.Poly(sp.expand(num), *X).coeffs()
    Mx = sp.Matrix([[sp.diff(e, s) for s in coeffs] for e in eqs])
    print(f"  D={Dv}: polynomial ansatz deg<=3 -> dimension of Killing solution space = {len(coeffs) - Mx.rank()}")

print("== (E) infinitesimal limit of (135) and boundary limit eta -> 0 ==")
Dv = 4; X, a, eta, g = bg(Dv); t, x1, x2, x3 = X; s_ = sp.Symbol('s')
b1, b2, b3 = sp.symbols('b1 b2 b3')
fin = [t/Om] + [(xv[i] - th[i]*xx)/Om for i in range(3)]
fin = [f.subs({th1: s_*b1, th2: s_*b2, th3: s_*b3}) for f in fin]
gen = [sp.simplify(sp.diff(f, s_).subs(s_, 0)) for f in fin]
bx = b1*x1 + b2*x2 + b3*x3; bv = [b1, b2, b3]
expect = [2*bx*t] + [2*bx*xv[i] - xx*bv[i] for i in range(3)]
print("  d/ds (135)|_{theta = s b} == (2 b.x eta, 2 (b.x) x^i - (x.x) b^i):", all(sp.simplify(gen[i] - expect[i]) == 0 for i in range(4)))
bnd = [sp.limit(f, t, 0) for f in [(xv[i] - th[i]*xx)/Om for i in range(3)]]
r2 = x1**2 + x2**2 + x3**2
print("  eta -> 0: x'^i = (x^i - theta^i |x|^2)/(1 - 2 theta.x + |theta|^2 |x|^2)  (Euclidean SCT of R^{D-1}):",
      all(sp.simplify(bnd[i] - (xv[i] - th[i]*r2)/(1 - 2*thx + th2n*r2)) == 0 for i in range(3)))

print("== (F) gauge functional: F_mu = G_mu + E_mu, G = covariant de Donder of psi = a^2 h, E = H(2 psi_mn u^n - u_m psi), u^m = d^m_0/a ==")
Dv = 4; X, a, eta, g = bg(Dv); n = Dv; gi = g.inv(); G = christoffel(g, X)
h = sp.zeros(n, n)
for i in range(n):
    for j in range(i, n):
        f = sp.Function(f'h{i}{j}')(*X); h[i, j] = f; h[j, i] = f
psi = a**2*h
etai = eta.inv()
htr = sum(etai[r, s]*h[r, s] for r in range(n) for s in range(n))
F = [sum(etai[r, s]*(sp.diff(h[m, r], X[s]) - R(1, 2)*sp.diff(h[r, s], X[m]) + (Dv - 2)*H*a*h[m, r]*int(s == 0)) for r in range(n) for s in range(n)) for m in range(n)]
def cov_psi(s_, m_, r_):   # nabla_s psi_mr
    return sp.diff(psi[m_, r_], X[s_]) - sum(G[l][s_][m_]*psi[l, r_] + G[l][s_][r_]*psi[m_, l] for l in range(n))
Gdd = [sum(gi[r, s]*(cov_psi(s, m, r) - R(1, 2)*cov_psi(m, r, s)) for r in range(n) for s in range(n)) for m in range(n)]
u_up = [1/a if m == 0 else 0 for m in range(n)]
u_dn = [sum(g[m, k]*u_up[k] for k in range(n)) for m in range(n)]
psitr = sum(gi[r, s]*psi[r, s] for r in range(n) for s in range(n))
E = [H*(2*sum(psi[m, k]*u_up[k] for k in range(n)) - u_dn[m]*psitr) for m in range(n)]
print("  F_mu - G_mu - E_mu = 0 :", all(sp.simplify(sp.expand(F[m] - Gdd[m] - E[m])) == 0 for m in range(n)))
print("  L_GF = -1/2 a^{D-2} eta^mn F F = -1/2 sqrt(-g) g^mn F F :", sp.simplify(a**(Dv - 2)*etai - a**Dv*gi) == sp.zeros(n, n))

print("== (G) [xi, u] for the general Killing vector, D=4 ==")
t, x1, x2, x3 = X; xs = X[1:]
c = sp.symbols('c1:4'); bb = sp.symbols('b1:4'); lam = sp.Symbol('lam'); w12, w13, w23 = sp.symbols('w12 w13 w23')
om = sp.Matrix([[0, w12, w13], [-w12, 0, w23], [-w13, -w23, 0]])
r2 = x1**2 + x2**2 + x3**2; bx = sum(bb[i]*xs[i] for i in range(3))
xi = [t*(lam + 2*bx)] + [c[i] + sum(om[i, j]*xs[j] for j in range(3)) + lam*xs[i] + 2*bx*xs[i] - (r2 - t**2)*bb[i] for i in range(3)]
comm = [sp.simplify(sum(xi[k]*sp.diff(u_up[m], X[k]) - u_up[k]*sp.diff(xi[m], X[k]) for k in range(n))) for m in range(n)]
print("  [xi,u] =", comm, "  (= (0, 2 H eta^2 b^i): vanishes iff b = 0)")

print("== (H) variation of L_GF under each generator, D=3 (generic h_mn(eta,x1,x2)) ==")
Dv = 3; X, a, eta, g = bg(Dv); n = Dv; etai = eta.inv(); t, x1, x2 = X
h = sp.zeros(n, n)
for i in range(n):
    for j in range(i, n):
        f = sp.Function(f'h{i}{j}')(*X); h[i, j] = f; h[j, i] = f
def Fof(hh):
    return [sum(etai[r, s]*(sp.diff(hh[m, r], X[s]) - R(1, 2)*sp.diff(hh[r, s], X[m]) + (Dv - 2)*H*a*hh[m, r]*int(s == 0)) for r in range(n) for s in range(n)) for m in range(n)]
Fh = Fof(h)
L = -R(1, 2)*a**(Dv - 2)*sum(etai[m, nn]*Fh[m]*Fh[nn] for m in range(n) for nn in range(n))
def delta_h(xi):   # active: delta psi = Lie_xi psi, psi = a^2 h
    psi = a**2*h
    lp = sp.Matrix(n, n, lambda m, nn: sum(xi[l]*sp.diff(psi[m, nn], X[l]) + psi[l, nn]*sp.diff(xi[l], X[m]) + psi[m, l]*sp.diff(xi[l], X[nn]) for l in range(n)))
    return lp/a**2
c1_, c2_, w_, lam_, b1_, b2_ = sp.symbols('c1 c2 w lam b1 b2')
bxv = b1_*x1 + b2_*x2; r2v = x1**2 + x2**2
gens = {'translation': [0, c1_, c2_], 'rotation': [0, w_*x2, -w_*x1], 'dilatation': [lam_*t, lam_*x1, lam_*x2],
        'spatial SCT': [2*bxv*t, 2*bxv*x1 - (r2v - t**2)*b1_, 2*bxv*x2 - (r2v - t**2)*b2_]}
htr = sum(etai[r, s]*h[r, s] for r in range(n) for s in range(n))
for name, xi in gens.items():
    dh = delta_h(xi)
    dF = Fof(dh)
    dL = -a**(Dv - 2)*sum(etai[m, nn]*Fh[m]*dF[nn] for m in range(n) for nn in range(n))
    anom = sp.simplify(sp.expand(dL - sum(sp.diff(xi[l]*L, X[l]) for l in range(n))))
    if anom == 0:
        print(f"  {name:12s}: delta L_GF - d_l(xi^l L_GF) = 0   -> symmetry of the gauge-fixing action")
    else:
        bvec = [0, b1_, b2_]
        Delta = [-(4*sum(h[m, i]*bvec[i] for i in range(1, n)) - 2*bvec[m]*htr) for m in range(n)]
        pred = -a**(Dv - 2)*sum(etai[m, nn]*Fh[m]*Delta[nn] for m in range(n) for nn in range(n))
        same = sp.simplify(sp.expand(anom - pred)) == 0
        # Euler-Lagrange derivative of the anomaly density (first-order Lagrangian pred) w.r.t. h_00
        ELs = []
        for (i, j) in [(0, 0), (1, 2)]:
            hf = h[i, j]
            dsym = sp.Symbol('HF'); ders = [sp.Symbol(f'HF{k}') for k in range(n)]
            P = pred.subs({sp.diff(hf, X[k]): ders[k] for k in range(n)}).subs(hf, dsym)
            el = sp.diff(P, dsym) - sum(sp.diff(sp.diff(P, ders[k]).subs({dsym: hf, **{ders[q]: sp.diff(hf, X[q]) for q in range(n)}}), X[k]) for k in range(n))
            el = el.subs({dsym: hf, **{ders[q]: sp.diff(hf, X[q]) for q in range(n)}})
            ELs.append(sp.simplify(sp.expand(el)) != 0)
        print(f"  {name:12s}: anomaly != 0;  anomaly == -a^(D-2) eta^mn F_m Delta_n with Delta_m = -(4 h_mi b^i - 2 b_m h): {same};"
              f"  Euler-Lagrange of anomaly w.r.t. h_00, h_12 nonzero: {ELs}")

print("== (I) Lie algebra of the Killing vectors, D=4: P_i, M_ij, Dil, K_i ==")
Dv = 4; X = coords(Dv); t = X[0]; xs = X[1:]; n = Dv
def P(i): return [0] + [1 if j == i else 0 for j in range(3)]
def M(i, j): return [0] + [xs[i] if k == j else (-xs[j] if k == i else 0) for k in range(3)]
Dil = [t] + list(xs)
r2 = sum(v**2 for v in xs)
def K(i): return [2*xs[i]*t] + [2*xs[i]*xs[k] - (r2 - t**2)*(1 if k == i else 0) for k in range(3)]
def br(A, B): return [sp.expand(sum(A[k]*sp.diff(B[m], X[k]) - B[k]*sp.diff(A[m], X[k]) for k in range(n))) for m in range(n)]
def eq(A, B): return all(sp.expand(A[m] - B[m]) == 0 for m in range(n))
lin = lambda *terms: [sp.expand(sum(c*v[m] for c, v in terms)) for m in range(n)]
print("  [Dil, P_1] = -P_1        :", eq(br(Dil, P(0)), lin((-1, P(0)))))
print("  [Dil, K_1] = +K_1        :", eq(br(Dil, K(0)), lin((1, K(0)))))
print("  [K_1, P_1] = -2 Dil      :", eq(br(K(0), P(0)), lin((-2, Dil))))
print("  [K_1, P_2] = -2 M_12     :", eq(br(K(0), P(1)), lin((-2, M(0, 1)))))
print("  [K_1, K_2] = 0, [P_1,P_2] = 0:", eq(br(K(0), K(1)), [0]*n), eq(br(P(0), P(1)), [0]*n))
print("  [M_12, K_1] = -K_2, [M_12, P_1] = -P_2:", eq(br(M(0, 1), K(0)), lin((-1, K(1)))), eq(br(M(0, 1), P(0)), lin((-1, P(1)))))
