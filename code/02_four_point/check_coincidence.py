"""check_coincidence.py -- Task 4: coincidence limits, Eqs. (I.34)-(I.40).

Purpose
    Supplementary checks of the analytic derivation of the coincidence limits of
    the scalar propagators iDelta_A, iDelta_B, iDelta_C and of the four
    contractions of the coincident graviton propagator used in Paper I
    (Kahya & Woodard, PRD 76, 124005 (2007)).  The analytic derivation is in notes
    chapter "The 4-point contribution: Eqs. (I.31)-(I.44)", section "Task 4".

Checks
    (1) Gauss's theorem 2F1(a,b;c;1) = G(c)G(c-a-b)/(G(c-a)G(c-b)) applied to the
        hypergeometric forms (I.29)-(I.30) reproduces (I.35) and (I.36) for symbolic D.
    (2) Numerical coincidence limits at D = 0.6 (Re D < 2, where the limit exists
        literally), H = 1, 30 digits: series (I.27)-(I.28) vs hypergeometric
        (I.29)-(I.30) at y = 0.5, and the approach of the series (I.14), (I.27),
        (I.28) to (I.34)-(I.36) for y = 1e-2 ... 1e-8 (differences ~ (y/4)^0.7).
    (3) Explicit-component construction of T^A, T^B, T^C (Eqs. (I.22)-(I.24)) for
        D = 4, 5, 6, 7 with iDelta_A, iDelta_B, iDelta_C symbolic: pair symmetries,
        the closed forms of (I.37)-(I.40) in every component, and the trace
        identities trace(I.40) = (I.38), eta.(I.39) = (I.37).
    (4) D = 4 reduction of the finite (B, C) parts (iDelta_B -> -h, iDelta_C -> +h,
        h = H^2/16pi^2) compared with the printed finite parts.

Method
    SymPy (symbolic D, explicit tensor components) and mpmath (30-digit series).

Expected result
    All True; finite parts (37): 4h, (38): -2h, (39): [2,2,0]h, (40): [-2,0,0]h
    (the paper prints +2h for the (alpha,sigma)=(0,0) entry of (I.40)).
    See expected_output/check_coincidence.txt.

Runtime
    about 2 s.

Provenance
    Reconstructed verbatim from the Task 4 command log (script
    build4/check_coincidence.py); only this docstring was added.
"""
# Supplementary checks for Eqs. (34)-(40): coincidence limits and graviton-propagator contractions.
import sympy as sp, mpmath as mp, itertools
D = sp.Symbol('D')
N = sp.Symbol('N')            # N = H^{D-2}/(4 pi)^{D/2}
G = sp.gamma

print("== (1) Gauss theorem on the hypergeometric forms (29)-(30): 2F1(a,b;c;1) = G(c)G(c-a-b)/(G(c-a)G(c-b)), Re(c-a-b)=1-D/2>0 ==")
gauss = lambda a, b, c: G(c)*G(c - a - b)/(G(c - a)*G(c - b))
B0 = N*G(D - 2)*G(1)/G(D/2)*gauss(D - 2, 1, D/2)
C0 = N*G(D - 3)*G(2)/G(D/2)*gauss(D - 3, 2, D/2)
print("  iDelta_B(x;x) via (29) == -N G(D-1)/((D-2) G(D/2)) :", sp.simplify(sp.gammasimp(B0 + N*G(D - 1)/((D - 2)*G(D/2)))) == 0)
print("  iDelta_C(x;x) via (30) ==  N G(D-1)/((D-2)(D-3) G(D/2)):", sp.simplify(sp.gammasimp(C0 - N*G(D - 1)/((D - 2)*(D - 3)*G(D/2)))) == 0)

print("== (2) numerical coincidence limits at D = 0.6 (Re D < 2: the limit exists literally), H = 1 ==")
mp.mp.dps = 30
Dn = mp.mpf('0.6'); Nn = 1/(4*mp.pi)**(Dn/2)
def A_series(y, nmax=400):
    s = mp.gamma(Dn/2 - 1)/(Dn/2 - 1)*(4/y)**(Dn/2 - 1) + mp.gamma(Dn/2 + 1)/(Dn/2 - 2)*(4/y)**(Dn/2 - 2) \
        - mp.pi*mp.cot(mp.pi*Dn/2)*mp.gamma(Dn - 1)/mp.gamma(Dn/2)
    for n in range(1, nmax):
        s += mp.gamma(n + Dn - 1)/(n*mp.gamma(n + Dn/2))*(y/4)**n \
             - mp.gamma(n + Dn/2 + 1)/((n - Dn/2 + 2)*mp.gamma(n + 2))*(y/4)**(n - Dn/2 + 2)
    return Nn*s
def cf(y): return Nn*mp.gamma(Dn/2 - 1)*(4/y)**(Dn/2 - 1)
def B_series(y, nmax=400):
    s = sum(mp.gamma(n + Dn - 2)/mp.gamma(n + Dn/2)*(y/4)**n - mp.gamma(n + Dn/2)/mp.gamma(n + 2)*(y/4)**(n - Dn/2 + 2) for n in range(nmax))
    return cf(y) - Nn*s
def C_series(y, nmax=400):
    s = sum((n + 1)*mp.gamma(n + Dn - 3)/mp.gamma(n + Dn/2)*(y/4)**n - (n - Dn/2 + 3)*mp.gamma(n + Dn/2 - 1)/mp.gamma(n + 2)*(y/4)**(n - Dn/2 + 2) for n in range(nmax))
    return cf(y) + Nn*s
B_hyp = lambda y: Nn*mp.gamma(Dn - 2)/mp.gamma(Dn/2)*mp.hyp2f1(Dn - 2, 1, Dn/2, 1 - y/4)
C_hyp = lambda y: Nn*mp.gamma(Dn - 3)*mp.gamma(2)/mp.gamma(Dn/2)*mp.hyp2f1(Dn - 3, 2, Dn/2, 1 - y/4)
k = Nn*mp.gamma(Dn - 1)/mp.gamma(Dn/2)
A34 = k*(-mp.pi*mp.cot(mp.pi*Dn/2)); B35 = -k/(Dn - 2); C36 = k/((Dn - 2)*(Dn - 3))
print("  series vs hypergeometric at y=0.5:  B:", mp.nstr(B_series(mp.mpf('0.5')) - B_hyp(mp.mpf('0.5')), 5),
      "  C:", mp.nstr(C_series(mp.mpf('0.5')) - C_hyp(mp.mpf('0.5')), 5))
for yv in ('1e-2', '1e-4', '1e-6', '1e-8'):
    y = mp.mpf(yv)
    print(f"  y={yv:>5}:  A(y)-(34) = {mp.nstr(A_series(y) - A34, 4):>11}   B(y)-(35) = {mp.nstr(B_series(y) - B35, 4):>11}   C(y)-(36) = {mp.nstr(C_series(y) - C36, 4):>11}")
print("  (differences scale like (y/4)^(1-D/2) = (y/4)^0.7 -> 0)")

print("== (3) explicit-component contractions of T^A, T^B, T^C, D = 4,5,6,7 ==")
a_, b_, c_ = sp.symbols('iDA iDB iDC')
def closed(Dv):
    d = sp.Integer(Dv)
    eb = lambda r, s: (1 if r == s and r > 0 else 0)           # eta-bar^{rs} = eta-bar_{rs} numerically
    db = lambda s, a: (1 if s == a and a > 0 else 0)          # delta-bar^s_a
    d0 = lambda r: (1 if r == 0 else 0)
    c37 = -4*(d - 1)/(d - 3)*a_ + 8/((d - 2)*(d - 3))*c_
    c38 = (d - 1)*(d**2 - 3*d - 2)/(d - 3)*a_ + 2*(d - 1)*b_ + 2*(d**2 - 5*d + 8)/((d - 2)*(d - 3))*c_
    c39 = lambda r, s: -4/(d - 3)*eb(r, s)*a_ + 4/((d - 2)*(d - 3))*((d - 3)*d0(r)*d0(s) + eb(r, s))*c_
    c40 = lambda a, s: (d**2 - 3*d - 2)/(d - 3)*db(s, a)*a_ + ((d - 1)*d0(a)*d0(s) + db(s, a))*b_ \
                       + 2/((d - 2)*(d - 3))*((d - 3)**2*d0(a)*d0(s) + db(s, a))*c_
    return c37, c38, c39, c40
for Dv in (4, 5, 6, 7):
    d = sp.Integer(Dv); R = range(Dv)
    eta = sp.diag(*([-1] + [1]*(Dv - 1))); ei = eta.inv()
    u = [1 if m == 0 else 0 for m in R]                       # u_mu = delta^0_mu
    eb = [[eta[m, n] + u[m]*u[n] for n in R] for m in R]      # eta-bar_{mu nu}
    sym = lambda f: sp.Rational(1, 4)*(f(0) + f(1) + f(2) + f(3))
    TA = lambda m, n, r, s: eb[m][r]*eb[s][n] + eb[m][s]*eb[r][n] - 2/(d - 3)*eb[m][n]*eb[r][s]
    TB = lambda m, n, r, s: -(u[m]*eb[n][r]*u[s] + u[n]*eb[m][r]*u[s] + u[m]*eb[n][s]*u[r] + u[n]*eb[m][s]*u[r])
    f = lambda m, n: (d - 3)*u[m]*u[n] + eb[m][n]
    TC = lambda m, n, r, s: 2/((d - 2)*(d - 3))*f(m, n)*f(r, s)
    X = {}
    for m, n, r, s in itertools.product(R, R, R, R):
        X[(m, n, r, s)] = TA(m, n, r, s)*a_ + TB(m, n, r, s)*b_ + TC(m, n, r, s)*c_
    pair_ok = all(sp.expand(X[(m, n, r, s)] - X[(r, s, m, n)]) == 0 and sp.expand(X[(m, n, r, s)] - X[(n, m, r, s)]) == 0
                  for m, n, r, s in itertools.product(R, R, R, R))
    x37 = sum(ei[m, n]*ei[r, s]*X[(m, n, r, s)] for m, n, r, s in itertools.product(R, R, R, R))
    x38 = sum(ei[m, r]*ei[n, s]*X[(m, n, r, s)] for m, n, r, s in itertools.product(R, R, R, R))
    x39 = lambda P, Q: sum(ei[m, n]*ei[P, r]*ei[Q, s]*X[(m, n, r, s)] for m, n, r, s in itertools.product(R, R, R, R))
    x40 = lambda A, S: sum(ei[r, l]*ei[S, kk]*X[(A, r, l, kk)] for r, l, kk in itertools.product(R, R, R))
    c37, c38, c39, c40 = closed(Dv)
    ok37 = sp.simplify(x37 - c37) == 0
    ok38 = sp.simplify(x38 - c38) == 0
    ok39 = all(sp.simplify(x39(P, Q) - c39(P, Q)) == 0 for P in R for Q in R)
    ok40 = all(sp.simplify(x40(A, S) - c40(A, S)) == 0 for A in R for S in R)
    tr40 = sp.simplify(sum(x40(A, A) for A in R) - x38) == 0
    tr39 = sp.simplify(sum(eta[P, Q]*x39(P, Q) for P in R for Q in R) - x37) == 0
    print(f"  D={Dv}: pair/within-pair symmetric {pair_ok}; (37) {ok37}; (38) {ok38}; (39) {ok39}; (40) {ok40}; "
          f"trace(40)=(38) {tr40}; eta.(39)=(37) {tr39}")

print("== (4) D = 4 reduction of the finite (B,C) parts, h = H^2/(16 pi^2): iDB -> -h, iDC -> +h ==")
h = sp.Symbol('h')
c37, c38, c39, c40 = closed(4)
sub = {b_: -h, c_: h}
print("  (37) finite:", sp.simplify(c37.subs(sub).coeff(h)), "* h      paper: +4 h")
print("  (38) finite:", sp.simplify(c38.subs(sub).coeff(h)), "* h      paper: -2 h")
print("  (39) finite [rho,sigma]=(0,0),(1,1),(0,1):", [sp.simplify(c39(P, Q).subs(sub).coeff(h)) for P, Q in ((0, 0), (1, 1), (0, 1))],
      "* h   paper: [2, 2, 0]")
print("  (40) finite [alpha,sigma]=(0,0),(1,1),(0,1):", [sp.simplify(c40(A, S).subs(sub).coeff(h)) for A, S in ((0, 0), (1, 1), (0, 1))],
      "* h   paper: [+2, 0, 0]")
print("  trace of paper's (40) finite part = +2h, but paper's (38) finite part = -2h  -> paper's (40) inconsistent with its (38)")
