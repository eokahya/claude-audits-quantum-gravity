"""check_41_44.py -- Task 5: the final 4-point contribution, Eqs. (I.41)-(I.44).

Purpose
    Supplementary checks of the analytic derivation of (I.41), (I.43) and (I.44)
    of Paper I (Kahya & Woodard, PRD 76, 124005 (2007)) from the derived (I.33')
    (Task 3) and the derived coincidence limits (Task 4).  The analytic derivation
    is in notes chapter "The 4-point contribution: Eqs. (I.31)-(I.44)",
    section "Task 5".

Checks
    (A) Symbolic-D bookkeeping of (I.33) in the basis P[g] = d_0(g d_0 delta),
        S[g] = g Lap delta: coefficients of the iDelta_A part and of the finite
        (B, C) part, which is O(D-4); reproduction of the paper's chain
        ((I.33) as printed + printed (I.40)) giving the extra bracket of (I.41).
    (B) Explicit D = 4 tensors T^A, T^B, T^C, a = -1/(H t), Lambda = k(pihat + 2 ln a)
        and a test function f: the kernel equals -i kappa^2 {3 d^mu(a^2 Lambda d_mu)
        - 4 a^2 Lambda Lap} f; the B/C part alone (c = -b) gives zero.
    (C) Eq. (I.43) as an identity of kernels acting on test functions, in D = 4 and
        general D; false if ln(aa') is replaced by ln(a^2);
        Box ln a = -(D-1) H^2.
    (D) (I.41) -> (I.44) on a test function in D = 4: derived (I.41') == derived
        (I.44'); printed (I.41) == printed (I.44).

Method
    SymPy.

Expected result
    Coefficients I*D*(D-1)/4 and -I*D*(D-5)/4; finite parts proportional to (D-4);
    paper chain: 4*I*h, 0; all identity checks True except the deliberate
    ln(a^2) counter-examples (False).  See expected_output/check_41_44.txt.

Runtime
    about 1 s.

Provenance
    Reconstructed verbatim from the Task 5 command log (script
    build5/check_41_44.py, including the two in-log clean-up edits: the
    lhs43 line and the removal of an unused helper); only this docstring was added.
"""
# Supplementary checks for Eqs. (41)-(44).
import sympy as sp, itertools
I, R = sp.I, sp.Rational
D, kap, b, H, k, pih = sp.symbols('D kappa b H k pihat')   # pihat = -pi*cot(pi D/2)

print("== (A) symbolic-D bookkeeping of (33) in the basis P[g]=d_0(g d_0 delta), S[g]=g Lap delta ==")
# exact-D coincidence contractions (Task 4): A-coefficients and finite parts, c = -b/(D-3)
c = -b/(D - 3)
A1, A2 = -4*(D - 1)/(D - 3), (D - 1)*(D**2 - 3*D - 2)/(D - 3)
A3, A4 = -4/(D - 3), (D**2 - 3*D - 2)/(D - 3)                 # times etabar^{rs} / deltabar^s_a
F1 = 8*c/((D - 2)*(D - 3))
F2 = 2*(D - 1)*b + 2*(D**2 - 5*D + 8)*c/((D - 2)*(D - 3))
F3_uu, F3_eb = 4*c/(D - 2), 4*c/((D - 2)*(D - 3))            # coefficients of d^r_0 d^s_0 and etabar^{rs}
F4_00, F4_db = (D - 1)*b + 2*(D - 3)*c/(D - 2), b + 2*c/((D - 2)*(D - 3))   # of d^0_a d^s_0 and deltabar^s_a
def kernel(s1, s2, s3, s4, F4_00v=F4_00):
    """Coefficients of P[a^{D-2}Lam], S[a^{D-2}Lam], P[a^{D-2}], S[a^{D-2}] for
       s1 d^mu{a X1 d_mu} + s2 d^mu{a X2 d_mu} + s3 d_r{a X3^{rs} d_s} + s4 d^a{a X4_a^s d_s}  (times kappa^2)."""
    PL = s1*A1*(-1) + s2*A2*(-1)               # d^mu(g d_mu) = -P + S
    SL = s1*A1 + s2*A2 + s3*A3 + s4*A4          # etabar -> S, deltabar -> S
    P1 = s1*F1*(-1) + s2*F2*(-1) + s3*F3_uu + s4*F4_00v*(-1)   # d^0_a d^s_0 -> d^0(...d_0) = -P
    S1 = s1*F1 + s2*F2 + s3*F3_eb + s4*F4_db
    return [sp.factor(sp.simplify(e)) for e in (PL, SL, P1, S1)]
paper = (-I/8, I/4, I/2, -I)                    # paper's (33) signs
derived = tuple(-s for s in paper)              # Task 3: derived (33) = -(paper's 33)
PL, SL, P1, S1 = kernel(*derived)
print("  derived: coeff of P[a^{D-2}Lam] =", PL, " ; S[a^{D-2}Lam] =", SL)
print("           [expected  +i D(D-1)/4  and  -i D(D-1)/4 + i D, i.e. -i{D(D-1)/4 d^mu(.d_mu) - D Lap}]")
print("  derived: finite part  P[a^{D-2}] =", P1, " ; S[a^{D-2}] =", S1, "   (both O(D-4))")
# paper's chain at D=4 with printed (40): F4_00 -> +2h instead of the derived value; h = -b at D=4
h = sp.Symbol('h')
PLp, SLp, P1p, S1p = kernel(*paper, F4_00v=2*h)
P1p4 = sp.simplify(P1p.subs(D, 4).subs(b, -h)); S1p4 = sp.simplify(S1p.subs(D, 4).subs(b, -h))
print("  paper chain (33)+printed(40), D=4: finite part P[a^2] =", P1p4, " ; S[a^2] =", S1p4,
      "   -> (i kappa^2 H^2/4pi^2){-d^mu(a^2 d_mu) + a^2 Lap} = 4 i kappa^2 h P[a^2]")
P1d4 = sp.simplify(kernel(*paper)[2].subs(D, 4)); S1d4 = sp.simplify(kernel(*paper)[3].subs(D, 4))
print("  paper signs + derived (40), D=4: finite part P =", P1d4, " ; S =", S1d4)

print("== (B) explicit D=4 tensors and functions: derived kernel vs derived (41) on a test function ==")
t, x1, x2, x3 = X = sp.symbols('t x1 x2 x3', real=True)
a = -1/(H*t)
Lam = k*(pih + 2*sp.log(a))
bb, cc = sp.symbols('bB cC')
f = sp.Function('f')(*X)
Dv = 4; d = sp.Integer(Dv); Rg = range(Dv)
eta = sp.diag(-1, 1, 1, 1); ei = eta.inv()
u = [1, 0, 0, 0]
eb = [[eta[m, n] + u[m]*u[n] for n in Rg] for m in Rg]
TA = lambda m, n, r, s: eb[m][r]*eb[s][n] + eb[m][s]*eb[r][n] - 2/(d - 3)*eb[m][n]*eb[r][s]
TB = lambda m, n, r, s: -(u[m]*eb[n][r]*u[s] + u[n]*eb[m][r]*u[s] + u[m]*eb[n][s]*u[r] + u[n]*eb[m][s]*u[r])
fT = lambda m, n: (d - 3)*u[m]*u[n] + eb[m][n]
TC = lambda m, n, r, s: 2/((d - 2)*(d - 3))*fT(m, n)*fT(r, s)
def Xc(m, n, r, s, LamV, bV, cV): return TA(m, n, r, s)*LamV + TB(m, n, r, s)*bV + TC(m, n, r, s)*cV
def Kt(LamV, bV, cV):
    K = sp.zeros(Dv, Dv)
    for m, n, r, s in itertools.product(Rg, Rg, Rg, Rg):
        Xv = Xc(m, n, r, s, LamV, bV, cV)
        if Xv == 0: continue
        for M, N in itertools.product(Rg, Rg):
            # U1..U4 (derived Table II, U2=+i/4) contracted with X, coefficient of d_{1M} d_{2N}
            K[M, N] += kap**2*a**(Dv - 2)*Xv*(-I/8*ei[m, n]*ei[r, s]*ei[M, N] + I/4*ei[m, r]*ei[n, s]*ei[M, N]
                                             + I/2*ei[m, n]*ei[r, M]*ei[s, N] - I*ei[m, M]*ei[n, r]*ei[s, N])
    return K
dd = lambda g, m: sp.diff(g, X[m])
kern = lambda K: -sum(dd((K[m, n] + K[n, m])/2*dd(f, n), m) for m in Rg for n in Rg)
up_div = lambda g: sum(ei[m, m]*dd(g*dd(f, m), m) for m in Rg)          # d^mu(g d_mu f)
lap = lambda g: sum(dd(dd(g, j), j) for j in (1, 2, 3))
lhs = sp.expand(kern(Kt(Lam, bb, -bb)))                                  # c = -b at D=4
rhs = sp.expand(-I*kap**2*(3*up_div(a**2*Lam) - 4*a**2*Lam*lap(f)))
print("  kernel(X with Lam, b, c=-b) == -i kappa^2 {3 d^mu(a^2 Lam d_mu) - 4 a^2 Lam Lap} f :", sp.simplify(lhs - rhs) == 0)
print("  B/C finite part alone (Lam=0, c=-b) gives zero kernel:", sp.simplify(kern(Kt(0, bb, -bb))) == 0)

print("== (C) Eq. (43) as an identity of kernels acting on test functions ==")
tt = sp.Symbol('t', negative=True)
af = sp.Function('a')(tt)
L = sp.log(af)
g = sp.Function('g')(tt, x1)            # test function (one spatial direction suffices: a depends on t only)
def ev(expr):                            # impose a' = H a^2 (from a = -1/(H t))
    expr = expr.doit()
    expr = expr.subs(sp.Derivative(af, (tt, 2)), 2*H**2*af**3).subs(sp.Derivative(af, tt), H*af**2)
    return sp.simplify(sp.expand(expr))
divD = lambda w, G: -sp.diff(af**w*sp.diff(G, tt), tt) + af**w*sp.diff(G, x1, 2)   # d^mu(a^w d_mu G)
for Dval, lab in ((4, 'D=4'), (D, 'general D')):
    lhs43 = -sp.diff(af**(Dval - 2)*L*sp.diff(g, tt), tt) + af**(Dval - 2)*L*sp.diff(g, x1, 2)   # d^mu(a^{D-2} ln a d_mu g)
    # (1/2) ln(a a') a^D Box delta acting on g:  (1/2)[ L * a^D Box g + a^D Box (L g) ],  a^D Box G = d^mu(a^{D-2} d_mu G)
    rhs43 = R(1, 2)*(L*divD(Dval - 2, g) + divD(Dval - 2, L*g)) + R(1, 2)*(Dval - 1)*H**2*af**Dval*g
    wrong = L*divD(Dval - 2, g) + R(1, 2)*(Dval - 1)*H**2*af**Dval*g      # with ln(a^2) instead of ln(aa')
    print(f"  {lab}: d^mu(a^(D-2) ln a d_mu delta) == (1/2)ln(aa') a^D Box delta + (D-1)/2 H^2 a^D delta :",
          ev(lhs43 - rhs43) == 0, "  | with ln(a^2) instead:", ev(lhs43 - wrong) == 0)
print("  Box ln a = -(D-1) H^2 :", ev(divD(D - 2, L)/af**D + (D - 1)*H**2) == 0)

print("== (D) (41) -> (44) on a test function, D=4 (pihat, k kept symbolic) ==")
L4 = sp.log(a)
def up_div_g(G): return sum(ei[m, m]*dd(a**2*dd(G, m), m) for m in Rg)  # a^4 Box G
lapf = lap(f)
def M41(sign, extra):    # sign*i kappa^2 k { pihat[3 a^4Box - 4 a^2 Lap] + 2[3 d(a^2 ln a d) - 4 ln a a^2 Lap] + extra*2[-a^4Box + a^2 Lap] } f
    return sign*I*kap**2*k*(pih*(3*up_div_g(f) - 4*a**2*lapf) + 2*(3*up_div(a**2*L4) - 4*L4*a**2*lapf)
                            + extra*2*(-up_div_g(f) + a**2*lapf))
def M44(sign, c0):       # sign*i kappa^2 k { [3 pihat + c0*(-2) + 3 ln(aa')] a^4Box + [-4 pihat + c0*2 - 4 ln(aa')] a^2 Lap + 9 H^2 a^4 } f
    lnaa_box = L4*up_div_g(f) + up_div_g(L4*f)          # ln(aa') a^4 Box delta acting on f
    lnaa_lap = L4*a**2*lapf + a**2*lap(L4*f)             # ln(aa') a^2 Lap delta acting on f
    return sign*I*kap**2*k*((3*pih - 2*c0)*up_div_g(f) + 3*lnaa_box + (-4*pih + 2*c0)*a**2*lapf - 4*lnaa_lap + 9*H**2*a**4*f)
print("  derived (41) == derived (44)  [sign -, no -2/+2]:", sp.simplify(sp.expand(M41(-1, 0) - M44(-1, 0))) == 0)
print("  paper   (41) == paper   (44)  [sign +, with -2/+2]:", sp.simplify(sp.expand(M41(+1, 1) - M44(+1, 1))) == 0)
print("  (k -> H^2/8pi^2 at D=4 makes 2k = H^2/4pi^2, the prefactor of the finite bracket in (41))")
