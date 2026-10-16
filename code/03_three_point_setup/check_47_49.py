"""Task 7 -- checks of Eqs. (I.47)-(I.49): the graviton propagator inserted into the 3-point contribution.

Paper I = E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

Purpose
    Supplementary code checks of the analytic derivation in Chapter 3, Section "Task 7" of the notes.
    They do not replace the derivation.

Method
    For D = 4, 5, 6 the tensor factors T^A, T^B, T^C of (I.22)-(I.24) are built as explicit component arrays.
    Outer derivatives (acting on everything) and inner derivatives (acting only on i Delta_A) are represented
    by formal vectors P, P' and Q, Q'; the scalar propagator multiplying the outer bracket is tracked by the
    symbols dA, dB, dC.  Two routes are compared:
      route 1: the derived (I.46) (four terms with traces),
      route 2: directly (I.45) with Table I and the rules d_1 -> -P, d_2 -> +Q, d'_1 -> -P', d'_2 -> +Q'.
    The result is compared with the derived nine terms (our (I.47')-(I.49')) plus the antisymmetric AA pair,
    and with the nine terms as printed in Paper I.  Homogeneity of i Delta_A (d'_i = -d_i on i Delta_A) is
    imposed where stated, and verified separately for F = A(y) + k ln(a a').

Expected result
    For each D: route1 == route2 True; route1 == derived + antisym. pair True; antisym. pair vanishes True;
    route1 == printed (I.47)-(I.49) False; route1 == -(printed (I.47)-(I.49)) True.  Homogeneity: True, True.

Runtime
    About 2 seconds.

Origin
    Run as check_47_49.py in Task 7 (15 Oct 2026).  The original one-line header comment was replaced by this
    docstring; the code is unchanged.
"""
import sympy as sp, itertools
I, R = sp.I, sp.Rational
dA, dB, dC = sp.symbols('dA dB dC')        # which scalar propagator sits in the outer bracket: (aa')^{D-2} i Delta_I

def setup(Dv):
    d = sp.Integer(Dv); Rg = range(Dv)
    eta = sp.diag(*([-1] + [1]*(Dv - 1))); ei = eta.inv()
    u = [1 if m == 0 else 0 for m in Rg]
    eb = [[eta[m, n] + u[m]*u[n] for n in Rg] for m in Rg]
    TA = lambda m, n, r, s: eb[m][r]*eb[s][n] + eb[m][s]*eb[r][n] - 2/(d - 3)*eb[m][n]*eb[r][s]
    TB = lambda m, n, r, s: -(u[m]*eb[n][r]*u[s] + u[n]*eb[m][r]*u[s] + u[m]*eb[n][s]*u[r] + u[n]*eb[m][s]*u[r])
    fT = lambda m, n: (d - 3)*u[m]*u[n] + eb[m][n]
    TC = lambda m, n, r, s: 2/((d - 2)*(d - 3))*fT(m, n)*fT(r, s)
    X = lambda m, n, r, s: TA(m, n, r, s)*dA + TB(m, n, r, s)*dB + TC(m, n, r, s)*dC
    P = sp.symbols(f'P0:{Dv}'); Q = sp.symbols(f'Q0:{Dv}'); Pp = sp.symbols(f'Pp0:{Dv}'); Qp = sp.symbols(f'Qp0:{Dv}')
    return d, Rg, eta, ei, X, P, Q, Pp, Qp

for Dv in (4, 5, 6):
    d, Rg, eta, ei, X, P, Q, Pp, Qp = setup(Dv)
    up = lambda v: [sum(ei[m, n]*v[n] for n in Rg) for m in Rg]
    dot = lambda a_, b_: sum(ei[m, n]*a_[m]*b_[n] for m in Rg for n in Rg)
    sdot = lambda a_, b_: sum(a_[i]*b_[i] for i in range(1, Dv))          # spatial dot
    Pu, Qu, Ppu, Qpu = up(P), up(Q), up(Pp), up(Qp)
    quad = list(itertools.product(Rg, Rg, Rg, Rg))
    # ---- route 1: from (46): outer P,P' outside, inner Q,Q' inside (commuting bookkeeping) ----
    t1 = -sum(X(a, b, r, s)*Pu[a]*Qu[b]*Ppu[r]*Qpu[s] for a, b, r, s in quad)
    t2 = R(1, 2)*dot(P, Q)*sum(ei[a, b]*X(a, b, r, s)*Ppu[r]*Qpu[s] for a, b, r, s in quad)
    t3 = R(1, 2)*dot(Pp, Qp)*sum(X(a, b, r, s)*ei[r, s]*Pu[a]*Qu[b] for a, b, r, s in quad)
    t4 = -R(1, 4)*dot(P, Q)*dot(Pp, Qp)*sum(ei[a, b]*ei[r, s]*X(a, b, r, s) for a, b, r, s in quad)
    from46 = sp.expand(t1 + t2 + t3 + t4)
    # ---- route 2: directly from (45) with Table I and rules d1->-P, d2->Q, d1'->-P', d2'->Q' ----
    V = lambda d1u, d2u, d1, d2: (lambda a, b: I*d1u[a]*d2u[b] - I/2*ei[a, b]*dot(d1, d2))
    Vx = V([-v for v in Pu], Qu, [-v for v in P], Q); Vxp = V([-v for v in Ppu], Qpu, [-v for v in Pp], Qp)
    from45 = sp.expand(sum(Vx(a, b)*Vxp(r, s)*X(a, b, r, s) for a, b, r, s in quad))
    # ---- derived (47)-(49), kappa^2 (aa')^{D-2} stripped ----
    P0, Q0, Pp0, Qp0 = P[0], Q[0], Pp[0], Qp[0]
    AA = (-sdot(P, Pp)*sdot(Q, Qp) + (d - 1)/(d - 3)*P0*Pp0*Q0*Qp0 - sdot(P, Q)*Pp0*Qp0 - P0*Q0*sdot(Pp, Qp))*dA
    AB = (P0*Pp0*sdot(Q, Qp) + sdot(P, Qp)*Pp0*Q0 + P0*sdot(Pp, Q)*Qp0 + sdot(P, Pp)*Q0*Qp0)*dB
    AC = -2*(d - 2)/(d - 3)*P0*Pp0*Q0*Qp0*dC
    antisym = (sdot(P, Q)*sdot(Pp, Qp) - sdot(P, Qp)*sdot(Pp, Q))*dA
    derived = sp.expand(AA + AB + AC)
    paper = sp.expand(-(AA + AB + AC))                    # paper's (47)-(49) as printed
    hom = {Qp[i]: -Q[i] for i in range(1, Dv)}           # inner derivatives on i Delta_A(x;x'): d'_i = -d_i
    print(f"D={Dv}: (46)-route == (45)-route                         :", sp.expand(from46 - from45) == 0)
    print(f"       (46)-route == derived (47)-(49) + antisym. pair      :", sp.expand(from46 - derived - antisym) == 0)
    print(f"       antisym. pair vanishes when d'_i = -d_i on iDelta_A  :", sp.expand(antisym.subs(hom)) == 0)
    print(f"       (46)-route == paper's (47)-(49)  (after d'_i=-d_i)   :", sp.expand((from46 - paper).subs(hom)) == 0)
    print(f"       (46)-route == -(paper's (47)-(49)) (after d'_i=-d_i) :", sp.expand((from46 + paper).subs(hom)) == 0)

print("== homogeneity of iDelta_A: d'_i F = -d_i F and d_i d'_j F symmetric, F = A(y) + k ln(a a') ==")
H, k = sp.symbols('H k', positive=True)
t, tp = sp.symbols('t tp', negative=True)
xs = sp.symbols('x1:4', real=True); xps = sp.symbols('y1:4', real=True)
a = -1/(H*t); ap = -1/(H*tp)
y = a*ap*H**2*(sum((xs[i] - xps[i])**2 for i in range(3)) - (t - tp)**2)
Af = sp.Function('A')
F = Af(y) + k*sp.log(a*ap)
print("  d'_i F + d_i F = 0           :", all(sp.simplify(sp.diff(F, xps[i]) + sp.diff(F, xs[i])) == 0 for i in range(3)))
print("  d_i d'_j F - d_j d'_i F = 0  :", all(sp.simplify(sp.diff(F, xs[i], xps[j]) - sp.diff(F, xs[j], xps[i])) == 0
                                           for i in range(3) for j in range(3)))
