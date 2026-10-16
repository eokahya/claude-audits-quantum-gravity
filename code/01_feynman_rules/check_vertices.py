"""check_vertices.py -- supplementary computer-algebra checks of Eqs. (I.7), (I.8),
Table I and Table II (Task 2).

Purpose
    Verify the 3-point vertex operators V_I^{ab} (Table I) and the 4-point vertex
    operators U_I^{abrs} (Table II) of E. O. Kahya and R. P. Woodard,
    Phys. Rev. D 76, 124005 (2007), obtained from Eq. (I.6) with the paper's
    normalization:
        kappa   L^(1) = (1/2)(1/i) sum_I [V_I^{ab}   phi_1 phi_2] h_ab,
        kappa^2 L^(2) = (1/2)(1/i) sum_I [U_I^{abrs} phi_1 phi_2] h_ab h_rs
    (a factor 1/2 extracted for the two identical scalars, none for the two
    gravitons, a factor i per vertex from exp(iS)).

Conventions
    p_mu stands for d_mu acting on a scalar leg (lower index); P^a = eta^{ab} p_b.
    Overall factors a^{D-2} are set to 1; powers of kappa are tracked by hand
    (V ~ kappa, U ~ kappa^2).  The sign of U_2 is a switch: s2 = +1 is the
    derived value +i/4, s2 = -1 the printed value -i/4.

Method
    (A) Full vertices from Eq. (I.6) by polarization with distinct legs
        (p1 != p2, h3 != h4), D = 3, 4 fully symbolic:
        i * coeff[phi1 phi2 h3] of kappa L^(1) versus V_1 + V_2, and
        i * coeff[phi1 phi2 h3 h4] of kappa^2 L^(2) versus U[3,4] + U[4,3],
        for the derived and for the printed Table II.
    (B) Tables versus the EXACT metric (no Eq. (I.6) input), coincident legs:
        D = 2, 3 symbolic h; D = 4, 5 random rational h (random.seed(7)),
        symbolic p.  Checks kappa L^(1) = V[h]/(2i), kappa^2 L^(2) = U[h,h]/(2i).
    (C) Determines the table coefficients from the exact metric alone: general
        ansatz V = i[x d1^a d2^b + y eta^{ab} d1.d2],
        U = i[u1 eta^ab eta^rs d1.d2 + u2 eta^ar eta^bs d1.d2 + u3 eta^ab d1^r d2^s
              + u4 d1^a eta^br d2^s],
        fitted to (i) the conformal case h = c eta for general D and
        (ii) the case h_00 = b only, D = 4, generic p; solves the linear system.

Expected result
    (A) True / 0 for the derived table; for the printed table a residual
        -I*(p.q)*(h3^{rs} h4_{rs}) (the sign of U_2), printed in components as
        I*(p0*q0 - p1*q1 - ...)*(a00*b00 - 2*a01*b01 + ...), where
        p0*q0 - p1*q1 - ... = -p.q.
    (B) all True; printed table minus exact = -(1/4)(d phi.d phi)(h^{rs}h_{rs})
        (shown in components; numbers for the random rational h).
    (C) unique solution x=1, y=-1/2, u=(-1/8, 1/4, 1/2, -1).

Runtime
    About 10 s (sympy 1.14, Apple silicon).

Run
    python3 check_vertices.py
"""
import random

import sympy as sp

I = sp.I
R = sp.Rational


def eta(D): return sp.diag(*([-1] + [1]*(D-1)))


def symm(D, name):
    M = sp.zeros(D, D)
    for i in range(D):
        for j in range(i, D):
            s = sp.Symbol(f'{name}{i}{j}'); M[i, j] = s; M[j, i] = s
    return M


def vec(D, name): return sp.Matrix([sp.Symbol(f'{name}{i}') for i in range(D)])


def ip(u, M, v): return (u.T*M*v)[0, 0]


# ---- Lagrangian pieces read off from Eq. (6)  (a = 1, explicit kappa removed) ----
def L1_eq6(p, h, D):
    E = eta(D); Ei = E.inv(); hu = Ei*h*Ei; tr = (Ei*h).trace()
    return -R(1, 2)*ip(p, -hu + R(1, 2)*Ei*tr, p)


def L2_eq6(p, h, D):
    E = eta(D); Ei = E.inv(); hu = Ei*h*Ei; tr = (Ei*h).trace()
    hh = sum(h[i, j]*hu[i, j] for i in range(D) for j in range(D))
    chain = hu*E*hu
    return -R(1, 2)*ip(p, chain - R(1, 2)*tr*hu + R(1, 8)*Ei*tr**2 - R(1, 4)*Ei*hh, p)


# ---- vertex tables, contracted: V(p1,p2)[h], U(p1,p2)[h3,h4]  (kappa, a stripped) ----
def V_table(p1, p2, h, D):
    Ei = eta(D).inv(); P1, P2 = Ei*p1, Ei*p2; tr = (Ei*h).trace()
    V1 = I*ip(P1, h, P2)                       # i  d1^a d2^b
    V2 = -I/2*tr*ip(p1, Ei, p2)                # -(i/2) eta^{ab} d1.d2
    return V1 + V2


def U_table(p1, p2, h3, h4, D, s2):
    Ei = eta(D).inv(); P1, P2 = Ei*p1, Ei*p2
    t3, t4 = (Ei*h3).trace(), (Ei*h4).trace(); d12 = ip(p1, Ei, p2)
    h3u = Ei*h3*Ei; h3h4 = sum(h3u[i, j]*h4[i, j] for i in range(D) for j in range(D))
    U1 = -I/8*t3*t4*d12                        # -(i/8) eta^{ab} eta^{rs} d1.d2
    U2 = s2*I/4*h3h4*d12                       # s2*(i/4) eta^{ar} eta^{bs} d1.d2  (s2=+1 derived, -1 paper)
    U3 = I/2*t3*ip(P1, h4, P2)                 # (i/2) eta^{ab} d1^r d2^s
    U4 = -I*ip(P1, h3*Ei*h4, P2)               # -i d1^a eta^{br} d2^s
    return U1 + U2 + U3 + U4


t1, t2, s3, s4, k = sp.symbols('t1 t2 s3 s4 kappa')

print("== (A) full vertices from Eq. (6) by polarization, distinct legs p1,p2,h3,h4 ==")
for D in (3, 4):
    p1, p2, h3, h4 = vec(D, 'p'), vec(D, 'q'), symm(D, 'a'), symm(D, 'b')
    L1 = sp.expand(L1_eq6(t1*p1 + t2*p2, h3, D))
    full3 = sp.expand(I*L1.coeff(t1, 1).coeff(t2, 1))
    print(f"D={D}: i*coeff[phi1 phi2 h3] of kappa L^(1)  == Table I :", sp.expand(full3 - V_table(p1, p2, h3, D)) == 0)
    L2 = sp.expand(L2_eq6(t1*p1 + t2*p2, s3*h3 + s4*h4, D))
    full4 = sp.expand(I*L2.coeff(t1, 1).coeff(t2, 1).coeff(s3, 1).coeff(s4, 1))
    for s2, lab in ((+1, 'derived (U2=+i/4)'), (-1, 'paper   (U2=-i/4)')):
        rhs = U_table(p1, p2, h3, h4, D, s2) + U_table(p1, p2, h4, h3, D, s2)
        diff = sp.factor(sp.expand(rhs - full4))
        print(f"D={D}: U[3,4]+U[4,3] - i*coeff[phi1 phi2 h3 h4], {lab}:", diff)

print("== (B) tables vs EXACT metric (no Eq. (6) input), coincident legs ==")


def exact_series(h, D):
    g = eta(D) + k*h; detg = g.det()
    M = sp.sqrt(-detg)*g.adjugate()/detg
    return M.applyfunc(lambda e: sp.series(e, k, 0, 3).removeO())


random.seed(7)
for D in (2, 3, 4, 5):
    if D <= 3: h = symm(D, 'h')
    else:
        h = sp.zeros(D, D)
        for i in range(D):
            for j in range(i, D):
                v = R(random.randint(-9, 9), random.randint(1, 7)); h[i, j] = v; h[j, i] = v
    p = vec(D, 'p')
    Lex = sp.expand(-R(1, 2)*ip(p, exact_series(h, D), p))
    ex1, ex2 = Lex.coeff(k, 1), Lex.coeff(k, 2)
    ok3 = sp.expand(V_table(p, p, h, D)/(2*I) - ex1) == 0
    ok4d = sp.expand(U_table(p, p, h, h, D, +1)/(2*I) - ex2) == 0
    d4p = sp.factor(sp.expand(U_table(p, p, h, h, D, -1)/(2*I) - ex2))
    tag = 'symbolic h' if D <= 3 else 'random rational h'
    print(f"D={D} ({tag}): L^(1) == (1/2i)V[h]: {ok3};  L^(2) == (1/2i)U[h,h] derived: {ok4d};  paper - exact = {d4p}")

print("== (C) determine the table coefficients from the exact metric alone ==")
Dg, c, b = sp.symbols('D c b')
x, y, u1, u2, u3, u4 = sp.symbols('x y u1 u2 u3 u4')
# conformal case h = c*eta, general D: structures evaluated with eta^{ab}eta_{ab} = D
S = sp.Symbol('S')    # S = d phi . d phi
ex = sp.series(-R(1, 2)*S*(1 + k*c)**((Dg - 2)/2), k, 0, 3).removeO()
rec3 = R(1, 2)*(x*c + y*c*Dg)*S                       # (1/2i)*i*(x P.h.P + y tr(h) S)
rec4 = R(1, 2)*(u1*Dg**2 + u2*Dg + u3*Dg + u4)*c**2*S
eqs = sp.Poly(sp.expand(rec3 - ex.coeff(k, 1)), Dg).coeffs() + sp.Poly(sp.expand(rec4 - ex.coeff(k, 2)), Dg).coeffs()
# h_00 = b case, explicit D = 4, generic p
D = 4; E = eta(D); Ei = E.inv(); p = vec(D, 'p'); P = Ei*p
h = sp.zeros(D, D); h[0, 0] = b
Lex = sp.expand(-R(1, 2)*ip(p, exact_series(h, D), p))
hu = Ei*h*Ei; tr = (Ei*h).trace(); hh = sum(h[i, j]*hu[i, j] for i in range(D) for j in range(D))
rec4b = R(1, 2)*(u1*tr**2*ip(p, Ei, p) + u2*hh*ip(p, Ei, p) + u3*tr*ip(P, h, P) + u4*ip(P, h*Ei*h, P))
eqs += sp.Poly(sp.expand(rec4b - Lex.coeff(k, 2)), *list(p)).coeffs()
print("solution:", sp.solve(eqs, [x, y, u1, u2, u3, u4], dict=True))
print("   (Table I : x=1 -> V1 = i kappa d1^a d2^b ; y=-1/2 -> V2 = -(i/2) kappa eta^{ab} d1.d2)")
print("   (Table II: U_I = i kappa^2 a^{D-2} u_I x structure_I)")
