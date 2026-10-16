#!/usr/bin/env python3
"""
check_54_59.py -- code checks for the local 3-point contributions, Eqs. (I.54)-(I.59)
of Kahya & Woodard, Phys. Rev. D 76, 124005 (2007) ("Paper I").

Purpose
    Supplementary (not a substitute for the analytic derivation in the notes,
    chapter "3-point local and logarithm terms: Eqs. (I.54)-(I.67)", Task 9).

What it checks
    (A) The local part of the 3-point contribution extracted DIRECTLY from (I.45)
        and Table I (the graviton tensor factors T^A, T^B, T^C built as explicit
        component arrays for D = 4, 5, 6), without using the intermediate step
        (I.47)-(I.49).  The inner derivative pair is replaced by the delta-function
        part of (I.52), i.e. Q_rho Q'_sigma -> delta^0_rho delta^0_sigma.
        Expected: local polynomial = (D-1)/(D-3) P0 P0' dA + (P.P')_spatial dB
                  - 2(D-2)/(D-3) P0 P0' dC   (dA, dB, dC label the outer propagator).
    (B) The outer-derivative rules (R2) of the notes (Task 9 section) and Paper I's
        (I.55)-(I.56), acting on a test function g, with explicit
        a(eta) = -1/(H eta) and symbolic D in the powers of a:
            d0 d0'[a^{D-2} delta]       = (a^D Box - a^{D-2} Lap) delta
            sum_i d_i d_i'[a^{D-2} delta] = - a^{D-2} Lap delta
    (C) Assembly of (I.59) in units of i kappa^2 k,
            k = H^{D-2}/(4 pi)^{D/2} Gamma(D-1)/Gamma(D/2),
        using the coincidence limits A -> k*pihat, b = -k/(D-2), c = k/((D-2)(D-3)).
        Prints the derived coefficients (definition sign), their finite parts at D = 4
        (Box: -2, Lap: +5/2), the sum of the paper's own (I.54)+(I.57)+(I.58)
        (Lap finite part -5/2, whereas (I.59) prints -7/2), and checks that the
        paper's (I.54), (I.57), (I.58) are -1 times the derived ones term by term.

Method
    SymPy exact algebra; the delta function is represented by its action on a
    test function g(t, x1, x2, x3).

Expected result
    Every True/False check prints True; finite parts Box: -2, Lap: 5/2.

Runtime
    About 1 s (SymPy 1.14, Python 3.10+).

Run
    cd code/04_three_point_local_log && python3 check_54_59.py
"""
# Supplementary checks for Eqs. (54)-(59): local 3-point contributions.
import sympy as sp, itertools
I, R = sp.I, sp.Rational
dA, dB, dC = sp.symbols('dA dB dC')
Dsym = sp.Symbol('D')

print("== (A) local part extracted directly from (45) + Table I (no use of (47)-(49)) ==")
for Dv in (4, 5, 6):
    d = sp.Integer(Dv); Rg = range(Dv)
    eta = sp.diag(*([-1] + [1]*(Dv - 1))); ei = eta.inv()
    u = [1 if m == 0 else 0 for m in Rg]
    eb = [[eta[m, n] + u[m]*u[n] for n in Rg] for m in Rg]
    TA = lambda m, n, r, s: eb[m][r]*eb[s][n] + eb[m][s]*eb[r][n] - 2/(d - 3)*eb[m][n]*eb[r][s]
    TB = lambda m, n, r, s: -(u[m]*eb[n][r]*u[s] + u[n]*eb[m][r]*u[s] + u[m]*eb[n][s]*u[r] + u[n]*eb[m][s]*u[r])
    fT = lambda m, n: (d - 3)*u[m]*u[n] + eb[m][n]
    TC = lambda m, n, r, s: 2/((d - 2)*(d - 3))*fT(m, n)*fT(r, s)
    X = lambda m, n, r, s: TA(m, n, r, s)*dA + TB(m, n, r, s)*dB + TC(m, n, r, s)*dC
    P = sp.symbols(f'P0:{Dv}'); Pp = sp.symbols(f'Pp0:{Dv}')
    # inner derivatives replaced by the delta-part of (52): Q_rho Q'_sigma -> delta^0_rho delta^0_sigma
    Q = [1] + [0]*(Dv - 1); Qp = [1] + [0]*(Dv - 1)
    up = lambda v: [sum(ei[m, n]*v[n] for n in Rg) for m in Rg]
    dot = lambda a_, b_: sum(ei[m, n]*a_[m]*b_[n] for m in Rg for n in Rg)
    Pu, Qu, Ppu, Qpu = up(P), up(Q), up(Pp), up(Qp)
    V = lambda d1u, d2u, d1, d2: (lambda a, b: I*d1u[a]*d2u[b] - I/2*ei[a, b]*dot(d1, d2))
    Vx = V([-v for v in Pu], Qu, [-v for v in P], Q); Vxp = V([-v for v in Ppu], Qpu, [-v for v in Pp], Qp)
    loc = sp.expand(sum(Vx(a, b)*Vxp(r, s)*X(a, b, r, s) for a, b, r, s in itertools.product(Rg, Rg, Rg, Rg)))
    sdot = sum(P[i]*Pp[i] for i in range(1, Dv))
    expected = sp.expand((d - 1)/(d - 3)*P[0]*Pp[0]*dA + sdot*dB - 2*(d - 2)/(d - 3)*P[0]*Pp[0]*dC)
    print(f"  D={Dv}: local polynomial == (D-1)/(D-3) P0P0' dA + (P.P')_spatial dB - 2(D-2)/(D-3) P0P0' dC :",
          sp.expand(loc - expected) == 0)

print("== (B) outer derivatives on a^{D-2} delta^D(x-x'), acting on a test function g (D=4 coordinates) ==")
H = sp.Symbol('H', positive=True); t = sp.Symbol('t', negative=True)
x1, x2, x3 = xs = sp.symbols('x1 x2 x3', real=True)
a = -1/(H*t); Dn = sp.Symbol('D')
g = sp.Function('g')(t, *xs)
lap = lambda f: sum(sp.diff(f, xi, 2) for xi in xs)
aD2 = a**(Dn - 2)
# int dx' d_mu d'_nu [a^{D-2} delta(x-x')] g(x') = d_mu [a^{D-2} (-d_nu g)]
op00 = sp.diff(aD2*(-sp.diff(g, t)), t)                         # d0 d0'
opii = sum(sp.diff(aD2*(-sp.diff(g, xi)), xi) for xi in xs)    # sum_i d_i d_i'
aDbox = -sp.diff(aD2*sp.diff(g, t), t) + aD2*lap(g)             # a^D Box g = d^mu(a^{D-2} d_mu g)
print("  d0 d0'[a^{D-2} delta] == (a^D Box - a^{D-2} Lap) delta :", sp.simplify(op00 - (aDbox - aD2*lap(g))) == 0)
print("  sum_i d_i d_i'[a^{D-2} delta] == - a^{D-2} Lap delta   :", sp.simplify(opii + aD2*lap(g)) == 0)
print("  (55)-(56): -d0(a^{D-2} d0' delta) == -a^D Box + a^{D-2} Lap :", sp.simplify(-op00 - (-aDbox + aD2*lap(g))) == 0)

print("== (C) assembling (59): units of i kappa^2 k, k = H^{D-2}/(4pi)^{D/2} Gamma(D-1)/Gamma(D/2) ==")
pihat = sp.Symbol('pihat')          # pihat = -pi cot(pi D/2) = -2/(D-4) + O(D-4)
D = Dsym
LamA, b, c = pihat, -1/(D - 2), 1/((D - 2)*(D - 3))         # coincidence limits / k (Task 4), A without the log
# local polynomial with P0P0' -> (Box_ - Lap_), (P.P')_s -> -Lap_
BOX, LAP = sp.symbols('BOX LAP')     # BOX = a^D Box delta, LAP = a^{D-2} Lap delta
derived = sp.expand(I*((D - 1)/(D - 3)*LamA*(BOX - LAP) + b*(-LAP) - 2*(D - 2)/(D - 3)*c*(BOX - LAP)))
cB, cL = derived.coeff(BOX), derived.coeff(LAP)
print("  derived coeff of a^D Box     :", sp.simplify(cB/I))
print("  derived coeff of a^{D-2} Lap :", sp.simplify(cL/I))
fin = lambda e: sp.simplify((e/I).subs(pihat, 0).subs(D, 4))
print("  finite parts at D=4 (pihat-terms kept exact):  Box:", fin(cB), "   Lap:", fin(cL))
# paper's own (54),(57),(58) in the same units, D=4 for finite parts: (57) -> -1/2 Lap ; (58) -> +2 (Box - Lap)
p54 = I*(D - 1)/(D - 3)*pihat*(-BOX + LAP)
paper_sum = sp.expand(p54 + I*(-R(1, 2))*LAP + I*2*(BOX - LAP))
print("  sum of the paper's own (54)+(57)+(58):  Box:", sp.simplify(paper_sum.coeff(BOX)/I), "   Lap:", sp.simplify(paper_sum.coeff(LAP)/I))
print("  paper's printed (59):                    Box: -(D-1)/(D-3)*pihat + 2   Lap: (D-1)/(D-3)*pihat - 7/2")
print("  paper's (54),(57),(58) = -(derived) term by term:",
      sp.simplify(p54 + I*(D - 1)/(D - 3)*LamA*(BOX - LAP)) == 0,
      sp.simplify((-R(1, 2)) - (-(b*(-1))).subs(D, 4)) == 0,
      sp.simplify(2 - (2*(D - 2)/(D - 3)*c).subs(D, 4)) == 0)
print("  check k -> H^2/(8 pi^2) at D=4:", sp.simplify((H**(D - 2)/(4*sp.pi)**(D/2)*sp.gamma(D - 1)/sp.gamma(D/2)).subs(D, 4)))
