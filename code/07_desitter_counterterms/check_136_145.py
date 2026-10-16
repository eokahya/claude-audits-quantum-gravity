"""
check_136_145.py -- code checks for Task 15 (Paper I, Eqs. (I.136)-(I.145)).

Paper I: E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

Purpose
    Independent symbolic (SymPy) checks of the reduction of the de Sitter
    noninvariant counterterms in Task 15.  Supplements, and does not replace,
    the analytic derivation in the notes.

Notation
    a = -1/(H eta), t stands for conformal time eta < 0, D = 4 spacetime
    coordinates (eta, x1, x2, x3) for the field, but the power D of a is kept
    symbolic wherever possible.  Basis of admissible counterterm densities:
      B1 = a^D (Box phi)^2                  (de Sitter invariant)
      B2 = H^2 a^{D-2} eta^mn d_m phi d_n phi (de Sitter invariant)
      B3 = a^{D-2} Box phi Lap phi
      B4 = a^{D-4} (Lap phi)^2
      B5 = H^2 a^{D-2} grad phi . grad phi
    with Box = a^{-D} d_m (a^{D-2} eta^mn d_n)  [Paper I, Eq. (I.42)].

What is checked
    * (I.136), (I.137): exact identities, symbolic D.
    * (I.138), (I.140), (I.141), (I.143): the difference across each "->" step
      equals exactly the divergence of the stated discarded current; (I.139)
      and (I.142) are exact equalities.  Symbolic D.
    * Completeness (D = 4, 5, 6 only): the Euler-Lagrange operators of B1..B5
      are linearly independent, and the Euler-Lagrange operator of each of
      15 admissible candidate monomials lies in span{B1..B5}; the coefficient
      vectors c(B1..B5) are printed.  Equality of Euler-Lagrange operators is
      equality of the densities up to total derivatives.
    * (I.145): the Euler-Lagrange operator of -1/2 alpha3 kappa^2 H^2 a^{D-2}
      (grad phi)^2 is alpha3 kappa^2 H^2 a^{D-2} Lap phi, and the two printed
      forms of the counterterm are identical.
    * Remark after (I.145): the Euler-Lagrange operator of
      c(eta) d0phi Lap phi is -c'(eta) Lap phi (no d0 survives).

Expected result
    All tests print True; the completeness table is reproduced exactly as in
    expected_output/check_136_145.txt.

Runtime
    About 5 s (Python 3.13, sympy 1.14).

Run
    python3 check_136_145.py
"""
import sympy as sp
from sympy.calculus.euler import euler_equations
R = sp.Rational
H = sp.Symbol('H', positive=True); D = sp.Symbol('D')
t = sp.Symbol('t', negative=True); x1, x2, x3 = xs = sp.symbols('x1 x2 x3', real=True)
X = (t, x1, x2, x3)
a = -1/(H*t)
phi = sp.Function('phi')(*X)
d0 = lambda f: sp.diff(f, t)
lap = lambda f: sum(sp.diff(f, xi, 2) for xi in xs)
gradsq = sum(sp.diff(phi, xi)**2 for xi in xs)
box = lambda f: (-sp.diff(a**(D - 2)*d0(f), t) + a**(D - 2)*lap(f))/a**D          # (42), D symbolic
etadd = -d0(phi)**2 + gradsq                                                         # eta^{mn} d_m phi d_n phi
div = lambda J: sum(sp.diff(J[m], X[m]) for m in range(4))
Z = lambda e: sp.simplify(sp.expand(e)) == 0

print("== (136)-(137) ==")
print("  (136) a^-2 d0^2 phi = -Box phi - (D-2)(H/a) d0 phi + a^-2 Lap phi :", Z(d0(d0(phi))/a**2 - (-box(phi) - (D - 2)*H/a*d0(phi) + lap(phi)/a**2)))
print("  (137) (d0 phi)^2/a^2 = -g^mn d phi d phi + a^-2 grad phi^2         :", Z(d0(phi)**2/a**2 - (-etadd/a**2 + gradsq/a**2)))

print("== (138)-(140): H a^{D-1} d0phi Box phi ==")
T0 = H*a**(D - 1)*d0(phi)*box(phi)
T138 = -H*a**(D - 3)*sum((-1 if m == 0 else 1)*sp.diff(phi, X[m], t)*sp.diff(phi, X[m]) for m in range(4)) - H**2*a**(D - 2)*d0(phi)**2
J138 = [H*a**(D - 3)*d0(phi)*(-1 if m == 0 else 1)*sp.diff(phi, X[m]) for m in range(4)]
print("  (138) T0 - [rhs] = d_m (H a^{D-3} d0phi d^m phi)                :", Z(T0 - T138 - div(J138)))
T139 = -R(1, 2)*H*a**(D - 3)*d0(etadd) + H**2*a**(D - 2)*etadd - H**2*a**(D - 2)*gradsq
print("  (139) exact equality                                            :", Z(T138 - T139))
T140 = R(1, 2)*(D - 1)*H**2*a**(D - 2)*etadd - H**2*a**(D - 2)*gradsq
print("  (140) (139) - [rhs] = d0(-1/2 H a^{D-3} eta^mn dphi dphi)        :", Z(T139 - T140 - d0(-R(1, 2)*H*a**(D - 3)*etadd)))

print("== (141)-(143): H a^{D-3} d0phi Lap phi ==")
U0 = H*a**(D - 3)*d0(phi)*lap(phi)
U141 = -H*a**(D - 3)*sum(sp.diff(phi, xi, t)*sp.diff(phi, xi) for xi in xs)
Jsp = [0] + [H*a**(D - 3)*d0(phi)*sp.diff(phi, xi) for xi in xs]
print("  (141) U0 - [rhs] = d_i(H a^{D-3} d0phi d_i phi)                 :", Z(U0 - U141 - div(Jsp)))
U142 = -R(1, 2)*H*a**(D - 3)*d0(gradsq)
print("  (142) exact                                                     :", Z(U141 - U142))
U143 = R(1, 2)*(D - 3)*H**2*a**(D - 2)*gradsq
print("  (143) (142) - [rhs] = d0(-1/2 H a^{D-3} grad phi^2)              :", Z(U142 - U143 - d0(-R(1, 2)*H*a**(D - 3)*gradsq)))

print("== completeness: Euler-Lagrange operators of all admissible monomials lie in span{B1..B5} (D = 4, 5, 6) ==")
def EL(L):
    return sp.expand(euler_equations(L, [phi], X)[0].lhs)
def run_complete(Dv):
    aD = lambda p: a**(Dv + p)
    boxN = lambda f: (-sp.diff(a**(Dv - 2)*d0(f), t) + a**(Dv - 2)*lap(f))/a**Dv
    basis = {'B1': aD(0)*boxN(phi)**2, 'B2': H**2*aD(-2)*etadd, 'B3': aD(-2)*boxN(phi)*lap(phi),
             'B4': aD(-4)*lap(phi)**2, 'B5': H**2*aD(-2)*gradsq}
    ELb = [EL(v) for v in basis.values()]
    cands = {
     '(a^-2 d0^2 phi)^2':           aD(0)*(d0(d0(phi))/a**2)**2,
     '(a^-2 d0^2)(a^-2 Lap)':       aD(0)*(d0(d0(phi))/a**2)*(lap(phi)/a**2),
     '(a^-2 Lap phi)^2':            aD(0)*(lap(phi)/a**2)**2,
     'a^-4 (d0 d_i phi)^2':         aD(-4)*sum(sp.diff(phi, t, xi)**2 for xi in xs),
     'a^-4 (d_i d_j phi)^2':        aD(-4)*sum(sp.diff(phi, xi, xj)**2 for xi in xs for xj in xs),
     'a^-4 d0phi d0^3 phi':         aD(-4)*d0(phi)*sp.diff(phi, t, 3),
     'a^-4 d0phi d0 Lap phi':       aD(-4)*d0(phi)*d0(lap(phi)),
     'a^-4 d_i phi d_i d0^2 phi':   aD(-4)*sum(sp.diff(phi, xi)*sp.diff(phi, xi, t, t) for xi in xs),
     'a^-4 d_i phi d_i Lap phi':    aD(-4)*sum(sp.diff(phi, xi)*lap(sp.diff(phi, xi)) for xi in xs),
     'a^D (a^-1 d0 (a^-1 d0 phi))^2': aD(0)*(d0(d0(phi)/a)/a)**2,
     'H a^{D-3} d0phi d0^2 phi':    H*aD(-3)*d0(phi)*d0(d0(phi)),
     'H a^{D-3} d0phi Lap phi':     H*aD(-3)*d0(phi)*lap(phi),
     'H a^{D-3} d_i phi d0 d_i phi': H*aD(-3)*sum(sp.diff(phi, xi)*sp.diff(phi, xi, t) for xi in xs),
     'H a^{D-1} d0phi Box phi':     H*aD(-1)*d0(phi)*boxN(phi),
     'H^2 a^{D-2} (d0 phi)^2':      H**2*aD(-2)*d0(phi)**2,
    }
    cs = sp.symbols('c1:6')
    def eqs_of(expr):
        expr = sp.expand(expr)
        ders = sorted(expr.atoms(sp.Derivative) | {phi}, key=lambda d: -sp.count_ops(d))
        qs = [sp.Symbol(f'Q{i}') for i in range(len(ders))]
        e2 = sp.expand(expr.xreplace(dict(zip(ders, qs)))*t**(Dv + 4))
        out = []
        for cf in sp.Poly(e2, *qs).coeffs():
            cf = sp.expand(cf)
            out += sp.Poly(cf, t).coeffs() if cf.has(t) else [cf]
        return out
    sol0 = sp.solve(eqs_of(sum(c*e for c, e in zip(cs, ELb))), cs, dict=True)
    indep = (sol0 == [] or sol0 == [{c: 0 for c in cs}])
    print(f"  D={Dv}: basis EL operators linearly independent: {indep}")
    for name, Lc in cands.items():
        sol = sp.solve(eqs_of(EL(Lc) - sum(c*e for c, e in zip(cs, ELb))), cs, dict=True)
        print(f"    {name:34s} -> c(B1..B5) = {[sol[0].get(c, 0) for c in cs] if sol else 'NOT in span'}")
for Dv in (4, 5, 6):
    run_complete(Dv)

print("== (145): vertex of -1/2 alpha3 kappa^2 H^2 a^{D-2} grad phi^2 ==")
al3, ka = sp.symbols('alpha3 kappa')
L145 = -R(1, 2)*al3*ka**2*H**2*a**(D - 2)*gradsq
print("  E-L operator (= delta^2 S acting on phi) == alpha3 kappa^2 H^2 a^{D-2} Lap phi :", Z(EL(L145) - al3*ka**2*H**2*a**(D - 2)*lap(phi)))
print("  equivalently -1/2 alpha3 kappa^2 H^2 a^D (grad/a phi)^2 :", Z(L145 - (-R(1, 2)*al3*ka**2*H**2*a**D*sum((sp.diff(phi, xi)/a)**2 for xi in xs))))

print("== remark after (145): S = int c(eta) d0phi Lap phi gives no d0 in its vertex ==")
c = sp.Function('c')(t)
print("  E-L of c d0phi Lap phi == -c' Lap phi :", Z(EL(c*d0(phi)*lap(phi)) - (-sp.diff(c, t)*lap(phi))))
