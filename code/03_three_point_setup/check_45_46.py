"""Task 6 -- checks of Eqs. (I.45)-(I.46): set-up of the 3-point contribution to the scalar self-mass.

Paper I = E. O. Kahya and R. P. Woodard, Phys. Rev. D 76, 124005 (2007).

Purpose
    Supplementary code checks of the analytic derivation in Chapter 3, Section "Task 6" of the notes.
    They do not replace the derivation.

What is checked
    (a) Index algebra (I.45) -> (I.46).  In D = 3 a generic tensor
        X_{ab rs}(x;x') that is symmetric within each index pair (36 independent symbols; no pair-exchange
        symmetry is assumed because x != x') stands for the graviton propagator.  Derivatives are formal
        vectors: P, P' = outer-leg derivatives at x, x' (act on everything), Q, Q' = inner-leg derivatives
        (act only on i Delta_A).  With the derivative rules d_1 -> -P, d_2 -> +Q, d'_1 -> -P', d'_2 -> +Q'
        the sum  Sum_IJ V_I(x) V_J(x') X  is compared with the four terms of (I.46).
        Also: with the sign of the printed left-hand side of (I.7) (V_1 -> -V_1) the identity fails and the
        difference is exactly the sign flip of the 2nd and 3rd terms.
    (b) Overall normalization  (1/2) * 2 * 2 * (1/2i)^2 * i^2 * 2 = 1.
    (c) Derivative rule with test functions: a 1D toy  i Gamma = (1/2) int int W(z,z') Phi'(z) Phi'(z')
        with explicit smooth symmetric kernels.  The Gateaux second derivative along test functions f1, f2
        is compared with the kernel obtained from the rule of (I.45) [(-d_z)(-d_z') on everything], and with
        a rule that takes one outer leg with the opposite sign (the type of rule used in (I.32)).

Expected result
    (a) True, and False with difference = -2*(t2+t3) True;  (b) factor = 1;
    (c) 0.779492709611 for the Gateaux derivative and the (I.45) rule, -0.779492709611 for the wrong rule.

Runtime
    About 3 minutes (dominated by the two-dimensional mpmath quadratures in (c)).

Origin
    Run as check_45_46.py in Task 6 (15 Oct 2026).  The original one-line header comment was replaced by this
    docstring; the code is unchanged.
"""
import sympy as sp, mpmath as mp, itertools
I, R = sp.I, sp.Rational

print("== (a) index algebra (45) -> (46), D=3, generic within-pair-symmetric X_{ab rs}(x;x') ==")
D = 3; Rg = range(D)
eta = sp.diag(-1, 1, 1); ei = eta.inv()
pairs = [(i, j) for i in Rg for j in range(i, D)]
Xs = {(P, Q): sp.Symbol(f'X{P}{Q}') for P in range(len(pairs)) for Q in range(len(pairs))}   # no pair-exchange symmetry assumed (x != x')
def X(a, b, r, s): return Xs[(pairs.index(tuple(sorted((a, b)))), pairs.index(tuple(sorted((r, s)))))]
P  = sp.symbols('P0:3');  Q  = sp.symbols('Q0:3')     # P: outer derivative at x (acts on everything); Q: inner derivative at x (on iDelta_A)
Pp = sp.symbols('Pp0:3'); Qp = sp.symbols('Qp0:3')   # same at x'
up = lambda v: [sum(ei[m, n]*v[n] for n in Rg) for m in Rg]
dot = lambda u, v: sum(ei[m, n]*u[m]*v[n] for m in Rg for n in Rg)
Pu, Qu, Ppu, Qpu = up(P), up(Q), up(Pp), up(Qp)
def V(s1, d1u, d2u, d1, d2):   # Table I (kappa a^{D-2} stripped); d1 = outer leg, d2 = inner leg
    V1 = lambda a, b: s1*I*d1u[a]*d2u[b]
    V2 = lambda a, b: -I/2*ei[a, b]*dot(d1, d2)
    return lambda a, b: V1(a, b) + V2(a, b)
# derivative rules derived in Sec. 2: d_1 -> -P (all), d_2 -> +Q (iDelta_A), d'_1 -> -P', d'_2 -> +Q'
mP, mPu, mPp, mPpu = [-v for v in P], [-v for v in Pu], [-v for v in Pp], [-v for v in Ppu]
def lhs45(s1):
    Vx  = V(s1, mPu, Qu, mP, Q)
    Vxp = V(s1, mPpu, Qpu, mPp, Qp)
    return sp.expand(sum(Vx(a, b)*Vxp(r, s)*X(a, b, r, s) for a, b, r, s in itertools.product(Rg, Rg, Rg, Rg)))
# (46) with kappa^2 (aa')^{D-2} stripped: outer derivatives P, P' written outside, inner Q, Q' inside
t1 = -sum(Pu[a]*Ppu[r]*X(a, b, r, s)*Qu[b]*Qpu[s] for a, b, r, s in itertools.product(Rg, Rg, Rg, Rg))
t2 = R(1, 2)*sum(dot(P, Q)*Ppu[r]*ei[a, b]*X(a, b, r, s)*Qpu[s] for a, b, r, s in itertools.product(Rg, Rg, Rg, Rg))
t3 = R(1, 2)*sum(Pu[a]*dot(Pp, Qp)*X(a, b, r, s)*ei[r, s]*Qu[b] for a, b, r, s in itertools.product(Rg, Rg, Rg, Rg))
t4 = -R(1, 4)*dot(P, Q)*dot(Pp, Qp)*sum(ei[a, b]*ei[r, s]*X(a, b, r, s) for a, b, r, s in itertools.product(Rg, Rg, Rg, Rg))
rhs46 = sp.expand(t1 + t2 + t3 + t4)
print("  Sum_IJ V_I(x) V_J(x') X  ==  (46)  with derived V_1 = +i kappa ...  :", sp.expand(lhs45(+1) - rhs46) == 0)
d = sp.expand(lhs45(-1) - rhs46)
print("  with the sign of (7)'s printed LHS (V_1 = -i kappa ...)            :", d == 0,
      "  (difference = -2*(t2+t3):", sp.expand(d + 2*(t2 + t3)) == 0, ")")

print("== (b) normalization: (1/2) * 2 * 2 * (1/2i)^2 * i^2 * 2 ==")
print("  factor =", R(1, 2)*2*2*(1/(2*I))**2*I**2*2)

print("== (c) derivative rules, 1D toy with explicit smooth kernels and test functions ==")
mp.mp.dps = 20
z, zp = sp.symbols('z zp', real=True)
C = lambda u: 1 + u**2
G = (1 + z*zp)*sp.exp(-(z - zp)**2/2)                  # stands for the graviton propagator (symmetric)
Dl = sp.exp(-(z - zp)**2)                              # stands for iDelta_A (symmetric)
W = C(z)*C(zp)*G*sp.diff(Dl, z, zp)                    # coefficient multiplying Phi'(z) Phi'(z')
f1 = sp.exp(-z**2); f2 = z*sp.exp(-(z - 1)**2)
sub = lambda e, u: e.subs(z, u)
# Gateaux: i Gamma = (1/2) int int W Phi'(z) Phi'(z')  ->  d^2/de1 de2 = (1/2) int int W [f1'(z) f2'(z') + f2'(z) f1'(z')]
g1 = sp.diff(f1, z); g2 = sp.diff(f2, z)
A_int = R(1, 2)*W*(g1*g2.subs(z, zp) + g2*g1.subs(z, zp))
# rule of (45): outer derivatives -> (-d_x)(-d_x') acting on everything
K_rule = sp.diff(W, z, zp)
K_4pt_style = -sp.diff(W, z, zp)                       # one outer derivative with the opposite sign (as in the paper's (32))
B_int = f1*K_rule*f2.subs(z, zp)
Bw_int = f1*K_4pt_style*f2.subs(z, zp)
num = lambda e: mp.quad(sp.lambdify((z, zp), e, 'mpmath'), [-mp.inf, mp.inf], [-mp.inf, mp.inf])
A, B, Bw = num(A_int), num(B_int), num(Bw_int)
print("  Gateaux  d^2(i Gamma)/de1 de2 =", mp.nstr(A, 12))
print("  (45) rule  int f1 K f2        =", mp.nstr(B, 12), "  equal:", abs(A - B) < mp.mpf('1e-12'))
print("  one outer leg with + sign     =", mp.nstr(Bw, 12), "  equal:", abs(A - Bw) < mp.mpf('1e-12'))
