"""check_4pt.py -- Task 3: the 4-point (seagull) contribution, Eqs. (I.31)-(I.33).

Purpose
    Supplementary computer-algebra checks of the analytic derivation of the
    kernel -iM^2_4pt(x;x') of Paper I (Kahya & Woodard, PRD 76, 124005 (2007)),
    in particular of the derivative-reading rule behind Eq. (I.32).
    The analytic derivation is in notes chapter "The 4-point contribution:
    Eqs. (I.31)-(I.44)", section "Task 3".

Checks
    (a) Lemma A (second functional derivative of F = int K^{mu nu} d_mu phi d_nu phi)
        with explicit Gaussian test functions in D=2, eta=diag(-1,1), and a
        deliberately non-symmetric polynomial K^{mu nu}.  The kernel
        -d_mu[(K+K^T)^{mu nu} d_nu delta] reproduces d^2F/de1 de2; the opposite
        sign gives minus that value.
    (b) Lattice Hessian of S = sum_z c_z (phi_{z+1}-phi_z)^2 on a periodic 1D
        lattice with N=7 sites and symbolic c_z: equal to -D^-(2c D^+) (the lattice
        version of -d(2c d delta)), not to +D^-(2c D^+); eigenvalues for c=1 >= 0.
    (c) Operator-level check of (I.33) in D=3 with a generic pair-symmetric
        coincident propagator X_{abrs}(x) (21 independent functions), generic a(eta)
        and test function f: the kernel (derived rule) equals the derived (I.33'),
        which is exactly -1 times the printed (I.33); the paper's rule plus the
        derived U2 = +i/4 reproduces the printed (I.33), the printed U2 = -i/4 does not.
    (d) Calibration with the paper's counterterm vertex (I.131), D=3: Lemma A
        reproduces (I.131); the paper's 4-point rule gives the opposite sign.

Method
    SymPy (exact symbolic integration / differentiation).

Expected result
    Every "expected" line prints True, the deliberate counter-examples print False,
    A = 97*sqrt(6)*pi*exp(-1/2)/54 ~ 8.384, eigenvalues [0.0, 1.506041, 4.890084, 7.603875].
    See expected_output/check_4pt.txt.

Runtime
    about 20-25 s (SymPy 1.14).

Provenance
    Reconstructed verbatim from the Task 3 command log (script build3/check_4pt.py,
    including the in-log edit of the eigenvalue print statement); only this
    docstring was added.
"""
# Supplementary checks for Eqs. (31)-(33): kernel of the 4-point (seagull) contribution.
import sympy as sp, itertools
I, R = sp.I, sp.Rational

print("== (a) Lemma A with explicit test functions, D=2, eta=diag(-1,1) ==")
t, x, e1, e2 = sp.symbols('t x e1 e2', real=True)
X = (t, x)
K = sp.Matrix([[1 + t**2, x], [t*x, 2 + x**2]])          # deliberately non-symmetric K^{mu nu}
f1 = sp.exp(-t**2 - x**2)
f2 = (1 + t + x**2)*sp.exp(-(t - 1)**2 - 2*x**2)
phi = e1*f1 + e2*f2
F_density = sum(K[m, n]*sp.diff(phi, X[m])*sp.diff(phi, X[n]) for m in range(2) for n in range(2))
A_density = sp.diff(F_density, e1, e2)                    # d^2 F / de1 de2 (integrand)
Ks = K + K.T
divK = sum(sp.diff(Ks[m, n]*sp.diff(f2, X[n]), X[m]) for m in range(2) for n in range(2))
B_ours_density = -f1*divK                                 # f1 * ( -d_mu[(K+K^T) d_nu f2] )
B_paper_density = +f1*divK                                # f1 * ( +d_mu[(K+K^T) d_nu f2] )
integ = lambda g: sp.simplify(sp.integrate(sp.expand(g), (t, -sp.oo, sp.oo), (x, -sp.oo, sp.oo)))
A, Bo, Bp = integ(A_density), integ(B_ours_density), integ(B_paper_density)
print("  A  = d2F/de1de2           =", A, "~", sp.N(A))
print("  B_ours  (kernel -d[K_s d delta]) =", sp.N(Bo), " equal to A:", sp.simplify(A - Bo) == 0)
print("  B_paper (kernel +d[K_s d delta]) =", sp.N(Bp), " equal to -A:", sp.simplify(A + Bp) == 0)

print("== (b) lattice Hessian, 1D periodic, N=7 ==")
N = 7
c = sp.symbols(f'c0:{N}'); p = sp.symbols(f'p0:{N}')
S = sum(c[z]*(p[(z + 1) % N] - p[z])**2 for z in range(N))
Hs = sp.hessian(S, p)
Dp = sp.zeros(N, N); Dm = sp.zeros(N, N)
for z in range(N):
    Dp[z, (z + 1) % N] += 1; Dp[z, z] -= 1               # forward difference
    Dm[z, z] += 1; Dm[z, (z - 1) % N] -= 1               # backward difference
C2 = sp.diag(*[2*ci for ci in c])
print("  Hessian == -D^-(2c D^+)  [lattice version of -d(2c d delta)]:", (Hs - (-Dm*C2*Dp)).expand() == sp.zeros(N, N))
print("  Hessian == +D^-(2c D^+)  [lattice version of +d(2c d delta)]:", (Hs - (Dm*C2*Dp)).expand() == sp.zeros(N, N))
ev = (Hs.subs({ci: 1 for ci in c})).eigenvals()
print("  eigenvalues for c=1 (must be >= 0 for a sum of squares):", sorted(round(float(sp.re(sp.N(k))), 6) for k in ev))

print("== (c) operator-level check of (33), D=3, generic pair-symmetric coincident propagator ==")
D = 3
xs = sp.symbols('x0:3', real=True)
eta = sp.diag(-1, 1, 1); ei = eta.inv()
a = sp.Function('a')(xs[0]); f = sp.Function('f')(*xs); kap = sp.Symbol('kappa')
pairs = [(i, j) for i in range(D) for j in range(i, D)]
Xf = {}
for P, Q in itertools.combinations_with_replacement(range(len(pairs)), 2):
    Xf[(P, Q)] = Xf[(Q, P)] = sp.Function(f'X{P}{Q}')(*xs)
def Xc(al, be, rh, si):                                   # X_{alpha beta rho sigma} = i[_{ab}Delta_{rs}](x;x)
    P = pairs.index(tuple(sorted((al, be)))); Q = pairs.index(tuple(sorted((rh, si))))
    return Xf[(P, Q)]
rng = range(D)
X1 = sum(ei[al, be]*ei[rh, si]*Xc(al, be, rh, si) for al in rng for be in rng for rh in rng for si in rng)
X2 = sum(ei[al, rh]*ei[be, si]*Xc(al, be, rh, si) for al in rng for be in rng for rh in rng for si in rng)
X3 = sp.Matrix(D, D, lambda m, n: sum(ei[al, be]*ei[rh, m]*ei[si, n]*Xc(al, be, rh, si) for al in rng for be in rng for rh in rng for si in rng))   # i[^a_a Delta^{mn}]
X4 = sp.Matrix(D, D, lambda m, n: sum(ei[al, m]*ei[be, rh]*ei[si, n]*Xc(al, be, rh, si) for al in rng for be in rng for rh in rng for si in rng))   # i[^m_r Delta^{r n}]
aD2 = a**(D - 2)
def Ktilde(u2):   # U_I contracted with X: coefficient of d_{1mu} d_{2nu}
    return (-I/8*kap**2*aD2*X1*ei + u2*kap**2*aD2*X2*ei + I/2*kap**2*aD2*X3 - I*kap**2*aD2*X4)
d = lambda g, m: sp.diff(g, xs[m])
def kernel_ours(Kt):   # -iM^2 acting on f:  -d_mu[ K_s^{mu nu} d_nu f ]
    Ks = (Kt + Kt.T)/2
    return -sum(d(Ks[m, n]*d(f, n), m) for m in rng for n in rng)
up = lambda g: [sum(ei[m, n]*d(g, n) for n in rng) for m in rng]           # d^mu g
def form33(s1, s2, s3, s4):   # s1 d^mu{a^{D-2} X1 d_mu f} + s2 d^mu{a X2 d_mu f} + s3 d_r{a X3^{rs} d_s f} + s4 d^a{a X4_a^s d_s f}
    T1 = sum(d(sum(ei[m, n]*aD2*X1*d(f, n) for n in rng), m) for m in rng)
    T2 = sum(d(sum(ei[m, n]*aD2*X2*d(f, n) for n in rng), m) for m in rng)
    T3 = sum(d(aD2*X3[m, n]*d(f, n), m) for m in rng for n in rng)
    T4 = sum(d(aD2*X4[m, n]*d(f, n), m) for m in rng for n in rng)          # X4[m,n] = i[^m_r Delta^{rn}] = d^alpha{i[_{alpha r}Delta^{r sigma}] d_sigma}
    return kap**2*(s1*T1 + s2*T2 + s3*T3 + s4*T4)
Kd = Ktilde(+I/4)
print("  K~_4 symmetric (pair symmetry of X):", sp.simplify(sp.expand((-I*kap**2*aD2*X4) - (-I*kap**2*aD2*X4).T)) == sp.zeros(D, D))
ours = sp.expand(kernel_ours(Kd))
der33 = sp.expand(form33(+I/8, -I/4, -I/2, +I))
pap33 = sp.expand(form33(-I/8, +I/4, +I/2, -I))
print("  derived kernel == derived (33)  [+i/8, -i/4, -i/2, +i]:", sp.expand(ours - der33) == 0)
print("  derived kernel == -(paper's 33) [-i/8, +i/4, +i/2, -i]:", sp.expand(ours + pap33) == 0)
paper_rule = lambda Kt: -kernel_ours(Kt)                  # d_1 -> -d (everything), d_2 -> +d' on delta
print("  paper's rule + derived U2(+i/4) reproduces paper's (33):", sp.expand(paper_rule(Kd) - pap33) == 0)
print("  paper's rule + printed U2(-i/4) reproduces paper's (33):", sp.expand(paper_rule(Ktilde(-I/4)) - pap33) == 0)

print("== (d) calibration with the paper's counterterm vertex (131), D=3 ==")
al2, Hh = sp.symbols('alpha2 H')
Kct = -R(1,2)*al2*kap**2*Hh**2*aD2*ei                     # Delta L = K^{mu nu} d_mu phi d_nu phi
box_f = sum(d(aD2*ei[m, n]*d(f, n), m) for m in rng for n in rng)/a**D
paper131 = I*al2*kap**2*Hh**2*a**D*box_f                   # paper: -> i alpha2 kappa^2 H^2 a^D Box delta
ours131 = I*(-sum(d((Kct + Kct.T)[m, n]*d(f, n), m) for m in rng for n in rng))   # i * delta^2 S
print("  i*delta^2 S (Lemma A) == paper's (131):", sp.simplify(ours131 - paper131) == 0)
Uct = 2*I*Kct                                              # vertex operator in the Table-II normalization
print("  paper's 4-pt rule applied to (131) == paper's (131):", sp.simplify(paper_rule(Uct) - paper131) == 0)
