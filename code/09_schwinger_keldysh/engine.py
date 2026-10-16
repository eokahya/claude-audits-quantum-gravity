"""Exact late-time engine for the Schwinger-Keldysh (SK) integrals of Paper II (helper module).

Purpose
    Closed-form evaluation of the integrals that appear when the corrected Table 1 of
    Paper II (= corrected Table X of Paper I) is inserted into the SK mode equation (II.29).
    Units H = 1 throughout; ln(a) is represented by the symbol L_a.

Method
    Inside the past light cone y++ = -|y| + i0 and y+- = -|y| - i0, so (Eqs. (II.39)-(II.41))
        Delta ln x       = 2 pi i theta,
        Delta ln^2 x     = 4 pi i theta ln|x|,
        Delta x^j ln x   = 2 pi i theta x^j,      x -> -X0 (1 - z^2),
    with X0 = a a' Delta_eta^2 / 4, r = Delta_eta z.  After the substitution a' = a s one has
    Delta_eta = (1-s)/(a s) and X0 = (1-s)^2/(4 s), independent of a.  Every integral then reduces to
        int_{1/a}^1 ds s^alpha (1-s)^beta {1, ln s, ln(1-s)},
    which is done exactly (harmonic numbers, -pi^2/6); the result is an exact expansion in eps = 1/a.

Main functions
    W(q, ftype, nu, m, j, order)
        = -(i/16 pi^2) int d eta' a'^(q-nu) int d^3r r^(2m) { f(y++/4) - f(y+-/4) }
        for ftype in {'ln' (ln x), 'ln2' (ln^2 x), 'xj' (x^j ln x)}; with nu = m = j = 0 this is
        S_q[f] of the notes (Tables 2-4 of Paper II for q = 3, 4, 5).
    Dop(f)       Box/H^2 at k = 0 on functions of a: -(a^2 d^2/da^2 + 4 a d/da)   (Eq. (II.37), k = 0)
    truncate(f, n)  keep terms of order a^(>= -n) (powers of ln a do not count).

Expected result / runtime
    No output when imported; used by all other scripts in this folder.  Import time < 1 s.

Origin
    Task 18 script engine.py (scratch build18), unchanged except for the removal of an unused
    helper (Bz) and of absolute paths.
"""
import sympy as sp

a, La = sp.symbols('a L_a', positive=True)        # La stands for ln(a)
eps = sp.Symbol('epsilon', positive=True)
LN2 = sp.log(2)


def H(n):
    """harmonic number H_n"""
    return sum(sp.Rational(1, k) for k in range(1, n + 1))


def J1(p, order):
    """int_eps^1 s^p ds, exact"""
    if p == -1:
        return La                      # -ln eps = ln a
    return sp.Rational(1, p + 1) - eps**(p + 1)/(p + 1)


def Jlns(p, order):
    """int_eps^1 s^p ln s ds, exact (ln eps = -La)"""
    if p == -1:
        return -La**2/2
    return -sp.Rational(1, (p + 1)**2) + eps**(p + 1)*La/(p + 1) + eps**(p + 1)/sp.Rational((p + 1)**2)


def Jln1ms(p, order):
    """int_eps^1 s^p ln(1-s) ds, exact up to eps^order.
    C_p = -H_{p+1}/(p+1) (p >= 0), -pi^2/6 (p = -1), -H_q/q (p <= -2, q = -p-1)."""
    q = -p - 1 if p <= -2 else 0
    if p >= 0:
        C = -H(p + 1)/(p + 1)
    elif p == -1:
        C = -sp.pi**2/6
    else:
        C = -H(q)/q
    res = C
    for n in range(1, q + 1):                       # - sum_{n<=q} (1/n) int_eps^1 s^{p+n}
        res -= sp.Rational(1, n)*J1(p + n, order)
    n = q + 1                                       # tail  sum_{n>q} eps^{p+n+1}/(n(p+n+1))
    while p + n + 1 <= order:
        res += eps**(p + n + 1)/(n*(p + n + 1))
        n += 1
    return res


def Sint(alpha, beta, kind, order):
    """int_eps^1 s^alpha (1-s)^beta L ds, L in {'1','lns','ln1ms'} (binomial expansion of (1-s)^beta)"""
    res = 0
    for i in range(beta + 1):
        c = sp.binomial(beta, i)*(-1)**i
        p = alpha + i
        res += c*{'1': J1, 'lns': Jlns, 'ln1ms': Jln1ms}[kind](p, order)
    return sp.expand(res)


def to_a(expr, apow):
    """multiply by a^apow and substitute eps = 1/a"""
    return sp.expand(a**apow*expr.subs(eps, 1/a))


def eta_int(Q, N, j, G, order=3):
    """int_{eta_i}^{eta} d eta' a'^Q Dleta^N X0^j G(X0),  X0 = a a' Dleta^2/4 = (1-s)^2/(4s);  G in {'1','lnX0'}
       = a^{Q-1-N} int ds s^{Q-2-N} (1-s)^N X0^j G"""
    apow = Q - 1 - N
    alpha = Q - 2 - N - j
    beta = N + 2*j
    pref = sp.Rational(1, 4**j)
    o = order + apow + 2
    if G == '1':
        I = Sint(alpha, beta, '1', o)
    else:   # ln X0 = 2 ln(1-s) - ln s - 2 ln 2
        I = 2*Sint(alpha, beta, 'ln1ms', o) - Sint(alpha, beta, 'lns', o) - 2*LN2*Sint(alpha, beta, '1', o)
    return to_a(pref*I, apow)


z = sp.Symbol('z', positive=True)


def W(q, ftype, nu=0, m=0, j=0, order=3):
    """W = -(i/16 pi^2) int d eta' a'^{q-nu} int d^3r r^{2m} { f(y++/4) - f(y+-/4) }
       ftype: 'ln' (ln x), 'ln2' (ln^2 x), 'xj' (x^j ln x)"""
    Q = q - nu
    N = 3 + 2*m
    if ftype == 'ln':
        return sp.Rational(1, 2)/(3 + 2*m)*eta_int(Q, N, 0, '1', order)
    if ftype == 'ln2':
        Cm = sp.integrate(z**(2 + 2*m)*sp.log(1 - z**2), (z, 0, 1))
        return sp.expand(eta_int(Q, N, 0, 'lnX0', order)/(3 + 2*m) + Cm*eta_int(Q, N, 0, '1', order))
    if ftype == 'xj':
        B = sp.integrate(z**(2 + 2*m)*(1 - z**2)**j, (z, 0, 1))
        return sp.expand(sp.Rational(1, 2)*(-1)**j*B*eta_int(Q, N, j, '1', order))
    raise ValueError(ftype)


# ---- operators on functions of a (with La = ln a) ----
def d_a(f):
    """total derivative wrt a, with La = ln a"""
    return sp.expand(sp.diff(f, a) + sp.diff(f, La)/a)


def Dop(f):
    """Box/H^2 at k = 0 :  -(a^2 d^2 + 4 a d)"""
    return sp.expand(-(a**2*d_a(d_a(f)) + 4*a*d_a(f)))


def apower(t):
    c, e = t.as_coeff_exponent(a)
    return e


def truncate(f, n):
    """keep terms of order a^{>= -n} (log powers do not count)"""
    out = 0
    for t in sp.Add.make_args(sp.expand(f)):
        if apower(t) >= -n:
            out += t
    return sp.expand(out)


def clean(f):
    return sp.expand(sp.simplify(sp.expand(f)))
