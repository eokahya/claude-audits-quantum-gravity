"""Independent numerical check of the closed-form engine (engine.W) by direct quadrature.

Purpose
    Verify the exact a-expansions of the SK integrals (Eqs. (II.39)-(II.43) at the level of the
    radial z-integral and the a'-integral) against direct mpmath quadrature at finite a.

Method
    For a sample of (q, f, nu, m, j) -- including ln^2 x, x ln x, x^2 ln x, an r^2 weight (m = 1) and a
    u-correction power (nu = 2) -- evaluate
        -(i/16 pi^2) int d eta' a'^(q-nu) int d^3r r^(2m) {f(y++/4) - f(y+-/4)}
    (i) from engine.W expanded through a^-6 and (ii) by direct quadrature over a' and z (30 digits),
    at a = 50 and a = 400.

Expected result
    |difference| <~ 1e-19 for every case (Task 18 log: max 1.02e-19).

Runtime: a few seconds.
Origin: Task 18 script numcheck.py (unchanged).
"""
import mpmath as mp
import sympy as sp
from engine import W, La, a

mp.mp.dps = 30


def direct(q, ftype, A, nu=0, m=0, j=0):
    # -(i/16pi^2) int deta' a'^{q-nu} int d^3r r^{2m} Delta f, with the same prefactors as W but integrated numerically
    def inner(ap):
        de = 1/ap - 1/A
        X0 = A*ap*de**2/4
        if ftype == 'ln':
            zi = mp.mpf(1)/(3 + 2*m)
            pref = mp.mpf(1)/2
            val = pref*zi*de**(3 + 2*m)
        elif ftype == 'ln2':
            zi = mp.quad(lambda zz: zz**(2 + 2*m)*(mp.log(X0) + mp.log(1 - zz**2)), [0, 1])
            val = zi*de**(3 + 2*m)
        else:
            zi = mp.quad(lambda zz: zz**(2 + 2*m)*(1 - zz**2)**j, [0, 1])
            val = mp.mpf(1)/2*(-1)**j*zi*X0**j*de**(3 + 2*m)
        return ap**(q - nu - 2)*val            # d eta' = d a'/a'^2
    return mp.quad(inner, [1, mp.sqrt(A), A])


if __name__ == '__main__':
    for (q, ft, nu, m, j) in [(3, 'ln2', 0, 0, 0), (5, 'ln2', 0, 0, 0), (4, 'ln', 0, 0, 0), (2, 'xj', 0, 0, 2),
                              (4, 'xj', 0, 0, 1), (3, 'ln2', 2, 1, 0), (5, 'ln', 0, 1, 0)]:
        ex = W(q, ft, nu=nu, m=m, j=j, order=6)
        for A in (50, 400):
            e = sp.N(ex.subs({La: sp.log(A), a: A}), 20)
            d = direct(q, ft, mp.mpf(A), nu, m, j)
            print(f"q={q} {ft:4s} nu={nu} m={m} j={j} a={A:4d}: exact-expansion={float(e): .12e}  "
                  f"quadrature={float(d): .12e}  diff={float(e-d): .2e}")
