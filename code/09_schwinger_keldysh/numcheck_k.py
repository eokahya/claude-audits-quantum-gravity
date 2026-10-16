"""Check of the k-expansion (kexp.py) against direct quadrature with the exact k dependence.

Purpose
    Validate the k-expanded SK integrals
        G_q[g](a;k) = -(i/16 pi^2) int d eta' a'^q [u(eta',k)/u(0,k)] int d^3r [sin(kr)/(kr)] Delta g
    truncated at k^8, against direct mpmath quadrature with the exact u(eta',k)/u(0,k) = (1 - i k/a') e^{ik/a'}
    and the exact sin(kr)/(kr), at a = 200 and k/H = 0.3, 1.0.

Expected result
    Differences of the size of the k^9 truncation error: ~1e-10 - 2e-8 at k = 0.3, ~5e-6 - 1e-3 at k = 1
    (Task 18 log).

Runtime: about 5 s.
Origin: Task 18 script numcheck_k.py (now builds KExpansion(8) instead of setting sys.argv).
"""
import mpmath as mp
import sympy as sp
from kexp import KExpansion, k
from engine import La, a

mp.mp.dps = 25
E = KExpansion(8)


def U(kk, ap):  # u(eta')/u(0), eta' = -1/a'
    w = kk/ap
    return (1 - 1j*w)*mp.exp(1j*w)


def direct(q, g, A, kk, j=0):
    def inner(ap):
        de = 1/ap - 1/A
        if g == 'ln':
            rad = (mp.sin(kk*de) - kk*de*mp.cos(kk*de))/kk**3            # int_0^de r^2 sinc(kr) dr
            val = mp.mpf(1)/2*rad
        elif g == 'ln2':
            X0 = A*ap*de**2/4
            rad = mp.quad(lambda r: r**2*mp.sin(kk*r)/(kk*r)*(mp.log(X0) + mp.log(1 - (r/de)**2)), [0, de])
            val = rad
        else:
            X0 = A*ap*de**2/4
            rad = mp.quad(lambda r: r**2*mp.sin(kk*r)/(kk*r)*(1 - (r/de)**2)**j, [0, de])
            val = mp.mpf(1)/2*(-1)**j*X0**j*rad
        return ap**(q - 2)*U(kk, ap)*val
    return mp.quad(inner, [1, 2, 10, A])


if __name__ == '__main__':
    for (q, g, j) in [(4, 'ln', 0), (3, 'ln2', 0), (2, 'xj', 2), (4, 'xj', 1)]:
        ser = E.G(q, g, j)
        for (A, kk) in [(200, 0.3), (200, 1.0)]:
            sv = complex(sp.N(ser.subs({La: sp.log(A), a: A, k: kk}), 20))
            dv = complex(direct(q, g, mp.mpf(A), mp.mpf(kk), j))
            print(f"q={q} {g:3s} j={j} a={A} k={kk}: series(k^8)={sv:.10e}  quadrature={dv:.10e}  |diff|={abs(sv-dv):.1e}")
