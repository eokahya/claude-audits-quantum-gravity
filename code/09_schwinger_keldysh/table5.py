"""Table 5 of Paper II: the 19 non-derivative-row contributions under simplification (II.44).

Purpose
    For every term c * a^p a'^q Box^n * f(x) of the corrected Table 1 without nabla^2, compute the
    non-falling part of
        c a^(p-4) Dop^n S_q[f]       (units H^4 u(0,k)/(16 pi^2); kappa^2 stripped off, as in Phi = u + kappa^2 Phi_1)
    i.e. -(u(0,k)/a^4) (i H^8/(4pi)^4) int d^4x' (Ext. Operator) x {f(y++/4) - f(y+-/4)},
    which is the definition of the entries of Table 5 of Paper II (Eq. (II.45)).

Convention
    TABLE-X (paper) CONVENTION, to compare with the printed Table 5.  In the definition sign of
    Paper I, Eq. (I.168), every entry is multiplied by -1 (the total stays 0).
    Analytic terms (constants, including the corrected 72 gamma_2 constant +18 of row 7) drop out of
    f(y++) - f(y+-) identically; rows with nabla^2 vanish under (II.44).

Expected result
    'paper-derived: 0' for all 19 entries; TOTAL (non-falling) = 0.

Runtime: about 1 s.
Origin: Task 18 script table5.py (now recomputes S_K[f] instead of loading a pickle).
"""
import sympy as sp
from engine import Dop, truncate, La, a
from rows import ROWS
from tables234 import S as S_K

L = La
_S = {}


def S(q, f):
    if (q, f) not in _S:
        _S[(q, f)] = sp.expand(sp.simplify(sp.expand(S_K(q, f))))
    return _S[(q, f)]


paper5 = [0, 0, 0, -18, 0, 0, 0, 2, 3, -9, 33 - 18*L,
          (-325 + 12*sp.pi**2)/27 + sp.Rational(14, 3)*L - sp.Rational(2, 3)*L**2,
          (85 - 12*sp.pi**2)/27 - 4*L + sp.Rational(2, 3)*L**2,
          sp.Rational(8, 9) - sp.Rational(2, 3)*L,
          (160 - 12*sp.pi**2)/27 - sp.Rational(10, 3)*L + sp.Rational(2, 3)*L**2,
          sp.Rational(55, 3) - 10*L, -15 + 18*L,
          (-91 + 12*sp.pi**2)/27 + sp.Rational(8, 3)*L - sp.Rational(2, 3)*L**2,
          -sp.Rational(80, 9) + sp.Rational(32, 3)*L]


def Dn(f, n):
    for _ in range(n):
        f = Dop(f)
    return f


if __name__ == '__main__':
    entries = []
    for (lab, pq, nb, nl, fs) in ROWS:
        if nl:
            continue
        for f, c in fs.items():
            if f == 'const':
                continue
            val = 0
            for (p, q) in pq:
                val += c*a**(p - 4)*Dn(S(q, f), nb)
            val = truncate(sp.expand(val), 0)
            entries.append((lab, f, c, sp.simplify(val)))
    tot = 0
    print("entry (Table-X convention) vs paper Table 5")
    for (e, pv) in zip(entries, paper5):
        lab, f, c, v = e
        tot += v
        print(f"  {lab:34s} {str(c)+'*'+f:12s}  derived: {str(sp.expand(v)):45s} paper-derived: {sp.simplify(pv - v)}")
    print("  TOTAL (non-falling):", sp.simplify(tot))
