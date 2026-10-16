"""Late-time coefficients of R(eta,k), by power of k and of 1/a, for the corrected and printed tables.

Purpose
    Decompose the full SK right-hand side (through k^6, falling terms through a^-4) into
    coefficients  k^n a^-m  for
        tot   = corrected Table 1,
        tot_p = printed Table 1.
    This shows that the corrected total has no a^-1 term (leading fall-off ln a / a^2), while the
    printed table leaves k^4/a and k^6/a terms coming from the zeta row.

Convention
    TABLE-X (paper) CONVENTION; the notes quote these numbers multiplied by -1 (definition sign).

Expected result (Table-X convention)
    tot   : k^0 a^-2: -12, k^2 a^-2: (24 L_a + 103)/3, k^3 a^-2: 8i, k^3 a^-3: 12i, k^4 a^-2: -8/3,
            k^5 a^-2: -32i/45, k^6 a^-2: 2/15 (plus a^-4 terms); no a^-1 terms.
    tot_p : the same plus k^4 a^-1: -4/45, k^6 a^-1: -4/175 (and modified a^-3, a^-4 terms).

Runtime: about 10 s.
Origin: Task 18 script analysis3.py (recomputes the totals via analysis2.compute() instead of loading
    a pickle).
"""
import sympy as sp
from engine import truncate, a
from kexp import k
from analysis2 import compute

if __name__ == '__main__':
    E, tot, tot_p, zr = compute(6)
    d = {'tot': tot, 'tot_p': tot_p}
    for name in ('tot', 'tot_p'):
        print(f"== {name} ==")
        for n in range(7):
            c = sp.expand(d[name].coeff(k, n))
            for m in range(0, 5):
                cm = sp.expand(truncate(c, m) - truncate(c, m - 1))
                if cm != 0:
                    print(f"   k^{n}  a^-{m}: {sp.factor(sp.simplify(cm*a**m))}")
