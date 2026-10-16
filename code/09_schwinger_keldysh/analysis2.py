"""Printed versus corrected Table 1 in the full (beyond-(II.44)) SK right-hand side; the zeta correction.

Purpose
    (i)   non-falling part of R(eta,k) for the PRINTED Table 1, order by order in k (through k^6);
    (ii)  the contribution of the zeta-row correction  (aa')^2 nabla^4 (8 ln x - 16 x ln x)
          (corrected minus printed row 12), keeping falling terms through a^-3;
    (iii) the total for the CORRECTED Table 1 through O(1/a).
    The 72 gamma_2 correction (a constant) drops out of f(y++) - f(y+-) identically and needs no code.

Convention
    TABLE-X (paper) CONVENTION; multiply by -1 for the definition sign of (I.168).

Expected result
    (i)  [0, 0, 0, 0, 0, 0, 0];
    (ii) k^4: 4(a^2 - 15)/(45 a^3) = 4/(45a) - 4/(3a^3),  k^6: 2(6a^2 - 35)/(525 a^3) = 4/(175a) - 2/(15a^3);
         i.e. only falling terms;
    (iii) 0 at every order k^0 ... k^6 through O(1/a)  (no 1/a term in the corrected total).

Runtime: about 10 s.
Origin: Task 18 script analysis2.py (pickle output replaced by the function compute(), reused by
    analysis3.py).
"""
import sympy as sp
from kexp import KExpansion, k
from engine import truncate
from rows import ROWS, paper_rows

NMAX = 6


def compute(nmax=NMAX):
    E = KExpansion(nmax)
    tot_p = sp.expand(sum(E.row_rhs(r) for r in paper_rows()))
    zrow = ("zeta-corr", [(2, 2)], 0, 2, {'lnx': 8, 'xlnx': -16})
    zr = E.row_rhs(zrow)
    tot = sp.expand(sum(E.row_rhs(r) for r in ROWS))
    return E, tot, tot_p, zr


if __name__ == '__main__':
    E, tot, tot_p, zr = compute()
    # (i) printed (paper) rows
    print("printed Table 1, non-falling part per k-order:",
          [sp.simplify(truncate(tot_p.coeff(k, n), 0)) for n in range(NMAX + 1)])
    # (ii) zeta-row correction  (aa')^2 nabla^4 (8 ln x - 16 x ln x), Table-X convention, keep falling terms
    print("zeta correction, Table-X convention:")
    for n in range(NMAX + 1):
        c = truncate(zr.coeff(k, n), 3)
        if c != 0:
            print(f"   k^{n}: {sp.simplify(c)}")
    # 72 gamma_2 (constant) : f(y++)-f(y+-) = 0 identically -> contributes nothing (by construction)
    # (iii) full total, corrected table, falling terms through 1/a, Table-X convention
    print("corrected Table X total, through O(1/a):")
    for n in range(NMAX + 1):
        c = truncate(tot.coeff(k, n), 1)
        print(f"   k^{n}: {sp.simplify(c)}")
