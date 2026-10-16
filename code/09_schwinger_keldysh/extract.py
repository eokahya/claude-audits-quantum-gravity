"""Coefficients quoted in the notes: Table 5'' (true sign) and the late-time form of R(eta,k), through k^10.

Purpose
    (i)   per-row non-falling coefficients of k^0, k^2, k^3, k^4, k^6 in the DEFINITION SIGN of (I.168)
          (= -1 x Table-X convention)  -> Table 5'' of the notes;
    (ii)  totals of the non-falling part, k^0 ... k^10 (must all vanish);
    (iii) true-sign a^-2 and a^-3 coefficients of the corrected total, k^0 ... k^10  -> Eq. for R in the notes;
    (iv)  printed table: non-falling part, and the true-sign difference (printed - corrected): its
          a^-1, a^-2, a^-3 coefficients per k^n (the 1/a tails removed by the zeta correction).

Expected result
    (ii) all zero; (iii) k^0: 12, k^2: -8 ln a - 103/3, k^3: -8i (a^-2) and -12i (a^-3), k^4: 8/3,
    k^5: 32i/45, k^6: -2/15, ...; (iv) non-falling all zero; a^-1 tails k^4: 4/45, k^6: 4/175,
    k^7: 2i/135, k^8: -8/1323, ...; a^-2 differences all zero.

Runtime: about 25-40 s.
Origin: Task 18 script extract.py (recomputes the N = 10 rows instead of loading a pickle).
"""
import sympy as sp
from engine import truncate, La, a
from rows import ROWS, paper_rows
from kexp import KExpansion, k


def nice(e):
    return sp.latex(sp.simplify(sp.expand(e).subs(La, sp.Symbol('\\ln a'))))


if __name__ == '__main__':
    E = KExpansion(10)
    rows = {i: E.row_rhs(r) for i, r in enumerate(ROWS)}
    print("=== per-row non-falling coefficients, TRUE sign (= - Table-X convention) ===")
    for i, r in rows.items():
        out = []
        for n in (0, 2, 3, 4, 6):
            c = -truncate(r.coeff(k, n), 0)
            out.append(f"k^{n}: {nice(c)}")
        print(i, ROWS[i][0], ' | '.join(out))
    tot = sp.expand(sum(rows.values()))
    print("totals non-falling k^0..k^10:", [sp.simplify(truncate(tot.coeff(k, n), 0)) for n in range(11)])
    print("=== TRUE-sign falling a^-2 coefficient of corrected total, k^0..k^10 ===")
    for n in range(11):
        c = sp.expand(-tot.coeff(k, n))
        c2 = sp.expand(truncate(c, 2) - truncate(c, 1))
        c3 = sp.expand(truncate(c, 3) - truncate(c, 2))
        print(f"  k^{n}: a^-2: {nice(c2*a**2)}   a^-3: {nice(c3*a**3)}")
    # printed table, full N10
    totp = sp.expand(sum(E.row_rhs(r) for r in paper_rows()))
    print("printed: non-falling", [sp.simplify(truncate(totp.coeff(k, n), 0)) for n in range(11)])
    print("printed minus corrected, TRUE sign, a^-1 coefficient per k^n:")
    for n in range(11):
        d = sp.expand(-(totp - tot).coeff(k, n))
        print(f"  k^{n}: a^-1: {nice(sp.expand(truncate(d, 1)-truncate(d, 0))*a)}  "
              f"a^-3: {nice(sp.expand(truncate(d, 3)-truncate(d, 2))*a**3)}  "
              f"a^-2: {nice(sp.expand(truncate(d, 2)-truncate(d, 1))*a**2)}")
