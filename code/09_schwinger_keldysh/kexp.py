"""Beyond simplification (II.44): the full SK right-hand side of (II.29), expanded exactly in k.

Purpose
    Keep everything that (II.44) throws away -- the exact tree mode function u(eta',k), the phase
    e^{i k.(x'-x)}, the k^2/a^2 part of Box, and all rows of Table 1 that contain nabla^2 -- and
    compute, row by row, the non-falling (a -> infinity) part of the right-hand side R(eta,k) of
    (II.29) through order k^NMAX.

Method (H = 1, khat = k/H written as k)
    u(eta',k)/u(0,k) = (1 - i w) e^{i w},  w = k/a'          (series in w)
    angular average of e^{i k.r} = sin(kr)/(kr)                (series in (kr)^2)
    Box/H^2 -> D_k = Dop - k^2/a^2 ,  nabla^2 -> -k^2
    1/x -> D_k ln x ,  ln x/x -> D_k(ln^2 x/2 - ln x) + 3 ln x   (Eqs. (II.35), (II.36); constants drop)
    analytic f (const, x, x^2) drop out of f(y++) - f(y+-).
    Each term is reduced to the closed-form engine integrals with shifted a'-power (nu) and r^{2m}
    weight (m).  The integrand is entire in k on a compact region, so the expansion converges for
    every fixed a and k.

Convention
    TABLE-X (paper) CONVENTION; the definition sign of (I.168) multiplies everything by -1.

Usage
    python3 kexp.py N        (N = maximal power of k; default 4)

Expected result
    'TOTAL (Table-X convention), non-falling part': 0 at every order k^0 ... k^N.
    Rows 10 and 11 separately grow like +-(4/15) k^2 a (they cancel each other); row 12 (corrected
    zeta row) gives +1/70 at k^4.

Runtime: N = 2: ~2 s, N = 4: ~5 s, N = 6: ~10 s, N = 10: ~25 s.
Origin: Task 18 script kexp.py, wrapped in a class so that other scripts can choose NMAX without
    manipulating sys.argv; the mathematics is unchanged.
"""
import sys
import sympy as sp
from engine import W, Dop, truncate, La, a
from rows import ROWS, paper_rows

k = sp.Symbol('k', positive=True)
w = sp.Symbol('w')
ORDER = 4


class KExpansion:
    def __init__(self, nmax):
        self.NMAX = nmax
        Useries = sp.series((1 - sp.I*w)*sp.exp(sp.I*w), w, 0, nmax + 1).removeO()   # u(eta')/u(0), w = k/a'
        self.cU = {n: Useries.coeff(w, n) for n in range(nmax + 1)}
        self.cS = {m: sp.Rational((-1)**m, sp.factorial(2*m + 1)) for m in range(nmax//2 + 1)}   # sin(kr)/(kr)
        self._G = {}

    def ktrunc(self, f):
        f = sp.expand(f)
        return sum(f.coeff(k, n)*k**n for n in range(self.NMAX + 1))

    def G(self, q, g, j=0):
        """k-expanded -(i/16pi^2) int d^4x' a'^q [u(eta')/u(0)] e^{ik.r} Delta g, g in {'ln','ln2','xj'}"""
        key = (q, g, j)
        if key in self._G:
            return self._G[key]
        res = 0
        for nu in range(self.NMAX + 1):
            for m in range(self.NMAX//2 + 1):
                if nu + 2*m > self.NMAX or self.cU[nu] == 0:
                    continue
                res += self.cU[nu]*self.cS[m]*k**(nu + 2*m)*W(q, g, nu=nu, m=m, j=j, order=ORDER)
        self._G[key] = sp.expand(res)
        return self._G[key]

    def Dk(self, f):
        """Box/H^2 acting on F(eta) e^{-ik.x}"""
        return self.ktrunc(Dop(f) - k**2*f/a**2)

    def F(self, q, f):
        if f == 'lnx':
            return self.G(q, 'ln')
        if f == 'ln2x':
            return self.G(q, 'ln2')
        if f == '1/x':
            return self.Dk(self.G(q, 'ln'))                                            # (35), constant drops
        if f == 'lnx/x':
            return self.ktrunc(self.Dk(self.G(q, 'ln2')/2 - self.G(q, 'ln')) + 3*self.G(q, 'ln'))   # (36)
        if f == 'xlnx':
            return self.G(q, 'xj', 1)
        if f == 'x2lnx':
            return self.G(q, 'xj', 2)
        return 0                                         # analytic: const, x, x^2 drop out of ++ minus +-

    def row_rhs(self, row, fsel=None):
        lab, pq, nb, nl, fs = row
        tot = 0
        for f, c in fs.items():
            if fsel and f not in fsel:
                continue
            for (p, q) in pq:
                v = self.F(q, f)
                for _ in range(nb):
                    v = self.Dk(v)
                v = self.ktrunc(v*(-k**2)**nl)
                tot += c*a**(p - 4)*v
        return sp.expand(tot)


def main(nmax):
    E = KExpansion(nmax)
    res = {}
    for i, row in enumerate(ROWS):
        r = E.row_rhs(row)
        res[i] = r
        print(f"[{i:2d}] {row[0]}")
        for n in range(nmax + 1):
            cn = truncate(r.coeff(k, n), 0)
            if cn != 0:
                print(f"      k^{n}: {sp.simplify(cn)}")
    tot = sp.expand(sum(res.values()))
    print("TOTAL (Table-X convention), non-falling part:")
    for n in range(nmax + 1):
        print(f"  k^{n}: {sp.simplify(truncate(tot.coeff(k, n), 0))}")
    return E, res


if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4)
