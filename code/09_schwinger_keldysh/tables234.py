"""Tables 2-4 of Paper II: the twelve SK integrals S_K[f], K = 3, 4, 5, f = 1/x, ln x/x, ln x, ln^2 x.

Purpose
    Re-derive, without using the tables as input, the entries
        S_K[f] = -(i H^4/16 pi^2) int d^4x' a'^K { f(y++/4) - f(y+-/4) }      (k = 0, H = 1)
    and compare them with the printed Tables 2 (K = 3), 3 (K = 4) and 4 (K = 5) of Paper II.

Method
    ln x, ln^2 x : closed-form engine (engine.W), Eqs. (II.39), (II.40), (II.42)-(II.43) at k = 0.
    1/x          : identity (II.35), 1/x = Box ln x + 3  ->  S_K[1/x] = Dop S_K[ln x]  (constant drops).
    ln x/x       : identity (II.36)  ->  S_K[ln x/x] = Dop(S_K[ln^2 x]/2 - S_K[ln x]) + 3 S_K[ln x].
    The printed entries are truncated where the paper truncates (through O(1/a) for K = 3, O(a^0) for
    K = 4, 5); the next order of our exact expansion is also printed.
    (The direct light-cone evaluation of S_K[1/x] with Delta(1/x) = -8 pi i delta(y) is done by hand in
    the notes, not in this script.)

Expected result
    'paper - derived = 0' for all 12 entries.

Runtime: about 1 s.
Origin: Task 18 script tables234.py (pickle output removed; S() is reused by table5.py).
"""
import sympy as sp
from engine import W, Dop, truncate, clean, La, a

L = La


def S(K, f, order=3):
    """S_K[f] = -(i/16pi^2) int d^4x' a'^K {f(y++/4)-f(y+-/4)} at k = 0  (paper Tables 2-4)"""
    if f == 'lnx':
        return W(K, 'ln', order=order)
    if f == 'ln2x':
        return W(K, 'ln2', order=order)
    if f == '1/x':
        return Dop(W(K, 'ln', order=order))                                      # (35)
    if f == 'lnx/x':
        return sp.expand(Dop(W(K, 'ln2', order=order)/2 - W(K, 'ln', order=order))
                         + 3*W(K, 'ln', order=order))                            # (36)
    raise ValueError(f)


paper = {
    3: {'1/x': -L/a + 1/a,
        'lnx/x': -L**2/(2*a) + 2*L/a - 3/a + sp.pi**2/(3*a),
        'lnx': sp.Rational(1, 6) - L/(2*a) + 1/(4*a),
        'ln2x': L/3 - sp.Rational(11, 9) - L**2/(2*a) + 2*L/a - sp.Rational(9, 4)/a + sp.pi**2/(3*a)},
    4: {'1/x': -sp.Rational(1, 2),
        'lnx/x': sp.Rational(3, 4),
        'lnx': L/6 - sp.Rational(11, 36),
        'ln2x': L**2/6 - sp.Rational(8, 9)*L + sp.Rational(7, 4) - sp.pi**2/9},
    5: {'1/x': -a/6,
        'lnx/x': sp.Rational(17, 36)*a,
        'lnx': a/24 - sp.Rational(1, 6),
        'ln2x': -sp.Rational(13, 144)*a - L/3 + sp.Rational(5, 9)}}
keep = {3: 1, 4: 0, 5: 0}

if __name__ == '__main__':
    for K in (3, 4, 5):
        for f in ('1/x', 'lnx/x', 'lnx', 'ln2x'):
            full = clean(S(K, f))
            mine = truncate(full, keep[K])
            pk = truncate(sp.expand(paper[K][f]), keep[K])
            print(f"K={K} {f:6s}: derived = {sp.simplify(mine)}\n            paper - derived = {sp.simplify(pk - mine)}")
            print(f"            next order: {truncate(full, keep[K]+1) - mine}")
