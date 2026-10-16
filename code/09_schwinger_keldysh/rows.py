"""Corrected Table 1 of Paper II (= corrected Table X of Paper I) as data (helper module).

Purpose
    Input data for the Schwinger-Keldysh scripts of this folder.

Content
    ROWS: the 13 rows of the corrected nonlocal table, in the TABLE-X (paper) CONVENTION, i.e. with the
    sign of the printed table; header kappa^2 H^4/(4 pi)^4.  The physical self-mass-squared in the
    DEFINITION SIGN of Paper I, Eq. (I.168), is  -i M^2 = (local) - (this table).
    Each row: (label, [(p,q), ...], n_box, n_lap, {f: coeff}) where
        a^p a'^q            scale-factor prefactor ((aa')^3(a^2+a'^2) -> (5,3) + (3,5), etc.),
        n_box, n_lap        powers of Box and nabla^2 (acting on f(y) only; scale factors stand on the left),
        f in {'lnx/x','1/x','ln2x','lnx','xlnx','x2lnx','x','x2','const'}.
    Corrections with respect to the printed table (Tasks 13 and 17):
        row 7  (aa')^3(a^2+a'^2)H^4 : constant -54 -> +18  (72 gamma_2; analytic, drops out of SK)
        row 12 (aa')^2 nabla^4      : 10/3 ln x - 24 x ln x -> 34/3 ln x - 40 x ln x  (zeta row)
    paper_rows() returns the printed table.

Runtime: import only.
Origin: Task 18 script rows.py (unchanged).
"""
import copy
import sympy as sp

R = sp.Rational
# each row: (label, [(p,q),...] from the scale-factor prefactor, n_box, n_lap, {f: coeff})
ROWS = [
    ("(aa')^4 \\Box^3/H^2",             [(4, 4)],          3, 0, {'lnx/x': R(-1, 3)}),
    ("(aa')^4 \\Box^2",                 [(4, 4)],          2, 0, {'lnx/x': R(26, 3), '1/x': R(38, 3), 'ln2x': -6, 'lnx': -18}),
    ("(aa')^4 H^2\\Box",                [(4, 4)],          1, 0, {'lnx/x': -6, '1/x': 4, 'lnx': -4}),
    ("(aa')^4 H^4",                     [(4, 4)],          0, 0, {'lnx/x': 4, '1/x': 18, 'lnx': -108, 'const': -120}),
    ("(aa')^3(a^2+a'^2)\\Box^3/H^2",    [(5, 3), (3, 5)],  3, 0, {'lnx/x': R(1, 6)}),
    ("(aa')^3(a^2+a'^2)\\Box^2",        [(5, 3), (3, 5)],  2, 0, {'lnx/x': R(-1, 3), '1/x': R(1, 6)}),
    ("(aa')^3(a^2+a'^2)H^2\\Box",       [(5, 3), (3, 5)],  1, 0, {'lnx/x': R(-2, 3), '1/x': 5, 'lnx': -18}),
    ("(aa')^3(a^2+a'^2)H^4",            [(5, 3), (3, 5)],  0, 0, {'lnx/x': R(4, 3), '1/x': R(-32, 3), 'const': 18}),
    ("(aa')^3 H^2\\nabla^2",            [(3, 3)],          0, 1, {'lnx/x': R(-2, 3), '1/x': -16, 'lnx': 84, 'xlnx': -48, 'x': 96}),
    ("(aa')^3 \\nabla^2\\Box",          [(3, 3)],          1, 1, {'lnx/x': R(1, 3)}),
    ("(aa')^2(a^2+a'^2)H^2\\nabla^2",   [(4, 2), (2, 4)],  0, 1, {'lnx/x': R(-7, 3), '1/x': R(11, 2), 'lnx': -12, 'xlnx': 48, 'x': 12}),
    ("(aa')^2(a^2+a'^2)\\nabla^2\\Box", [(4, 2), (2, 4)],  1, 1, {'lnx/x': R(-17, 6), '1/x': R(-49, 6), 'ln2x': 4, 'lnx': 10, 'xlnx': 12}),
    ("(aa')^2\\nabla^4",                [(2, 2)],          0, 2, {'lnx': R(34, 3), 'xlnx': -40, 'x2lnx': 24, 'x2': -36}),
]
# printed values where they differ
PAPER_DIFF = {7: {'const': -54}, 12: {'lnx': R(10, 3), 'xlnx': -24}}


def paper_rows():
    rows = copy.deepcopy(ROWS)
    for i, d in PAPER_DIFF.items():
        rows[i][4].update(d)
    return rows
