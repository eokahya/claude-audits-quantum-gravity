"""check_table4.py -- Table IV of Paper I in general D  [Task 13, check (T4)]

Purpose
  Derive mechanically, in general D, the two columns of Table IV (coefficients of
  Gamma^2(D/2) kappa^2 H^(2D-4)/(4 pi)^D x (4/y)^(D-1) and x (4/y)^(D-2)) from the totals of the
  appendix Tables XI-XX, using only the exact-D leading coefficients of the eight products
  (I.116)-(I.123) and the rules  I[z^q] = -4 z^(q-1)/(q-1),  y z^q = 4 z^(q-1),
  d/dy z^q = -(q/4) z^(q+1)   (z = 4/y).  Compares with the printed Table IV.

Expected result
  19 of 20 entries agree; gamma2, second column: derived +(D^3-3D^2-4D+8)/(4(D-1)),
  printed -(D^3-3D^2-4D+8)/(4(D-1))  (printed as "True False" in the gamma2 line).
Runtime: a few seconds.
Run:     python3 check_table4.py
"""
# Table IV in general D from the appendix totals, using only the exact-D leading coefficients of the products.
import sympy as sp
D = sp.Symbol('D'); R = sp.Rational
# product -> {exponent of z=4/y : coefficient} (units Gamma^2(D/2) N^2/16)
P = {'A2': {D: 1, D - 1: D}, 'AApp': {D: D/(D - 2)}, 'ApBp': {D: 1, D - 1: D - 2}, 'BApp': {D: D/(D - 2)},
     'ApdBp': {D - 1: -2}, 'dBApp': {}, 'ApdCp': {D - 1: -8}, 'dCApp': {}}
def add(*fs):
    out = {}
    for c, f in fs:
        for q, v in f.items(): out[q] = sp.simplify(out.get(q, 0) + c*v)
    return out
I = lambda f: {sp.simplify(q - 1): sp.simplify(-4*v/(q - 1)) for q, v in f.items()}
Yx = lambda f: {sp.simplify(q - 1): 4*v for q, v in f.items()}
dY = lambda f: {sp.simplify(q + 1): sp.simplify(-q*v/4) for q, v in f.items()}
In = lambda n, f: f if n == 0 else In(n - 1, I(f))
T = {}
T['alpha'] = add((1, In(2, P['AApp'])), (1, In(2, P['dCApp'])))
T['beta'] = add((D, I(P['ApBp'])), (-1, I(P['ApdBp'])), (-R(1, 2), Yx(I(P['A2']))), (-(D - 1)/2, In(2, P['A2'])),
                (1, I(P['ApdCp'])), (-R(1, 2), Yx(I(P['ApdCp']))), (-(D - 1)/2, In(2, P['ApdCp'])))
T['gamma1'] = {q: sp.simplify((D - 1)*(D - 2)*v) for q, v in T['beta'].items()}
T['gamma2'] = add((D, Yx(P['ApBp'])), (-1, Yx(P['ApdBp'])), (-R(1, 2), Yx(Yx(P['A2']))), (-D/2, Yx(I(P['A2']))),
                  (1, Yx(P['ApdCp'])), (-D/2, Yx(I(P['ApdCp']))), (-R(1, 2), Yx(Yx(P['ApdCp']))))
T['gamma3'] = add(((D - 1)*D, Yx(P['ApBp'])), (D, Yx(Yx(dY(P['ApBp'])))), (-(D - 1), Yx(P['ApdBp'])), (-1, Yx(Yx(dY(P['ApdBp'])))),
                  (-D, Yx(Yx(P['A2']))), (-R(1, 2), Yx(Yx(Yx(dY(P['A2']))))), (-R(1, 2)*(D - 1)*D, Yx(I(P['A2']))),
                  (D - 1, Yx(P['ApdCp'])), (1, Yx(Yx(dY(P['ApdCp'])))), (-D, Yx(Yx(P['ApdCp']))), (-R(1, 2), Yx(Yx(Yx(dY(P['ApdCp']))))),
                  (-R(1, 2)*(D - 1)*D, Yx(I(P['ApdCp']))))
T['delta'] = add((2, In(2, P['dBApp'])), (R(1, 4), In(3, P['A2'])), (-1, In(2, P['dCApp'])), (R(1, 4), In(3, P['ApdCp'])))
T['eps1'] = add((-2*(D - 2), I(P['ApdBp'])), (-1, Yx(I(P['ApdBp']))), (-(5*D - 4)/2, In(2, P['ApdBp'])), (D**2/2, In(2, P['ApBp'])),
                (-D/4, Yx(In(2, P['A2']))), ((D - 2)**2, In(2, P['BApp'])), (-R(3, 4)*(D - 2), In(3, P['A2'])), (-2, I(P['ApdCp'])),
                (1, Yx(I(P['ApdCp']))), ((3*D - 2)/2, In(2, P['ApdCp'])), (-D/4, Yx(In(2, P['ApdCp']))), (-R(3, 4)*(D - 2), In(3, P['ApdCp'])))
T['eps2'] = add((2, Yx(I(P['BApp']))), (1, In(2, P['ApBp'])), (-(D - 1), In(2, P['ApdBp'])), (R(1, 4), Yx(In(2, P['A2']))), (R(1, 4), Yx(In(2, P['ApdCp']))))
T['eps3'] = add((2, Yx(Yx(P['BApp']))), (2*(D - 1), Yx(I(P['BApp']))), (D - 2, In(2, P['ApBp'])), (1, Yx(I(P['A2']))), ((D - 1)/4, Yx(In(2, P['A2']))),
                (R(1, 4), Yx(Yx(I(P['A2'])))), ((D - 1)/4, Yx(In(2, P['ApdCp']))), (R(1, 4), Yx(Yx(I(P['ApdCp'])))))
T['zeta'] = add((-4, In(2, P['dBApp'])), (1, In(3, P['ApdBp'])), ((D - 2)/8, In(4, P['A2'])), (1, In(2, P['dCApp'])), (-R(1, 2), In(3, P['ApdCp'])), ((D - 2)/8, In(4, P['ApdCp'])))
paper = {'alpha': (0, D/((D - 1)*(D - 2)**2)), 'beta': (-D/(4*(D - 1)), -(D**3 - 3*D**2 - 4*D + 8)/(4*(D - 1)*(D - 2))),
         'gamma1': (-D*(D - 2)/4, -(D**3 - 3*D**2 - 4*D + 8)/4), 'gamma2': (D/4, -(D**3 - 3*D**2 - 4*D + 8)/(4*(D - 1))),
         'gamma3': (-D/4, 0), 'delta': (0, 0), 'eps1': (0, (D**2 - 6*D + 4)/(2*(D - 1)*(D - 2))), 'eps2': (0, (1 - 2*D)/((D - 1)*(D - 2))),
         'eps3': (0, 0), 'zeta': (0, 0)}
print(f"{'op':7s} {'derived col1 (4/y)^(D-1)':32s} {'derived col2 (4/y)^(D-2)':46s}  paper-match")
for k, f in T.items():
    c1 = sp.factor(sp.simplify(sum(v for q, v in f.items() if sp.simplify(q - (D - 1)) == 0)/16))
    c2 = sp.factor(sp.simplify(sum(v for q, v in f.items() if sp.simplify(q - (D - 2)) == 0)/16))
    ok1 = sp.simplify(c1 - paper[k][0]) == 0; ok2 = sp.simplify(c2 - paper[k][1]) == 0
    print(f"{k:7s} {str(c1):32s} {str(c2):46s}  {ok1} {ok2}")
