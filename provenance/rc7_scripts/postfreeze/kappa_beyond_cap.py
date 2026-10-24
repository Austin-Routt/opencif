"""Dinh figure-aligned: how the row count -> cutoff inference moves if the cutoff law's prefactor kappa differs beyond
Gifford's calibration cap (f > 5.76e-4) by a factor (1 + e).  Uses the shipped record's own linear section: the centre
width as a function of the centre flow P (the constant-shear taper maps P to w_c whatever D is), the frozen cubic
f = 3x^2 - 2x^3 with x = delta/w_c, D = kappa delta.  Marches from the first linear row's P to the last's and finds the
D that needs as many rows as the record draws.  usage: python kappa_beyond_cap.py TREE"""
import pickle, sys, math
import numpy as np
T = sys.argv[1]
sys.path[:0] = [T, T + '/opencif', T + '/opencif/engine']
a = pickle.load(open(T + '/runs/dinh2024_v1_figure_aligned/dinh2024_v1_figure_aligned.pkl', 'rb'))['arrays']
k = np.asarray(a['kind']); L = np.where(k == 'L')[0]
P = np.asarray(a['P'], float); wc = np.asarray(a['wc'], float); f = np.asarray(a['f'], float)
KG, CAP = 4.9138, 5.76e-4
Pin = P[L[0] - 1] if L[0] > 0 else 1.0          # centre flow entering the first linear row
Pend = P[L[-1]]
order = np.argsort(P[L]); Pk, wk = P[L][order], wc[L][order]
wc_of = lambda p: float(np.interp(p, Pk, wk))
def g_inv(fv):                                    # x from f = 3x^2 - 2x^3 (x < 1/2)
    lo, hi = 0.0, 0.5
    for _ in range(60):
        m = 0.5 * (lo + hi)
        (lo, hi) = (m, hi) if 3 * m * m - 2 * m ** 3 < fv else (lo, m)
    return 0.5 * (lo + hi)
def rows(D, e):
    p, n = Pin, 0
    while p > Pend and n < 5000:
        w = wc_of(p)
        x = (D / KG) / w; fv = 3 * x * x - 2 * x ** 3
        if fv > CAP:                                # beyond the cap the prefactor is KG (1 + e): delta = D / (KG (1 + e))
            x = (D / (KG * (1 + e))) / w; fv = max(3 * x * x - 2 * x ** 3, CAP) if 3 * x * x - 2 * x ** 3 < CAP else 3 * x * x - 2 * x ** 3
        p *= (1 - 2 * fv); n += 1
    return n
def D_for(nrows, e):
    lo, hi = 2.0, 8.0
    for _ in range(40):
        m = 0.5 * (lo + hi)
        (lo, hi) = (lo, m) if rows(m, e) < nrows else (m, hi)
    return 0.5 * (lo + hi)
import json
D_rec = float(json.load(open(T + '/runs/dinh2024_v1_figure_aligned/dinh2024_v1_figure_aligned.json'))['design']['D_star_um'])
n_rec = len(L)
print('record: %d linear rows, D %.3f um; P %.4f -> %.4f; w_c %.1f -> %.1f um; f %.2e -> %.2e' % (
    n_rec, D_rec, Pin, Pend, wc[L[0]], wc[L[-1]], f[L[0]], f[L[-1]]))
print('check, e = 0: rows(%.2f) = %d; D for %d rows = %.3f um' % (D_rec, rows(D_rec, 0.0), n_rec, D_for(n_rec, 0.0)))
for e in (-0.10, -0.05, -0.034, 0.0, 0.034, 0.05, 0.10):
    print('kappa beyond the cap x %.3f: D for %d rows = %.3f um' % (1 + e, n_rec, D_for(n_rec, e)))
