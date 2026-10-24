"""Dinh two-cutoff device: the depicted head (as the figure-aligned preset draws it) and a linear section at cutoff D:
linear rows, linear length, device length, for D from 2.2 to 4.8 um (the 1-D outlet identity's F; the record's own
taper map w_c(P); frozen cubic and kappa_g).  usage: python length_vs_D.py TREE"""
import pickle, sys, json
import numpy as np
T = sys.argv[1]
sys.path[:0] = [T, T + '/opencif', T + '/opencif/engine']
a = pickle.load(open(T + '/runs/dinh2024_v1_figure_aligned/dinh2024_v1_figure_aligned.pkl', 'rb'))['arrays']
r = json.load(open(T + '/runs/dinh2024_v1_figure_aligned/dinh2024_v1_figure_aligned.json'))
k = np.asarray(a['kind']); L = np.where(k == 'L')[0]
P = np.asarray(a['P'], float); wc = np.asarray(a['wc'], float)
Pin, Pend = P[L[0] - 1], P[L[-1]]
o = np.argsort(P[L]); Pk, wk = P[L][o], wc[L][o]
KG, PITCH = 4.9138, 45.1
fp = r['design'].get('footprint_mm') or None
def rows(D):
    p, n = Pin, 0
    while p > Pend and n < 20000:
        x = (D / KG) / float(np.interp(p, Pk, wk)); f = 3 * x * x - 2 * x ** 3
        p *= (1 - 2 * f); n += 1
    return n
head_mm, other_mm = 25.566, 65.204 - 53.112      # the record: array head 25.6 mm; footprint minus array
print('rows check at 4.29 um: %d (record 638 L rows)' % rows(4.29))
for D in (2.17, 2.5, 3.0, 3.3, 3.5, 3.6, 3.7, 4.0, 4.29, 4.5, 4.79):
    n = rows(D); lin = n * PITCH / 1e3
    print('D %.2f um: %4d linear rows, %5.1f mm linear; device about %5.1f mm (head %.1f + other %.1f); printed 75.8 mm' % (
        D, n, lin, lin + head_mm + other_mm, head_mm, other_mm))
