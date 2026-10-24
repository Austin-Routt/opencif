"""Sensitivity of the last rows' cutoffs to an extra side-outlet resistance (an unmodelled sudden expansion), on the
design network of a preset, with the drawn outlet resistances of its record.
usage: step_sens.py ROOT PRESET.yaml RECORD.json   (e.g. runs/dinh2024_v1_depicted/dinh2024_v1_depicted.json)"""
import sys, math
import numpy as np
root, preset = sys.argv[1], sys.argv[2]
sys.path.insert(0, root)
import opencif.pipeline as PL, opencif.spec as SP
S = sys.modules['cifsched']; M2 = sys.modules.get('measure2') or __import__('measure2')
d = SP.load(preset)
spec, opts, sc, g, env, notes = PL.prepare(d)
N = sc.N
print('attrs', [a for a in ('Rc','Rs','Rg','Rs_out','Rc_out','Qg') if hasattr(sc, a)])
Rc, Rs, Rg = np.asarray(sc.Rc, float), np.asarray(sc.Rs, float), np.asarray(sc.Rg, float)
Rso = float(sc.Rs_out)
P, Sd, Qg = np.asarray(sc.P, float), np.asarray(sc.S, float), np.asarray(sc.Qg, float)
Rco = (Rso * Sd[N - 1] + Rg[N - 1] * Qg[N - 1]) / P[N - 1]
wc = np.asarray(sc.wc, float)
D = S.KAPPA_G * spec.delta_star
def cut(Rso_):
    r = M2.solve3(Rc[:-1], Rs, Rs, Rg, Rg, Rco, Rso_, Rso_)
    dd = np.array([S.KAPPA_G * S.x_of(fe) * w if 0 < fe < 0.5 else np.nan for fe, w in zip(r['fA'], wc)])
    return dd / D, r['F']
c0, F0 = cut(Rso)
print('base: last rows D/D*', np.round(c0[-5:], 5), 'F', F0, 'rows |c-1|>1e-3:', int(np.sum(np.abs(c0 - 1) > 1e-3)))
h = spec.h
wN = float(sc.ws[N - 2])
for cfac in (0.1, 0.25, 0.5, 1.0):
    dR = S.rr(wN, h) * cfac * wN * S.UM            # an excess of cfac * w of the narrow channel
    c1, F1 = cut(Rso + dR)
    dc = c1 - c0
    print('excess %.2f w (%.1f um): dR/Rs_out %.3e; max |dD/D| %.2e at row %d (of %d); rows moved > 1e-3: %d, > 1e-2: %d; dF %.2e'
          % (cfac, cfac * wN, dR / Rso, np.nanmax(np.abs(dc)), int(np.nanargmax(np.abs(dc))) + 1, N,
             int(np.sum(np.abs(dc) > 1e-3)), int(np.sum(np.abs(dc) > 1e-2)), F1 - F0))
# --- with the drawn outlet resistances from the record
import json, os
rec = json.load(open(sys.argv[3]))
om = rec['validation']['outlet_measured']
Rco_d, Rso_d = om['Rc_out'], om['RsA_out']
def cut2(Rso_, Rco_):
    r = M2.solve3(Rc[:-1], Rs, Rs, Rg, Rg, Rco_, Rso_, Rso_)
    dd = np.array([S.KAPPA_G * S.x_of(fe) * w if 0 < fe < 0.5 else np.nan for fe, w in zip(r['fA'], wc)])
    return dd / D, r['F']
c0, F0 = cut2(Rso_d, Rco_d)
print('drawn outlet: Rso %.4e Rco %.4e; base last rows' % (Rso_d, Rco_d), np.round(c0[-6:], 5), 'F %.6f' % F0,
      'rows |c-1|>1e-3:', int(np.sum(np.abs(c0 - 1) > 1e-3)))
for cfac in (0.05, 0.1, 0.25, 0.5):
    dR = S.rr(wN, h) * cfac * wN * S.UM
    c1, F1 = cut2(Rso_d + dR, Rco_d)
    dc = c1 - c0
    big = np.where(np.abs(dc) > 1e-2)[0]
    print('excess %.2f w (%.1f um): dR/Rs_out %.3e; max |dD/D| %.3e at row %d; rows moved >1e-3: %d, >1e-2: %d (rows %s); dF %.2e'
          % (cfac, cfac * wN, dR / Rso_d, np.nanmax(np.abs(dc)), int(np.nanargmax(np.abs(dc))) + 1,
             int(np.sum(np.abs(dc) > 1e-3)), len(big), (big[0] + 1, big[-1] + 1) if len(big) else '-', F1 - F0))
