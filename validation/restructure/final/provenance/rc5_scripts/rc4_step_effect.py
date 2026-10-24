"""rc4's Dinh V1 (depicted): what the computed right-angle step excess (rc5's c_exp, 0.513 of the 44 um side
channel, i.e. 3.49 % of R_s,out) does to the last rows when it is added to the side outlets and NOT to the closure
(rc4's situation), solved on the design network with the drawn outlet resistances of the rc4 record (as rc4's
step_sens.py).  Nonlinear: the full ladder is re-solved.
usage: python3 rc4_step_effect.py RC4_ROOT PRESET.yaml RECORD.json C_EXP"""
import sys, json
import numpy as np
root, preset, recp, cexp = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
sys.path.insert(0, root)
import opencif.pipeline as PL, opencif.spec as SP
S = sys.modules['cifsched']; M2 = sys.modules.get('measure2') or __import__('measure2')
spec, opts, sc, g, env, notes = PL.prepare(SP.load(preset))
N = sc.N
Rc, Rs, Rg = np.asarray(sc.Rc, float), np.asarray(sc.Rs, float), np.asarray(sc.Rg, float)
wc = np.asarray(sc.wc, float)
D = S.KAPPA_G * spec.delta_star
om = json.load(open(recp))['validation']['outlet_measured']
Rco, Rso = om['Rc_out'], 0.5 * (om['RsA_out'] + om['RsB_out'])
def cut(Rso_):
    r = M2.solve3(Rc[:-1], Rs, Rs, Rg, Rg, Rco, Rso_, Rso_)
    f = np.asarray(r['fA'], float)
    dd = np.array([S.KAPPA_G * S.x_of(fe) * w if 0 < fe < 0.5 else np.nan for fe, w in zip(f, wc)])
    return dd / D, f, float(r['F'])
c0, f0, F0 = cut(Rso)
wN = float(sc.ws[N - 2])
dR = S.rr(wN, spec.h) * cexp * wN * S.UM
c1, f1, F1 = cut(Rso + dR)
ok0 = np.isfinite(c0) & (np.abs(c0 - 1) <= 0.01)
out1 = ~(np.isfinite(c1) & (np.abs(c1 - 1) <= 0.01))
newly = ok0 & out1
lost = ok0 & ~np.isfinite(c1)
rows = np.where(newly)[0] + 1
mv = np.where(ok0 & np.isfinite(c1), np.abs(c1 - c0), -1)
res = dict(c_exp=cexp, w_s_N_um=wN, dR_over_Rs_out=dR / Rso, rows_newly_outside_1pct=int(newly.sum()),
           first_last=[int(rows[0]), int(rows[-1])] if len(rows) else None, rows_losing_cutoff=int(lost.sum()),
           lost_rows=[int(r) + 1 for r in np.where(lost)[0]], max_change_D_finite=float(mv.max()),
           last_rows_f_over_fstar_before=[round(float(x), 4) for x in (f0[-6:] / f0[N // 2])],
           last_rows_f_over_fstar_after=[round(float(x), 4) for x in (f1[-6:] / f0[N // 2])],
           F_before=F0, F_after=F1, N=N)
print(json.dumps(res, indent=1))
