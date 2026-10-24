"""Dinh 2024 V1: the least number of gap rows, and so of device length, that the 1-D model needs.
Premises: the frozen cutoff law D = kappa_g x(f) w_c with the cubic kernel; constant D (Dinh's f* = 1.04e-4 at the
75 um inlet width gives D); every centre width >= 52 um, so f(i) <= f_max = f at 52 um; the march's own recursion
P(i+1) = P(i) (1 - 2 f(i)) (cifsched.plan: two gaps per row); every gap row at least the measured post pitch long.
  N >= ln(P_end / P_start) / ln(1 - 2 f_max),  P_end = 1 - F
P_start = 1 counts every row, head and linear section alike, PROVIDED no row extracts more than f_max: true of a head
at constant D (its rows extract f* at 75 um), not of a head drawn with larger fractions.  The bound is also computed
with the preset's own head, crediting every row whose f exceeds f_max (Dinh V1's head as depicted: its first 18 rows,
up to 6.4 f_max; as specified: gap 1 at 10.1 f_max).  P_start = P after the reconstructed head counts the linear
section only.

rc6 port of provenance/rc5_scripts/dinh_length_bound.py: (1) the pitch is the preset's own row pitch G + g_p (rc6:
19 + 26.1 = 45.1 um, the 50 um bars read with their anti-aliased ends; rc5 hard-coded 46.0), its tolerance an argument
(rc6 0.2 um: the six column pitches' SD 0.044 um and the bar's reading, measure_rc6.py; rc5 1.0); (2) the head's rows
come from cifsched.plan of the preset itself, not from a run pickle, so the bound needs no build; (3) the length key
is length_mm_at_pitch.  With --pitch-tol 1.0 on the rc5 tree it reproduces rc5's dinh_length_bound_*.json apart from
that key (provenance/rc6_scripts/out/length_regression_rc5.out).
Usage: python3 dinh_length_bound.py OPENCIF_ROOT PRESET.yaml OUT.json [--pitch-tol UM]"""
import json
import math
import os
import sys

argv, tol = sys.argv[1:], 0.2
if '--pitch-tol' in argv:
    i_ = argv.index('--pitch-tol')
    tol = float(argv[i_ + 1])
    del argv[i_:i_ + 2]
R, PRESET, OUT = argv[:3]
sys.path[:0] = [R]
from opencif import spec as SP
from opencif.physics import schedule as S
import numpy as np

sp, _ = SP.to_engine(SP.load(PRESET))
pitch, dpitch = float(sp.G + sp.g_p), tol
fstar, w_in, w_min = 1.04e-4, 75.0, 52.0
x_star = S.x_of(fstar)                    # cubic kernel inverse
D = S.KAPPA_G * x_star * w_in
x_max = x_star * w_in / w_min
f_max = 3 * x_max ** 2 - 2 * x_max ** 3
per_row = math.log(1.0 - 2.0 * f_max)
out = dict(preset=os.path.basename(PRESET), pitch_um=pitch, pitch_tolerance_um=dpitch, D_um=D, kappa_g=S.KAPPA_G,
           x_star=x_star, f_max=f_max, ln_one_minus_2fmax=per_row)
F_identity = 0.5592819159228324           # the outlet identity at the printed 44 / 52 um (the presets' F)
F_measured = 0.505                        # Dinh Fig. 4, mean +- SD 0.033
sc = S.plan(sp)
wc, P, NM, fa = np.asarray(sc.wc, float), np.asarray(sc.P, float), int(sc.NM), np.asarray(sc.f, float)
P_head = float(P[NM - 1])                                         # P after the head's last row
i52 = int(np.argmax(wc <= w_min + 1e-9))
F_taper = 1.0 - float(P[i52])                                     # constant-shear taper reaches 52 um here
out.update(P_after_reconstructed_head=P_head, row_where_wc_reaches_52=i52 + 1, F_at_taper_end_width=F_taper,
           P_over_Pstart_at_taper_end=float(P[i52]) / P_head)
over = np.nonzero(fa > f_max)[0]
head_credit = float(sum(math.log(1.0 - 2.0 * fa[i]) / per_row - 1.0 for i in over))
out.update(rows_above_f_max=[int(i) + 1 for i in over], largest_f_over_f_max=float(fa.max() / f_max),
           rows_credited_to_this_head=head_credit)
rows = []
for Flab, F in (('outlet identity (printed 44/52 um)', F_identity),
                ('constant-shear taper reaching 52 um (assumed taper)', F_taper),
                ('measured recovery (Dinh Fig. 4)', F_measured)):
    for Plab, P0, credit in (('every row, a head at constant D (no row above f_max)', 1.0, 0.0),
                             ('every row, this design\'s own head', 1.0, head_credit),
                             ('linear rows after the reconstructed head', P_head, 0.0)):
        n = math.log((1.0 - F) / P0) / per_row - credit
        N = math.ceil(n)
        rows.append(dict(F_reading=Flab, F=round(F, 4), counted=Plab, P_start=round(P0, 5), N_continuous=round(n, 2), N_min=N,
                         length_mm_at_pitch=round(N * pitch / 1000, 2),
                         length_mm_range=[round(N * (pitch - dpitch) / 1000, 2), round(N * (pitch + dpitch) / 1000, 2)],
                         pitch_for_75p8_mm_um=round(75800.0 / N, 2)))
out['bounds'] = rows
json.dump(out, open(OUT, 'w'), indent=1)
print('%s: pitch %.2f +- %.2f um; D %.6f um, x* %.6g, f_max %.6g' % (out['preset'], pitch, dpitch, D, x_star, f_max))
for k in ('P_after_reconstructed_head', 'row_where_wc_reaches_52', 'F_at_taper_end_width', 'P_over_Pstart_at_taper_end',
          'rows_credited_to_this_head'):
    print(k, out[k])
for r in rows:
    print('%-52s F %.4f | %-52s P0 %.5f | N >= %.2f -> %d rows | %.2f mm at %.1f um (%.2f-%.2f) | 75.8 mm needs pitch <= %.2f um'
          % (r['F_reading'], r['F'], r['counted'], r['P_start'], r['N_continuous'], r['N_min'], r['length_mm_at_pitch'], pitch,
             r['length_mm_range'][0], r['length_mm_range'][1], r['pitch_for_75p8_mm_um']))
