"""Figure 1 of the rc5 verification report: the stretch-end coefficients c (in units of the narrow width) against the width
ratio, deep (planar Stokes) and shallow (Hele-Shaw) limits, from opencif's tables, with rc4's 21 changes of width.
Usage: python3 fig_coeffs.py OPENCIF_ROOT OUT.png rc4_step_excess.json"""
import sys, os, json, math
ROOT = sys.argv[1]; OUT = sys.argv[2]; SURV = sys.argv[3]
sys.path[:0] = [ROOT, os.path.join(ROOT, 'opencif'), os.path.join(ROOT, 'opencif', 'engine')]
import numpy as np
import cifsched as S
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
R = np.exp(np.linspace(0, math.log(8.0), 300))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
a = ax[0]
a.plot(R, [S._interp_logR(S.STEP_ONE_DEEP, r) for r in R], color='#1f4e79', lw=2, label='one-sided step, deep')
a.plot(R, [S.step_shallow(r, True) for r in R], color='#1f4e79', lw=2, ls='--', label='one-sided step, shallow (closed form)')
a.plot(R, [S._interp_logR(S.STEP_SYM_DEEP, r) for r in R], color='#b35900', lw=2, label='centred step, deep')
a.plot(R, [S.step_shallow(r, False) for r in R], color='#b35900', lw=2, ls='--', label='centred step, shallow (closed form)')
a.plot([x for x, _ in S.STEP_ONE_DEEP], [y for _, y in S.STEP_ONE_DEEP], 'o', color='#1f4e79', ms=4)
a.plot([x for x, _ in S.STEP_SYM_DEEP], [y for _, y in S.STEP_SYM_DEEP], 'o', color='#b35900', ms=4)
Rf = [r for r in R if r <= 4.0]
a.plot(Rf, [S._flare_tab(S.FLARE_DEEP, 14.5, r) for r in Rf], color='#2e7d32', lw=2, label='flare 14.5°, deep')
a.plot(Rf, [S._flare_tab(S.FLARE_SHALLOW, 14.5, r) for r in Rf], color='#2e7d32', lw=2, ls='--', label='flare 14.5°, shallow')
sv = json.load(open(SURV))
rows = sv if isinstance(sv, list) else sv.get('rows', sv)
for d in [x for x in rows if x['kind'] != 'none']:
    a.plot(d['ratio'], d['c'], marker='x', color='k', ms=6, ls='none')
a.plot([], [], marker='x', color='k', ls='none', label="rc4's 21 changes of width (c used, depth-weighted)")
a.set_xscale('log'); a.set_xticks([1, 1.5, 2, 3, 4, 6, 8]); a.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
a.set_xlabel('width ratio  w_wide / w_narrow'); a.set_ylabel('excess c  (lengths of the narrow width)')
a.set_title('Stretch-end excess: 2-D solves in the two limits', fontsize=10); a.grid(alpha=0.3); a.legend(fontsize=7, loc='upper left')
b = ax[1]
degs = np.linspace(0, 45, 100)
for r, col in ((1.5, '#9e9e9e'), (2.4, '#2e7d32'), (4.0, '#1b5e20')):
    b.plot(degs, [S._flare_tab(S.FLARE_DEEP, g, r) for g in degs], color=col, lw=2, label='ratio %.1f, deep' % r)
    b.plot(degs, [S._flare_tab(S.FLARE_SHALLOW, g, r) for g in degs], color=col, lw=2, ls='--', label='ratio %.1f, shallow' % r)
b.axvline(14.5, color='k', lw=0.8, ls=':'); b.text(15, 0.27, 'Dinh 14.5°', fontsize=8)
b.axhline(S._interp_logR(S.STEP_ONE_DEEP, 2.4), color='#1f4e79', lw=1, ls='-.'); b.text(1, S._interp_logR(S.STEP_ONE_DEEP, 2.4) + 0.01, 'right-angle step at ratio 2.4, deep: 0.488', fontsize=8, color='#1f4e79')
b.set_xlabel('flare angle per side (deg)'); b.set_ylabel('excess c (lengths of the narrow width)'); b.set_ylim(0, 0.55)
b.set_title('Flare excess over its integrated widening', fontsize=10); b.grid(alpha=0.3); b.legend(fontsize=7, loc='upper left', bbox_to_anchor=(0.0, 0.86))
fig.tight_layout(); fig.savefig(OUT, dpi=150); print('wrote', OUT)
