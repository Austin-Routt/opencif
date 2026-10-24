"""Every scale bar the traces use, read two ways (full pixels only, as rc5 and both reviewers did; and with the
anti-aliased end pixels counted by their coverage), the Xia printed-width scale, and the Dinh Fig. 3C consistency check.

usage: python3 measure_scales.py RASTERS OUT.json"""
import json
import os
import sys

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trace_lib as T
import panels_xia as PX

RAS, OUTP = sys.argv[1:3]
im = {f: T.Raster(os.path.join(RAS, f)).im for f in ('dinh5-000.png', 'dinh6-000.png', 'gif5-001.png', 'xia5-000.png', 'xia3-000.png')}
bars = {}
for i, r0 in enumerate((470, 952, 1436, 1920)):
    b = T.bar_length_aa(im['dinh5-000.png'], r0, r0 + 20, 870, 975)
    bars['dinh_fig2B_panel%d' % (i + 1)] = dict(um=50, px_aa=b[0], sd_px=b[1], full_px=b[3], rows=b[2])
for lab, (r0, r1) in (('dinh_fig3C_i', (540, 560)), ('dinh_fig3C_ii', (1215, 1240))):
    b = T.bar_length_aa(im['dinh6-000.png'], r0, r1, 1300, 1400)
    bars[lab] = dict(um=100, px_aa=b[0], sd_px=b[1], full_px=b[3], rows=b[2])
b = T.bar_length_aa(im['gif5-001.png'], 550, 600, 690, 829)
bars['gifford_fig4_inset'] = dict(um=300, px_aa=b[0], sd_px=b[1], full_px=b[3], rows=b[2])
b = T.bar_length_aa(im['xia5-000.png'], 3085, 3120, 1860, 2025, dark=False)
bars['xia_fig3_panel_i'] = dict(um=250, px_aa=b[0], sd_px=b[1], full_px=b[3], rows=b[2])
b = T.bar_length_aa(im['xia5-000.png'], 2925, 2960, 1000, 1230)
bars['xia_fig3a'] = dict(um=5000, px_aa=b[0], sd_px=b[1], full_px=b[3], rows=b[2])
for lab, (c0, c1) in (('xia_fig1_e_i', (1030, 1210)), ('xia_fig1_e_ii', (1830, 2010))):
    b = T.bar_length_aa(im['xia3-000.png'], 2900, 2925, c0, c1)
    bars[lab] = dict(um=250 if lab.endswith('_i') else 50, px_aa=b[0], sd_px=b[1], full_px=b[3], rows=b[2])
for v in bars.values():
    v['um_per_px_aa'] = v['um'] / v['px_aa']
    v['um_per_px_full'] = v['um'] / v['full_px']
    v['difference_pct'] = 100.0 * (v['px_aa'] / v['full_px'] - 1.0)

# Dinh Fig. 3C: the post pitch by autocorrelation of the post columns' profiles, against its 100 µm bar
def autocorr_period(p, guess):
    p = p - p.mean()
    ac = np.correlate(p, p, 'full')[len(p) - 1:]
    pk = []
    for m in range(1, 9):
        lo, hi = int(m * guess * 0.9), int(m * guess * 1.1) + 1
        j = lo + int(np.argmax(ac[lo:hi]))
        a, b_, c = ac[j - 1], ac[j], ac[j + 1]
        pk.append((m, j + 0.5 * (a - c) / (a - 2 * b_ + c)))
    return float(np.polyfit([m for m, _ in pk], [q for _, q in pk], 1)[0])
fig3c = {}
for lab, (r0, r1, cols) in (('C_i', (110, 535, (1177, 1248))), ('C_ii', (715, 1215, (1179, 1240)))):
    per = [autocorr_period(T.col_profile(im['dinh6-000.png'], c, r0, r1, 2), 27.0) for c in cols]
    k = bars['dinh_fig3C_i']['px_aa'] / 100.0
    fig3c[lab] = dict(period_px=per, pitch_um_at_bar=[p / k for p in per])
# C(i) at the start of the linear section (caption: "the beginning of the linear section"): widths in px
rows = []
p = ndi.gaussian_filter1d(T.col_profile(im['dinh6-000.png'], 1177, 115, 530, 2), 1.5)
mx = [115 + i for i in range(3, len(p) - 3) if p[i] == p[i - 3:i + 4].max() and p[i] > np.percentile(p, 60)]
W = []
for r in mx:
    L = T.dark_lines(T.row_profile(im['dinh6-000.png'], r, 1120, 1310, 1), 8, off=1120)
    a = [T.line_near(L, x, 3.0) for x in (1142.0, 1157.0, 1187.8, 1233.0, 1263.8, 1278.6)]
    if None not in a:
        W.append(a)
W = np.array(W)
w = dict(walls=float(np.mean(W[:, 5] - W[:, 0])), side_L=float(np.mean(W[:, 1] - W[:, 0])), post_L=float(np.mean(W[:, 2] - W[:, 1])),
         centre=float(np.mean(W[:, 3] - W[:, 2])), post_R=float(np.mean(W[:, 4] - W[:, 3])), side_R=float(np.mean(W[:, 5] - W[:, 4])), rows=len(W))
k = bars['dinh_fig3C_i']['px_aa'] / 100.0
fig3c['C_i_widths_px'] = w
fig3c['C_i_widths_um_at_bar'] = {a: (b / k if a != 'rows' else b) for a, b in w.items()}
kp = float(np.mean(fig3c['C_i']['period_px'])) / 45.1
fig3c['C_i_widths_um_at_pitch_45.1'] = {a: (b / kp if a != 'rows' else b) for a, b in w.items()}
fig3c['note'] = ('at its 100 µm bar Fig. 3C reads a 51.2 µm post pitch and an 86 µm centre at "the beginning of the linear section" '
                 '(Fig. 2 caption: 75 µm there); at the Fig. 2B pitch the centre reads 75 µm: the bar does not calibrate the panel')
out = dict(bars=bars, xia_printed_width_scale=PX.scale_evidence(im['xia5-000.png']), dinh_fig3C=fig3c)
json.dump(out, open(OUTP, 'w'), indent=1, default=float)
for a, v in bars.items():
    print('%-22s %5d µm  full %4d px  aa %8.3f px  (%+.2f %%)  %.4f µm/px' % (a, v['um'], v['full_px'], v['px_aa'], v['difference_pct'], v['um_per_px_aa']))
print('Xia printed-width scale %.4f µm/px' % out['xia_printed_width_scale']['um_per_px_printed'])
print('Dinh Fig. 3C', json.dumps(fig3c, default=float)[:900])
