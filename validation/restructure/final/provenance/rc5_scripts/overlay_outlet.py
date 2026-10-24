"""Overlay the drawn Dinh V1 outlet end (rc4: step; rc5: flare) on Dinh 2024 Fig. 2B, fourth panel (native raster,
1.087 um/px from the 50 um bar).  Registration: the last gap's centre row (mean of the two pill columns' measured last
gap centres) and the channel axis (mid-point of the two outer walls' straight fits); the walls' own tilt (0.35 deg left,
0.70 deg right) is not removed, so the overlay is good to about 1 px over the 300 um shown.  No fitting to the flare.
Usage: OPENCIF_ROOT=<tree> python3 -I overlay_outlet.py RASTER FLARE_JSON PITCH_JSON OUT_PNG LABEL=PKL [LABEL=PKL ...]"""
import json
import pickle
import sys

import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

raster, fj, pj, out = sys.argv[1:5]
runs = [a.split('=', 1) for a in sys.argv[5:]]
import os
ROOT = os.environ.get('OPENCIF_ROOT', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
for p in (ROOT, os.path.join(ROOT, 'opencif'), os.path.join(ROOT, 'opencif', 'engine')):
    sys.path.append(p)               # pickles reference opencif's classes; appended so stdlib wins (OPENCIF_ROOT)

F = json.load(open(fj))
P = json.load(open(pj))
um = F['scale_bar']['um_per_px']
row_gap = 0.5 * (P['left_column']['last_gap_centre_row'] + P['right_column']['last_gap_centre_row'])
lp = np.array(F['left_picks'], float).reshape(-1, 2)
rp = np.array(F['right_picks'], float).reshape(-1, 2)
sel = lambda a: a[(a[:, 0] >= 1590) & (a[:, 0] <= 1750)][:, 1].mean()
col_axis = 0.5 * (sel(lp) + sel(rp))
img = np.asarray(Image.open(raster).convert('L'), float)
r0, r1, c0, c1 = 1480, 1935, 560, 980

fig, axes = plt.subplots(1, len(runs), figsize=(6.2 * len(runs), 7.2), squeeze=False)
colors = ['#d62728', '#1f77b4']
for k, (label, pk) in enumerate(runs):
    ax = axes[0][k]
    ax.imshow(img[r0:r1, c0:c1], cmap='gray', extent=(c0, c1, r1, r0))
    d = pickle.load(open(pk, 'rb'))
    a, fluid = d['arrays'], d['fluid']
    N = int(a['N'])
    leg = a['legs'][-1]
    gx0 = float(a['gx'][-1])
    x_last = gx0 + (float(a['x'][N - 1]) - float(leg['p0'] or 0.0)) * leg['dirn']
    Y = float(leg['Y'])
    from shapely.geometry import box
    win = box(x_last - (row_gap - r0) * um - 5, Y - 260, x_last + (r1 - row_gap) * um + 5, Y + 260)
    clip = fluid.intersection(win)
    geoms = getattr(clip, 'geoms', [clip])
    for g in geoms:
        for ring in [g.exterior] + list(g.interiors):
            xy = np.asarray(ring.coords)
            rows = row_gap + (xy[:, 0] - x_last) * leg['dirn'] / um
            cols = col_axis + (xy[:, 1] - Y) / um
            ax.plot(cols, rows, '-', color=colors[k % 2], lw=0.9)
    ax.set_xlim(c0, c1)
    ax.set_ylim(r1, r0)
    ax.axhline(row_gap, color='#ff7f0e', lw=0.6, ls=':')
    ax.set_title('%s drawn outline on Dinh 2024 Fig. 2B (4th panel)' % label, fontsize=9)
    ax.set_xlabel('raster column (1 px = %.3f um)' % um, fontsize=8)
    ax.set_ylabel('raster row (flow downward)', fontsize=8)
    ax.tick_params(labelsize=7)
fig.text(0.5, 0.01, 'Registered on the last gap centre (dotted) and the channel axis only; no fit to the flare.',
         ha='center', fontsize=8)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(out, dpi=150)
print('wrote', out, 'row_gap %.2f col_axis %.2f' % (row_gap, col_axis))
