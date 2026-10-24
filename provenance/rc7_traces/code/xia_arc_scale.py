"""Xia Fig. 3d: the scale from the loops' curvature.  Each arc wall of the two loop ends (a dark line followed down
232 rows) is fitted with a circle; each loop's arcs share a centre (concentric fit); their radii in px against the
mask's arc radii in µm (concentric circles of the shipped DXF about the loop centres, set by Fig. 3a's leg axes at its
own 5 mm bar) give µm/px.  Curvature does not depend on the outline bias (a 1 px bias moves a 1300 px radius by
0.08 %) or on the printed widths.  usage: xia_arc_scale.py RASTER DXF OUT_STEM"""
import json
import os
import sys

import numpy as np
from scipy import optimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trace_lib as T

ras, dxf, stem = sys.argv[1:4]
im = T.Raster(ras, 'e9fb7fbb348567d3').im
c0, r0, c1, r1 = 1, 800, 747, 1538
rows = list(range(r0 + 6, r0 + 470, 2))
pts = {r: T.dark_lines(T.row_profile(im, r, c0, c1, 1), 8, halfwin=5, off=c0) for r in rows}
curves = []
for x in pts[rows[0]]:
    cur = [(x, rows[0])]
    for r in rows[1:]:
        v = T.line_near(pts[r], cur[-1][0], 3.0)
        if v is None:
            break
        cur.append((v, r))
    if len(cur) > 150:
        curves.append(np.array(cur))


def concentric(Ps, guess):
    def f(p):
        return np.concatenate([np.hypot(P[:, 0] - p[0], P[:, 1] - p[1]) - p[2 + i] for i, P in enumerate(Ps)])
    x0 = list(guess) + [float(np.hypot(P[0, 0] - guess[0], P[0, 1] - guess[1])) for P in Ps]
    s = optimize.least_squares(f, x0)
    J = s.jac
    cov = np.linalg.inv(J.T @ J) * np.var(s.fun)
    return s.x, np.sqrt(np.diag(cov)), float(np.std(s.fun))


left = [P for P in curves if P[-1, 0] < 300]
right = [P for P in curves if P[-1, 0] > 450]
pl, el, sl = concentric(left, (1490.0, 1350.0))
pr, er, sr = concentric(right, (1490.0, 1350.0))
# the mask's arcs at the left end: concentric circles about the two loop centres
mask = T.Mask(dxf)


def mask_arcs():
    out = []
    for (L, p, x0, x1, y0, y1) in mask.polys:
        if x0 > -4100:
            continue
        q = np.vstack([p, p[:1]])
        idx = np.where(q[:, 0] < -4150)[0]
        if len(idx) < 20:
            continue
        for run in np.split(idx, np.where(np.diff(idx) != 1)[0] + 1):
            P = q[run]
            if len(P) < 20:
                continue
            x, y = P[:, 0], P[:, 1]
            A = np.c_[x, y, np.ones_like(x)]
            c = np.linalg.lstsq(A, x ** 2 + y ** 2, rcond=None)[0]
            xc, yc = c[0] / 2, c[1] / 2
            R = np.sqrt(c[2] + xc ** 2 + yc ** 2)
            if np.max(np.abs(np.hypot(x - xc, y - yc) - R)) < 0.5 and 300 < R < 4000:
                out.append((float(R), float(xc), float(yc)))
    return sorted(set((round(a, 2), round(b, 2), round(c, 2)) for a, b, c in out))


MA = mask_arcs()
# the loops at the device's top (x < -4100): loop 2 (legs 2 -> 3) and loop 4 (legs 4 -> 5), centred half way between
# their legs' axes (the DXF's record: design.leg_Y; rc5's hard-coded -206.0 / -209.5 belonged to rc5's axes)
_Y = json.load(open(dxf[:-4] + '.json'))['design']['leg_Y']
yc2, yc4 = 0.5 * (_Y[1] + _Y[2]), 0.5 * (_Y[3] + _Y[4])
loop4 = [a[0] for a in MA if abs(a[2] - yc4) < 1.0]
loop2 = [a[0] for a in MA if abs(a[2] - yc2) < 1.0]
# pairing by the pattern of radial gaps (islands 78.1 µm read about 42 px): the left set's five arcs are loop 4's
# second to sixth... found by least squares over every consecutive assignment
def best_assign(Rpx, Rum):
    best = None
    for j in range(len(Rum) - len(Rpx) + 1):
        s = np.array(Rum[j:j + len(Rpx)]) / np.array(Rpx)
        sc = float(np.std(s))
        if best is None or sc < best[0]:
            best = (sc, j, s)
    return best
Lpx = sorted(pl[2:])
Rpx = sorted(pr[2:])
# joint assignment: every consecutive placement of each set's arcs among its loop's mask arcs; the pair of placements
# whose seven ratios agree best (one scale for the panel)
L4, L2 = sorted(loop4), sorted(loop2)
best = None
for jl in range(len(L4) - len(Lpx) + 1):
    for jr in range(len(L2) - len(Rpx) + 1):
        s = np.r_[np.array(L4[jl:jl + len(Lpx)]) / np.array(Lpx), np.array(L2[jr:jr + len(Rpx)]) / np.array(Rpx)]
        if best is None or np.std(s) < best[0]:
            best = (float(np.std(s)), jl, jr, s)
bl = (None, best[1], best[3][:len(Lpx)])
br = (None, best[2], best[3][len(Lpx):])
res = dict(curves=len(curves),
           left_set=dict(centre_px=pl[:2].tolist(), radii_px=Lpx, sd_px=sl, mask_radii_um=sorted(loop4)[bl[1]:bl[1] + len(Lpx)], scale_um_per_px=bl[2].tolist()),
           right_set=dict(centre_px=pr[:2].tolist(), radii_px=Rpx, sd_px=sr, mask_radii_um=sorted(loop2)[br[1]:br[1] + len(Rpx)], scale_um_per_px=br[2].tolist()),
           mask_loop2_radii_um=sorted(loop2), mask_loop4_radii_um=sorted(loop4),
           mask=os.path.basename(dxf), mask_sha256=mask.sha256, loop_centres_y_um=[yc2, yc4])
allS = np.r_[bl[2], br[2]]
res['scale_um_per_px'] = dict(mean=float(allS.mean()), sd=float(allS.std()), n=int(len(allS)),
                              fig3a_leg_axes_uncertainty_pct=100 * 20.0 / 2430.5)
json.dump(res, open(stem + '.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
# figure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(3.4, 3.6), dpi=600)
ax = fig.add_axes([0, 0.1, 1, 0.9])
ax.imshow(im[r0:r1, c0:c1], cmap='gray', vmin=0, vmax=255, interpolation='nearest', extent=(c0 - .5, c1 - .5, r1 - .5, r0 - .5))
t = np.linspace(0, 2 * np.pi, 4000)
for (p, Rs) in ((pl, pl[2:]), (pr, pr[2:])):
    for R in Rs:
        ax.plot(p[0] + R * np.cos(t), p[1] + R * np.sin(t), '--', color=T.TRACE_COLOUR, lw=0.35, dashes=(3, 2))
for P in curves:
    ax.plot(P[:, 0], P[:, 1], '.', color='#ffb000', ms=0.4)
ax.set_xlim(c0 - .5, c1 - .5)
ax.set_ylim(r1 - .5, r0 - .5)
ax.axis('off')
fig.text(0.02, 0.01, 'Xia Fig. 3d: arc walls (amber, dark-line centres) and their concentric circle fits (cyan dashed). Radii against the\n'
         'mask\'s loop arcs (Fig. 3a\'s leg axes): %.3f ± %.3f µm/px (n %d); the panel-i bar: 1.936; the printed widths: 1.803.' % (
             allS.mean(), allS.std(), len(allS)), fontsize=3.6, va='bottom')
fig.savefig(stem + '.png')
fig.savefig(stem + '.pdf')
