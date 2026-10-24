"""xia_fig3a_device.py -- the whole Xia RS-7 mask on the paper's whole-device drawing, Fig. 3a (rc6).

Fig. 3a (Xia et al. 2016, xia5-000.png columns 765-1203) is a line drawing of the whole device at a 5 mm bar; its
post columns appear as dotted lines, so it is a rendering of the device's own layout.  This script registers a
generated mask on it and asks two questions of the registration:

  V11  where does each part of the device sit on the drawing?  The mask's polygon edges are sampled every drawing pixel
       and matched to the drawing's ink by a truncated chamfer distance.  Scale: Fig. 3a's bar read with its
       anti-aliased ends (178.32 px per 5 mm).  Pose: the rotation and translation that fit the inlet port, the head
       and the first loops (the device's top, where the drawing and every schedule agree on the layout), then, with
       that pose held, a translation fitted to each further group alone: each leg's straight part, the loops at the
       far end, the outlet ring and its ports.  A group's shift along the legs is the drawing's length minus the
       mask's up to that group (the rows the schedule places against the rows the drawing shows); across, the leg
       axes (which every preset takes from Fig. 3a) must agree to a drawing pixel.  Also reported: one pose fitted to
       the whole device, and with the scale free (a check on the bar).
  V12  do the micrographs sit where Fig. 3a's boxes put them?  Each box b-i (its centre in xia5-000 px) is mapped
       into the mask by the top-anchored pose; for each traced panel (the trace's own map from mask µm to panel px)
       the panel's centre is mapped into the mask; the two must be within the box's half-size.

usage: python3 -I xia_fig3a_device.py XIA5_000_PNG DXF RECORD OUT_JSON OUT_PNG [TRACE_RESULTS.json:KEY:PANEL ...]
       (KEY '-' for a results file holding one panel; PANEL a letter b-i)"""
import json
import math
import os
import sys

import numpy as np
from scipy import ndimage as ndi
from scipy import optimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trace_lib as T

RAS, DXF, REC, OUT, PNG = sys.argv[1:6]
TRACED = sys.argv[6:]
FIG3_SHA = 'e9fb7fbb348567d3'
A_COLS = (765, 1203)                     # panel a
A_ROWS = (40, 2920)                      # below the letter, above the bar
BOXES = dict(b=(1072, 430, 21), c=(1063, 851, 20.5), d=(988, 277, 20.5), e=(972, 1447, 20.5), f=(921, 967, 20.5),
             g=(908, 2454, 21.5), h=(918, 2658, 21.5), i=(972, 2872, 21.5))
PANEL_CROPS = dict(b=(1, 1, 748, 748), c=(1277, 1, 2025, 748), d=(1, 796, 748, 1542), e=(1277, 797, 2025, 1543),
                   f=(1, 1594, 748, 2340), g=(1277, 1594, 2025, 2340), h=(1, 2390, 748, 3137), i=(1277, 2390, 2025, 3137))
TRUNC = 3.0                               # px: a drawing line is 2-3 px wide

r = T.Raster(RAS, FIG3_SHA)
im = r.im.astype(float)
bar = T.bar_length_aa(r.im, 2925, 2960, 1000, 1230)
k = bar[0] / 5000.0                       # px per µm
crop = (A_COLS[0], A_ROWS[0], A_COLS[1], A_ROWS[1])
sub = im[A_ROWS[0]:A_ROWS[1] + 1, A_COLS[0]:A_COLS[1] + 1]
ink = sub < 150.0
# the lines' centres: ink pixels at a local maximum of their distance to the paper (a line 2-3 px wide has a 1-px
# ridge); the chamfer distance is measured to these, so that a point on a line's edge scores about 1 px, not 0
d_in = ndi.distance_transform_edt(ink)
ridge = ink & (d_in >= ndi.maximum_filter(d_in, size=3) - 1e-9)
dt = ndi.distance_transform_edt(~ridge)
rec = json.load(open(REC))
mask = T.Mask(DXF)
assert mask.sha256 == rec['files']['dxf_sha256'], 'the DXF is not the record\'s'
legY = [float(v) for v in rec['design']['leg_Y']]

# ---- the mask's edges, sampled every drawing pixel, each point tagged with its group
step_um = 1.0 / k
P, G = [], []
posts = mask.posts(-1e9, 1e9, -1e9, 1e9, max_len=60.0, max_wid=100.0)
leg_x = []
for Y in legY:
    q = posts[np.abs(posts[:, 1] - Y) < 400.0]
    leg_x.append((float(q[:, 0].min()), float(q[:, 0].max())) if len(q) else (np.nan, np.nan))
x_head_end = leg_x[0][0]
x_top = min(v[0] for v in leg_x[1:])
x_bot = max(v[1] for v in leg_x)


def group_of(x, y):
    j = int(np.argmin([abs(y - Y) for Y in legY]))
    if x < x_head_end and abs(y - legY[0]) < 1500.0:
        return 'inlet and head'
    if x < x_top:
        return 'loops at the top'
    if x > x_bot + 200.0:
        return 'outlet ring and ports' if x > x_bot + 2500.0 else 'loops at the bottom'
    return 'leg %d' % (j + 1)


for (L, p, x0, x1, y0, y1) in mask.polys:
    q = np.vstack([p, p[:1]])
    d = np.diff(q, axis=0)
    Ls = np.hypot(d[:, 0], d[:, 1])
    cum = np.r_[0, np.cumsum(Ls)]
    if cum[-1] < step_um:
        pts = q[:1]
    else:
        t = np.arange(0.5 * step_um, cum[-1], step_um)
        i = np.clip(np.searchsorted(cum, t) - 1, 0, len(Ls) - 1)
        f = (t - cum[i]) / np.maximum(Ls[i], 1e-12)
        pts = q[:-1][i] + f[:, None] * d[i]
    P.append(pts)
    G += [group_of(x, y) for x, y in pts]
P = np.vstack(P)
G = np.array(G)
groups = sorted(set(G), key=lambda g: (0 if g == 'inlet and head' else 1 if g == 'loops at the top' else
                                       2 if g.startswith('leg') else 3 if g == 'loops at the bottom' else 4, g))


def to_px(pts, prm, kk=None):
    """mask µm -> image px: x (along the legs) down the drawing, y (across) to the right, rotated by theta"""
    c, rr, th = prm[:3]
    kk = k if kk is None else kk
    s, co = math.sin(th), math.cos(th)
    col = c + kk * (pts[:, 1] * co + pts[:, 0] * s)
    row = rr + kk * (pts[:, 0] * co - pts[:, 1] * s)
    return np.c_[col, row]


def score(pts_px):
    x = pts_px[:, 0] - crop[0]
    y = pts_px[:, 1] - crop[1]
    ok = (x >= 0) & (y >= 0) & (x <= dt.shape[1] - 1) & (y <= dt.shape[0] - 1)
    v = np.full(len(x), TRUNC)
    v[ok] = np.minimum(ndi.map_coordinates(dt, [y[ok], x[ok]], order=1), TRUNC)
    return float(v.mean())


# initial pose: the inlet port's ring and the drawing's top; the mask's extent against the ink's
ys, xs = np.nonzero(ink)
P_top = P[np.isin(G, ['inlet and head', 'loops at the top'])]
c0 = crop[0] + 0.5 * (xs.min() + xs.max()) - k * 0.5 * (P[:, 1].min() + P[:, 1].max())
r0 = crop[1] + ys.min() - k * P[:, 0].min()
best = None
for dc in np.arange(-40, 40.01, 1.0):                   # a grid fine enough for the 1-px ridges (the top's shape fixes it)
    for dr in np.arange(-200, 200.01, 1.0):
        s_ = score(to_px(P_top, (c0 + dc, r0 + dr, 0.0)))
        if best is None or s_ < best[0]:
            best = (s_, c0 + dc, r0 + dr)
# the rotation from the whole device (its seven straight legs, 55-62 mm long, fix it); then the translation from the
# top alone, the rotation held; every later group is measured against that pose
res_all = None
for th0 in (-0.2, 0.0, 0.2):
    r_ = optimize.minimize(lambda q: score(to_px(P, (q[0], q[1], math.radians(q[2])))), [best[1], best[2], th0],
                           method='Nelder-Mead', options=dict(xatol=0.005, fatol=1e-6, maxiter=6000,
                                                              initial_simplex=[[best[1], best[2], th0],
                                                                               [best[1] + 1.0, best[2], th0],
                                                                               [best[1], best[2] + 1.0, th0],
                                                                               [best[1], best[2], th0 + 0.05]]))
    if res_all is None or r_.fun < res_all.fun:
        res_all = r_
pose_all = (res_all.x[0], res_all.x[1], math.radians(res_all.x[2]))
res = optimize.minimize(lambda q: score(to_px(P_top, (q[0], q[1], pose_all[2]))), [best[1], best[2]],
                        method='Nelder-Mead', options=dict(xatol=0.005, fatol=1e-6, maxiter=4000,
                                                           initial_simplex=[[best[1], best[2]], [best[1] + 1.0, best[2]],
                                                                            [best[1], best[2] + 1.0]]))
pose_top = (res.x[0], res.x[1], pose_all[2])
res_s = optimize.minimize(lambda q: score(to_px(P, (q[0], q[1], math.radians(q[2])), k * q[3])),
                          list(res_all.x) + [1.0], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-6, maxiter=6000))


def group_shift(g, pose):
    """the translation (px) that fits group g alone, the pose otherwise held: a grid, then Nelder-Mead.  A leg's
    straight part is invariant along its axis (rows of identical posts), so for a leg only the shift across is fitted"""
    Q = to_px(P[G == g], pose)
    Sg = {}
    best_ = None
    leg = g.startswith('leg')
    for dc in np.arange(-12, 12.01, 0.5):
        for dr in ([0.0] if leg else np.arange(-260, 260.01, 1.0)):
            v = score(Q + (dc, dr))
            Sg[(dc, dr)] = v
            if best_ is None or v < best_[0]:
                best_ = (v, dc, dr)
    if leg:
        r2 = optimize.minimize(lambda q: score(Q + (q[0], 0.0)), [best_[1]], method='Nelder-Mead',
                               options=dict(xatol=0.02, fatol=1e-6))
        return float(r2.x[0]), 0.0, float(r2.fun), float(score(Q)), int(len(Q))
    r2 = optimize.minimize(lambda q: score(Q + (q[0], q[1])), [best_[1], best_[2]], method='Nelder-Mead',
                           options=dict(xatol=0.02, fatol=1e-6))
    return float(r2.x[0]), float(r2.x[1]), float(r2.fun), float(score(Q)), int(len(Q))


shifts = {}
for g in groups:
    dc, dr, s1, s0, n = group_shift(g, pose_top)
    th = pose_top[2]
    # image px -> mask µm (along x, across y)
    d_across = (dc * math.cos(th) - dr * math.sin(th)) / k
    d_along = (dc * math.sin(th) + dr * math.cos(th)) / k
    shifts[g] = dict(points=n, shift_px=[dc, dr], shift_along_um=(None if g.startswith('leg') else d_along),
                     shift_across_um=d_across, chamfer_px_before=s0, chamfer_px_after=s1)


def inv_px(px, pose):
    c, rr, th = pose
    dc, dr = px[0] - c, px[1] - rr
    s, co = math.sin(th), math.cos(th)
    y = (dc * co - dr * s) / k
    x = (dc * s + dr * co) / k
    return np.array([x, y])


# ---- V12: the boxes against the traced panels
v12 = {}
for letter, (bc, br, hs) in BOXES.items():
    v12[letter] = dict(box_centre_px=[bc, br], half_size_px=hs, box_centre_mask_um=inv_px((bc, br), pose_top).tolist(),
                       half_size_um=hs / k)
for spec_ in TRACED:
    path, key, letter = spec_.split(':')
    R_ = json.load(open(path))
    rec_p = R_[key] if key not in ('-', '') else R_
    m = rec_p['map_chamfer']
    Mp = T.Map(m['origin_px'], m['ex_px_per_um'], m['ey_px_per_um'])
    c0_, r0_, c1_, r1_ = PANEL_CROPS[letter]
    centre_um = Mp.inverse(np.array([[0.5 * (c0_ + c1_), 0.5 * (r0_ + r1_)]]))[0]
    b = np.array(v12[letter]['box_centre_mask_um'])
    d = centre_um - b
    v12[letter].update(traced=os.path.basename(path), panel_centre_mask_um=centre_um.tolist(), offset_um=d.tolist(),
                       offset_over_half_size=float(np.max(np.abs(d)) / v12[letter]['half_size_um']),
                       ok=bool(np.max(np.abs(d)) <= v12[letter]['half_size_um']))

out = dict(raster=os.path.basename(RAS), raster_sha256=r.sha256, dxf=os.path.basename(DXF), dxf_sha256=mask.sha256,
           record=os.path.basename(REC), bar_px=bar[0], um_per_px=1.0 / k, truncation_px=TRUNC,
           leg_extent_um=leg_x, pose_top=dict(col=pose_top[0], row=pose_top[1], rotation_deg=math.degrees(pose_top[2]),
                                              chamfer_px=float(res.fun)),
           pose_whole=dict(col=pose_all[0], row=pose_all[1], rotation_deg=math.degrees(pose_all[2]),
                           chamfer_px=float(res_all.fun)),
           scale_free=dict(scale_over_bar=float(res_s.x[3]), chamfer_px=float(res_s.fun)),
           V11_group_shifts=shifts, V12_boxes=v12,
           mask_extent_um=dict(x=[float(P[:, 0].min()), float(P[:, 0].max())], y=[float(P[:, 1].min()), float(P[:, 1].max())]))
# the drawing as an underlay for full_device.py (its Placement reads map_chamfer, px_per_um and crop): the top-anchored
# pose, so the underlay shows where the drawing's device ends against the mask's
th_ = pose_top[2]
out['underlay'] = dict(name='xia2016_fig3a_whole_device', raster_sha256=r.sha256, crop=list(crop), px_per_um=k,
                       map_chamfer=dict(origin_px=[pose_top[0], pose_top[1]],
                                        ex_px_per_um=[k * math.sin(th_), k * math.cos(th_)],
                                        ey_px_per_um=[k * math.cos(th_), -k * math.sin(th_)]))
json.dump(out, open(OUT, 'w'), indent=1, default=float)

# ---- the figure: the drawing, the mask's edges at the top-anchored pose (cyan), each group's fitted shift (orange)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H, W = sub.shape
fig = plt.figure(figsize=(W / 150.0 * 2.2, H / 150.0 * 1.0), dpi=150)
for j, (pose, title) in enumerate(((pose_top, 'top-anchored pose'), (None, 'each group at its own fit'))):
    ax = fig.add_axes([j * 0.5, 0.0, 0.5, 0.97])
    ax.imshow(sub, cmap='gray', vmin=0, vmax=255, extent=(crop[0] - 0.5, crop[2] + 0.5, crop[3] + 0.5, crop[1] - 0.5),
              interpolation='none')
    for g in groups:
        Q = to_px(P[G == g], pose_top)
        if pose is None:
            Q = Q + shifts[g]['shift_px']
        ax.plot(Q[:, 0], Q[:, 1], ',', color=T.TRACE_COLOUR if pose is not None else '#ff6a00', alpha=0.6)
    for letter, (bc, br, hs) in BOXES.items():
        ax.add_patch(plt.Rectangle((bc - hs, br - hs), 2 * hs, 2 * hs, fill=False, ec='#ff00ff', lw=0.6))
        ax.text(bc + hs + 2, br, letter, color='#ff00ff', fontsize=5)
    ax.set_xlim(crop[0], crop[2])
    ax.set_ylim(crop[3], crop[1])
    ax.set_axis_off()
    ax.set_title(title, fontsize=6)
fig.savefig(PNG, dpi=150)
print('bar %.2f px (%.4f um/px); top pose: rot %.3f deg, chamfer %.3f px; whole: rot %.3f, chamfer %.3f; scale free %.4f'
      % (bar[0], 1 / k, math.degrees(pose_top[2]), res.fun, math.degrees(pose_all[2]), res_all.fun, res_s.x[3]))
for g in groups:
    s_ = shifts[g]
    print('%-24s n %6d  along %s  across %+6.1f um  chamfer %.3f -> %.3f px' % (
        g, s_['points'], ('%+8.1f um' % s_['shift_along_um']) if s_['shift_along_um'] is not None else '(invariant)',
        s_['shift_across_um'], s_['chamfer_px_before'], s_['chamfer_px_after']))
for letter, v in v12.items():
    if 'ok' in v:
        print('V12 %s: offset %s um, %.2f of the half-size -> %s' % (letter, np.round(v['offset_um'], 0),
                                                                        v['offset_over_half_size'], 'ok' if v['ok'] else 'FAIL'))
