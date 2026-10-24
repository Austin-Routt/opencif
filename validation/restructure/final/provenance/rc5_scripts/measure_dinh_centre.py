"""rc5 (after the rc5 reviews): two more measurements on Dinh 2024 Fig. 2B's fourth (outlet) panel, on the same
publisher raster as measure_dinh_flare.py (pdfimages page 5, image 0: 982 x 1940 px), with every pick printed:
  (1) the centre channel after the last gap: its two edges (the facing outlines of the two long walls that follow the
      last posts) per pixel row, from the long walls' rounded ends to the panel's foot, as a width against the distance
      from the last gap's centre; a plateau-then-linear fit gives where it starts to widen.  Calibration: the same
      outline-to-outline width between the two post columns at three post rows clear of the caption's arrows, against
      the printed 52 um ("The middle channel narrowed from 75 um near the inlet to 52 um near the outlet");
  (2) the outer side-channel walls' straight-part tilt over sub-ranges of rows (Reviewer 1's E3: the right wall's
      tilt is the least certain input of the flare angle), and the flare angle from each wall's own direction for
      each sub-range.
Edges are the dark outlines: per row, the intensity minimum within +-win px of the previous row's pick (or a guide),
refined to sub-pixel by a parabola; a minimum pinned at its window's edge is dropped.  Per-edge pick uncertainty is
taken as +-0.5 px; the 50 um bar is 46 +- 1 px (+-2.2 % on every length).
Usage: python3 measure_dinh_centre.py <dinh_fig2_raster.png> [OUT.json]"""
import json, math, sys
import numpy as np
from PIL import Image

IM = np.asarray(Image.open(sys.argv[1]).convert('L')).astype(float)
OUT = sys.argv[2] if len(sys.argv) > 2 else 'dinh_centre_measurement.json'
bar = [int(np.ptp(np.nonzero(IM[y, 860:982] < 60)[0])) + 1 for y in range(1922, 1940)
       if len(np.nonzero(IM[y, 860:982] < 60)[0]) > 20]
BAR_PX = float(np.median(bar))
UM = 50.0 / BAR_PX


def pick(y, g, win):
    lo = int(round(g)) - win
    row = IM[y, lo:lo + 2 * win + 1]
    i = int(np.argmin(row))
    if i == 0 or i == len(row) - 1:
        return None
    a, b, c = row[i - 1], row[i], row[i + 1]
    den = a - 2 * b + c
    return lo + i + (0.5 * (a - c) / den if den > 0 else 0.0)


def track(x0, y0, y1, win=4, guide=None):
    ys, xs, g = [], [], x0
    for y in range(y0, y1 + 1):
        gg = guide(y) if guide else g
        x = pick(y, gg, win)
        if x is None:
            continue
        ys.append(y); xs.append(x); g = x
    return np.array(ys, float), np.array(xs)


out = dict(raster=sys.argv[1], raster_shape=list(IM.shape), scale_bar_px=BAR_PX, um_per_px=UM,
           pick_uncertainty_px_per_edge=0.5, scale_uncertainty_rel=1.0 / BAR_PX)
# ---- the last gap: between the last small post's downstream outline and the long wall's upstream outline -----------
col = IM[1560:1760, 700:740].mean(axis=1)
dark = list(np.nonzero(col < np.percentile(col, 15))[0] + 1560)
clusters, cur = [], [dark[0]]
for v in dark[1:]:
    if v - cur[-1] <= 2:
        cur.append(v)
    else:
        clusters.append(cur); cur = [v]
clusters.append(cur)
post_end, wall_start = float(np.mean(clusters[-2])), float(np.mean(clusters[-1]))
y_gap = 0.5 * (post_end + wall_start)
out.update(outline_row_clusters=[[int(c[0]), int(c[-1])] for c in clusters], last_post_downstream_outline_row=post_end,
           long_wall_upstream_outline_row=wall_start, last_gap_centre_row=y_gap)

# ---- calibration: the centre width between the post columns at three post rows clear of the arrows ---------------
cal = []
for yc in (1518, 1560, 1603):                    # mid-rows of posts 1-3 of the panel (the arrows cross posts 4 and 5)
    for y in range(yc - 3, yc + 4):
        xl, xr = pick(y, 746.5, 4), pick(y, 796.5, 4)
        if xl is not None and xr is not None:
            cal.append([y, round(xl, 3), round(xr, 3), round((xr - xl) * UM, 3)])
cw = np.array([c[3] for c in cal])
out['calibration_between_post_columns'] = dict(picks=cal, width_um_mean=float(cw.mean()), width_um_sd=float(cw.std(ddof=1)),
                                               printed_um=52.0, offset_um=float(cw.mean() - 52.0))

# ---- the centre channel after the last gap -----------------------------------------------------------------------
y0 = int(math.ceil(wall_start)) + 12               # below the long walls' rounded ends
ysl, xsl = track(749.0, y0, 1922)
ysr, xsr = track(795.5, y0, 1922)
common = sorted(set(ysl.astype(int)) & set(ysr.astype(int)))
L = dict(zip(ysl.astype(int), xsl)); Rr = dict(zip(ysr.astype(int), xsr))
rows = [[int(y), round(L[y], 3), round(Rr[y], 3), round((y - y_gap) * UM, 2), round((Rr[y] - L[y]) * UM, 3)] for y in common]
s = np.array([r[3] for r in rows]); w = np.array([r[4] for r in rows])
best = None
for k in range(10, len(s) - 5):                    # plateau then linear, continuous at the break s[k]
    sb = s[k]
    X = np.column_stack([np.ones_like(s), np.maximum(s - sb, 0.0)])
    coef, *_ = np.linalg.lstsq(X, w, rcond=None)
    sse = float(np.sum((X @ coef - w) ** 2))
    if best is None or sse < best[0]:
        best = (sse, sb, coef)
sse, sb, coef = best
sig_w = math.sqrt(2.0) * 0.5 * UM
plateau = w[s < sb]
out['centre_channel_after_last_gap'] = dict(
    rows_picked=len(rows), first_row=int(common[0]), last_row=int(common[-1]),
    s_um_range=[float(s.min()), float(s.max())],
    plateau_width_um_mean=float(plateau.mean()), plateau_width_um_sd=float(plateau.std(ddof=1)),
    widening_starts_um_after_last_gap_centre=float(sb), widening_slope_um_per_um=float(coef[1]),
    widening_angle_per_side_deg=float(math.degrees(math.atan(coef[1] / 2.0))),
    width_um_at_last_row=float(w[-1]), width_um_last_5_rows_mean=float(w[-5:].mean()),
    left_edge_shift_px_last_vs_plateau=float(xsl[-1] - np.mean([L[y] for y in common if (y - y_gap) * UM < sb])),
    right_edge_shift_px_last_vs_plateau=float(xsr[-1] - np.mean([Rr[y] for y in common if (y - y_gap) * UM < sb])),
    width_pick_uncertainty_um=sig_w, rms_fit_um=math.sqrt(sse / len(s)))
out['centre_channel_picks_[row, x_left, x_right, s_um, width_um]'] = rows

# ---- the outer walls' straight-part tilt over sub-ranges, and the flare angle from each ---------------------------
GUIDE = {'left': lambda y: 653.5 if y < 1770 else 653.5 - (y - 1770) * (36.5 / 148.0),
         'right': lambda y: 886.2 + (y - 1585) * 0.0118 if y < 1768 else 888.4 + (y - 1768) * (39.8 / 147.0)}
EXCL = {'left': (1636, 1656), 'right': (1690, 1702)}  # the caption's arrowhead / "52 um" label touch the walls there
tilt = {}
for side in ('left', 'right'):
    ys, xs = track(None, 1500, 1918, win=5, guide=GUIDE[side])
    keep = ~((ys >= EXCL[side][0]) & (ys <= EXCL[side][1]))
    ys, xs = ys[keep], xs[keep]
    mf = (ys >= 1790) & (ys <= 1910)
    pf = np.polyfit(ys[mf], xs[mf], 1)
    flare_img = math.degrees(math.atan(pf[0]))
    subs = []
    for a, b in ((1500, 1585), (1590, 1640), (1645, 1689), (1703, 1750), (1590, 1750), (1500, 1750)):
        m = (ys >= a) & (ys <= b)
        p = np.polyfit(ys[m], xs[m], 1)
        t = math.degrees(math.atan(p[0]))
        subs.append(dict(rows=[a, b], n=int(m.sum()), tilt_deg=round(t, 3),
                         rms_px=round(float(np.sqrt(np.mean((xs[m] - np.polyval(p, ys[m])) ** 2))), 3),
                         flare_from_this_wall_deg=round(abs(flare_img - t), 3)))
    tilts = [d['tilt_deg'] for d in subs]
    angs = [d['flare_from_this_wall_deg'] for d in subs]
    tilt[side] = dict(flare_from_image_vertical_deg=round(flare_img, 3), sub_ranges=subs,
                      tilt_deg_range=[min(tilts), max(tilts)], flare_from_wall_deg_range=[min(angs), max(angs)])
out['outer_wall_tilt'] = tilt
json.dump(out, open(OUT, 'w'), indent=1)
c = out['centre_channel_after_last_gap']
k = out['calibration_between_post_columns']
print('scale %.1f px / 50 um -> %.4f um/px; last gap centre at row %.1f (post outline %.1f, long wall outline %.1f)'
      % (BAR_PX, UM, y_gap, post_end, wall_start))
print('calibration: between the post columns %.2f +- %.2f um (sd, %d rows) against the printed 52 um (offset %+.2f um)'
      % (k['width_um_mean'], k['width_um_sd'], len(cal), k['offset_um']))
print('centre channel after the last gap: rows %d-%d (s %.1f-%.1f um): %.2f +- %.2f um (sd) to s = %.1f um, then widening'
      ' %.3f um/um (%.2f deg per side) to %.2f um at the last row (last 5 rows %.2f); pick uncertainty +-%.2f um, scale +-2.2 %%'
      % (c['first_row'], c['last_row'], c['s_um_range'][0], c['s_um_range'][1], c['plateau_width_um_mean'],
         c['plateau_width_um_sd'], c['widening_starts_um_after_last_gap_centre'], c['widening_slope_um_per_um'],
         c['widening_angle_per_side_deg'], c['width_um_at_last_row'], c['width_um_last_5_rows_mean'], c['width_pick_uncertainty_um']))
print('  edge shifts at the last row against the plateau: left %+.2f px, right %+.2f px'
      % (c['left_edge_shift_px_last_vs_plateau'], c['right_edge_shift_px_last_vs_plateau']))
for side in ('left', 'right'):
    t = tilt[side]
    print('%-5s wall: flare %.2f deg from the image vertical; straight-part tilt %+.3f to %+.3f deg over %d sub-ranges -> '
          'flare %.2f-%.2f deg from the wall' % (side, t['flare_from_image_vertical_deg'], t['tilt_deg_range'][0],
                                                  t['tilt_deg_range'][1], len(t['sub_ranges']),
                                                  t['flare_from_wall_deg_range'][0], t['flare_from_wall_deg_range'][1]))
    for d in t['sub_ranges']:
        print('      rows %d-%d (n %d): tilt %+.3f deg, rms %.3f px, flare from it %.2f deg' % (
            d['rows'][0], d['rows'][1], d['n'], d['tilt_deg'], d['rms_px'], d['flare_from_this_wall_deg']))
