"""Dinh 2024 Fig. 2B, fourth (outlet) panel: the outer side-channel walls' flare, measured on the publisher-embedded
raster (pdfimages page 5, image 0: 982 x 1940 px, identical to the audit package's dinh_fig2_original_embedded.png).
Every pick is printed.  Edge = the dark wall outline: per pixel row, the intensity minimum within +-6 px of the
previous row's pick, refined to sub-pixel by a parabola through the minimum and its two neighbours.
Scale: the 50 um bar (caption: "The scale bars are 50 μm").  Angles are measured from the straight wall's own
direction (the flow axis), not from the image's vertical; both walls are fitted separately.
Usage: python3 measure_dinh_flare.py <dinh_fig2_raster.png>"""
import sys, json, math
import numpy as np
from PIL import Image

IM = np.asarray(Image.open(sys.argv[1]).convert('L')).astype(float)

# ---- scale bar: the black bar under the panel (rows 1922-1940, columns 860-982) --------------------------------
bar_rows = []
for y in range(1922, 1940):
    d = np.nonzero(IM[y, 860:982] < 60)[0]
    if len(d) > 20:
        bar_rows.append((y, int(d.min() + 860), int(d.max() + 860), int(d.max() - d.min() + 1)))
bar_px = float(np.median([b[3] for b in bar_rows]))
UM = 50.0 / bar_px


def track(guess, y0, y1, win=5):
    """per row: the intensity minimum within +-win px of the guide line guess(y), refined to sub-pixel"""
    xs, ys = [], []
    for y in range(y0, y1 + 1):
        g = guess(y)
        lo, hi = int(round(g)) - win, int(round(g)) + win + 1
        row = IM[y, lo:hi]
        i = int(np.argmin(row))
        off = 0.0
        if 0 < i < len(row) - 1:
            a, b, c = row[i - 1], row[i], row[i + 1]
            den = a - 2 * b + c
            off = 0.5 * (a - c) / den if den > 0 else 0.0
        if i == 0 or i == len(row) - 1:              # a minimum pinned at the window's edge is not on the wall: dropped
            continue
        xs.append(lo + i + off)
        ys.append(y)
    return np.array(ys, float), np.array(xs)


def fit(ys, xs, a, b):
    m = (ys >= a) & (ys <= b)
    p = np.polyfit(ys[m], xs[m], 1)
    res = xs[m] - np.polyval(p, ys[m])
    return p, float(np.sqrt(np.mean(res ** 2))), int(m.sum())


out = dict(raster=sys.argv[1], raster_shape=list(IM.shape), scale_bar=dict(rows=bar_rows, length_px=bar_px,
                                                                              um_per_px=UM))
# ---- pills: the last post's downstream edge and the long walls' upstream end, on the left post column ------------
col = IM[1560:1760, 700:740].mean(axis=1)          # dark outline rows of the left pill column (x 700-740)
dark = np.nonzero(col < np.percentile(col, 15))[0] + 1560
out['pill_column_dark_rows'] = [int(v) for v in dark]
# guide lines (read off the raster by eye: straight wall, kink near row 1770, wall at the panel's foot), used only to
# keep each row's search window on the wall; the picks themselves are the intensity minima.  Rows 1636-1656 on the
# left are excluded from every fit: the caption's w_s(out) arrowhead touches that wall there; rows 1690-1702 on the
# right likewise: the caption's "52 um" label starts beside that wall there (v1 of this script kept two picks pinned on
# the label's "5" at rows 1696-1697; an independent check found them).  A pick at its window's edge is dropped.
GUIDE = {'left': lambda y: 653.5 if y < 1770 else 653.5 - (y - 1770) * (36.5 / 148.0),
         'right': lambda y: 886.2 + (y - 1585) * 0.0118 if y < 1768 else 888.4 + (y - 1768) * (39.8 / 147.0)}
EXCL = {'left': (1636, 1656), 'right': (1690, 1702)}
for side in ('left', 'right'):
    ys, xs = track(GUIDE[side], 1585, 1918)
    keep = ~((ys >= EXCL[side][0]) & (ys <= EXCL[side][1]))
    ys, xs = ys[keep], xs[keep]
    rec = dict(picks=[[int(y), round(float(x), 3)] for y, x in zip(ys, xs)])
    ps, rms_s, n_s = fit(ys, xs, 1590, 1750)          # the straight wall, well above the kink
    axis = math.atan(ps[0])
    results = []
    for a in (1785, 1790, 1795, 1800, 1810, 1820):
        for b in (1880, 1895, 1910, 1916):
            pf, rms_f, n_f = fit(ys, xs, a, b)
            ang = math.degrees(abs(math.atan(pf[0]) - axis))
            ang_img = math.degrees(abs(math.atan(pf[0])))
            yk = (ps[1] - pf[1]) / (pf[0] - ps[0])
            results.append(dict(fit_rows=[a, b], n=n_f, rms_px=round(rms_f, 3), angle_deg=round(ang, 2),
                                angle_from_image_vertical_deg=round(ang_img, 2), onset_row=round(float(yk), 1)))
    angs = [r['angle_deg'] for r in results]
    onsets = [r['onset_row'] for r in results]
    rec.update(straight=dict(fit_rows=[1590, 1750], slope=float(ps[0]), tilt_deg=round(math.degrees(axis), 3),
                             rms_px=round(rms_s, 3), n=n_s),
               flare_fits=results, angle_deg_range=[min(angs), max(angs)], angle_deg_median=float(np.median(angs)),
               onset_row_range=[min(onsets), max(onsets)], onset_row_median=float(np.median(onsets)))
    out[side] = {k: v for k, v in rec.items() if k != 'picks'}
    out[side + '_picks'] = rec['picks']
json.dump(out, open('dinh_flare_measurement.json', 'w'), indent=1)
print('scale bar %.1f px -> %.4f um/px' % (bar_px, UM))
print('pill column dark rows (outline rows of the last posts and the long wall top):', out['pill_column_dark_rows'])
for side in ('left', 'right'):
    r = out[side]
    print('%-5s straight wall tilt %+.3f deg (rms %.2f px); flare angle from the wall axis %.2f-%.2f deg (median %.2f) '
          'over %d fit intervals; onset row %.1f-%.1f (median %.1f)'
          % (side, r['straight']['tilt_deg'], r['straight']['rms_px'], r['angle_deg_range'][0], r['angle_deg_range'][1],
             r['angle_deg_median'], len(r['flare_fits']), r['onset_row_range'][0], r['onset_row_range'][1],
             r['onset_row_median']))
