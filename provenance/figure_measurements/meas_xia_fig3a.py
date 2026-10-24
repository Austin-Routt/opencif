"""meas_xia_fig3a.py -- Xia 2016 Fig. 3a (the RS-7 layout) measured at its own 5 mm scale bar.

Raster: pdfimages page xia-002.png, panel a = columns 760:1270 (saved as r3/x3a.png, 510 x 3138 px).
Scale: the 5 mm bar is the black run at rows 2940-2950, columns 267-445 of the crop = 179 px -> 27.93 um/px.
(An earlier figure script searched rows 2600-2900 for the bar, found none and fell back to 160 px; that scale,
31.25 um/px, is withdrawn.)  Check: leg 7's outer width at its outlet end reads 1508 um against the design's 1502.

Lines in this CAD drawing are 3-8 px thick; every position below is a line CENTRE (a wall boundary), +-1 px."""
import json
import numpy as np
from PIL import Image
from scipy.optimize import least_squares

IM = np.asarray(Image.open(__file__.rsplit('/', 2)[0] + '/r3/x3a.png').convert('L')).astype(float)
DARK = IM < 110
out = {}

# ---- scale bar
bar = 0
for y in range(2900, 3000):                      # the longest contiguous dark run (the bar is 15 rows thick)
    xs = np.where(DARK[y])[0]
    if len(xs):
        bar = max([bar] + [len(r) for r in np.split(xs, np.where(np.diff(xs) > 1)[0] + 1)])
UM = 5000.0 / bar
out['bar_px'] = int(bar)
out['um_per_px'] = UM


def runs_at(y, x0=120, x1=470):
    xs = np.where(DARK[y, x0:x1])[0] + x0
    if not len(xs):
        return []
    return [(float(r[0]), float(r[-1])) for r in np.split(xs, np.where(np.diff(xs) > 1)[0] + 1)]


def fitc(pts, g):
    pts = np.asarray(pts, float)
    r = least_squares(lambda p: np.hypot(pts[:, 0] - p[0], pts[:, 1] - p[1]) - p[2], g)
    return r.x, float(np.sqrt(np.mean(r.fun ** 2)))


# ---- leg axes: each leg's two post-row lines (or merged wall+post lines on the narrow early legs), median over rows
# 600-2300 (the straight legs).  Windows (px) bracket each leg's pair of inner lines.
WIN = {7: (150, 175), 5: (200, 218), 3: (235, 249), 1: (290, 318), 2: (338, 365), 4: (377, 390), 6: (412, 428)}
cent = {}
for leg, (a, b) in WIN.items():
    cs = []
    for y in range(600, 2300, 10):
        rr = [((r0 + r1) / 2.0) for r0, r1 in runs_at(y, a, b)]
        if len(rr) >= 2:
            cs.append(0.5 * (rr[0] + rr[-1]))
    cent[leg] = float(np.median(cs))
out['leg_axis_px'] = cent
out['leg_Y_um'] = {leg: (cent[leg] - cent[1]) * UM for leg in sorted(cent)}

# ---- outer widths of the legs with distinct outer walls, near their two ends (rows 440-460 and 2390-2410).  A wall
# line is 4 px thick where it stands alone; where a wall line merges with a post-row line the run's CENTRE is biased
# inward, so each wall is located from the run's OUTER edge minus half a line (2 px): width = outer extent - 4 px
LINE_PX = 4.0
OUTER = {7: (125, 195), 5: (192, 229), 3: (228, 256), 4: (368, 399), 6: (398, 444)}
ow = {}
for leg, (a, b) in OUTER.items():
    ww = {}
    for y0 in (450, 2400):
        vals = []
        for y in range(y0 - 10, y0 + 11):
            rr = runs_at(y, a, b)
            if rr:
                vals.append((rr[-1][1] - rr[0][0] + 1.0 - LINE_PX) * UM)
        ww[y0] = float(np.median(vals))
    ow[leg] = ww
out['leg_outer_width_um_rows_450_2400'] = ow
out['leg_outer_width_note'] = 'median over rows y0 +- 10; wall located at the run outer edge minus 2 px (line 4 px)'

# ---- loops: centreline radius = half the axis spacing of the two legs a loop joins; centre rows from the loop arcs
out['loops'] = {}
for k, (a, b) in {1: (1, 2), 2: (2, 3), 3: (3, 4), 4: (4, 5), 5: (5, 6), 6: (6, 7)}.items():
    out['loops'][k] = dict(legs=(a, b), centre_x_px=0.5 * (cent[a] + cent[b]),
                           r_um=0.5 * abs(cent[b] - cent[a]) * UM, end='bottom' if k % 2 else 'top')
# loop-centre rows: the top loops' arcs cross column 296 at rows 126-231, the bottom loops' at column 320 at rows
# 2482-2592; each loop's four lines are averaged to its centreline and the radius subtracted
col_top = [((r0 + r1) / 2.0) for r0, r1 in [(126, 128), (141, 143), (149, 151), (162, 165), (176, 179), (183, 185),
                                             (191, 194), (197, 200), (209, 212), (214, 216), (224, 231)]]
top_cl = {6: np.mean(col_top[0:4]), 4: np.mean(col_top[4:8]), 2: np.mean(col_top[8:11])}
bot = [((r0 + r1) / 2.0) for r0, r1 in [(2482, 2489), (2499, 2507), (2532, 2535), (2537, 2540), (2545, 2547),
                                         (2550, 2553), (2562, 2566), (2572, 2574), (2580, 2582), (2589, 2592)]]
bot_cl = {1: np.mean(bot[0:2]), 3: np.mean(bot[2:6]), 5: np.mean(bot[6:10])}
rows_top = {k: float(top_cl[k] + out['loops'][k]['r_um'] / UM) for k in top_cl}
rows_bot = {k: float(bot_cl[k] - out['loops'][k]['r_um'] / UM) for k in bot_cl}
out['loop_centre_row_top'] = rows_top
out['loop_centre_row_bottom'] = rows_bot
x_lo, x_hi = float(np.mean(list(rows_top.values()))), float(np.mean(list(rows_bot.values())))
out['leg_straight_length_um'] = (x_hi - x_lo) * UM

# ---- circles: inlet port, centre port, filtrate port, ring walls
def centres_in(y0, y1, x0, x1, maxw=8):
    P = []
    for y in range(y0, y1):
        for r0, r1 in runs_at(y, x0, x1):
            if r1 - r0 + 1 <= maxw:
                P.append(((r0 + r1) / 2.0, y))
    return np.array(P)


pin = centres_in(334, 400, 262, 338)
c, rms = fitc(pin, [305, 366, 31])
out['inlet_port'] = dict(centre_px=c[:2].tolist(), r_um=c[2] * UM, rms_px=rms)
P = centres_in(2600, 2950, 0, 420)
for name, g in (('centre_port', [162, 2700, 61]), ('filtrate_port', [162, 2872, 61]),
                ('ring_outer_wall', [162, 2745, 155]), ('ring_inner_wall', [162, 2744, 144])):
    sel = P[np.abs(np.hypot(P[:, 0] - g[0], P[:, 1] - g[1]) - g[2]) < 5]
    c, rms = fitc(sel, g)
    out[name] = dict(centre_px=c[:2].tolist(), r_um=c[2] * UM, rms_px=rms, n=int(len(sel)))
cp, fp = np.array(out['centre_port']['centre_px']), np.array(out['filtrate_port']['centre_px'])
ro, ri = out['ring_outer_wall'], out['ring_inner_wall']
rc_ = 0.5 * (np.array(ro['centre_px']) + np.array(ri['centre_px']))
out['ring'] = dict(R_centreline_um=0.5 * (ro['r_um'] + ri['r_um']), width_um=ro['r_um'] - ri['r_um'],
                   centre_beyond_centre_port_um=(rc_[1] - cp[1]) * UM,
                   filtrate_port_beyond_centre_port_um=(fp[1] - cp[1]) * UM,
                   top_beyond_bottom_loops_um=(rc_[1] - 0.5 * (ro['r_um'] + ri['r_um']) / UM - x_hi) * UM)
out['inlet_beyond_top_loops_um'] = (out['inlet_port']['centre_px'][1] - x_lo) * UM
out['inlet_lateral_from_leg1_um'] = (out['inlet_port']['centre_px'][0] - cent[1]) * UM
# ---- overall extents: top loop's outer wall to the filtrate port's bottom; ring's left wall to leg 6's right wall
top = 126.0 + 1.0
bottom = fp[1] + out['filtrate_port']['r_um'] / UM
left = np.array(ro['centre_px'])[0] - ro['r_um'] / UM
right = max(np.mean([(r0 + r1) / 2.0 for r0, r1 in runs_at(y, 425, 445)]) for y in (450,))
out['extent_um'] = dict(length=(bottom - top) * UM, width=(right - left) * UM)

# ---- the head: its band's outer edge along leg 1, its top, and the first pill (box c located on the Fig. 3c field)
ax_ = out['inlet_port']['centre_px'][0]
yin = out['inlet_port']['centre_px'][1]


def band_edge(y):
    """half-width (px) of the dark band round leg 1's trunk at row y: outer edges of the dark runs within 48 px"""
    prof = IM[y - 2:y + 3].mean(axis=0)
    xl = int(ax_) - 6
    while prof[xl] > 128 and xl > 256:
        xl -= 1
    while prof[xl] <= 128 and xl > 256:
        xl -= 1
    xr = int(ax_) + 6
    while prof[xr] > 128 and xr < 345:
        xr += 1
    while prof[xr] <= 128 and xr < 345:
        xr += 1
    return 0.5 * (xr - xl)


prof = [(float((y - yin) * UM / 1000.0), float(band_edge(y) * UM)) for y in range(430, 900, 12)]
out['head_band_profile_mm_um'] = prof
top_row = next(y for y in range(395, 440) if DARK[y, int(ax_) - 32:int(ax_) - 18].any())
out['head_band_top_row_px'] = int(top_row)
out['head_top_below_inlet_mm'] = (top_row - yin) * UM / 1000.0 + 2.0 * UM / 1000.0    # + half a drawn line
box = [y for y in range(800, 900) if (DARK[y, 280:327].sum() >= 40)]
bc_top, bc_bot = min(box), max(box)
out['box_c_rows_px'] = [int(bc_top), int(bc_bot)]
# Fig. 3c (r3/x3c.png, 1.953 um/px): the first pill's centre row, from the outlines of the left post column
IMC = np.asarray(Image.open(__file__.rsplit('/', 2)[0] + '/r3/x3c.png').convert('L')).astype(float)
pc = IMC[:, 228:252].mean(axis=1)
dd = np.where(pc < np.percentile(pc, 25))[0]
cent = [float(r.mean()) for r in np.split(dd, np.where(np.diff(dd) > 1)[0] + 1)]
gaps = np.diff(cent)
# blocks: outline pairs ~10 px apart separated by ~55 px; pills: pairs ~6-7 px apart separated by ~22 px
j = next(i for i in range(1, len(gaps) - 2) if gaps[i] < 30 and gaps[i - 1] < 15 and gaps[i + 1] < 15 and gaps[i + 2] < 30)
first_pill_c = 0.5 * (cent[j] + cent[j + 1])                  # the first pill: top outline cent[j], bottom cent[j+1]
SC3 = 250.0 / 128.0
out['fig3c_first_pill_row_px'] = first_pill_c
pill_mm = ((bc_top + bc_bot) / 2.0 - yin) * UM / 1000.0 + (first_pill_c - IMC.shape[0] / 2.0) * SC3 / 1000.0
out['first_pill_below_inlet_mm'] = pill_mm
out['head_length_mm'] = pill_mm - out['head_top_below_inlet_mm']
out['head_note'] = ('head top: the band\'s first dark row beside the stem (rendered edge) plus half a drawn line (2 px); '
                    'first pill: box c\'s centre in Fig. 3a mapped onto the centre of the Fig. 3c field, plus the pill\'s '
                    'offset from that centre in Fig. 3c; uncertainty +-0.2 mm (the box is 1.23 mm tall, the field 1.47 mm)')
json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1, default=float)
for k, v in out.items():
    print(k, v)
