# rc7 (code review R11): shipped with the tree; copied unchanged from the rc6 device-package tools (devpkg/tools/measure_p1.py).
# usage: python3 -I measure_p1.py DINH5_RASTER DINH_RESULTS.json  (the raster is the paper's image: not in the repo)
"""measure the Dinh Fig. 2B panel-1 features the rc6 head does not draw (maintenance tool, reads the paper image)"""
import sys, json
import numpy as np
from PIL import Image
im = np.asarray(Image.open(sys.argv[1]).convert('L'), float)
rec = json.load(open(sys.argv[2]))['p1']
mc = rec['map_chamfer'] if 'map_chamfer' in rec else rec['maps']['chamfer']
O = np.array(mc['origin_px']); ex = np.array(mc['ex_px_per_um']); ey = np.array(mc['ey_px_per_um'])
A = np.column_stack([ex, ey])
def to_um(c, r):
    return np.linalg.solve(A, np.array([c, r], float) - O)
def dark_min(profile, lo):
    """sub-pixel positions of local minima darker than the median - 25"""
    p = np.asarray(profile); out = []
    thr = np.median(p) - 25
    for i in range(1, len(p) - 1):
        if p[i] < thr and p[i] <= p[i - 1] and p[i] < p[i + 1]:
            a, b, c = p[i - 1], p[i], p[i + 1]
            d = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0.0
            out.append(lo + i + d)
    return out
print('frame', rec.get('frame') or rec.get('frame_px_c0_r0_c1_r1'))
# rectangles: top/bottom edges along columns inside them; inner end along rows
for side, cols, rows in (('left', range(614, 664, 6), range(185, 250, 6)), ('right', range(876, 936, 6), range(185, 250, 6))):
    tops, bots = [], []
    for c in cols:
        L = dark_min(im[150:290, c], 150)
        L = [v for v in L if 165 < v < 265]
        if len(L) >= 2:
            tops.append(min(L)); bots.append(max(L))
    ends = []
    for r in rows:
        if side == 'left':
            L = [v for v in dark_min(im[r, 650:720], 650) if 655 < v < 700]
            if L: ends.append(max(L))
        else:
            L = [v for v in dark_min(im[r, 840:900], 840) if 850 < v < 895]
            if L: ends.append(min(L))
    t, b, e = np.median(tops), np.median(bots), np.median(ends)
    c_mid = np.median(list(cols))
    pt, pb = to_um(c_mid, t), to_um(c_mid, b)
    pe = to_um(e, 0.5 * (t + b))
    print(side, 'edges rows %.2f / %.2f (n %d), along x %.1f .. %.1f um (%.1f um); inner end col %.2f (n %d) -> |y| %.1f um' % (
        t, b, len(tops), pt[0], pb[0], pb[0] - pt[0], e, len(ends), abs(pe[1])))


# ---- corner radii: the outline's distance from the apex of its two straight edges, along the bisector
from scipy.ndimage import map_coordinates
def sample(c, r):
    return map_coordinates(im, [np.atleast_1d(r), np.atleast_1d(c)], order=1)
def darkest(p0, n, half=5.0, step=0.1):
    t = np.arange(-half, half + step, step)
    v = sample(p0[0] + t * n[0], p0[1] + t * n[1])
    i = int(np.argmin(v))
    return p0 + t[i] * n, v[i]
def corner(apex, d1, d2, k=1.066):
    apex = np.array(apex, float); d1 = np.array(d1, float); d2 = np.array(d2, float)
    lines = []
    for d in (d1, d2):
        n = np.array([-d[1], d[0]])
        pts = np.array([darkest(apex + t * d, n)[0] for t in np.arange(9, 24, 1.0)])
        # least-squares line through pts: point + direction
        m = pts.mean(0); u = np.linalg.svd(pts - m)[2][0]
        lines.append((m, u))
    (m1, u1), (m2, u2) = lines
    s = np.linalg.solve(np.column_stack([u1, -u2]), m2 - m1)
    A = m1 + s[0] * u1
    b = (d1 + d2) / np.linalg.norm(d1 + d2)
    t = np.arange(0, 14, 0.05)
    v = sample(A[0] + t * b[0], A[1] + t * b[1])
    i = int(np.argmin(v))
    d = t[i]
    return A, d, d * k / (np.sqrt(2) - 1)
for name, apex, d1, d2 in (
        ('left: the turn into branch 2, outer corner (C1 -> B2)', (685, 332), (-1, 0), (0, 1)),
        ('right: the turn into branch 2, outer corner', (855, 332), (1, 0), (0, 1)),
        ('left: the meander wall\'s end, upstream corner', (656, 358), (-1, 0), (0, 1)),
        ('left: the meander wall\'s end, downstream corner', (656, 383), (-1, 0), (0, -1)),
        ('right: the meander wall\'s end, upstream corner', (884, 358), (1, 0), (0, 1)),
        ('right: the meander wall\'s end, downstream corner', (884, 383), (1, 0), (0, -1)),
        ('left: rectangle inner end, upstream corner', (671, 179), (-1, 0), (0, 1)),
        ('left: rectangle inner end, downstream corner', (671, 257), (-1, 0), (0, -1)),
        ('right: rectangle inner end, upstream corner', (868, 177), (1, 0), (0, 1)),
        ('right: rectangle inner end, downstream corner', (868, 255), (1, 0), (0, -1))):
    A, d, r = corner(apex, d1, d2)
    print('%-58s apex (%.1f, %.1f) px; outline %.2f px from the apex -> radius %.1f um' % (name, A[0], A[1], d, r))


# ---- the meander walls' ends: the end outline x(y) across the wall, against a semicircle and rounded squares
def end_profile(c_guess, rows, sign):
    out = []
    for r in rows:
        p, v = darkest(np.array([c_guess, r], float), np.array([1.0, 0.0]), half=9.0)
        out.append((r, p[0]))
    return np.array(out)
def model(y, y0, y1, tip, r, sign):
    """end outline column at row y for a wall spanning rows y0..y1 whose end is at column tip, corner radius r (px)"""
    yy = np.minimum(y - y0, y1 - y)            # distance from the nearer long edge
    dx = np.where(yy < r, r - np.sqrt(np.maximum(r * r - (r - yy) ** 2, 0)), 0.0)
    return tip - sign * dx
for name, c_guess, sign in (('left wall end', 657, 1), ('right wall end', 882, -1)):
    rows = np.arange(357.0, 385.0, 0.5)
    prof = end_profile(c_guess, rows, sign)
    # the wall's long edges: dark lines along columns a few px back from the end
    cs = c_guess - sign * np.arange(8.0, 15.0)
    top = np.median([darkest(np.array([c, 358.0]), np.array([0.0, 1.0]), half=4)[0][1] for c in cs])
    bot = np.median([darkest(np.array([c, 383.0]), np.array([0.0, 1.0]), half=4)[0][1] for c in cs])
    best = None
    for rpx in np.arange(1.0, (bot - top) / 2 + 0.01, 0.25):
        for tip in np.arange(c_guess - 4, c_guess + 4, 0.05):
            m = model(prof[:, 0], top, bot, tip, rpx, sign)
            inside = (prof[:, 0] > top + 0.5) & (prof[:, 0] < bot - 0.5)
            e = np.sqrt(np.mean((prof[inside, 1] - m[inside]) ** 2))
            if best is None or e < best[0]:
                best = (e, rpx, tip)
    semi = None
    for tip in np.arange(c_guess - 4, c_guess + 4, 0.05):
        m = model(prof[:, 0], top, bot, tip, (bot - top) / 2, sign)
        inside = (prof[:, 0] > top + 0.5) & (prof[:, 0] < bot - 0.5)
        e = np.sqrt(np.mean((prof[inside, 1] - m[inside]) ** 2))
        semi = e if semi is None else min(semi, e)
    print('%s: wall %.1f um across (rows %.1f..%.1f); best corner radius %.1f um (rms %.2f px); a full semicircle (r %.1f um): rms %.2f px' % (
        name, (bot - top) * 1.066, top, bot, best[1] * 1.066, best[0], (bot - top) / 2 * 1.066, semi))
