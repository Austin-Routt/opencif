"""Dinh 2024 Fig. 2A inset: the inlet tube, chevron filter, bulb with cross row, funnel and stem.
Native raster dinh-001.png; scale 9.89 um/px (inset width = 3.2 mm chip, cross-checked by the row pitch 2(GS+wall))."""
import json
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
IM = np.asarray(Image.open('../r3/papers/dinh-001.png').convert('L')).astype(float)
SC = 9.89
XC = 192.6                     # trunk axis (native px), from the trunk profile at y = 640
Y_ROW1 = 599.5                 # row-1 mouth (native px)
y0, y1 = 185, 600
prof = []
for y in range(y0, y1):
    row = IM[y, 120:270]
    blk = row < 128
    # contiguous black run containing the axis (or nearest to it)
    xs = np.nonzero(blk)[0] + 120
    if len(xs) == 0:
        continue
    # sub-pixel edges by linear interpolation at level 128 around the outermost black pixels connected to the axis
    lab, n = ndi.label(blk)
    ia = int(round(XC)) - 120
    l = lab[ia] if lab[ia] else None
    if l is None:
        continue
    idx = np.nonzero(lab == l)[0]
    a, b = idx[0], idx[-1]
    def edge(i, j):          # i black, j white neighbour
        vi, vj = row[i], row[j]
        return i + (128 - vi) / (vj - vi) * (j - i)
    xl = edge(a, a - 1) + 120 if a > 0 else a + 120
    xr = edge(b, b + 1) + 120 if b < len(row) - 1 else b + 120
    prof.append((y, xl, xr))
prof = np.array(prof)
half = 0.5 * (prof[:, 2] - prof[:, 1])
xs_um = (prof[:, 0] - Y_ROW1) * SC
out = dict(scale_um_px=SC, y_row1_px=Y_ROW1, axis_px=XC,
           profile=[(float(x), float(h * SC), float(((l + r) / 2 - XC) * SC)) for x, h, l, r in
                    zip(xs_um, half, prof[:, 1], prof[:, 2])])
# white islands (posts) inside the black inlet structure
sub = IM[y0:y1, 120:270]
white = sub > 128
lab, n = ndi.label(white)
posts = []
for i, sl in enumerate(ndi.find_objects(lab)):
    ys, xs = sl
    if ys.start == 0 or xs.start == 0 or ys.stop == sub.shape[0] or xs.stop == sub.shape[1]:
        continue
    m = lab[sl] == i + 1
    area = m.sum()
    if area < 2:
        continue
    yy, xx = np.nonzero(m)
    cy, cx = yy.mean() + ys.start + y0, xx.mean() + xs.start + 120
    # orientation from second moments
    cov = np.cov(np.vstack([xx, yy])) if area > 2 else np.eye(2)
    ev, evec = np.linalg.eigh(cov)
    ang = float(np.degrees(np.arctan2(evec[1, 1], evec[0, 1])))
    posts.append(dict(x_um=float((cy - Y_ROW1) * SC), y_um=float((cx - XC) * SC), area_um2=float(area * SC * SC),
                      len_um=float(4 * np.sqrt(max(ev[1], 0)) * SC), wid_um=float(4 * np.sqrt(max(ev[0], 0)) * SC),
                      ang_deg=ang))
out['posts'] = posts
json.dump(out, open('dinh_filter_meas.json', 'w'), indent=1)
# summary
P = np.array([(p['x_um'], p['y_um'], p['len_um'], p['wid_um'], p['ang_deg']) for p in posts])
print('profile rows', len(prof), 'posts', len(posts))
for x, h, c in out['profile'][::12]:
    print('  x %8.1f um  half %6.1f um  centre %5.1f' % (x, h, c))
print('posts (x, y, len, wid, ang):')
for r in P[np.argsort(P[:, 0])]:
    print('  %8.1f %7.1f  %5.1f %5.1f %6.1f' % tuple(r))
