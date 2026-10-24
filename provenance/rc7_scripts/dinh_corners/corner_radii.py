"""rc7 (A. H. Routt, 00:23: "all inner corners that impact particle flow are necessarily still rounded"): the rounding of
every right-angle corner of the drawn head on Dinh Fig. 2B's first panel, on the device (the image) and on the mask, on
the trace record's own registration (its chamfer map; nothing re-fitted).

For each corner of the mask's edges in the panel's region (a turn of 75-105 degrees within 30 µm of path, with straight
runs on both sides): the two straight runs' lines meet at the corner's virtual apex; the rounding r follows from the
apex's distance to the outline along the bisector, d = (sqrt 2 - 1) r.  On the image, the outline (the darkest point
across each run, as tools/measure_p1.py) is sampled along both runs and fitted with two lines, which give the image's
apex; the darkest point along the bisector from it gives d.  A corner is a WALL corner (convex wall, the fluid wrapping
around it: the inner corner of a turn, a wall's end, a junction's corner at the trunk) when its apex lies in the fluid,
and a FLUID corner (the outer corner of a turn) when its apex lies in the solid.
usage: python -I corner_radii.py RASTER DXF RECORD_JSON [panel key p1] OUT.json"""
import json
import math
import os
import sys

import numpy as np
from scipy.ndimage import map_coordinates
from shapely.geometry import Point
from shapely.ops import unary_union

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "rc7_traces", "code"))
import trace_lib as T                                         # noqa: E402

RAS, DXF, REC = sys.argv[1:4]
KEY = sys.argv[4] if len(sys.argv) > 5 else 'p1'
OUT = sys.argv[-1]
rec = json.load(open(REC))[KEY]
mp = rec['map_chamfer']
M = T.Map(mp['origin_px'], mp['ex_px_per_um'], mp['ey_px_per_um'])
A = np.column_stack([M.ex, M.ey])
k = float(np.hypot(*M.ex))                                    # px per µm
im = T.Raster(RAS).im
mask = T.Mask(DXF)
x0, x1, y0, y1 = rec['region_um']
from shapely.geometry import Polygon
chan = unary_union([Polygon(p) for L, p, *_ in mask.polys if L == 'CHANNELS' and p[:, 0].max() > x0 - 500 and p[:, 0].min() < x1 + 500])
isl = unary_union([Polygon(p) for L, p, *_ in mask.polys if L == 'ISLANDS' and p[:, 0].max() > x0 - 500 and p[:, 0].min() < x1 + 500])
fluid = chan.difference(isl)
c0, r0, c1, r1 = rec.get('micrograph_box', rec['crop'])


def to_px(p):
    return M(np.atleast_2d(p))[0]


def to_um(q):
    return np.linalg.solve(A, np.asarray(q, float) - M.origin)


def sample(q):
    q = np.atleast_2d(q)
    return map_coordinates(im, [q[:, 1], q[:, 0]], order=1)


def darkest(q0, n, half=5.0, step=0.05):
    t = np.arange(-half, half + step, step)
    v = sample(q0[None, :] + t[:, None] * n[None, :])
    i = int(np.argmin(v))
    if i in (0, len(t) - 1):
        return None
    return q0 + t[i] * n


def resample(P, ds=0.25):
    seg = np.hypot(*np.diff(P, axis=0).T)
    s = np.concatenate([[0], np.cumsum(seg)])
    ss = np.arange(0, s[-1], ds)
    return np.c_[np.interp(ss, s, P[:, 0]), np.interp(ss, s, P[:, 1])], ss


def line_fit(Q):
    c = Q.mean(axis=0)
    u, sv, vt = np.linalg.svd(Q - c)
    return c, vt[0]


def intersect(c1_, d1_, c2_, d2_):
    Mx = np.column_stack([d1_, -d2_])
    t = np.linalg.solve(Mx, c2_ - c1_)
    return c1_ + t[0] * d1_


L_RUN, GAP = 12.0, 22.0                                       # straight-run length and its offset from the corner (µm)
res = []
seen = []
for e in mask.edges(x0, x1, y0, y1):
    P, s = resample(np.asarray(e, float))
    if len(P) < 40:
        continue
    th = np.unwrap(np.arctan2(*np.diff(P, axis=0).T[::-1]))
    n = len(th)
    i_ = int(GAP / 0.25) + int(L_RUN / 0.25)
    for i in range(i_, n - i_, 4):
        a = th[i - int(GAP / 0.25) - int(L_RUN / 0.25):i - int(GAP / 0.25)]
        b = th[i + int(GAP / 0.25):i + int(GAP / 0.25) + int(L_RUN / 0.25)]
        if np.ptp(a) > 0.02 or np.ptp(b) > 0.02:
            continue
        dth = float(np.mean(b) - np.mean(a))
        if not (math.radians(75) < abs(dth) < math.radians(105)):
            continue
        Pa = P[i - int(GAP / 0.25) - int(L_RUN / 0.25):i - int(GAP / 0.25)]
        Pb = P[i + int(GAP / 0.25):i + int(GAP / 0.25) + int(L_RUN / 0.25)]
        ca, da = line_fit(Pa)
        cb, db = line_fit(Pb)
        apex = intersect(ca, da, cb, db)
        if any(np.hypot(*(apex - q)) < 8.0 for q in seen):
            continue
        # the mask's rounding: the polyline's nearest point to the apex
        dm = float(np.min(np.hypot(*(P - apex).T)))
        if dm > 25.0:
            continue
        seen.append(apex)
        r_mask = dm / (math.sqrt(2) - 1)
        kind = 'wall' if fluid.contains(Point(*apex)) else 'fluid'
        # the image: the outline across each run, two lines, the apex, the bisector
        pts = {}
        ok = True
        for nm, Q, d in (('a', Pa, da), ('b', Pb, db)):
            nrm = np.array([-d[1], d[0]])
            qs = []
            for q in Q[::8]:
                qp = to_px(q)
                npx = to_px(q + nrm) - qp
                npx /= np.hypot(*npx)
                if not (c0 + 8 < qp[0] < c1 - 8 and r0 + 8 < qp[1] < r1 - 8):
                    ok = False
                    break
                hit = darkest(qp, npx)
                if hit is not None:
                    qs.append(to_um(hit))
            if not ok or len(qs) < 4:
                ok = False
                break
            pts[nm] = np.array(qs)
        if not ok:
            continue
        ia, ida = line_fit(pts['a'])
        ib, idb = line_fit(pts['b'])
        apex_i = intersect(ia, ida, ib, idb)
        bis = (P[int(np.argmin(np.hypot(*(P - apex).T)))] - apex)
        bis /= np.hypot(*bis)
        q0 = to_px(apex_i)
        bpx = to_px(apex_i + bis) - q0
        bpx /= np.hypot(*bpx)
        t = np.arange(0.5, 25.0 * k, 0.05)
        v = sample(q0[None, :] + t[:, None] * bpx[None, :])
        j = int(np.argmin(v))
        di = float(t[j] / k)
        res.append(dict(apex_um=[float(apex[0]), float(apex[1])], kind=kind, turn_deg=float(math.degrees(dth)),
                        r_mask_um=r_mask, d_mask_um=dm, r_image_um=di / (math.sqrt(2) - 1), d_image_um=di,
                        image_apex_minus_mask_apex_um=float(np.hypot(*(apex_i - apex))),
                        edge_tag=e.tag))
res.sort(key=lambda r: (r['kind'], r['apex_um'][0]))
json.dump(dict(method=__doc__.split('\n\n')[0], panel=KEY, record=REC, dxf=DXF, corners=res), open(OUT, 'w'), indent=1)
for r in res:
    print('%-5s at (%7.1f, %7.1f) turn %+5.0f  r mask %5.1f  r image %5.1f  (apex shift %.1f um)  %s' % (
        r['kind'], r['apex_um'][0], r['apex_um'][1], r['turn_deg'], r['r_mask_um'], r['r_image_um'],
        r['image_apex_minus_mask_apex_um'], r['edge_tag']))
for kd in ('wall', 'fluid'):
    rr = [r for r in res if r['kind'] == kd]
    if rr:
        print('%s corners: n %d, r mask median %.1f, r image median %.1f (range %.1f-%.1f) um' % (
            kd, len(rr), np.median([r['r_mask_um'] for r in rr]), np.median([r['r_image_um'] for r in rr]),
            min(r['r_image_um'] for r in rr), max(r['r_image_um'] for r in rr)))
