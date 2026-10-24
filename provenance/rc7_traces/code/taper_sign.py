"""V13 -- the taper-direction SIGN test on a micrograph (RC8 trace code, not hashed; written for opencif, sharing no
code with review 1's validator).

In a CIF leg the side channels widen along the flow (Gifford 2014, Fig. 1 caption), so the sign of a leg's outer
wall-to-wall width slope along the leg is its flow direction.  Given a trace record (its feature map from mask
micrometres to image pixels, its micrograph box and its annotation boxes), the photomask (DXF, layer CHANNELS: the fluid
outline) and a leg's centre line Y, this measures, at positions X along the leg inside the micrograph:
  - the mask's outer walls (the line x = X cut with the fluid outline; the segment that contains Y);
  - the image's walls: for each mask wall, the darkest point of the image profile across the leg within +-window_um
    (the profile averaged over +-avg_um along the wall, sampled every 0.25 um), accepted only if it lies inside the
    window (not within 1 um of its ends) and is at least min_depth grey levels below the profile's median;
then fits the image's width and the mask's width against X (least squares) and compares their signs.

A position is ELIGIBLE when both walls' search segments (+-window, +-avg) lie inside the micrograph box and clear of the
annotation boxes; it is FOUND when both walls' minima are accepted; coverage = found / eligible (as the A7 coverage of
the outline metric).

Verdict (a sign test; the magnitude ratio is reported, never judged):
  PASS          the image slope has the mask's sign and |slope| >= 3 SE;
  FAIL          the opposite sign at >= 3 SE;
  INCONCLUSIVE  coverage below 0.5 (the window does not hold the walls: it must exceed the largest expected offset
                between the image's walls and the mask's -- review 3, E6: on the nested Fig. 4 reading the neighbour's
                walls lie ~18 um outside the mask's, up to about 20 um, so its window is +-30 um), fewer than 30 found
                positions, |slope| < 3 SE, or a mask slope below 1e-4 um/um (a mask that does not taper in the window
                has no sign to compare; RC8 code review, n2).
The SE is the slope's Newey-West (HAC) standard error, Bartlett kernel, lag 100 um of positions (20 at the 5 um step,
beyond two 45.8 um row pitches, whose period the residuals repeat), and never less than the independent-sample SE:
neighbouring positions share pixels (5 um apart on a 2.79 um/px raster, each averaged over +-3 um) and the residuals are
strongly autocorrelated (lag-1 r 0.64-0.77), so the independent-sample SE alone understates the uncertainty (on the
Fig. 4 neighbours the HAC SE is about 2 x the iid one; on leg 4 the iid SE is the larger) (RC8 code review, M1, N2); both are recorded (se_image_hac, se_image_iid).  The 30 positions are a floor on the fit's support, not 30 independent samples.
The DXF must be the photomask the trace registered (its sha256 against the record's mask_sha256; m4).

usage: python -I taper_sign.py RASTER TRACE_JSON DXF LEG LABEL OUT_JSON [window_um=30] [step_um=5] [avg_um=3]
       [x_min] [x_max]
  LEG: a leg's centre line Y (um), or 'neighbour' / 'leg4' (from the trace record's measured.neighbour_leg and
       mask_record.leg_Y);  x_min: a number, or '@array' = leg 4's array part only: X >= mask_record.last_row_x_um
       + 25 um, i.e. upstream of its last row, since leg 4 flows toward lower X on the Gifford Fig. 4 masks (the
       outlet leads downstream of the last row are not a taper and are left out)."""
import json, math, sys
import numpy as np
import ezdxf
from PIL import Image
from scipy import ndimage
from shapely.geometry import Polygon, LineString

MIN_FOUND = 30
MIN_COVERAGE = 0.5
MIN_DEPTH = 6.0          # grey levels below the profile median
HAC_LAG_UM = 100.0       # the Newey-West lag, in um along the leg (20 positions at the 5 um step): beyond two of the
                         # inset's 45.8 um row pitches, whose period the residuals repeat (RC8 code review, N2)
MIN_MASK_SLOPE = 1e-4    # um/um: below it the mask has no taper direction in the window


def hac_slope_se(X, resid, lag):
    """the least-squares slope's Newey-West standard error (Bartlett kernel): var = sum_l w_l sum_i u_i u_{i-l} /
    Sxx^2, u_i = (X_i - mean X) resid_i, w_l = 1 - l / (lag + 1), both directions"""
    x = X - X.mean()
    u = x * resid
    S = float(u @ u)
    for l_ in range(1, min(lag, len(u) - 1) + 1):
        S += 2.0 * (1.0 - l_ / (lag + 1.0)) * float(u[l_:] @ u[:-l_])
    n = len(u)
    S *= n / max(n - 2, 1)                       # the small-sample factor of the iid formula
    return math.sqrt(max(S, 0.0)) / float(x @ x)


def load_outline(dxf):
    doc = ezdxf.readfile(dxf)
    polys = [Polygon([(float(p[0]), float(p[1])) for p in e.get_points()])
             for e in doc.modelspace() if e.dxf.layer == 'CHANNELS']
    return max(polys, key=lambda q: q.area)


def mask_walls(poly, X, Y):
    """the outer walls (y_low, y_high) of the leg whose fluid contains (X, Y), or None"""
    b = poly.bounds
    s = poly.intersection(LineString([(X, b[1] - 10), (X, b[3] + 10)]))
    parts = [s] if s.geom_type == 'LineString' else list(getattr(s, 'geoms', []))
    for g in parts:
        ys = [c[1] for c in g.coords]
        if min(ys) - 1e-6 <= Y <= max(ys) + 1e-6:
            return min(ys), max(ys)
    return None


def main():
    import hashlib
    ras, trj, dxf, legY, label, out = sys.argv[1:7]
    window = float(sys.argv[7]) if len(sys.argv) > 7 else 30.0
    step = float(sys.argv[8]) if len(sys.argv) > 8 else 5.0
    avg = float(sys.argv[9]) if len(sys.argv) > 9 else 3.0
    x_min_arg = sys.argv[10] if len(sys.argv) > 10 else -1e18
    x_max = float(sys.argv[11]) if len(sys.argv) > 11 else 1e18
    g = np.asarray(Image.open(ras).convert('L'), float)
    t = json.load(open(trj))
    h = hashlib.sha256(open(dxf, 'rb').read()).hexdigest()
    want = t.get('mask_sha256') or (t.get('mask_record') or {}).get('dxf_sha256')
    if want and h != want:
        raise SystemExit('%s: sha256 %s is not the traced photomask %s (the trace record\'s mask_sha256)' % (
            dxf, h[:16], want[:16]))
    mr = t.get('mask_record') or {}
    if legY == 'neighbour':
        legY = float(mr['leg_Y'][int(t['measured']['neighbour_leg']) - 1])
    elif legY == 'leg4':
        legY = float(mr['leg_Y'][3])
    else:
        legY = float(legY)
    x_min = (float(mr['last_row_x_um']) + 25.0) if x_min_arg == '@array' else float(x_min_arg)
    m = t['map_feature']
    o = np.array(m['origin_px']); ex = np.array(m['ex_px_per_um']); ey = np.array(m['ey_px_per_um'])
    x0, y0, x1, y1 = t['micrograph_box']
    excl = t.get('exclude_px') or []
    poly = load_outline(dxf)

    def px(X, Y):
        return o + ex * X + ey * Y

    def ok_px(p, margin=6):
        if not (x0 + margin <= p[0] <= x1 - margin and y0 + margin <= p[1] <= y1 - margin):
            return False
        return not any(a <= p[0] <= c and b <= p[1] <= d for a, b, c, d in excl)

    # X range: every X whose leg centre point lies inside the micrograph box
    Xs = []
    lo, hi = poly.bounds[0], poly.bounds[2]
    for X in np.arange(math.floor(lo), math.ceil(hi), step):
        if x_min <= X <= x_max and ok_px(px(X, legY), margin=0):
            Xs.append(float(X))
    s_grid = np.arange(-window, window + 1e-9, 0.25)
    a_grid = np.arange(-avg, avg + 1e-9, 0.5)
    rows = []
    eligible = 0
    for X in Xs:
        w = mask_walls(poly, X, legY)
        if w is None:
            continue
        corners = [px(X + a, yw + s) for yw in w for a in (-avg, avg) for s in (-window, window)]
        if not all(ok_px(c) for c in corners):
            continue
        eligible += 1
        found = []
        for yw in w:
            P = np.array([[px(X + a, yw + s) for s in s_grid] for a in a_grid])      # (along, across, 2)
            if not all(ok_px(P[i, j]) for i in (0, -1) for j in (0, -1)):
                found.append(None)
                continue
            vals = ndimage.map_coordinates(g, [P[..., 1].ravel(), P[..., 0].ravel()], order=1).reshape(P.shape[:2])
            prof = vals.mean(axis=0)
            k = int(np.argmin(prof))
            if k < 4 or k > len(prof) - 5 or np.median(prof) - prof[k] < MIN_DEPTH:
                found.append(None)
                continue
            # sub-sample minimum (parabola through the three samples)
            y_m1, y_0, y_p1 = prof[k - 1], prof[k], prof[k + 1]
            den = y_m1 - 2 * y_0 + y_p1
            dk = 0.5 * (y_m1 - y_p1) / den if den > 0 else 0.0
            found.append(float(s_grid[k] + dk * 0.25))
        if None in found:
            continue
        rows.append(dict(X=X, mask_low=w[0], mask_high=w[1], off_low=found[0], off_high=found[1],
                         width_mask=w[1] - w[0], width_image=(w[1] + found[1]) - (w[0] + found[0])))
    n = len(rows)
    res = dict(label=label, raster=ras.split('/')[-1], dxf=dxf.split('/')[-1], dxf_sha256=h, leg_Y=legY, window_um=window, step_um=step,
               avg_um=avg, x_limits=[x_min if x_min > -1e17 else None, x_max if x_max < 1e17 else None], positions_tried=len(Xs), eligible=eligible, positions=n,
               coverage=(n / eligible) if eligible else None, min_found=MIN_FOUND, min_coverage=MIN_COVERAGE)
    if n >= 3:
        X = np.array([r['X'] for r in rows]); wi = np.array([r['width_image'] for r in rows])
        wm = np.array([r['width_mask'] for r in rows])
        A = np.vstack([X - X.mean(), np.ones(n)]).T
        ci, rs, _, _ = np.linalg.lstsq(A, wi, rcond=None)
        resid = wi - A @ ci
        se_iid = math.sqrt((resid @ resid) / max(n - 2, 1) / ((X - X.mean()) @ (X - X.mean())))
        lag = max(1, int(round(HAC_LAG_UM / step)))
        se_hac = hac_slope_se(X, resid, lag)
        se = max(se_hac, se_iid)                 # never below the independent-sample SE (N2: conservative)
        r1 = float(np.corrcoef(resid[:-1], resid[1:])[0, 1]) if n > 3 else None
        cm = np.linalg.lstsq(A, wm, rcond=None)[0]
        z = ci[0] / se if se > 0 else float('inf')
        same = math.copysign(1, ci[0]) == math.copysign(1, cm[0])
        if n < MIN_FOUND or (n / eligible) < MIN_COVERAGE or abs(z) < 3 or abs(cm[0]) < MIN_MASK_SLOPE:
            verdict = 'INCONCLUSIVE'
        else:
            verdict = 'PASS' if same else 'FAIL'
        res.update(X_range=[float(X.min()), float(X.max())], slope_image=float(ci[0]), se_image=float(se),
                   se_method='max(Newey-West (Bartlett, lag %d positions, %g um), the independent-sample SE)' % (
                       lag, lag * step), se_image_hac=float(se_hac),
                   se_image_iid=float(se_iid), z_image_iid=float(ci[0] / se_iid) if se_iid > 0 else None,
                   residual_lag1_r=r1, min_mask_slope=MIN_MASK_SLOPE,
                   z_image=float(z), slope_mask=float(cm[0]), magnitude_ratio=float(ci[0] / cm[0]) if cm[0] else None,
                   offset_low_mean_um=float(np.mean([r['off_low'] for r in rows])),
                   offset_low_sd_um=float(np.std([r['off_low'] for r in rows])),
                   offset_high_mean_um=float(np.mean([r['off_high'] for r in rows])),
                   offset_high_sd_um=float(np.std([r['off_high'] for r in rows])),
                   width_image_mean_um=float(wi.mean()), width_mask_mean_um=float(wm.mean()), verdict=verdict,
                   test='sign test: PASS when the image slope has the mask slope\'s sign beyond 3 SE (the larger of the Newey-West SE and the independent-sample SE), '
                        'with coverage >= %g, >= %d positions and a mask slope >= %g um/um; the magnitude ratio is '
                        'reported, not judged' % (MIN_COVERAGE, MIN_FOUND, MIN_MASK_SLOPE))
    else:
        res.update(verdict='INCONCLUSIVE')
    res['rows'] = rows
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != 'rows'}))


if __name__ == '__main__':
    main()
