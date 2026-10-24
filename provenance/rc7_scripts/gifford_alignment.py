# rc7 (code review R11): shipped with the tree; copied from the rc6 Gifford Fig. 6 package (devpkg/opencif_rc6_gifford2014_fig6_v2/gifford_alignment.py), a library its measure step imports; rc7: fig4_profiles takes the leg beside leg 4 from the record (path order: leg 3).
"""gifford_alignment.py -- the two Gifford packages' measurements beyond rc6's trace: the start width a micrograph implies
(the frozen physics' row-fraction law), and the checks that the aligned reproduction lies on the micrograph.

The start width.  Gifford's design (constant w_c, constant f, App. A eq. (g4) = his eqn (4)) gives the side channel
after gap i as r(w_s(i)) = 2 r(w_c) q^i / (1 - q^i), q = 1 - 2f (rc6's exact closure agrees to 0.1 um).  So the row
fraction at which the schedule reaches a width w is

    phi(w) = ln Q(w) / ln Q(w_end),    Q(w) = r(w) / (2 r(w_c) + r(w))

(f, N and the gap G drop out; r is rc6's duct law, cifsched.rr).  A start width w0 > 0 (spec v5 section 8, DIRECT mode;
App. A: "the design code's other start mode"; the glossary's sigma0, a supplied input) enters the same schedule at the
fraction delta with Q(w0) = Q(w_end)^delta, so a place at the chip's row fraction phi_c shows the width of the fraction
delta + (1 - delta) phi_c.  The sealed reference build (rc6's design preset) gives phi_c at each place a micrograph shows,
from its drawn width there: phi_c = phi(w_drawn) at the design's own final width.  The micrograph's width there gives
delta, and delta the start width.

Nothing here changes the rc6 code; it reads rc6's duct law and the rc6 trace code's records and raster tools."""
import math

import numpy as np


def _rr():
    from opencif.physics import schedule as S
    return S.rr


def Q(w, wc, h):
    rr = _rr()
    r = rr(w, h)
    return r / (2.0 * rr(wc, h) + r)


def phi(w, w_end, wc, h):
    return math.log(Q(w, wc, h)) / math.log(Q(w_end, wc, h))


def width_at(ph, w_end, wc, h):
    from scipy.optimize import brentq
    qq = Q(w_end, wc, h) ** ph
    return brentq(lambda w: Q(w, wc, h) - qq, 1e-3, 5.0 * w_end)


def start_from(phi_chip, w_img, w_end, wc, h):
    """(delta, start width um, App. A's mass defect 2 sigma0) from one place: its chip row fraction and the image's
    width there, at the reproduction's final width w_end"""
    d = (phi(w_img, w_end, wc, h) - phi_chip) / (1.0 - phi_chip)
    q0 = Q(w_end, wc, h) ** d
    from scipy.optimize import brentq
    w0 = brentq(lambda w: Q(w, wc, h) - q0, 1e-3, 5.0 * w_end)
    return d, w0, (1.0 - q0) / q0


def predict(phi_chip, delta, w_end, wc, h):
    """the width the start offset delta puts at the chip row fraction phi_chip"""
    return width_at(delta + (1.0 - delta) * phi_chip, w_end, wc, h)


# --------------------------------------------------------------------------- Fig. 4: leg 2 and leg 4 read on profiles
def fig4_profiles(T, im, rec, mask=None, neighbour=None):
    """the inset's walls and post rows across leg 4 and the leg beside it (leg 2 with nested legs; rc7: leg 3 with the
    legs in path order, the record's measured.neighbour_leg, or `neighbour`), as dark-line centres on profiles across
    the legs every 20 um along them, at the registration rec (the method panels_gifford.py uses for leg 4's widths);
    with a mask, the mask's own walls on the same profiles (the drawn side channel = |wall| - (w_c/2 + l_g)).
    Returns {leg: {'image': {name: (median um, n)}, 'sides_image': (lo, hi), 'sides_mask': (lo, hi) or None}}"""
    m = rec['map_chamfer']
    O, ex, ey = (np.array(m[k], float) for k in ('origin_px', 'ex_px_per_um', 'ey_px_per_um'))
    M = np.c_[ex, ey]
    k = float(np.sqrt(abs(np.linalg.det(M))))
    fmap = lambda xy: O + (M @ np.asarray(xy, float).T).T
    legY = rec['mask_record']['leg_Y']
    vimg = ey / np.linalg.norm(ey)
    x0, x1 = rec['region_um'][0], rec['region_um'][1]
    out = {}
    nb = int(neighbour or (rec.get('measured') or {}).get('neighbour_leg') or 2)
    for Y, leg, half, wall in ((legY[nb - 1], nb, 300.0, 214.0), (legY[3], 4, 480.0, 409.0)):
        names = (('wall_lo', -wall), ('post_out_lo', -109.5), ('post_in_lo', -50.0), ('post_in_hi', 50.0),
                 ('post_out_hi', 109.5), ('wall_hi', wall))
        acc, xs, mw = {}, [], {'lo': [], 'hi': []}
        for x in np.arange(x0 + 100, x1 - 100, 20.0):
            p0 = fmap([[x, Y]])[0]
            if not (205 < p0[0] < 820 and 60 < p0[1] < 545):
                continue
            s, pr = T.profile_line(im, p0, vimg, half * k, step=0.25, width=1.0)
            L = T.dark_lines(pr, 10, halfwin=12, off=0) * 0.25 / k - half
            xs.append(float(x))
            for n, c in names:
                v = T.line_near(np.asarray(L), c, 25.0 if 'wall' in n else 6.0)
                if v is not None:
                    acc.setdefault(n, []).append(v)
            if mask is not None:
                for side, sgn in (('lo', -1), ('hi', 1)):
                    yw = _mask_wall(mask, x, Y + sgn * wall, 40.0)
                    if yw is not None:
                        mw[side].append(abs(yw - Y) - 109.5)
        a = {n: (float(np.median(v)), len(v)) for n, v in acc.items()}
        sides = (a['post_out_lo'][0] - a['wall_lo'][0], a['wall_hi'][0] - a['post_out_hi'][0])
        out[leg] = dict(image=a, profiles=len(xs), x_um=[min(xs), max(xs)], sides_image=sides,
                        sides_mask=(float(np.mean(mw['lo'])), float(np.mean(mw['hi']))) if mask is not None and mw['lo'] and mw['hi'] else None)
    return out


def _mask_wall(mask, x, y_near, tol):
    """the y where the vertical line at x crosses the mask's edges nearest y_near (within tol)"""
    from shapely.geometry import LineString
    ln = LineString([(x, y_near - tol), (x, y_near + tol)])
    ys = []
    for e in mask.edges(x - 1.0, x + 1.0, y_near - tol, y_near + tol):
        g = LineString(np.asarray(e)).intersection(ln)
        for part in (g.geoms if hasattr(g, 'geoms') else [g]):
            if not part.is_empty and part.geom_type == 'Point':
                ys.append(part.y)
    return min(ys, key=lambda y: abs(y - y_near)) if ys else None


# --------------------------------------------------------------------------- the outline offsets at narrow windows
def narrow_window(T, raster, mask, rec, post_max, claim, windows=(8.0, 3.5)):
    """the claim's nearest-outline offsets at the rc6 window and at a window below half the gap (register GIF4-2):
    the posts' flanks along the flow (signed median, mean, p90 |d|), and the whole claim (median, p90 |d|).
    claim(edge) -> bool selects the claim's edges.  Returns {window_px: {...}}"""
    m = rec['map_chamfer']
    O, ex, ey = (np.array(m[k], float) for k in ('origin_px', 'ex_px_per_um', 'ey_px_per_um'))
    M = np.c_[ex, ey]
    k = float(np.sqrt(abs(np.linalg.det(M))))
    det = float(np.linalg.det(M))
    crop = rec.get('micrograph_box') or rec['crop']
    x0, x1, y0, y1 = rec['region_um']
    along = np.array(rec['along_px'], float)
    along = along / np.linalg.norm(along)
    excl = [tuple(v) for v in rec.get('exclude_px', [])]
    ps, pn, cs, cn = [], [], [], []
    for e in mask.edges(x0, x1, y0, y1, post_max=post_max):
        if not claim(e):
            continue
        seg = O + (M @ np.asarray(e, float).T).T
        sg = int(e.orient * np.sign(det))
        cs.append(seg); cn.append(sg)
        if e.tag == 'post':
            ps.append(seg); pn.append(sg)
    out = {}
    for spx in windows:
        r = {}
        P, N = T.sample_edges(ps, 1.0, pn)
        off, _ = T.nearest_outline(raster.im, P, N, k, spx)
        off[T.in_rects(P, excl, pad=spx)] = np.nan
        off[~T.in_rects(P, [crop], pad=-spx)] = np.nan
        f = off[(np.abs(N @ along) >= 0.7071) & np.isfinite(off)]
        r['post_flanks_along'] = dict(n=int(len(f)), median_um=float(np.median(f)), mean_um=float(np.mean(f)),
                                      p90_abs_um=float(np.percentile(np.abs(f), 90)))
        P, N = T.sample_edges(cs, 0.5, cn)
        off, _ = T.nearest_outline(raster.im, P, N, k, spx)
        off[T.in_rects(P, excl, pad=spx)] = np.nan
        off[~T.in_rects(P, [crop], pad=-spx)] = np.nan
        f = np.abs(off[np.isfinite(off)])
        r['claim'] = dict(n=int(len(f)), of=int(len(off)), median_abs_um=float(np.median(f)),
                          p90_abs_um=float(np.percentile(f, 90)))
        r['window_um'] = spx / k
        out['%.1f' % spx] = r
    return out
