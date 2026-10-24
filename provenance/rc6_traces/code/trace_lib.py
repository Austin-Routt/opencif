"""trace_lib.py -- the shared machinery of the micrograph traces (provenance/rc6_traces).

A trace places the edges of a shipped opencif DXF mask on a paper's micrograph.  The map from mask coordinates (µm)
to image pixels is an affine map  img = origin + x * ex + y * ey  whose scale is the panel's printed scale bar and
whose rotation is the layout's (or, for a rotated inset, the direction measured on the image), never fitted to the
mask.  Only the translation (origin) is fitted, and it is fitted twice, by two independent methods that must agree:

  * chamfer registration: the mask's edges, sampled every half pixel, are scored against the distance transform of
    the image's dark outlines (truncated mean distance, px); the translation is the minimum of that score over a
    declared search window (grid search, then a parabolic refinement);
  * feature registration: the same features (post centres, wall lines) measured on the image as dark-line centres
    and on the mask from its polygons, matched by a declared anchor rule, least squares.

Validators computed for every panel (trace_results.json):
  V1  the two registrations agree (|difference| <= 1.5 px in each direction);
  V2  the fitted score is far below the null: the median score over the whole search window, and the score one
      post pitch away along the flow (the periodic ambiguity), are reported with their ratios to the best score;
  V3  per-class signed residuals (post / block / island / channel x across / along the flow) of the nearest-outline
      metric, with a flag when a class is offset by more than 1 px on average;
  V4  the residual map (every trace segment coloured by its signed offset).
The nearest-outline metric (median |d|, 90th percentile, signed mean) is reported as before (search +-8 px along the
local normal), together with its value at the null registrations, so that a small median can be read against what
the same metric gives for a trace that is NOT registered.
"""
import hashlib
import json
import math
import os
import textwrap

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.ndimage import map_coordinates

TRACE_COLOUR = '#00d5ff'
NOCLAIM_COLOUR = '#ffb000'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------------------ mask
class Mask:
    """the closed polygons of a shipped DXF (layers CHANNELS: fluid; ISLANDS: posts, blocks and walls in the fluid)"""

    def __init__(self, path):
        import ezdxf
        self.path, self.sha256 = path, sha256(path)
        doc = ezdxf.readfile(path)
        self.polys = []
        self.orient = []
        for e in doc.modelspace():
            if e.dxftype() != 'LWPOLYLINE':
                continue
            p = np.array([(q[0], q[1]) for q in e.get_points()], float)
            if len(p) < 3:
                continue
            area2 = float(np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1]))
            self.polys.append((e.dxf.layer, p, p[:, 0].min(), p[:, 0].max(), p[:, 1].min(), p[:, 1].max()))
            self.orient.append(-1 if area2 > 0 else 1)     # +1: the left normal points out of the polygon (clockwise)
        isl = [(x0, x1, y0, y1) for L, p, x0, x1, y0, y1 in self.polys if L == 'ISLANDS']
        self.islands = np.array(isl) if isl else np.zeros((0, 4))

    @staticmethod
    def tag(layer, x0, x1, y0, y1, post_max=(70.0, 130.0)):
        if layer != 'ISLANDS':
            return 'channel'
        lx, ly = x1 - x0, y1 - y0
        if lx <= post_max[0] and ly <= post_max[1]:
            return 'post'
        if lx <= 260 and ly <= 260:
            return 'block'
        return 'island'

    def edges(self, x0, x1, y0, y1, post_max=(70.0, 130.0)):
        """edge polylines clipped to the box (mask µm), each an (N, 2) array with a .tag attribute in a list"""
        from shapely.geometry import LineString, box
        clip = box(x0, y0, x1, y1)
        out = []
        for (L, p, bx0, bx1, by0, by1), o in zip(self.polys, self.orient):
            if bx1 < x0 or bx0 > x1 or by1 < y0 or by0 > y1:
                continue
            g = LineString(np.vstack([p, p[:1]])).intersection(clip)
            tg = self.tag(L, bx0, bx1, by0, by1, post_max)
            for part in (g.geoms if hasattr(g, 'geoms') else [g]):
                if not part.is_empty and part.geom_type == 'LineString':
                    out.append(Edge(np.array(part.coords), tg, o))
        return out

    def posts(self, x0, x1, y0, y1, max_len=70.0, max_wid=130.0):
        """centres (x, y) of the islands of post size whose bbox centre lies in the box"""
        I = self.islands
        if not len(I):
            return np.zeros((0, 2))
        cx, cy = (I[:, 0] + I[:, 1]) / 2, (I[:, 2] + I[:, 3]) / 2
        s = (cx >= x0) & (cx <= x1) & (cy >= y0) & (cy <= y1) & ((I[:, 1] - I[:, 0]) <= max_len) & ((I[:, 3] - I[:, 2]) <= max_wid)
        c = np.c_[cx[s], cy[s]]
        return c[np.argsort(c[:, 0])]


class Edge(np.ndarray):
    """a mask edge polyline (µm) with its class tag and orientation (+1: its left normal points out of its polygon)"""
    def __new__(cls, a, tag='?', orient=1):
        o = np.asarray(a, float).view(cls)
        o.tag = tag
        o.orient = orient
        return o

    def __array_finalize__(self, obj):
        self.tag = getattr(obj, 'tag', '?')
        self.orient = getattr(obj, 'orient', 1)


# ------------------------------------------------------------------------------------------------------- raster
class Raster:
    def __init__(self, path, expect_sha=None):
        self.path, self.sha256 = path, sha256(path)
        if expect_sha and not self.sha256.startswith(expect_sha):
            raise ValueError('%s: sha256 %s is not the expected %s' % (path, self.sha256[:16], expect_sha))
        self.im = np.asarray(Image.open(path).convert('L')).astype(float)

    def outline_distance(self, crop, sigma_small=0.7, sigma_large=4.0, thr=None, exclude=(), min_component=0):
        """distance (px) from every pixel of the crop to the nearest dark-outline pixel.  Outline pixels are dark
        ridges: darker than their neighbourhood (difference of Gaussians above a threshold) and a local minimum across
        x or y.  `exclude`: rectangles (c0, r0, c1, r1) of caption ink, arrows and labels removed from the outline set."""
        c0, r0, c1, r1 = crop
        a = self.im[r0:r1, c0:c1]
        d = ndi.gaussian_filter(a, sigma_large) - ndi.gaussian_filter(a, sigma_small)
        if thr is None:
            thr = max(4.0, 1.5 * float(np.median(np.abs(d - np.median(d)))) * 1.4826 * 2)
        s = ndi.gaussian_filter(a, sigma_small)
        minx = (s <= np.roll(s, 1, 1)) & (s <= np.roll(s, -1, 1))
        miny = (s <= np.roll(s, 1, 0)) & (s <= np.roll(s, -1, 0))
        ridge = (d > thr) & (minx | miny)
        for (e0, f0, e1, f1) in exclude:
            ridge[max(f0 - r0, 0):max(f1 - r0, 0), max(e0 - c0, 0):max(e1 - c0, 0)] = False
        if min_component:            # isolated dark specks (beads, dust) are not outlines
            lab, n = ndi.label(ridge, structure=np.ones((3, 3)))
            if n:
                size = ndi.sum(ridge, lab, index=np.arange(1, n + 1))
                ridge &= np.r_[False, size >= min_component][lab]
        dt = ndi.distance_transform_edt(~ridge)
        return dt, ridge, thr


# ------------------------------------------------------------------------------------------------------- map
class Map:
    """image px = origin + x * ex + y * ey (x, y in mask µm)"""

    def __init__(self, origin, ex, ey):
        self.origin = np.asarray(origin, float)
        self.ex = np.asarray(ex, float)
        self.ey = np.asarray(ey, float)

    def __call__(self, pts):
        pts = np.asarray(pts, float)
        return self.origin + np.outer(pts[:, 0], self.ex) + np.outer(pts[:, 1], self.ey)

    def shifted(self, d):
        return Map(self.origin + np.asarray(d, float), self.ex, self.ey)

    def inverse(self, px):
        A = np.c_[self.ex, self.ey]
        return np.linalg.solve(A, (np.asarray(px, float) - self.origin).T).T

    def to_dict(self):
        return dict(origin_px=self.origin.tolist(), ex_px_per_um=self.ex.tolist(), ey_px_per_um=self.ey.tolist())


def vertical_map(k, origin=(0.0, 0.0)):
    """flow (mask +x) down the image, mask +y to the right (a non-mirrored placement)"""
    return Map(origin, (0.0, k), (k, 0.0))


def rotated_map(k, theta_rad, origin=(0.0, 0.0)):
    """mask +x along (cos t, sin t) in image px (y down), mask +y along (sin t, -cos t): non-mirrored"""
    return Map(origin, (k * math.cos(theta_rad), k * math.sin(theta_rad)), (k * math.sin(theta_rad), -k * math.cos(theta_rad)))


# ------------------------------------------------------------------------------------------------ sampling
def sample_edges(segs_px, step=0.5, signs=None):
    """points every `step` px along each polyline, and their unit normals (the left normal times the polyline's sign:
    pass signs so that every normal points out of its polygon)"""
    P, N = [], []
    signs = [1] * len(segs_px) if signs is None else signs
    for s, sg in zip(segs_px, signs):
        s = np.asarray(s, float)
        if len(s) < 2:
            continue
        d = np.diff(s, axis=0)
        L = np.hypot(d[:, 0], d[:, 1])
        ok = L > 1e-9
        d, L, s0 = d[ok], L[ok], s[:-1][ok]
        if not len(L):
            continue
        cum = np.r_[0, np.cumsum(L)]
        t = np.arange(step / 2, cum[-1], step)
        i = np.clip(np.searchsorted(cum, t) - 1, 0, len(L) - 1)
        f = (t - cum[i]) / L[i]
        P.append(s0[i] + f[:, None] * d[i])
        N.append(sg * np.c_[-d[i, 1], d[i, 0]] / L[i][:, None])
    if not P:
        return np.zeros((0, 2)), np.zeros((0, 2))
    return np.vstack(P), np.vstack(N)


def chamfer(dt, crop, pts, trunc=6.0):
    """truncated mean distance (px) from the trace points to the image's outlines; points outside the crop ignored"""
    c0, r0, c1, r1 = crop
    x, y = pts[:, 0] - c0, pts[:, 1] - r0
    ok = (x >= 0) & (y >= 0) & (x <= c1 - c0 - 1) & (y <= r1 - r0 - 1)
    if ok.sum() < 20:
        return float('nan'), int(ok.sum())
    v = map_coordinates(dt, [y[ok], x[ok]], order=1, mode='nearest')
    return float(np.mean(np.minimum(v, trunc))), int(ok.sum())


def register_chamfer(dt, crop, pts_px0, search, step=0.5, trunc=6.0):
    """grid search of the translation (dx, dy) in px over search = (dx0, dx1, dy0, dy1) minimising the chamfer score of
    pts_px0 + (dx, dy); returns (best (dx, dy), best score, the landscape (dxs, dys, S))"""
    dxs = np.arange(search[0], search[1] + 1e-9, step)
    dys = np.arange(search[2], search[3] + 1e-9, step)
    S = np.full((len(dys), len(dxs)), np.nan)
    for j, dy in enumerate(dys):
        for i, dx in enumerate(dxs):
            S[j, i] = chamfer(dt, crop, pts_px0 + (dx, dy), trunc)[0]
    j, i = np.unravel_index(np.nanargmin(S), S.shape)

    def para(v, a, b, c):                       # parabolic refinement along one axis
        den = a - 2 * b + c
        return v + (0.5 * (a - c) / den) * step if den > 0 else v
    bx = para(dxs[i], S[j, i - 1], S[j, i], S[j, i + 1]) if 0 < i < len(dxs) - 1 else dxs[i]
    by = para(dys[j], S[j - 1, i], S[j, i], S[j + 1, i]) if 0 < j < len(dys) - 1 else dys[j]
    best = chamfer(dt, crop, pts_px0 + (bx, by), trunc)[0]
    return (float(bx), float(by)), float(best), (dxs, dys, S)


# ------------------------------------------------------------------------------------------------ profiles
def dark_lines(prof, min_depth=8.0, halfwin=6, off=0.0):
    """centres of the dark lines of a 1-D profile: local minima deeper than min_depth relative to the lower of the two
    maxima within +-halfwin samples, refined by a parabola through the three samples"""
    p = np.asarray(prof, float)
    out = []
    for i in range(1, len(p) - 1):
        if p[i] < p[i - 1] and p[i] <= p[i + 1]:
            lo, hi = max(0, i - halfwin), min(len(p), i + halfwin + 1)
            depth = min(p[lo:i].max(), p[i + 1:hi].max()) - p[i]
            if depth >= min_depth:
                den = p[i - 1] - 2 * p[i] + p[i + 1]
                d = 0.5 * (p[i - 1] - p[i + 1]) / den if den != 0 else 0.0
                out.append(off + i + d)
    return np.array(out)


def row_profile(im, row, c0, c1, half=1):
    """the profile along a row (x), averaged over 2*half+1 rows"""
    r = int(round(row))
    return im[r - half:r + half + 1, c0:c1].mean(axis=0)


def col_profile(im, col, r0, r1, half=2):
    c = int(round(col))
    return im[r0:r1, c - half:c + half + 1].mean(axis=1)


def pair_centres(lines, lo, hi, skip=1):
    """line pairs whose separation lies in [lo, hi]: their midpoints (a post between its two outlines).  Up to `skip`
    lines between the two outlines are passed over (an annotation arrow drawn across a post)."""
    c, i = [], 0
    v = list(lines)
    while i < len(v) - 1:
        for j in range(i + 1, min(len(v), i + 2 + skip)):
            if lo <= v[j] - v[i] <= hi:
                c.append((v[i] + v[j]) / 2.0)
                i = j + 1
                break
        else:
            i += 1
    return np.array(c)


def in_rects(P, rects, pad=0.0):
    """True for the points inside any of the rectangles (c0, r0, c1, r1), grown by pad px"""
    m = np.zeros(len(P), bool)
    for (a, b, c, d) in rects:
        m |= (P[:, 0] >= a - pad) & (P[:, 0] <= c + pad) & (P[:, 1] >= b - pad) & (P[:, 1] <= d + pad)
    return m


def line_near(lines, x, tol):
    if not len(lines):
        return None
    j = int(np.argmin(np.abs(lines - x)))
    return float(lines[j]) if abs(lines[j] - x) <= tol else None


# ------------------------------------------------------------------------------------------------ metric
def nearest_outline(im, pts, nrm, k, search_px=8.0, min_contrast=18.0, ink=40.0):
    """for each trace point: the signed distance (µm) along its normal to the nearest dark minimum within +-search_px
    (parabolic), or NaN when the profile has no outline (contrast below min_contrast, or caption ink darker than ink,
    or the minimum at the window's end).  The same definition as the earlier rounds' metric."""
    H, W = im.shape
    s = np.arange(-search_px, search_px + 0.25, 0.5)
    xs = pts[:, 0:1] + s[None, :] * nrm[:, 0:1]
    ys = pts[:, 1:2] + s[None, :] * nrm[:, 1:2]
    inside = (xs.min(1) >= 1) & (ys.min(1) >= 1) & (xs.max(1) <= W - 2) & (ys.max(1) <= H - 2)
    out = np.full(len(pts), np.nan)
    if not inside.any():
        return out, inside
    prof = map_coordinates(im, [ys[inside].ravel(), xs[inside].ravel()], order=1, mode='nearest').reshape(-1, len(s))
    j = np.argmin(prof, axis=1)
    ok = (prof.max(1) - prof.min(1) >= min_contrast) & (prof.min(1) >= ink) & (j > 0) & (j < len(s) - 1)
    jj = np.clip(j, 1, len(s) - 2)
    a, b, c = prof[np.arange(len(j)), jj - 1], prof[np.arange(len(j)), jj], prof[np.arange(len(j)), jj + 1]
    den = a - 2 * b + c
    dj = np.where(den != 0, 0.5 * (a - c) / np.where(den != 0, den, 1), 0.0)
    val = (s[jj] + dj * 0.5) / k
    idx = np.where(inside)[0]
    out[idx[ok]] = val[ok]
    return out, inside


def summarise(d):
    d = d[np.isfinite(d)]
    if not len(d):
        return dict(n=0)
    a = np.abs(d)
    return dict(n=int(len(d)), median_abs_um=float(np.median(a)), p90_abs_um=float(np.percentile(a, 90)),
                mean_signed_um=float(np.mean(d)))


def by_class(edges_px, tags, im, k, crop, along_dir, search_px=8.0, exclude=(), signs=None):
    """V3: the signed nearest-outline offsets (µm; + = the image's outline lies outside the mask's polygon: an island
    bigger, a channel wider) per class: tag x the offset's direction (across / along the flow) x the third of the crop
    (left / centre / right).  Classified point by point, so a closed post outline contributes to both directions."""
    c0, r0, c1, r1 = crop
    signs = [1] * len(edges_px) if signs is None else signs
    along_dir = np.asarray(along_dir, float) / np.hypot(*along_dir)
    res = {}
    for seg, tag, sg in zip(edges_px, tags, signs):
        P, Nn = sample_edges([seg], 1.0, [sg])
        if not len(P):
            continue
        off, _ = nearest_outline(im, P, Nn, k, search_px)
        off[in_rects(P, exclude, pad=search_px)] = np.nan
        off[~in_rects(P, [crop], pad=-search_px)] = np.nan
        across = np.abs(Nn @ along_dir) < 0.7071
        third = np.clip((3 * (P[:, 0] - c0) / max(c1 - c0, 1)).astype(int), 0, 2)
        for a in (True, False):
            for th in range(3):
                sel = (across == a) & (third == th)
                if sel.any():
                    res.setdefault('%s_%s_%s' % (tag, 'across' if a else 'along', 'LCR'[th]), []).extend(off[sel].tolist())
    out = {}
    for cls, v in sorted(res.items()):
        v = np.array(v)
        f = v[np.isfinite(v)]
        out[cls] = dict(n=int(len(f)), total=int(len(v)), mean_signed_um=float(f.mean()) if len(f) else None,
                        sd_um=float(f.std()) if len(f) else None, median_abs_um=float(np.median(np.abs(f))) if len(f) else None,
                        systematic=bool(len(f) >= 20 and abs(f.mean()) > 1.0 / k))
    return out


def micrograph_box(im, crop, frac=0.9, white=250):
    """the micrograph's own rectangle inside the crop: the rows and columns of which at least `frac` are not
    paper-white (the figure's margins, captions and bars lie outside it)"""
    c0, r0, c1, r1 = crop
    a = im[r0:r1, c0:c1] < white
    rows = np.where(a.mean(axis=1) >= frac)[0]
    cols = np.where(a.mean(axis=0) >= frac)[0]
    if not len(rows) or not len(cols):
        return crop
    return (int(c0 + cols.min()), int(r0 + rows.min()), int(c0 + cols.max() + 1), int(r0 + rows.max() + 1))


# ------------------------------------------------------------------------------------------------ drawing
def draw_panel(stem, im, crop, segs_px, k, bar_um, caption=None, colours=None, lw=0.55, dpi=600, legend='opencif mask edges',
               clip=None, own_bar=False):
    """PNG (600 dpi) and vector PDF; the micrograph at its native pixels (nearest-neighbour), the trace as vectors;
    a caption-free `_clean.png` for composites"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    c0, r0, c1, r1 = crop
    img = im[r0:r1, c0:c1]
    H, W = img.shape
    colours = colours or [TRACE_COLOUR] * len(segs_px)
    for with_caption in (True, False):
        cap = textwrap.wrap(caption or '', 92) if with_caption else []
        fw = 3.3
        cap_h = (0.08 + len(cap) * 4.6 * 1.25 / 72.0 + 0.06) if cap else 0.0
        fh = fw * H / W + cap_h
        fig = plt.figure(figsize=(fw, fh), dpi=dpi)
        ax = fig.add_axes([0.0, cap_h / fh, 1.0, 1.0 - cap_h / fh])
        ax.imshow(img, cmap='gray', vmin=0, vmax=255, interpolation='nearest', extent=(c0 - 0.5, c1 - 0.5, r1 - 0.5, r0 - 0.5))
        lc = LineCollection([np.asarray(s) for s in segs_px], colors=colours, linewidths=lw, alpha=0.95)
        ax.add_collection(lc)
        if clip is not None:
            from matplotlib.patches import Rectangle
            rect = Rectangle((clip[0] - 0.5, clip[1] - 0.5), clip[2] - clip[0], clip[3] - clip[1], transform=ax.transData,
                             facecolor='none', edgecolor='none')
            ax.add_patch(rect)
            lc.set_clip_path(rect)
        if own_bar:
            L = bar_um * k
            ax.plot([c1 - 10 - L, c1 - 10], [r1 - 12, r1 - 12], '-', color='white', lw=2.2, solid_capstyle='butt')
            ax.text(c1 - 10 - L / 2, r1 - 16, '%g µm' % bar_um, color='white', ha='center', va='bottom', fontsize=5.5)
        ax.set_xlim(c0 - 0.5, c1 - 0.5)
        ax.set_ylim(r1 - 0.5, r0 - 0.5)
        ax.axis('off')
        if cap:
            fig.text(0.015, 0.04 / fh, '\n'.join(cap), fontsize=4.6, va='bottom', ha='left', linespacing=1.25)
            fig.savefig(stem + '.png', dpi=dpi)
            fig.savefig(stem + '.pdf', dpi=dpi)
        else:
            fig.savefig(stem + '_clean.png', dpi=dpi)
        plt.close(fig)


def draw_residual_map(stem, im, crop, segs_px, k, search_px=8.0, dpi=300, signs=None, clip=None):
    """V4: every trace segment coloured by its signed offset to the nearest outline (blue inside, red outside,
    saturating at +-6 µm); grey: no outline within the window"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.collections import LineCollection
    c0, r0, c1, r1 = crop
    img = im[r0:r1, c0:c1]
    H, W = img.shape
    segs, cols = [], []
    norm = mcolors.Normalize(-6, 6)
    cmap = plt.get_cmap('coolwarm')
    signs = [1] * len(segs_px) if signs is None else signs
    cb0, rb0, cb1, rb1 = clip if clip is not None else (-1e9, -1e9, 1e9, 1e9)
    for s, sg in zip(segs_px, signs):
        s = np.asarray(s, float)
        if len(s) < 2:
            continue
        mid = 0.5 * (s[1:] + s[:-1])
        d = np.diff(s, axis=0)
        L = np.hypot(d[:, 0], d[:, 1])
        ok = (L > 1e-9) & (mid[:, 0] >= cb0) & (mid[:, 0] <= cb1) & (mid[:, 1] >= rb0) & (mid[:, 1] <= rb1)
        nrm = np.zeros_like(d)
        nrm[ok] = sg * np.c_[-d[ok, 1], d[ok, 0]] / L[ok][:, None]
        off, _ = nearest_outline(im, mid, nrm, k, search_px)
        for i in range(len(mid)):
            if not ok[i]:
                continue
            segs.append(s[i:i + 2])
            cols.append(cmap(norm(off[i])) if np.isfinite(off[i]) else (0.55, 0.55, 0.55, 0.9))
    fh = 3.3 * H / W + 0.35
    fig = plt.figure(figsize=(3.3, fh), dpi=dpi)
    ax = fig.add_axes([0, 0.35 / fh, 1, 1 - 0.35 / fh])
    ax.imshow(img, cmap='gray', vmin=0, vmax=255, interpolation='nearest', extent=(c0 - 0.5, c1 - 0.5, r1 - 0.5, r0 - 0.5), alpha=0.55)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=0.9))
    ax.set_xlim(c0 - 0.5, c1 - 0.5)
    ax.set_ylim(r1 - 0.5, r0 - 0.5)
    ax.axis('off')
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cax = fig.add_axes([0.1, 0.2 / fh, 0.8, 0.05 / fh])
    cb = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cb.ax.tick_params(labelsize=4.5)
    cb.set_label('signed offset to the nearest outline, µm (+: outside the mask polygon; grey: none within ±%g px)' % search_px, fontsize=4.5)
    fig.savefig(stem + '_residuals.png', dpi=dpi)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ scale bars
def bar_length_aa(im, r0, r1, c0, c1, dark=True, full_tol=30.0):
    """the length (px) of a horizontal scale bar inside the window, counting its anti-aliased end pixels by their
    coverage: a full pixel counts 1, an end pixel (v - bg) / (bar - bg).  bar: the bar's own level (0 or 255);
    bg: the median of 3-8 px beyond each end.  Returned: (median length over the bar's full rows, sd, rows, the
    count of full pixels alone -- what an integer run-length reading gives)"""
    lvl = 0.0 if dark else 255.0
    out = []
    for y in range(r0, r1):
        row = im[y, c0:c1]
        full = np.abs(row - lvl) <= full_tol
        if full.sum() < 10:
            continue
        idx = np.where(full)[0]
        a, b = idx.min(), idx.max() + 1
        if (b - a) != full[a:b].sum():          # not one solid run (text, an arrow)
            continue
        if a < 9:
            continue
        bgL = np.median(row[a - 8:a - 2])
        # a bar that ends within the raster's last pixels (Gifford's inset): its far end is read against the near
        # end's background (the page's white); a bar cut by the raster's edge is not measured (its end pixel is full)
        if b >= len(row):
            continue
        bgR = np.median(row[b + 2:b + 8]) if b <= len(row) - 9 else bgL
        cl = (row[a - 1] - bgL) / (lvl - bgL)
        cr = (row[b] - bgR) / (lvl - bgR)
        out.append((y, (b - a) + float(np.clip(cl, 0, 1)) + float(np.clip(cr, 0, 1)), b - a))
    if not out:
        return None
    L = np.array([o[1] for o in out])
    # the bar's interior rows (drop the first and last, which are partly covered vertically)
    core = L[1:-1] if len(L) > 3 else L
    return float(np.median(core)), float(np.std(core)), [o[0] for o in out], int(np.median([o[2] for o in out]))


def profile_line(im, p0, dirn, half_len, step=0.5, width=2.0):
    """the image sampled along p0 + s * dirn, s in [-half_len, half_len] px (step px), averaged over +-width px
    across; returns (s, profile)"""
    dirn = np.asarray(dirn, float) / np.hypot(*dirn)
    nrm = np.array([-dirn[1], dirn[0]])
    s = np.arange(-half_len, half_len + 1e-9, step)
    acc = np.zeros(len(s))
    ws = np.arange(-width, width + 1e-9, 1.0)
    for w in ws:
        P = np.asarray(p0, float) + np.outer(s, dirn) + w * nrm
        acc += map_coordinates(im, [P[:, 1], P[:, 0]], order=1, mode='nearest')
    return s, acc / len(ws)
