"""full_device.py -- the whole generated device, every polygon of the shipped DXF, with the papers' micrographs placed at
the positions their traces registered, so that one figure goes from the whole device down to a micrograph's own pixels.

Each micrograph is placed by the inverse of the map its trace stored (`map_chamfer` in the panel's results JSON: image
px = origin + x ex + y ey, x and y in mask µm), nothing re-fitted.  A whole-device drawing (Xia's Fig. 3a) can be laid
under the device by the map `xia_fig3a_device.py` fitted (key "underlay").  Outputs, per device:

  <name>_full.pdf                 vector edges (fills rasterised), the micrographs embedded at their native pixels
  <name>_strip<i>.png             the device in strips at one declared scale (µm per px), micrographs resampled to it
  <name>_zoom_<letter>.png        each micrograph's neighbourhood at half its native µm/px: the micrograph, the mask
                                  edges over it (cyan) and the generated device continuing beyond its frame
and, for all devices together, full_device_composite.{png,pdf} at one scale, and full_device_index.json: every frame's
corners in mask µm, every check below, and each panel's trace validators copied from its results.

Checks (all must hold, else the script exits 2 and says which):
  V9  placement round trip: the placement is the exact inverse of the stored registration (|p - B(M p + o - o)| over the
      frame and 400 random points, < 1e-6 µm); the map is a rotation and one scale, not mirrored as displayed (det < 0
      for image rows down, mask y up), its scale the trace's declared px_per_um (1e-9 relative); the raster's sha256 and
      the DXF's are those the trace recorded (the micrograph is placed on the very mask it was traced on).
  V10 one scale: every strip, zoom and composite axis has the declared µm per px (transData, 1e-6 relative).

usage: python3 -I full_device.py CONFIG.json OUT_DIR
CONFIG: {"strip_um_per_px": 2.0, "composite_um_per_px": 5.0, "strip_max_mm": 40, "pdf_dpi": 600 (optional),
         "devices": [{"name", "dxf", "title", "panels": [{"letter", "results", "key", "raster", "title"}],
                      "underlay": {"results": FIG3A_JSON, "key": "underlay", "raster": PATH} | null, "note": "..."}]}"""
import gc
import json
import math
import os
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.transforms import Affine2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trace_lib as T

FLUID, ISLAND, EDGE, ORANGE = '#cfe3f7', '#ffffff', '#1b4f9c', '#ff6a00'
CYAN = T.TRACE_COLOUR
FAIL = []


def check(ok, what):
    if not ok:
        FAIL.append(what)
    return bool(ok)


def load_rec(p):
    R = json.load(open(p['results']))
    return R[p['key']] if p.get('key', '-') not in ('-', None) else R


class Placement:
    """the inverse of a trace's stored map: mask µm <- image px"""

    def __init__(self, rec):
        m = rec['map_chamfer']
        self.o = np.array(m['origin_px'], float)
        self.M = np.c_[np.array(m['ex_px_per_um'], float), np.array(m['ey_px_per_um'], float)]   # px = o + M (x, y)
        self.B = np.linalg.inv(self.M)
        self.k = float(rec['px_per_um'])

    def to_mask(self, px):
        return (self.B @ (np.asarray(px, float) - self.o).T).T

    def to_px(self, xy):
        return self.o + (self.M @ np.asarray(xy, float).T).T

    def affine(self):
        """Affine2D taking image (col, row) to mask (x, y)"""
        t = -self.B @ self.o
        return Affine2D(np.array([[self.B[0, 0], self.B[0, 1], t[0]], [self.B[1, 0], self.B[1, 1], t[1]], [0, 0, 1]]))

    def v9(self, frame_px, rng):
        pts = np.vstack([frame_px, rng.uniform(frame_px.min(0), frame_px.max(0), (400, 2))])
        xy = self.to_mask(pts)
        err_px = float(np.abs(self.to_px(xy) - pts).max())
        back = self.to_mask(self.to_px(xy))
        err_um = float(np.abs(back - xy).max())
        ex, ey = self.M[:, 0], self.M[:, 1]
        det = float(np.linalg.det(self.M))
        sc = math.sqrt(abs(det))
        return dict(roundtrip_err_px=err_px, roundtrip_err_um=err_um, det=det, not_mirrored=det < 0,
                    orthogonal=abs(float(np.dot(ex, ey))) / (sc * sc) < 1e-9,
                    equal_axes=abs(np.hypot(*ex) - np.hypot(*ey)) / sc < 1e-9,
                    scale_px_per_um=sc, declared_px_per_um=self.k, scale_rel_err=abs(sc - self.k) / self.k)


def mask_draw_lists(mask, x0, x1, y0, y1):
    """polygons for the fills (whole polygons whose bbox meets the box) and clipped edge polylines (vectors)"""
    ch, isl = [], []
    for (L, p, bx0, bx1, by0, by1) in mask.polys:
        if bx1 < x0 or bx0 > x1 or by1 < y0 or by0 > y1:
            continue
        (ch if L == 'CHANNELS' else isl).append(p)
    edges = [np.asarray(e) for e in mask.edges(x0, x1, y0, y1)]
    return ch, isl, edges


def draw_device(ax, mask, box, lw=0.25, raster_fill=True, alpha=1.0):
    x0, x1, y0, y1 = box
    ch, isl, edges = mask_draw_lists(mask, x0, x1, y0, y1)
    pc = PolyCollection(ch, facecolors=FLUID, edgecolors='none', zorder=0, alpha=alpha)
    pi = PolyCollection(isl, facecolors=ISLAND, edgecolors='none', zorder=1, alpha=alpha)
    pc.set_rasterized(raster_fill)
    pi.set_rasterized(raster_fill)
    ax.add_collection(pc)
    ax.add_collection(pi)
    ax.add_collection(LineCollection(edges, colors=EDGE, linewidths=lw, zorder=3))
    return edges


def view_box_px(pl, rec, box, pad_px=4):
    """the part of a record's image box that the mask-µm view box (x0, x1, y0, y1) shows (plus a few px): the underlay
    cropped to the view draws the same picture with less memory (Xia's Fig. 3a under a 0.9 µm/px zoom)"""
    c0, r0, c1, r1 = rec.get('micrograph_box') or rec['crop']
    P = pl.to_px(np.array([[box[0], box[2]], [box[1], box[2]], [box[1], box[3]], [box[0], box[3]]], float))
    a0, b0 = np.floor(P.min(axis=0)) - pad_px
    a1, b1 = np.ceil(P.max(axis=0)) + pad_px
    cb = (int(max(c0, a0)), int(max(r0, b0)), int(min(c1, a1)), int(min(r1, b1)))
    return cb if cb[2] > cb[0] and cb[3] > cb[1] else None


def place_micrograph(ax, rec, raster_im, pl, letter, mask, lw_trace=0.35, label=True, interp='none', alpha=1.0, z=4, sub=None):
    c0, r0, c1, r1 = sub or rec.get('micrograph_box') or rec['crop']
    img = raster_im[r0:r1, c0:c1]
    im = ax.imshow(img, cmap='gray', vmin=0, vmax=255, interpolation=interp, extent=(c0 - 0.5, c1 - 0.5, r1 - 0.5, r0 - 0.5),
                   zorder=z, alpha=alpha)
    im.set_transform(pl.affine() + ax.transData)
    frame_px = np.array([[c0, r0], [c1, r0], [c1, r1], [c0, r1]], float) - 0.5
    corners = pl.to_mask(frame_px)
    if letter is None:
        return corners, frame_px
    frame = MplPolygon(corners, closed=True, fill=False, edgecolor=ORANGE, linewidth=0.8, zorder=z + 2)
    ax.add_patch(frame)
    xs, ys = corners[:, 0], corners[:, 1]
    edges = [np.asarray(e) for e in mask.edges(xs.min(), xs.max(), ys.min(), ys.max())]
    lc = LineCollection(edges, colors=CYAN, linewidths=lw_trace, zorder=z + 1)
    ax.add_collection(lc)
    lc.set_clip_path(frame)
    if label:
        i = int(np.argmax(ys - 1e-3 * xs))
        ax.text(xs[i], ys[i], ' ' + letter, color=ORANGE, fontsize=7, fontweight='bold', va='bottom', ha='left', zorder=z + 3)
    return corners, frame_px


def axis_scale_um_per_px(fig, ax):
    p = ax.transData.transform([[0.0, 0.0], [1000.0, 0.0], [0.0, 1000.0]])
    return 1000.0 / abs(p[1, 0] - p[0, 0]), 1000.0 / abs(p[2, 1] - p[0, 1])


def exact_axes(fig_w_px, fig_h_px, ax_box_px, dpi=100):
    """a figure of exactly fig_w_px x fig_h_px and an axes at the pixel box (left, bottom, width, height)"""
    fig = plt.figure(figsize=(fig_w_px / dpi, fig_h_px / dpi), dpi=dpi)
    l, b, w, h = ax_box_px
    ax = fig.add_axes([l / fig_w_px, b / fig_h_px, w / fig_w_px, h / fig_h_px])
    ax.set_axis_off()
    return fig, ax


def scale_bar(ax, x, y, um, text_size=6, colour='k', lw=2.0):
    ax.plot([x, x + um], [y, y], '-', color=colour, lw=lw, solid_capstyle='butt', zorder=20)
    ax.text(x + um / 2, y, ('%g mm' % (um / 1000)) if um >= 1000 else ('%g µm' % um), ha='center', va='bottom',
            fontsize=text_size, color=colour, zorder=20)


def device_bbox(mask):
    a = np.array([[p[2], p[3], p[4], p[5]] for p in mask.polys])
    return float(a[:, 0].min()), float(a[:, 1].max()), float(a[:, 2].min()), float(a[:, 3].max())


def nice_bar(length_um):
    for v in (100, 200, 250, 500, 1000, 2000, 5000, 10000):
        if v >= length_um / 12:
            return v
    return 10000


def render_device(dev, cfg, out, index, rng):
    name = dev['name']
    mask = T.Mask(dev['dxf'])
    X0, X1, Y0, Y1 = device_bbox(mask)
    pad = 0.04 * (Y1 - Y0) + 150.0
    Lx, Ly = X1 - X0, Y1 - Y0
    n = max(1, int(math.ceil(Lx / (cfg.get('strip_max_mm', 40.0) * 1000.0))))
    seg = Lx / n
    s = float(cfg.get('strip_um_per_px', 2.0))
    rasters = {}
    panels = []
    for p in dev.get('panels', []):
        rec = load_rec(p)
        if p['raster'] not in rasters:
            rasters[p['raster']] = T.Raster(p['raster'])
        R = rasters[p['raster']]
        pl = Placement(rec)
        c0, r0, c1, r1 = rec.get('micrograph_box') or rec['crop']
        frame_px = np.array([[c0, r0], [c1, r0], [c1, r1], [c0, r1]], float) - 0.5
        v9 = pl.v9(frame_px, rng)
        v9['raster_sha256_matches'] = (R.sha256 == rec['raster_sha256'])
        v9['dxf_sha256_matches'] = (mask.sha256 == rec['mask_sha256'])
        ok = (v9['roundtrip_err_um'] < 1e-6 and v9['roundtrip_err_px'] < 1e-6 and v9['not_mirrored'] and v9['orthogonal']
              and v9['equal_axes'] and v9['scale_rel_err'] < 1e-9 and v9['raster_sha256_matches'] and v9['dxf_sha256_matches'])
        v9['ok'] = check(ok, '%s %s: V9 placement %s' % (name, p['letter'], {k: v for k, v in v9.items() if k != 'ok'}))
        corners = pl.to_mask(frame_px)
        tv = dict(V1_px=rec['V1_feature_vs_chamfer_px'], median_um=rec['metric'].get('median_abs_um'),
                  p90_um=rec['metric'].get('p90_abs_um'), null_half_pitch_um=rec['metric_null']['half_pitch_along'].get('median_abs_um'),
                  null_5um_across_um=rec['metric_null']['5um_across'].get('median_abs_um'),
                  V5=dict(along=(rec.get('V5_scale_scan') or {}).get('best_along'), across=(rec.get('V5_scale_scan') or {}).get('best_across')),
                  um_per_px=1.0 / rec['px_per_um'], trace=rec['name'])
        panels.append(dict(p=p, rec=rec, pl=pl, R=R, corners=corners, v9=v9, trace=tv))
    under = None
    if dev.get('underlay'):
        u = dev['underlay']
        urec = json.load(open(u['results']))
        if u.get('key'):                                  # xia_fig3a_device.py nests it under "underlay"
            urec = urec[u['key']]
        if u['raster'] not in rasters:
            rasters[u['raster']] = T.Raster(u['raster'])
        under = (urec, Placement(urec), rasters[u['raster']])
    D = dict(name=name, dxf=os.path.basename(dev['dxf']), dxf_sha256=mask.sha256, bbox_um=[X0, X1, Y0, Y1],
             footprint_mm=[Lx / 1000.0, Ly / 1000.0], strips=n, strip_um_per_px=s, panels={}, files=[])

    def draw_all(ax, box, lw=0.25):
        if under is not None:
            urec, upl, uR = under
            sb = view_box_px(upl, urec, box)
            if sb is not None:
                place_micrograph(ax, urec, uR.im, upl, None, mask, interp='bilinear', alpha=0.55, z=-1, sub=sb)
        draw_device(ax, mask, box, lw=lw)
        for q in panels:
            place_micrograph(ax, q['rec'], q['R'].im, q['pl'], q['p']['letter'], mask)

    # --- strips (PNG at s µm/px) and the vector PDF (one page, the strips stacked)
    head_px = 70
    strip_boxes = [(X0 + i * seg - 50.0, X0 + (i + 1) * seg + 50.0, Y0 - pad, Y1 + pad) for i in range(n)]
    strip_boxes = [(x0, x0 + math.ceil((x1 - x0) / s) * s, y0, y0 + math.ceil((y1 - y0) / s) * s) for x0, x1, y0, y1 in strip_boxes]
    for i, (x0, x1, y0, y1) in enumerate(strip_boxes):
        W, H = int(round((x1 - x0) / s)), int(round((y1 - y0) / s))      # exact: the box is a whole number of px
        fig, ax = exact_axes(W, H + head_px, (0, 0, W, H))
        draw_all(ax, (x0, x1, y0, y1))
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        sx, sy = axis_scale_um_per_px(fig, ax)
        check(abs(sx - s) / s < 1e-6 and abs(sy - s) / s < 1e-6, '%s strip %d: V10 scale %.6f/%.6f != %.3f' % (name, i + 1, sx, sy, s))
        b = nice_bar(x1 - x0)
        scale_bar(ax, x1 - b - 300.0, y0 + 0.25 * pad, b, text_size=8)
        fig.text(0.005, 1.0 - 8.0 / (H + head_px), '%s — strip %d of %d, x %.2f to %.2f mm, %.1f µm per px (every polygon of %s, sha256 %s…); '
                 'orange frames: micrographs at their registered positions; cyan: the mask edges over them'
                 % (dev.get('title', name), i + 1, n, x0 / 1000, x1 / 1000, s, D['dxf'], mask.sha256[:12]),
                 fontsize=9, va='top')
        f = os.path.join(out, '%s_strip%d.png' % (name, i + 1))
        fig.savefig(f, dpi=100)
        plt.close(fig)
        gc.collect()                                   # free the canvas before the next figure (Xia: GBs)
        D['files'].append(os.path.basename(f))
        D.setdefault('V10', []).append(dict(file=os.path.basename(f), um_per_px=[sx, sy]))
    # PDF: one page, strips stacked, physical width 16 in
    fw = 16.0
    hs = [(b[3] - b[2]) / (b[1] - b[0]) * fw for b in strip_boxes]
    cap = textwrap.wrap(dev.get('note', ''), 220)
    head_in = 0.5 + 0.13 * len(cap)
    fh = sum(hs) + 0.25 * n + head_in
    fig = plt.figure(figsize=(fw, fh))
    ytop = fh - head_in
    for i, b in enumerate(strip_boxes):
        ytop -= hs[i]
        ax = fig.add_axes([0.0, ytop / fh, 1.0, hs[i] / fh])
        ax.set_axis_off()
        draw_all(ax, b, lw=0.12)
        ax.set_xlim(b[0], b[1])
        ax.set_ylim(b[2], b[3])
        bb = nice_bar(b[1] - b[0])
        scale_bar(ax, b[1] - bb - 300.0, b[2] + 0.25 * pad, bb, text_size=6)
        ytop -= 0.25
    fig.text(0.005, 1.0 - 0.1 / fh, dev.get('title', name), fontsize=9, fontweight='bold', va='top')
    fig.text(0.005, 1.0 - 0.32 / fh, '\n'.join(cap), fontsize=5.5, va='top')
    f = os.path.join(out, '%s_full.pdf' % name)
    fig.savefig(f, dpi=float(cfg.get('pdf_dpi', 600)))      # the rasterised fills and images
    plt.close(fig)
    gc.collect()
    D['files'].append(os.path.basename(f))
    # --- zooms: each micrograph's neighbourhood at half its native µm/px
    for q in panels:
        letter = q['p']['letter']
        cx, cy = q['corners'].mean(axis=0)
        hw = 1.0 * np.abs(q['corners'][:, 0] - cx).max()
        hh = 1.0 * np.abs(q['corners'][:, 1] - cy).max()
        sz = 0.5 / q['rec']['px_per_um']
        W = H = int(math.ceil(4.0 * max(hw, hh) / sz))
        half = 0.5 * W * sz
        box = (cx - half, cx + half, cy - half, cy + half)
        fig, ax = exact_axes(W, H + 90, (0, 0, W, H))
        if under is not None:
            urec, upl, uR = under
            sb = view_box_px(upl, urec, box)
            if sb is not None:
                place_micrograph(ax, urec, uR.im, upl, None, mask, interp='bilinear', alpha=0.55, z=-1, sub=sb)
        draw_device(ax, mask, box, lw=0.6)
        place_micrograph(ax, q['rec'], q['R'].im, q['pl'], letter, mask, lw_trace=0.8, interp='nearest')
        ax.set_xlim(box[0], box[1])
        ax.set_ylim(box[2], box[3])
        sx, sy = axis_scale_um_per_px(fig, ax)
        check(abs(sx - sz) / sz < 1e-6 and abs(sy - sz) / sz < 1e-6, '%s zoom %s: V10 scale' % (name, letter))
        b = nice_bar(2 * half)
        scale_bar(ax, box[1] - b - 0.05 * half, box[2] + 0.05 * half, b, text_size=10, lw=3)
        tv = q['trace']
        fig.text(0.005, 1.0 - 6.0 / (H + 90), '%s, panel %s (%s) on %s at %.3f µm per px (half the micrograph\'s %.3f): the micrograph at its '
                 'registered position (orange frame), the mask edges over it (cyan), the generated device beyond it.  Trace: '
                 'registrations agree %.2f / %.2f px; median offset %.2f µm (null %.2f half a pitch along, %.2f 5 µm across); V5 %s / %s'
                 % (dev.get('title', name), letter, q['p'].get('title', ''), D['dxf'], sz, 1.0 / q['rec']['px_per_um'],
                    abs(tv['V1_px']['along']), abs(tv['V1_px']['across']), tv['median_um'], tv['null_half_pitch_um'],
                    tv['null_5um_across_um'], ('%.3f' % tv['V5']['along']) if tv['V5']['along'] else '—',
                    ('%.3f' % tv['V5']['across']) if tv['V5']['across'] else '—'), fontsize=10, va='top', wrap=True)
        f = os.path.join(out, '%s_zoom_%s.png' % (name, letter))
        fig.savefig(f, dpi=100)
        plt.close(fig)
        gc.collect()                                   # free the canvas before the next figure (Xia: GBs)
        D['files'].append(os.path.basename(f))
        D['panels'][letter] = dict(title=q['p'].get('title', ''), results=os.path.basename(q['p']['results']), key=q['p'].get('key', '-'),
                                   corners_um=q['corners'].tolist(), centre_um=[float(cx), float(cy)], V9=q['v9'], trace=tv,
                                   zoom_um_per_px=sz)
    index['devices'][name] = D
    return mask, strip_boxes, panels, under


def composite(cfg, out, rendered, index):
    s = float(cfg.get('composite_um_per_px', 5.0))
    rows = []
    for name, (mask, boxes, panels, under) in rendered.items():
        for b in boxes:
            b = (b[0], b[0] + math.ceil((b[1] - b[0]) / s) * s, b[2], b[2] + math.ceil((b[3] - b[2]) / s) * s)
            rows.append((name, mask, b, panels, under))
    W = int(round(max((b[1] - b[0]) for _, _, b, _, _ in rows) / s))
    gap = 60
    hs = [int(round((b[3] - b[2]) / s)) for _, _, b, _, _ in rows]
    H = sum(hs) + gap * (len(rows) + 1)
    fig = plt.figure(figsize=(W / 100.0, H / 100.0), dpi=100)
    y = H - gap
    out_scales = []
    for (name, mask, b, panels, under), h in zip(rows, hs):
        y -= h
        w = int(round((b[1] - b[0]) / s))
        ax = fig.add_axes([0.0, y / H, w / W, h / H])
        ax.set_axis_off()
        if under is not None:
            urec, upl, uR = under
            sb = view_box_px(upl, urec, b)
            if sb is not None:
                place_micrograph(ax, urec, uR.im, upl, None, mask, interp='bilinear', alpha=0.55, z=-1, sub=sb)
        draw_device(ax, mask, b, lw=0.2)
        for q in panels:
            place_micrograph(ax, q['rec'], q['R'].im, q['pl'], q['p']['letter'], mask, label=True)
        ax.set_xlim(b[0], b[1])
        ax.set_ylim(b[2], b[3])
        sx, sy = axis_scale_um_per_px(fig, ax)
        check(abs(sx - s) / s < 1e-6 and abs(sy - s) / s < 1e-6, 'composite %s: V10 scale %.6f' % (name, sx))
        out_scales.append([name, sx, sy])
        fig.text(2.0 / W, (y + h + 8) / H, '%s  x %.1f-%.1f mm' % (name, b[0] / 1000, b[1] / 1000), fontsize=10, va='bottom')
        y -= gap
    for ext in ('png', 'pdf'):
        f = os.path.join(out, 'full_device_composite.' + ext)
        fig.savefig(f, dpi=100)
    plt.close(fig)
    gc.collect()
    index['composite'] = dict(um_per_px=s, rows=[[r[0], r[2]] for r in rows], V10=out_scales,
                              files=['full_device_composite.png', 'full_device_composite.pdf'])


def main(cfg_path, out):
    cfg = json.load(open(cfg_path))
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(1)
    index = dict(config=os.path.basename(cfg_path), devices={})
    rendered = {}
    for dev in cfg['devices']:
        rendered[dev['name']] = render_device(dev, cfg, out, index, rng)
        print('rendered', dev['name'], flush=True)
    if cfg.get('composite', True):
        composite(cfg, out, rendered, index)
    index['failures'] = FAIL
    index['ok'] = not FAIL
    json.dump(index, open(os.path.join(out, 'full_device_index.json'), 'w'), indent=1, default=float)
    for f in FAIL:
        print('FAIL', f)
    print('full_device: %s (%d failures)' % ('ok' if not FAIL else 'FAILED', len(FAIL)))
    sys.exit(0 if not FAIL else 2)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
