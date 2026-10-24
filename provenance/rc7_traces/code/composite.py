"""Dissertation figures from the trace records: each panel is redrawn from its raster, its DXF and its registered map
(the record's map_chamfer), so the figure is regenerated from data, not pasted from earlier renders.

usage: python3 composite.py SPEC.json OUT_STEM
SPEC: {"title": ..., "columns": n, "width_in": 7.0, "panels": [{"results": path, "key": "p4" | "-", "raster": path,
       "dxf": path, "letter": "a", "title": "...", "claim_y_above": null | float, "post_max": [27, 55],
       "lw": 0.45, "extra_lines": [...], "claim_tags_any": ["channel"]}, ...],
       "footer": "..."}"""
import json
import os
import sys
import textwrap

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trace_lib as T

CYAN, AMBER = T.TRACE_COLOUR, T.NOCLAIM_COLOUR


def load(p):
    R = json.load(open(p['results']))
    return R[p['key']] if p.get('key', '-') != '-' else R


def panel_segments(rec, p, mask):
    mp = rec['map_chamfer']
    M = T.Map(mp['origin_px'], mp['ex_px_per_um'], mp['ey_px_per_um'])
    pm = tuple(p.get('post_max', (70.0, 130.0)))
    E = mask.edges(*rec['region_um'], post_max=pm)
    yc = p.get('claim_y_above')
    segs, cols = [], []
    for e in E:
        claim = True if yc is None else (float(np.mean(e[:, 1])) > yc or e.tag in p.get('claim_tags_any', ()))
        segs.append(M(e))
        cols.append(CYAN if claim else AMBER)
    return segs, cols


def main(spec_path, stem):
    S = json.load(open(spec_path))
    P = S['panels']
    ncol = S.get('columns', len(P))
    nrow = int(np.ceil(len(P) / ncol))
    recs = [load(p) for p in P]
    rasters = {}
    masks = {}
    boxes = []
    for p, r in zip(P, recs):
        if p['raster'] not in rasters:
            rasters[p['raster']] = T.Raster(p['raster']).im
        if p['dxf'] not in masks:
            masks[p['dxf']] = T.Mask(p['dxf'])
        boxes.append(tuple(r.get('micrograph_box', r['crop'])))
    W = S.get('width_in', 7.0)
    gap = 0.06
    colw = (W - gap * (ncol + 1)) / ncol
    heights = []
    for i in range(nrow):
        hs = [colw * (b[3] - b[1]) / (b[2] - b[0]) for b in boxes[i * ncol:(i + 1) * ncol]]
        heights.append(max(hs))
    title_h, metric_h = 0.24, 0.30 + 0.075 * max(len(p.get('extra_lines', [])) for p in P)
    foot = textwrap.wrap(S.get('footer', ''), int(150 * W / 7.0))      # 150 characters at 4.6 pt fill a 7 in figure
    foot_h = 0.08 + 0.083 * len(foot)
    H = sum(h + title_h + metric_h for h in heights) + foot_h + 0.05
    fig = plt.figure(figsize=(W, H), dpi=600)
    y = H
    for i, (p, r, b) in enumerate(zip(P, recs, boxes)):
        row, col = divmod(i, ncol)
        if col == 0:
            y -= title_h + heights[row]
            ytop = y
        x0 = gap + col * (colw + gap)
        h = colw * (b[3] - b[1]) / (b[2] - b[0])
        ax = fig.add_axes([x0 / W, (ytop + heights[row] - h) / H, colw / W, h / H])
        im = rasters[p['raster']]
        c0, r0, c1, r1 = b
        ax.imshow(im[r0:r1, c0:c1], cmap='gray', vmin=0, vmax=255, interpolation='nearest', extent=(c0 - .5, c1 - .5, r1 - .5, r0 - .5))
        segs, cols = panel_segments(r, p, masks[p['dxf']])
        lc = LineCollection(segs, colors=cols, linewidths=p.get('lw', 0.45), alpha=0.95)
        ax.add_collection(lc)
        rect = Rectangle((c0 - .5, r0 - .5), c1 - c0, r1 - r0, transform=ax.transData, facecolor='none', edgecolor='none')
        ax.add_patch(rect)
        lc.set_clip_path(rect)
        ax.set_xlim(c0 - .5, c1 - .5)
        ax.set_ylim(r1 - .5, r0 - .5)
        ax.axis('off')
        fig.text(x0 / W, (ytop + heights[row] + 0.05) / H, '%s  %s' % (p['letter'], p['title']), fontsize=6.0, ha='left', va='bottom',
                 fontweight='bold')
        m = r['metric']
        v5 = r.get('V5_scale_scan', {})
        lines = ['median |d| %.2f µm · p90 %.2f µm' % (m['median_abs_um'], m['p90_abs_um']),
                 'null: ½ pitch %.2f · 5 µm across %.2f µm' % (r['metric_null']['half_pitch_along']['median_abs_um'],
                                                              r['metric_null']['5um_across']['median_abs_um'])]
        if v5:
            lines.append('scale fit along %.3f · across %.3f' % (v5['best_along'], v5['best_across']))
        lines += p.get('extra_lines', [])
        fig.text(x0 / W, (ytop + heights[row] - h - 0.03) / H, '\n'.join(lines), fontsize=4.6, ha='left', va='top', linespacing=1.2)
        if col == ncol - 1 or i == len(P) - 1:
            y -= metric_h
    if foot:
        fig.text(gap / W, 0.05 / H, '\n'.join(foot), fontsize=4.6, ha='left', va='bottom', linespacing=1.25)
    fig.savefig(stem + '.pdf')
    fig.savefig(stem + '.png', dpi=600)
    plt.close(fig)
    print(stem)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
