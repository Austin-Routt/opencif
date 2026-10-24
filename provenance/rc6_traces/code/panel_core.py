"""panel_core.py -- the common steps of one traced panel: chamfer registration beside the feature registration,
the validators V1-V4, the nearest-outline metric with its null values, the figure files and the result record."""
import json
import math
import os

import numpy as np

import trace_lib as T

DRAW = True          # tests set this False: the validators without the figure files


def run_panel(name, out_dir, raster, crop, k, mask, region, fmap, along, pitch_um, exclude=(), search_across_px=8.0,
              search_along=None, claim_filter=None, caption_head='', caption_tail='', bar_um=50, extra=None,
              post_max=(70.0, 130.0), noclaim_filter=None, draw=True, min_component=0, scale_src='scale from the bar',
              metric_search_px=8.0):
    """name: the panel's file stem; raster: trace_lib.Raster; crop: (c0, r0, c1, r1) px; k: px per µm (the bar);
    mask: trace_lib.Mask; region: (x0, x1, y0, y1) mask µm whose edges are traced; fmap: trace_lib.Map from the feature
    registration; along: unit vector of the flow in image px (for the class split and the periodic null);
    pitch_um: the post pitch (the periodic ambiguity along the flow); scale_src: the caption's words for where the
    scale k came from (rc6: Xia's panels are at the scale of their printed widths, not panel i's bar); search_along: (lo, hi) px of the chamfer search
    along the flow around fmap (default +-0.6 pitch); claim_filter(edge) -> bool: the edges that are a reproduction
    claim (default all); the others are drawn in the no-claim colour and left out of the metric; metric_search_px: the
    nearest-outline metric's window along each normal (+-8 px, 8.5-22 µm on the panels at 1.07-2.79 µm/px; Gifford's
    Fig. 6 at 4.52 µm/px uses +-2.5 px, 11.3 µm, so that the window does not reach the next post across a 17 µm gap)."""
    c0, r0, c1, r1 = crop
    edges = mask.edges(*region, post_max=post_max)
    claim = [e for e in edges if (claim_filter(e) if claim_filter else True)]
    other = [e for e in edges if not (claim_filter(e) if claim_filter else True)]
    dt, ridge, thr = raster.outline_distance(crop, exclude=exclude, min_component=min_component)
    along = np.asarray(along, float) / np.hypot(*along)
    across = np.array([-along[1], along[0]])
    pitch_px = pitch_um * k
    sa = search_along or (-0.6 * pitch_px, 0.6 * pitch_px)
    det = float(fmap.ex[0] * fmap.ey[1] - fmap.ex[1] * fmap.ey[0])
    signs = [int(e.orient * np.sign(det)) for e in claim]          # every normal points out of its polygon on the image
    P0, N0 = T.sample_edges([fmap(e) for e in claim], 0.5, signs)
    # the chamfer search is done in (along, across) coordinates, then mapped to image px
    best, best_s = None, None
    grid = []
    for u in np.arange(sa[0], sa[1] + 1e-9, 0.5):
        row = []
        for v in np.arange(-search_across_px, search_across_px + 1e-9, 0.5):
            s = T.chamfer(dt, crop, P0 + u * along + v * across)[0]
            row.append(s)
            if np.isfinite(s) and (best_s is None or s < best_s):
                best_s, best = s, (u, v)
        grid.append(row)
    grid = np.array(grid)
    us = np.arange(sa[0], sa[1] + 1e-9, 0.5)
    vs = np.arange(-search_across_px, search_across_px + 1e-9, 0.5)
    iu, iv = int(np.where(np.isclose(us, best[0]))[0][0]), int(np.where(np.isclose(vs, best[1]))[0][0])

    def para(arr, i, step=0.5):
        if 0 < i < len(arr) - 1:
            a, b, c = arr[i - 1], arr[i], arr[i + 1]
            den = a - 2 * b + c
            if den > 0:
                return 0.5 * (a - c) / den * step
        return 0.0
    u_best = best[0] + para(grid[:, iv], iu)
    v_best = best[1] + para(grid[iu, :], iv)
    cmap = fmap.shifted(u_best * along + v_best * across)
    P1, N1 = T.sample_edges([cmap(e) for e in claim], 0.5, signs)
    # V2: the null -- the whole landscape's median, and one pitch, half a pitch and 5 µm away.  Every score of the
    # comparison is taken on the same trace points: those that stay inside the crop (8 px in) under every shift.
    shifts = [s * pitch_px * along for s in (-1, 1)] + [s * 0.5 * pitch_px * along for s in (-1, 1)] + \
             [s * 5.0 * k * across for s in (-1, 1)]
    common = T.in_rects(P1, [crop], pad=-8.0)
    for dv in shifts:
        common &= T.in_rects(P1 + dv, [crop], pad=-8.0)
    Pc = P1[common]
    s_best = T.chamfer(dt, crop, Pc)[0]
    s_pm = [T.chamfer(dt, crop, Pc + dv)[0] for dv in shifts[0:2]]
    s_half = [T.chamfer(dt, crop, Pc + dv)[0] for dv in shifts[2:4]]
    s_acc = [T.chamfer(dt, crop, Pc + dv)[0] for dv in shifts[4:6]]
    # the nearest-outline metric at the registration and at the null shifts
    # (points whose +-8 px window touches an excluded annotation rectangle, or that lie outside the crop, are left out)
    def metric_at(P, Nn, keep):
        d, _ = T.nearest_outline(raster.im, P, Nn, k, metric_search_px)
        d[T.in_rects(P, exclude, pad=8.0)] = np.nan
        d[~keep] = np.nan
        return d
    d_best = metric_at(P1, N1, T.in_rects(P1, [crop], pad=-8.0))
    m_best = T.summarise(d_best)
    nulls = {}
    Nc = N1[common]
    m_common = T.summarise(metric_at(Pc, Nc, np.ones(len(Pc), bool)))
    for lab, dv in (('half_pitch_along', 0.5 * pitch_px * along), ('5um_across', 5.0 * k * across),
                    ('one_pitch_along', pitch_px * along)):
        nulls[lab] = T.summarise(metric_at(Pc + dv, Nc, np.ones(len(Pc), bool)))
    nulls['registered_same_points'] = m_common
    # V1: the two registrations agree
    dshift = cmap.origin - fmap.origin
    d_along, d_across = float(np.dot(dshift, along)), float(np.dot(dshift, across))
    v1 = abs(d_along) <= 1.5 and abs(d_across) <= 1.5
    cls = T.by_class([cmap(e) for e in claim], [e.tag for e in claim], raster.im, k, crop, along, search_px=metric_search_px,
                     exclude=exclude, signs=signs)
    rec = dict(
        name=name, raster=os.path.basename(raster.path), raster_sha256=raster.sha256, mask=os.path.basename(mask.path),
        mask_sha256=mask.sha256, crop=list(crop), px_per_um=k, pitch_um=pitch_um,
        map_feature=fmap.to_dict(), map_chamfer=cmap.to_dict(),
        V1_feature_vs_chamfer_px=dict(along=d_along, across=d_across, agree=bool(v1)),
        V2_chamfer=dict(points=int(len(Pc)), best_px=s_best, one_pitch_along_px=s_pm, half_pitch_along_px=s_half, five_um_across_px=s_acc,
                        landscape_median_px=float(np.nanmedian(grid)),
                        ratio_landscape_median_to_best=float(np.nanmedian(grid) / s_best) if s_best else None,
                        ratio_one_pitch_to_best=float(min(s_pm) / s_best) if s_best else None),
        metric=m_best, metric_null=nulls, V3_by_class=cls, outline_threshold=thr,
        claim_edges=len(claim), other_edges=len(other), region_um=[float(v) for v in region], along_px=[float(v) for v in along],
        exclude_px=[list(map(int, e)) for e in exclude], metric_search_px=float(metric_search_px))
    if extra:
        rec.update(extra)
    if draw and DRAW:
        segs = [cmap(e) for e in claim] + [cmap(e) for e in other]
        cols = [T.TRACE_COLOUR] * len(claim) + [T.NOCLAIM_COLOUR] * len(other)
        cap = (caption_head + '  Registration: ' + scale_src + ', %s; translation fitted on features and, independently, by chamfer '
               'matching (they agree within %.1f / %.1f px along / across).  Nearest-outline offset: median %.2f µm, 90th percentile %.2f µm '
               '(n %d); for comparison %.2f µm with the mask shifted half a pitch along and %.2f µm with it shifted 5 µm across.  ' % (
                   'no rotation' if abs(fmap.ex[0]) < 1e-12 else 'rotation measured on the image',
                   abs(d_along), abs(d_across), m_best.get('median_abs_um', float('nan')), m_best.get('p90_abs_um', float('nan')),
                   m_best.get('n', 0), nulls['half_pitch_along'].get('median_abs_um', float('nan')),
                   nulls['5um_across'].get('median_abs_um', float('nan'))) + caption_tail)
        rec['caption'] = cap
        stem = os.path.join(out_dir, name)
        box = T.micrograph_box(raster.im, crop)
        rec['micrograph_box'] = list(box)
        T.draw_panel(stem, raster.im, crop, segs, k, bar_um, cap, cols, clip=box)
        T.draw_residual_map(stem, raster.im, crop, [cmap(e) for e in claim], k, search_px=metric_search_px, signs=signs, clip=box)
    return rec, cmap


def scale_scan(raster, crop, mask, region, cmap, along, exclude=(), claim_filter=None, post_max=(70.0, 130.0),
               factors=None, shift_px=3.0, step=0.5, min_component=0):
    """V5: is the bar's scale the one the image prefers, separately along and across the flow?  The trace is
    rescaled about the region's centre by (s_along, s_across) -- mask x is the flow direction in every traced
    panel -- and, at each pair, re-registered by chamfer over +-shift_px; reported: the pair with the lowest score
    (parabolic refinement on each axis), the score at (1, 1), and the scale-only profiles through (1, 1).  The
    implied pitch is the mask pitch times s_along: a check of the mask's along-flow geometry against the image at the
    bar's scale."""
    factors = np.round(np.arange(0.95, 1.0501, 0.005), 4) if factors is None else np.asarray(factors)
    edges = mask.edges(*region, post_max=post_max)
    claim = [e for e in edges if (claim_filter(e) if claim_filter else True)]
    dt, _, _ = raster.outline_distance(crop, exclude=exclude, min_component=min_component)
    xc, yc = 0.5 * (region[0] + region[1]), 0.5 * (region[2] + region[3])
    c_img = cmap(np.array([[xc, yc]]))[0]
    along = np.asarray(along, float) / np.hypot(*along)
    across = np.array([-along[1], along[0]])
    pts_um = np.vstack([np.asarray(e, float) for e in claim])
    # sample the mask edges in µm finely enough for any scale in the range (0.5 px at the largest scale)
    k = float(np.hypot(*cmap.ex))
    segs_um = [np.asarray(e, float) for e in claim]
    P_um, _ = T.sample_edges(segs_um, 0.5 / k)
    rel = P_um - (xc, yc)
    shifts = np.arange(-shift_px, shift_px + 1e-9, step)
    S = np.full((len(factors), len(factors)), np.nan)
    for i, sa in enumerate(factors):
        for j, sc in enumerate(factors):
            base = c_img + np.outer(rel[:, 0] * sa, cmap.ex) + np.outer(rel[:, 1] * sc, cmap.ey)
            best = np.inf
            for u in shifts:
                for v in shifts:
                    s = T.chamfer(dt, crop, base + u * along + v * across)[0]
                    if np.isfinite(s) and s < best:
                        best = s
            S[i, j] = best
    i, j = np.unravel_index(np.nanargmin(S), S.shape)

    def para(arr, m):
        if 0 < m < len(arr) - 1:
            a, b, c = arr[m - 1], arr[m], arr[m + 1]
            den = a - 2 * b + c
            if den > 0:
                return float(factors[m] + 0.5 * (a - c) / den * (factors[1] - factors[0]))
        return float(factors[m])
    i1 = int(np.argmin(np.abs(factors - 1.0)))
    return dict(best_along=para(S[:, j], i), best_across=para(S[i, :], j), score_best=float(S[i, j]),
                score_at_bar=float(S[i1, i1]), factors=factors.tolist(),
                profile_along_at_across_1=S[:, i1].tolist(), profile_across_at_along_1=S[i1, :].tolist(),
                on_grid_edge=bool(i in (0, len(factors) - 1) or j in (0, len(factors) - 1)))
