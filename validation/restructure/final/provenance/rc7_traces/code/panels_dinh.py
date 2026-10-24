"""Dinh et al. 2024 Lab Chip 24, 913-923, Fig. 2B (four panels, each with its own 50 µm bar): the outlet end (4th),
the linear section (3rd) and the head-to-array transition (2nd), traced with an opencif mask.

Scale: each panel's own bar, read with its anti-aliased ends (46 full pixels + 0.33 + 0.57 = 46.91 px; the integer
run length, 46 px, is what rc5's anchors and both reviewers used).  Rotation: none (the flow runs down the page).
Translation: by features (the posts' pill centres along the flow, the post columns' axis across it) and, independently,
by chamfer matching (panel_core.run_panel).  The mask is the shipped dinh2024_v1_depicted DXF or a variant built from a
preset that differs in declared fields (ctx['label'])."""
import json
import os

import numpy as np

import trace_lib as T
from panel_core import run_panel, scale_scan   # rc7: the Dinh run uses panels_dinh_v2 (window rule there)

SRC = 'Dinh et al. 2024 Lab Chip 24, 913–923'
FIG2_SHA = 'a21e2e3ce0710e32'
POST_ALONG_UM = (21.5, 29.0)          # a pill's two outlines (dark-line centres) are this far apart along the flow;
                                      # a gap's are at most 20.5 µm apart (19.2 px), so the two never mix


def setup(root, rasters, preset=None, dxf=None, record=None, label=None, tag='shipped'):
    import sys
    for p in (root,):
        if p not in sys.path:
            sys.path.insert(0, p)
    from opencif import spec as SP
    from opencif.physics import schedule as CS
    name = 'dinh2024_v1_depicted'
    preset = preset or os.path.join(root, 'presets', name + '.yaml')
    dxf = dxf or os.path.join(root, 'runs', name, name + '.dxf')
    record = record or os.path.join(root, 'runs', name, name + '.json')
    rec = json.load(open(record))
    mask = T.Mask(dxf)
    if mask.sha256 != rec['files']['dxf_sha256']:
        raise ValueError('the DXF is not its record\'s: %s' % dxf)
    spec = SP.load(preset)
    sp, op = SP.to_engine(spec)
    sc = CS.plan(sp)
    if sc.N != rec['design']['N'] or sc.NM != rec['design']['NM']:
        raise ValueError('the preset\'s plan (N %d, NM %d) is not the record\'s (N %d, NM %d)' % (sc.N, sc.NM, rec['design']['N'], rec['design']['NM']))
    r2 = T.Raster(os.path.join(rasters, 'dinh5-000.png'), FIG2_SHA)
    allposts = mask.posts(-1e9, 1e9, -1e9, 1e9, max_len=sp.g_p + 1.0, max_wid=sp.l_g + 1.0)
    return dict(mask=mask, rec=rec, sp=sp, sc=sc, r2=r2, posts=allposts, label=label or os.path.basename(dxf)[:-4],
                preset=os.path.basename(preset), dxf_sha256=mask.sha256, tag=tag, outlet=spec.get('outlet', {}))


def pills_along(im, col, r0, r1, k, half=3, depth=10):
    """pill centres along a post column: the mid-points of outline pairs a post length apart (an annotation arrow
    between the two outlines is passed over)"""
    L = T.dark_lines(T.col_profile(im, col, r0, r1, half), depth, off=r0)
    return T.pair_centres(L, POST_ALONG_UM[0] * k, POST_ALONG_UM[1] * k, skip=1), L


def across_lines(im, rows, c0, c1, cL, cR, half_post_px, tol=5.0, depth=10):
    """the six across-flow lines (outer walls, the two posts' outer and inner ends) averaged over pill rows: the
    walls are the first and last dark lines of the row; a post's ends are the lines nearest its column centre
    +- half a post width"""
    names = ('wl', 'pLo', 'pLi', 'pRi', 'pRo', 'wr')
    acc = {n: [] for n in names}
    for r in rows:
        L = T.dark_lines(T.row_profile(im, r, c0, c1, half=1), depth, off=c0)
        if len(L) < 2:
            continue
        acc['wl'].append(L[0])
        acc['wr'].append(L[-1])
        for n, x in (('pLo', cL - half_post_px), ('pLi', cL + half_post_px), ('pRi', cR - half_post_px), ('pRo', cR + half_post_px)):
            v = T.line_near(L, x, tol)
            if v is not None:
                acc[n].append(v)
    mean = {n: (float(np.mean(v)) if v else None) for n, v in acc.items()}
    sd = {n: (float(np.std(v)) if len(v) > 1 else None) for n, v in acc.items()}
    return mean, sd, {n: len(v) for n, v in acc.items()}


def widths(lines, k):
    g = lambda a, b: (lines[b] - lines[a]) / k
    return dict(outer_walls=g('wl', 'wr'), side_L=g('wl', 'pLo'), side_R=g('pRo', 'wr'), post_L=g('pLo', 'pLi'),
                post_R=g('pRi', 'pRo'), centre=g('pLi', 'pRi'),
                post_column_separation=(0.5 * (lines['pRi'] + lines['pRo']) - 0.5 * (lines['pLi'] + lines['pLo'])) / k)


def pitch_fit(centres_lists, k):
    """the post pitch (µm) from pill centres: a straight-line fit of centre against post index in each column"""
    out = []
    for c in centres_lists:
        c = np.asarray(c, float)
        if len(c) < 3:
            continue
        d = np.diff(c)
        idx = np.r_[0, np.cumsum(np.round(d / np.median(d)))]
        out.append(float(np.polyfit(idx, c, 1)[0] / k))
    return float(np.mean(out)), out


def side_posts(ctx, side):
    P = ctx['posts']
    return np.sort(P[P[:, 1] < 0, 0] if side == 'L' else P[P[:, 1] > 0, 0])


def post_centres(im, pill_rows, guess_col, half_post_px, c0=620, c1=940, tol=6.0):
    """image post centres (col, row) at the given pill rows: the mid-point of the post's two ends on the row's
    profile, found near guess_col +- half a post width (the guess follows the column down the panel)"""
    pts = []
    g = float(guess_col)
    for r in pill_rows:
        cs = []
        for dr in (-6, 6):        # 6 px above and below the pill's mid-row: clear of the caption arrows drawn along it;
            L = T.dark_lines(T.row_profile(im, r + dr, c0, c1, half=1), 10, off=c0)   # the two ends move in equally
            a, b = T.line_near(L, g - half_post_px, tol), T.line_near(L, g + half_post_px, tol)
            if a is not None and b is not None:
                cs.append((0.5 * (a + b), b - a))
        if cs:
            c = float(np.mean([q[0] for q in cs]))
            pts.append((c, float(r), float(np.mean([q[1] for q in cs]))))
            g = c
    return np.array(pts)


def fit_rotation(cols):
    """the common slope d(col)/d(row) of the post columns (a separate intercept per column): the image's rotation"""
    A, y = [], []
    for j, P in enumerate(cols):
        for c, r, _ in P:
            row = [r] + [1.0 if q == j else 0.0 for q in range(len(cols))]
            A.append(row)
            y.append(c)
    A, y = np.array(A), np.array(y)
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ sol
    m = float(sol[0])
    return m, float(np.degrees(np.arctan(m))), float(np.std(res))


def wall_angle(im, rows, guess_col, tol=4.0, c0=620, c1=940):
    """the outer wall's slope d(col)/d(row): its dark line followed down the rows from guess_col (each row's line
    nearest the previous one, within tol px)"""
    pts, g = [], float(guess_col)
    for r in rows:
        L = T.dark_lines(T.row_profile(im, r, c0, c1, half=1), 10, off=c0)
        v = T.line_near(L, g, tol)
        if v is not None:
            pts.append((v, r))
            g = v
    P = np.array(pts)
    m = np.polyfit(P[:, 1], P[:, 0], 1)[0]
    return float(m), float(np.degrees(np.arctan(m))), len(P)


def feature_map(k, m, img_pts, mask_pts):
    """the map with the bar's scale and the measured rotation (flow along (m, 1) in image px), its origin by least
    squares over matched (image, mask) post centres; returns the map and the residuals (µm, along / across)"""
    theta = float(np.arctan2(1.0, m))
    M0 = T.rotated_map(k, theta)
    O = np.mean(np.asarray(img_pts) - M0(np.asarray(mask_pts)), axis=0)
    M = T.rotated_map(k, theta, O)
    d = np.asarray(img_pts) - M(np.asarray(mask_pts))
    u = np.array([np.cos(theta), np.sin(theta)])
    v = np.array([-u[1], u[0]])
    return M, theta, (d @ u) / k, (d @ v) / k


def _mask_widths(ctx, i):
    sp, sc = ctx['sp'], ctx['sc']
    ws = float(sc.ws[i])
    return dict(centre=float(sc.wc[i]), side=ws, post=float(sp.l_g), outer_walls=float(sc.wc[i] + 2 * sp.l_g + 2 * ws),
                post_column_separation=float(sc.wc[i] + sp.l_g))


def _pm(sp):
    return (sp.g_p + 1.0, sp.l_g + 1.0)


def panel4(ctx, out, scan=True):
    r, mask, sp, sc = ctx['r2'], ctx['mask'], ctx['sp'], ctx['sc']
    im = r.im
    bar = T.bar_length_aa(im, 1920, 1940, 870, 975)
    k = ctx.get('k_override') or bar[0] / 50.0
    crop = (610, 1486, 945, 1935)
    pitch = sp.g_p + sp.G
    pl, _ = pills_along(im, 721, 1495, 1725, k)
    pr, _ = pills_along(im, 820, 1495, 1725, k)
    mL, mR = side_posts(ctx, 'L')[-5:], side_posts(ctx, 'R')[-5:]
    tt = np.r_[pl[-5:] - mL * k, pr[-5:] - mR * k]
    t = float(tt.mean())
    rows = [int(round(v)) for v in pr[-5:-2]]
    lines, sds, ns = across_lines(im, rows, 620, 940, 720.0, 819.5, 0.5 * sp.l_g * k)
    axis = 0.5 * (lines['pLo'] + lines['pRo'])
    fmap = T.vertical_map(k, (axis, t))
    W = widths(lines, k)
    Mw = _mask_widths(ctx, sc.N - 2)
    p_um, p_cols = pitch_fit([pl, pr], k)
    meas = dict(bar_px_aa=bar[0], bar_px_full_pixels=bar[3], px_per_um=k, pitch_um=p_um, pitch_um_by_column=p_cols,
                mask_pitch_um=float(pitch), image_widths_um=W, mask_widths_um_last_row=Mw,
                printed_um=dict(centre_out=52.0, side_out=44.0, gap=19.0),
                feature_fit_residual_um=((tt - t) / k).tolist(), lines_sd_px=sds, lines_n=ns)
    x_last = float(max(mL.max(), mR.max()))
    region = (x_last - 4 * pitch - 70, x_last + 300, -200, 200)
    o = ctx['outlet']
    head = ('%s, Fig. 2B, 4th panel (the outlet end; the panel\'s 50 µm bar, %.2f px) with the edges of %s: the last five post rows '
            'and the outlet stretch (%s µm straight from the last gap\'s centre, then %s° per side).' % (
                SRC, bar[0], ctx['label'], o.get('flare_start_um'), o.get('flare_deg')))
    tail = ('On the image (dark-line centres, µm at the bar): post pitch %.2f (mask %.1f); post columns %.1f apart (mask %.1f); '
            'outer walls %.1f (mask %.1f); centre %.1f (printed 52); sides %.1f / %.1f (printed 44).' % (
                p_um, pitch, W['post_column_separation'], Mw['post_column_separation'], W['outer_walls'], Mw['outer_walls'],
                W['centre'], W['side_L'], W['side_R']))
    exclude = ((610, 1622, 660, 1660), (690, 1638, 765, 1654), (676, 1680, 752, 1697), (792, 1680, 864, 1697),
               (862, 1668, 945, 1705), (885, 1915, 945, 1935))
    rec, cmap = run_panel('dinh_fig2B_p4_outlet__' + ctx['tag'], out, r, crop, k, mask, region, fmap, (0, 1), pitch,
                          exclude=exclude, caption_head=head, caption_tail=tail, bar_um=50, extra=dict(measured=meas),
                          post_max=_pm(sp))
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, (0, 1), exclude=exclude, post_max=_pm(sp))
    return rec


def panel3(ctx, out, p4, scan=True):
    r, mask, sp, sc = ctx['r2'], ctx['mask'], ctx['sp'], ctx['sc']
    im = r.im
    bar = T.bar_length_aa(im, 1436, 1456, 870, 975)
    k = ctx.get('k_override') or bar[0] / 50.0
    crop = (610, 1000, 945, 1452)
    pitch = sp.g_p + sp.G
    pl, _ = pills_along(im, 727, 1010, 1445, k)          # col 727: clear of the vertical caption arrow at col ~717
    pr, _ = pills_along(im, 820, 1010, 1445, k)
    rows = [int(round(v)) for v in pr[[0, 1, 6, 7, 8]]]
    lines, sds, ns = across_lines(im, rows, 620, 940, 720.5, 822.0, 0.5 * sp.l_g * k, tol=7.0)
    W = widths(lines, k)
    axis = 0.5 * (lines['pLo'] + lines['pRo'])
    x = np.asarray(sc.x, float)
    wc = np.asarray(sc.wc, float)
    ws = np.asarray(sc.ws, float)
    lin = np.arange(sc.NM + sp.n_rect + 1, sc.N - 1)
    # the row on the taper: from the post columns' separation (centre to centre: free of the outline bias), i.e.
    # w_c = separation - l_g; and, as a check, from the outer-wall width corrected by the 4th panel's image-minus-mask
    i_sep = int(lin[np.argmin(np.abs(wc[lin] - (W['post_column_separation'] - sp.l_g)))])
    biasW = p4['measured']['image_widths_um']['outer_walls'] - p4['measured']['mask_widths_um_last_row']['outer_walls']
    Wm = wc + 2 * sp.l_g + 2 * ws
    i_W = int(lin[np.argmin(np.abs(Wm[lin] - (W['outer_walls'] - biasW)))])
    # along: the image's pills against consecutive mask posts centred on row i_sep (the phase, by least squares)
    mR = side_posts(ctx, 'R')
    n = len(pr)
    j0 = int(np.argmin(np.abs(mR - x[i_sep])))
    j = int(np.clip(j0 - n // 2, 0, len(mR) - n))
    used = mR[j:j + n]
    tt = pr - used * k
    t = float(tt.mean())
    res = tt - t
    i_first = int(np.argmin(np.abs(x - used[0])))
    fmap = T.vertical_map(k, (axis, t))
    Mw = _mask_widths(ctx, i_sep)
    p_um, p_cols = pitch_fit([pl, pr], k)
    meas = dict(bar_px_aa=bar[0], bar_px_full_pixels=bar[3], px_per_um=k, pitch_um=p_um, pitch_um_by_column=p_cols,
                mask_pitch_um=float(pitch), image_widths_um=W, row_from_post_columns=i_sep, row_from_outer_walls=i_W,
                outer_wall_bias_from_panel4_um=biasW, fraction_along_linear=float((i_sep - sc.NM) / (sc.N - sc.NM)),
                mask_widths_um_at_row=Mw, rows_traced=[i_first, i_first + n - 1],
                feature_fit_residual_um=(res / k).tolist(), lines_sd_px=sds, lines_n=ns)
    head = ('%s, Fig. 2B, 3rd panel (the linear section; the panel\'s bar, %.2f px per 50 µm) with %s, rows %d–%d.  The row is found '
            'from the post columns\' separation (centre to centre, %.1f µm: w_c = %.1f µm, row %d, %.0f %% along the linear section); '
            'the outer walls, corrected by the 4th panel\'s image-minus-mask offset, give row %d.' % (
                SRC, bar[0], ctx['label'], i_first, i_first + n - 1, W['post_column_separation'], W['post_column_separation'] - sp.l_g,
                i_sep, 100 * meas['fraction_along_linear'], i_W))
    tail = ('At row %d the mask has w_c %.1f, w_s %.1f, outer walls %.1f µm; the image: centre %.1f, sides %.1f / %.1f, outer walls %.1f, '
            'pitch %.2f µm (mask %.1f).' % (i_sep, Mw['centre'], Mw['side'], Mw['outer_walls'], W['centre'], W['side_L'], W['side_R'],
                                             W['outer_walls'], p_um, pitch))
    region = (used[0] - pitch, used[-1] + pitch, -200, 200)
    exclude = ((610, 1168, 662, 1202), (712, 1100, 724, 1265), (885, 1435, 945, 1452))
    rec, cmap = run_panel('dinh_fig2B_p3_linear__' + ctx['tag'], out, r, crop, k, mask, region, fmap, (0, 1), pitch,
                          exclude=exclude, caption_head=head, caption_tail=tail, bar_um=50, extra=dict(measured=meas),
                          post_max=_pm(sp))
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, (0, 1), exclude=exclude, post_max=_pm(sp))
    return rec


def panel2(ctx, out, scan=True):
    r, mask, sp, sc = ctx['r2'], ctx['mask'], ctx['sp'], ctx['sc']
    im = r.im
    bar = T.bar_length_aa(im, 952, 972, 870, 975)
    k = ctx.get('k_override') or bar[0] / 50.0
    crop = (610, 522, 945, 985)
    pitch = sp.g_p + sp.G
    pr, _ = pills_along(im, 820, 735, 980, k)
    pl, _ = pills_along(im, 721, 735, 980, k)
    x = np.asarray(sc.x, float)
    mR = side_posts(ctx, 'R')
    mR = mR[mR > x[sc.NM] - 1.0]
    n = min(len(pr), 4)
    tt = pr[:n] - mR[:n] * k
    t = float(tt.mean())
    rows = [int(round(v)) for v in pr[2:5]]
    lines, sds, ns = across_lines(im, rows, 620, 940, 721.0, 820.0, 0.5 * sp.l_g * k, tol=7.0)
    axis = 0.5 * (lines['pLo'] + lines['pRo'])
    fmap = T.vertical_map(k, (axis, t))
    W = widths(lines, k)
    i1 = int(np.argmin(np.abs(x - mR[0])))
    Mw = _mask_widths(ctx, i1 + 2)
    p_um, p_cols = pitch_fit([pl, pr], k)
    meas = dict(bar_px_aa=bar[0], bar_px_full_pixels=bar[3], px_per_um=k, pitch_um=p_um, pitch_um_by_column=p_cols,
                first_post_row=i1, image_widths_um=W, mask_widths_um_third_post_row=Mw,
                printed_um=dict(centre_in=75.0, side_in=21.0), feature_fit_residual_um=((tt - t) / k).tolist(),
                lines_sd_px=sds, lines_n=ns)
    head = ('%s, Fig. 2B, 2nd panel (the head-to-array transition; the panel\'s bar, %.2f px per 50 µm) with %s around the seam: '
            'the rectilinear rows and the first post rows (registered on the first %d posts of the right column).' % (SRC, bar[0], ctx['label'], n))
    tail = ('On the image at the 3rd–5th post rows: centre %.1f µm (printed 75), sides %.1f / %.1f (printed 21), outer walls %.1f '
            '(mask %.1f), post columns %.1f apart (mask %.1f).' % (W['centre'], W['side_L'], W['side_R'], W['outer_walls'],
                                                                    Mw['outer_walls'], W['post_column_separation'], Mw['post_column_separation']))
    region = (mR[0] - 300, mR[0] + 300, -200, 200)
    exclude = ((610, 776, 672, 796), (688, 776, 762, 796), (672, 748, 744, 768), (806, 748, 882, 768), (885, 738, 945, 776),
               (885, 965, 945, 985))
    rec, cmap = run_panel('dinh_fig2B_p2_transition__' + ctx['tag'], out, r, crop, k, mask, region, fmap, (0, 1), pitch,
                          exclude=exclude, caption_head=head, caption_tail=tail, bar_um=50, extra=dict(measured=meas),
                          post_max=_pm(sp))
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, (0, 1), exclude=exclude, post_max=_pm(sp))
    return rec
