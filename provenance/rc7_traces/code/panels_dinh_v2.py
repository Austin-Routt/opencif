"""Dinh 2024 Fig. 2B panels 4, 3, 2 with the rotation measured on the image (the post columns' common slope) and the
origin fitted on matched post centres (both columns).  The helpers are panels_dinh's."""
import numpy as np

import trace_lib as T
from panel_core import run_panel, scale_scan, window_px
from panels_dinh import (SRC, pills_along, across_lines, widths, pitch_fit, side_posts, _mask_widths, _pm,
                         post_centres, fit_rotation, wall_angle, feature_map)


def _common(ctx, bar_rows, crop, colL, colR, prow0, prow1, guessL, guessR):
    r, sp = ctx['r2'], ctx['sp']
    im = r.im
    bar = T.bar_length_aa(im, bar_rows[0], bar_rows[1], 870, 975)
    k = ctx.get('k_override') or bar[0] / 50.0
    pl, _ = pills_along(im, colL, prow0, prow1, k)
    pr, _ = pills_along(im, colR, prow0, prow1, k)
    hp = 0.5 * sp.l_g * k
    PL = post_centres(im, [int(round(v)) for v in pl], guessL, hp)
    PR = post_centres(im, [int(round(v)) for v in pr], guessR, hp)
    m, deg, sd = fit_rotation([PL, PR])
    return im, bar, k, pl, pr, PL, PR, m, deg, sd


def _rows_lines(im, PL, PR, idx, k, sp):
    rows = [int(round(PR[i, 1])) for i in idx if i < len(PR)]
    cL = float(np.mean(PL[:, 0])) if len(PL) else 720.0
    cR = float(np.mean(PR[:, 0]))
    return across_lines(im, rows, 620, 940, cL, cR, 0.5 * sp.l_g * k, tol=6.0)


def _register(k, m, PL, PR, mL, mR, used_rows_L, used_rows_R):
    img = np.r_[PL[used_rows_L, :2], PR[used_rows_R, :2]]
    msk = np.r_[np.c_[mL, np.full(len(mL), -1.0)], np.c_[mR, np.full(len(mR), 1.0)]]
    return img, msk


def _yfix(ctx, msk):
    """the mask post centres' y (the post columns at +-(w_c + l_g)/2 of their row)"""
    P = ctx['posts']
    out = msk.copy()
    for i, (x, s) in enumerate(msk):
        sel = P[(np.abs(P[:, 0] - x) < 1.0) & (np.sign(P[:, 1]) == s)]
        out[i, 1] = sel[0, 1]
    return out


def panel4(ctx, out, scan=True):
    sp, sc, mask, r = ctx['sp'], ctx['sc'], ctx['mask'], ctx['r2']
    crop = (610, 1486, 945, 1935)
    im, bar, k, pl, pr, PL, PR, m, deg, sd = _common(ctx, (1920, 1940), crop, 721, 820, 1495, 1725, 720.0, 819.5)
    pitch = sp.g_p + sp.G
    n = min(5, len(PL), len(PR))
    mL, mR = side_posts(ctx, 'L')[-n:], side_posts(ctx, 'R')[-n:]
    img, msk = _register(k, m, PL, PR, mL, mR, slice(-n, None), slice(-n, None))
    msk = _yfix(ctx, msk)
    fmap, theta, r_al, r_ac = feature_map(k, m, img, msk)
    lines, sds, ns = _rows_lines(im, PL, PR, range(len(PR) - n, len(PR) - 2), k, sp)
    W = widths(lines, k)
    Mw = _mask_widths(ctx, sc.N - 2)
    p_um, p_cols = pitch_fit([pl, pr], k)
    wl = wall_angle(im, range(1500, 1620, 3), 654.0)
    wr = wall_angle(im, range(1500, 1660, 3), 886.0)
    meas = dict(bar_px_aa=bar[0], bar_px_full_pixels=bar[3], px_per_um=k, rotation_deg=deg, rotation_fit_sd_px=sd,
                wall_angles_deg=dict(left=wl[1], right=wr[1]), pitch_um=p_um, pitch_um_by_column=p_cols,
                mask_pitch_um=float(pitch), image_widths_um=W, mask_widths_um_last_row=Mw,
                printed_um=dict(centre_out=52.0, side_out=44.0, gap=19.0),
                feature_fit_residual_um=dict(along=r_al.tolist(), across=r_ac.tolist()), lines_sd_px=sds, lines_n=ns)
    x_last = float(max(mL.max(), mR.max()))
    region = (x_last - 4 * pitch - 70, x_last + 300, -200, 200)
    o = ctx['outlet']
    head = ('%s, Fig. 2B, 4th panel (the outlet end; the panel\'s 50 µm bar, %.2f px; the image\'s rotation %.2f°, from the post '
            'columns) with the edges of %s: the last five post rows and the outlet stretch (%s µm straight from the last gap\'s '
            'centre, then %s° per side).' % (SRC, bar[0], deg, ctx['label'], o.get('flare_start_um'), o.get('flare_deg')))
    tail = ('On the image (dark-line centres, µm at the bar): post pitch %.2f (mask %.1f); post columns %.1f apart (mask %.1f); '
            'outer walls %.1f (mask %.1f); centre %.1f (printed 52); sides %.1f / %.1f (printed 44).' % (
                p_um, pitch, W['post_column_separation'], Mw['post_column_separation'], W['outer_walls'], Mw['outer_walls'],
                W['centre'], W['side_L'], W['side_R']))
    exclude = ((610, 1622, 660, 1660), (690, 1638, 765, 1654), (676, 1680, 752, 1697), (792, 1680, 864, 1697),
               (862, 1668, 945, 1705), (885, 1915, 945, 1935))
    along = (np.cos(theta), np.sin(theta))
    rec, cmap = run_panel('dinh_fig2B_p4_outlet__' + ctx['tag'], out, r, crop, k, mask, region, fmap, along, pitch,
                          exclude=exclude, caption_head=head, caption_tail=tail, bar_um=50, extra=dict(measured=meas),
                          post_max=_pm(sp), metric_search_px=window_px(k, sp.G), gap_um=sp.G)
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, along, exclude=exclude, post_max=_pm(sp))
    return rec


def panel3(ctx, out, p4, scan=True):
    sp, sc, mask, r = ctx['sp'], ctx['sc'], ctx['mask'], ctx['r2']
    crop = (610, 1000, 945, 1452)
    im, bar, k, pl, pr, PL, PR, m, deg, sd = _common(ctx, (1436, 1456), crop, 727, 834, 1010, 1445, 716.5, 831.0)
    pitch = sp.g_p + sp.G
    lines, sds, ns = _rows_lines(im, PL, PR, [0, 1, 6, 7, 8], k, sp)
    W = widths(lines, k)
    x = np.asarray(sc.x, float)
    wc = np.asarray(sc.wc, float)
    ws = np.asarray(sc.ws, float)
    lin = np.arange(sc.NM + sp.n_rect + 1, sc.N - 1)
    # the post columns' separation (centre to centre; across the flow, i.e. corrected for the rotation)
    sep = float(np.mean(PR[:, 0]) - np.mean(PL[:, 0])) * np.cos(np.arctan(m)) / k
    i_sep = int(lin[np.argmin(np.abs(wc[lin] - (sep - sp.l_g)))])
    biasW = p4['measured']['image_widths_um']['outer_walls'] - p4['measured']['mask_widths_um_last_row']['outer_walls']
    Wm = wc + 2 * sp.l_g + 2 * ws
    i_W = int(lin[np.argmin(np.abs(Wm[lin] - (W['outer_walls'] - biasW)))])
    mLa, mRa = side_posts(ctx, 'L'), side_posts(ctx, 'R')
    n = min(len(PL), len(PR))
    j0 = int(np.argmin(np.abs(mRa - x[i_sep])))
    j = int(np.clip(j0 - n // 2, 0, len(mRa) - n))
    mR = mRa[j:j + n]
    mL = np.array([mLa[np.argmin(np.abs(mLa - v))] for v in mR])
    img, msk = _register(k, m, PL, PR, mL, mR, slice(0, n), slice(0, n))
    msk = _yfix(ctx, msk)
    fmap, theta, r_al, r_ac = feature_map(k, m, img, msk)
    i_first = int(np.argmin(np.abs(x - mR[0])))
    Mw = _mask_widths(ctx, i_sep)
    p_um, p_cols = pitch_fit([pl, pr], k)
    wl = wall_angle(im, range(1010, 1440, 3), 661.5)
    wr = wall_angle(im, range(1010, 1440, 3), 885.5)
    meas = dict(bar_px_aa=bar[0], bar_px_full_pixels=bar[3], px_per_um=k, rotation_deg=deg, rotation_fit_sd_px=sd,
                wall_angles_deg=dict(left=wl[1], right=wr[1]),
                pitch_um=p_um, pitch_um_by_column=p_cols, mask_pitch_um=float(pitch), image_widths_um=W,
                post_column_separation_um_from_centres=sep, row_from_post_columns=i_sep, row_from_outer_walls=i_W,
                outer_wall_bias_from_panel4_um=biasW, fraction_along_linear=float((i_sep - sc.NM) / (sc.N - sc.NM)),
                mask_widths_um_at_row=Mw, rows_traced=[i_first, i_first + n - 1],
                feature_fit_residual_um=dict(along=r_al.tolist(), across=r_ac.tolist()), lines_sd_px=sds, lines_n=ns)
    head = ('%s, Fig. 2B, 3rd panel (the linear section; the panel\'s bar, %.2f px per 50 µm; rotation %.2f°) with %s, rows %d–%d.  '
            'The row is found from the post columns\' separation (centre to centre, %.1f µm: w_c = %.1f µm, row %d, %.0f %% along the '
            'linear section); the outer walls, corrected by the 4th panel\'s image-minus-mask offset, give row %d.' % (
                SRC, bar[0], deg, ctx['label'], i_first, i_first + n - 1, sep, sep - sp.l_g, i_sep, 100 * meas['fraction_along_linear'], i_W))
    tail = ('At row %d the mask has w_c %.1f, w_s %.1f, outer walls %.1f µm; the image: centre %.1f, sides %.1f / %.1f, outer walls %.1f, '
            'pitch %.2f µm (mask %.1f).' % (i_sep, Mw['centre'], Mw['side'], Mw['outer_walls'], W['centre'], W['side_L'], W['side_R'],
                                             W['outer_walls'], p_um, pitch))
    region = (mR[0] - pitch, mR[-1] + pitch, -200, 200)
    exclude = ((610, 1168, 662, 1202), (712, 1100, 724, 1265), (885, 1435, 945, 1452))
    along = (np.cos(theta), np.sin(theta))
    rec, cmap = run_panel('dinh_fig2B_p3_linear__' + ctx['tag'], out, r, crop, k, mask, region, fmap, along, pitch,
                          exclude=exclude, caption_head=head, caption_tail=tail, bar_um=50, extra=dict(measured=meas),
                          post_max=_pm(sp), metric_search_px=window_px(k, sp.G), gap_um=sp.G)
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, along, exclude=exclude, post_max=_pm(sp))
    return rec


def panel2(ctx, out, scan=True):
    sp, sc, mask, r = ctx['sp'], ctx['sc'], ctx['mask'], ctx['r2']
    crop = (610, 522, 945, 985)
    im, bar, k, pl, pr, PL, PR, m, deg, sd = _common(ctx, (952, 972), crop, 717, 837, 735, 980, 717.0, 837.0)
    pitch = sp.g_p + sp.G
    x = np.asarray(sc.x, float)
    mLa, mRa = side_posts(ctx, 'L'), side_posts(ctx, 'R')
    mRa = mRa[mRa > x[sc.NM] - 1.0]
    mLa = mLa[mLa > x[sc.NM] - 1.0]
    n = min(len(PL), len(PR), 5)
    img, msk = _register(k, m, PL, PR, mLa[:n], mRa[:n], slice(0, n), slice(0, n))
    msk = _yfix(ctx, msk)
    fmap, theta, r_al, r_ac = feature_map(k, m, img, msk)
    lines, sds, ns = _rows_lines(im, PL, PR, [2, 3, 4], k, sp)
    W = widths(lines, k)
    i1 = int(np.argmin(np.abs(x - mRa[0])))
    Mw = _mask_widths(ctx, i1 + 2)
    p_um, p_cols = pitch_fit([pl, pr], k)
    meas = dict(bar_px_aa=bar[0], bar_px_full_pixels=bar[3], px_per_um=k, rotation_deg=deg, rotation_fit_sd_px=sd,
                pitch_um=p_um, pitch_um_by_column=p_cols, first_post_row=i1, image_widths_um=W,
                mask_widths_um_third_post_row=Mw, printed_um=dict(centre_in=75.0, side_in=21.0),
                feature_fit_residual_um=dict(along=r_al.tolist(), across=r_ac.tolist()), lines_sd_px=sds, lines_n=ns)
    head = ('%s, Fig. 2B, 2nd panel (the head-to-array transition; the panel\'s bar, %.2f px per 50 µm; rotation %.2f°) with %s '
            'around the seam: the rectilinear rows and the first post rows (registered on the first %d posts of each column).' % (
                SRC, bar[0], deg, ctx['label'], n))
    tail = ('On the image at the 3rd–5th post rows: centre %.1f µm (printed 75), sides %.1f / %.1f (printed 21), outer walls %.1f '
            '(mask %.1f), post columns %.1f apart (mask %.1f).' % (W['centre'], W['side_L'], W['side_R'], W['outer_walls'],
                                                                    Mw['outer_walls'], W['post_column_separation'], Mw['post_column_separation']))
    region = (mRa[0] - 300, mRa[0] + 300, -200, 200)
    exclude = ((610, 776, 672, 796), (688, 776, 762, 796), (672, 748, 744, 768), (806, 748, 882, 768), (885, 738, 945, 776),
               (885, 965, 945, 985))
    along = (np.cos(theta), np.sin(theta))
    rec, cmap = run_panel('dinh_fig2B_p2_transition__' + ctx['tag'], out, r, crop, k, mask, region, fmap, along, pitch,
                          exclude=exclude, caption_head=head, caption_tail=tail, bar_um=50, extra=dict(measured=meas),
                          post_max=_pm(sp), metric_search_px=window_px(k, sp.G), gap_um=sp.G)
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, along, exclude=exclude, post_max=_pm(sp))
    return rec


def panel1(ctx, out, scan=True):
    """the head's start: the inlet funnel's end, the stem and the first arm junctions.  Across: the centre channel's
    two walls (rotation from their slope); along: the first side channel's upstream wall (image row of the dark line
    at the junction, both sides) against the mask's first channel edge across the flow; chamfer over +-0.6 unit pitch"""
    sp, sc, mask, r = ctx['sp'], ctx['sc'], ctx['mask'], ctx['r2']
    im = r.im
    bar = T.bar_length_aa(im, 470, 490, 870, 975)
    k = ctx.get('k_override') or bar[0] / 50.0
    crop = (610, 34, 945, 495)
    wl = wall_angle(im, range(140, 266, 2), 736.5, c0=700, c1=840)
    wr = wall_angle(im, range(140, 266, 2), 802.5, c0=700, c1=840)
    m = 0.5 * (wl[0] + wr[0])
    rows = range(140, 266, 2)
    Lw, Rw = [], []
    for rr in rows:
        L = T.dark_lines(T.row_profile(im, rr, 700, 840, 1), 10, off=700)
        a, b = T.line_near(L, 736.5, 3.0), T.line_near(L, 802.5, 3.0)
        if a is not None and b is not None:
            Lw.append((a, rr)); Rw.append((b, rr))
    Lw, Rw = np.array(Lw), np.array(Rw)
    rmid = float(np.mean(Lw[:, 1]))
    axis = float(0.5 * (np.mean(Lw[:, 0]) + np.mean(Rw[:, 0])))
    stem_um = float(np.mean(Rw[:, 0] - Lw[:, 0]) * np.cos(np.arctan(m)) / k)
    # the first junction's upstream wall: the dark line at rows ~270-292 in four columns beside the centre channel
    jr = []
    for c in (670, 690, 850, 870):
        L = T.dark_lines(T.col_profile(im, c, 260, 300, 2), 10, off=260)
        if len(L):
            jr.append(float(L[0]))
    row_j = float(np.mean(jr))
    # the mask's first channel edge across the flow beside the centre channel (the channel polygon's edge at the
    # smallest x > -100 that spans y 50..200)
    xj = None
    for e in mask.edges(-150, 150, 40, 260):
        d = np.diff(e, axis=0)
        for p, q, dd in zip(e[:-1], e[1:], d):
            if e.tag == 'channel' and abs(dd[0]) < 1e-6 and abs(dd[1]) > 50:
                xj = p[0] if xj is None else min(xj, p[0])
    theta = float(np.arctan2(1.0, m))
    M0 = T.rotated_map(k, theta)
    # origin: the image point (axis at row_j, corrected for the rotation) <- the mask point (xj, 0)
    img = np.array([axis + m * (row_j - rmid), row_j])
    O = img - M0(np.array([[xj, 0.0]]))[0]
    fmap = T.rotated_map(k, theta, O)
    meas = dict(bar_px_aa=bar[0], bar_px_full_pixels=bar[3], px_per_um=k, rotation_deg=float(np.degrees(np.arctan(m))),
                wall_angles_deg=dict(left=wl[1], right=wr[1]), stem_width_um=stem_um, first_junction_row_px=row_j,
                first_junction_rows_px=jr, mask_first_channel_edge_x_um=float(xj),
                note="along: one anchor (the first junction's upstream wall); the funnel, the stem length and the junction "
                     "geometry are the depicted head's (Fig. 2A inset unit lengths) and the inlet style's")
    head = ('%s, Fig. 2B, 1st panel (the head\'s start; the panel\'s bar, %.2f px per 50 µm; rotation %.2f°, from the centre '
            'channel\'s walls) with %s: the inlet funnel\'s end, the stem and the first arm junctions of the depicted head.' % (
                SRC, bar[0], meas['rotation_deg'], ctx['label']))
    tail = ('Registered across on the stem\'s walls (%.1f µm apart on the image) and along on the first junction\'s upstream wall; '
            'the funnel and the stem length are the inlet\'s (dinh_filter), the junctions the depicted head\'s.' % stem_um)
    region = (xj - 300, xj + 300, -200, 200)
    exclude = ((872, 130, 945, 170), (900, 205, 915, 282), (900, 300, 915, 378), (664, 316, 738, 330), (804, 316, 878, 330),
               (885, 470, 945, 495))
    along = (np.cos(theta), np.sin(theta))
    rec, cmap = run_panel('dinh_fig2B_p1_head__' + ctx['tag'], out, r, crop, k, mask, region, fmap, along, 110.0,
                          exclude=exclude, caption_head=head, caption_tail=tail, bar_um=50, extra=dict(measured=meas),
                          post_max=_pm(sp), metric_search_px=window_px(k, sp.G), gap_um=sp.G)
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, along, exclude=exclude, post_max=_pm(sp))
    return rec
