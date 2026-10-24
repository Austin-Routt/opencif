"""Xia et al. 2016 Sci. Rep. 6:35943, Fig. 3g (the device's outlet: leg 7's last post rows and the leads), traced with
an opencif Xia mask at a declared scale.

Scale: there is one bar for panels b-i, drawn on panel i (250 µm: 128 full pixels + its anti-aliased ends = 129.11 px,
1.936 µm/px; rc5 and Reviewer A read 128 px, 1.953).  The panels' printed widths give another: the centre channel
reads 164.3 px at w_c(0) = 300 µm (panel c) and 80.3 px at the outlet's 150 µm (panel g), so k = (164.3 - 80.3) /
150 (a difference of two wide channels of the same kind: free of the outline bias), and Fig. 3a's leg spacing (its own
5 mm bar) against panel d's agrees with it (scale_evidence()).  The trace is run at the scale it declares.
Rotation: from the post columns.  Translation: features (leg 7's last posts) and chamfer."""
import os

import numpy as np
from scipy import ndimage as ndi

import trace_lib as T
from panel_core import run_panel, scale_scan
from panels_dinh import fit_rotation, feature_map

SRC = 'Xia et al. 2016 Sci. Rep. 6:35943'
FIG3_SHA = 'e9fb7fbb348567d3'
PANELS = {'c': (1277, 1, 2025, 748), 'd': (1, 800, 747, 1538), 'g': (1277, 1594, 2025, 2338)}


def bar_scale(im):
    b = T.bar_length_aa(im, 3085, 3120, 1860, 2025, dark=False)
    return b[0] / 250.0, b


def centre_rows(im, colL, colR, r0, r1, half_post, c0, c1, offs=(-3, 3)):
    """post-end lines at the pill rows of the right column: (row, Lo, Li, Ri, Ro) in px"""
    p = ndi.gaussian_filter1d(T.col_profile(im, colR, r0, r1, 3), 2.0)
    rows = [r0 + i for i in range(4, len(p) - 4) if p[i] == p[i - 4:i + 5].max()]
    out = []
    for r in rows:
        acc = []
        for dr in offs:
            L = T.dark_lines(T.row_profile(im, r + dr, c0, c1, 1), 6, halfwin=5, off=c0)
            v = [T.line_near(L, x, 6.0) for x in (colL - half_post, colL + half_post, colR - half_post, colR + half_post)]
            if None not in v:
                acc.append(v)
        if acc:
            out.append([r] + list(np.mean(acc, axis=0)))
    return np.array(out)


def scale_evidence(im):
    """the printed-width scale: the centre channel (inner post ends) at the transition's first post rows (panel c,
    w_c(0) = 300 µm, Xia Fig. 1b(i) and p. 3) and at the outlet (panel g, 150 µm, Fig. 1d(i)); and panel d's leg axes"""
    c = centre_rows(im, 1519.5, 1725.5, 450, 745, 21, 1400, 1850)
    g = centre_rows(im, 1769.4, 1891.6, 1600, 1975, 21, 1700, 1950)
    Wc = float(np.mean(c[:, 3] - c[:, 2]))
    Wg = float(np.mean(g[:, 3] - g[:, 2]))
    k = (Wc - Wg) / 150.0
    b = 0.5 * (Wg - 150.0 * k)
    # panel d: the two legs' axes (mid-points of their post columns) at the bottom rows
    axes = []
    for (cl, cr) in ((70.0, 190.0), (552.0, 715.0)):
        d = centre_rows(im, cl, cr, 1400, 1535, 21, int(max(1, cl - 60)), int(min(747, cr + 60)))
        axes.append(float(np.mean(0.25 * (d[:, 1] + d[:, 2] + d[:, 3] + d[:, 4]))) if len(d) else None)
    return dict(centre_px_panel_c=Wc, centre_px_panel_g=Wg, rows_c=len(c), rows_g=len(g),
                px_per_um_printed=k, um_per_px_printed=1.0 / k, outline_bias_px=b,
                panel_d_leg_axes_px=axes, panel_d_axis_separation_px=(axes[1] - axes[0]) if None not in axes else None,
                fig3a_leg_spacing_um_legs_3_5=901.0, fig3a_leg_spacing_uncertainty_um=28.0)


def panel_g(mask, sp, rec, out, k, tag, label, scale_note, scan=True):
    r = T.Raster(os.path.join(os.environ['XIA_RASTERS'], 'xia5-000.png'), FIG3_SHA)
    im = r.im
    crop = PANELS['g']
    pitch = sp.g_p + sp.G
    # pill centres along the two post columns
    def pills(col):
        p = ndi.gaussian_filter1d(T.col_profile(im, col, 1600, 1990, 2), 0.8)
        L = T.dark_lines(p, 6, halfwin=5, off=1600)
        return T.pair_centres(L, 0.75 * sp.g_p * k, 1.25 * sp.g_p * k, skip=1)
    pl, pr = pills(1769.4), pills(1891.6)
    hp = 0.5 * sp.l_g * k
    def centres(rows, guess):
        pts = []
        g = guess
        for rr in rows:
            cs = []
            for dr in (-3, 3):
                L = T.dark_lines(T.row_profile(im, int(round(rr)) + dr, 1700, 1950, 1), 6, halfwin=5, off=1700)
                a, b = T.line_near(L, g - hp, 6.0), T.line_near(L, g + hp, 6.0)
                if a is not None and b is not None:
                    cs.append(0.5 * (a + b))
            if cs:
                pts.append((float(np.mean(cs)), float(rr), 0.0))
                g = float(np.mean(cs))
        return np.array(pts)
    PL, PR = centres(pl, 1769.4), centres(pr, 1891.6)
    m, deg, sd = fit_rotation([PL, PR])
    Y7 = rec['design']['leg_Y'][6]
    P = mask.posts(-1e9, 1e9, Y7 - 300, Y7 + 300, max_len=sp.g_p + 1.0, max_wid=sp.l_g + 1.0)
    yL, yR = Y7 - 0.5 * (sp.wc_min + sp.l_g), Y7 + 0.5 * (sp.wc_min + sp.l_g)
    mL = np.sort(P[np.abs(P[:, 1] - yL) < 5, 0])
    mR = np.sort(P[np.abs(P[:, 1] - yR) < 5, 0])
    # index the pills from the last one by the image's own period (a pill whose outline the fluorescence hides is
    # skipped, not mis-assigned), then match them to the mask's posts counted back from its last row
    d = np.diff(np.r_[pl, pr]) if False else np.r_[np.diff(pl), np.diff(pr)]
    per = float(np.median(d[(d > 20) & (d < 36)]))
    img, msk = [], []
    for Pc, ms, yy in ((PL, mL, yL), (PR, mR, yR)):
        idx = np.round((Pc[:, 1] - Pc[-1, 1]) / per).astype(int)
        for (c, rr, _), j in zip(Pc, idx):
            if -12 < j <= 0:
                img.append((c, rr))
                msk.append((ms[len(ms) - 1 + j], yy))
    img, msk = np.array(img), np.array(msk)
    fmap, theta, r_al, r_ac = feature_map(k, m, img, msk)
    E = centre_rows(im, 1769.4, 1891.6, 1600, 1975, hp, 1700, 1950)
    # the left side channel: from the left outer wall (the dark line near col 1402) to the left post's outer end
    wl = []
    for row in E[:, 0]:
        L = T.dark_lines(T.row_profile(im, int(row), 1380, 1430, 1), 6, halfwin=5, off=1380)
        v = T.line_near(L, 1402.0, 6.0)
        if v is not None:
            wl.append(v)
    cosr = np.cos(np.arctan(m))
    meas = dict(px_per_um=k, scale_note=scale_note, rotation_deg=deg, rotation_fit_sd_px=sd,
                pitch_um=dict(image=per / k, image_px=per, mask=float(pitch)),
                centre_um=float(np.mean(E[:, 3] - E[:, 2]) * cosr / k), posts_um=[float(np.mean(E[:, 2] - E[:, 1]) * cosr / k), float(np.mean(E[:, 4] - E[:, 3]) * cosr / k)],
                column_separation_um=float(np.mean(PR[:, 0]) - np.mean(PL[:, 0])) * cosr / k,
                left_side_channel_um=float((np.mean(E[:, 1]) - np.mean(wl)) * cosr / k) if wl else None,
                mask=dict(centre=float(sp.wc_min), post=float(sp.l_g), column_separation=float(sp.wc_min + sp.l_g),
                          side_last_row=float(rec['design'].get('ws_last_um', np.nan)) if 'ws_last_um' in rec['design'] else None),
                feature_fit_residual_um=dict(along=r_al.tolist(), across=r_ac.tolist()))
    x_last = float(max(mL.max(), mR.max()))
    region = (x_last - 760, x_last + 700, Y7 - 820, Y7 + 380)
    head = ('%s, Fig. 3g (the outlet: leg 7\'s last post rows and the leads) at %.3f µm/px (%s; rotation %.2f°) with %s.' % (
        SRC, 1.0 / k, scale_note, deg, label))
    tail = ('On the image at this scale: centre %.1f µm (printed 150), posts %.1f / %.1f (mask %.1f), post columns %.1f apart (mask %.1f), '
            'pitch %.2f (mask %.2f), left side channel %s µm.' % (
                meas['centre_um'], meas['posts_um'][0], meas['posts_um'][1], sp.l_g, meas['column_separation_um'], sp.wc_min + sp.l_g,
                meas['pitch_um']['image'], pitch, ('%.0f' % meas['left_side_channel_um']) if meas['left_side_channel_um'] else 'n/a'))
    exclude = ((1285, 1600, 1345, 1680), (1285, 1870, 1345, 2010))
    along = (np.cos(theta), np.sin(theta))
    recp, cmap = run_panel('xia_fig3g_outlet__' + tag, out, r, crop, k, mask, region, fmap, along, pitch, exclude=exclude,
                           caption_head=head, caption_tail=tail, bar_um=250, extra=dict(measured=meas),
                           scale_src=('scale from the printed widths, column separations and loop curvature'
                                      if 'rc6 scale' in scale_note else 'scale from the printed widths'
                                      if 'printed' in scale_note else 'scale from the bar'),
                           post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=10)
    if scan:
        recp['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, along, exclude=exclude, post_max=(sp.g_p + 1.0, sp.l_g + 1.0),
                                           min_component=10)
    return recp
