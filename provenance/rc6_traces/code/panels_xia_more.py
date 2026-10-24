"""Xia et al. 2016, Fig. 3 panels c, d, e and f traced with an opencif Xia mask (rc6): the post rows of the legs at
four places of the device, beside panel g (the outlet, panels_xia.panel_g).

  c  leg 1 where the papers' rectilinear transition (square wall blocks beside a GS-wide side channel) gives way to the
     first post rows; the first pill below the last block is the mask's first post of leg 1 (a non-periodic feature
     fixes the along-flow phase: the blocks' pitch is 117 µm, the posts' 51.8 µm);
  d  legs 5 (left) and 3 (right) where they leave the loops at the device's top; the first pill below each loop is that
     leg's first post (one map for both legs: their spacing, 905 µm in the layout, is tested, not fitted);
  e  legs 5 and 3 half way down: a periodic stretch (rows of identical posts), so no feature of the panel fixes which
     rows it shows.  It is placed where Fig. 3a's box e puts it (xia_fig3a_device.py: the box's centre mapped into the
     mask by the drawing's registration at the device's top); the image's middle pill on each side is matched to the
     post nearest that place and the others by index.  The along-flow phase is therefore the drawing's, to within its
     registration; the panel tests the local geometry (pitch, posts, centre, side channels, the two legs' spacing).
     The centre width read on the image is compared with the mask's at that place, and the place where the mask's
     taper has the image's width is reported (the taper's form is assumed, A);
  f  leg 7 a quarter of the way down: the same as e, on one leg.  Panel f is reproduced at its own scale: its post
     columns' period is 29.19 px against 28.75-28.81 px on panels e and g (a Fourier peak over each column's full
     height), so the figure enlarged it about 1.5 % more than the others; it is traced at the pitch measured on the
     other panels (51.78 µm, data/r6/rc6_values.json) over its own period, which leaves its pitch untested and its
     widths tested.
Scale: the rc6 scale (data/r6/rc6_values.json, 1.798 µm/px).  Rotation: the post columns' common slope.  Translation:
the matched posts (features) and, independently, chamfer matching (panel_core.run_panel), with V1-V5 as for every
panel.  Where each panel lies on the mask is also compared with Fig. 3a's box for it (xia_fig3a_device.py, V12)."""
import os

import numpy as np
from scipy import ndimage as ndi

import trace_lib as T
from panel_core import run_panel, scale_scan
from panels_dinh import fit_rotation, feature_map

SRC = 'Xia et al. 2016 Sci. Rep. 6:35943'
FIG3_SHA = 'e9fb7fbb348567d3'
CROPS = dict(c=(1277, 1, 2025, 748), d=(1, 796, 748, 1542), e=(1277, 797, 2025, 1543), f=(1, 1594, 748, 2340))
# per panel: the legs shown (1-based), each with its two post columns' image columns (approximate, raster px), the rows
# holding pills, how the along-flow phase is fixed, the ink to leave out (panel letter, arrow), and what it shows
PANELS = dict(
    c=dict(legs=[(1, (1519.5, 1727.0))], rows=(430, 745), anchor='first',
           exclude=((1280, 4, 1345, 80), (1280, 290, 1345, 430)),
           what='leg 1 where the rectilinear transition\'s wall blocks give way to the first post rows'),
    d=dict(legs=[(5, (70.0, 190.0)), (3, (552.0, 712.0))], rows=(1340, 1540), anchor='first',
           exclude=((2, 800, 65, 870), (355, 1170, 420, 1325)),
           what='legs 5 and 3 leaving the loops at the device\'s top'),
    e=dict(legs=[(5, (1337.0, 1462.0)), (3, (1830.0, 1972.0))], rows=(800, 1540), anchor='box',
           exclude=((1280, 800, 1345, 872), (1630, 1150, 1700, 1295)),
           what='legs 5 and 3 half way down: a periodic stretch, placed at Fig. 3a\'s box'),
    f=dict(legs=[(7, (311.0, 436.0))], rows=(1597, 2338), anchor='box', scale='own_pitch',
           exclude=((2, 1597, 60, 1670), (2, 1870, 62, 2010)),
           what='leg 7 a quarter of the way down: a periodic stretch, placed at Fig. 3a\'s box'))


def fourier_period(im, col, r0, r1, lo=26.0, hi=32.0):
    """the post columns' period (px): the Fourier peak of the column's profile over its full height"""
    p = im[r0:r1, int(col) - 3:int(col) + 4].mean(axis=1)
    p = p - p.mean()
    n = np.arange(len(p))
    pers = np.arange(lo, hi, 0.005)
    c = [abs(np.sum(p * np.exp(-2j * np.pi * n / q))) for q in pers]
    return float(pers[int(np.argmax(c))])


def pills(im, col, r0, r1, k, sp):
    """pill centres (rows) along one post column: the pairs of dark outline lines g_p apart"""
    p = ndi.gaussian_filter1d(T.col_profile(im, col, r0, r1, 2), 0.8)
    L = T.dark_lines(p, 6, halfwin=5, off=r0)
    return T.pair_centres(L, 0.75 * sp.g_p * k, 1.25 * sp.g_p * k, skip=1)


def centres(im, rows, guess, hp, c0, c1):
    """each pill's centre column (the mid-point of its two outer ends, 3 rows above and below its centre row)"""
    pts, g = [], guess
    for rr in rows:
        cs = []
        for dr in (-3, 3):
            L = T.dark_lines(T.row_profile(im, int(round(rr)) + dr, c0, c1, 1), 6, halfwin=5, off=c0)
            a, b = T.line_near(L, g - hp, 6.0), T.line_near(L, g + hp, 6.0)
            if a is not None and b is not None:
                cs.append(0.5 * (a + b))
        if cs:
            pts.append((float(np.mean(cs)), float(rr), 0.0))
            g = float(np.mean(cs))
    return np.array(pts)


def leg_posts(mask, sp, Y):
    """the leg's posts on each side of its axis: (left, right) arrays of (x, y), sorted along the flow"""
    P = mask.posts(-1e9, 1e9, Y - 450.0, Y + 450.0, max_len=sp.g_p + 1.0, max_wid=sp.l_g + 1.0)
    return P[P[:, 1] < Y], P[P[:, 1] > Y]


def mask_centre_width(Lp, Rp, sp):
    """the drawn centre width along the leg: each left post paired with the right post nearest in x"""
    xs, ws = [], []
    for x, y in Lp:
        j = int(np.argmin(np.abs(Rp[:, 0] - x)))
        if abs(Rp[j, 0] - x) < 1.0:
            xs.append(x)
            ws.append(Rp[j, 1] - y - sp.l_g)
    return np.array(xs), np.array(ws)


def panel(letter, mask, sp, rec, out, k, tag, label, scale_note, rasters, scan=True, fig3a=None):
    """fig3a: the xia_fig3a_device.py result for this mask (needed by the panels placed at Fig. 3a's box)"""
    cfg = PANELS[letter]
    r = T.Raster(os.path.join(rasters, 'xia5-000.png'), FIG3_SHA)
    im = r.im
    crop = CROPS[letter]
    pitch = sp.g_p + sp.G
    hp = 0.5 * sp.l_g * k
    r0, r1 = cfg['rows']
    legY = rec['design']['leg_Y']
    if cfg.get('scale') == 'own_pitch':
        import json as _json
        lam = _json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))), 'opencif', 'data', 'r6', 'rc6_values.json')))['xia']['row_pitch_um']
        per_f = float(np.mean([fourier_period(im, c_, r0, r1) for _, cc in cfg['legs'] for c_ in cc]))
        k_rc6 = k
        k = per_f / lam
        scale_note = ('its own scale %.4f µm/px: the pitch measured on panels c and g, %.2f µm, over its own period '
                      '%.3f px (the rc6 scale of the other panels is %.3f)' % (1.0 / k, lam, per_f, 1.0 / k_rc6))
        hp = 0.5 * sp.l_g * k
    cols, img, msk, per_leg, tips = [], [], [], [], []
    for leg, (cl, cr) in cfg['legs']:
        Y = float(legY[leg - 1])
        Lp, Rp = leg_posts(mask, sp, Y)
        rows_l, rows_r = pills(im, cl, r0, r1, k, sp), pills(im, cr, r0, r1, k, sp)
        PL = centres(im, rows_l, cl, hp, int(cl - 3 * hp), int(cl + 3 * hp))
        PR = centres(im, rows_r, cr, hp, int(cr - 3 * hp), int(cr + 3 * hp))
        cols += [PL, PR]
        d = np.r_[np.diff(PL[:, 1]), np.diff(PR[:, 1])]
        per = float(np.median(d[(d > 0.7 * pitch * k) & (d < 1.3 * pitch * k)]))
        sep_px = float(np.median(PR[:, 0]) - np.median(PL[:, 0]))
        w_img = sep_px / k - sp.l_g                     # the centre width at the rc6 scale (column separation - l_g)
        xs_w, ws_w = mask_centre_width(Lp, Rp, sp)
        if cfg['anchor'] == 'first':
            jL = np.round((PL[:, 1] - PL[0, 1]) / per).astype(int)
            jR = np.round((PR[:, 1] - PR[0, 1]) / per).astype(int)
            pairsL = [((c, rr), (Lp[j, 0], Lp[j, 1])) for (c, rr, _), j in zip(PL, jL) if j < len(Lp)]
            pairsR = [((c, rr), (Rp[j, 0], Rp[j, 1])) for (c, rr, _), j in zip(PR, jR) if j < len(Rp)]
            located = dict(by='the first pill below the non-periodic feature = the leg\'s first post')
        else:
            x_box = float(fig3a['V12_boxes'][letter]['box_centre_mask_um'][0])
            # where the mask's taper has the image's width (reported, not used)
            o = np.argsort(xs_w)
            xs_o, ws_o = xs_w[o], ws_w[o]
            dirn = np.sign(ws_o[-1] - ws_o[0]) or 1.0
            x_mid = float(np.interp(w_img * dirn, ws_o * dirn, xs_o)) if (ws_o * dirn).min() <= w_img * dirn <= (ws_o * dirn).max() \
                else float(xs_o[int(np.argmin(np.abs(ws_o - w_img)))])
            # the image's pill nearest the panel's centre row on each side <-> the post nearest Fig. 3a's box centre
            pairsL, pairsR = [], []
            rc = 0.5 * (crop[1] + crop[3])
            for PC, MP, pr in ((PL, Lp, pairsL), (PR, Rp, pairsR)):
                m = int(np.argmin(np.abs(PC[:, 1] - rc)))
                j0 = int(np.argmin(np.abs(MP[:, 0] - x_box)))
                jj = np.round((PC[:, 1] - PC[m, 1]) / per).astype(int) + j0
                pr += [((c, rr), (MP[j, 0], MP[j, 1])) for (c, rr, _), j in zip(PC, jj) if 0 <= j < len(MP)]
            w_box = float(np.interp(x_box, xs_o, ws_o))
            located = dict(by='Fig. 3a\'s box centre (mask x %.0f µm); the mask\'s centre there %.1f µm, the image\'s %.1f µm'
                              % (x_box, w_box, w_img), x_box_um=x_box, w_mask_at_box_um=w_box, w_image_um=w_img,
                           x_where_mask_taper_has_image_width_um=x_mid,
                           w_range_on_leg_um=[float(ws_w.min()), float(ws_w.max())])
        for a, b in pairsL + pairsR:
            img.append(a)
            msk.append(b)
        per_leg.append(dict(leg=leg, Y_um=Y, period_px=per, pitch_image_um=per / k, column_separation_px=sep_px,
                            centre_um_image=w_img, located=located, n_pills=[len(PL), len(PR)],
                            first_post_x_um=[float(Lp[0, 0]), float(Rp[0, 0])]))
        tips.append(Y)
    img, msk = np.array(img), np.array(msk)
    m, deg, sd = fit_rotation(cols)
    fmap, theta, r_al, r_ac = feature_map(k, m, img, msk)
    # the traced region: the whole panel mapped into the mask (its four corners through the feature map), plus 60 µm
    corners = fmap.inverse(np.array([[crop[0], crop[1]], [crop[2], crop[1]], [crop[0], crop[3]], [crop[2], crop[3]]], float))
    region = (corners[:, 0].min() - 60.0, corners[:, 0].max() + 60.0, corners[:, 1].min() - 60.0, corners[:, 1].max() + 60.0)
    meas = dict(px_per_um=k, scale_note=scale_note, rotation_deg=deg, rotation_fit_sd_px=sd, legs=per_leg,
                pitch_um=dict(image=float(np.mean([q['pitch_image_um'] for q in per_leg])), mask=float(pitch)),
                feature_fit_residual_um=dict(along=r_al.tolist(), across=r_ac.tolist()),
                mask=dict(post_across=float(sp.l_g), post_along=float(sp.g_p)))
    if len(per_leg) == 2:
        meas['leg_axis_separation_um'] = dict(
            image=float((np.mean(cols[2][:, 0]) + np.mean(cols[3][:, 0]) - np.mean(cols[0][:, 0]) - np.mean(cols[1][:, 0]))
                        / 2.0 * np.cos(np.arctan(m)) / k),
            mask=float(abs(per_leg[1]['Y_um'] - per_leg[0]['Y_um'])))
    head = '%s, Fig. 3%s (%s) at %.3f µm/px (%s; rotation %.2f°) with %s.' % (SRC, letter, cfg['what'], 1.0 / k,
                                                                           scale_note, deg, label)
    tail = 'On the image at this scale: pitch %.2f µm (mask %.2f); ' % (meas['pitch_um']['image'], pitch) + '; '.join(
        'leg %d centre %.1f µm (%s)' % (q['leg'], q['centre_um_image'], 'first post below the feature' if cfg['anchor'] == 'first'
                                       else 'mask %.1f µm at Fig. 3a\'s box' % q['located']['w_mask_at_box_um'])
        for q in per_leg)
    if 'leg_axis_separation_um' in meas:
        tail += '; leg axes %.1f µm apart (mask %.1f)' % (meas['leg_axis_separation_um']['image'],
                                                          meas['leg_axis_separation_um']['mask'])
    tail += '.'
    along = (np.cos(theta), np.sin(theta))
    recp, cmap = run_panel('xia_fig3%s__%s' % (letter, tag), out, r, crop, k, mask, region, fmap, along, pitch,
                           exclude=cfg['exclude'], caption_head=head, caption_tail=tail, bar_um=250,
                           extra=dict(measured=meas, panel=letter),
                           scale_src=('scale from its own post period at the other panels\' pitch'
                                      if cfg.get('scale') == 'own_pitch' else
                                      'scale from the printed widths, column separations and loop curvature'),
                           post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=10)
    if scan:
        recp['V5_scale_scan'] = scale_scan(r, crop, mask, region, cmap, along, exclude=cfg['exclude'],
                                           post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=10)
    return recp
