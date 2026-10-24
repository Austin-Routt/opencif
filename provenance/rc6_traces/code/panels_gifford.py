"""Gifford et al. 2014 Lab Chip 14, 4496-4505, Fig. 4 inset (300 µm bar): the end of leg 4 (its last post rows, side
channels, outer walls and the leads beyond) beside the neighbouring leg, traced with an opencif Gifford mask.

Scale: the bar read with its anti-aliased ends (107 full pixels + 0.11 + 0.42 = 107.52 px; rc5 and Reviewer A: 107).
Rotation: the leg direction measured on the image (the principal direction of each post column's pill centroids,
averaged over the four columns).  Translation: by features (leg 4's last post rows matched to the mask's last rows)
and, independently, by chamfer matching.  Leg 4 is the reproduction claim.  The neighbouring leg (leg 2) is placed by
the mask's leg spacing: Gifford prints none (his Fig. 4 schematic is 'not drawn to scale', p. 4502), so rc6's presets
take the spacing read on this inset (legs 2 and 4 817.9 um apart; provenance/rc6_scripts/apply_rc6_layout_presets.py).
Leg 2 stays in the no-claim colour either way: with the inset's spacing its centre line is where the image's is (by
construction), but its side channels' widths there depend on which of leg 2's rows pass the inset, which the inset does
not fix; their outer walls' offsets are measured and recorded (measured.leg2_side_walls), not claimed."""
import json
import os

import numpy as np
from scipy import ndimage as ndi

import trace_lib as T
from panel_core import run_panel, scale_scan

SRC = 'Gifford et al. 2014 Lab Chip 14, 4496–4505'
INSET_SHA = 'e70a8b38b87f9639'


def pill_centroids(im, box, smooth=8.0, thr=-12.0, area=(25, 250), size=(5, 25)):
    """the centroids of the closed dark rings' interiors (a pill's bright inside), in image px"""
    c0, r0, c1, r1 = box
    a = im[r0:r1, c0:c1]
    d = a - ndi.gaussian_filter(a, smooth)
    dark = ndi.binary_closing(d < thr, iterations=1)
    inter = ndi.binary_fill_holes(dark) & ~dark
    lab, n = ndi.label(inter)
    out = []
    for i, sl in enumerate(ndi.find_objects(lab)):
        m = lab[sl] == i + 1
        A = int(m.sum())
        h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        if area[0] <= A <= area[1] and size[0] <= h <= size[1] and size[0] <= w <= size[1]:
            yy, xx = np.nonzero(m)
            out.append((c0 + sl[1].start + xx.mean(), r0 + sl[0].start + yy.mean(), A))
    return np.array(out)


def columns(C, u, gap=20.0):
    """split the centroids into post columns by their coordinate across u (unit vector along the legs)"""
    v = np.array([-u[1], u[0]])
    pc = C[:, :2] @ v
    o = np.argsort(pc)
    groups, cur = [], [o[0]]
    for a, b in zip(o[:-1], o[1:]):
        if pc[b] - pc[a] > gap:
            groups.append(cur)
            cur = []
        cur.append(b)
    groups.append(cur)
    return [C[g] for g in groups if len(g) >= 5], v


def direction(cols, u0):
    """the leg direction: each column's principal direction (oriented like u0), weighted by its number of posts"""
    acc = np.zeros(2)
    for P in cols:
        X = P[:, :2] - P[:, :2].mean(axis=0)
        w, V = np.linalg.eigh(X.T @ X)
        d = V[:, -1]
        if d @ u0 < 0:
            d = -d
        acc += len(P) * d
    return acc / np.hypot(*acc)


def run(root, rasters, out, dxf=None, record=None, preset=None, tag='design', label=None, scan=True):
    import sys
    for p in (root,):
        if p not in sys.path:
            sys.path.insert(0, p)
    from opencif import spec as SP
    from opencif.physics import schedule as CS
    name = 'gifford2014_fig4_design'
    dxf = dxf or os.path.join(root, 'runs', name, name + '.dxf')
    record = record or os.path.join(root, 'runs', name, name + '.json')
    preset = preset or os.path.join(root, 'presets', name + '.yaml')
    rec0 = json.load(open(record))
    mask = T.Mask(dxf)
    if mask.sha256 != rec0['files']['dxf_sha256']:
        raise ValueError('the DXF is not its record\'s')
    sp, op = SP.to_engine(SP.load(preset))
    sc = CS.plan(sp)
    d = rec0['design']
    legY = d['leg_Y']
    Y4, Y2 = legY[3], legY[1]
    r = T.Raster(os.path.join(rasters, 'gif5-001.png'), INSET_SHA)
    im = r.im
    bar = T.bar_length_aa(im, 550, 600, 690, 829)
    k = bar[0] / 300.0
    box = T.micrograph_box(im, (190, 40, 829, 560))
    C = pill_centroids(im, box)
    u0 = np.array([1.0, 1.0]) / np.sqrt(2)                   # mask +x (against leg 4's flow) points down-right
    cols, v = columns(C, u0)
    u = direction(cols, u0)
    cols, v = columns(C, u)
    theta = float(np.arctan2(u[1], u[0]))
    # the four columns in order across, along the image direction of mask +y: leg 2's two, then leg 4's two
    vimg = np.array([u[1], -u[0]])
    order = np.argsort([np.mean(P[:, :2] @ vimg) for P in cols])
    cols = [cols[i] for i in order]
    if len(cols) != 4:
        raise ValueError('expected four post columns in the inset, found %d' % len(cols))
    pitch = sp.g_p + sp.G
    pitch_px = pitch * k
    # leg 4's columns: index each pill along the leg from the last (upper-left, smallest along-coordinate) one
    P = mask.posts(-1e9, 1e9, Y4 - 150, Y4 + 150, max_len=sp.g_p + 1.0, max_wid=sp.l_g + 1.0)
    img_pts, msk_pts, info = [], [], {}
    for side, Pc in (('inner', cols[2]), ('outer', cols[3])):
        pa = Pc[:, :2] @ u
        o = np.argsort(pa)
        Pc, pa = Pc[o], pa[o]
        idx = np.round((pa - pa[0]) / pitch_px).astype(int)
        ym = Y4 - 0.5 * (sp.wc0 + sp.l_g) if side == 'inner' else Y4 + 0.5 * (sp.wc0 + sp.l_g)
        Ms = np.sort(P[np.abs(P[:, 1] - ym) < 5.0, 0])
        n = min(12, len(Pc))
        for j in range(n):
            img_pts.append(Pc[j, :2])
            msk_pts.append((Ms[idx[j]], ym))
        info[side] = dict(n_image=len(Pc), index_gaps=int(idx[-1] + 1 - len(idx)),
                          pitch_px_fit=float(np.polyfit(idx, pa, 1)[0]))
    img_pts, msk_pts = np.array(img_pts), np.array(msk_pts)
    M0 = T.rotated_map(k, theta)
    O = np.mean(img_pts - M0(msk_pts), axis=0)
    fmap = T.rotated_map(k, theta, O)
    dres = img_pts - fmap(msk_pts)
    r_al, r_ac = (dres @ u) / k, (dres @ np.array([-u[1], u[0]])) / k
    # widths across leg 4 and leg 2, from profiles across the legs at gap rows and post rows (dark-line centres)
    meas = dict(bar_px_aa=bar[0], bar_px_full_pixels=bar[3], px_per_um=k, um_per_px=1.0 / k, leg_direction_deg=float(np.degrees(theta)),
                columns_across_px=[float(np.mean(Pc[:, :2] @ vimg)) for Pc in cols], posts_per_column=[len(Pc) for Pc in cols],
                feature_fit=dict(n=len(img_pts), along_um=r_al.tolist(), across_um=r_ac.tolist(), rms_um=float(np.sqrt(np.mean(r_al ** 2 + r_ac ** 2)))),
                leg4=info)
    mid4 = np.mean([np.mean(cols[2][:, :2], axis=0), np.mean(cols[3][:, :2], axis=0)], axis=0)
    mid2 = np.mean([np.mean(cols[0][:, :2], axis=0), np.mean(cols[1][:, :2], axis=0)], axis=0)
    meas['post_column_separation_um'] = dict(leg4=float((np.mean(cols[3][:, :2] @ vimg) - np.mean(cols[2][:, :2] @ vimg)) / k),
                                             leg2=float((np.mean(cols[1][:, :2] @ vimg) - np.mean(cols[0][:, :2] @ vimg)) / k),
                                             mask=float(sp.wc0 + sp.l_g))
    meas['leg_centreline_separation_um'] = dict(image=float(((mid4 - mid2) @ vimg) / k), mask=float(Y4 - Y2))
    meas['pitch_um'] = dict(inner=info['inner']['pitch_px_fit'] / k, outer=info['outer']['pitch_px_fit'] / k, mask=float(pitch))
    # profile lines across leg 4 at its posts (the across-flow lines: walls and post ends)
    prof = []
    for j in range(2, 10):
        p0 = fmap(np.array([[msk_pts[j][0], Y4]]))[0]
        s, pr = T.profile_line(im, p0, vimg, 480.0 * k, step=0.25, width=1.0)
        L = T.dark_lines(pr, 10, halfwin=12, off=0) * 0.25 - 480.0 * k
        prof.append(L / k)
    meas['leg4_profile_lines_um'] = [np.round(L, 2).tolist() for L in prof]
    # the six across-flow lines of leg 4 (outer walls, posts' outer and inner ends), averaged over the profiles
    names = (('wall_in', -409.0), ('post_out_in', -109.5), ('post_in_in', -50.0), ('post_in_out', 50.0), ('post_out_out', 109.5), ('wall_out', 409.0))
    acc = {n: [] for n, _ in names}
    for L in prof:
        for n, x in names:
            v = T.line_near(np.asarray(L), x, 6.0)
            if v is not None:
                acc[n].append(v)
    a = {n: float(np.mean(v)) for n, v in acc.items() if v}
    if len(a) == 6:
        meas['leg4_widths_um'] = dict(centre=a['post_in_out'] - a['post_in_in'], side_inner=a['post_out_in'] - a['wall_in'],
                                      side_outer=a['wall_out'] - a['post_out_out'], posts=[a['post_in_in'] - a['post_out_in'], a['post_out_out'] - a['post_in_out']],
                                      outer_walls=a['wall_out'] - a['wall_in'], n_profiles={n: len(v) for n, v in acc.items()},
                                      printed=dict(centre=100.0, side_final=300.0), mask=dict(centre=float(sp.wc0), post=float(sp.l_g),
                                      outer_walls=float(sp.wc0 + 2 * sp.l_g + 600.0)))
    Mx = float(np.sort(P[:, 0])[0])
    region = (Mx - 450, Mx + 1900, Y2 - 420, Y4 + 420)
    mid = 0.5 * (Y2 + Y4)
    sep_img = meas['leg_centreline_separation_um']['image']
    spacing_from_image = abs((Y4 - Y2) - sep_img) <= 0.02 * sep_img
    # which leg an edge belongs to: by its fluid side, not by the mid line alone -- with the inset's spacing (818 µm)
    # below leg 4's full width at its end (819 µm), leg 4's inner wall lies only ~3 µm above the mid line, so a slightly
    # wider leg would cross it.  A straight edge along the legs with its fluid above it is the lower wall of the leg
    # above (leg 4), with its fluid below it the upper wall of the leg below (leg 2); any other edge (posts, leads,
    # divider ends) goes to the nearer axis.  The mid-line rule is recorded alongside, and the composite (which colours
    # by the mid line) is valid only while the two agree, which run() checks.
    def leg_of(e):
        a = np.asarray(e, float)
        d = np.diff(a, axis=0)
        n = np.c_[-d[:, 1], d[:, 0]].sum(axis=0)                     # the left normal, summed over the polyline
        fluid = n * (-e.orient if e.tag == 'channel' else e.orient)  # CHANNELS: fluid inside; ISLANDS: fluid outside
        ym = float(np.mean(a[:, 1]))
        if np.ptp(a[:, 0]) > 4.0 * max(np.ptp(a[:, 1]), 1e-9) and abs(fluid[1]) > 0.7 * float(np.hypot(*fluid)):
            return 4 if fluid[1] > 0 and Y2 < ym < Y4 + 1.0 else (2 if fluid[1] < 0 and Y2 - 1.0 < ym < Y4 else (4 if ym > mid else 2))
        return 4 if ym > mid else 2
    claim = lambda e: leg_of(e) == 4
    split = [(leg_of(e) == 4, float(np.mean(np.asarray(e)[:, 1])) > mid)
             for e in mask.edges(Mx - 450, Mx + 1900, Y2 - 420, Y4 + 420, post_max=(sp.g_p + 1.0, sp.l_g + 1.0))]
    meas['claim_split'] = dict(rule='leg 4 (the claim) by the fluid side of each straight edge along the legs; else the nearer leg axis',
                               edges=len(split), disagree_with_mid_line=int(sum(a_ != b_ for a_, b_ in split)), mid_line_um=float(mid))
    if meas['claim_split']['disagree_with_mid_line']:
        raise SystemExit('Fig. 4 claim split: %d edge(s) change leg between the fluid-side rule and the mid line; the composite '
                         '(claim_y_above) would colour them wrongly' % meas['claim_split']['disagree_with_mid_line'])
    meas['spacing_from_this_inset'] = bool(spacing_from_image)
    lab = label or os.path.basename(dxf)[:-4]
    if spacing_from_image:
        head = ('%s, Fig. 4 inset (the 300 µm bar, %.2f px; the legs at %.2f° on the image, from the post columns) with %s: leg 4\'s last '
                'post rows and the leads beyond (cyan: the reproduction claim), and the neighbouring leg at the preset\'s leg spacing, '
                'read on this inset (amber: its centre line is the image\'s by construction; its side channels\' widths there depend '
                'on which of its rows pass the inset, which the inset does not fix).' % (SRC, bar[0], np.degrees(theta), lab))
    else:
        head = ('%s, Fig. 4 inset (the 300 µm bar, %.2f px; the legs at %.2f° on the image, from the post columns) with %s: leg 4\'s last '
                'post rows and the leads beyond (cyan: the reproduction claim) and the neighbouring leg where the mask puts it (amber: '
                'the leg spacing and the legs\' lengths are the generator\'s, not printed).' % (SRC, bar[0], np.degrees(theta), lab))
    tail = ('On the image: leg 4\'s post columns %.1f µm apart (mask w_c + l_g = %.1f), pitch %.2f / %.2f µm (mask %.1f); the legs\' '
            'centrelines %.0f µm apart (mask %.0f).' % (meas['post_column_separation_um']['leg4'], sp.wc0 + sp.l_g, meas['pitch_um']['inner'],
                                                       meas['pitch_um']['outer'], pitch, meas['leg_centreline_separation_um']['image'], Y4 - Y2))
    exclude = ((735, 60, 810, 165), (390, 240, 455, 315))
    rec, cmap = run_panel('gifford_fig4_inset_leg4__' + tag, out, r, box, k, mask, region, fmap, tuple(u), pitch,
                          exclude=exclude, claim_filter=claim, caption_head=head, caption_tail=tail, bar_um=300,
                          extra=dict(measured=meas), post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=12)
    leg2_walls = [e for e in mask.edges(*region, post_max=(sp.g_p + 1.0, sp.l_g + 1.0))
                  if e.tag == 'channel' and leg_of(e) == 2 and np.ptp(np.asarray(e)[:, 0]) > 300.0]
    if leg2_walls:
        det = float(np.linalg.det(np.c_[cmap.ex, cmap.ey]))
        cls = T.by_class([cmap(e) for e in leg2_walls], ['wall_y%.0f' % float(np.mean(e[:, 1])) for e in leg2_walls], r.im, k, box,
                         tuple(u), search_px=14.0, exclude=exclude, signs=[int(e.orient * np.sign(det)) for e in leg2_walls])
        rec['measured']['leg2_side_walls'] = dict(
            rule='leg 2\'s straight channel walls (y < the legs\' mid line), their signed nearest-outline offsets (µm; + = the '
                 'image\'s outline outside the mask\'s channel: wider) within +-14 px (39 µm), by wall',
            walls={c: dict(n=v['n'], mean_signed_um=v['mean_signed_um'], sd_um=v['sd_um']) for c, v in cls.items() if v['n']})
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, box, mask, region, cmap, tuple(u), exclude=exclude, claim_filter=claim,
                                          post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=12)
    rec['mask_record'] = dict(dxf_sha256=mask.sha256, N=d['N'], leg_Y=legY, last_row_x_um=Mx)
    return rec
