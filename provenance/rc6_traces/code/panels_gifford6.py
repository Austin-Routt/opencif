"""Gifford et al. 2014, Fig. 6 (right panel: the platelet device in operation) traced with an opencif Fig. 6 mask (rc6).

What the panel shows: the downstream ends of four legs side by side (the top one's outer row is cut by the panel's
edge), every leg a pair of post rows whose last post stands one gap before a long wall -- the start of the divider
that carries the centre and side channels apart through the turn -- and the flow (the caption's arrows) running along
the rows towards the walls.  Beyond the walls three legs' channels curve away (their U-turns, nested) and the top
leg's run straight on.

Which legs: in an opencif paperclip leg j lies at Y = 0, +d, -d, +2d, ... and odd legs flow +x (layout.py), so the
downstream ends of legs 1, 3, 5, 7 lie side by side at the far end, legs 1, 3, 5 turning and leg 7 running on to the
outlet.  The identification: (1) the flow runs towards the walls (the caption's arrows point along -u: read on the
figure, not computed), so the four legs share one direction, which leaves the windows (7, 5, 3, 1) and its mirror
(1, 3, 5, 7); (2) computed: the side channels' widths at the four ends, read on the image without the mask (wall line
to the posts' outer tip line), against the drawn widths at every leg's end for both orders -- only (7, 5, 3, 1) has a
constant offset; (3) seen on the figure, not computed: the top leg's channels run straight on (the outlet), the
others' curve (the turns).  The map is then non-mirrored (mask +x along the image's -u, mask +y along +v).

Scale and rotation are the measurement's (provenance/rc6_scripts/gifford_fig6_measurement.json): the printed 125 µm
centre between the post rows' tip lines, 4.5248 µm/px, and the rows' direction on the image, 49.806°; the caption's
250 µm bar runs off the raster and only bounds the scale.  Consistency, not tests (fixed by how this image was
measured): the centre channel (it sets the scale), the pitch (V5 along) and the post length l_g (the posts' across
positions), and the layout the presets take from this image (rc6: reproduction.leg_Y and pc_last_short, read on it
without the mask by gifford6_turns.py and this module) -- the legs' axes, so the four legs' spacing and the turns' radii
and centre across the rows, and the last leg's early end.  What the trace tests: the leg identity; the side channels at
the four ends (the schedule's physics at its own scale); the leg ends -- the image's last pill stands one gap before
the divider wall, as the mask's last post does (END_CHECK on both, row by row), which fixes the along-flow pairing that
a periodic array cannot; that the generator, given those axes, draws four legs that one map registers with their
walls and divider walls; and, in the turns run, where along the flow the turns start (their centre's u, which nothing
read on the turns sets: the one independent turn result).

Translation: by features (each leg's pills indexed from its last one, which stands beside the wall: a non-periodic
feature, so the along-flow phase is fixed) and, independently across the flow, by chamfer; along the flow the chamfer
scores one pitch away about as well as the registration (a periodic array), so the along-flow position rests on the
last-pill pairing, which END_CHECK verifies on the image and on the mask.  Two runs: each leg registered on its own
(the offsets between the four maps reported), and all four on ONE map fitted on the ends of legs 5, 3 and 1
(--one-map; leg 7 drawn by that map and its own offset reported) -- the shipped presets' legs lie where the image's do
(Gifford prints no spacing; the schematic is not to scale, p. 4502), so one map holds them all; the generator's own
rule (723 µm apart for this device) would not.  Claims: each leg's posts, channels and walls up to 60 µm into the
divider wall, and with --turns the turns of legs 5, 3 and 1 as far as the frame shows them; the rest (the outlet
stretch, leg 7's end) is drawn in the no-claim colour."""
import json
import os

import numpy as np
from scipy import ndimage as ndi

import trace_lib as T
from panel_core import run_panel, scale_scan

SRC = 'Gifford et al. 2014 Lab Chip 14, 4496–4505'
FIG6_SHA = '182826e7aaccf14c'
CROP = (199, 59, 828, 541)
DU = 0.25
ARROWS = ((333, 91, 390, 152), (551, 261, 608, 324), (652, 345, 710, 408), (744, 433, 801, 496))   # the caption's flow arrows
CLAIM_PAST_WALL_UM = 60.0
LEG7_END_ROWS = 6             # one-map run: leg 7's claims end this many rows before its drawn end.  Fixed, not derived:
                              # the shipped presets draw leg 7 4 rows beyond the image's end, so the claims stop 2 rows
                              # inside it; a mask whose leg 7 ends elsewhere needs it re-set (the caption reports the end)
METRIC_PX = 2.5               # the nearest-outline window, +-11.3 µm here (panel_core.run_panel)
# the image's legs in v order: the measurement's row indices (lower-v row, higher-v row); the top leg's outer row is
# cut by the panel's edge (the measurement's rows_of finds only its inner row)
IMAGE_LEGS = (('top', (None, 0)), ('second', (1, 2)), ('third', (3, 4)), ('bottom', (5, 6)))


def frame(theta_deg):
    th = np.radians(theta_deg)
    return np.array([np.cos(th), -np.sin(th)]), np.array([np.sin(th), np.cos(th)])


def resample(im, theta):
    eu, ev = frame(theta)
    x0, y0, x1, y1 = CROP
    C = np.array([[x0, y0], [x1, y0], [x0, y1], [x1, y1]], float)
    U, V = C @ eu, C @ ev
    us = np.arange(np.floor(U.min()), U.max(), DU)
    vs = np.arange(np.floor(V.min()), V.max(), DU)
    UU, VV = np.meshgrid(us, vs)
    X, Y = UU * eu[0] + VV * ev[0], UU * eu[1] + VV * ev[1]
    ok = (X >= x0 + 1) & (X <= x1 - 1) & (Y >= y0 + 1) & (Y <= y1 - 1)
    R = ndi.map_coordinates(im, [Y, X], order=3, mode='nearest')
    bg = float(np.median(im[y0:y1, x0:x1]))
    return np.where(ok, np.clip(bg - R, 0.0, None), np.nan), us, vs


def row_pills(D, us, vs, row, a_px):
    """pill centres (u) of one row from its side outlines (the measurement's row record gives the band, u range and
    tip lines), and the row's middle line v at each"""
    v_lo, v_hi = row['band']
    m_ = (vs >= v_lo + 0.3 * (v_hi - v_lo)) & (vs <= v_hi - 0.3 * (v_hi - v_lo))
    s = np.nanmean(D[m_], axis=0)
    u0, u1 = row['u_range']
    m = (us >= u0) & (us <= u1) & ~np.isnan(s)
    uu, ss = us[m], s[m]
    pk = [uu[i] for i in range(1, len(ss) - 1) if ss[i] >= ss[i - 1] and ss[i] > ss[i + 1] and ss[i] > 0.35 * ss.max()]
    pk = np.array(pk)
    d = np.diff(pk)
    cent = np.array([0.5 * (pk[i] + pk[i + 1]) for i in range(len(d)) if abs(d[i] - a_px) < 0.25 * a_px])
    lo, hi = np.poly1d(row['apex_lo']), np.poly1d(row['apex_hi'])
    vmid = 0.5 * (lo(cent) + hi(cent))
    return cent, vmid


def to_raster(u, v, theta):
    eu, ev = frame(theta)
    u, v = np.asarray(u, float), np.asarray(v, float)
    return np.c_[u * eu[0] + v * ev[0], u * eu[1] + v * ev[1]]


def dark_peaks(prof, vv, thr):
    """centres (v) of the darkness peaks of a profile above thr, refined by a parabola"""
    out = []
    for i in range(1, len(prof) - 1):
        if prof[i] >= prof[i - 1] and prof[i] > prof[i + 1] and prof[i] > thr:
            den = prof[i - 1] - 2 * prof[i] + prof[i + 1]
            out.append(vv[i] + (0.5 * (prof[i - 1] - prof[i + 1]) / den * (vv[1] - vv[0]) if den < 0 else 0.0))
    return np.array(out)


def image_end_check(D, us, vs, row, c0, a_px, pitch_px, g_px):
    """END_CHECK on the image: downstream (-u) of the row's last pill (centre c0) there must be one outline within two
    gaps (the divider wall's end) and no further outline for at least a pitch beyond it (the wall's bright inside, not
    another pill).  Returns the gap (px) from the pill's edge to the wall's end and whether the pattern holds."""
    v_lo, v_hi = row['band']
    m_ = (vs >= v_lo + 0.3 * (v_hi - v_lo)) & (vs <= v_hi - 0.3 * (v_hi - v_lo))
    s = np.nanmean(D[m_], axis=0)
    ref = (us >= c0) & (us <= c0 + 5 * pitch_px) & ~np.isnan(s)
    thr = 0.35 * float(np.nanmax(s[ref]))
    edge = c0 - 0.5 * a_px

    def peaks_in(lo, hi):
        m = (us >= lo) & (us <= hi) & ~np.isnan(s)
        uu, ss = us[m], s[m]
        return [float(uu[i]) for i in range(1, len(ss) - 1) if ss[i] >= ss[i - 1] and ss[i] > ss[i + 1] and ss[i] > thr]
    near = peaks_in(edge - 2.0 * g_px, edge - 0.5)
    wall = max(near) if near else None
    beyond = peaks_in(edge - 2.0 * g_px - pitch_px, (wall if wall is not None else edge) - 1.0)
    return dict(pill_edge_u=float(edge), wall_end_u=wall, gap_px=float(edge - wall) if wall is not None else None,
                outlines_within_two_gaps=len(near), outlines_in_the_next_pitch=len(beyond),
                ok=bool(len(near) == 1 and not beyond))


def mask_end_check(mask, sp, y_row, x_last, direction=1):
    """END_CHECK on the mask: the divider wall (an island spanning the row's line) starts one gap G after the row's last
    post (x_last its centre; direction +1: the leg flows +x)"""
    x_edge = x_last + direction * 0.5 * sp.g_p
    best = None
    for L, p, x0, x1, y0, y1 in mask.polys:
        if L != 'ISLANDS' or not (y0 - 1e-6 <= y_row <= y1 + 1e-6) or (x1 - x0) <= sp.g_p + 1.0:
            continue
        start = x0 if direction > 0 else x1
        gap = direction * (start - x_edge)
        if 0.0 < gap <= 3.0 * sp.G and (best is None or gap < best):
            best = gap
    return dict(gap_um=best, ok=bool(best is not None and abs(best - sp.G) <= 1.0))


def image_side_widths(D, us, vs, row, side, tip_bias_px, u_win):
    """the side channel beside one row, on the image: from the row's outer tip line (corrected outwards by the tip
    estimator's bias) to the first dark line beyond it (the wall), on the mean profile over u_win; px"""
    tip_f = np.poly1d(row['apex_lo'] if side < 0 else row['apex_hi'])
    um = 0.5 * (u_win[0] + u_win[1])
    tip = float(tip_f(um)) + side * tip_bias_px
    mu = (us >= u_win[0]) & (us <= u_win[1])
    if side < 0:
        mv = (vs >= tip - 60.0) & (vs <= tip - 1.0)
    else:
        mv = (vs >= tip + 1.0) & (vs <= tip + 60.0)
    prof = np.nanmean(D[np.ix_(mv, mu)], axis=1)
    vv = vs[mv]
    if np.all(np.isnan(prof)):
        return None
    prof = np.nan_to_num(prof)
    pk = dark_peaks(prof, vv, 0.6 * prof.max())             # the wall's straight line; cells in the channel are fainter
    if not len(pk):
        return None
    wall = pk.max() if side < 0 else pk.min()              # the first line beyond the tip
    return dict(tip_v=tip, wall_v=float(wall), width_px=float(abs(wall - tip)), u_window=list(map(float, u_win)))


def leg_posts(mask, sp, Y):
    """the mask's post centres of leg Y's two rows (sorted by x): the lower row (y < Y) and the upper"""
    off = 0.5 * (float(sp.wc0) + sp.l_g)
    P = mask.posts(-1e9, 1e9, Y - off - 0.5 * sp.l_g - 5, Y + off + 0.5 * sp.l_g + 5, max_len=sp.g_p + 1.0, max_wid=sp.l_g + 1.0)
    lo = np.sort(P[np.abs(P[:, 1] - (Y - off)) < 5.0, 0])
    hi = np.sort(P[np.abs(P[:, 1] - (Y + off)) < 5.0, 0])
    return lo, hi, off


def mask_side_widths(mask, sp, Y, x_last, pitch, direction=1):
    """the drawn side channels at a leg's downstream end (µm): from the posts' outer tips to the nearest straight
    channel edge beyond, over the leg's last four rows (direction +1: the leg flows +x and ends at x_last)"""
    off = 0.5 * (float(sp.wc0) + sp.l_g)
    xa, xb = sorted((x_last - direction * 4.5 * pitch, x_last - direction * 0.5 * pitch))
    I = mask.islands
    cx, cy = (I[:, 0] + I[:, 1]) / 2, (I[:, 2] + I[:, 3]) / 2
    sel = (cx >= xa) & (cx <= xb) & ((I[:, 1] - I[:, 0]) <= sp.g_p + 1.0)
    lo_tip = float(np.mean(I[sel & (np.abs(cy - (Y - off)) < 5.0), 2]))
    hi_tip = float(np.mean(I[sel & (np.abs(cy - (Y + off)) < 5.0), 3]))
    ys = []
    for e in mask.edges(xa, xb, Y - 600.0, Y + 600.0):
        if e.tag != 'channel':
            continue
        a = np.asarray(e, float)
        if np.ptp(a[:, 0]) > 0.5 * (xb - xa) and np.ptp(a[:, 1]) < 1.0:
            ys.append(float(np.mean(a[:, 1])))
    ys = np.array(ys)
    below, above = ys[ys < lo_tip], ys[ys > hi_tip]
    return dict(lower=float(lo_tip - below.max()) if len(below) else None, upper=float(above.min() - hi_tip) if len(above) else None,
                x_window_um=[float(xa), float(xb)])


class SplitMask:
    """a trace_lib.Mask whose edges are cut at each leg's claim limit: the pieces up to it are claims, the pieces
    beyond (the turns, the outlet stretch) carry .beyond = True (drawn in the no-claim colour, left out of the metric)"""

    def __init__(self, mask, bands, outside_claim=False):
        self.mask, self.bands = mask, bands          # bands: [(y0, y1, x_cut)]
        self.outside_claim = outside_claim           # edges outside every band: claims (the turns run) or not.  In the
        # turns run that is every edge outside leg 7's band, not only the turns of legs 5, 3, 1: the panel's crop and
        # the region (to TURN_X_UM beyond the ends) keep only those turns in view, which is what makes it the turns' claim
        self.path, self.sha256 = mask.path, mask.sha256

    def edges(self, x0, x1, y0, y1, post_max=(70.0, 130.0)):
        out = []
        cuts = sorted(set([y0, y1] + [b for a in self.bands for b in a[:2] if y0 < b < y1]))
        for ya, yb in zip(cuts[:-1], cuts[1:]):
            yc = 0.5 * (ya + yb)
            xc = next((b[2] for b in self.bands if b[0] <= yc <= b[1]), None)
            if xc is None:
                parts = [(x0, x1, not self.outside_claim)]
            else:
                parts = [(x0, min(x1, xc), False), (max(x0, xc), x1, True)]
            for a, b, beyond in parts:
                if b <= a:
                    continue
                for e in self.mask.edges(a, b, ya, yb, post_max=post_max):
                    e.beyond = beyond
                    out.append(e)
        return out


def _setup(root, rasters, dxf, record, preset):
    import sys
    for p in (root,):
        if p not in sys.path:
            sys.path.insert(0, p)
    from opencif import spec as SP
    rec0 = json.load(open(record))
    mask = T.Mask(dxf)
    if mask.sha256 != rec0['files']['dxf_sha256']:
        raise ValueError('the DXF is not its record\'s')
    sp, op = SP.to_engine(SP.load(preset))
    M = json.load(open(os.path.join(root, 'provenance', 'rc6_scripts', 'gifford_fig6_measurement.json')))
    r = T.Raster(os.path.join(rasters, 'gif7-000.png'), FIG6_SHA)
    if not M['raster_sha256'].startswith(FIG6_SHA):
        raise ValueError('the measurement is not of this raster')
    return rec0, mask, sp, M, r


def image_legs(D, us, vs, M):
    """per image leg: the pills of each visible row (u, v), the side widths (px), the leg's centre line v"""
    rows = M['rows']
    bias = float(np.mean([L['bias_px']['wc'] for L in M['per_leg']])) / 2.0      # per tip, outwards
    g_px = float(M['G_read_um']) / float(M['um_per_px'])                          # the gap as the image reads it
    wc_px = float(np.mean([L['corrected_px']['wc'] for L in M['per_leg']]))
    out = []
    for name, (ilo, ihi) in IMAGE_LEGS:
        rl = rows[ilo] if ilo is not None else None
        rh = rows[ihi]
        q = dict(name=name, rows=[ilo, ihi], pills={})
        q['end_check'] = {}
        for side, rr in (('lower', rl), ('upper', rh)):
            if rr is None:
                continue
            c, vm = row_pills(D, us, vs, rr, rr['a_px'])
            q['pills'][side] = (c, vm)
            q['end_check'][side] = image_end_check(D, us, vs, rr, c[0], rr['a_px'], rr['pitch_px'], g_px)
        q['last_pill_u'] = float(np.mean([p[0][0] for p in q['pills'].values()]))
        u_hi = min(rr_['u_range'][1] for rr_ in (rl, rh) if rr_ is not None)
        u_win = (q['last_pill_u'] + 5.0, u_hi - 10.0)           # the leg's straight part, its last rows onwards
        q['side_px'] = {}
        if rl is not None:
            q['side_px']['lower'] = image_side_widths(D, us, vs, rl, -1, bias, u_win)
        q['side_px']['upper'] = image_side_widths(D, us, vs, rh, +1, bias, u_win)
        um = 0.5 * (u_win[0] + u_win[1])
        if rl is not None:
            q['centre_v'] = float(0.5 * (np.polyval(rl['apex_hi'], um) + np.polyval(rh['apex_lo'], um)))
        else:                                               # only the inner (upper) row: its inner tip line (corrected
            q['centre_v'] = float(np.polyval(rh['apex_lo'], um) - bias - 0.5 * wc_px)   # by the tip bias) less half the centre
        out.append(q)
    return out


def identify(img, ws_end, legY):
    """the image's legs (v order) against every window of four neighbouring mask legs (ordered by Y, both handednesses):
    flow direction (odd legs +x) and the side widths' offsets image - mask"""
    order = list(np.argsort(legY))                            # leg indices by Y
    cands = []
    for mirrored in (False, True):
        seq = order[::-1] if mirrored else order
        for s in range(len(seq) - len(img) + 1):
            legs = [int(seq[s + i]) for i in range(len(img))]
            same_dir = len({j % 2 for j in legs}) == 1
            offs = []
            for q, j in zip(img, legs):
                for side, w in q['side_um'].items():
                    if w is None:
                        continue
                    # mirrored: the image's lower side is the mask's upper side
                    ms = ws_end[j][('upper' if side == 'lower' else 'lower') if mirrored else side]
                    if ms is not None:
                        offs.append(w - ms)
            offs = np.array(offs)
            cands.append(dict(legs=[j + 1 for j in legs], mirrored=mirrored, same_flow_direction=same_dir,
                              offset_mean_um=float(offs.mean()), offset_sd_um=float(offs.std()), n=int(len(offs))))
    ok = [c for c in cands if c['same_flow_direction']]
    best = min(ok, key=lambda c: c['offset_sd_um'])
    return best, cands


def leg_map(q, mask_rows, theta, k, Y, off, n_max=15):
    """feature registration of one leg: its pills indexed from the last (beside the divider wall) against the mask's
    posts indexed from the leg's last; returns the map and the matched points"""
    eu, ev = frame(theta)
    img_pts, msk_pts = [], []
    for side, (c, vm) in q['pills'].items():
        Ms = mask_rows[side]
        ym = Y - off if side == 'lower' else Y + off
        d = np.diff(c)
        pp = float(np.median(d))
        idx = np.round((c - c[0]) / pp).astype(int)
        for j in range(min(n_max, len(c))):
            if idx[j] >= len(Ms):
                break
            img_pts.append(to_raster(c[j], vm[j], theta)[0])
            msk_pts.append((Ms[-1 - idx[j]], ym))
    img_pts, msk_pts = np.array(img_pts), np.array(msk_pts)
    ex, ey = -k * eu, k * ev                                   # mask +x along -u, +y along +v (non-mirrored)
    M0 = T.Map((0.0, 0.0), ex, ey)
    O = np.mean(img_pts - M0(msk_pts), axis=0)
    fmap = T.Map(O, ex, ey)
    res = img_pts - fmap(msk_pts)
    return fmap, img_pts, msk_pts, res


TURN_X_UM = 2000.0            # turns run: the region reaches this far beyond the legs' ends (the frame shows <= ~1.3 mm)


def turn_compare(D, us, vs, img, mask, fmap, theta, k, x_from, region, crop):
    """the turn walls on the image (gifford6_turns.follow_walls: legs 5, 3, 1, u <= U_ARC) against the mask's edges
    beyond the turning legs' ends mapped by fmap onto the same frame and window: (1) the nearest mask edge to every
    image wall point (µm; no window); (2) each set fitted with circles about one common centre (mask points grouped
    into walls by their radius) -- centre offset mask - image and the radii, paired wall by wall"""
    from scipy.spatial import cKDTree
    import gifford6_turns as G
    eu, ev = frame(theta)
    _, arcs = G.follow_walls(D, us, vs, img[0]['centre_v'])
    Ps = [w['points'][w['points'][:, 0] <= G.U_ARC] for w in arcs]
    (cu, cv), csd, radii, fun = G.fit_common(Ps)
    segs = [np.asarray(e, float) for e in mask.edges(x_from, region[1], region[2], region[3])]
    c0, r0, c1, r1 = crop
    u_min = min(float(P_[:, 0].min()) for P_ in Ps)
    pts, gid = [], []
    for i, sg in enumerate(segs):                                        # every edge sampled each 0.5 µm, kept by edge
        Pm, _ = T.sample_edges([sg], 0.5)
        if not len(Pm):
            continue
        Q = fmap(Pm)
        uv = np.c_[Q @ eu, Q @ ev]
        keep = ((Q[:, 0] > c0 + 2) & (Q[:, 0] < c1 - 2) & (Q[:, 1] > r0 + 2) & (Q[:, 1] < r1 - 2) &
                (uv[:, 0] <= G.U_ARC) & (uv[:, 0] >= u_min) & (uv[:, 1] >= img[0]['centre_v'] + 70.0))
        pts.append(uv[keep])
        gid.append(np.full(int(keep.sum()), i))
    uvm, gid = np.vstack(pts), np.concatenate(gid)
    tree = cKDTree(uvm)
    per_wall = []
    alld = []
    for P_, R_ in zip(Ps, radii):
        d_, _ = tree.query(P_)
        alld.append(d_)
        per_wall.append(dict(image_radius_um=float(R_ / k), n=int(len(P_)), median_um=float(np.median(d_) / k),
                             max_um=float(np.max(d_) / k)))
    alld = np.concatenate(alld) / k
    # the mask's turn walls: its edges (each a polyline; >= 40 samples, 20 µm or ~4.4 px, of it in the window), fitted
    # about one common centre;
    # pieces of one wall (equal radii within 1 µm) are then counted once
    groups = [uvm[gid == i] for i in np.unique(gid) if np.sum(gid == i) >= 40]
    cm, msd, mrad_all, mfun = G.fit_common(groups, c0=(cu, cv))
    mrad = []
    for R_ in sorted(mrad_all):
        if not mrad or R_ - mrad[-1] > 1.0 * k:
            mrad.append(R_)
    mrad = np.array(mrad)
    pairs = []
    for R_ in radii:
        i = int(np.argmin(np.abs(mrad - R_)))
        pairs.append(dict(image_um=float(R_ / k), mask_um=float(mrad[i] / k), diff_um=float((mrad[i] - R_) / k)))
    diffs = np.array([p_['diff_um'] for p_ in pairs])
    # on arcs of 8-23 degrees a common centre and the radii trade off against each other, so the walls are also compared
    # about the IMAGE's centre: each mask edge's mean distance from it, against the image wall at that radius
    rim = []
    for G_ in groups:
        rim.append(float(np.mean(np.hypot(G_[:, 0] - cu, G_[:, 1] - cv))))
    rim_u = []
    for R_ in sorted(rim):
        if not rim_u or R_ - rim_u[-1] > 1.0 * k:
            rim_u.append(R_)
    rim_u = np.array(rim_u)
    pairs_c = []
    for R_ in radii:
        i = int(np.argmin(np.abs(rim_u - R_)))
        pairs_c.append(dict(image_um=float(R_ / k), mask_um=float(rim_u[i] / k), diff_um=float((rim_u[i] - R_) / k)))
    dc = np.array([p_['diff_um'] for p_ in pairs_c])
    return dict(
        rule='the image\'s turn walls (cyan in fig_gifford2014_fig6_turns.png) against the mask\'s edges beyond the legs\' ends on '
             'the one map, on the same u window (u <= %.0f px) and inside the frame' % G.U_ARC,
        image=dict(walls=len(Ps), points=int(sum(len(P_) for P_ in Ps)), centre_uv_px=[cu, cv], centre_sd_px=csd.tolist(),
                   residual_sd_um=float(np.std(fun) / k)),
        mask=dict(walls=int(len(mrad)), edge_pieces=len(groups), points=int(len(uvm)), centre_uv_px=[float(cm[0]), float(cm[1])], centre_sd_px=msd.tolist(),
                  residual_sd_um=float(np.std(mfun) / k), radii_um=sorted((mrad / k).tolist(), reverse=True)),
        centre_offset_mask_minus_image_um=dict(u=float((cm[0] - cu) / k), v=float((cm[1] - cv) / k)),
        nearest_mask_edge_um=dict(n=int(len(alld)), median=float(np.median(alld)), p90=float(np.percentile(alld, 90)),
                                  max=float(np.max(alld))),
        per_wall=per_wall, radii_pairs=pairs,
        radii_diff_um=dict(mean=float(diffs.mean()), sd=float(diffs.std()), max_abs=float(np.abs(diffs).max()),
                           note='each set about its own fitted centre (centre and radius correlate on a short arc)'),
        radii_about_image_centre=dict(pairs=pairs_c, mean_um=float(dc.mean()), sd_um=float(dc.std()), max_abs_um=float(np.abs(dc).max()),
                                      note='the mask\'s edges\' mean distance from the IMAGE\'s turn centre, paired with the image walls'))


def run(root, rasters, out, dxf=None, record=None, preset=None, tag='design', label=None, scan=True, one_map=False,
        name='gifford2014_fig6_design', draw=True, turns=False):
    dxf = dxf or os.path.join(root, 'runs', name, name + '.dxf')
    record = record or os.path.join(root, 'runs', name, name + '.json')
    preset = preset or os.path.join(root, 'presets', name + '.yaml')
    rec0, mask, sp, M, r = _setup(root, rasters, dxf, record, preset)
    im = r.im
    d = rec0['design']
    legY = [float(v) for v in d['leg_Y']]
    theta = float(M['rows_direction_deg'])
    k = 1.0 / float(M['um_per_px'])
    eu, ev = frame(theta)
    pitch = float(sp.g_p + sp.G)
    D, us, vs = resample(im, theta)
    img = image_legs(D, us, vs, M)
    for q in img:
        q['side_um'] = {s: (w['width_px'] / k if w else None) for s, w in q['side_px'].items()}
    # the mask: every leg's rows and its end (odd legs flow +x and end at their largest x; even legs at their smallest)
    legs_m = []
    ws_end = []
    for j, Y in enumerate(legY):
        lo, hi, off = leg_posts(mask, sp, Y)
        x_last = float(max(lo.max(), hi.max())) if j % 2 == 0 else float(min(lo.min(), hi.min()))
        legs_m.append(dict(leg=j + 1, Y=Y, rows=dict(lower=lo, upper=hi), off=off, x_last=x_last))
        ws_end.append(mask_side_widths(mask, sp, Y, x_last, pitch, 1 if j % 2 == 0 else -1))
    best, cands = identify(img, ws_end, legY)
    if best['mirrored'] or best['legs'] != [7, 5, 3, 1]:
        raise ValueError('leg identification gave %s (mirrored %s), not the expected 7-5-3-1' % (best['legs'], best['mirrored']))
    # the claim limit of each leg: CLAIM_PAST_WALL_UM into the divider wall, which starts at the last post + g_p/2 + G
    bands = []
    ys4 = sorted(legs_m[j - 1]['Y'] for j in best['legs'])
    half = 0.5 * float(np.min(np.diff(ys4)))                   # half the spacing of the image's four legs
    for q, j in zip(img, best['legs']):
        L = legs_m[j - 1]
        x_cut = L['x_last'] + 0.5 * sp.g_p + sp.G + CLAIM_PAST_WALL_UM
        if one_map and j == 7:
            # on the common map the last leg's end is not a claim: the generator draws each middle leg a whole number of
            # rows and gives the last leg the remainder, so reproduction.pc_last_short moves its end only in steps of
            # about n - 1 rows, every earlier leg gaining a row (Fig. 6: 257 um, the image's 5 rows, gives 1 row); its rows up to LEG7_END_ROWS before
            # the drawn end are claims (a periodic array: their across position is the test)
            x_cut = L['x_last'] - LEG7_END_ROWS * pitch
        if turns and j != 7:
            continue                                            # turns run: legs 5, 3, 1 and their turns are claims
        bands.append((L['Y'] - abs(half), L['Y'] + abs(half), x_cut))
    smask = SplitMask(mask, bands, outside_claim=turns)
    claim = lambda e: not getattr(e, 'beyond', False)
    lab = label or os.path.basename(dxf)[:-4]
    excl = ARROWS
    recs = []
    meas_common = dict(px_per_um=k, um_per_px=1.0 / k, rows_direction_deg=theta,
                       scale_note='the printed 125 µm centre between the post rows\' tip lines (gifford_fig6_measurement.json)',
                       identification=dict(best=best, candidates=cands),
                       mask_leg_spacing_um=float(np.median(np.abs(np.diff(sorted(legY))))),
                       image_leg_spacing_um=[float((img[i + 1]['centre_v'] - img[i]['centre_v']) / k) for i in range(len(img) - 1)],
                       image_last_pill_u_px=[q['last_pill_u'] for q in img],
                       last_leg_end_offset_um=float((img[0]['last_pill_u'] - np.mean([q['last_pill_u'] for q in img[1:]])) / k),
                       mask_last_post_x_um={str(L['leg']): L['x_last'] for L in legs_m if L['leg'] % 2 == 1})
    side_tab = []
    for q, j in zip(img, best['legs']):
        for s in ('lower', 'upper'):
            if q['side_um'].get(s) is not None:
                side_tab.append(dict(image_leg=q['name'], leg=j, side=s, image_um=q['side_um'][s], mask_um=ws_end[j - 1][s]))
    meas_common['side_channels_at_leg_ends'] = side_tab
    # END_CHECK: on every visible row, the image's last pill and the mask's last post each stand one gap before a wall
    ends = []
    for q, j in zip(img, best['legs']):
        L = legs_m[j - 1]
        for side, e in q['end_check'].items():
            y_row = L['Y'] - L['off'] if side == 'lower' else L['Y'] + L['off']
            x_last_row = float(L['rows'][side].max())
            mk = mask_end_check(mask, sp, y_row, x_last_row, 1)
            ends.append(dict(leg=j, side=side, image_gap_um=(e['gap_px'] / k) if e['gap_px'] is not None else None,
                             image_ok=e['ok'], image_outlines=[e['outlines_within_two_gaps'], e['outlines_in_the_next_pitch']],
                             mask_gap_um=mk['gap_um'], mask_ok=mk['ok']))
    meas_common['end_check'] = dict(rows=ends, ok=bool(all(r['image_ok'] and r['mask_ok'] for r in ends)),
                                    rule='the last pill (image) and the last post (mask) are each followed, one gap downstream, '
                                         'by the divider wall and by no further post')
    if not one_map:
        for q, j in zip(img, best['legs']):
            L = legs_m[j - 1]
            fmap, ip, mp, res = leg_map(q, L['rows'], theta, k, L['Y'], L['off'])
            Y = L['Y']
            hw = 0.5 * float(sp.wc0) + sp.l_g + max(ws_end[j - 1][s_] or 0.0 for s_ in ('lower', 'upper')) + 40.0
            region = (L['x_last'] - 2600.0, L['x_last'] + 700.0, Y - hw, Y + hw)
            along = tuple(-eu)
            r_al, r_ac = (res @ -eu) / k, (res @ ev) / k
            meas = dict(meas_common, leg=j, image_leg=q['name'], feature_fit=dict(n=len(ip), along_um=r_al.tolist(), across_um=r_ac.tolist(),
                        rms_um=float(np.sqrt(np.mean(r_al ** 2 + r_ac ** 2)))),
                        pills_per_row={s: int(len(c)) for s, (c, _) in q['pills'].items()})
            head = ('%s, Fig. 6 (the micrograph; scale from the printed 125 µm centre, %.4f µm/px; the rows at %.2f° on the '
                    'image) with %s: the end of leg %d (the image\'s %s leg; cyan: the reproduction claim, its posts, channels and '
                    'walls to %.0f µm into the divider wall that starts one gap after the last post; amber: the turn or '
                    'outlet stretch beyond).' % (SRC, 1.0 / k, theta, lab, j, q['name'], CLAIM_PAST_WALL_UM))
            sides = '; '.join('%s side %.1f µm (mask %.1f)' % (t['side'], t['image_um'], t['mask_um']) for t in side_tab if t['leg'] == j)
            ec = [r for r in meas_common['end_check']['rows'] if r['leg'] == j]
            tail = ('On the image: side channels at the leg\'s end %s.  END_CHECK (the along-flow pairing): on %s the last pill '
                    'stands %s µm and the mask\'s last post %.1f µm before the divider wall, no further post beyond (%s).  The '
                    'centre (the scale), the pitch and the post length were measured on this micrograph: their agreement is '
                    'consistency.' % (sides, ' and '.join('the %s row' % r['side'] for r in ec),
                                      ' / '.join('%.1f' % r['image_gap_um'] for r in ec if r['image_gap_um'] is not None),
                                      ec[0]['mask_gap_um'] or float('nan'),
                                      'holds' if all(r['image_ok'] and r['mask_ok'] for r in ec) else 'FAILS'))
            rec, cmap = run_panel('gifford_fig6_leg%d__%s' % (j, tag), out, r, CROP, k, smask, region, fmap, along, pitch,
                                  exclude=excl, claim_filter=claim, caption_head=head, caption_tail=tail, bar_um=250,
                                  extra=dict(measured=meas), post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=12,
                                  scale_src='scale and rotation from the printed centre width and the rows on the image', draw=draw,
                                  metric_search_px=METRIC_PX)
            if scan:
                rec['V5_scale_scan'] = scale_scan(r, CROP, smask, region, cmap, along, exclude=excl, claim_filter=claim,
                                                  post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=12)
            rec['mask_record'] = dict(dxf_sha256=mask.sha256, N=d['N'], leg_Y=legY, turns=d.get('turns'))
            recs.append((rec, cmap))
        # the four maps against one another: where leg j's map puts the mask point (x_last_j, Y_j), against leg 5's map
        ref = recs[1][1]
        for (rec, cmap), j in zip(recs, best['legs']):
            L = legs_m[j - 1]
            p = np.array([[L['x_last'], L['Y']]])
            dp = cmap(p)[0] - ref(p)[0]
            rec['measured']['offset_from_leg5_map_um'] = dict(along=float(dp @ -eu / k), across=float(dp @ ev / k))
        if draw:
            # the whole panel: every leg on its own map (cyan claims, amber beyond)
            segs, cols = [], []
            for (rec, cmap), j in zip(recs, best['legs']):
                for e in smask.edges(*rec['region_um']):
                    segs.append(cmap(e))
                    cols.append(T.TRACE_COLOUR if claim(e) else T.NOCLAIM_COLOUR)
            cap = ('%s, Fig. 6 (the micrograph; %.4f µm/px from the printed 125 µm centre; rows at %.2f°) with %s: the ends '
                   'of legs 7, 5, 3, 1, each registered on its own (the mask\'s legs are %.0f µm apart, the image\'s %.0f; '
                   'Gifford prints no spacing).  Cyan: the claims; amber: the turns and the outlet stretch.  Nearest-outline '
                   'median per leg: %s µm; side channels at the ends, image / mask: %s µm.' % (
                       SRC, 1.0 / k, theta, lab, meas_common['mask_leg_spacing_um'], float(np.mean(meas_common['image_leg_spacing_um'])),
                       ' / '.join('%.2f' % r_['metric']['median_abs_um'] for r_, _ in recs),
                       ', '.join('%.0f / %.0f' % (t['image_um'], t['mask_um']) for t in side_tab)))
            T.draw_panel(os.path.join(out, 'gifford_fig6_legs__%s' % tag), im, CROP, segs, k, 250, cap, cols,
                         clip=T.micrograph_box(im, CROP))
        return [r_ for r_, _ in recs]
    # one map for the four legs: features of legs 5, 3, 1 (their ends side by side on the image) fixed together; leg 7
    # drawn by the same map and its offset reported
    ips, mps = [], []
    for q, j in zip(img, best['legs']):
        if j == 7:
            continue
        L = legs_m[j - 1]
        _, ip, mp, _ = leg_map(q, L['rows'], theta, k, L['Y'], L['off'])
        ips.append(ip)
        mps.append(mp)
    ip, mp = np.vstack(ips), np.vstack(mps)
    ex, ey = -k * eu, k * ev
    O = np.mean(ip - T.Map((0.0, 0.0), ex, ey)(mp), axis=0)
    fmap = T.Map(O, ex, ey)
    res = ip - fmap(mp)
    r_al, r_ac = (res @ -eu) / k, (res @ ev) / k
    q7, L7 = img[0], legs_m[6]
    f7, ip7, mp7, _ = leg_map(q7, L7['rows'], theta, k, L7['Y'], L7['off'])
    d7 = f7(np.array([[L7['x_last'], L7['Y']]]))[0] - fmap(np.array([[L7['x_last'], L7['Y']]]))[0]
    ys = [legs_m[j - 1]['Y'] for j in best['legs']]
    hw = 0.5 * float(sp.wc0) + sp.l_g + 200.0
    xl = min(legs_m[j - 1]['x_last'] for j in best['legs'])
    region = (xl - 2600.0, xl + (TURN_X_UM if turns else 700.0), min(ys) - hw, max(ys) + hw)
    along = tuple(-eu)
    tc = None
    if turns:
        x_from = max(legs_m[j - 1]['x_last'] for j in best['legs'] if j != 7) + 0.5 * sp.g_p + sp.G
        tc = turn_compare(D, us, vs, img, mask, fmap, theta, k, x_from, region, CROP)
    meas = dict(meas_common, legs=best['legs'], feature_fit=dict(n=len(ip), legs=[5, 3, 1], along_um=r_al.tolist(), across_um=r_ac.tolist(),
                rms_um=float(np.sqrt(np.mean(r_al ** 2 + r_ac ** 2)))),
                leg7_own_map_minus_common_um=dict(along=float(d7 @ -eu / k), across=float(d7 @ ev / k)))
    if tc is not None:
        meas['turns'] = tc
    head = ('%s, Fig. 6 (the micrograph; scale from the printed 125 µm centre, %.4f µm/px; the rows at %.2f° on the image) '
            'with %s, all four legs on ONE map, fitted on the ends of legs 5, 3 and 1 (cyan: the claim, %s leg 7 to %d rows '
            'before its drawn end; amber: %s).' % (
                SRC, 1.0 / k, theta, lab,
                'legs 5, 3 and 1 with their divider walls and their turns as far as the frame shows them,' if turns else
                'each of those legs to %.0f µm into its divider wall,' % CLAIM_PAST_WALL_UM,
                LEG7_END_ROWS, 'the outlet stretch and leg 7\'s end' if turns else 'the turns, the outlet stretch and leg 7\'s end'))
    meas_axes = None
    tj = os.path.join(root, 'provenance', 'rc6_traces', 'measurements', 'gifford_fig6_turns.json')
    if os.path.exists(tj):
        meas_axes = json.load(open(tj))['leg_Y_um']
    axes_from_turns = meas_axes is not None and len(meas_axes) == len(legY) and max(abs(a_ - b_) for a_, b_ in zip(meas_axes, legY)) < 1.0
    if turns and not axes_from_turns:
        # a mask whose legs 2, 4, 6 are not where the image's turns put them (the first rc6 form: uniform legs): the turns
        # are then a test of its layout, and its turns are not concentric with the image's
        co = tc['centre_offset_mask_minus_image_um']
        tail = ('This mask\'s legs 2, 4, 6 lie where the generator\'s uniform paperclip puts them (%.0f µm apart), not where the '
                'image\'s turns put them, so the turns test that layout: the image\'s %d turn walls lie a median %.1f µm (90th '
                'percentile %.1f, max %.0f) from the mask\'s nearest edge, and the mask\'s own turn centre lies %.0f µm across from the '
                'image\'s: the turns are not concentric with the image\'s (a radius comparison about one centre does not apply).  '
                'Leg 7, not in the fit, lies %.2f µm across from where its own end registers it.' % (
                    meas_common['mask_leg_spacing_um'], tc['image']['walls'], tc['nearest_mask_edge_um']['median'],
                    tc['nearest_mask_edge_um']['p90'], tc['nearest_mask_edge_um']['max'], co['v'],
                    meas['leg7_own_map_minus_common_um']['across']))
    elif turns:
        du = tc['centre_offset_mask_minus_image_um']['u']            # the flow runs along -u: mask - image < 0 is downstream
        tail = ('The leg axes are the preset\'s (reproduction.leg_Y), read on this image without the mask (legs 1, 3, 5, 7 on the '
                'post rows; legs 2, 4, 6 from the turns\' common centre, a semicircle assumed), so the turns\' radii and their '
                'centre across the rows agree with the image by construction (consistency, not a test): the image\'s %d walls lie '
                'a median %.1f µm (90th percentile %.1f) from the mask\'s nearest edge, and about the image\'s turn centre the '
                'mask\'s walls lie at radii %+.1f ± %.1f µm from the image\'s (the side channels\' own offset included).  The '
                'independent turn result is along the flow, which nothing read on the turns sets: fitted alone, the mask\'s turn '
                'centre lies %.1f µm %s of the image\'s (and %.1f µm across).  Leg 7, not in the fit, lies %.2f µm '
                'across from where its own end registers it; its drawn end is %.0f µm (%d rows) beyond the image\'s (pc_last_short '
                'moves it in steps of about n - 1 rows: LIMITATIONS 25).' % (
                    tc['image']['walls'], tc['nearest_mask_edge_um']['median'], tc['nearest_mask_edge_um']['p90'],
                    tc['radii_about_image_centre']['mean_um'], tc['radii_about_image_centre']['sd_um'],
                    abs(du), 'downstream' if du < 0 else 'upstream', tc['centre_offset_mask_minus_image_um']['v'],
                    meas['leg7_own_map_minus_common_um']['across'], -meas['leg7_own_map_minus_common_um']['along'],
                    int(round(-meas['leg7_own_map_minus_common_um']['along'] / pitch))))
    else:
        tail = ('This mask\'s legs 1, 3, 5, 7 lie %.0f µm apart, the spacing read on this image without the mask; leg 7, not in '
                'the fit, lies %.2f µm across from where its own end registers it.  Its drawn end is %.0f µm (%d rows) beyond the image\'s: '
                'with pc_last_short 257 µm (the image\'s 5 rows) the generator ends it 1 row short of the other legs, not 5 (each middle '
                'leg is drawn a whole number of rows and the last leg takes the remainder).' % (
                    meas_common['mask_leg_spacing_um'], meas['leg7_own_map_minus_common_um']['across'], -meas['leg7_own_map_minus_common_um']['along'],
                    int(round(-meas['leg7_own_map_minus_common_um']['along'] / pitch))))
    rec, cmap = run_panel('gifford_fig6_four_legs__%s' % tag, out, r, CROP, k, smask, region, fmap, along, pitch,
                          exclude=excl, claim_filter=claim, caption_head=head, caption_tail=tail, bar_um=250,
                          extra=dict(measured=meas), post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=12,
                          scale_src='scale and rotation from the printed centre width and the rows on the image', draw=draw,
                                  metric_search_px=METRIC_PX)
    if scan:
        rec['V5_scale_scan'] = scale_scan(r, CROP, smask, region, cmap, along, exclude=excl, claim_filter=claim,
                                          post_max=(sp.g_p + 1.0, sp.l_g + 1.0), min_component=12)
    rec['mask_record'] = dict(dxf_sha256=mask.sha256, N=d['N'], leg_Y=legY, turns=d.get('turns'))
    return [rec]
